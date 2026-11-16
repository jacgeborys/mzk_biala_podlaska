"""Convert the MZK Biała Podlaska timetable export (data/raw/1488, semicolon CSV,
CP1250) into a GTFS feed in data/gtfs/ (+ data/gtfs.zip).

Raw -> GTFS:
  przystanki.csv            -> stops.txt (posts = stops, NR_PRZYST groups = parent stations)
  rozklad.csv               -> trips.txt, stop_times.txt (stop order inferred from times)
  warianty_tras.csv         -> trip headsigns, direction_id, route long names
  typy_dni + kalendarz.csv  -> calendar_dates.txt (service_id = day type)
  OSM road network          -> shapes.txt (bus route per variant, map-matched with Viterbi)

Usage:  python convert_to_gtfs.py [--no-shapes]
"""
import argparse
import itertools
import re
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw" / "1488"
OUT = ROOT / "data" / "gtfs"
CACHE = ROOT / "data" / "cache"
QA = ROOT / "data" / "qa"

AGENCY = dict(agency_id="MZK", agency_name="MZK Biała Podlaska",
              agency_url="https://mzkbp.pl/", agency_timezone="Europe/Warsaw",
              agency_lang="pl")
METRIC_CRS = "EPSG:2180"

# Routing cost = length x class weight: mild preference for main roads, so shapes avoid
# shortcuts through side streets that buses do not take. Living streets are never used.
ROAD_WEIGHTS = {
    "motorway": 1.0, "motorway_link": 1.0, "trunk": 1.0, "trunk_link": 1.0,
    "primary": 1.0, "primary_link": 1.0, "secondary": 1.0, "secondary_link": 1.0,
    "tertiary": 1.05, "tertiary_link": 1.05,
    "unclassified": 1.1, "residential": 1.1,
    "service": 1.5,
}
DEFAULT_ROAD_WEIGHT = 1.1
EXCLUDED_ROADS = {"living_street"}

# One colour per line, distinct hues (text colour picked for contrast below)
LINE_COLORS = {
    "A": "E6194B", "B": "3CB44B", "C": "4363D8", "D": "F58231", "E": "911EB4",
    "F": "42D4F4", "G": "F032E6", "H": "9A6324", "I": "469990", "K": "800000",
    "L": "808000", "M": "000075", "T": "A9A9A9",
}


def read(name):
    df = pd.read_csv(RAW / name, sep=";", encoding="cp1250", dtype=str, keep_default_na=False)
    return df.apply(lambda c: c.str.strip())


def clean_name(s):
    return re.sub(r"\s+", " ", s).strip()


def text_color(hex_):
    r, g, b = (int(hex_[i:i + 2], 16) for i in (0, 2, 4))
    return "000000" if 0.299 * r + 0.587 * g + 0.114 * b > 150 else "FFFFFF"


# ---------------------------------------------------------------- stops
def build_stops(prz):
    prz = prz.copy()
    prz["stop_lat"] = prz.POZ_Y.astype(int) / 1e5
    prz["stop_lon"] = prz.POZ_X.astype(int) / 1e5
    prz["stop_name"] = prz.NAZWA.map(clean_name)

    children = pd.DataFrame({
        "stop_id": prz.ID_PRZYST,
        "stop_code": prz.NR_PRZYST + prz.NR_SLUPKA.str.zfill(2),
        "stop_name": prz.stop_name,
        "stop_lat": prz.stop_lat.round(6),
        "stop_lon": prz.stop_lon.round(6),
        "location_type": 0,
        "parent_station": "P" + prz.NR_PRZYST,
    })
    parents = (prz.groupby("NR_PRZYST")
               .agg(stop_name=("stop_name", lambda s: s.mode().iloc[0]),
                    stop_lat=("stop_lat", "mean"), stop_lon=("stop_lon", "mean"))
               .reset_index())
    parents = pd.DataFrame({
        "stop_id": "P" + parents.NR_PRZYST, "stop_code": parents.NR_PRZYST,
        "stop_name": parents.stop_name, "stop_lat": parents.stop_lat.round(6),
        "stop_lon": parents.stop_lon.round(6), "location_type": 1, "parent_station": "",
    })

    # sanity: posts far from their station centroid
    m = children.merge(parents[["stop_id", "stop_lat", "stop_lon"]],
                       left_on="parent_station", right_on="stop_id", suffixes=("", "_p"))
    d = haversine(m.stop_lat, m.stop_lon, m.stop_lat_p, m.stop_lon_p)
    far = m[d > 300]
    if len(far):
        print(f"  warning: {len(far)} posts > 300 m from their station centroid:",
              ", ".join(far.stop_name.unique()))
    return pd.concat([parents, children], ignore_index=True)


def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371000 * 2 * np.arcsin(np.sqrt(a))


# ---------------------------------------------------------------- stop order
def variant_orders(roz, xy):
    """Stop sequence per (line, variant). Running times are identical for all trips
    of a variant, so stops are sorted by their offset from the trip start; stops that
    share a minute are ordered by the permutation giving the shortest path between
    the neighbouring stops."""
    roz = roz.copy()
    roz["off"] = roz.m - roz.groupby("ID_KRS").m.transform("min")
    off = roz.groupby(["LINIA", "WARIANT", "ID_PRZYST"]).off.agg(["min", "max"])
    assert (off["min"] == off["max"]).all(), "running times differ between trips of a variant"
    orders, resolved = {}, 0
    for key, g in off["min"].reset_index().groupby(["LINIA", "WARIANT"]):
        g = g.sort_values(["min", "ID_PRZYST"])
        seq = []
        groups = [list(x.ID_PRZYST) for _, x in g.groupby("min", sort=True)]
        for i, grp in enumerate(groups):
            if len(grp) > 1:
                prev = seq[-1] if seq else None
                nxt = groups[i + 1][0] if i + 1 < len(groups) and len(groups[i + 1]) == 1 else None

                def cost(p):
                    pts = ([prev] if prev else []) + list(p) + ([nxt] if nxt else [])
                    return sum(np.hypot(*(xy[a] - xy[b])) for a, b in zip(pts, pts[1:]))
                grp = list(min(itertools.permutations(grp), key=cost))
                resolved += 1
            seq.extend(grp)
        orders[key] = seq
    print(f"  {len(orders)} variants, {resolved} same-minute groups ordered by geometry")
    return orders


# ---------------------------------------------------------------- shapes
class Router:
    """Map-matches stop sequences onto the OSM drive network (Viterbi over the
    candidate edges near each stop, transition cost = network distance)."""

    def __init__(self, bbox):
        import osmnx as ox
        from shapely.strtree import STRtree
        CACHE.mkdir(parents=True, exist_ok=True)
        ox.settings.cache_folder = str(CACHE / "overpass")
        path = CACHE / "osm_drive_service.graphml"
        if path.exists():
            G = ox.load_graphml(path)
        else:
            print("  downloading OSM road network ...")
            G = ox.graph_from_bbox(bbox=bbox, network_type="drive_service", simplify=True)
            ox.save_graphml(G, path)
        G = ox.project_graph(G, to_crs=METRIC_CRS)
        edges = ox.graph_to_gdfs(G, nodes=False).reset_index()
        nodes = ox.graph_to_gdfs(G, edges=False)
        self.G = G
        from shapely.geometry import LineString
        geoms = []
        for u, v, geom in zip(edges.u, edges.v, edges.geometry):
            if geom is None:
                geom = LineString([nodes.loc[u].geometry, nodes.loc[v].geometry])
            # orient u -> v
            pu = nodes.loc[u].geometry
            if geom.interpolate(0).distance(pu) > geom.interpolate(geom.length).distance(pu):
                geom = LineString(list(geom.coords)[::-1])
            geoms.append(geom)
        edges["geometry"] = geoms
        # merged (simplified) edges carry a list of classes: drop if any part is excluded,
        # otherwise weight by the worst class
        classes = edges.highway.map(lambda h: h if isinstance(h, list) else [h])
        edges = edges[~classes.map(lambda c: bool(EXCLUDED_ROADS & set(c)))].reset_index(drop=True)
        classes = edges.highway.map(lambda h: h if isinstance(h, list) else [h])
        penalty = classes.map(lambda c: max(ROAD_WEIGHTS.get(x, DEFAULT_ROAD_WEIGHT) for x in c))
        edges["w"] = edges.geometry.length * penalty.values
        self.edges = edges
        # best parallel edge per (u, v) for path reconstruction
        best = edges.sort_values("w").drop_duplicates(["u", "v"])
        self.uv = {(u, v): (w, g) for u, v, w, g in zip(best.u, best.v, best.w, best.geometry)}
        import networkx as nx
        self.D = nx.DiGraph()
        for (u, v), (w, _) in self.uv.items():
            self.D.add_edge(u, v, w=w)
        self.tree = STRtree(list(edges.geometry))
        self.sp_cache = {}

    def candidates(self, pt, k=6, radius=60):
        idx = self.tree.query(pt.buffer(radius))
        if len(idx) == 0:
            idx = [self.tree.nearest(pt)]
        c = []
        for i in idx:
            g = self.edges.geometry.iat[i]
            c.append((g.distance(pt), i, g.project(pt)))
        return sorted(c)[:k]

    def sp(self, src):
        import networkx as nx
        if src not in self.sp_cache:
            self.sp_cache[src] = nx.dijkstra_predecessor_and_distance(self.D, src, cutoff=10000, weight="w")
        return self.sp_cache[src]

    def path_geom(self, src, dst):
        pred, _ = self.sp(src)
        nodes = [dst]
        while nodes[-1] != src:
            nodes.append(pred[nodes[-1]][0])
        nodes.reverse()
        coords = []
        for a, b in zip(nodes, nodes[1:]):
            coords.extend(self.uv[(a, b)][1].coords)
        return coords

    def match(self, pts):
        """pts: list of shapely Points (metric). Returns (coords, dist_at_stop, n_fallback)."""
        from shapely.ops import substring
        E = self.edges
        cands = [self.candidates(p) for p in pts]
        SNAP = 3.0
        score = [d * SNAP for d, _, _ in cands[0]]
        back = []
        for i in range(1, len(pts)):
            new, bk = [], []
            for d, ei, pos in cands[i]:
                best, arg = np.inf, None
                for j, (dp, ej, pp) in enumerate(cands[i - 1]):
                    if ej == ei and pos >= pp:
                        t = pos - pp
                    else:
                        _, dist = self.sp(E.v.iat[ej])
                        net = dist.get(E.u.iat[ei], np.inf)
                        t = (E.geometry.iat[ej].length - pp) + net + pos
                    if score[j] + t < best:
                        best, arg = score[j] + t, j
                if arg is None:  # unreachable from every previous candidate
                    j = int(np.argmin(score))
                    best, arg = score[j] + 5 * pts[i].distance(pts[i - 1]) + 1e6, j
                new.append(best + d * SNAP)
                bk.append(arg)
            score, back = new, back + [bk]
        # backtrack
        sel = [int(np.argmin(score))]
        for bk in reversed(back):
            sel.append(bk[sel[-1]])
        sel.reverse()
        chosen = [cands[i][s] for i, s in enumerate(sel)]

        coords, dists, fallback = [], [0.0], 0
        for (_, ej, pp), (_, ei, pos), pa, pb in zip(chosen, chosen[1:], pts, pts[1:]):
            ga, gb = E.geometry.iat[ej], E.geometry.iat[ei]
            if ej == ei and pos >= pp:
                leg = list(substring(ga, pp, pos).coords)
            else:
                src, dst = E.v.iat[ej], E.u.iat[ei]
                if dst in self.sp(src)[1]:
                    leg = (list(substring(ga, pp, ga.length).coords) + self.path_geom(src, dst)
                           + list(substring(gb, 0, pos).coords))
                else:
                    fallback += 1
                    leg = [ga.interpolate(pp).coords[0], gb.interpolate(pos).coords[0]]
            leg = dedupe(leg)
            if coords and leg and coords[-1] == leg[0]:
                leg = leg[1:]
            seg_len = polyline_len(([coords[-1]] if coords else []) + leg)
            coords.extend(leg)
            dists.append(dists[-1] + seg_len)
        return coords, dists, fallback


def dedupe(c):
    out = []
    for p in c:
        if not out or p != out[-1]:
            out.append(p)
    return out


def polyline_len(c):
    if len(c) < 2:
        return 0.0
    a = np.asarray(c)
    return float(np.hypot(*np.diff(a, axis=0).T).sum())


def build_shapes(orders, prz):
    from pyproj import Transformer
    from shapely.geometry import Point
    to_m = Transformer.from_crs("EPSG:4326", METRIC_CRS, always_xy=True)
    to_ll = Transformer.from_crs(METRIC_CRS, "EPSG:4326", always_xy=True)
    lon = prz.POZ_X.astype(int) / 1e5
    lat = prz.POZ_Y.astype(int) / 1e5
    x, y = to_m.transform(lon.values, lat.values)
    pt = {sid: Point(a, b) for sid, a, b in zip(prz.ID_PRZYST, x, y)}
    pad = 0.02
    router = Router((lat.max() + pad, lat.min() - pad, lon.max() + pad, lon.min() - pad))

    shapes, stop_dist, qa = [], {}, []
    for (line, var), seq in orders.items():
        sid = f"{line}_{var}"
        coords, dists, fb = router.match([pt[s] for s in seq])
        straight = sum(pt[a].distance(pt[b]) for a, b in zip(seq, seq[1:]))
        for i, (a, b) in enumerate(zip(seq, seq[1:])):
            leg = dists[i + 1] - dists[i]
            crow = pt[a].distance(pt[b])
            qa.append(dict(shape_id=sid, from_stop=a, to_stop=b, crow_m=round(crow),
                           route_m=round(leg), ratio=round(leg / crow, 2) if crow > 0 else np.nan))
        lo, la = to_ll.transform(*np.asarray(coords).T)
        cum = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(np.asarray(coords), axis=0).T))])
        shapes.append(pd.DataFrame({
            "shape_id": sid, "shape_pt_lat": np.round(la, 6), "shape_pt_lon": np.round(lo, 6),
            "shape_pt_sequence": np.arange(len(coords)), "shape_dist_traveled": np.round(cum, 1)}))
        stop_dist[(line, var)] = np.round(dists, 1)
        if fb:
            print(f"  warning: {sid}: {fb} legs not routable, straight line used")
        if straight and dists[-1] / straight > 1.6:
            print(f"  note: {sid} routed {dists[-1]/straight:.2f}x the stop-to-stop straight line")
    QA.mkdir(parents=True, exist_ok=True)
    qa = pd.DataFrame(qa).sort_values("ratio", ascending=False)
    qa.to_csv(QA / "shape_legs.csv", index=False)
    sus = qa[(qa.ratio > 2.5) & (qa.route_m > 400)]
    print(f"  shapes: {len(shapes)}, suspicious legs (ratio>2.5, >400 m): {len(sus)} -> data/qa/shape_legs.csv")
    return pd.concat(shapes, ignore_index=True), stop_dist


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-shapes", action="store_true")
    args = ap.parse_args()

    prz, roz, war = read("przystanki.csv"), read("rozklad.csv"), read("warianty_tras.csv")
    kal, typy = read("kalendarz.csv"), read("typy_dni.csv")
    roz["m"] = roz.ODJAZD.str[:2].astype(int) * 60 + roz.ODJAZD.str[3:5].astype(int)

    assert set(roz.ID_PRZYST) <= set(prz.ID_PRZYST)
    tr = roz.drop_duplicates("ID_KRS")
    assert roz.groupby("ID_KRS")[["LINIA", "WARIANT", "TYP_DNIA"]].nunique().max().max() == 1
    var_keys = set(zip(war.LINIA, war.WARIANT))
    missing = set(zip(tr.LINIA, tr.WARIANT)) - var_keys
    assert not missing, f"variants missing in warianty_tras: {missing}"

    print("stops")
    stops = build_stops(prz)

    print("stop order")
    xy = {s: np.array([float(a) * 0.6157, float(b)]) for s, a, b in zip(prz.ID_PRZYST, prz.POZ_X, prz.POZ_Y)}
    orders = variant_orders(roz, xy)  # x scaled by cos(52°) so distances are isotropic

    # routes
    routes = []
    for line in sorted(roz.LINIA.unique()):
        w = war[war.LINIA == line]
        used = tr[tr.LINIA == line].WARIANT.value_counts()
        base = w[w.PODSTAWOWY == "Y"]
        # basic variant if flagged, else the variant with most trips
        main_var = base.iloc[0] if len(base) else w[w.WARIANT == used.index[0]].iloc[0]
        long_name = clean_name(main_var.OPIS_KIER.replace("->", " – "))
        col = LINE_COLORS.get(line, "777777")
        routes.append(dict(route_id=line, agency_id=AGENCY["agency_id"], route_short_name=line,
                           route_long_name=long_name, route_type=3,
                           route_color=col, route_text_color=text_color(col)))
    routes = pd.DataFrame(routes)

    # trips
    wv = war.set_index(["LINIA", "WARIANT"])
    trips = pd.DataFrame({
        "route_id": tr.LINIA.values,
        "service_id": tr.TYP_DNIA.values,
        "trip_id": tr.ID_KRS.values,
        "trip_headsign": [wv.loc[(l, v)].OPIS_TABL for l, v in zip(tr.LINIA, tr.WARIANT)],
        "direction_id": [0 if wv.loc[(l, v)].KIERUNEK == "T" else 1 for l, v in zip(tr.LINIA, tr.WARIANT)],
        "shape_id": (tr.LINIA + "_" + tr.WARIANT).values,
    }).sort_values(["route_id", "service_id", "trip_id"])

    # stop_times
    seq_idx = {(l, v, s): i for (l, v), seq in orders.items() for i, s in enumerate(seq)}
    roz["stop_sequence"] = [seq_idx[k] for k in zip(roz.LINIA, roz.WARIANT, roz.ID_PRZYST)]
    roz = roz.sort_values(["ID_KRS", "stop_sequence"])
    assert (roz.groupby("ID_KRS").m.diff().fillna(0) >= 0).all()
    t = roz.ODJAZD + ":00"
    st = pd.DataFrame({"trip_id": roz.ID_KRS, "arrival_time": t, "departure_time": t,
                       "stop_id": roz.ID_PRZYST, "stop_sequence": roz.stop_sequence, "timepoint": 1})

    # calendar_dates
    cal = pd.DataFrame({"service_id": kal.TYP_DNIA, "date": kal.DATA.str.replace("-", ""),
                        "exception_type": 1})
    unknown = set(cal.service_id) - set(typy.TYP_DNIA)
    assert not unknown, unknown

    feed_info = pd.DataFrame([dict(
        feed_publisher_name="Unofficial conversion of MZK Biała Podlaska timetable export",
        feed_publisher_url=AGENCY["agency_url"], feed_lang="pl", feed_contact_url=AGENCY["agency_url"],
        feed_start_date=cal.date.min(), feed_end_date=cal.date.max(),
        feed_version=f"1488_{cal.date.min()}")])

    shapes = None
    if not args.no_shapes:
        print("shapes")
        shapes, stop_dist = build_shapes(orders, prz)
        st["shape_dist_traveled"] = [stop_dist[(l, v)][i] for l, v, i in
                                     zip(roz.LINIA, roz.WARIANT, roz.stop_sequence)]
    else:
        trips = trips.drop(columns="shape_id")

    # write
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    files = {"agency.txt": pd.DataFrame([AGENCY]), "stops.txt": stops, "routes.txt": routes,
             "trips.txt": trips, "stop_times.txt": st, "calendar_dates.txt": cal,
             "feed_info.txt": feed_info}
    if shapes is not None:
        files["shapes.txt"] = shapes
    for name, df in files.items():
        df.to_csv(OUT / name, index=False, encoding="utf-8", lineterminator="\n")
    zpath = OUT.parent / "gtfs.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for name in files:
            z.write(OUT / name, name)
    print(f"written {OUT} and {zpath}")
    for name, df in files.items():
        print(f"  {name:20s} {len(df):7d} rows")


if __name__ == "__main__":
    main()
