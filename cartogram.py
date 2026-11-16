"""Flow cartogram: every street piece drawn as thick as the number of buses passing it in one day.

All shapes are noded into one network; each piece counts every pass of a shape running parallel
to it within MERGE_M (either direction), or within DUAL_M on the same named OSM street (a loop
shape using both carriageways, or an out-and-back spur, counts twice), so both
directions are summed and both carriageways of a divided road / both halves of a split around a
traffic island carry the street's full count; such pieces get a filler copy shifted halfway towards
the other carriageway, so a divided road reads as one band. Beyond 2 m, parallel passes count only
on the same street (OSM ref, else name), so close parallel streets (Chłodna / Twarda) stay apart. Layout follows the TransitFrequency posters (Bahnschrift header, frame, scale
bar, legend box, credit line).

Output: png/cartogram_<date>.png, layer `segments` in data/gis/biala_podlaska_gtfs.gpkg
        (one LineString per run of pieces with the same count: trips, lines)

Usage:  python cartogram.py [--date 20261118] [--city]   (--city: frame on the city limits, drawn dashed)
"""
import argparse
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
from pyproj import Transformer
from shapely.geometry import LineString, Point, box
from shapely.affinity import translate
from shapely.ops import linemerge

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager, patheffects  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

ROOT = Path(__file__).resolve().parent
GTFS = ROOT / "data" / "gtfs"
BASEMAP = ROOT / "data" / "basemap"
GPKG = ROOT / "data" / "gis" / "biala_podlaska_gtfs.gpkg"
LIMITS = ROOT / "data" / "gis" / "city_limits.gpkg"  # OSM relation 3016433 (Nominatim)
CRS = "EPSG:2180"
FRAME_ALL = (773970, 466100, 790310, 476690)   # all stops + 300 m (as animate_buses.py)
MERGE_M = 20          # shapes this close and parallel count on the same piece (traffic islands)
DUAL_M = 45           # ... and up to this far when both lie on the same named OSM street (dual carriageways)
PARALLEL_DEG = 30
SAMPLE_M = 5           # shapes sampled this often to count their passes
SCALE_PCTL = 99.5      # widths / colours top out at this length percentile of counts (clipped above)
MAX_WIDTH_PT = 16
DAY_NAMES = ["poniedziałek", "wtorek", "środa", "czwartek", "piątek", "sobota", "niedziela"]

# TransitFrequency poster look
FONT = Path("C:/Windows/Fonts/bahnschrift.ttf")
INK = "#e1e1e8"
INK_DIM = "#9a9aa6"
PAGE = "#000000"
PAGE_W_IN = 19.2
MARGIN_IN = 0.30


def bahn(size):
    if FONT.exists():
        return font_manager.FontProperties(fname=str(FONT), size=size)
    return font_manager.FontProperties(family="sans-serif", size=size)


def heading(g, d):
    a = g.interpolate(max(0.0, d - 5))
    b = g.interpolate(min(g.length, d + 5))
    return np.degrees(np.arctan2(b.y - a.y, b.x - a.x))


def piece_counts(date):
    """Pieces of the noded union of all shapes, each with the trips of every shape running on it
    or parallel to it within MERGE_M."""
    from shapely.ops import unary_union
    from shapely.strtree import STRtree
    trips = pd.read_csv(GTFS / "trips.txt", dtype=str)
    cal = pd.read_csv(GTFS / "calendar_dates.txt", dtype=str)
    shapes = pd.read_csv(GTFS / "shapes.txt", dtype={"shape_id": str})
    services = set(cal[(cal.date == date) & (cal.exception_type == "1")].service_id)
    day = trips[trips.service_id.isin(services)]
    n_shape = day.groupby("shape_id").size()
    tr = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)

    sids, geoms = [], []
    for sid, g in shapes[shapes.shape_id.isin(n_shape.index)].sort_values(
            ["shape_id", "shape_pt_sequence"]).groupby("shape_id"):
        x, y = tr.transform(g.shape_pt_lon.values, g.shape_pt_lat.values)
        sids.append(sid)
        geoms.append(LineString(zip(x, y)))
    tree = STRtree(geoms)
    street = street_namer()
    # each shape sampled every SAMPLE_M with its heading, to count its separate passes near a piece
    samples = []
    for g in geoms:
        d = np.arange(0, g.length, SAMPLE_M)
        xy = np.array([g.interpolate(t).coords[0] for t in d])
        hd = np.array([heading(g, t) for t in d])
        samples.append((xy, hd))
    noded = unary_union(geoms)
    pieces = []
    for piece in (noded.geoms if hasattr(noded, "geoms") else [noded]):
        mid = piece.interpolate(0.5, normalized=True)
        h = heading(piece, piece.length / 2)
        key = street(mid)
        total, lines = 0, set()
        other = None  # nearest pass on the other carriageway of the same street (offset from the piece)
        for i in tree.query(mid.buffer(DUAL_M)):
            xy, hd = samples[i]
            dist = np.hypot(xy[:, 0] - mid.x, xy[:, 1] - mid.y)
            near = np.flatnonzero(dist <= DUAL_M)
            if len(near) == 0:
                continue
            # one pass = a run of consecutive samples near the piece (a loop or an out-and-back shape passes twice)
            for run in np.split(near, np.flatnonzero(np.diff(near) > 1) + 1):
                k = run[np.argmin(dist[run])]
                if dist[k] > 1.0 + SAMPLE_M / 2:  # not the piece's own course: a parallel carriageway ...
                    if abs((hd[k] - h + 90) % 180 - 90) > PARALLEL_DEG:
                        continue
                    k_key = street(Point(xy[k]))
                    same = bool(key) and k_key == key            # ... of the same street (ref, else name)
                    unnamed_island = not key and not k_key and dist[k] <= MERGE_M
                    if not (same or unnamed_island):
                        continue
                    if other is None or dist[k] < np.hypot(*other):
                        other = (xy[k][0] - mid.x, xy[k][1] - mid.y)
                total += int(n_shape[sids[i]])
                lines.add(sids[i].split("_")[0])
        if total:
            pieces.append((piece, total, ",".join(sorted(lines)), other))
    return pieces, len(day)


def street_namer():
    """Point -> street of the nearest drivable edge within 8 m: its road number (ref) if it has one, else its
    name ('' if neither). DK2's carriageways are partly named Białka, partly unnamed, all ref 2."""
    import osmnx as ox
    from shapely.strtree import STRtree
    e = ox.graph_to_gdfs(ox.load_graphml(ROOT / "data" / "cache" / "osm_drive_service.graphml"), nodes=False)
    def key(ref, name):
        v = ref if isinstance(ref, (str, list)) else name
        return v if isinstance(v, str) else "/".join(v) if isinstance(v, list) else ""
    e = e[e["name"].notna() | e["ref"].notna()].to_crs(CRS)
    names = [key(r, n) for r, n in zip(e["ref"], e["name"])]
    lines = list(e.geometry)
    tree = STRtree(lines)

    def name_at(p):
        idx = tree.query(p.buffer(8))
        if len(idx) == 0:
            return ""
        return names[min(idx, key=lambda i: lines[i].distance(p))]
    return name_at


def merged_segments(pieces):
    """Join pieces with the same count and line set into longer LineStrings (for the gpkg)."""
    groups = defaultdict(list)
    for g, n, ls, _ in pieces:
        groups[(n, ls)].append(g)
    rows = []
    for (n, ls), segs in groups.items():
        m = linemerge(segs)
        for g in (m.geoms if hasattr(m, "geoms") else [m]):
            rows.append(dict(trips=n, lines=ls, n_lines=len(ls.split(",")), geometry=g))
    return gpd.GeoDataFrame(rows, crs=CRS)


def width_of(v, vmax):
    return 0.6 + (MAX_WIDTH_PT - 0.6) * np.minimum(v, vmax) / vmax


def color_of(v, vmax):
    return plt.get_cmap("plasma")(0.25 + 0.75 * np.sqrt(min(v, vmax) / vmax))


def scale_max(pieces):
    """Count at which SCALE_PCTL % of the network length is reached: short extremes (the Dworzec PKP spur,
    passed in and out by every bus going there) are drawn at full width without thinning everything else."""
    n = np.array([p[1] for p in pieces])
    L = np.array([p[0].length for p in pieces])
    o = np.argsort(n)
    cum = np.cumsum(L[o]) / L.sum()
    return int(n[o][np.searchsorted(cum, SCALE_PCTL / 100)])


def draw(pieces, n_trips, date, frame, out, limits=None, caps="round", cap=True):
    w, h = frame[2] - frame[0], frame[3] - frame[1]
    map_w = PAGE_W_IN - 2 * MARGIN_IN
    map_h = map_w * h / w
    page_h = map_h + 2 * MARGIN_IN
    fig = plt.figure(figsize=(PAGE_W_IN, page_h), facecolor=PAGE)
    ax = fig.add_axes([MARGIN_IN / PAGE_W_IN, MARGIN_IN / page_h, map_w / PAGE_W_IN, map_h / page_h])
    ax.set_facecolor("#0b0b0e")
    ax.set_xlim(frame[0], frame[2])
    ax.set_ylim(frame[1], frame[3])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#c8c8d0")
        s.set_linewidth(0.8)

    clip = gpd.GeoSeries([box(*frame).buffer(500)], crs=CRS).to_crs(4326).total_bounds.tolist()
    for name, style in [("water", dict(facecolor="#1a2e4a", edgecolor="none")),
                        ("forests", dict(facecolor="#0e130f", edgecolor="none")),
                        ("buildings", dict(facecolor="#1c1c22", edgecolor="none")),
                        ("roads", dict(color="#000000", linewidth=0.9)),
                        ("railways", dict(color="#2c2c34", linewidth=0.8))]:
        f = BASEMAP / f"{name}.gpkg"
        if f.exists():
            g = gpd.read_file(f, bbox=clip)
            g = g[g.geometry.notna()].to_crs(CRS)
            if len(g):
                g.plot(ax=ax, zorder=1, aspect=None, **style)
    if limits is not None:
        limits.boundary.plot(ax=ax, color="#8a8a96", linewidth=1.0, linestyle=(0, (4, 3)), zorder=2, aspect=None)

    pieces = sorted(pieces, key=lambda p: p[1])  # busiest drawn last, on top
    n = np.array([p[1] for p in pieces])
    vmax = scale_max(pieces) if cap else int(max(p[1] for p in pieces))
    top = int(n.max())
    if caps == "flat":
        # runs of equal count merged into continuous lines, flat ends: width steps meet squarely and branches
        # leave the trunk without the round bulbs every short piece gets with round caps
        runs = defaultdict(list)
        for g, v, _, _ in pieces:
            runs[v].append(g)
        segs, vals = [], []
        for v in sorted(runs):
            m = linemerge(runs[v])
            for g in (m.geoms if hasattr(m, "geoms") else [m]):
                segs.append(list(g.coords))
                vals.append(v)
        vals = np.array(vals)
        ax.add_collection(LineCollection(segs, linewidths=width_of(vals, vmax),
                                         colors=[color_of(v, vmax) for v in vals], capstyle="butt",
                                         joinstyle="round", zorder=3))
    else:
        # every piece where it is (stays connected), flat ends: round caps on hundreds of short pieces bulge
        # wherever the width changes. Bends are closed by a dot at each joint, as wide as the thinnest piece
        # meeting there (fills the notch, never sticks out of a thicker line); dead ends stay flat.
        # A carriageway of a divided road gets a filler copy shifted halfway towards the other one, so both
        # carriageways and the gap between them read as one street.
        segs = [list(p[0].coords) for p in pieces]
        segs += [list(translate(p[0], p[3][0] / 2, p[3][1] / 2).coords) for p in pieces if p[3] is not None]
        vals = np.concatenate([n, [p[1] for p in pieces if p[3] is not None]])
        joint = defaultdict(list)
        for c, v in zip(segs, vals):
            for end in (c[0], c[-1]):
                joint[(round(end[0], 1), round(end[1], 1))].append(v)
        jx, jv = [], []
        for xy, vs in joint.items():
            if len(vs) > 1:
                jx.append(xy)
                jv.append(min(vs))
        jx, jv = np.array(jx), np.array(jv)
        ax.scatter(jx[:, 0], jx[:, 1], s=width_of(jv, vmax) ** 2, c=[color_of(v, vmax) for v in jv],
                   linewidths=0, zorder=3)
        order = np.argsort(vals, kind="stable")
        ax.add_collection(LineCollection([segs[i] for i in order], linewidths=width_of(vals[order], vmax),
                                         colors=[color_of(v, vmax) for v in vals[order]], capstyle="butt",
                                         joinstyle="round", zorder=3.1))

    # header (as TransitFrequency's draw_headline: city name two rows tall, question + date beside it)
    halo = [patheffects.withStroke(linewidth=6, foreground=PAGE, alpha=0.55)]
    title_pt, sub_pt = 22.0, 13.6
    city = ax.text(0.022, 0.975, "BIAŁA PODLASKA", transform=ax.transAxes, ha="left", va="top",
                   fontproperties=bahn(58), color=INK, path_effects=halo, zorder=10)
    fig.canvas.draw()
    bb = city.get_window_extent().transformed(ax.transAxes.inverted())
    x2 = bb.x1 + 0.015
    ax.text(x2, 0.972, "Ile autobusów MZK przejeżdża ulicą w ciągu doby?", transform=ax.transAxes,
            ha="left", va="top", fontproperties=bahn(title_pt), color=INK, path_effects=halo, zorder=10)
    d = pd.Timestamp(date)
    ax.text(x2, bb.y0 + 0.004, f"dzień powszedni · {DAY_NAMES[d.weekday()]} {d:%d.%m.%Y} · {n_trips} kursów, "
            "oba kierunki", transform=ax.transAxes, ha="left", va="bottom", fontproperties=bahn(sub_pt),
            color=INK_DIM, path_effects=halo, zorder=10)

    # legend box, bottom right
    ticks = [v for v in (10, 50, 100, 200, 300) if v < vmax * 0.9] + [int(vmax)]
    lw_in, row_in = 2.3, 0.33
    lh_in = 0.55 + row_in * len(ticks)
    lx = 1 - (lw_in + 0.15) / map_w
    ly = 0.15 / map_h
    lax = fig.add_axes([ax.get_position().x0 + lx * ax.get_position().width,
                        ax.get_position().y0 + ly * ax.get_position().height,
                        lw_in / PAGE_W_IN, lh_in / page_h])
    lax.set_xlim(0, lw_in)
    lax.set_ylim(0, lh_in)
    lax.set_facecolor(PAGE)
    lax.set_xticks([])
    lax.set_yticks([])
    for s in lax.spines.values():
        s.set_color("#c8c8d0")
        s.set_linewidth(0.8)
    lax.text(0.15, lh_in - 0.12, "Autobusy/doba", ha="left", va="top", fontproperties=bahn(13), color=INK)
    for i, v in enumerate(ticks):
        y = lh_in - 0.62 - i * row_in
        lax.plot([0.18, 1.05], [y, y], lw=width_of(v, vmax), color=color_of(v, vmax), solid_capstyle="butt")
        label = f"{v}" if v < vmax else (f"{v}+" if top > vmax else f"{v}  (maks.)")
        lax.text(1.3, y, label, ha="left", va="center",
                 fontproperties=bahn(12), color=INK)

    # scale bar + north arrow, bottom left
    km = 3 if w > 9000 else 2
    sx, sy = frame[0] + w * 0.025, frame[1] + h * 0.045
    for k in range(km * 2):
        ax.add_patch(plt.Rectangle((sx + k * 500, sy), 500, h * 0.006, facecolor=INK if k % 2 == 0 else PAGE,
                                   edgecolor=INK, linewidth=0.6, zorder=10))
    ax.text(sx, sy + h * 0.014, "0", ha="center", va="bottom", fontproperties=bahn(10), color=INK, zorder=10)
    ax.text(sx + km * 1000, sy + h * 0.014, f"{km} km", ha="center", va="bottom", fontproperties=bahn(10),
            color=INK, zorder=10)
    ax.annotate("", xy=(sx + 120, sy + h * 0.16), xytext=(sx + 120, sy + h * 0.08),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2, mutation_scale=18), zorder=10)
    ax.text(sx + 120, sy + h * 0.17, "N", ha="center", va="bottom", fontproperties=bahn(11), color=INK, zorder=10)

    ax.text(0.5, 0.012, "Dane: rozkład MZK Biała Podlaska (GTFS nieoficjalny) · © OpenStreetMap contributors, "
            "2026 Jacek Gęborys", transform=ax.transAxes, ha="center", va="bottom", fontproperties=bahn(12),
            color=INK, path_effects=halo, zorder=10)
    fig.savefig(out, dpi=150, facecolor=PAGE)
    plt.close(fig)
    print(f"written {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="20261118")
    ap.add_argument("--city", action="store_true", help="frame on the city limits (+400 m) instead of all stops")
    ap.add_argument("--caps", choices=["round", "flat"], default="round", help="line ends: round pieces or flat runs")
    ap.add_argument("--no-cap", action="store_true", help="scale widths to the real maximum (the PKP spur) instead of the 99.5th percentile")
    ap.add_argument("--out", help="output PNG (default png/cartogram_<date>[_city].png)")
    args = ap.parse_args()

    pieces, n_trips = piece_counts(args.date)
    seg = merged_segments(pieces)
    seg.to_file(GPKG, layer="segments", driver="GPKG")
    print(f"segments: {len(seg)} lines, busiest {seg.trips.max()} buses/day "
          f"({seg.sort_values('trips').iloc[-1].lines}) -> {GPKG.name}:segments")
    out = Path(args.out) if args.out else ROOT / "png" / f"cartogram_{args.date}{'_city' if args.city else ''}.png"
    if args.city:
        limits = gpd.read_file(LIMITS).to_crs(CRS)
        x0, y0, x1, y1 = limits.total_bounds
        draw(pieces, n_trips, args.date, (x0 - 400, y0 - 400, x1 + 400, y1 + 400), out, limits, args.caps, not args.no_cap)
    else:
        draw(pieces, n_trips, args.date, FRAME_ALL, out, caps=args.caps, cap=not args.no_cap)


if __name__ == "__main__":
    main()
