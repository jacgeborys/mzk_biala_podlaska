"""Flow cartogram on a simplified street network: every edge drawn as thick as the buses passing it in a day.

1. corridors: all shapes of the day folded into one centreline per street (directions, carriageways, variants
   within TOL m share it; cut ends re-attached where the route goes on)
2. counts: each shape sampled every STEP_M m, each sample snapped to the nearest parallel corridor piece; a pass =
   a run of samples on one piece in one direction (an out-and-back spur counts twice); trips added per 10 m bin
   along the stretch the pass covers
3. graph: the centreline noded at crossings and wherever the count steps (a branch folded into the corridor, a
   route end); edges = runs between nodes with one value each (length-weighted median), so the count only
   changes at nodes; short edges (roundabouts, staggered crossings) contracted into one node

Output: png/cartogram_<date>[_city].png, layer `segments` in data/gis/biala_podlaska_gtfs.gpkg
        (one LineString per edge: trips, lines)

Usage:  python cartogram.py [--date 20261118] [--city]   (--city: frame on the city limits, drawn dashed)
"""
import argparse
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import matplotlib
import networkx as nx
import numpy as np
import pandas as pd
from pyproj import Transformer
from shapely.geometry import LineString, Point, box
from shapely.ops import nearest_points, substring, unary_union
from shapely.strtree import STRtree

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
TOL = 35              # m: shapes this close share one corridor (both directions, dual carriageways, variants)
MIN_PART = 40         # m: shorter leftovers of the corridor folding are dropped
DUP_M = 10            # m: corridor pieces overlapping another one within this distance are cut back
PARALLEL_DEG = 35     # a sample snaps only to a corridor piece within this angle of its heading
STEP_M = 5            # shapes sampled this often
CONTRACT_M = 40       # nodes joined by shorter edges become one node (roundabouts, staggered crossings)
SHORT_M = 60          # count steps closer than this are merged (corridor links and junctions: noisy counts)
BIN_M = 10            # passes counted in bins of this length along each corridor piece
SMOOTH_BINS = 5       # median filter over the bins (50 m)
SIMPLIFY_M = 6        # drawn edge geometry simplified to this tolerance
RB_MARGIN_M = 10      # roundabouts: edges cut this far outside the ring and joined straight to one node
RB_MAX_R = 60         # rings with a larger radius are left as they are
LOOP_M = 250          # edges starting and ending at the same node shorter than this are turning loops (dropped)
SPLIT_SUM = 1.15      # a node is a split / merge when the other forward edges' trips add up to at most this x the trunk
TAPER_M = 180         # a stacked branch eases from its slot in the trunk to its own centreline over this length
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


def parts(g):
    if g.is_empty:
        return []
    return [p for p in (g.geoms if hasattr(g, "geoms") else [g]) if isinstance(p, LineString)]


def heading(g, d):
    """Heading of line g at distance d (radians, over +-5 m)."""
    a = g.interpolate(max(0.0, d - 5))
    b = g.interpolate(min(g.length, d + 5))
    return np.arctan2(b.y - a.y, b.x - a.x)


def ang_diff(a, b):
    """Smallest difference of two undirected headings, degrees."""
    return abs((np.degrees(a - b) + 90) % 180 - 90)


def key(c):
    return round(c[0], 1), round(c[1], 1)


def load(date):
    trips = pd.read_csv(GTFS / "trips.txt", dtype=str)
    cal = pd.read_csv(GTFS / "calendar_dates.txt", dtype=str)
    shapes = pd.read_csv(GTFS / "shapes.txt", dtype={"shape_id": str})
    day = trips[trips.service_id.isin(set(cal[(cal.date == date) & (cal.exception_type == "1")].service_id))]
    n_shape = day.groupby("shape_id").size()
    line_of = day.drop_duplicates("shape_id").set_index("shape_id").route_id
    tr = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    geoms = {}
    for sid, g in shapes[shapes.shape_id.isin(n_shape.index)].sort_values(
            ["shape_id", "shape_pt_sequence"]).groupby("shape_id"):
        x, y = tr.transform(g.shape_pt_lon.values, g.shape_pt_lat.values)
        geoms[sid] = LineString(zip(x, y)).simplify(2)
    return geoms, n_shape, line_of, len(day)


# ---------------------------------------------------------------- 1-2: corridors and graph
def self_fold(p, step=5.0):
    """Line p without the stretches where it runs back along itself (within DUP_M of its own course at least
    2 * DUP_M earlier, for more than 3 * DUP_M; no street comes that close to itself, a hairpin does at once): an
    out-and-back or loop shape keeps one copy of the street, up to its turning point."""
    ds = np.arange(0, p.length, step)
    dup = np.array([s > 2 * DUP_M and substring(p, 0, s - 2 * DUP_M).distance(p.interpolate(s)) <= DUP_M for s in ds])
    if not dup.any():
        return [p]
    keep, a = [], 0.0
    for run in np.split(np.flatnonzero(dup), np.flatnonzero(np.diff(np.flatnonzero(dup)) > 1) + 1):
        if len(run) * step > 3 * DUP_M:
            keep.append((a, ds[run[0]]))
            a = min(ds[run[-1]] + step, p.length)
    keep.append((a, p.length))
    return [substring(p, x, y) for x, y in keep if y - x > 0.5]


def corridors(geoms, n_shape):
    """One centreline per street: shapes added busiest first, each keeping only what lies farther than TOL from
    the corridors so far; every cut end linked to the corridor its route continues on (40 m past the cut)."""
    kept, src = [], []
    for s in sorted(geoms, key=lambda s: -n_shape[s]):
        g = geoms[s]
        rest = g.difference(unary_union(kept).buffer(TOL)) if kept else g
        for p in [q for r in parts(rest) for q in self_fold(r)]:
            if p.length >= MIN_PART:
                kept.append(p)
                src.append(g)
    links = []
    tree = STRtree(kept)
    for i, (p, g) in enumerate(zip(kept, src)):
        for at_start, end in ((True, Point(p.coords[0])), (False, Point(p.coords[-1]))):
            d = g.project(end)
            if d < 1 or d > g.length - 1:            # a real route end, not a cut
                continue
            ahead = g.interpolate(d - 40 if at_start else d + 40)
            cand = [j for j in tree.query(ahead.buffer(TOL + 30)) if j != i]
            if not cand:
                continue
            j = min(cand, key=lambda j: kept[j].distance(ahead))
            if kept[j].distance(ahead) > TOL + 30:
                continue
            q = kept[j].interpolate(kept[j].project(ahead))
            if end.distance(q) > 0.5:
                links.append(LineString([end, q]))
    # a shape is not folded against itself: a loop or out-and-back running back on its own street leaves a second
    # copy of that street (the Dworzec PKP spur); such overlaps are cut from the shorter piece and its ends relinked
    acc, relinks = [], []
    for p in sorted(parts(unary_union(kept + links)), key=lambda p: -p.length):
        near = unary_union(acc).buffer(DUP_M) if acc else None
        if near is None or p.intersection(near).length <= 3 * DUP_M:
            acc.append(p)
            continue
        base = unary_union(acc)
        for r in parts(p.difference(near)):
            if r.length < MIN_PART:
                continue
            acc.append(r)
            for end in (Point(r.coords[0]), Point(r.coords[-1])):
                if base.distance(end) <= DUP_M + 1:
                    relinks.append(LineString([end, nearest_points(base, end)[0]]))
    return parts(unary_union(acc + relinks))


# ---------------------------------------------------------------- 3: passes, counted along each piece
def piece_profiles(pieces, geoms, n_shape, line_of):
    """Trips per BIN_M bin along every piece. Samples snap to the nearest parallel piece within TOL + 8 m; a pass
    is a run of consecutive samples on one piece in one direction and adds its trips only to the bins it runs
    along, so a piece carrying a folded-in branch (no junction there) still shows where the count steps."""
    tree = STRtree(pieces)
    names = sorted(set(line_of))
    prof = [np.zeros(int(np.ceil(p.length / BIN_M))) for p in pieces]
    on = [np.zeros((len(t), len(names)), bool) for t in prof]
    for sid, g in geoms.items():
        seq, pos = [], []  # (piece, direction) and position along the piece per sample
        for d in np.arange(0, g.length, STEP_M):
            pt = g.interpolate(d)
            h = heading(g, d)
            best = None
            for i in tree.query(pt.buffer(TOL + 8)):
                dist = pieces[i].distance(pt)
                if dist > TOL + 8:
                    continue
                t = pieces[i].project(pt)
                hp = heading(pieces[i], t)
                if ang_diff(hp, h) <= PARALLEL_DEG and (best is None or dist < best[0]):
                    best = (dist, i, bool(np.cos(hp - h) > 0), t)
            seq.append(best[1:3] if best else None)
            pos.append(best[3] if best else 0.0)
        k = 0
        while k < len(seq):
            j = k
            while j + 1 < len(seq) and seq[j + 1] == seq[k]:
                j += 1
            if seq[k] is not None:
                i = seq[k][0]
                lo, hi = min(pos[k:j + 1]) - STEP_M / 2, max(pos[k:j + 1]) + STEP_M / 2
                b0, b1 = max(0, int(lo // BIN_M)), min(len(prof[i]), int(hi // BIN_M) + 1)
                prof[i][b0:b1] += n_shape[sid]
                on[i][b0:b1, names.index(line_of[sid])] = True
            k = j + 1
    return prof, on, names


def split_by_profile(pieces, prof, on, names):
    """Pieces cut where their count steps: profile median-filtered over SMOOTH_BINS bins (junction noise), runs
    shorter than SHORT_M merged into the neighbour with the closer value. Each part: trips = median of its bins."""
    items = []
    for p, t, o in zip(pieces, prof, on):
        h = SMOOTH_BINS // 2
        pad = np.concatenate([np.repeat(t[:1], h), t, np.repeat(t[-1:], h)])
        sm = np.array([np.median(pad[k:k + SMOOTH_BINS]) for k in range(len(t))])
        cuts = [0] + [k for k in range(1, len(t)) if sm[k] != sm[k - 1]] + [len(t)]
        runs = [[a, b, float(np.median(t[a:b]))] for a, b in zip(cuts, cuts[1:])]
        while len(runs) > 1:
            k = min(range(len(runs)), key=lambda k: runs[k][1] - runs[k][0])
            if (runs[k][1] - runs[k][0]) * BIN_M >= SHORT_M:
                break
            nb = [m for m in (k - 1, k + 1) if 0 <= m < len(runs)]
            m = min(nb, key=lambda m: abs(runs[m][2] - runs[k][2]))
            a, b = min(runs[k][0], runs[m][0]), max(runs[k][1], runs[m][1])
            runs[min(k, m)] = [a, b, float(np.median(t[a:b]))]
            del runs[max(k, m)]
        for a, b, v in runs:
            g = substring(p, a * BIN_M, min(b * BIN_M, p.length)) if len(runs) > 1 else p
            items.append(dict(geom=g, pt=v, u=key(g.coords[0]), v=key(g.coords[-1]),
                              pl={l: o[a:b, n].sum() * BIN_M for n, l in enumerate(names) if o[a:b, n].any()}))
    return items


# ---------------------------------------------------------------- 3b: edges with one value each
def chain(items):
    """Join items {geom, u, v, ...} through nodes of degree 2 when their values are equal or one of them is
    shorter than SHORT_M (a different value on two longer pieces is a real step: a branch folded into the corridor
    leaves no junction, e.g. the Dworzec PKP spur). The joined item keeps its source pieces in `parts`."""
    items = [dict(it, parts=it.get("parts", [it])) for it in items]
    for it in items:
        it["trips"], it["lines"] = edge_value(it)
    changed = True
    while changed:
        changed = False
        inc = defaultdict(list)
        for i, e in enumerate(items):
            inc[e["u"]].append(i)
            inc[e["v"]].append(i)
        for n, es in inc.items():
            if len(es) != 2 or es[0] == es[1]:
                continue
            a, b = items[es[0]], items[es[1]]
            if a["trips"] != b["trips"] and min(a["geom"].length, b["geom"].length) >= SHORT_M:
                continue
            ca = list(a["geom"].coords) if a["v"] == n else list(a["geom"].coords)[::-1]
            cb = list(b["geom"].coords) if b["u"] == n else list(b["geom"].coords)[::-1]
            g = LineString(ca + cb[1:])
            merged = dict(geom=g, u=key(g.coords[0]), v=key(g.coords[-1]), parts=a["parts"] + b["parts"])
            merged["trips"], merged["lines"] = edge_value(merged)
            items = [it for k, it in enumerate(items) if k not in es] + [merged]
            changed = True
            break
    return items


def edge_value(e):
    """Length-weighted median of the source pieces' trips; lines covering at least half of the edge."""
    ps = e["parts"]
    t = np.array([p["pt"] for p in ps])
    L = np.array([p["geom"].length for p in ps])
    o = np.argsort(t)
    trips = int(t[o][np.searchsorted(np.cumsum(L[o]), L.sum() / 2)])
    cover = defaultdict(float)
    for p in ps:
        for l, v in p["pl"].items():
            cover[l] += min(v, p["geom"].length)
    return trips, frozenset(l for l, v in cover.items() if v >= 0.5 * L.sum())


def contract_short(edges):
    """Nodes joined by edges shorter than CONTRACT_M, and dead ends closer than that to each other, become one node
    at their centre (short dead-end edges are kept)."""
    deg = defaultdict(int)
    for e in edges:
        deg[e["u"]] += 1
        deg[e["v"]] += 1
    G = nx.Graph()
    G.add_edges_from((e["u"], e["v"]) for e in edges
                     if e["geom"].length < CONTRACT_M and e["u"] != e["v"] and deg[e["u"]] > 1 and deg[e["v"]] > 1)
    ends = [n for n, k in deg.items() if k == 1]               # two loose ends this close: a corridor cut not relinked
    for x in range(len(ends)):
        for y in range(x + 1, len(ends)):
            if np.hypot(ends[x][0] - ends[y][0], ends[x][1] - ends[y][1]) < CONTRACT_M:
                G.add_edge(ends[x], ends[y])
    rep = {}
    for comp in nx.connected_components(G):
        c = tuple(np.round(np.array(list(comp)).mean(axis=0), 1))
        for n in comp:
            rep[n] = c
    out = []
    for e in edges:
        u, v = rep.get(e["u"], e["u"]), rep.get(e["v"], e["v"])
        if u == v and e["geom"].length < CONTRACT_M:
            continue
        coords = list(e["geom"].coords)
        coords[0], coords[-1] = u, v
        out.append(dict(e, geom=LineString(coords), u=u, v=v))
    return out


def roundabouts():
    """(centre, radius) of every OSM roundabout (junction=roundabout rings), EPSG:2180."""
    import osmnx as ox
    e = ox.graph_to_gdfs(ox.load_graphml(ROOT / "data" / "cache" / "osm_drive_service.graphml"), nodes=False)
    e = e[e["junction"] == "roundabout"].to_crs(CRS)
    G = nx.Graph()
    for (u, v, _), g in zip(e.index, e.geometry):
        G.add_edge(u, v, g=g)
    out = []
    for comp in nx.connected_components(G):
        ring = unary_union([d["g"] for _, _, d in G.subgraph(comp).edges(data=True)])
        c = ring.convex_hull.centroid
        r = max(Point(xy).distance(c) for g in parts(ring) for xy in g.coords)
        if r <= RB_MAX_R:
            out.append((c, r))
    return out


def straighten_roundabouts(edges):
    """At every roundabout the network touches, edges are cut RB_MARGIN_M outside the ring and joined by straight
    spokes to one node (trip-weighted mean of the cut points): no arcs around the ring. Geometry only (counts
    unchanged); an edge passing through is split there and joined again by chain()."""
    n = 0
    for c, r in roundabouts():
        R = r + RB_MARGIN_M
        if not any(e["geom"].distance(c) < R for e in edges):
            continue
        n += 1
        disk = c.buffer(R)
        out, cut = [], []
        for e in edges:
            if e["geom"].distance(c) >= R:
                out.append(e)
                continue
            for g in parts(e["geom"].difference(disk)):
                xy = list(g.coords)
                ends = [Point(xy[0]).distance(c) < R + 0.5, Point(xy[-1]).distance(c) < R + 0.5]
                cut += [(xy[0] if k == 0 else xy[-1], e["trips"]) for k in (0, 1) if ends[k]]
                out.append((e, xy, ends))
        if not cut:
            continue
        # the node: trip-weighted mean of the cut points, so the busiest line through stays (nearly) straight
        w = np.array([t for _, t in cut], float)
        m = tuple(np.round((np.array([xy for xy, _ in cut]) * w[:, None]).sum(axis=0) / w.sum(), 1))
        res = []
        for o in out:
            if isinstance(o, dict):
                res.append(o)
                continue
            e, xy, ends = o
            xy = ([m] if ends[0] else []) + xy + ([m] if ends[1] else [])
            res.append(dict(e, geom=LineString(xy), u=key(xy[0]), v=key(xy[-1])))
        out = res
        edges = out
    print(f"roundabouts straightened: {n}")
    return edges


def network(date):
    geoms, n_shape, line_of, n_trips = load(date)
    pieces = corridors(geoms, n_shape)
    items = split_by_profile(pieces, *piece_profiles(pieces, geoms, n_shape, line_of))
    edges = [e for e in chain(items) if e["trips"] > 0]   # junction to junction, route end or real step
    edges = [e for e in contract_short(edges)          # roundabouts to one node; turning loops (their buses are
             if not (e["u"] == e["v"] and e["geom"].length < LOOP_M)]   # counted on the way in and out) dropped
    edges = chain(straighten_roundabouts(edges))   # roundabouts as one node, neighbours joined again
    print(f"network: {len(pieces)} corridor pieces, {len(items)} after cuts at steps -> {len(edges)} edges, "
          f"{len({n for e in edges for n in (e['u'], e['v'])})} nodes")
    return edges, n_trips


def network_gdfs(edges):
    """The simplified network for reuse: edges (trips, lines, from_node, to_node) and nodes (node_id, degree,
    trips = busiest edge there)."""
    nodes = sorted({n for e in edges for n in (e["u"], e["v"])})
    nid = {n: k for k, n in enumerate(nodes)}
    deg, top = defaultdict(int), defaultdict(int)
    for e in edges:
        for n in (e["u"], e["v"]):
            deg[n] += 1
            top[n] = max(top[n], e["trips"])
    seg = gpd.GeoDataFrame([dict(trips=e["trips"], lines=",".join(sorted(e["lines"])), n_lines=len(e["lines"]),
                                 from_node=nid[e["u"]], to_node=nid[e["v"]], length_m=round(e["geom"].length, 1),
                                 geometry=e["geom"]) for e in edges], crs=CRS)
    nod = gpd.GeoDataFrame([dict(node_id=nid[n], degree=deg[n], trips=top[n], geometry=Point(n)) for n in nodes],
                           crs=CRS)
    return seg, nod


def branch_offsets(edges, m_per_pt, vmax):
    """Per edge the lateral offset (m, left of its own direction) at its start and end. At a node where a trunk
    (busiest edge) divides into edges whose trips add up to at most SPLIT_SUM of it, those edges are stacked side by
    side across the trunk in the order they turn (leftmost on the left, slot by share of trips, kept inside the
    trunk), so branches peel off its sides; at a plain
    crossing or a simple width step everything stays centred. Also returns the stacked nodes and edge ends."""
    W = {i: width_of(e["trips"], vmax) * m_per_pt for i, e in enumerate(edges)}
    inc = defaultdict(list)
    for i, e in enumerate(edges):
        c = list(e["geom"].coords)
        for at_u, n in ((True, e["u"]), (False, e["v"])):
            g = e["geom"]
            q = g.interpolate(min(30.0, g.length / 2) if at_u else max(0.0, g.length - min(30.0, g.length / 2)))
            o = c[0] if at_u else c[-1]
            inc[n].append((i, at_u, np.arctan2(q.y - o[1], q.x - o[0])))
    off = defaultdict(lambda: [0.0, 0.0])
    stacked, slots = set(), set()
    for n, es in inc.items():
        if len(es) < 3:
            continue
        t, _, ht = max(es, key=lambda x: W[x[0]])
        h_in = ht + np.pi
        cand = []
        for i, at_u, h in es:
            if i == t:
                continue
            rel = (h - h_in + np.pi) % (2 * np.pi) - np.pi
            if abs(rel) < np.radians(110):
                cand.append((rel, i, at_u))
        T = sum(edges[i]["trips"] for _, i, _ in cand)
        if len(cand) < 2 or T > SPLIT_SUM * edges[t]["trips"]:
            continue
        stacked.add(n)
        cum = 0.0
        for rel, i, at_u in sorted(cand, reverse=True):          # leftmost turn first
            share = edges[i]["trips"] / T
            room = max(0.0, (W[t] - W[i]) / 2)                     # the branch stays inside the trunk
            lat = float(np.clip(W[t] * (0.5 - cum - share / 2), -room, room)) * np.cos(rel)
            cum += share
            off[i][0 if at_u else 1] = lat if at_u else -lat       # left of the direction into the node
            slots.add((i, 0 if at_u else 1))
    return off, stacked, slots


def offset_line(g, o_start, o_end, step=4.0):
    """g shifted sideways by o_start at its start and o_end at its end, each easing to 0 over TAPER_M (cosine)."""
    if abs(o_start) < 0.01 and abs(o_end) < 0.01:
        return list(g.coords)
    L = g.length
    d = np.append(np.arange(0, L, step), L)
    xy = np.array([g.interpolate(t).coords[0] for t in d])
    tx, ty = np.gradient(xy[:, 0]), np.gradient(xy[:, 1])
    nrm = np.hypot(tx, ty)
    nrm[nrm == 0] = 1
    T = min(TAPER_M, L / 2)
    f = lambda s: np.where(s < T, 0.5 * (1 + np.cos(np.pi * np.clip(s / T, 0, 1))), 0.0)  # noqa: E731
    o = o_start * f(d) + o_end * f(L - d)
    return list(zip(xy[:, 0] - ty / nrm * o, xy[:, 1] + tx / nrm * o))


def width_of(v, vmax):
    return 0.6 + (MAX_WIDTH_PT - 0.6) * np.minimum(v, vmax) / vmax


def color_of(v, vmax):
    return plt.get_cmap("plasma")(0.25 + 0.75 * np.sqrt(min(v, vmax) / vmax))


def draw(edges, n_trips, date, frame, out, limits=None):
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

    # one line per edge, flat ends; where a trunk divides, its branches start side by side inside it (branch_offsets).
    # Each node gets a dot as wide as the second-widest edge there, which closes the bend of a line passing through
    # without sticking out where a thick edge simply continues as a thin one
    edges = sorted(edges, key=lambda e: e["trips"])  # busiest drawn last, on top
    vals = np.array([e["trips"] for e in edges])
    vmax = top = int(vals.max())
    off, stacked, slots = branch_offsets(edges, w / (map_w * 72), vmax)
    at = defaultdict(list)
    for e in edges:
        at[e["u"]].append(e["trips"])
        at[e["v"]].append(e["trips"])
    jx = np.array([n for n, vs in at.items() if len(vs) > 1 and n not in stacked])
    jv = np.array([sorted(vs)[-2] for n, vs in at.items() if len(vs) > 1 and n not in stacked])
    ax.scatter(jx[:, 0], jx[:, 1], s=width_of(jv, vmax) ** 2, c=[color_of(v, vmax) for v in jv],
               linewidths=0, zorder=3)
    lines = [offset_line(e["geom"].simplify(SIMPLIFY_M), *off[i]) for i, e in enumerate(edges)]
    # a stacked branch gets a round cap at its slot in the trunk: closes the corner where it bends away
    caps = [(lines[i][0 if end == 0 else -1], edges[i]["trips"]) for i, end in slots]
    if caps:
        ax.scatter([c[0][0] for c in caps], [c[0][1] for c in caps], s=width_of(np.array([c[1] for c in caps]), vmax) ** 2,
                   c=[color_of(c[1], vmax) for c in caps], linewidths=0, zorder=3)
    ax.add_collection(LineCollection(lines,
                                     linewidths=width_of(vals, vmax), colors=[color_of(v, vmax) for v in vals],
                                     capstyle="butt", joinstyle="round", zorder=3.1))

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

    # scale bar, bottom left (no north arrow: north is up, as on the other maps)
    km = 3 if w > 9000 else 2
    sx, sy = frame[0] + w * 0.025, frame[1] + h * 0.045
    for k in range(km * 2):
        ax.add_patch(plt.Rectangle((sx + k * 500, sy), 500, h * 0.006, facecolor=INK if k % 2 == 0 else PAGE,
                                   edgecolor=INK, linewidth=0.6, zorder=10))
    ax.text(sx, sy + h * 0.014, "0", ha="center", va="bottom", fontproperties=bahn(10), color=INK, zorder=10)
    ax.text(sx + km * 1000, sy + h * 0.014, f"{km} km", ha="center", va="bottom", fontproperties=bahn(10),
            color=INK, zorder=10)

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
    ap.add_argument("--out", help="output PNG (default png/cartogram_<date>[_city].png)")
    args = ap.parse_args()

    edges, n_trips = network(args.date)
    seg, nod = network_gdfs(edges)
    seg.to_file(GPKG, layer="segments", driver="GPKG")
    nod.to_file(GPKG, layer="segment_nodes", driver="GPKG")
    print(f"segments: {len(seg)} edges, busiest {seg.trips.max()} buses/day "
          f"({seg.sort_values('trips').iloc[-1].lines}) -> {GPKG.name}:segments")
    out = Path(args.out) if args.out else ROOT / "png" / f"cartogram_{args.date}{'_city' if args.city else ''}.png"
    if args.city:
        limits = gpd.read_file(LIMITS).to_crs(CRS)
        x0, y0, x1, y1 = limits.total_bounds
        draw(edges, n_trips, args.date, (x0 - 400, y0 - 400, x1 + 400, y1 + 400), out, limits)
    else:
        draw(edges, n_trips, args.date, FRAME_ALL, out)


if __name__ == "__main__":
    main()
