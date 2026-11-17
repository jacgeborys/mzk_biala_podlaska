"""Bus network diagram for Biała Podlaska: lines bundled side by side on the real streets, one colour per line.

Built for a small town (10 km across) straight from the feed in data/gtfs/, no tram-network machinery:
  1. corridors: every shape of the day, busiest lines first, minus whatever already lies within TOL m of the
     corridors kept so far (opposite directions, divided carriageways and parallel variants fold into one);
     cut-off branch ends are joined back to the corridor they left
  2. the corridors are noded into pieces; a line sits on a piece when its own shapes run along it (TOL m, parallel)
  3. pieces with the same lines are chained into edges between junctions
  4. one global line order for all bundles, searched (annealing) to minimise lines crossing at junctions
  5. each line is drawn as continuous chains along its lane (offset curves, joined at junctions, corners
     rounded over a few metres); stops are black dots on the lines that call there; interchanges, termini and
     Dworzec PKP get a white capsule across the bundle and a name
Colours / legend names: config/diagram_lines.json. Output: png/bus_network_city.png + .pdf (vector, for print),
zoomed to the city limits.

Usage:  python diagram.py [--date 20261118] [--out png/bus_network_city.png] [--seed 1]
"""
import argparse
import json
import math
import random
import re
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import shapely
import pandas as pd
from pyproj import Transformer
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import linemerge, unary_union
from shapely.strtree import STRtree

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager, patheffects  # noqa: E402

ROOT = Path(__file__).resolve().parent
GTFS = ROOT / "data" / "gtfs"
BASEMAP = ROOT / "data" / "basemap"
LIMITS = ROOT / "data" / "gis" / "city_limits.gpkg"
CONFIG = ROOT / "config" / "diagram_lines.json"
CRS = "EPSG:2180"

TOL = 35            # m: routes this close and parallel are one corridor (directions, carriageways, variants)
MIN_PART = 40       # m: shorter leftovers of the corridor dedupe are dropped
PARALLEL_DEG = 35
SPACING = 26        # m between lanes of a bundle
STEP_M = 6          # drawn lines are sampled this often
SIMPLIFY_M = 30     # sections are simplified this much before their corners are rounded (street wiggles go)
CORNER_R = 150      # m: corners become circular arcs of this radius (smaller where the legs are short)
RAMP_M = 120        # a line changing lane at a junction eases over this length (no jag)
EXT_TOL = 45        # m: shapes this close to a drawn point pass it (extension test)
EXT_MIN_M = 150     # m: shorter solid / dashed runs are merged into their neighbours
DASH_ON, DASH_OFF = 4.0, 0.55   # extensions: dash / gap length in line widths
EXT_SHARE = 0.25    # where a line runs with under this share of its trips of the day it is an extension (dashed)
U_DEPTH_M = 690      # Dworzec PKP: the U hangs this far below the Sidorska / Al. 1000-lecia crossing
U_HALF_M = 130       # ... its legs this far either side of Al. 1000-lecia (also the radius of its bottom)
STATION_CROSS_W = 0.3   # ... crossings / lane shifts at the station (under its bar) cost this much, so they go there
ABOVE_CROSS_W = 3.0     # ... and at the crossing above it this much (as few crossings above the station as possible)
U_TOP_M = 300        # ... the legs first run straight down side by side this far (arms turn into them on full arcs)
U_S_R = 160          # ... then move apart to their place on an S-bend of two arcs of this radius (room for every lane)
CONTRACT_M = 40     # junctions joined by shorter edges are merged into one
PATH_FOLD_EDGES, PATH_FOLD_M = 4, 80   # an edge beside a path of this many edges, this close, is the same link
TURN_STUB_M = 120   # dead-end stubs this short with their stop by the junction are terminal turning loops ...
TURN_DUMMY_M = 220  # ... drawn as a straight dummy street this long (which lines turn in stays readable)
TRIANGLE_M = 80     # ... and by shorter edges closing a triangle (a split rejoining right away)
SHIFT_COST = 0.15   # order search: cost per lane a line shifts through a junction (vs 1 per crossing)
KEY_STATIONS = {"Dworzec PKP"}   # always drawn as interchange with a name
MIN_INTERCHANGE = 5               # stations served by this many lines get a capsule and a name
MAIN_TERMINUS = 0.5               # a station is a line's terminus when this share of its trips start or end there
CONNECT_M = 500     # a line broken into pieces is joined along the corridors if the gap is at most this long
DROP_M = 150        # pieces of a line shorter than this that cannot be joined are dropped (strays)
END_JOIN_M = 400    # a loose end of a line inside the network continues this far at most to reach the line again
PAIR_M = 80         # two loose ends of one line this close are bridged
LEAD_FOLLOW_M = 120  # a line's continuation at a junction: its routes are followed this far past the node
SNAP_M = 10         # corridor ends this close are linked before lines are repaired (gaps closed in place)
RB_MARGIN_M = 10    # roundabouts: edges cut this far outside the ring and joined straight to one node
STUB_M = 150        # dead-end bits of a line this short that do not end where its routes end are dropped
FONT = Path("C:/Windows/Fonts/bahnschrift.ttf")
LAND, WATER, GREEN, INK = "#f5f2ea", "#a9cde6", "#d9e8cf", "#1d1d1f"
DAY_NAMES = ["poniedziałek", "wtorek", "środa", "czwartek", "piątek", "sobota", "niedziela"]


def font(size):
    if FONT.exists():
        return font_manager.FontProperties(fname=str(FONT), size=size)
    return font_manager.FontProperties(family="sans-serif", size=size)


def parts(g):
    if g.is_empty:
        return []
    return [p for p in (g.geoms if hasattr(g, "geoms") else [g]) if isinstance(p, LineString)]


def heading_at(g, d):
    a, b = g.interpolate(max(0.0, d - 6)), g.interpolate(min(g.length, d + 6))
    return math.atan2(b.y - a.y, b.x - a.x)


def ang_diff(a, b):
    """Smallest difference of two undirected headings, degrees."""
    return abs((math.degrees(a - b) + 90) % 180 - 90)


# ---------------------------------------------------------------- data
def load(date):
    trips = pd.read_csv(GTFS / "trips.txt", dtype=str)
    cal = pd.read_csv(GTFS / "calendar_dates.txt", dtype=str)
    st = pd.read_csv(GTFS / "stop_times.txt", dtype={"trip_id": str, "stop_id": str})
    stops = pd.read_csv(GTFS / "stops.txt", dtype={"stop_id": str, "parent_station": str})
    shapes = pd.read_csv(GTFS / "shapes.txt", dtype={"shape_id": str})
    day = trips[trips.service_id.isin(set(cal[(cal.date == date) & (cal.exception_type == "1")].service_id))]
    n_shape = day.groupby("shape_id").size()
    line_of = day.drop_duplicates("shape_id").set_index("shape_id").route_id
    tr = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    geoms = {}
    for sid, g in shapes[shapes.shape_id.isin(n_shape.index)].sort_values(["shape_id", "shape_pt_sequence"]).groupby("shape_id"):
        x, y = tr.transform(g.shape_pt_lon.values, g.shape_pt_lat.values)
        geoms[sid] = LineString(zip(x, y)).simplify(2)

    # stations (parent stops): position, lines calling, lines starting / ending there
    st = st[st.trip_id.isin(day.trip_id)].merge(day[["trip_id", "route_id"]], on="trip_id")
    parent = stops.set_index("stop_id").parent_station
    st["station"] = st.stop_id.map(parent)
    st["seq"] = st.stop_sequence.astype(int)
    first_last = pd.concat([st.loc[st.groupby("trip_id").seq.idxmin()], st.loc[st.groupby("trip_id").seq.idxmax()]])
    ends = first_last.groupby(["station", "route_id"]).size()
    line_trips = day.groupby("route_id").size()
    st_tab = stops[stops.location_type == 1].set_index("stop_id")
    sx, sy = tr.transform(st_tab.stop_lon.values, st_tab.stop_lat.values)
    stations = {}
    for sid_, name, x, y in zip(st_tab.index, st_tab.stop_name, sx, sy):
        lines = set(st[st.station == sid_].route_id)
        if not lines:
            continue
        term = {l for l in lines if ends.get((sid_, l), 0) >= MAIN_TERMINUS * line_trips[l]}
        stations[sid_] = dict(name=name, xy=np.array([x, y]), lines=lines, term=term)
    return geoms, n_shape, line_of, stations, line_trips


# ---------------------------------------------------------------- 1-3: corridors, pieces, edges
def corridors(geoms, n_shape, line_of, line_trips):
    order = sorted(geoms, key=lambda s: (-line_trips[line_of[s]], line_of[s], -n_shape[s]))
    kept, src = [], []
    for s in order:
        g = geoms[s]
        rest = g.difference(unary_union(kept).buffer(TOL)) if kept else g
        for p in parts(rest):
            if p.length >= MIN_PART:
                kept.append(p)
                src.append(g)
    # join every cut end back to the corridor its route continues on: follow the original route 40 m past the
    # cut (into the folded stretch) and link to the kept corridor nearest that point
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
    return parts(unary_union(kept + links))


def lines_on(pieces, geoms, line_of, step=10.0):
    """Lines per corridor piece. Every shape is sampled each `step` m and each sample votes only for the NEAREST
    parallel piece within TOL + 8 m, so a route lands on one corridor and never on a parallel street next to it.
    A piece carries a line when the line's votes cover at least half of it (short pieces: one vote)."""
    tree = STRtree(pieces)
    votes = defaultdict(lambda: defaultdict(float))
    for sid, g in geoms.items():
        l = line_of[sid]
        seen = defaultdict(float)
        for d in np.arange(0, g.length, step):
            pt = g.interpolate(d)
            h = heading_at(g, d)
            cand = [i for i in tree.query(pt.buffer(TOL + 8))
                    if pieces[i].distance(pt) <= TOL + 8
                    and ang_diff(heading_at(pieces[i], pieces[i].project(pt)), h) <= PARALLEL_DEG]
            if cand:
                i = min(cand, key=lambda i: pieces[i].distance(pt))
                seen[i] += step
        for i, v in seen.items():                       # a line counts once per piece (max over its shapes)
            votes[i][l] = max(votes[i][l], v)
    out = []
    for i, p in enumerate(pieces):
        need = min(0.5 * p.length, step)
        out.append({l for l, v in votes[i].items() if v >= max(need, 0.5 * p.length if p.length > 2 * step else step)})
    return out


def build_edges(pieces, sets):
    key = lambda c: (round(c[0], 1), round(c[1], 1))  # noqa: E731
    return chain_edges([dict(geom=p, lines=frozenset(s), u=key(p.coords[0]), v=key(p.coords[-1]))
                        for p, s in zip(pieces, sets) if s])


def chain_edges(edge_list):
    """Edges with the same lines chained through nodes of degree 2."""
    key = lambda c: (round(c[0], 1), round(c[1], 1))  # noqa: E731
    edges = dict(enumerate(edge_list))
    changed = True
    while changed:  # chain pieces with the same lines through nodes of degree 2
        changed = False
        inc = defaultdict(list)
        for i, e in edges.items():
            inc[e["u"]].append(i)
            inc[e["v"]].append(i)
        for n, es in inc.items():
            if len(es) == 2 and es[0] != es[1] and edges[es[0]]["lines"] == edges[es[1]]["lines"]:
                a, b = edges[es[0]], edges[es[1]]
                ca = list(a["geom"].coords) if a["v"] == n else list(a["geom"].coords)[::-1]
                cb = list(b["geom"].coords) if b["u"] == n else list(b["geom"].coords)[::-1]
                g = LineString(ca + cb[1:])
                edges[es[0]] = dict(geom=g, lines=a["lines"], u=key(g.coords[0]), v=key(g.coords[-1]))
                del edges[es[1]]
                changed = True
                break
    return list(edges.values())


# ---------------------------------------------------------------- 3b: continuous lines
def repair_lines(edges, geoms, line_of):
    """Every line drawn as one connected network ending where its routes end: pieces of a line separated by a
    corridor the membership test missed are joined (shortest corridor path, CONNECT_M), short unjoinable strays
    and dead-end stubs that are not route ends are dropped."""
    import networkx as nx
    # corridor ends a few metres apart that do not touch (a cut-end link landing next to a node) are linked
    # first, so a gap there is closed in place and not by a detour along another street
    edges = list(edges)
    ends = sorted({n for e in edges for n in (e["u"], e["v"])})
    tree = STRtree([Point(n) for n in ends])
    linked = {frozenset((e["u"], e["v"])) for e in edges}
    for a in ends:
        for k in tree.query(Point(a).buffer(SNAP_M)):
            b = ends[k]
            if b != a and frozenset((a, b)) not in linked and Point(a).distance(Point(b)) <= SNAP_M:
                linked.add(frozenset((a, b)))
                edges.append(dict(geom=LineString([a, b]), lines=frozenset(), u=a, v=b))
    G = nx.Graph()
    for i, e in enumerate(edges):
        if not G.has_edge(e["u"], e["v"]) or G[e["u"]][e["v"]]["w"] > e["geom"].length:
            G.add_edge(e["u"], e["v"], i=i, w=e["geom"].length)
    route_ends = defaultdict(list)
    for sid, g in geoms.items():
        route_ends[line_of[sid]] += [Point(g.coords[0]), Point(g.coords[-1])]
    edges = list(edges)
    lines = [set(e["lines"]) for e in edges]
    joined, dropped = 0, 0

    def sub(l):
        H = nx.Graph()
        for i, ls in enumerate(lines):
            if l in ls:
                H.add_edge(edges[i]["u"], edges[i]["v"], i=i)
        return H

    def comp_len(H, c):
        return sum(edges[H[a][b]["i"]]["geom"].length for a, b in H.subgraph(c).edges)

    for l in sorted(route_ends):
        for _ in range(50):
            H = sub(l)
            comps = sorted(nx.connected_components(H), key=lambda c: comp_len(H, c))
            if len(comps) <= 1:
                break
            c, others = comps[0], set().union(*comps[1:])
            dist, path = nx.multi_source_dijkstra(G, c, weight="w", cutoff=CONNECT_M)
            targets = [n for n in others if n in dist]
            if targets:
                p = path[min(targets, key=dist.get)]
                for a, b in zip(p, p[1:]):
                    lines[G[a][b]["i"]].add(l)
                joined += 1
            elif comp_len(H, c) < DROP_M:
                for a, b in H.subgraph(c).edges:
                    lines[H[a][b]["i"]].discard(l)
                dropped += 1
            else:
                break
        # loose ends inside the network (the line goes on along a corridor piece that starts a few metres away,
        # where carriageways / one-way pairs were folded): continue to the nearest other part of the line
        H = sub(l)
        for n in [n for n in H.nodes if H.degree(n) == 1]:
            if any(Point(n).distance(q) <= TOL + 40 for q in route_ends[l]):
                continue
            (m,) = list(H.neighbors(n))
            own = H[n][m]["i"]
            G2 = G.copy()
            if G2.has_edge(n, m) and G2[n][m]["i"] == own:
                G2.remove_edge(n, m)
            dist, path = nx.single_source_dijkstra(G2, n, weight="w", cutoff=END_JOIN_M)
            targets = [t for t in H.nodes if t not in (n, m) and t in dist]
            if targets:
                p = path[min(targets, key=dist.get)]
                for a, b in zip(p, p[1:]):
                    lines[G2[a][b]["i"]].add(l)
                joined += 1
        # two loose ends of the line close together (it arrives on one corridor strand and leaves on another that
        # does not touch it): bridge them
        H = sub(l)
        loose = [n for n in H.nodes if H.degree(n) == 1
                 and not any(Point(n).distance(q) <= TOL + 40 for q in route_ends[l])]
        done = set()
        for x in loose:
            if x in done:
                continue
            near = [y for y in loose if y != x and y not in done and Point(x).distance(Point(y)) <= PAIR_M]
            if near:
                y = min(near, key=lambda y: Point(x).distance(Point(y)))
                edges.append(dict(geom=LineString([x, y]), lines=frozenset(), u=x, v=y))
                lines.append({l})
                done |= {x, y}
                joined += 1
        changed = True
        while changed:
            changed = False
            H = sub(l)
            for n in [n for n in H.nodes if H.degree(n) == 1]:
                if any(Point(n).distance(q) <= TOL + 40 for q in route_ends[l]):
                    continue
                (m,) = list(H.neighbors(n))
                i = H[n][m]["i"]
                if edges[i]["geom"].length < STUB_M:
                    lines[i].discard(l)
                    dropped += 1
                    changed = True
                    break
    out = [dict(e, lines=frozenset(ls)) for e, ls in zip(edges, lines) if ls]
    print(f"line repair: {joined} gaps joined, {dropped} strays / stubs dropped")
    return out


def contract_short(edges, min_len=None, only=None, quiet=False):
    """Junction clusters joined by edges shorter than CONTRACT_M m become one junction (at their centre): lines
    then change order once, with room for the lane ramps, instead of curling over a few metres of corridor."""
    import networkx as nx
    min_len = min_len or CONTRACT_M
    G = nx.Graph()
    G.add_edges_from((e["u"], e["v"]) for k, e in enumerate(edges)
                     if e["geom"].length < min_len and e["u"] != e["v"] and (only is None or k in only))
    rep = {}
    for comp in nx.connected_components(G):
        c = tuple(np.round(np.array(list(comp)).mean(axis=0), 1))
        for n in comp:
            rep[n] = c
    out = []
    for e in edges:
        u, v = rep.get(e["u"], e["u"]), rep.get(e["v"], e["v"])
        if u == v and e["geom"].length < min_len:
            continue
        coords = list(e["geom"].coords)
        coords[0], coords[-1] = u, v
        out.append(dict(geom=LineString(coords), lines=e["lines"], u=u, v=v))
    if not quiet:
        print(f"junction clusters merged: {len(set(rep.values()))} ({len(edges) - len(out)} short edges removed)")
    return out


def straighten_roundabouts(edges, keep_near=()):
    """Every OSM roundabout the network touches becomes one node (as in cartogram.py): edges are cut RB_MARGIN_M
    outside the ring and joined by straight spokes to the line-weighted mean of the cut points; the ring itself
    goes. Then edges left running between the same two nodes (a frontage road beside the main carriageway, two
    carriageways into the ring, e.g. at Białka) are folded into one (fold_parallel). Roundabouts within
    800 m of a point in keep_near are left alone (Dworzec PKP: its schematic needs the real junctions)."""
    from cartogram import roundabouts
    key = lambda c: (round(c[0], 1), round(c[1], 1))  # noqa: E731
    n = 0
    for c, r in roundabouts():
        R = r + RB_MARGIN_M
        if not any(e["geom"].distance(c) < R for e in edges):
            continue
        if any(Point(s).distance(c) < 800 for s in keep_near):
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
                cut += [(xy[0] if k == 0 else xy[-1], len(e["lines"])) for k in (0, 1) if ends[k]]
                out.append((e, xy, ends))
        if not cut:
            edges = [o for o in out if isinstance(o, dict)]
            continue
        w = np.array([t for _, t in cut], float)
        m = tuple(np.round((np.array([xy for xy, _ in cut]) * w[:, None]).sum(axis=0) / w.sum(), 1))
        res = []
        for o in out:
            if isinstance(o, dict):
                res.append(o)
                continue
            e, xy, ends = o
            xy = ([m] if ends[0] else []) + xy + ([m] if ends[1] else [])
            res.append(dict(geom=LineString(xy), lines=e["lines"], u=key(xy[0]), v=key(xy[-1])))
        edges = res
    print(f"roundabouts straightened: {n}")
    return fold_parallel(edges)


def fold_parallel(edges):
    """Edges between the same two junctions with lengths within 1.5x (+50 m) of each other (two carriageways, a
    frontage road beside the main road, both approaches to a straightened roundabout) become one corridor: the
    one with more lines (then the shorter) carrying the lines of both. A detour much longer than the direct edge
    (the Dworzec PKP spur vs Sidorska) is kept apart. Repeated until stable: chaining after a fold can line up
    new pairs. Short edges closing a triangle (two similar edges from one junction to both its ends: a split that
    rejoins 50-80 m on, e.g. at the Białka roundabout) are contracted first, which turns the triangle into a pair."""
    total = 0
    while True:
        edges = collapse_triangles(edges)
        edges, n = _fold_once(edges)
        if not n:
            edges, n = fold_paths(edges)
        total += n
        if not n:
            break
    if total:
        print(f"parallel edges folded: {total}")
    return edges


def fold_paths(edges):
    """An edge running beside a path of up to PATH_FOLD_EDGES other edges between the same two junctions, each
    already carrying all its lines, of similar length and never more than PATH_FOLD_M away (Hausdorff), is dropped:
    the same link drawn twice, e.g. E's far-east variant out along the north road and up the spur past the
    cemetery, back down the spur and along DK2 50 m south. A loop round a block (A at Sławacinek) is farther apart
    and stays. Removes one edge per call: (edges, 0 / 1)."""
    import networkx as nx
    for k, e in sorted(enumerate(edges), key=lambda t: -t[1]["geom"].length):
        if e.get("raw") or e["u"] == e["v"]:
            continue
        G = nx.Graph()
        for i, f in enumerate(edges):
            if i != k and not f.get("raw") and f["u"] != f["v"] and e["lines"] <= f["lines"]:
                if not G.has_edge(f["u"], f["v"]) or G[f["u"]][f["v"]]["w"] > f["geom"].length:
                    G.add_edge(f["u"], f["v"], w=f["geom"].length, i=i)
        if e["u"] not in G or e["v"] not in G:
            continue
        try:
            path = nx.shortest_path(G, e["u"], e["v"], weight="w")
        except nx.NetworkXNoPath:
            continue
        if len(path) - 1 > PATH_FOLD_EDGES:
            continue
        geoms = [edges[G[a][b]["i"]]["geom"] for a, b in zip(path, path[1:])]
        if not similar(e["geom"].length, sum(g.length for g in geoms)):
            continue
        if e["geom"].hausdorff_distance(unary_union(geoms)) > PATH_FOLD_M:
            continue
        return chain_edges([f for i, f in enumerate(edges) if i != k]), 1
    return edges, 0


def similar(a, b, same_lines=False):
    """Lengths close enough for two edges between the same junctions to be one corridor. Edges carrying exactly the
    same lines (one link, its two directions on different roads, e.g. C on the DK2 / Janowska ramps: 289 vs
    556 m) fold up to 2.5x."""
    a, b = sorted((a, b))
    return b < (2.5 * a if same_lines else 1.5 * a + 50)


def collapse_triangles(edges):
    ends = defaultdict(dict)                         # node -> {neighbour: length}
    for e in edges:
        if e["u"] != e["v"]:
            ends[e["u"]][e["v"]] = ends[e["v"]][e["u"]] = e["geom"].length
    only = set()
    for k, e in enumerate(edges):
        u, v = e["u"], e["v"]
        if u == v or e["geom"].length >= TRIANGLE_M:
            continue
        if any(w not in (u, v) and w in ends[v] and similar(lu, ends[v][w]) for w, lu in ends[u].items()):
            only.add(k)
    return contract_short(edges, TRIANGLE_M, only, quiet=True) if only else edges


def _fold_once(edges):
    by_ends = defaultdict(list)
    out = []
    for e in edges:
        if e["u"] == e["v"]:
            out.append(e)
        else:
            by_ends[frozenset((e["u"], e["v"]))].append(e)
    folded = 0
    for es in by_ends.values():
        es = sorted(es, key=lambda e: (-len(e["lines"]), e["geom"].length))
        keep = [es[0]]
        for e in es[1:]:
            if similar(e["geom"].length, keep[0]["geom"].length, e["lines"] == keep[0]["lines"]):
                keep[0] = dict(keep[0], lines=keep[0]["lines"] | e["lines"])
                folded += 1
            else:
                keep.append(e)
        out += keep
    return chain_edges(out), folded


def drop_turning_stubs(edges, stations):
    """Short dead-end stubs (< TURN_STUB_M) whose stop lies nearer their junction than their dead end are the
    terminal turning loop (Terebelska Pętla: 55 m, stop 23 m from the junction, all six lines ending there turn
    into it). Drawn as a dummy straight street of TURN_DUMMY_M in the stub's direction with the stop moved near
    its end, so the lines that turn in are seen to turn and end flush at the capsule there, instead of a stub too
    short to read. A spur with its stop out at the end (Cmentarz Komunalny, 101 m) stays as it is."""
    deg = defaultdict(int)
    for e in edges:
        deg[e["u"]] += 1
        deg[e["v"]] += 1
    out, n = [], 0
    for e in edges:
        if e["geom"].length < TURN_STUB_M and (deg[e["u"]] == 1) != (deg[e["v"]] == 1):
            dead, junc = (e["u"], e["v"]) if deg[e["u"]] == 1 else (e["v"], e["u"])
            near = [k for k, s_ in stations.items() if Point(s_["xy"]).distance(e["geom"]) < 60]
            if near and all(Point(stations[k]["xy"]).distance(Point(junc)) < Point(stations[k]["xy"]).distance(Point(dead))
                            for k in near):
                d_ = np.array(dead) - np.array(junc)
                d_ /= np.hypot(*d_)
                end = tuple(np.round(np.array(junc) + d_ * TURN_DUMMY_M, 1))
                for k in near:
                    stations[k]["xy"] = np.array(junc) + d_ * (TURN_DUMMY_M - 25)
                out.append(dict(geom=LineString([junc, end]), lines=e["lines"], u=junc, v=end, dummy=True))
                n += 1
                continue
        out.append(e)
    print(f"terminal turning loops drawn as dummy streets: {n}")
    return chain_edges(out)


def edge_trips(edges, geoms, n_shape, line_of):
    """Trips of the day per line along every edge (shapes running along it, TOL m): e["trips"] = {line: n}."""
    for e in edges:
        g = e["geom"]
        probe = [g.interpolate(f, normalized=True) for f in (0.25, 0.5, 0.75)]
        t = defaultdict(int)
        for sid, sh in geoms.items():
            if line_of[sid] in e["lines"] and all(sh.distance(q) < TOL for q in probe):
                t[line_of[sid]] += int(n_shape[sid])
        e["trips"] = dict(t)
    return edges


def station_access(edges, stations):
    """Dworzec PKP drawn schematically as a U. In reality the buses run down Al. 1000-lecia from its crossing with
    Sidorska to a dead end at the station and turn back; D, F and H leave along Kościuszki instead. On the map a
    U hangs below the crossing: lines go down one leg, along the bottom through the station (a vertical bar
    across them) and up the other leg back to the crossing, where they split north / west / east as their routes
    do. Lines that leave along Kościuszki go straight on east from the U's bottom corner; lines that turn back at
    the station without passing (B) end at the bar."""
    key = lambda c: (round(c[0], 1), round(c[1], 1))  # noqa: E731
    out = [e for e in edges if e["u"] != e["v"] and e["geom"].length >= 1]
    for name in KEY_STATIONS:
        sts = [k for k, s in stations.items() if s["name"] == name]
        if not sts:
            continue
        xy = np.mean([stations[k]["xy"] for k in sts], axis=0)
        cand = [i for i, e in enumerate(out) if e["geom"].distance(Point(xy)) < 200]
        if not cand:
            continue
        i = min(cand, key=lambda i: out[i]["geom"].distance(Point(xy)))
        e, g = out[i], out[i]["geom"]
        d = g.project(Point(xy))
        if d < 50 or d > g.length - 50:
            continue
        # the shorter leg from the station is the access road (Al. 1000-lecia), the longer one Kościuszki
        legs = [(e["u"], [g.interpolate(x) for x in np.linspace(d, 0, 40)]),
                (e["v"], [g.interpolate(x) for x in np.linspace(d, g.length, 80)])]
        legs.sort(key=lambda t: LineString(t[1]).length)
        (N2, acc_pts), (N1, thr_pts) = legs
        thr_c = np.array([p.coords[0] for p in thr_pts])[::-1]      # roundabout -> station
        P0 = np.array(N2)
        dv = np.array(acc_pts[0].coords[0]) - P0
        dv /= np.hypot(*dv)                                           # down the access road
        perp = np.array([-dv[1], dv[0]])
        side = lambda q: np.sign(dv[0] * (q[1] - P0[1]) - dv[1] * (q[0] - P0[0]))  # noqa: E731
        if side(P0 + perp) != side(np.array(N1)):
            perp = -perp                                              # towards Kościuszki (east)
        # lines: those reaching Kościuszki (on the roundabout's other edges and not back on Sidorska) go on east
        at = lambda n: [f for k, f in enumerate(out) if k != i and n in (f["u"], f["v"])]  # noqa: E731
        L = e["lines"]
        far = lambda f, n: f["v"] if f["u"] == n else f["u"]  # noqa: E731
        direct = [f for f in at(N2) if far(f, N2) == N1]             # Sidorska straight to the roundabout
        on_direct = frozenset().union(*[f["lines"] for f in direct]) if direct else frozenset()
        n1_lines = frozenset().union(*[f["lines"] for f in at(N1) if f not in direct])
        kos_lines = (L & n1_lines) - on_direct                        # leave along Kościuszki
        end_lines = frozenset(l for l in L - kos_lines
                              if sum(1 for f in at(N2) if l in f["lines"]) == 1)  # in and out the same street: end
        u_lines = L - kos_lines - end_lines                           # round the U and back up
        trips = e.get("trips", {})
        tr = lambda ls: {l: trips.get(l, 0) for l in ls}  # noqa: E731
        # the U, drawn as is (no simplification): legs leave the crossing diagonally, run straight down, and
        # meet in a half circle whose lowest point is the station
        C = P0 + dv * (U_DEPTH_M - U_HALF_M)
        Sx = P0 + dv * U_DEPTH_M
        ss = lambda t: t * t * (3 - 2 * t)  # noqa: E731
        down = u_lines | kos_lines | end_lines

        def leg(sign, sh):
            # one long straight piece down from the crossing (so the arms turning into it get full-size arcs), a
            # smooth move outwards to the leg's place, straight down, a half circle at the bottom. Returns the points
            # and, per point, the bundle shift: at the top the two legs run side by side as one bundle continuing
            # Al. 1000-lecia (each pushed outwards by half the other's width), eased out as the legs move apart
            ease, lat = [], []
            for a in np.linspace(0, phi, 10)[1:]:              # first arc, bending outwards
                ease.append((U_S_R * math.sin(a), U_S_R * (1 - math.cos(a))))
            for a in np.linspace(phi, 0, 10)[1:]:              # second arc, back to straight down
                ease.append((2 * U_S_R * math.sin(phi) - U_S_R * math.sin(a),
                             2 * U_S_R * (1 - math.cos(phi)) - U_S_R * (1 - math.cos(a))))
            ease_p = [P0 + dv * (U_TOP_M + x) + sign * perp * y for x, y in ease]
            arc = [C + sign * perp * U_HALF_M * math.cos(a) + dv * U_HALF_M * math.sin(a)
                   for a in np.linspace(0, math.pi / 2, 19)]
            prof = [sh, sh] + [sh * (1 - y / U_HALF_M) for _, y in ease] + [0.0] * len(arc)
            return [tuple(p) for p in [P0, P0 + dv * U_TOP_M] + ease_p + arc], prof
        keep = [f for k, f in enumerate(out) if k != i]
        nL, nR = len(down), len(u_lines)
        phi = math.acos(1 - U_HALF_M / (2 * U_S_R))                   # S-bend: two arcs, U_HALF_M sideways
        s_len = 2 * U_S_R * math.sin(phi)                             # ... over this length
        # shift: metres left of the edge's direction (left leg: down, so west is negative; right leg: up, so east
        # is negative)
        lg, lprof = leg(-1, -nR * SPACING / 2)
        probe = U_TOP_M + s_len                                   # headings measured where the legs are apart
        keep.append(dict(geom=LineString([N2] + lg[1:-1] + [tuple(Sx)]), lines=down, u=N2, v=key(Sx), raw=True,
                         probe_m=probe, trips=tr(down), shift_prof=lprof, schematic=Point(xy)))
        if u_lines:
            rg, rprof = leg(1, -nL * SPACING / 2)
            keep.append(dict(geom=LineString([tuple(Sx)] + rg[::-1][1:-1] + [N2]), lines=u_lines, u=key(Sx), v=N2,
                             raw=True, probe_m=probe, trips=tr(u_lines), shift_prof=rprof[::-1], schematic=Point(xy)))
        if kos_lines:
            kc = int(np.argmax([LineString([thr_c[0], thr_c[-1]]).distance(Point(p)) for p in thr_c]))
            corner = thr_c[kc]
            turn = Sx + perp * max(float(np.dot(corner - Sx, perp)), U_HALF_M + 200)  # east along the bottom, then up
            keep.append(dict(geom=LineString([tuple(Sx), tuple(turn)] + [tuple(p) for p in thr_c[kc:0:-1]] + [N1]),
                             lines=kos_lines, u=key(Sx), v=N1, trips=tr(kos_lines),
                             schematic=LineString(thr_c).interpolate(0.5, normalized=True)))
            kos_g = keep[-1]["geom"]
            for st in stations.values():                              # stops on Kościuszki follow it
                if st["name"] != name and Point(st["xy"]).distance(LineString(thr_c)) < 40:
                    st["xy"] = np.array(kos_g.interpolate(kos_g.project(Point(st["xy"]))).coords[0])
        for k_ in sts:
            stations[k_]["xy"] = Sx - perp * 1.0                     # on the leg every line uses
        out = keep
        # how A / C / I pair their four sections at the crossing: by side, the more westerly arm with the west leg,
        # the other with the east leg, so no line crosses itself there (C, I: Sidorska west + west leg, Al.
        # 1000-lecia north + east leg; A: Al. 1000-lecia north + west leg, Sidorska east + east leg)
        idx = {id(f): k for k, f in enumerate(out)}
        left_k = idx[id(next(f for f in out if f.get("raw") and f["u"] == N2))]
        right_k = next((idx[id(f)] for f in out if f.get("raw") and f["v"] == N2), None)

        def east(k):
            g_ = out[k]["geom"]
            d_ = min(60, g_.length / 2)
            q = g_.interpolate(d_ if out[k]["u"] == N2 else g_.length - d_)
            return float(np.dot(np.array(q.coords[0]) - P0, perp))
        for l in u_lines:
            arms = [k for k, f in enumerate(out) if l in f["lines"] and N2 in (f["u"], f["v"])
                    and k not in (left_k, right_k)]
            if len(arms) != 2 or right_k is None:
                continue
            w_arm, e_arm = sorted(arms, key=east)
            PAIRS[(l, N2)] = {left_k: w_arm, w_arm: left_k, right_k: e_arm, e_arm: right_k}
        for k, f in enumerate(out):                                   # arms: shifted where they join the legs
            if N2 in (f["u"], f["v"]) and k not in (left_k, right_k):
                if f["geom"].length >= 2 * U_TOP_M:
                    f["split_at"] = (N2, U_TOP_M)
                else:                 # too short to shift back before the next junction (A on Sidorska east,
                    f["shift_all"] = True       # 445 m): shifted all along, so A runs straight between two arcs
        U_INFO.update(node=N2, left=left_k, right=right_k)
        NODE_W[key(Sx)] = STATION_CROSS_W
        NODE_W[N2] = ABOVE_CROSS_W
        print(f"{name}: U {sorted(u_lines)}, on along Kościuszki {sorted(kos_lines)}, ending at the station "
              f"{sorted(end_lines)}")
    return out


# ---------------------------------------------------------------- 4: line order (per edge)
def node_passes(edges):
    """For every node: incident edges with their heading away from the node and +1 / -1 (edge leaves / enters)."""
    nodes = defaultdict(list)
    for i, e in enumerate(edges):
        g = e["geom"]
        k = min(e.get("probe_m", 25), g.length / 2)        # the U's legs: measured where they have peeled apart
        if "probe_m" in e:
            pu, pv, qu, qv = g.interpolate(0), g.interpolate(g.length), g.interpolate(k), g.interpolate(g.length - k)
            nodes[e["u"]].append((i, math.atan2(qu.y - pu.y, qu.x - pu.x), +1))
            nodes[e["v"]].append((i, math.atan2(qv.y - pv.y, qv.x - pv.x), -1))
            continue
        nodes[e["u"]].append((i, heading_at(g, k), +1))
        nodes[e["v"]].append((i, heading_at(g, g.length - k) + math.pi, -1))
    return nodes


def left_of(orders, i, l, sgn):
    """Lateral lane position of line l on edge i, left positive, looking away from the node (sgn: +1 at u)."""
    ls = orders[i]
    return sgn * ((len(ls) - 1) / 2 - ls.index(l))


def turn(h_from_entry, h_out):
    """Signed turn angle (left positive, -pi..pi) from the direction of travel into the node to an exit."""
    return (h_out - h_from_entry + math.pi) % (2 * math.pi) - math.pi


def node_cost(inc, orders, edges, node=None):
    """Crossings at one node (pairs of passes whose left/right order is wrong) + SHIFT_COST per lane shifted,
    times the node's weight (NODE_W).
    A pass is a line going through the node from one section to another: a line with two sections here, or the
    pairs set in PAIRS (the Dworzec PKP U). Lines touching a node 3+ times otherwise (terminal loops) are left
    out: every order "crosses" somewhere there, which would drown the real junctions."""
    per_line = defaultdict(dict)
    for i, h, sgn in inc:
        for l in orders[i]:
            per_line[l][i] = (h, sgn)
    passes = []
    for l, ea in per_line.items():
        if len(ea) == 2:
            passes.append((l, *ea))
        elif (l, node) in PAIRS:
            pr = PAIRS[(l, node)]
            passes += [(l, p, q) for p, q in pr.items() if p < q and p in ea and q in ea]
        elif (l, node) in JOINS:             # 3+ sections: the joins the drawing makes (the routes' main relations)
            passes += [(l, p, q) for p, q in JOINS[(l, node)] if p in ea and q in ea]
    cost = 0.0
    for a, p, q in passes:                               # lane shifts of a line passing through
        ea = per_line[a]
        cost += SHIFT_COST * abs(-left_of(orders, p, a, ea[p][1]) - left_of(orders, q, a, ea[q][1]))
    crossings = 0
    for x in range(len(passes)):
        for y in range(x + 1, len(passes)):
            a, pa, qa_ = passes[x]
            b, pb, qb_ = passes[y]
            if a == b:
                continue
            ea, eb = per_line[a], per_line[b]
            if not {pa, qa_} & {pb, qb_}:                 # no section shared: cross if their ends interleave
                ang = lambda e_, ed: ed[e_][0] % (2 * math.pi)  # noqa: E731
                a0, a1 = sorted((ang(pa, ea), ang(qa_, ea)))
                inside = [a0 < ang(e_, eb) < a1 for e_ in (pb, qb_)]
                if inside[0] != inside[1]:
                    crossings += 1
                continue
            for p in {pa, qa_} & {pb, qb_}:
                qa = qa_ if p == pa else pa
                qb = qb_ if p == pb else pb
                hp, sp = ea[p]
                entry = hp + math.pi                     # direction of travel into the node along p
                la, lb = -left_of(orders, p, a, sp), -left_of(orders, p, b, sp)
                (ha, sa), (hb, sb) = ea[qa], eb[qb]
                if qa == qb:                             # on together: left/right must stay
                    if (la > lb) != (left_of(orders, qa, a, sa) > left_of(orders, qa, b, sb)):
                        crossings += 1
                else:                                    # split: the line turning further left keeps left
                    ta, tb = turn(entry, ha), turn(entry, hb)
                    if abs(ta - tb) > 0.05 and (la > lb) != (ta > tb):
                        crossings += 1
    return (cost + crossings) * NODE_W.get(node, 1.0)


def global_start(lines, edges, nodes, seed, iters=3000):
    """A good common order to start from (annealing over one order for all edges)."""
    rnd = random.Random(seed)
    def orders_for(order):
        rank = {l: k for k, l in enumerate(order)}
        return [sorted(e["lines"], key=rank.get) for e in edges]
    def total(order):
        o = orders_for(order)
        return sum(node_cost(inc, o, edges, n) for n, inc in nodes.items())
    cur = sorted(lines)
    cur_c = total(cur)
    best, best_c = cur[:], cur_c
    for it in range(iters):
        cand = cur[:]
        i, j = rnd.sample(range(len(cand)), 2)
        cand.insert(j, cand.pop(i))
        c = total(cand)
        T = 2.0 * (1 - it / iters) + 1e-3
        if c <= cur_c or rnd.random() < math.exp((cur_c - c) / T):
            cur, cur_c = cand, c
            if c < best_c:
                best, best_c = cand[:], c
    return orders_for(best), best_c


def best_orders(lines, edges, seed, iters=60000):
    """Per-edge line orders: start from the best common order, then anneal single edges (swap two lines on one
    edge; only the costs at its two end nodes change)."""
    nodes = node_passes(edges)
    orders, start_c = global_start(lines, edges, nodes, seed)
    ends = [(e["u"], e["v"]) for e in edges]
    cost_at = {n: node_cost(inc, orders, edges, n) for n, inc in nodes.items()}
    multi = [i for i, o in enumerate(orders) if len(o) > 1]
    rnd = random.Random(seed + 7)
    cur_total = sum(cost_at.values())
    best, best_total = [o[:] for o in orders], cur_total
    for it in range(iters):
        i = rnd.choice(multi)
        a_, b_ = rnd.sample(range(len(orders[i])), 2)
        orders[i][a_], orders[i][b_] = orders[i][b_], orders[i][a_]
        touched = set(ends[i])
        new = {n: node_cost(nodes[n], orders, edges, n) for n in touched}
        delta = sum(new[n] - cost_at[n] for n in touched)
        T = 1.0 * (1 - it / iters) + 1e-3
        if delta <= 0 or rnd.random() < math.exp(-delta / T):
            cost_at.update(new)
            cur_total += delta
            if cur_total < best_total - 1e-9:
                best, best_total = [o[:] for o in orders], cur_total
        else:
            orders[i][a_], orders[i][b_] = orders[i][b_], orders[i][a_]
    print(f"line order: common order cost {start_c:.1f} -> per-edge orders {best_total:.1f}")
    return best


# ---------------------------------------------------------------- 5: drawing geometry
def lane(orders, i, l):
    ls = orders[i]
    return ((len(ls) - 1) / 2 - ls.index(l)) * SPACING


def fillet(c, R):
    """Every corner of a polyline replaced by a circular arc (tangent length up to half of each leg), as in
    subway_like_map: a bundle offset from it bends as concentric arcs. Returns the points and, per point, the
    index of the input segment it belongs to (an arc's first half to the leg before the corner)."""
    c = np.asarray(c, float)
    if len(c) < 3:
        return c, list(range(len(c) - 1)) + [max(len(c) - 2, 0)]
    seg = np.linalg.norm(np.diff(c, axis=0), axis=1)
    out, sid = [c[0]], [0]
    for i in range(1, len(c) - 1):
        v1 = (c[i - 1] - c[i]) / (seg[i - 1] + 1e-9)
        v2 = (c[i + 1] - c[i]) / (seg[i] + 1e-9)
        th = math.pi - math.acos(float(np.clip(np.dot(v1, v2), -1, 1)))      # turning angle
        if th < math.radians(1) or th > math.radians(179):
            out.append(c[i])
            sid.append(i)
            continue
        tl = min(R * math.tan(th / 2), seg[i - 1] * 0.5, seg[i] * 0.5)
        r = tl / math.tan(th / 2)
        T1, T2 = c[i] + v1 * tl, c[i] + v2 * tl
        b = (v1 + v2) / (np.linalg.norm(v1 + v2) + 1e-9)
        C = c[i] + b * math.hypot(r, tl)
        a1, a2 = math.atan2(*(T1 - C)[::-1]), math.atan2(*(T2 - C)[::-1])
        da = (a2 - a1 + math.pi) % (2 * math.pi) - math.pi
        for f in np.linspace(0, 1, max(3, int(abs(da) / math.radians(5)) + 2)):
            out.append(C + r * np.array([math.cos(a1 + da * f), math.sin(a1 + da * f)]))
            sid.append(i - 1 if f < 0.5 else i)
    out.append(c[-1])
    sid.append(len(c) - 2)
    return np.array(out), sid


def section_coords(e, width):
    """An edge's centreline as drawn: simplified (SIMPLIFY_M), then bends inside it that are too tight for the
    bundle's width (a smaller arc would fold the inner lanes) removed, as subway_like_map does."""
    if e.get("raw"):                              # drawn schematically already (Dworzec PKP)
        return np.array(e["geom"].coords)
    c = np.array(e["geom"].simplify(SIMPLIFY_M).coords)
    rmin = (width - 1) / 2 * SPACING + 0.8 * SPACING
    while len(c) >= 3:
        seg = np.linalg.norm(np.diff(c, axis=0), axis=1)
        worst, wi = 1.0, None
        for i in range(1, len(c) - 1):
            v1 = (c[i - 1] - c[i]) / (seg[i - 1] + 1e-9)
            v2 = (c[i + 1] - c[i]) / (seg[i] + 1e-9)
            th = math.pi - math.acos(float(np.clip(np.dot(v1, v2), -1, 1)))
            if th < math.radians(1):
                continue
            ratio = rmin * math.tan(th / 2) / (min(seg[i - 1], seg[i]) * 0.5 + 1e-9)
            if ratio > worst:
                worst, wi = ratio, i
        if wi is None:
            break
        c = np.delete(c, wi, axis=0)
    if e.get("split_at"):                       # an extra vertex this far from that end (the arm's shifted piece)
        node, dist = e["split_at"]
        rev = tuple(np.round(c[0], 1)) != tuple(np.round(node, 1)) and             np.hypot(*(c[-1] - np.array(node))) < np.hypot(*(c[0] - np.array(node)))
        cc = c[::-1] if rev else c
        sd = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(cc, axis=0), axis=1))])
        if sd[-1] > dist + 60:
            k = int(np.searchsorted(sd, dist))
            f = (dist - sd[k - 1]) / (sd[k] - sd[k - 1])
            cc = np.insert(cc, k, cc[k - 1] + f * (cc[k] - cc[k - 1]), axis=0)
            c = cc[::-1] if rev else cc
    return c


def smooth(v, window):
    """Moving average along the first axis (ends held), applied twice: a soft, symmetric kernel."""
    k = max(1, int(round(window / STEP_M)))
    if k < 2 or len(v) < 3:
        return v
    pad = k
    for _ in range(2):
        vp = np.concatenate([np.repeat(v[:1], pad, axis=0), v, np.repeat(v[-1:], pad, axis=0)])
        c = np.cumsum(vp, axis=0)
        c = np.concatenate([np.zeros((1,) + v.shape[1:]), c])
        v = (c[k:] - c[:-k])[pad - k // 2: pad - k // 2 + len(v)] / k
    return v


def chain_seqs(l, edges):
    """Edges of line l chained through nodes where it has exactly two edges: [(edge, forward), ...]."""
    mine = [i for i, e in enumerate(edges) if l in e["lines"]]
    inc = defaultdict(list)
    for i in mine:
        inc[edges[i]["u"]].append(i)
        inc[edges[i]["v"]].append(i)
    used, chains = set(), []

    def walk(start_node, i):
        seq, node = [], start_node
        while i is not None and i not in used:
            used.add(i)
            e = edges[i]
            fwd = e["u"] == node
            seq.append((i, fwd))
            node = e["v"] if fwd else e["u"]
            nxt = [j for j in inc[node] if j != i and j not in used]
            pr = PAIRS.get((l, node), {})
            if i in pr:
                i = pr[i] if pr[i] not in used else None
            else:
                i = nxt[0] if len(inc[node]) == 2 and nxt else None
        return seq
    for n in [n for n, es in inc.items() if len(es) != 2 and not set(es) <= set(PAIRS.get((l, n), {}))]:
        for i in inc[n]:
            if i not in used:
                chains.append(walk(n, i))
    for i in mine:  # closed loops
        if i not in used:
            chains.append(walk(edges[i]["u"], i))
    return chains


def chain_node(step, edges, start):
    """The node a chain starts (start=True) or ends at, given its first / last (edge, forward) step."""
    i, fwd = step
    e = edges[i]
    return (e["u"] if fwd else e["v"]) if start else (e["v"] if fwd else e["u"])


def continuations(node, i, others, edges, shapes):
    """Trips (of the given shapes) running from / to section i through node on to each of the sections `others`
    (see lead_in): {section: trips}."""
    g_i = edges[i]["geom"]
    at = lambda d: g_i.interpolate(min(d, g_i.length / 2) if edges[i]["u"] == node  # noqa: E731
                                   else g_i.length - min(d, g_i.length / 2))
    qi, qi2 = at(40.0), at(150.0)
    score = {}
    for g, n in shapes:
        # the route must run ALONG section i (near both probes), not just sweep past its end at a merged junction
        # (K's Białka -> cemetery trips passed 40 m into Al. Solidarności at the Białka roundabout)
        if g.distance(qi) >= EXT_TOL or g.distance(qi2) >= EXT_TOL:
            continue
        p, pn = g.project(qi2), g.project(Point(node))
        if g.interpolate(pn).distance(Point(node)) > EXT_TOL:
            continue
        pt = g.interpolate(pn + (LEAD_FOLLOW_M if pn > p else -LEAD_FOLLOW_M))
        dist = {j: edges[j]["geom"].distance(pt) for j in others}
        j = min(dist, key=dist.get)
        if dist[j] < EXT_TOL and dist[j] < g_i.distance(pt):
            score[j] = score.get(j, 0) + n
    return score


def lead_in(l, node, i, edges, inc_l, before, shapes=()):
    """Line l's continuation through node from section i, as a chain step (None if l ends there): the section its
    buses really drive on to / come from, else the straightest. Every shape of the line running along section i
    (within EXT_TOL of a point 40 m out) is followed to where it passes the node and LEAD_FOLLOW_M beyond it (on
    the far side from i); the other section of the line it is nearest to there gets the shape's trips (measured
    from the node, not from the point on i: near a straightened roundabout that left K's Białka -> cemetery
    route landing on Al. Solidarności). Nearest, not "within a
    tolerance": next to a merged junction a route can run within 31 m of an arm it does not take (H at Sidorska /
    Al. Solidarności), and a tolerance either ties there or misses the routes elsewhere (L at Francuska).
    E.g. E from Daleka north into Wolska, H from Al. Solidarności west along Sidorska, L north to Francuska."""
    others = [j for j in inc_l[(l, node)] if j != i]
    if not others:
        return None
    if i in PAIRS.get((l, node), {}):
        j = PAIRS[(l, node)][i]
        return (j, edges[j]["v"] == node) if before else (j, edges[j]["u"] == node)

    def away(j):
        g = edges[j]["geom"]
        return heading_at(g, 0) if edges[j]["u"] == node else heading_at(g, g.length) + math.pi
    score = continuations(node, i, others, edges, shapes)
    h = away(i)
    straight = lambda j: abs((away(j) - h - math.pi + math.pi) % (2 * math.pi) - math.pi)  # noqa: E731
    j = min(others, key=lambda j: (-score.get(j, 0), straight(j)))
    return (j, edges[j]["v"] == node) if before else (j, edges[j]["u"] == node)


def extension_mask(Q, shapes, total, forced=None):
    """Per drawn point (centreline Q): True where the line runs with under EXT_SHARE of its trips of the day (an
    extension, drawn dashed), from the shapes passing within EXT_TOL m. Points no shape passes (schematic parts,
    rounded corners) take their neighbours' value; runs shorter than EXT_MIN_M m are merged into their neighbours."""
    pts = shapely.points(Q)
    cov = np.zeros(len(Q))
    for g, n in shapes:
        cov += n * (shapely.distance(pts, g) < EXT_TOL)
    val = np.where(cov > 0, (cov < EXT_SHARE * total).astype(float), np.nan)
    if forced is not None:                             # schematic parts (the Dworzec PKP U): from their own trips
        val = np.where(np.isnan(forced), val, forced)
    if np.isnan(val).all():
        return np.zeros(len(Q), bool)
    idx = np.where(~np.isnan(val))[0]                  # fill the gaps from the nearest covered point
    val = val[idx[np.clip(np.searchsorted(idx, np.arange(len(Q))), 0, len(idx) - 1)]]
    ext = val > 0.5
    sd = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(Q, axis=0).T))])
    for _ in range(100):
        cut = [0] + [k for k in range(1, len(ext)) if ext[k] != ext[k - 1]] + [len(ext)]
        runs = [(a, b) for a, b in zip(cut, cut[1:])]
        if len(runs) < 2:
            break
        a, b = min(runs, key=lambda r: sd[r[1] - 1] - sd[r[0]])
        if sd[b - 1] - sd[a] >= EXT_MIN_M:
            break
        ext[a:b] = ~ext[a]
    return ext


_SECTIONS = {}
PAIRS = {}      # (line, node) -> {edge: edge}: how a line touching a node 4 times pairs its sections (the U)
JOINS = {}      # (line, node) -> [(edge, edge)]: how a line with 3+ sections at a node is drawn through it
U_INFO = {}     # the Dworzec PKP U: crossing node, left / right leg edges
NODE_W = {}     # node -> weight of crossings there (crossings under the station bar are cheap)


def relax_junction_corners(edges, orders):
    """Corners AT junctions too tight for the bundles turning there (an arc with room for the innermost lane needs
    tangents of r * tan(angle / 2) on both sides, at most half of each piece): the piece's next vertex inside the
    edge is dropped, so the piece gets longer and the arc wider (e.g. the 9-line bundle from the south-west turning
    ~100 degrees onto Al. 1000-lecia, whose first piece was 195 m). Decided per edge, before drawing, so every
    line of a bundle uses the same centreline and the lanes stay concentric."""
    for i, e in enumerate(edges):
        if i not in _SECTIONS:
            _SECTIONS[i] = section_coords(e, len(orders[i]))
    at = defaultdict(list)
    for i, e in enumerate(edges):
        if e["u"] != e["v"]:
            at[e["u"]].append(i)
            at[e["v"]].append(i)

    def end(i, n):                                  # the section oriented to end at node n
        c = _SECTIONS[i]
        return c if np.hypot(*(c[-1] - np.array(n))) <= np.hypot(*(c[0] - np.array(n))) else c[::-1]
    dropped = 0
    for _ in range(10):
        changed = False
        for n, es in at.items():
            for a in es:
                for b in es:
                    if a >= b or not edges[a]["lines"] & edges[b]["lines"]:
                        continue
                    ca, cb = end(a, n), end(b, n)[::-1]          # into n along a, out of n along b
                    if len(ca) < 2 or len(cb) < 2:
                        continue
                    v1, v2 = ca[-1] - ca[-2], cb[1] - cb[0]
                    l1, l2 = np.hypot(*v1), np.hypot(*v2)
                    th = math.acos(float(np.clip(np.dot(v1, v2) / (l1 * l2 + 1e-9), -1, 1)))
                    if th < math.radians(10) or th > math.radians(160):
                        continue
                    w = max(len(orders[a]), len(orders[b]))
                    need = 2 * ((w - 1) / 2 * SPACING + 0.8 * SPACING) * math.tan(th / 2)
                    for k, c_, l_ in ((a, ca, l1), (b, cb[::-1], l2)):
                        e_ = edges[k]
                        if l_ >= need or len(c_) < 3 or e_.get("raw") or (e_.get("split_at") or (None,))[0] == n:
                            continue
                        c2 = np.delete(c_, len(c_) - 2, axis=0)          # the vertex next to n
                        _SECTIONS[k] = c2 if np.array_equal(c_[0], _SECTIONS[k][0]) else c2[::-1]
                        dropped += 1
                        changed = True
        if not changed:
            break
    print(f"junction corners relaxed: {dropped} vertices dropped")


def set_joins(edges, shapes):
    """JOINS for the order search: where a line has 3+ sections at a node (a loop end, a branch), the pairs the
    drawing joins (lead_in: the routes' main relations) count as passes, so crossings there are counted. They were
    left out, which let K and E swap for free at the Białka roundabout."""
    inc_l = defaultdict(list)
    for i, e in enumerate(edges):
        for l in e["lines"]:
            inc_l[(l, e["u"])].append(i)
            inc_l[(l, e["v"])].append(i)
    for (l, n), secs in inc_l.items():
        if len(secs) < 3 or (l, n) in PAIRS:
            continue
        pairs = set()
        for i in secs:
            s_ = lead_in(l, n, i, edges, inc_l, True, shapes.get(l, []))
            if s_:
                pairs.add(tuple(sorted((i, s_[0]))))
        JOINS[(l, n)] = sorted(pairs)


def arm_shifts(edges, orders):
    """Arms meeting the U at the crossing approach it already shifted to where their lines go on the legs (lane
    offsets left of travel kept through the turn), so every lane turns on a concentric arc instead of folding
    round the crossing point: e["end_shift"] = (node, metres left of the edge's direction)."""
    if not U_INFO:
        return
    N2, lk, rk = U_INFO["node"], U_INFO["left"], U_INFO["right"]

    def leg_x(k, l):                            # lane position at the crossing, across Al. 1000-lecia (east +)
        if k == lk:
            return lane(orders, k, l) + edges[k]["shift_prof"][0]
        return -(lane(orders, k, l) + edges[k]["shift_prof"][-1])
    for k, e in enumerate(edges):
        if "split_at" not in e and not e.get("shift_all"):
            continue
        d = []
        for l in e["lines"]:
            pr = PAIRS.get((l, N2))
            leg = pr[k] if pr and k in pr else lk if l in edges[lk]["lines"] else rk if l in edges[rk]["lines"] else None
            if leg is None:
                continue
            lo = lane(orders, k, l) * (1 if e["v"] == N2 else -1)   # left of travel towards the crossing
            d.append(leg_x(leg, l) - lo)
        if d:
            sh = float(np.mean(d))
            e["end_shift"] = (N2, sh if e["v"] == N2 else -sh)


def draw_path(seq, l, edges, orders):
    """Drawn coordinates of one chain (subway_like_map style): each section simplified, corners (inside sections
    and at junctions) replaced by circular arcs of radius CORNER_R, the lane offset eased over RAMP_M m wherever
    it changes, then applied along the normals, so a bundle bends as concentric arcs."""
    pts, seg_off, seg_edge = [], [], []
    for i, fwd in seq:
        if i not in _SECTIONS:
            _SECTIONS[i] = section_coords(edges[i], len(orders[i]))
        c = _SECTIONS[i] if fwd else _SECTIONS[i][::-1]
        o = lane(orders, i, l) * (1 if fwd else -1)
        offs = np.full(len(c), o)                   # per point: offset of the segment ending there
        if "end_shift" in edges[i] and len(c) > 1:  # an arm of the U: its last piece shifted (short arms: all)
            node, sh = edges[i]["end_shift"]
            at_start = (edges[i]["u"] if fwd else edges[i]["v"]) == node
            if edges[i].get("shift_all"):
                offs = offs + sh * (1 if fwd else -1)
            else:
                offs[1 if at_start else len(c) - 1] += sh * (1 if fwd else -1)
        if "shift_prof" in edges[i]:                # whole bundle shifted (the U's legs, side by side at the top)
            pr = np.array(edges[i]["shift_prof"], float)
            offs = offs + (pr if fwd else -pr[::-1])
        if pts:
            c, offs = c[1:], offs[1:]
        pts.extend(c)
        seg_off.extend(offs)
        seg_edge.extend([i] * len(c))
    P = np.array(pts)
    seg_off = np.array(seg_off[1:] if len(seg_off) > 1 else seg_off, dtype=float)  # offset of segment k: P[k]-P[k+1]
    seg_edge = np.array(seg_edge[1:] if len(seg_edge) > 1 else seg_edge)
    if len(P) < 2:
        return P, np.full(len(P), seq[0][0])
    F, sid = fillet(P, CORNER_R)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(F, axis=0), axis=1))])
    keep = np.concatenate([[True], np.diff(s) > 1e-6])
    F, s, sid = F[keep], s[keep], np.array(sid)[keep]
    if s[-1] < 1:
        return F, np.full(len(F), seq[0][0])
    ds = np.append(np.arange(0, s[-1], STEP_M), s[-1])
    Q = np.column_stack([np.interp(ds, s, F[:, 0]), np.interp(ds, s, F[:, 1])])
    k_ = np.clip(sid[np.clip(np.searchsorted(s, ds, side="right") - 1, 0, len(sid) - 1)], 0, len(seg_off) - 1)
    O, E = seg_off[k_], seg_edge[k_]
    O = smooth(O[:, None].astype(float), RAMP_M)[:, 0]
    t = np.gradient(Q, axis=0)
    t = smooth(t, STEP_M * 3)
    t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-9)[:, None]
    nrm = np.stack([-t[:, 1], t[:, 0]], axis=1)
    return Q + nrm * O[:, None], E, Q


def station_spots(stations, edges):
    """Station -> (edge index, point on the corridor, unit tangent, lines of the station on that edge)."""
    geoms = [e["geom"] for e in edges]
    tree = STRtree(geoms)
    deg = defaultdict(int)
    for e in edges:
        deg[e["u"]] += 1
        deg[e["v"]] += 1
    # short dead-end stubs (a terminal turning loop, e.g. the 80 m one at Terebelska Pętla, which runs ACROSS the
    # bundle) are no place for a station capsule when a through edge is near - unless the stop is out at the
    # stub's dead end (Cmentarz Komunalny, at the end of its spur)
    def stub(i, pt):
        e = edges[i]
        if geoms[i].length >= 200 or (deg[e["u"]] != 1 and deg[e["v"]] != 1):
            return False
        dead, junc = (e["u"], e["v"]) if deg[e["u"]] == 1 else (e["v"], e["u"])
        return pt.distance(Point(junc)) < pt.distance(Point(dead))
    spots = {}
    for sid, s in stations.items():
        pt = Point(s["xy"])
        cand = [i for i in tree.query(pt.buffer(90)) if edges[i]["lines"] & s["lines"]]
        cand = [i for i in cand if not stub(i, pt)] or cand
        if not cand:
            continue
        i = min(cand, key=lambda i: geoms[i].distance(pt))
        g = geoms[i]
        d = g.project(pt)
        h = heading_at(g, d)
        spots[sid] = (i, np.array(g.interpolate(d).coords[0]), np.array([math.cos(h), math.sin(h)]),
                      s["lines"] & edges[i]["lines"])
    return spots


# ---------------------------------------------------------------- render
# grey street names along the main bus streets: OSM name -> name on the map
STREET_LABELS = {
    "Sidorska": "Sidorska", "Aleja Tysiąclecia": "al. Tysiąclecia", "Gabriela Narutowicza": "Narutowicza",
    "Zamkowa": "Zamkowa", "Terebelska": "Terebelska", "Janowska": "Janowska", "Warszawska": "Warszawska",
    "Brzeska": "Brzeska", "Grzybowa": "Grzybowa", "Aleja Jana Pawła II": "al. Jana Pawła II",
    "Orzechowa": "Orzechowa", "Aleja Solidarności": "al. Solidarności", "Akademicka": "Akademicka",
    "Łomaska": "Łomaska", "Francuska": "Francuska", "Białka": "Białka", "Królowej Jadwigi": "Królowej Jadwigi",
    "Podleśna": "Podleśna", "Lubelska": "Lubelska",
}
STREET_SIZE, STREET_COLOR = 10.5, "#8e8e8e"


def street_labels(ax, edges, orders, occ, occ_tree, placed, frame, m_per_pt):
    """Small grey names of the main bus streets (STREET_LABELS), along a straight stretch of the street's bundle,
    rotated with it, just outside its outermost lane: below a (near) horizontal street, right of a (near) vertical
    one where that is free, else on the other side; clear of every line, capsule, stop name and other street name
    (a street with no free straight stretch gets no name). Streets are found from the OSM road names along each
    drawn section."""
    import osmnx as ox
    G = ox.project_graph(ox.load_graphml(ROOT / "data" / "cache" / "osm_drive_service.graphml"), to_crs=CRS)
    o = ox.graph_to_gdfs(G, nodes=False).reset_index()
    o = o[o["name"].notna()]
    names = [n[0] if isinstance(n, list) else n for n in o["name"]]
    keep = [k for k, n in enumerate(names) if n in STREET_LABELS]
    otree = STRtree([o.geometry.iat[k] for k in keep])
    # candidate spots: every 30 m along each drawn section where the street under it is a labelled one
    cands = defaultdict(list)                 # street -> [(edge, point on the section, unit tangent, half-width)]
    for i, e in enumerate(edges):
        if e.get("raw") or i not in _SECTIONS:
            continue
        sec = LineString(_SECTIONS[i])
        if sec.length < 200:
            continue
        half = len(orders[i]) * SPACING / 2
        for dd in np.arange(80, sec.length - 80, 30):
            p = sec.interpolate(dd)
            k = otree.nearest(p)
            if otree.geometries[k].distance(p) > 30:
                continue
            a, b = sec.interpolate(dd - 70), sec.interpolate(dd + 70)
            a2, b2 = sec.interpolate(dd - 70), sec.interpolate(dd)
            t = np.array([b.x - a.x, b.y - a.y])
            t2 = np.array([b2.x - a2.x, b2.y - a2.y])
            if np.hypot(*t) < 1:
                continue
            t, t2 = t / np.hypot(*t), t2 / max(np.hypot(*t2), 1e-9)
            if np.dot(t, t2) < math.cos(math.radians(6)):        # a bend here: not straight enough for text
                continue
            cands[names[keep[k]]].append((i, np.array(p.coords[0]), t, half, dd, sec.length))
    m = m_per_pt
    n_done = 0
    for street, cs in sorted(cands.items(), key=lambda kv: -len(kv[1])):
        text = STREET_LABELS[street]
        w = len(text) * STREET_SIZE * 0.55 * m
        h = STREET_SIZE * 1.15 * m
        # spots nearest the middle of the street's longest section first; the preferred side (below a horizontal
        # street, right of a vertical one) anywhere along the street before the other side
        cs = sorted(cs, key=lambda c_: (-c_[5], abs(c_[4] - c_[5] / 2)))
        done = False
        for pref in (True, False):
          for i, p, t, half, _, _ in cs:  # noqa: E111
            if t[0] < 0 or (abs(t[0]) < 1e-6 and t[1] < 0):
                t = -t                                           # text reads left to right (or bottom-up)
            nrm = np.array([-t[1], t[0]])                        # left of the reading direction (= up)
            horizontal = abs(t[0]) >= abs(t[1])
            sides = [-1, 1] if horizontal else ([-1, 1] if t[1] > 0 else [1, -1])
            for sd_ in sides[:1] if pref else sides[1:]:
                cen = p + nrm * sd_ * (half + SPACING * 0.5 + h / 2)
                poly = Polygon([cen + t * w / 2 + nrm * h / 2, cen - t * w / 2 + nrm * h / 2,
                                cen - t * w / 2 - nrm * h / 2, cen + t * w / 2 - nrm * h / 2])
                if not frame.contains(poly):
                    continue
                if any(occ[k].intersects(poly) for k in occ_tree.query(poly)):
                    continue
                if any(q.buffer(12).intersects(poly) for q in placed):
                    continue
                ang = math.degrees(math.atan2(t[1], t[0]))
                ax.text(cen[0], cen[1], text, fontproperties=font(STREET_SIZE), color=STREET_COLOR, ha="center",
                        va="center", rotation=ang, rotation_mode="anchor", zorder=6.5,
                        path_effects=[patheffects.withStroke(linewidth=2.5, foreground=LAND)])
                placed.append(poly)
                done = True
                break
            if done:
                break
          if done:  # noqa: E111
            break
        n_done += done
    print(f"street names: {n_done} of {len(STREET_LABELS)} placed")


EXT_STYLE = "pale"     # how extensions are drawn: dashed / thin / pale / hollow (--ext-style)


def draw_line(ax, xy, color, lw, ext, zorder=3):
    """One run of a line: solid, or as an extension in EXT_STYLE (dashed; thin: solid at 45 % width; pale: full
    width in the colour washed 55 % towards the paper; hollow: coloured outline with a paper-coloured centre)."""
    kw = dict(solid_capstyle="butt", solid_joinstyle="round", dash_capstyle="butt")
    x, y = xy[:, 0], xy[:, 1]
    if ext and zorder == 3:
        zorder = 2.9                                # extensions under the full-service lines
    if not ext:
        ax.plot(x, y, color=color, lw=lw, zorder=zorder, **kw)
    elif EXT_STYLE == "dashed":
        ax.plot(x, y, color=color, lw=lw, linestyle=(0, (DASH_ON, DASH_OFF)), zorder=zorder, **kw)
    elif EXT_STYLE == "thin":
        ax.plot(x, y, color=color, lw=lw * 0.45, zorder=zorder, **kw)
    elif EXT_STYLE == "pale":
        c = np.array(matplotlib.colors.to_rgb(color))
        ax.plot(x, y, color=c + (np.array(matplotlib.colors.to_rgb(LAND)) - c) * 0.55, lw=lw, zorder=zorder, **kw)
    elif EXT_STYLE == "hollow":
        ax.plot(x, y, color=color, lw=lw, zorder=zorder, **kw)
        ax.plot(x, y, color=LAND, lw=lw * 0.45, zorder=zorder + 0.01, **kw)


def render(edges, orders, stations, line_cfg, date, out, shapes, line_trips):
    order = sorted({l for e in edges for l in e["lines"]})
    limits = gpd.read_file(LIMITS).to_crs(CRS)
    x0, y0, x1, y1 = limits.total_bounds
    x0, y0, x1, y1 = x0 - 300, y0 - 300, x1 + 300, y1 + 300
    panel = (x1 - x0) * 0.36                      # legend panel on the east
    X0, X1 = x0, x1 + panel
    W_IN = 24.0
    H_IN = W_IN * (y1 - y0) / (X1 - X0)
    fig = plt.figure(figsize=(W_IN, H_IN), facecolor=LAND)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(X0, X1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.axis("off")
    pt_per_m = W_IN * 72 / (X1 - X0)
    lw = SPACING * 0.82 * pt_per_m

    clip = gpd.GeoSeries([box(X0, y0, X1, y1).buffer(500)], crs=CRS).to_crs(4326).total_bounds.tolist()
    for name, col in [("forests", GREEN), ("parks", GREEN), ("water", WATER), ("waterways", WATER)]:
        f = BASEMAP / f"{name}.gpkg"
        if f.exists():
            g = gpd.read_file(f, bbox=clip).to_crs(CRS)
            g = g[g.geometry.notna()]
            poly = g[g.geom_type.isin(["Polygon", "MultiPolygon"])]
            if len(poly):
                poly.plot(ax=ax, color=col, linewidth=0, zorder=1, aspect=None)
            ln = g[g.geom_type.isin(["LineString", "MultiLineString"])]
            if len(ln) and name == "waterways":
                ln.plot(ax=ax, color=col, linewidth=1.2, zorder=1, aspect=None)
    limits.boundary.plot(ax=ax, color="#9a9a9a", linewidth=1.0, linestyle=(0, (5, 4)), zorder=1.5, aspect=None)

    colors = {l: c for l, c in zip([g[0] for g in line_cfg["groups"]], line_cfg["colors"])}
    inc_l = defaultdict(list)                     # (line, node) -> sections of the line there
    for i, e in enumerate(edges):
        for l in e["lines"]:
            inc_l[(l, e["u"])].append(i)
            inc_l[(l, e["v"])].append(i)
    drawn = []                                    # drawn lines / capsules (map metres), kept clear of labels
    paths = defaultdict(list)                     # line -> drawn LineStrings
    line_ends = []                                # (line, end point, direction pointing out of the line)
    spots = station_spots(stations, edges)
    term_cut = defaultdict(list)                  # terminus station -> points where lines ending there were cut

    def trim_at_terminus(l, c, E, Q, at_end):
        """A line ending at one of its main termini is cut where it crosses the station (the line through the
        stop across the bundle), so every line ending there stops flush at the station capsule instead of
        wherever its drawing happens to end (Terebelska Pętla had six end bars scattered round a corner).
        Returns (c, E, Q, station or None)."""
        xy = c[-1] if at_end else c[0]
        cand = [k for k, s_ in stations.items() if l in s_["term"] and k in spots
                and np.hypot(*(s_["xy"] - xy)) < 300]
        if not cand:
            return c, E, Q, None
        k = min(cand, key=lambda k: np.hypot(*(stations[k]["xy"] - xy)))
        _, base, t, _ = spots[k]
        nrm = np.array([-t[1], t[0]])
        sgn = np.sign((c - base) @ t)
        idx = range(len(c) - 1, 0, -1) if at_end else range(0, len(c) - 1)
        sd = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(c, axis=0).T))])
        for j in idx:
            a, b = (j - 1, j) if at_end else (j, j + 1)
            if sgn[a] != sgn[b] and sgn[a] != 0:
                sa, sb = (c[a] - base) @ t, (c[b] - base) @ t
                p = c[a] + (c[b] - c[a]) * (sa / (sa - sb))
                far = (sd[-1] - sd[a]) if at_end else sd[b]
                if abs((p - base) @ nrm) > 4 * SPACING + 40 or far > 400:
                    return c, E, Q, None
                if at_end:
                    c, E, Q = np.vstack([c[:b], p]), np.append(E[:b], E[a]), np.vstack([Q[:b], Q[a]])
                else:
                    c, E, Q = np.vstack([p, c[b:]]), np.append(E[b], E[b:]), np.vstack([Q[b], Q[b:]])
                term_cut[k].append(p)
                return c, E, Q, k
        # ends just short of the station (I at Chłodna 1: inbound ends 17 m before it, outbound starts 17 m after):
        # carried on along the end's direction to the station line, so both meet at the capsule
        tip, prev = (c[-1], c[-2]) if at_end else (c[0], c[1])
        dv_ = tip - prev
        dv_ = dv_ / max(np.hypot(*dv_), 1e-9)
        st_ = -((tip - base) @ t) / (dv_ @ t) if abs(dv_ @ t) > 0.3 else -1
        if 0 <= st_ <= 60:
            p = tip + dv_ * st_
            if abs((p - base) @ nrm) <= 4 * SPACING + 40:
                if at_end:
                    c, E, Q = np.vstack([c, p]), np.append(E, E[-1]), np.vstack([Q, Q[-1]])
                else:
                    c, E, Q = np.vstack([p, c]), np.append(E[0], E), np.vstack([Q[0], Q])
                term_cut[k].append(p)
                return c, E, Q, k
        return c, E, Q, None

    for l in order:
        for seq in chain_seqs(l, edges):
            # a chain ending at a junction where the line goes on along other sections (a branch, a crossing of its
            # own) is drawn with a lead-in of its straightest continuation: it eases out of that lane, no gap
            pre = lead_in(l, chain_node(seq[0], edges, True), seq[0][0], edges, inc_l, True, shapes.get(l, []))
            post = lead_in(l, chain_node(seq[-1], edges, False), seq[-1][0], edges, inc_l, False, shapes.get(l, []))
            own = {i for i, _ in seq}
            pre = pre if pre and pre[0] not in own else None
            post = post if post and post[0] not in own else None    # may equal pre (a loop leaving and
            # returning through the same junction, L at Francuska): drawing that short lead-in twice is harmless
            c, E, Q = draw_path(([pre] if pre else []) + seq + ([post] if post else []), l, edges, orders)
            if len(c) < 2:
                continue
            sd = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(c, axis=0).T))])
            mask = np.ones(len(c), bool)
            if pre:
                k0 = int(np.argmax(E != pre[0]))
                mask &= sd >= sd[k0] - RAMP_M
            if post:
                k1 = len(E) - 1 - int(np.argmax(E[::-1] != post[0]))
                mask &= sd <= sd[k1] + RAMP_M
            c, E, Q = c[mask], E[mask], Q[mask]
            if len(c) < 2:
                continue
            trimmed = {}
            for end_i, step in ((0, seq[0]), (-1, seq[-1])):
                node = chain_node(step, edges, end_i == 0)
                if len(inc_l[(l, node)]) == 1 and not (end_i == 0 and pre) and not (end_i == -1 and post):
                    c, E, Q, trimmed[end_i] = trim_at_terminus(l, c, E, Q, end_i == -1)
            if len(c) < 2:
                continue
            # schematic parts (the Dworzec PKP U, Kościuszki): solid / pale from the trips passing the real place
            # they stand for (the station stop; the middle of Kościuszki), not from shapes near the drawing
            sch = {}
            for i in set(E.tolist()):
                if edges[i].get("schematic") is not None:
                    n_ = sum(n for g, n in shapes.get(l, []) if g.distance(edges[i]["schematic"]) < EXT_TOL)
                    sch[i] = float(n_ < EXT_SHARE * line_trips[l])
            forced = np.array([sch.get(i, np.nan) for i in E])
            ext = extension_mask(Q, shapes.get(l, []), line_trips[l], forced)
            cut = [0] + [k for k in range(1, len(c)) if ext[k] != ext[k - 1]] + [len(c)]
            for a_, b_ in zip(cut, cut[1:]):
                part = c[max(a_ - 1, 0):b_]
                if len(part) < 2:
                    continue
                draw_line(ax, part, colors.get(l, "#888"), lw, ext[a_])
            g = LineString(c)
            drawn.append(g)
            paths[l].append(g)
            # chain ends at nodes where the line has a single edge are line ends
            for end_i, (i, fwd), node_key in ((0, seq[0], "u"), (-1, seq[-1], "v")):
                e = edges[i]
                node = e[node_key] if fwd else e["v" if node_key == "u" else "u"]
                if sum(1 for e2 in edges if l in e2["lines"] and node in (e2["u"], e2["v"])) == 1:
                    k = 0 if end_i == 0 else len(c) - 1
                    k2 = min(len(c) - 1, 4) if end_i == 0 else max(0, len(c) - 5)
                    d = c[k] - c[k2]
                    line_ends.append((l, c[k], d / max(np.hypot(*d), 1e-9), trimmed.get(end_i)))

    # stops, interchanges, labels (one name per stop name without the post number)
    halo = [patheffects.withStroke(linewidth=4, foreground=LAND)]
    m_per_pt = 1 / pt_per_m
    base_name = lambda n: re.sub(r"\s+\d+[A-Za-z]?(\s*\(.*\))?$", "", n).strip()  # noqa: E731

    def capsule(p, q):
        ax.plot([p[0], q[0]], [p[1], q[1]], color=INK, lw=lw * 1.3, solid_capstyle="round", zorder=5)
        ax.plot([p[0], q[0]], [p[1], q[1]], color="white", lw=lw * 0.85, solid_capstyle="round", zorder=5.1)
        drawn.append(LineString([p, q]))

    def on_line(l, xy):
        """Point of line l's drawing nearest xy (None if the line is not drawn near)."""
        best = None
        for g in paths.get(l, []):
            d = g.distance(Point(xy))
            if d < 3 * SPACING + 30 and (best is None or d < best[0]):
                best = (d, np.array(g.interpolate(g.project(Point(xy))).coords[0]))
        return None if best is None else best[1]

    def drawn_tangent(l, xy):
        """Unit tangent of line l's drawing at its point nearest xy (sign as drawn)."""
        g = min(paths[l], key=lambda g: g.distance(Point(xy)))
        d = g.project(Point(xy))
        a, b = g.interpolate(max(d - 8, 0)), g.interpolate(min(d + 8, g.length))
        v = np.array([b.x - a.x, b.y - a.y])
        return v / max(np.hypot(*v), 1e-9)

    def one_way(l, base, t):
        """+1 / -1 if line l passes base only along / against the corridor tangent t (all its shapes of the day
        within EXT_TOL m run that way: a loop, a one-way pair of streets), 0 if both ways."""
        fw = bw = 0
        zone = Point(base).buffer(EXT_TOL)
        for g, n in shapes.get(l, []):
            if g.distance(Point(base)) >= EXT_TOL:
                continue
            for part in parts(g.intersection(zone)):             # every pass (a shape may pass twice)
                v = np.array(part.coords[-1]) - np.array(part.coords[0])
                if np.hypot(*v) < 10:
                    continue
                c = float(np.dot(v / np.hypot(*v), t))
                fw += n * (c > 0.5)
                bw += n * (c < -0.5)
        return 1 if fw and not bw else -1 if bw and not fw else 0

    # line ends: at a main terminus of the line an end bar, at the end of an extended variant a white diamond
    ends_at = defaultdict(set)                    # station -> lines whose drawing ends there
    var_end = defaultdict(set)                    # station -> lines whose extended variant ends there
    for l, xy, d, cut_at in line_ends:
        if cut_at is not None:                    # cut flush at its terminus: the station capsule marks the end
            ends_at[cut_at].add(l)
            continue
        cand = [k for k, s_ in stations.items() if l in s_["lines"] and np.hypot(*(s_["xy"] - xy)) < 250]
        k = min(cand, key=lambda k: np.hypot(*(stations[k]["xy"] - xy))) if cand else None
        if k is not None:
            ends_at[k].add(l)
        # no end mark where the line turns at Dworzec PKP (it gets its own capsule)
        if any(s_["name"] in KEY_STATIONS and np.hypot(*(s_["xy"] - xy)) < 80 for s_ in stations.values()):
            continue
        if k is not None and l in stations[k]["term"]:
            nrm = np.array([-d[1], d[0]])
            capsule(xy - nrm * SPACING * 0.5, xy + nrm * SPACING * 0.5)
        else:
            if k is not None:
                var_end[k].add(l)
            ax.scatter([xy[0]], [xy[1]], marker="D", s=(lw * 1.45) ** 2, facecolor="white",
                       edgecolor=colors.get(l, "#888"), linewidths=lw * 0.38, zorder=5.2)
            drawn.append(Point(xy).buffer(SPACING * 0.6))

    # names: termini only (main termini of a line with round bullets, ends of extended variants with diamond
    # bullets, and Dworzec PKP); every other stop is a dot on its lines
    requests = []                                 # labels to place
    n_arrows = 0
    cap_r = lw * 1.3 / 2 * m_per_pt
    capsuled = set()       # one capsule per stop name and place (Twarda 1, Chłodna 1: two stations side by side)
    for sid, (i, base, t, ls) in sorted(spots.items(), key=lambda kv: -len(stations[kv[0]]["lines"])):
        st_ = stations[sid]
        nrm = np.array([-t[1], t[0]])
        key_st = st_["name"] in KEY_STATIONS
        if not (st_["term"] or key_st or var_end.get(sid)):
            for l in ls:
                p = on_line(l, base)
                if p is None:
                    continue
                w = one_way(l, base, t)
                if w:                                    # served one way only: an arrowhead in the direction of travel
                    td = drawn_tangent(l, p)
                    dirv = td * w * (1 if np.dot(td, t) >= 0 else -1)
                    ang = math.degrees(math.atan2(dirv[1], dirv[0])) - 90
                    ax.scatter([p[0]], [p[1]], marker=(3, 0, ang), s=(lw * 1.35) ** 2, color=INK, zorder=4,
                               linewidths=0)
                    n_arrows += 1
                else:
                    ax.scatter([p[0]], [p[1]], s=(lw * 0.7) ** 2, color=INK, zorder=4, linewidths=0)
            continue
        pts_ = [p for p in (on_line(l, base) for l in edges[i]["lines"]) if p is not None] + term_cut.get(sid, [])
        if not pts_:
            continue
        side = [float(np.dot(p - base, nrm)) for p in pts_]
        lo, hi = min(side) - SPACING * 0.55, max(side) + SPACING * 0.55
        if not any(n_ == base_name(st_["name"]) and np.hypot(*(b_ - base)) < 150 for n_, b_ in capsuled):
            capsuled.add((base_name(st_["name"]), tuple(base)))
            capsule(base + nrm * lo, base + nrm * hi)
        zone = LineString([base + nrm * lo, base + nrm * hi]).buffer(cap_r + 5)
        main = st_["term"]                           # real termini only (B only visits Dworzec PKP, in and out)
        bullets = [(l, "o") for l in sorted(main)] + [(l, "D") for l in sorted(var_end.get(sid, set()) - main)]
        requests.append((0 if key_st else 1, -len(st_["lines"]), sid, base, zone, bullets))

    # label placement: as close to the stop as possible (16 directions, sliding outwards from the capsule),
    # off other names at any cost, off the lines where that costs little distance
    occ = [g.buffer(lw * m_per_pt / 2 + 4) for g in drawn]
    occ_tree = STRtree(occ)
    frame = box(x0 + 60, y0 + 60, x1 - 60, y1 - 60)  # names stay on the map (stops outside it get none)
    placed = []
    labelled = set()
    for prio, _, sid, base, zone, bullets in sorted(requests, key=lambda r: (r[0], r[1])):
        name = base_name(stations[sid]["name"])
        if name in labelled:
            continue
        labelled.add(name)
        size = 14 if prio == 0 else 12.5
        w_m = len(name) * size * 0.52 * m_per_pt
        h_m = size * 1.25 * m_per_pt
        if bullets:
            w_m = max(w_m, len(bullets) * 15 * m_per_pt)
            h_m += 17 * m_per_pt
        best = None
        for ang in range(0, 360, 22):
            v = np.array([math.cos(math.radians(ang)), math.sin(math.radians(ang))])
            first = None
            for dist in np.arange(0, 200, 5.0):
                if first is not None and dist > first + 40:
                    break                                  # at most 40 m past the closest clear spot
                anchor = base + v * dist
                bx0 = anchor[0] if v[0] > 0.3 else anchor[0] - w_m if v[0] < -0.3 else anchor[0] - w_m / 2
                top = anchor[1] + (h_m if v[1] > 0.3 else h_m / 2 if abs(v[1]) <= 0.3 else 0)
                bb = box(bx0, top - h_m, bx0 + w_m, top)
                if bb.intersects(zone) or not frame.contains(bb):
                    continue
                first = dist if first is None else first
                hit = sum(occ[k].intersection(bb).area for k in occ_tree.query(bb))
                hit += sum(q.buffer(10).intersection(bb).area for q in placed) * 100
                score = hit + dist * 25
                if best is None or score < best[0]:
                    best = (score, bb, top, "left" if v[0] > 0.3 else "right" if v[0] < -0.3 else "center")
        if best is None:
            continue
        _, bb, top, ha = best
        placed.append(bb)
        ty = top - size * 0.62 * m_per_pt
        # text aligned towards the stop: a name left of its stop ends at the stop's side of the box
        tx = {"left": bb.bounds[0], "right": bb.bounds[2], "center": (bb.bounds[0] + bb.bounds[2]) / 2}[ha]
        ax.text(tx, ty, name, fontproperties=font(size), color=INK, ha=ha, va="center", path_effects=halo, zorder=7)
        row_w = len(bullets) * 15 * m_per_pt
        bx0 = {"left": tx, "right": tx - row_w, "center": tx - row_w / 2}[ha]
        for k, (l, mk) in enumerate(bullets):  # lines ending here: bullets under the name (diamond: variant)
            bx = bx0 + (k * 15 + 7) * m_per_pt
            by = ty - (size * 0.62 + 9) * m_per_pt
            ax.scatter([bx], [by], marker=mk, s=(13 if mk == "o" else 11.5) ** 2, color=colors.get(l, "#888"),
                       zorder=7, linewidths=0)
            ax.text(bx, by, l, fontproperties=font(8), color="white", ha="center", va="center", zorder=7.1)

    street_labels(ax, edges, orders, occ, occ_tree, placed, frame, m_per_pt)

    # title + legend panel
    px0 = x1 + panel * 0.06
    pw = panel * 0.88
    ax.add_patch(plt.Rectangle((px0, y1 - (y1 - y0) * 0.17), pw, (y1 - y0) * 0.14, color=INK, zorder=8))
    ax.text(px0 + pw * 0.05, y1 - (y1 - y0) * 0.075, "Biała Podlaska", fontproperties=font(46), color="white",
            va="center", zorder=9)
    d = pd.Timestamp(date)
    ax.text(px0 + pw * 0.05, y1 - (y1 - y0) * 0.135, f"Autobusy MZK · {DAY_NAMES[d.weekday()]} {d:%d.%m.%Y}",
            fontproperties=font(18), color="#d8d8dc", va="center", zorder=9)
    names = {g[0]: n for g, n in zip(line_cfg["groups"], line_cfg["names"])}
    row = (y1 - y0) * 0.052
    top = y1 - (y1 - y0) * 0.22
    ax.add_patch(plt.Rectangle((px0, top - row * (len(order) - 0.4)), pw, row * len(order), facecolor="white",
                               edgecolor=INK, linewidth=1.2, zorder=8))
    for k, l in enumerate(sorted(order)):
        y = top - k * row
        ax.scatter([px0 + pw * 0.06], [y], s=24 ** 2, color=colors.get(l, "#888"), zorder=9, linewidths=0)
        ax.text(px0 + pw * 0.06, y, l, fontproperties=font(14), color="white", ha="center", va="center", zorder=10)
        ax.text(px0 + pw * 0.13, y, names.get(l, l).split("  ", 1)[-1], fontproperties=font(15), color=INK,
                va="center", zorder=9)
    # key: stop, one-way stop, extension
    ky = top - row * (len(order) - 0.4) - row * 0.9
    kx0, kx1 = px0 + pw * 0.02, px0 + pw * 0.10
    key_rows = [("dot", "przystanek"), ("arrow", "przystanek w jednym kierunku jazdy"),
                ("dash", "wydłużenia: mniej niż 1/4 kursów linii")]
    for k, (kind, text) in enumerate(key_rows):
        y = ky - k * row * 0.75
        draw_line(ax, np.array([[kx0, y], [kx1, y]]), "#8a8a8a", lw, kind == "dash", zorder=9)
        if kind == "dot":
            ax.scatter([(kx0 + kx1) / 2], [y], s=(lw * 0.7) ** 2, color=INK, zorder=9.5, linewidths=0)
        elif kind == "arrow":
            ax.scatter([(kx0 + kx1) / 2], [y], marker=(3, 0, -90), s=(lw * 1.35) ** 2, color=INK, zorder=9.5,
                       linewidths=0)
        ax.text(px0 + pw * 0.13, y, text, fontproperties=font(13), color=INK, va="center", zorder=9)
    print(f"one-way stop arrows: {n_arrows}")
    ax.text(x0 + (x1 - x0) * 0.01, y0 + (y1 - y0) * 0.012,
            "Schemat poglądowy · przebiegi linii wg rozkładu jazdy MZK Biała Podlaska na środę 18.11.2026 · "
            "podkład © OpenStreetMap contributors · opracowanie: Jacek Gęborys, 2026", fontproperties=font(10),
            color="#555", zorder=9)
    fig.savefig(out, dpi=150, facecolor=LAND)
    pdf = Path(out).with_suffix(".pdf")             # vector copy for print
    fig.savefig(pdf, facecolor=LAND)
    plt.close(fig)
    print(f"written {out} and {pdf}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="20261118")
    ap.add_argument("--out", default=str(ROOT / "png" / "bus_network_city.png"))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ext-style", choices=["dashed", "thin", "pale", "hollow"], default=None)
    args = ap.parse_args()
    global EXT_STYLE
    EXT_STYLE = args.ext_style or EXT_STYLE
    geoms, n_shape, line_of, stations, line_trips = load(args.date)
    pieces = corridors(geoms, n_shape, line_of, line_trips)
    sets = lines_on(pieces, geoms, line_of)
    edges = straighten_roundabouts(build_edges(pieces, sets),
                                   [s["xy"] for s in stations.values() if s["name"] in KEY_STATIONS])
    edges = fold_parallel(contract_short(repair_lines(edges, geoms, line_of)))
    edges = drop_turning_stubs(edges, stations)
    edges = station_access(edge_trips(edges, geoms, n_shape, line_of), stations)
    print(f"corridor pieces {len(pieces)}, edges {len(edges)}")
    shapes = defaultdict(list)                   # line -> [(shape, trips of the day)]
    for sid, g in geoms.items():
        shapes[line_of[sid]].append((g, int(n_shape[sid])))
    set_joins(edges, shapes)
    orders = best_orders(sorted(line_trips.index), edges, args.seed)
    arm_shifts(edges, orders)
    relax_junction_corners(edges, orders)
    render(edges, orders, stations, json.load(open(CONFIG, encoding="utf-8")), args.date, args.out, shapes, line_trips)


if __name__ == "__main__":
    main()
