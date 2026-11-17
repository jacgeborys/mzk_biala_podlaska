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
import copy
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
EXACT_TURNS = True    # a lane that would fold at a tight junction corner is redrawn lane to lane (exact_turns)
LANE_CHANGE_SIDE = False  # lane changes at a junction on the side with fewer lines (off: centred; it pushed G under E on al. Tysiąclecia)
SPUR_RAMP_M = 300   # lane changes onto / off the Dworzec PKP spur (--station spur) ease over this length
FIX = {"e_bridge", "l_bridge"}   # preview switches for the overlap fixes (user judges them): spur_east, e_bridge, l_bridge, diamonds
RAMP_M = 120        # a line changing lane at a junction eases over this length (no jag)
EXT_TOL = 45        # m: shapes this close to a drawn point pass it (extension test)
EXT_MIN_M = 150     # m: shorter solid / dashed runs are merged into their neighbours
DASH_ON, DASH_OFF = 4.0, 0.55   # extensions: dash / gap length in line widths
EXT_SHARE = 0.5     # where a line runs with under this share of its trips of the day it is an extension (pale)
# Dworzec PKP: the U hangs below the Sidorska / Al. 1000-lecia crossing as deep as its parts need
# (U_TOP_M + S-bend + U_HALF_M, 582 m; was 690 with a 285-300 m straight top)
U_HALF_M = 130       # ... its legs this far either side of Al. 1000-lecia (also the radius of its bottom)
STATION_STYLE = "spur"  # Dworzec PKP: u / untangle / under / spur / box (--station; user chose spur, 2026-11-17)
STATION_CROSS_W = 0.3   # ... crossings / lane shifts at the station (under its bar) cost this much, so they go there
ABOVE_CROSS_W = 3.0     # ... and at the crossing above it this much (as few crossings above the station as possible)
U_TOP_M = 260        # ... the legs first run straight down side by side this far (arms turn into them on full arcs)
U_ARM_ALL_M = 700    # ... arms shorter than this are shifted along their whole length (straight: Sidorska west, 604 m)
U_S_R = 175          # ... then move apart to their place on an S-bend of two arcs of this radius (room for every lane)
CONTRACT_M = 40     # junctions joined by shorter edges are merged into one
PATH_FOLD_EDGES, PATH_FOLD_M = 4, 80   # an edge beside a path of this many edges, this close, is the same link
PRECISE_LINES = frozenset("G")   # lines drawn along their real streets where they run alone (user, 2026-11-17)
PRECISE_SIMPLIFY_M, PRECISE_R = 10, 90   # ... simplified this little, corners rounded this tight
CORNER_EASE = 0.0   # a far street corner next to a too-tight junction turn is pulled this far towards its chord
TURN_STUB_M = 120   # dead-end stubs this short with their stop by the junction are turning manoeuvres (dropped)
TRIANGLE_M = 80     # ... and by shorter edges closing a triangle (a split rejoining right away)
END_SPLIT_COST = 0.5  # order search: a line passing between lines that end at the junction (see node_cost)
SHIFT_COST = 0.15   # order search: cost per lane a line shifts through a junction (vs 1 per crossing)
KEY_STATIONS = {"Dworzec PKP"}   # always drawn as interchange with a name
MIN_INTERCHANGE = 5               # stations served by this many lines get a capsule and a name
MAIN_TERMINUS = 0.25              # a station is a line's terminus when this share of its trips start or end there
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


SHORT_WORKINGS = set()   # (line, end): short turns at Dworzec PKP on the line's main route, not in the legend
EXT_SIDE = {}       # (line, variant end) -> the main terminus it lies beyond
LOOP_STREETS = {"L": "Francuska – Armii Krajowej – Glinki – Narutowicza – Jana Pawła II – Francuska"}
LINE_ENDS = {}      # line -> {terminus name: trip starts + ends}, for the legend
END_STATIONS = defaultdict(dict)   # (line, terminus name) -> {station: trip starts + ends there}
EXT_MIN_ENDS = 3    # legend: variant ends with fewer trip starts + ends are not listed


def clean_stop(n):
    """Stop name for the legend: no post number, no 'Pętla', no notes in brackets."""
    return re.sub(r"\s+", " ", re.sub(r"\bPętla\b|\(.*?\)|\s\d+\b", "", n)).strip(" -")


def legend_name(l, fallback):
    """(left terminus, right terminus, variant ends beyond the left one, beyond the right one); a circular line
    (one terminus): ("Linia okrężna", None, its streets, "")."""
    ends = sorted(LINE_ENDS.get(l, {}).items(), key=lambda kv: -kv[1])
    if len(ends) < 2:
        return "Linia okrężna", None, LOOP_STREETS.get(l, fallback.split("  ", 1)[-1]), ""
    t1, t2 = ends[0][0], ends[1][0]
    ext = [n for n, c in ends[2:] if c >= EXT_MIN_ENDS and (l, n) not in SHORT_WORKINGS]
    left = [n for n in ext if EXT_SIDE.get((l, n)) == t1]
    right = [n for n in ext if n not in left]
    if len(left) > len(right):                      # the side with more variant ends on the right (more room)
        t1, t2, left, right = t2, t1, right, left
    return t1, t2, ", ".join(left), ", ".join(right)


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
    calls = st.drop_duplicates(["trip_id", "station"]).groupby(["station", "route_id"]).size()   # trips per line
    stations = {}
    for sid_, name, x, y in zip(st_tab.index, st_tab.stop_name, sx, sy):
        lines = set(st[st.station == sid_].route_id)
        if not lines:
            continue
        term = {l for l in lines if ends.get((sid_, l), 0) >= MAIN_TERMINUS * line_trips[l]}
        stations[sid_] = dict(name=name, xy=np.array([x, y]), lines=lines, term=term,
                              trips={l: int(calls.get((sid_, l), 0)) for l in lines})
    first_last["nm"] = first_last.stop_id.map(stops.set_index("stop_id").stop_name).map(clean_stop)
    for l, g in first_last.groupby("route_id"):
        LINE_ENDS[l] = g.nm.value_counts().to_dict()        # trip starts + ends per terminus name
        for (nm, sta), c in g.groupby(["nm", "station"]).size().items():
            END_STATIONS[(l, nm)][sta] = int(c)
    # which main terminus each other variant end lies beyond: the main terminus its trips pass (nearest to it in
    # the stop sequence); a trip passing neither (a short working) goes on the side away from its other end
    # (full stop names here: once cleaned, "Sidorska 5" and the terminus "Sidorska Pętla" would be one)
    st["raw"] = st.stop_id.map(stops.set_index("stop_id").stop_name)
    first_last["raw"] = first_last.stop_id.map(stops.set_index("stop_id").stop_name)
    seqs = st.sort_values("seq").groupby("trip_id").raw.apply(list)
    trip_line = day.set_index("trip_id").route_id
    lx, ly = tr.transform(stops.stop_lon.astype(float).values, stops.stop_lat.astype(float).values)
    stop_xy = {}
    for n_, x_, y_ in zip(stops.stop_name, lx, ly):
        stop_xy.setdefault(n_, np.array([x_, y_]))
    main_stops = defaultdict(set)                         # stops of the trips between a line's two main termini
    for tid, seq in seqs.items():
        l = trip_line[tid]
        main = [n for n, _ in sorted(LINE_ENDS[l].items(), key=lambda kv: -kv[1])[:2]]
        if {clean_stop(seq[0]), clean_stop(seq[-1])} == set(main):
            main_stops[l] |= set(seq[1:-1])
    votes = defaultdict(lambda: defaultdict(int))
    for tid, seq in seqs.items():
        l = trip_line[tid]
        main = [n for n, _ in sorted(LINE_ENDS[l].items(), key=lambda kv: -kv[1])[:2]]
        if len(main) < 2:
            continue
        fl_l = first_last[first_last.route_id == l]
        main_raw = {m: set(fl_l[fl_l.nm == m].raw) for m in main}
        for at_start in (True, False):
            x, other = (seq[0], seq[-1]) if at_start else (seq[-1], seq[0])
            if clean_stop(x) in main:
                continue
            pos, opos = (0, len(seq) - 1) if at_start else (len(seq) - 1, 0)
            passed = [(abs(i - pos), m) for i, n in enumerate(seq) if i not in (pos, opos)
                      for m in main if n in main_raw[m]]
            if x in main_stops.get(l, set()) and x in KEY_STATIONS:
                SHORT_WORKINGS.add((l, clean_stop(x)))         # a short turn at Dworzec PKP: not an extension
                continue
            if passed:
                side = min(passed)[1]
            else:                                              # else the main terminus nearest to it
                side = min(main, key=lambda m: min(np.hypot(*(stop_xy[r] - stop_xy[x])) for r in main_raw[m]))
            votes[(l, clean_stop(x))][side] += 1
    for (l, x), v in votes.items():
        EXT_SIDE[(l, x)] = max(v, key=v.get)
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


HAIRPINS = {}      # line -> stop at the tip of its hairpin next to a junction (drawn tight beside the bundle); G
#                   was here until 2026-11-17, now follows its real route
HAIRPIN_ALONG_M = 330  # ... it runs this far beside the bundle it comes back to


def hairpin_tips(edges, stations):
    """A single line that runs past a junction on a parallel street and comes back into it in a hairpin, with a stop
    on the way back (G at Narutowicza: north up the street east of al. Tysiąclecia, back WSW through the stop into
    the Narutowicza junction): the hairpin is drawn tight, its tip at the stop, so the stop reads as a short walk
    from the junction. Geometry fixed here (raw), not simplified or eased later."""
    for l, name in HAIRPINS.items():
        sts = [s_ for s_ in stations.values() if s_["name"] == name and l in s_["lines"]]
        es = [k for k, e in enumerate(edges) if e["lines"] == frozenset(l) and not e.get("raw")]
        if not sts or not es:
            continue
        S = Point(np.mean([s_["xy"] for s_ in sts], axis=0))
        k = min(es, key=lambda k: edges[k]["geom"].distance(S))
        e = edges[k]
        if e["geom"].distance(S) > 60:
            continue
        c = section_coords(e, 1)
        node = min((e["u"], e["v"]), key=lambda n: S.distance(Point(n)))
        if np.hypot(*(c[0] - np.array(node))) < np.hypot(*(c[-1] - np.array(node))):
            c = c[::-1]                                       # ... ending at the junction by the stop
        turn = []
        for j in range(1, len(c) - 1):
            v1, v2 = c[j] - c[j - 1], c[j + 1] - c[j]
            turn.append(math.acos(float(np.clip(v1 @ v2 / (np.hypot(*v1) * np.hypot(*v2) + 1e-9), -1, 1))))
        j = int(np.argmax(turn)) + 1                          # the hairpin vertex
        if turn[j - 1] < math.radians(90):
            continue
        # the bundle it comes back to (al. Tysiąclecia): G runs beside it as the next lane out, touching it, up to
        # a tip by the junction where the stop is, then into the junction
        nd = np.array(node)
        others = [f for q, f in enumerate(edges) if q != k and node in (f["u"], f["v"]) and len(f["lines"]) >= 3]
        if not others:
            continue
        bnd = min(others, key=lambda f: f["geom"].distance(Point(c[j - 1])))
        bc = section_coords(bnd, len(bnd["lines"]))
        if np.hypot(*(bc[0] - nd)) > np.hypot(*(bc[-1] - nd)):
            bc = bc[::-1]                                     # from the junction outwards
        bl = LineString(bc)
        off = (len(bnd["lines"]) - 1) / 2 * SPACING + SPACING
        def along(d):
            p_ = np.array(bl.interpolate(d).coords[0])
            q_ = np.array(bl.interpolate(d + 5).coords[0])
            tg = (q_ - p_) / max(np.hypot(*(q_ - p_)), 1e-9)
            return p_, np.array([-tg[1], tg[0]])
        p_, n_ = along(HAIRPIN_ALONG_M)
        side = 1.0 if (c[j - 1] - p_) @ n_ > 0 else -1.0
        d0 = off + SPACING                                    # the tip: just clear of the junction's turn
        run = [along(d)[0] + along(d)[1] * side * off for d in np.linspace(HAIRPIN_ALONG_M, d0, 8)]
        T = run[-1]
        new = [tuple(p) for p in c[:j - 1]] + [tuple(q) for q in run] + [tuple(c[-1])]
        g = LineString(new)
        edges[k] = dict(e, geom=g if tuple(np.round(new[0], 1)) == e["u"] or np.hypot(*(np.array(new[0]) - e["u"])) < 1
                        else LineString(new[::-1]), raw=True, schematic=S,
                        hug=(next(q for q, f in enumerate(edges) if f is bnd), side, node, len(c[:j - 1])))
        for s_ in sts:
            s_["xy"] = T                                      # the stop dot at the tip
        print(f"hairpin {l} at {name}: tip {np.hypot(*(T - c[-1])):.0f} m from the junction")


def drop_turning_stubs(edges, stations):
    """Short dead-end stubs (< TURN_STUB_M) whose stop lies nearer their junction than their dead end are the
    terminal turning manoeuvre, not a route: dropped, so the lines end on the through road at the junction and are
    cut flush at the station capsule there (Terebelska Pętla: a 55 m stub across the end of the trunk, stop 23 m
    from the junction). A spur with its stop out at the end (Cmentarz Komunalny, 101 m) stays."""
    deg = defaultdict(int)
    for e in edges:
        deg[e["u"]] += 1
        deg[e["v"]] += 1
    pts = [Point(s["xy"]) for s in stations.values()]
    out, n = [], 0
    for e in edges:
        if e["geom"].length < TURN_STUB_M and (deg[e["u"]] == 1) != (deg[e["v"]] == 1):
            dead, junc = (e["u"], e["v"]) if deg[e["u"]] == 1 else (e["v"], e["u"])
            near = [p for p in pts if p.distance(e["geom"]) < 60]
            if near and all(p.distance(Point(junc)) < p.distance(Point(dead)) for p in near):
                n += 1
                continue
        out.append(e)
    print(f"terminal turning stubs dropped: {n}")
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
        # a line in and out along one arm from a junction where it comes from one street and leaves along another
        # (B: north from Narutowicza, south to Łomaska, Dworzec PKP a detour along Sidorska) is drawn twice along
        # that arm, its second pass a lane outside the bundle, and goes round the U: in on one leg, out on the
        # other (user, 2026-11-17: B from the north on the upper lane, into the U from the east, out to the west)
        twins = {}                                                    # line -> (arm, its far node, the two others)
        for l in sorted(end_lines):
            arm = next(f for f in at(N2) if l in f["lines"])
            W = far(arm, N2)
            ws = [f for f in at(W) if l in f["lines"] and f is not arm]
            if len(ws) == 2:
                twins[l] = (arm, W, ws)
        end_lines = end_lines - frozenset(twins)
        u_lines = L - kos_lines - end_lines                           # round the U and back up
        trips = e.get("trips", {})
        tr = lambda ls: {l: trips.get(l, 0) for l in ls}  # noqa: E731
        # the U, drawn as is (no simplification): legs leave the crossing diagonally, run straight down, and
        # meet in a half circle whose lowest point is the station
        phi = math.acos(1 - U_HALF_M / (2 * U_S_R))                   # S-bend: two arcs, U_HALF_M sideways
        s_len = 2 * U_S_R * math.sin(phi)                             # ... over this length
        depth = U_TOP_M + s_len + U_HALF_M          # no straight below the S-bend: the station as high as it fits
        C = P0 + dv * (depth - U_HALF_M)
        Sx = P0 + dv * depth
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
        if STATION_STYLE == "spur":
            # a short dead-end spur, as in reality: both legs straight down side by side (west: all lines in, the
            # lines going on along Kościuszki outermost; east: the lines coming back up), each line turning back on
            # a nested hairpin round one point on the boundary between the legs (innermost lane: half a lane
            # radius); the lines to Kościuszki go round the outside of the hairpins and on east. Centrelines run
            # RHO either side of the boundary (a jog of RHO at the top), the lanes put in place by shift_prof.
            RHO = SPACING / 2
            GAP = SPACING / 2                      # innermost hairpin: radius one lane (half a lane looked filled)
            O_ = P0 + dv * U_TOP_M                                       # hairpin centre, on the boundary
            Sx = O_ + dv * RHO
            C = O_
            jog = 60.0

            def spur_leg(sign):
                # one straight from the crossing to the hairpin (RHO sideways over its length: no jog)
                top = [P0]
                straight = [O_ + sign * perp * RHO]
                arc = [O_ + sign * perp * RHO * math.cos(a) + dv * RHO * math.sin(a)
                       for a in np.linspace(0, math.pi / 2, 7)[1:]]
                return [tuple(p) for p in top + straight + arc]
            lg = spur_leg(-1)
            sh_top_L = -nL * SPACING / 2                  # left of travel (down: east +): lanes west of the boundary
            lprof = [sh_top_L + RHO - GAP] * len(lg)
            keep.append(dict(geom=LineString([N2] + lg[1:-1] + [tuple(Sx)]), lines=down, u=N2, v=key(Sx), raw=True,
                             twin_side=-1, trips=tr(down), shift_prof=lprof, schematic=Point(xy), down_dir=dv))
            if u_lines:
                rg = spur_leg(1)
                sh_top_R = -nR * SPACING / 2              # left of travel (up: west +): lanes east of the boundary
                rprof = [sh_top_R + RHO - GAP] * len(rg)
                keep.append(dict(geom=LineString([tuple(Sx)] + rg[::-1][1:-1] + [N2]), lines=u_lines, u=key(Sx),
                                 v=N2, raw=True, twin_side=-1, trips=tr(u_lines), shift_prof=rprof[::-1],
                                 schematic=Point(xy), down_dir=dv))
            if kos_lines:
                kc = int(np.argmax([LineString([thr_c[0], thr_c[-1]]).distance(Point(p)) for p in thr_c]))
                corner = thr_c[kc]
                Lk = max(float(np.dot(corner - Sx, perp)), 520.0)
                turn = Sx + perp * Lk
                # outside the hairpins: the outermost (u_lines) lanes reach (nL - 0.5) lanes below the centre; past
                # the hairpins (RAMP_A on) they rise gently to the street, before the corner north
                ksh = -(len(u_lines) + len(kos_lines) / 2) * SPACING + RHO - GAP
                ra, rb = (len(u_lines) + 2) * SPACING, Lk - 220
                ramp = [Sx + perp * (ra + (rb - ra) * f) for f in np.linspace(0, 1, 9)]
                rest = list(LineString([tuple(turn)] + [tuple(p) for p in thr_c[kc:0:-1]] + [N1])
                            .simplify(SIMPLIFY_M).coords)
                ss_ = lambda t: t * t * (3 - 2 * t)  # noqa: E731
                # D / F / H keep hugging the spur (user, 2026-11-17): round its bottom on their lanes outside the
                # hairpins, north along the outside of the east leg, then a right turn onto Kościuszki at the
                # height of its real corner, east to that corner and on along the real street (unchanged). P is the
                # path of their middle lane; the centreline is P moved by the lane shift, which fades to 0 along
                # the horizontal so they arrive on the street's centreline
                north = -dv
                Rb = RHO - ksh                                      # radius of their middle lane round the hairpins
                E1 = O_ + perp * Rb
                h_up = float(np.dot(corner - E1, north))
                Rt = float(np.clip(0.75 * h_up, 60.0, 140.0))
                lift = max(h_up - Rt, 0.0)
                P_ = [O_ + Rb * (dv * math.cos(a) + perp * math.sin(a)) for a in np.linspace(0, math.pi / 2, 24)]
                P_ += [E1 + north * t for t in np.linspace(0, lift, max(2, int(lift / 8)))[1:]]
                T0 = E1 + north * lift
                cR = T0 + perp * Rt
                P_ += [cR - perp * Rt * math.cos(a) + north * Rt * math.sin(a) for a in np.linspace(0, math.pi / 2, 14)[1:]]
                H0 = P_[-1]
                corner = corner + perp * float(np.dot(np.array(N1) - corner, perp))   # straight under the junction
                hl = float(np.dot(corner - H0, perp))
                n_h = max(2, int(hl / 8))
                # the horizontal turns onto the real street before its corner, on an arc of radius Rc
                sdir = np.array(N1) - corner                    # straight on to the junction after the turn
                sdir = sdir / max(np.hypot(*sdir), 1e-9)
                th_c = math.acos(float(np.clip(perp @ sdir, -1, 1)))
                Rc = 70.0
                tl_c = Rc * math.tan(th_c / 2)
                hl -= tl_c
                n_h = max(2, int(hl / 8))
                P_ += [H0 + perp * t for t in np.linspace(0, hl, n_h)[1:]]
                A0 = P_[-1]
                lc = np.array([-perp[1], perp[0]])
                if lc @ sdir < 0:
                    lc = -lc
                cc = A0 + lc * Rc
                a_s = math.atan2(*(A0 - cc)[::-1])
                sg_ = 1 if (perp[0] * sdir[1] - perp[1] * sdir[0]) > 0 else -1
                arc_c = [cc + Rc * np.array([math.cos(a_s + sg_ * th_c * f), math.sin(a_s + sg_ * th_c * f)])
                         for f in np.linspace(0, 1, 10)[1:]]
                P_ += arc_c
                P_ = np.array(P_)
                # the shift fades gently along the whole way from halfway round the bottom to the corner (a short
                # fade bends the lanes: the shift is ~5 lanes), so the lanes run on P and reach the street centred
                sP = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(P_, axis=0).T))])
                s0 = sP[12]                                  # halfway round the bottom ... to the last turn
                s1 = sP[len(P_) - len(arc_c) - 1]
                fz_ = np.clip((sP - s0) / max(s1 - s0, 1.0), 0, 1)
                w_ = ksh * (1 - np.array([ss_(float(f)) for f in fz_]))
                tg = np.gradient(P_, axis=0)
                tg /= np.maximum(np.hypot(*tg.T), 1e-9)[:, None]
                C_ = P_ - np.stack([-tg[:, 1], tg[:, 0]], axis=1) * w_[:, None]
                C_[0] = Sx
                rest = [tuple(P_[-1]), tuple(N1)]
                kg = [tuple(p) for p in C_] + rest[1:]
                kprof = list(w_) + [0.0] * (len(rest) - 1)
                keep.append(dict(geom=LineString(kg), lines=kos_lines, u=key(Sx), v=N1, trips=tr(kos_lines), raw=True,
                                 probe_m=6.0,      # its heading at the station: east, as it leaves (not north)
                                 shift_prof=kprof,
                                 schematic=LineString(thr_c).interpolate(0.5, normalized=True)))
                kos_g = keep[-1]["geom"]
                for st in stations.values():
                    if st["name"] != name and Point(st["xy"]).distance(LineString(thr_c)) < 40:
                        st["xy"] = np.array(kos_g.interpolate(kos_g.project(Point(st["xy"]))).coords[0])
            for k_ in sts:
                stations[k_]["xy"] = Sx + dv * 1.0
            U_INFO.update(spur=True)
        else:
            Sx = None
        # shift: metres left of the edge's direction (left leg: down, so west is negative; right leg: up, so east
        # is negative)
        lg, lprof = leg(-1, -nR * SPACING / 2) if Sx is None else (None, None)
        probe = U_TOP_M + s_len                                   # headings measured where the legs are apart
        Sx = P0 + dv * depth if Sx is None else Sx
        spur = STATION_STYLE == "spur"
        if not spur:
            keep.append(dict(geom=LineString([N2] + lg[1:-1] + [tuple(Sx)]), lines=down, u=N2, v=key(Sx), raw=True,
                             probe_m=probe, trips=tr(down), shift_prof=lprof, schematic=Point(xy)))
        if u_lines and not spur:
            rg, rprof = leg(1, -nL * SPACING / 2)
            keep.append(dict(geom=LineString([tuple(Sx)] + rg[::-1][1:-1] + [N2]), lines=u_lines, u=key(Sx), v=N2,
                             raw=True, probe_m=probe, trips=tr(u_lines), shift_prof=rprof[::-1], schematic=Point(xy)))
        if kos_lines and not spur:
            kc = int(np.argmax([LineString([thr_c[0], thr_c[-1]]).distance(Point(p)) for p in thr_c]))
            corner = thr_c[kc]
            turn = Sx + perp * max(float(np.dot(corner - Sx, perp)), U_HALF_M + 200)  # east along the bottom, then up
            kg = [tuple(Sx), tuple(turn)] + [tuple(p) for p in thr_c[kc:0:-1]] + [N1]
            kd = dict(lines=kos_lines, u=key(Sx), v=N1, trips=tr(kos_lines),
                      schematic=LineString(thr_c).interpolate(0.5, normalized=True))
            if STATION_STYLE in ("untangle", "under"):
                # no crossings at the station: these lines are the outermost lanes of the west leg's bottom, so
                # they leave on those lanes (shifted south by the U lines' half width), drawn as is
                kg = [tuple(Sx)] + list(LineString(kg[1:]).simplify(SIMPLIFY_M).coords)
                ksh = -(len(down) - len(kos_lines)) / 2 * SPACING
                kd.update(raw=True, shift_prof=[ksh, ksh] + [0.0] * (len(kg) - 2))
            keep.append(dict(geom=LineString(kg), **kd))
            kos_g = keep[-1]["geom"]
            for st in stations.values():                              # stops on Kościuszki follow it
                if st["name"] != name and Point(st["xy"]).distance(LineString(thr_c)) < 40:
                    st["xy"] = np.array(kos_g.interpolate(kos_g.project(Point(st["xy"]))).coords[0])
        if not spur:
            for k_ in sts:
                stations[k_]["xy"] = Sx - perp * 1.0                 # on the leg every line uses
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
        for l, (arm, W, ws) in twins.items():
            b_k = idx[id(arm)]
            g_ = arm["geom"]
            s_ = 1 if arm["u"] == W else -1                           # the twin: left of travel from W to N2
            tw = dict(arm, lines=frozenset({l}), trips={l: arm.get("trips", {}).get(l, 0)}, twin_of=b_k,
                      twin_side=s_)
            out.append(tw)
            t_k = len(out) - 1
            LANE_SHIFT[b_k] = -s_ * SPACING / 2                       # the bundle half a lane over, the twin outside
            LANE_SHIFT[t_k] = s_ * len(arm["lines"]) * SPACING / 2
            # at W: the section coming in on the left of travel towards N2 joins the twin, the other the arm
            d_ = np.array(g_.interpolate(30 if arm["u"] == W else g_.length - 30).coords[0]) - np.array(W)
            def leftness(f):
                h = f["geom"]
                q = np.array(h.interpolate(30 if f["u"] == W else h.length - 30).coords[0]) - np.array(W)
                return d_[0] * q[1] - d_[1] * q[0]
            w_tw, w_arm = sorted(ws, key=leftness, reverse=True)
            PAIRS[(l, W)] = {t_k: idx[id(w_tw)], idx[id(w_tw)]: t_k, b_k: idx[id(w_arm)], idx[id(w_arm)]: b_k}
            # at N2: lanes keep their side through the turn down, so the twin takes the leg left of travel down
            nl = np.array([-dv[1], dv[0]])
            tw_leg, arm_leg = (right_k, left_k) if float(np.dot(perp, nl)) > 0 else (left_k, right_k)
            PAIRS[(l, N2)] = {t_k: tw_leg, tw_leg: t_k, b_k: arm_leg, arm_leg: b_k}
        for l in u_lines - frozenset(twins):
            arms = [k for k, f in enumerate(out) if l in f["lines"] and N2 in (f["u"], f["v"])
                    and k not in (left_k, right_k)]
            if len(arms) != 2 or right_k is None:
                continue
            w_arm, e_arm = sorted(arms, key=east)
            PAIRS[(l, N2)] = {left_k: w_arm, w_arm: left_k, right_k: e_arm, e_arm: right_k}
        for k, f in enumerate(out):                                   # arms: shifted where they join the legs
            if N2 in (f["u"], f["v"]) and k not in (left_k, right_k):
                if f["geom"].length >= U_ARM_ALL_M:
                    f["split_at"] = (N2, U_TOP_M)
                else:                 # too short to shift back before the next junction (A on Sidorska east,
                    f["shift_all"] = True       # 445 m): shifted all along, so A runs straight between two arcs
        U_INFO.update(node=N2, left=left_k, right=right_k, dv=dv, perp=perp, Sx=np.array(Sx), name=name)
        NODE_W[key(Sx)] = STATION_CROSS_W
        NODE_W[N2] = ABOVE_CROSS_W
        print(f"{name}: U {sorted(u_lines)}, on along Kościuszki {sorted(kos_lines)}, ending at the station "
              f"{sorted(end_lines)}, arm drawn twice {sorted(twins)}")
    return out


# ---------------------------------------------------------------- 4: line order (per edge)
def drop_duplicate_sections(edges):
    """A section between the same two junctions as another that carries all its lines and either runs within 45 m
    of it or is a longer way round of at most 600 m (F's one-way terminal loop at Akademicka, whose corners the
    relaxation flattened onto the K / F section beside it) is emptied (kept as an edge so indices stay valid): it
    was drawn as a second line on top of the bundle."""
    n = 0
    for i, e in enumerate(edges):
        if not e["lines"] or e.get("raw") or "twin_of" in e:
            continue                      # (a twin is B's second pass along Sidorska, meant to be there)
        for j, f in enumerate(edges):
            if j == i or not f["lines"] or f.get("raw") or {e["u"], e["v"]} != {f["u"], f["v"]}:
                continue
            if e["lines"] <= f["lines"] and (e["geom"].hausdorff_distance(f["geom"]) < 45 or
                                             f["geom"].length < e["geom"].length <= 600) and                     (e["lines"] < f["lines"] or e["geom"].length > f["geom"].length):
                e["lines"] = frozenset()
                n += 1
                break
    print(f"duplicate sections dropped: {n}")


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
        tw = 0.08 * e.get("twin_side", 0)                  # a twin lies beside its arm: just left of it (u -> v)
        nodes[e["u"]].append((i, heading_at(g, k) + tw, +1))
        nodes[e["v"]].append((i, heading_at(g, g.length - k) + math.pi - tw, -1))
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
    # lines ending here kept together: a line passing through with lines ending here on both sides of it on the
    # same section costs END_SPLIT_COST (H between I and M at Terebelska Pętla, its capsule cut across them)
    ends_here = {l for l, ea in per_line.items() if len(ea) == 1}
    if ends_here:
        for i, h, sgn in inc:
            ls = orders[i]
            for k, l in enumerate(ls):
                if l not in ends_here and any(x in ends_here for x in ls[:k]) and any(x in ends_here for x in ls[k + 1:]):
                    cost += END_SPLIT_COST
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
                # a line joining from a street of its own (a precise section, G) crossing there counts double:
                # it joins the bundle on its own side and the crossing moves along the bundle
                wx = 2 if any(is_precise(edges[x_]) for x_ in (qa, qb)) else 1
                if qa == qb:                             # on together: left/right must stay
                    if (la > lb) != (left_of(orders, qa, a, sa) > left_of(orders, qa, b, sb)):
                        crossings += wx
                else:                                    # split: the line turning further left keeps left
                    ta, tb = turn(entry, ha), turn(entry, hb)
                    if abs(ta - tb) > 0.05 and (la > lb) != (ta > tb):
                        crossings += wx
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
STATION_REORDER_M = 900   # --station other than u: sections this close to the station crossing are re-ordered


def refine_orders(edges, orders, only, seed, iters=20000):
    """Annealing of the line orders of the sections `only` (the rest kept), starting from `orders`."""
    nodes = node_passes(edges)
    orders = [o[:] for o in orders]
    ends = [(e["u"], e["v"]) for e in edges]
    cost_at = {n: node_cost(inc, orders, edges, n) for n, inc in nodes.items()}
    multi = [i for i in only if len(orders[i]) > 1]
    if not multi:
        return orders
    rnd = random.Random(seed + 11)
    cur = sum(cost_at.values())
    best, best_c = [o[:] for o in orders], cur
    for it in range(iters):
        i = rnd.choice(multi)
        a_, b_ = rnd.sample(range(len(orders[i])), 2)
        orders[i][a_], orders[i][b_] = orders[i][b_], orders[i][a_]
        touched = set(ends[i])
        new = {n: node_cost(nodes[n], orders, edges, n) for n in touched}
        delta = sum(new[n] - cost_at[n] for n in touched)
        T = 0.5 * (1 - it / iters) + 1e-3
        if delta <= 0 or rnd.random() < math.exp(-delta / T):
            cost_at.update(new)
            cur += delta
            if cur < best_c - 1e-9:
                best, best_c = [o[:] for o in orders], cur
        else:
            orders[i][a_], orders[i][b_] = orders[i][b_], orders[i][a_]
    print(f"line order near the station re-searched ({len(multi)} sections): cost {best_c:.1f}")
    return best


ORDER_NEXT_TO = [("E", "B")]   # (line, neighbour): in big bundles (8+ lines) the line runs right after it


def through_lines_east(edges, orders):
    """A terminus on a plain street (a node with two sections, lines ending there from both sides, some passing
    through: H at Al. Solidarności, where D / G end from the south and E / K from the north): the lines passing
    through run on the east side of both sections, straight through, the ending lines on the west side, so they
    read as ending there (user, 2026-11-17). Other lines keep their order."""
    at = defaultdict(list)
    for i, e in enumerate(edges):
        if e["lines"] and not e.get("raw"):
            at[e["u"]].append(i)
            at[e["v"]].append(i)
    n_ = 0
    for n, es in at.items():
        if len(es) != 2:
            continue
        a, b = es
        thru = edges[a]["lines"] & edges[b]["lines"]
        if not thru or len((edges[a]["lines"] | edges[b]["lines"]) - thru) < 2:
            continue
        if any(len(at.get(m, [])) == 1 for m in ()):
            continue
        for i in (a, b):
            g = edges[i]["geom"]
            p0, p1 = np.array(g.coords[0]), np.array(g.interpolate(min(60, g.length)).coords[0])
            v = p1 - p0
            left_east = -v[1] > 0                 # left normal (-vy, vx) points east
            rest = [l for l in orders[i] if l not in thru]
            th = [l for l in orders[i] if l in thru]
            orders[i] = th + rest if left_east else rest + th
        n_ += 1
    if n_:
        print(f"termini with through lines east: {n_}")


def lane(orders, i, l):
    ls = orders[i]
    return ((len(ls) - 1) / 2 - ls.index(l)) * SPACING + LANE_SHIFT.get(i, 0.0)


def vertex_R(edges, orders, adj):
    """Corner radius at a vertex between sections adj: PRECISE_R inside a precise line, else CORNER_R, but never
    tighter than the widest bundle there + 2 lanes (lanes and shifts of the U included), the same for every line
    through it (concentric lanes). Corners inside the U's legs keep CORNER_R (drawn as designed)."""
    if adj and all(is_precise(edges[k]) for k in adj):
        return PRECISE_R
    if adj and all(edges[k].get("raw") for k in adj):
        return CORNER_R
    return max([CORNER_R] + [edge_extent(edges, orders, k) + 2 * SPACING for k in adj])


TLCAP = {}      # (edge, segment of its section) -> metres an arc may use at its start / end


def compute_tl_caps(edges, orders):
    """For every segment of every drawn section: how much of it the arcs at its two ends may use. Each end asks
    for R * tan(turn / 2) (the largest over the ways lines turn there; at a junction every pair of sections sharing
    a line); when both do not fit, the segment is shared in proportion (0.95 of it), so a corner next to a short
    piece whose far end is straight can use nearly all of it (L into Zamkowa)."""
    TLCAP.clear()
    at = defaultdict(list)
    for i, e in enumerate(edges):
        if i in _SECTIONS and len(_SECTIONS[i]) >= 2:
            at[e["u"]].append(i)
            at[e["v"]].append(i)

    def ends(i, n):                       # unit direction of section i's end segment pointing away from node n
        c = _SECTIONS[i]
        if np.hypot(*(c[0] - np.array(n))) <= np.hypot(*(c[-1] - np.array(n))):
            v = c[1] - c[0]
        else:
            v = c[-2] - c[-1]
        return v / max(np.hypot(*v), 1e-9)
    node_need = {}
    for n, es in at.items():
        for a in es:
            best = 0.0
            for b in es:
                if a == b or not edges[a]["lines"] & edges[b]["lines"]:
                    continue
                th = math.pi - math.acos(float(np.clip(ends(a, n) @ ends(b, n), -1, 1)))
                if th < math.radians(1) or th > math.radians(179):
                    continue
                best = max(best, vertex_R(edges, orders, [a, b]) * math.tan(th / 2))
            node_need[(a, n)] = best
    for i, c in _SECTIONS.items():
        if i >= len(edges) or len(c) < 2:
            continue
        e = edges[i]
        need = np.zeros(len(c))
        for k in range(1, len(c) - 1):
            v1, v2 = c[k] - c[k - 1], c[k + 1] - c[k]
            th = math.acos(float(np.clip(v1 @ v2 / (np.hypot(*v1) * np.hypot(*v2) + 1e-9), -1, 1)))
            if th > math.radians(1):
                need[k] = vertex_R(edges, orders, [i]) * math.tan(min(th, math.radians(178)) / 2)
        for k_, idx in ((0, 0), (len(c) - 1, -1)):
            n = min((e["u"], e["v"]), key=lambda n_: np.hypot(*(np.array(n_) - c[idx])))
            need[k_] = node_need.get((i, n), 0.0)
        for k in range(len(c) - 1):
            L_ = float(np.hypot(*(c[k + 1] - c[k])))
            a, b = need[k], need[k + 1]
            if a + b <= 0.95 * L_ or a + b <= 0:
                TLCAP[(i, k)] = (0.95 * L_ - b, 0.95 * L_ - a)     # each end: whatever the other leaves
            else:
                TLCAP[(i, k)] = (0.95 * L_ * a / (a + b), 0.95 * L_ * b / (a + b))


def edge_extent(edges, orders, i):
    """How far the farthest lane of edge i lies from its centreline (lane shifts of an arm of the U included)."""
    e = edges[i]
    m = max((abs(lane(orders, i, l)) for l in e["lines"]), default=0.0)
    if e.get("shift_all") and "end_shift" in e:
        m += abs(e["end_shift"][1])
    if "shift_prof" in e:                       # a leg of the U: at the crossing, pushed outwards
        m += max(abs(x) for x in e["shift_prof"])
    if i in HOLD:
        m = max([m] + [abs(o) for o in HOLD[i][2].values()])
    return m


def is_precise(e):
    """A section drawn precisely along its streets: carries only PRECISE_LINES (G where it runs alone)."""
    return bool(e["lines"]) and e["lines"] <= PRECISE_LINES and not e.get("raw")


def fillet(c, R, cap=None):
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
        Ri = R[i] if np.ndim(R) else R               # per-vertex radius (precise lines: small)
        if cap is None:
            tl = min(Ri * math.tan(th / 2), seg[i - 1] * 0.5, seg[i] * 0.5)
        else:                                        # metres of each leg the arc may use (TLCAP)
            tl = min(Ri * math.tan(th / 2), cap[1][i - 1], cap[0][i])
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
    if is_precise(e):                             # along its real streets: only street wiggles under 5 m go
        return np.array(e["geom"].simplify(PRECISE_SIMPLIFY_M).coords)
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
LANE_SHIFT = {}  # edge -> metres its lanes are shifted left (an arm drawn twice for one line: B on Sidorska)
PAIRS = {}      # (line, node) -> {edge: edge}: how a line touching a node 4 times pairs its sections (the U)
HUG_TIPS = {}   # edge -> its hairpin tip as drawn (the stop dot goes there)
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
                    # sized for how far the lines turning from a onto b sit from the centre of either street;
                    # a vertex that is itself a sharp street corner is never dropped (sizing for the whole 9-line
                    # bundle used to erase G's corner onto al. Jana Pawła II and its S-bend)
                    shared = edges[a]["lines"] & edges[b]["lines"]
                    for k, c_, l_ in ((a, ca, l1), (b, cb[::-1], l2)):
                        r_k = max(max(abs(lane(orders, a, l)), abs(lane(orders, b, l))) for l in shared)                             + 0.8 * SPACING
                        need = 2 * r_k * math.tan(th / 2)
                        e_ = edges[k]
                        if l_ >= need or len(c_) < 3 or e_.get("raw") or is_precise(e_) or (e_.get("split_at") or (None,))[0] == n:
                            continue
                        v = len(c_) - 2                                    # the vertex next to n ...
                        w1, w2 = c_[v] - c_[v - 1], c_[v + 1] - c_[v]
                        cv = np.dot(w1, w2) / (np.hypot(*w1) * np.hypot(*w2) + 1e-9)
                        if cv < math.cos(math.radians(45)) and (np.hypot(*(c_[v] - c_[-1])) > 250 or len(e_["lines"]) <= 4):
                            # ... unless it is a real street corner away from the junction (G onto al. Jana Pawła II,
                            # or any street corner of a bundle of up to 4 lines: B, E, G, L north on al. Tysiąclecia
                            # then into Brzeska, so they do not cut through the neighbourhood;
                            # 456 m on): that one is only eased, pulled CORNER_EASE of the way towards the chord of its
                            # neighbours, once
                            if not e_.get("eased"):
                                c2 = c_.copy()
                                c2[v] = c_[v] + CORNER_EASE * ((c_[v - 1] + c_[v + 1]) / 2 - c_[v])
                                _SECTIONS[k] = c2 if np.array_equal(c_[0], _SECTIONS[k][0]) else c2[::-1]
                                e_["eased"] = True
                                changed = True
                            continue
                        c2 = np.delete(c_, v, axis=0)
                        _SECTIONS[k] = c2 if np.array_equal(c_[0], _SECTIONS[k][0]) else c2[::-1]
                        dropped += 1
                        changed = True
        if not changed:
            break
    # a single line leaving a bundle at a junction: its first sharp street corner within 500 m is pulled
    # CORNER_EASE of the way towards the chord of its neighbours, so it stays close to the bundle it left (G along
    # Narutowicza: its corner onto al. Jana Pawła II, 456 m on, drew it far away from al. Tysiąclecia)
    eased = 0
    for k, e in enumerate(edges):
        if len(e["lines"]) != 1 or e.get("raw") or k not in _SECTIONS or is_precise(e):
            continue
        if np.hypot(*(np.array(e["u"]) - np.array(e["v"]))) < 0.5 * e["geom"].length:
            continue                       # a loop out and back to the bundle (I via Os. Sitnickie): keep its shape
        for n in (e["u"], e["v"]):
            if len(at.get(n, [])) < 3 or max(len(edges[j]["lines"]) for j in at[n]) < 3:
                continue
            c_ = end(k, n)[::-1]                                   # starting at the junction
            if len(c_) < 3:
                continue
            w1, w2 = c_[1] - c_[0], c_[2] - c_[1]
            cv = np.dot(w1, w2) / (np.hypot(*w1) * np.hypot(*w2) + 1e-9)
            if cv < math.cos(math.radians(45)) and np.hypot(*w1) < 500:
                c2 = c_.copy()
                c2[1] = c_[1] + CORNER_EASE * ((c_[0] + c_[2]) / 2 - c_[1])
                c2 = c2[::-1]
                _SECTIONS[k] = c2 if np.hypot(*(c2[-1] - np.array(n))) <= np.hypot(*(c2[0] - np.array(n)))                     and np.array_equal(end(k, n)[-1], _SECTIONS[k][-1]) else c2[::-1]
                eased += 1
    for k, e in enumerate(edges):                   # a twin runs exactly beside its arm
        if "twin_of" in e:
            _SECTIONS[k] = _SECTIONS[e["twin_of"]]
    print(f"junction corners relaxed: {dropped} vertices dropped, {eased} single-line corners eased")


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


def rebuild_hugs(edges, orders):
    """The run of a hairpin line beside its bundle (hairpin_tips) redone from the bundle's final drawn centreline
    (after relax_junction_corners) and its real outermost lane, so it touches the bundle exactly one lane out."""
    for k, e in enumerate(edges):
        if "hug" not in e or k not in _SECTIONS:
            continue
        bi, side, node, n_keep = e["hug"]
        nd = np.array(node)
        bc = _SECTIONS[bi]
        if np.hypot(*(bc[0] - nd)) > np.hypot(*(bc[-1] - nd)):
            bc = bc[::-1]
        bl = LineString(bc)
        sgn = 1 if np.hypot(*(np.array(edges[bi]["geom"].coords[0]) - nd)) < 1 else -1   # lane sign from the node
        lanes = [lane(orders, bi, l) * sgn for l in edges[bi]["lines"]]
        off = (max(lanes) if side > 0 else -min(lanes)) + SPACING
        def along(d):
            p_ = np.array(bl.interpolate(d).coords[0])
            q_ = np.array(bl.interpolate(d + 5).coords[0])
            tg = (q_ - p_) / max(np.hypot(*(q_ - p_)), 1e-9)
            return p_ + np.array([-tg[1], tg[0]]) * side * off
        c = np.array(e["geom"].coords)
        if np.hypot(*(c[0] - nd)) < np.hypot(*(c[-1] - nd)):
            c = c[::-1]
        run = [along(d) for d in np.linspace(HAIRPIN_ALONG_M, off + SPACING, 8)]
        new = np.array(list(c[:n_keep]) + run + [nd])
        _SECTIONS[k] = new if np.hypot(*(new[0] - np.array(e["u"]))) < 1 else new[::-1]
        HUG_TIPS[k] = run[-1]


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
    g0 = edges[lk]["geom"]
    dvec = np.array(g0.coords[1]) - np.array(g0.coords[0])
    dvec /= np.hypot(*dvec)                     # down the legs
    for k, e in enumerate(edges):
        if "split_at" not in e and not e.get("shift_all"):
            continue
        g_ = e["geom"]
        q_ = np.array(g_.interpolate(g_.length - 60 if e["v"] == N2 else 60).coords[0])
        into = np.array(N2) - q_
        if np.dot(into / max(np.hypot(*into), 1e-9), dvec) < math.cos(math.radians(30)):
            # an arm turning into the legs (Sidorska west / east): not shifted at all, it stays on its street; its
            # lanes turn into theirs on the legs on concentric arcs (u_turns), since every lane's offset changes
            # by about the same amount there
            e["u_turn_arm"] = True
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
    for e in edges:                             # a twin moves with its arm (beside it), changing lane in the turn
        if "twin_of" in e and "end_shift" in edges[e["twin_of"]]:
            e["end_shift"] = edges[e["twin_of"]]["end_shift"]


HOLD = {}       # edge -> (node, metres from it, {line: lane left of u -> v}): lanes held past a terminus
HOLD_M = 200    # ... this far past the junction (the capsule slid along the road sits within it)


def hold_lanes(edges, orders, stations, spots):
    """Lines passing a terminus where other lines end (H, D at Terebelska Pętla, whose capsule is slid north along
    the road past the junction) keep their lanes HOLD_M past the junction: they cross the capsule parallel to the
    lines ending there and only then close up (H onto D). Lines joining there (C) run outside them."""
    for k, (i, base, t, _) in spots.items():
        if not (stations[k]["term"] - edges[i]["lines"]) or not (stations[k]["term"] & edges[i]["lines"]):
            continue
        node = min((edges[i]["u"], edges[i]["v"]), key=lambda n: np.hypot(*(np.array(n) - base)))
        held = {l: lane(orders, i, l) * (1 if edges[i]["v"] == node else -1) for l in edges[i]["lines"]}
        cur, left = i, float(HOLD_M)
        while left > 0:
            nxt = [j for j, e in enumerate(edges) if j != cur and node in (e["u"], e["v"]) and e["lines"] & held.keys()
                   and not e.get("raw")]
            if not nxt:
                break
            j = max(nxt, key=lambda j: len(edges[j]["lines"] & held.keys()))
            e = edges[j]
            sj = 1 if e["u"] == node else -1                      # left of travel away from the node = sj * left
            cont = {l: held[l] for l in e["lines"] if l in held}  # (u -> v)
            if len(cont) < 2:
                break
            offs = dict(cont)
            others = [l for l in orders[j] if l not in cont]
            ci = [orders[j].index(l) for l in cont]
            for l in others:          # joining lines: outside the whole incoming bundle, on their side (C next to G)
                o_ = orders[j].index(l)
                if o_ > max(ci):
                    offs[l] = min(list(offs.values()) + list(held.values())) - SPACING
                elif o_ < min(ci):
                    offs[l] = max(list(offs.values()) + list(held.values())) + SPACING
            c = _SECTIONS[j] if e["u"] == node else _SECTIONS[j][::-1]
            sd = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(c, axis=0), axis=1))])
            upto = left if sd[-1] > left + 30 else math.inf
            if upto < math.inf:                                    # a vertex where the hold ends
                q = int(np.searchsorted(sd, upto))
                f = (upto - sd[q - 1]) / (sd[q] - sd[q - 1])
                c = np.insert(c, q, c[q - 1] + f * (c[q] - c[q - 1]), axis=0)
                _SECTIONS[j] = c if e["u"] == node else c[::-1]
            HOLD[j] = (node, upto, {l: o * sj for l, o in offs.items()})
            left -= sd[-1]
            held = {l: cont[l] for l in cont}
            cur, node = j, (e["v"] if e["u"] == node else e["u"])


BRIDGE = {}     # (edge in, edge out, node) -> lines bridged from their lane before to theirs after the next corner


def bridge_lanes(seq, l, edges, orders, c, E, Q):
    """See hold_inner_lanes: a cubic curve (tangent to both lanes) from line l's lane on the section before the
    junction to its lane on the second piece of the section after it."""
    for (ia, fa), (ib, fb) in zip(seq, seq[1:]):
        n = edges[ia]["v"] if fa else edges[ia]["u"]
        key = (ia, ib, n) if (ia, ib, n) in BRIDGE else (ib, ia, n) if (ib, ia, n) in BRIDGE else None
        if key is None or l not in BRIDGE[key]:
            continue
        rev = key[0] != ia                        # chain runs from the short section's side
        if rev:
            c, E, Q = c[::-1], E[::-1], Q[::-1]
            ia, ib = key[0], key[1]
        sa = _SECTIONS[ia] if edges[ia]["v"] == n else _SECTIONS[ia][::-1]
        sb = _SECTIONS[ib] if edges[ib]["u"] == n else _SECTIONS[ib][::-1]
        if len(sb) >= 3:
            d1 = sa[-1] - sa[-2]
            la = float(np.hypot(*d1))
            d1 = d1 / la
            d2 = sb[2] - sb[1]
            lb = float(np.hypot(*d2))
            d2 = d2 / lb
            o_a = lane(orders, ia, l) * (1 if edges[ia]["v"] == n else -1)
            o_b = lane(orders, ib, l) * (1 if edges[ib]["u"] == n else -1)
            p1 = sa[-1] + np.array([-d1[1], d1[0]]) * o_a - d1 * min(150.0, 0.5 * la)
            p2 = sb[1] + np.array([-d2[1], d2[0]]) * o_b + d2 * min(180.0, 0.5 * lb)
            if "l_bridge" in FIX:       # end on the short piece, before the next corner (B's arc there kept clear)
                d2 = sb[1] - sb[0]
                d2 = d2 / max(np.hypot(*d2), 1e-9)
                p2 = sb[1] + np.array([-d2[1], d2[0]]) * o_b - d2 * 60
            ka = [k for k in range(len(c)) if E[k] == ia]
            kb = [k for k in range(len(c)) if E[k] == ib]
            if ka and kb:
                i1 = min(ka, key=lambda k: np.hypot(*(c[k] - p1)))
                i2 = min(kb, key=lambda k: np.hypot(*(c[k] - p2)))
                if i1 < i2 and np.hypot(*(c[i1] - p1)) < 15 and np.hypot(*(c[i2] - p2)) < 15:
                    L_ = float(np.hypot(*(p2 - p1)))
                    b1, b2 = p1 + d1 * L_ * 0.4, p2 - d2 * L_ * 0.4
                    m = max(8, int(L_ / STEP_M))
                    tt = np.linspace(0, 1, m)[:, None]
                    new = ((1 - tt) ** 3) * p1 + 3 * ((1 - tt) ** 2) * tt * b1 + 3 * (1 - tt) * tt ** 2 * b2 + tt ** 3 * p2
                    newE = np.array([ia] * (m // 2) + [ib] * (m - m // 2))
                    c = np.vstack([c[:i1], new, c[i2 + 1:]])
                    E = np.concatenate([E[:i1], newE, E[i2 + 1:]])
                    Q = np.vstack([Q[:i1], new, Q[i2 + 1:]])
        if rev:
            c, E, Q = c[::-1], E[::-1], Q[::-1]
    return c, E, Q


def hold_inner_lanes(edges, orders):
    """A line on the inside of a sharp junction turn (over 45 degrees) that would have to move more than 1.5 lanes
    towards the centre onto a short first piece (under 250 m, the next corner close: L from the 9-line bundle on
    Narutowicza north onto al. Tysiąclecia, 180 m before Brzeska) cannot get from its lane to the short piece's
    lane in between (a stair-step, or a hook if held): it is bridged (BRIDGE) straight from its lane before the
    junction to its lane after the next corner by one smooth curve (draw_path, bridge_lanes)."""
    BRIDGE.clear()
    at = defaultdict(list)
    for i, e in enumerate(edges):
        if e["lines"] and i in _SECTIONS and not e.get("raw"):
            at[e["u"]].append(i)
            at[e["v"]].append(i)

    def away(i, n):
        c = _SECTIONS[i] if edges[i]["u"] == n else _SECTIONS[i][::-1]
        v = c[1] - c[0]
        return c, v / max(np.hypot(*v), 1e-9)
    for n, es in at.items():
        for a in es:
            for b in es:
                if a == b or b in HOLD:
                    continue
                ca, va = away(a, n)
                cb, vb = away(b, n)
                d1 = -va                                           # travel into n along a, out along b
                th = math.acos(float(np.clip(d1 @ vb, -1, 1)))
                if th < math.radians(45) or len(cb) < 3:
                    continue
                turn = 1 if d1[0] * vb[1] - d1[1] * vb[0] > 0 else -1
                seg = float(np.hypot(*(cb[1] - cb[0])))
                if seg >= 250:
                    continue
                for l in edges[a]["lines"] & edges[b]["lines"]:
                    o_a = lane(orders, a, l) * (1 if edges[a]["v"] == n else -1)   # left of travel into n
                    o_b = lane(orders, b, l) * (1 if edges[b]["u"] == n else -1)   # left of travel out of n
                    if turn * o_a > 0 and turn * o_a - turn * o_b > 1.5 * SPACING:
                        if "e_bridge" in FIX and l == "E":
                            continue            # E: drawn as the bundle (no bridge pressing it against B)
                        BRIDGE.setdefault((a, b, n), set()).add(l)


def draw_path(seq, l, edges, orders):
    """Drawn coordinates of one chain (subway_like_map style): each section simplified, corners (inside sections
    and at junctions) replaced by circular arcs of radius CORNER_R, the lane offset eased over RAMP_M m wherever
    it changes, then applied along the normals, so a bundle bends as concentric arcs."""
    pts, seg_off, seg_edge, seg_k = [], [], [], []
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
        if i in HOLD and l in HOLD[i][2] and len(c) > 1:   # lanes held past a terminus (hold_lanes)
            node, upto, offd = HOLD[i]
            sd = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(c, axis=0), axis=1))])
            if (edges[i]["u"] if fwd else edges[i]["v"]) == node:
                inside = sd <= upto + 1e-6                          # segment ending at the point
            else:
                inside = np.concatenate([[False], (sd[-1] - sd[:-1]) <= upto + 1e-6])
            inside[0] = inside[1]
            offs = np.where(inside, offd[l] * (1 if fwd else -1), offs)
        if "shift_prof" in edges[i]:                # whole bundle shifted (the U's legs, side by side at the top)
            pr = np.array(edges[i]["shift_prof"], float)
            offs = offs + (pr if fwd else -pr[::-1])
        nseg = len(c) - 1
        seg_k.extend([(i, (t if fwd else nseg - 1 - t), fwd) for t in range(nseg)])
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
    # vertices inside a precise section (a line drawn along its real streets, PRECISE_LINES) get small arcs
    ns = len(seg_edge)
    Rv = np.array([PRECISE_R if all(is_precise(edges[seg_edge[q]]) for q in (v - 1, v) if 0 <= q < ns)
                   else CORNER_R for v in range(len(P))], float)
    # a corner is never tighter than the widest bundle meeting there (lanes and lane shifts included), so its inner
    # lanes do not fold where a shifted arm eases back (B, C, I from Łomaska into Sidorska; A, H at al. Jana Pawła
    # II); same radius for every line through it, so the lanes stay concentric
    # (at the Dworzec PKP crossing: an arm turning into a leg of the U, whose first 300 m are straight, may use 3/4
    # of that straight, so the shifted Sidorska arm's inner lanes get room)
    for v in range(1, len(P) - 1):
        adj = [seg_edge[q] for q in (v - 1, v) if 0 <= q < ns]
        Rv[v] = vertex_R(edges, orders, adj)
    # how much of each segment an arc may use: shared between the corners at its two ends by what they need
    # (TLCAP), not half each, so a sharp turn next to a short piece whose other end is straight gets room
    capA = np.array([TLCAP.get((i, k), (None, None))[0 if fw else 1] for i, k, fw in seg_k], dtype=object)
    capB = np.array([TLCAP.get((i, k), (None, None))[1 if fw else 0] for i, k, fw in seg_k], dtype=object)
    segl = np.linalg.norm(np.diff(P, axis=0), axis=1)
    capA = np.array([segl[q] * 0.5 if capA[q] is None else capA[q] for q in range(len(segl))], float)
    capB = np.array([segl[q] * 0.5 if capB[q] is None else capB[q] for q in range(len(segl))], float)
    F, sid = fillet(P, Rv, (capA, capB))
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(F, axis=0), axis=1))])
    keep = np.concatenate([[True], np.diff(s) > 1e-6])
    F, s, sid = F[keep], s[keep], np.array(sid)[keep]
    if s[-1] < 1:
        return F, np.full(len(F), seq[0][0])
    ds = np.append(np.arange(0, s[-1], STEP_M), s[-1])
    Q = np.column_stack([np.interp(ds, s, F[:, 0]), np.interp(ds, s, F[:, 1])])
    k_ = np.clip(sid[np.clip(np.searchsorted(s, ds, side="right") - 1, 0, len(sid) - 1)], 0, len(seg_off) - 1)
    O, E = seg_off[k_].astype(float), seg_edge[k_]
    # a lane change at a junction happens on the side with fewer lines: a line leaving a bundle moves over once it
    # has left it, one joining a bundle before it joins (the ease is centred on the junction otherwise, and A slid
    # over H for 60 m before turning off al. Jana Pawła II)
    half = RAMP_M / 2
    O2 = O.copy()
    bnd = [j for j in range(1, len(E)) if E[j] != E[j - 1]]
    # an arm turning into a leg of the U (u_turns redraws that turn): its lane change happens on the arm, inside the
    # redrawn stretch, so the leg is settled right below the crossing (its straight top is short)
    if U_INFO:
        legs = (U_INFO["left"], U_INFO["right"])
        for j in bnd:
            a_, b_ = E[j - 1], E[j]
            if (a_ in legs) == (b_ in legs):
                continue
            arm = b_ if a_ in legs else a_
            if not edges[arm].get("u_turn_arm"):
                continue
            ca_ = _SECTIONS[arm]
            seg_a = float(np.hypot(*(ca_[-1] - ca_[-2]))) if edges[arm]["v"] == U_INFO["node"] else                 float(np.hypot(*(ca_[1] - ca_[0])))
            reach = max(0.0, min(U_ARM_LEAD_M, 0.5 * seg_a) - 2 * RAMP_M)
            sel = (E == arm) & (np.abs(ds - ds[j]) < reach)
            O2[sel] = O[j] if b_ in legs else O[j - 1]
    for n_, j in enumerate(bnd):
        na, nb = len(edges[E[j - 1]]["lines"]), len(edges[E[j]]["lines"])
        if not LANE_CHANGE_SIDE or na == nb or O[j] == O[j - 1]:
            continue
        lo_ = ds[bnd[n_ - 1]] if n_ else -math.inf
        hi_ = ds[bnd[n_ + 1]] if n_ + 1 < len(bnd) else math.inf
        if na > nb:                                   # leaving: hold the old lane a little into the new section
            sel = (ds >= ds[j]) & (ds < min(ds[j] + half, hi_))
            O2[sel] = O[j - 1]
        else:                                         # joining: take the new lane a little before it
            sel = (ds < ds[j]) & (ds >= max(ds[j] - half, lo_))
            O2[sel] = O[j]
    # lanes spreading from the incoming bundle onto the dead-end spur's legs ease over a longer stretch
    O = smooth(O2[:, None], RAMP_M)[:, 0]
    if any(edges[i].get("raw") and "shift_prof" in edges[i] and "down_dir" not in edges[i] for i, _ in seq):
        # a schematic section with its own shift profile (Kościuszki on the spur: its centreline moves exactly as its
        # shift fades): that profile is followed closely, not eased over RAMP_M (which bent the lanes)
        own = np.array([1.0 if edges[e_].get("raw") and "shift_prof" in edges[e_] and "down_dir" not in edges[e_]
                        else 0.0 for e_ in E])
        for e_ in set(E[own > 0].tolist()):             # back to the normal easing near its far junction
            vend = np.array(edges[e_]["v"])
            own = own * np.where(E == e_, np.clip(np.hypot(*(Q - vend).T) / 150.0, 0, 1), 1.0)
        own = smooth(own[:, None], 40)[:, 0]
        O = O * (1 - own) + smooth(O2[:, None], 24)[:, 0] * own
    if U_INFO.get("spur") and any("down_dir" in edges[i] for i, _ in seq):    # only near the spur's top
        d_ = np.hypot(*(Q - np.array(U_INFO["node"])).T)
        w_ = np.clip((450 - d_) / 150, 0, 1)
        w_ = w_ * np.array([0.0 if edges[e_].get("raw") and "down_dir" not in edges[e_] else 1.0 for e_ in E])
        w_ = smooth(w_[:, None], 30)[:, 0]          # (not on Kościuszki: drawn with exact lanes, raw)
        if "spur_east" in FIX:      # not where a line turns off the east leg (A onto Sidorska east): its turn is exact
            w_ = w_ * np.array([0.0 if (e_ == U_INFO["right"] or edges[e_].get("u_turn_arm") and
                                        edges[e_]["geom"].distance(Point(U_INFO["node"])) < 1 and
                                        np.dot(np.array(edges[e_]["geom"].centroid.coords[0]) - np.array(U_INFO["node"]),
                                               U_INFO["perp"]) > 0) else 1.0 for e_ in E])
            w_ = smooth(w_[:, None], 60)[:, 0]
        O = O * (1 - w_) + smooth(O2[:, None], SPUR_RAMP_M)[:, 0] * w_
    t = np.gradient(Q, axis=0)
    t = smooth(t, STEP_M * 3)
    t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-9)[:, None]
    nrm = np.stack([-t[:, 1], t[:, 0]], axis=1)
    c = Q + nrm * O[:, None]
    if BRIDGE:
        c, E, Q = bridge_lanes(seq, l, edges, orders, c, E, Q)
    if EXACT_TURNS:
        c, E, Q = exact_turns(seq, l, edges, orders, c, E, Q)
    if U_INFO:
        c, E, Q = u_turns(seq, l, edges, orders, c, E, Q)
    return c, E, Q


def exact_turns(seq, l, edges, orders, c, E, Q):
    """An inner lane that would fold at a tight junction corner (the corner's arc smaller than the lane's distance
    inside it: L on the inside of the 9-line bundle at Zamkowa / Narutowicza and al. Tysiąclecia, where the next
    junction leaves too little road for a wide arc) is redrawn from its lane on the street before to its lane on
    the street after as straight - circular arc - straight (radius: the corner's centreline radius minus the line's
    mean offset inside, at least 0.8 lane), instead of looping back on itself (a V, a stair-step). Only such
    lanes: the others keep the bundle's concentric arcs."""
    for (ia, fa), (ib, fb) in zip(seq, seq[1:]):
        ea, eb = edges[ia], edges[ib]
        if any(e_.get("raw") or "end_shift" in e_ or "shift_prof" in e_ for e_ in (ea, eb)):
            continue
        if any(k_ in HOLD and l in HOLD[k_][2] for k_ in (ia, ib)):
            continue
        o_in = lane(orders, ia, l) * (1 if fa else -1)
        o_out = lane(orders, ib, l) * (1 if fb else -1)
        sa = _SECTIONS[ia] if fa else _SECTIONS[ia][::-1]
        sb = _SECTIONS[ib] if fb else _SECTIONS[ib][::-1]
        if len(sa) < 2 or len(sb) < 2:
            continue
        Nn = sa[-1]
        d1, d2 = sa[-1] - sa[-2], sb[1] - sb[0]
        seg_in, seg_out = float(np.hypot(*d1)), float(np.hypot(*d2))
        d1, d2 = d1 / seg_in, d2 / seg_out
        th = math.acos(float(np.clip(d1 @ d2, -1, 1)))
        if th < math.radians(20) or th > math.radians(170):
            continue
        turn = 1 if d1[0] * d2[1] - d1[1] * d2[0] > 0 else -1
        n1, n2 = np.array([-d1[1], d1[0]]), np.array([-d2[1], d2[0]])
        kin = len(_SECTIONS[ia]) - 2 if fa else 0
        kout = 0 if fb else len(_SECTIONS[ib]) - 2
        cap_in = TLCAP.get((ia, kin), (0.5 * seg_in, 0.5 * seg_in))[1 if fa else 0]
        cap_out = TLCAP.get((ib, kout), (0.5 * seg_out, 0.5 * seg_out))[0 if fb else 1]
        tan_ = math.tan(th / 2)
        tl_c = min(vertex_R(edges, orders, [ia, ib]) * tan_, cap_in, cap_out)
        M = np.array([d1, -d2]).T
        s_, t_ = np.linalg.solve(M, (Nn + n2 * o_out) - (Nn + n1 * o_in))
        C = Nn + n1 * o_in + d1 * s_                 # where the two lanes meet; s_ along d1 from the node
        t_ = -t_                                     # ... and C = Nn + n2 * o_out - d2 * t_
        R_eff = tl_c / tan_
        if min(R_eff - turn * o_in, R_eff - turn * o_out) >= SPACING:
            continue                                 # the lane bends on its own arc with room: drawn as is
        r = max(R_eff - turn * (o_in + o_out) / 2, 0.8 * SPACING)
        tl = r * tan_
        back = min(max(tl_c, 90.0, tl - s_), cap_in)          # how far before the node the redrawn piece starts
        fwd_ = min(max(tl_c, 90.0, tl - t_), cap_out)
        tl = min(tl, s_ + back, t_ + fwd_)
        if tl <= 1:
            continue
        r = tl / tan_
        p1 = Nn + n1 * o_in - d1 * back
        p2 = Nn + n2 * o_out + d2 * fwd_
        ka = [k for k in range(len(c)) if E[k] == ia]
        kb = [k for k in range(len(c)) if E[k] == ib]
        if not ka or not kb:
            continue
        i1 = min(ka, key=lambda k: np.hypot(*(c[k] - p1)))
        i2 = min(kb, key=lambda k: np.hypot(*(c[k] - p2)))
        if not (i1 < i2) or np.hypot(*(c[i1] - p1)) > 12 or np.hypot(*(c[i2] - p2)) > 12:
            continue
        T1, T2 = C - d1 * tl, C + d2 * tl
        ctr = T1 + n1 * turn * r
        a1 = math.atan2(*(T1 - ctr)[::-1])
        m1 = max(2, int(np.hypot(*(T1 - c[i1])) / STEP_M))
        m3 = max(2, int(np.hypot(*(c[i2] - T2)) / STEP_M))
        m2 = max(3, int(r * th / STEP_M))
        new = np.array([c[i1] + (T1 - c[i1]) * f for f in np.linspace(0, 1, m1, endpoint=False)]
                       + [ctr + r * np.array([math.cos(a1 + turn * th * f), math.sin(a1 + turn * th * f)])
                          for f in np.linspace(0, 1, m2, endpoint=False)]
                       + [T2 + (c[i2] - T2) * f for f in np.linspace(0, 1, m3)])
        newE = np.array([ia] * (m1 + m2 // 2) + [ib] * (len(new) - m1 - m2 // 2))
        c = np.vstack([c[:i1], new, c[i2 + 1:]])
        E = np.concatenate([E[:i1], newE, E[i2 + 1:]])
        Q = np.vstack([Q[:i1], new, Q[i2 + 1:]])
    return c, E, Q


U_TURN_R_IN = 60    # m: radius of the innermost lane of an arm turning into a leg of the U
U_ARM_LEAD_M = 350  # ... the turn is redrawn from this far back along the arm


def u_turns(seq, l, edges, orders, c, E, Q):
    """An arm turning into a leg of the Dworzec PKP U (Sidorska west / east): line l's lane on the arm (centred on
    its street) and its lane on the leg are joined by a straight, a circular arc and a straight, the arc's radius
    growing with the lane's distance from the inside of the turn, so the arm's lanes turn on concentric arcs
    (their offsets change by about the same amount from street to leg, so the arcs share a centre) and nothing
    is shifted along the street."""
    N2, legs = U_INFO["node"], (U_INFO["left"], U_INFO["right"])
    for (ia, fa), (ib, fb) in zip(seq, seq[1:]):
        if (edges[ia]["v"] if fa else edges[ia]["u"]) != N2:
            continue
        if U_INFO.get("spur") and U_INFO["right"] in (ia, ib):
            continue                  # off the spur's east leg: eased with its neighbours (A next to B, C)
        if ia in legs and edges[ib].get("u_turn_arm"):
            arm, leg, flip = ib, ia, True
        elif ib in legs and edges[ia].get("u_turn_arm"):
            arm, leg, flip = ia, ib, False
        else:
            continue
        if flip:
            c, E, Q = c[::-1], E[::-1], Q[::-1]
        sd = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(c, axis=0).T))])
        ia_ = np.where(E == arm)[0]
        jb = [k for k in np.where(E == leg)[0] if ia_.size and k > ia_.min()]
        done = False
        if ia_.size and jb:
            j = jb[0]
            # the arm's lane, exactly: its last straight piece into the crossing, line l's lane on it, from at most
            # half that piece back (the other half may be the corner at its far end)
            ca = _SECTIONS[arm] if edges[arm]["v"] == N2 else _SECTIONS[arm][::-1]
            d1 = ca[-1] - ca[-2]
            seg_a = np.hypot(*d1)
            d1 = d1 / seg_a
            a_l = lane(orders, arm, l) * (1 if edges[arm]["v"] == N2 else -1)
            p1 = np.array(N2) - d1 * min(U_ARM_LEAD_M, 0.5 * seg_a) + np.array([-d1[1], d1[0]]) * a_l
            cand1 = [k for k in range(0, j) if E[k] == arm]
            i1 = min(cand1, key=lambda k: np.hypot(*(c[k] - p1))) if cand1 else -1
            if 4 <= i1 and np.hypot(*(c[i1] - p1)) < 15:
                p1 = c[i1].copy()                          # join the arm exactly as drawn there
                d1 = c[i1] - c[i1 - 4]
                d1 = d1 / max(np.hypot(*d1), 1e-9)
            # the leg's lane, exactly: its straight top (U_TOP_M from the crossing), line l's lane on it
            gl = np.array(edges[leg]["geom"].coords)
            down_u = edges[leg]["u"] == N2                         # left leg: drawn from the crossing down
            d2 = edges[leg]["down_dir"] if "down_dir" in edges[leg] else                 ((gl[1] - gl[0]) if down_u else (gl[-2] - gl[-1]))
            d2 = d2 / np.hypot(*d2)
            pr_ = edges[leg]["shift_prof"]
            b_ = (lane(orders, leg, l) + pr_[0]) if down_u else -(lane(orders, leg, l) + pr_[-1])
            p2 = np.array(N2) + d2 * (U_TOP_M - 75) + np.array([-d2[1], d2[0]]) * b_   # before the S-bend's ease
            cand = [k for k in range(j, len(c)) if E[k] == leg]
            i2 = min(cand, key=lambda k: np.hypot(*(c[k] - p2))) if cand else -1
            if 0 <= i2 < len(c) - 4 and np.hypot(*(c[i2] - p2)) < 15:
                p2 = c[i2].copy()                          # join the leg exactly as drawn there
                d2 = c[i2 + 4] - c[i2]
                d2 = d2 / max(np.hypot(*d2), 1e-9)
            if 0 <= i1 < j < i2 < len(c) - 1 and np.hypot(*(c[i2] - p2)) < 15 and np.hypot(*(c[i1] - p1)) < 15:
                M = np.array([d1, -d2]).T
                if abs(np.linalg.det(M)) > 1e-6:
                    s_, t_ = np.linalg.solve(M, p2 - p1)
                    t_ = -t_                                       # C = p2 - t_ * d2: before p2
                    C = p1 + s_ * d1
                    th = math.acos(float(np.clip(d1 @ d2, -1, 1)))
                    turn = 1 if d1[0] * d2[1] - d1[1] * d2[0] > 0 else -1     # +1: a left turn
                    sg = 1 if edges[arm]["v"] == N2 else -1
                    out = lambda a: -a * turn  # noqa: E731  (left lanes are outside a right turn)
                    a_all = [lane(orders, arm, x) * sg for x in edges[arm]["lines"]]
                    r = U_TURN_R_IN + out(lane(orders, arm, l) * sg) - min(out(a) for a in a_all)
                    tl = r * math.tan(th / 2)
                    if s_ > 0 and t_ > 0 and th > math.radians(20):
                        if tl > 0.95 * min(s_, t_):
                            tl = 0.95 * min(s_, t_)
                            r = tl / math.tan(th / 2)
                        T1, T2 = C - d1 * tl, C + d2 * tl
                        nl = np.array([-d1[1], d1[0]]) * turn
                        ctr = T1 + nl * r
                        a1 = math.atan2(*(T1 - ctr)[::-1])
                        sweep = turn * th
                        n1 = max(2, int(np.hypot(*(T1 - p1)) / STEP_M))
                        n3 = max(2, int(np.hypot(*(p2 - T2)) / STEP_M))
                        n2 = max(3, int(r * th / STEP_M))
                        seg1 = [p1 + (T1 - p1) * f for f in np.linspace(0, 1, n1, endpoint=False)]
                        arc = [ctr + r * np.array([math.cos(a1 + sweep * f), math.sin(a1 + sweep * f)])
                               for f in np.linspace(0, 1, n2, endpoint=False)]
                        seg3 = [T2 + (p2 - T2) * f for f in np.linspace(0, 1, n3)]
                        new = np.array(seg1 + arc + seg3)
                        newE = np.array([arm] * (n1 + n2 // 2) + [leg] * (len(new) - n1 - n2 // 2))
                        c = np.vstack([c[:i1], new, c[i2 + 1:]])
                        E = np.concatenate([E[:i1], newE, E[i2 + 1:]])
                        Q = np.vstack([Q[:i1], new, Q[i2 + 1:]])
                        done = True
        if flip:
            c, E, Q = c[::-1], E[::-1], Q[::-1]
        _ = done
    return c, E, Q


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
        if not cand:      # a stop beyond the drawn end of its line (Ekologiczna, 217 m past E's end): still placed
            cand = [i for i in tree.query(pt.buffer(250)) if edges[i]["lines"] & s["lines"]
                    and geoms[i].distance(pt) < 250]
        cand = [i for i in cand if not stub(i, pt)] or cand
        if not cand:
            continue
        # a terminus sits on the section carrying most of the lines ending there (Terebelska Pętla: the trunk with
        # I, M, D, H, G, not the short northern piece with D, H, C), else the nearest
        dmin = min(geoms[i].distance(pt) for i in cand)        # (only as a tie-break among the nearest)
        i = min(cand, key=lambda i: (geoms[i].distance(pt) > dmin + 5, -len(edges[i]["lines"] & s["term"]),
                                     geoms[i].distance(pt)))
        g = geoms[i]
        d = g.project(pt)
        h = heading_at(g, d)
        if i in _SECTIONS and len(_SECTIONS[i]) >= 2 and min(d, g.length - d) < 80:
            # near a section's end (a dead end with a turning loop, Brzegowa - Baza): the direction of the section
            # as drawn there, not of the raw street's last wiggle
            ds_ = LineString(_SECTIONS[i])
            q_ = ds_.project(pt)
            h = heading_at(ds_, min(max(q_, 40), ds_.length - 40))
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
    "Sitnicka": "Sitnicka", "Robotnicza": "Robotnicza",
    "Armii Krajowej": "Armii Krajowej", "Witoroska": "Witoroska", "Aleja Riada Haidara": "al. Haidara",
}
STREET_SIZE, STREET_COLOR = 10.5, "#8e8e8e"
STREET_NEAR_CENTRE = {"Sidorska", "Francuska", "Lubelska", "Janowska"}
OFFSET_OK = set()   # street names allowed further out from their bundle
# names placed by hand (map x, y, angle): where the drawn lines run off the street (C / K on Sobieskiego)
STREET_MANUAL = {}
STREET_PIN = {"Zamkowa": (782262, 471235), "Gabriela Narutowicza": (782680, 470833)}  # names kept here (user)


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
    keep = list(range(len(names)))            # every named street: the one really under the line decides
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
            if otree.geometries[k].distance(p) > 30 or names[keep[k]] not in STREET_LABELS:
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
    centre = np.mean([np.array(e["geom"].interpolate(0.5, normalized=True).coords[0]) * len(e["lines"])
                      for e in edges], axis=0) / np.mean([len(e["lines"]) for e in edges])
    for street, (mx, my, ma) in STREET_MANUAL.items():
        text = STREET_LABELS[street]
        ax.text(mx, my, text, fontproperties=font(STREET_SIZE), color=STREET_COLOR, ha="center", va="center",
                rotation=ma, rotation_mode="anchor", zorder=6.5,
                path_effects=[patheffects.withStroke(linewidth=2.5, foreground=LAND)])
        w = len(text) * STREET_SIZE * 0.55 * m_per_pt
        placed.append(Point(mx, my).buffer(w / 2))
        n_done += 1
    for street, cs in sorted(cands.items(), key=lambda kv: -len(kv[1])):
        if street in STREET_MANUAL:
            continue
        text = STREET_LABELS[street]
        w = len(text) * STREET_SIZE * 0.55 * m
        h = STREET_SIZE * 1.15 * m
        # spots nearest the middle of the street's longest section first; the preferred side (below a horizontal
        # street, right of a vertical one) anywhere along the street before the other side
        if street in STREET_PIN:                       # placed where the user liked them
            cs = sorted(cs, key=lambda c_: np.hypot(*(c_[1] - np.array(STREET_PIN[street]))))
        elif street in STREET_NEAR_CENTRE:             # long radial streets: named close to the city
            cs = sorted(cs, key=lambda c_: np.hypot(*(c_[1] - centre)))
        else:
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
                ok_ = False
                # lines drawn off their street's centre (held lanes: C at Terebelska Pętla) push the name out
                for extra in ((0.0,) if street not in OFFSET_OK else (0.0, 1.5, 3.0, 4.5)):
                    cen = p + nrm * sd_ * (half + SPACING * (0.5 + extra) + h / 2)
                    poly = Polygon([cen + t * w / 2 + nrm * h / 2, cen - t * w / 2 + nrm * h / 2,
                                    cen - t * w / 2 - nrm * h / 2, cen + t * w / 2 - nrm * h / 2])
                    if not frame.contains(poly):
                        continue
                    tol = 120 if street in STREET_PIN else 0   # a pinned name may graze a line's edge (box corner)
                    if sum(occ[k].intersection(poly).area for k in occ_tree.query(poly)) > tol:
                        continue
                    if any(q.buffer(12).intersects(poly) for q in placed):
                        continue
                    ok_ = True
                    break
                if not ok_:
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


CASED_LINES = {"L"}    # lines drawn with a fine black outline (L's yellow on the cream paper)
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
    if U_INFO.get("spur"):                       # the station line across the spur's bottom: radially, down
        for k_, (i_, b_, t_, ls_) in list(spots.items()):
            if stations[k_]["name"] in KEY_STATIONS:
                spots[k_] = (U_INFO["left"], b_, U_INFO["perp"], ls_)
    hold_lanes(edges, orders, stations, spots)
    hold_inner_lanes(edges, orders)
    compute_tl_caps(edges, orders)
    ending_cut = defaultdict(set)                 # terminus station -> lines cut there
    term_cut = defaultdict(list)                  # terminus station -> points where lines ending there were cut
    shifted = {}      # terminus station -> (station point moved along the road, tip: the peak of a side line's arc)
    def pale_end(l, k):
        """True if line l ending at station k ends on an extension: under EXT_SHARE of its trips call there (counted
        at the stop, both directions: near a terminal loop one point sees only one direction)."""
        return stations[k]["trips"].get(l, 0) < EXT_SHARE * line_trips[l]

    def trim_at_terminus(l, c, E, Q, at_end):
        """A line ending at one of its main termini is cut where it crosses the station (the line through the
        stop across the bundle), so every line ending there stops flush at the station capsule instead of
        wherever its drawing happens to end (Terebelska Pętla had six end bars scattered round a corner).
        Returns (c, E, Q, station or None)."""
        xy = c[-1] if at_end else c[0]
        cand = [k for k, s_ in stations.items() if k in spots and np.hypot(*(s_["xy"] - xy)) < 300
                and (l in s_["term"] or k in END_STATIONS.get((l, clean_stop(s_["name"])), {}))]
        if not cand:
            return c, E, Q, None
        k = min(cand, key=lambda k: np.hypot(*(stations[k]["xy"] - xy)))
        _, base, t, _ = spots[k]
        if k in shifted:                          # the station line moved to the arc's peak, level (see below)
            base, _, t = shifted[k]
        elif k in tilted:                         # cut along the tilted station line instead
            pa, pb, u = tilted[k]
            t2 = np.array([u[1], -u[0]])
            t, base = (t2 if t2 @ t >= 0 else -t2), (pa + pb) / 2
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
                ending_cut[k].add(l)
                return c, E, Q, k
        # ends just short of the station (I at Chłodna 1: inbound ends 17 m before it, outbound starts 17 m after):
        # carried on along the end's direction to the station line, so both meet at the capsule
        tip, prev = (c[-1], c[-2]) if at_end else (c[0], c[1])
        dv_ = tip - prev
        dv_ = dv_ / max(np.hypot(*dv_), 1e-9)
        st_ = -((tip - base) @ t) / (dv_ @ t) if abs(dv_ @ t) > 0.3 else -1
        if 0 <= st_ <= (250 if k in shifted else 60):
            p = tip + dv_ * st_
            if abs((p - base) @ nrm) <= 4 * SPACING + 40:
                if at_end:
                    c, E, Q = np.vstack([c, p]), np.append(E, E[-1]), np.vstack([Q, Q[-1]])
                else:
                    c, E, Q = np.vstack([p, c]), np.append(E[0], E), np.vstack([Q[0], Q])
                term_cut[k].append(p)
                ending_cut[k].add(l)
                return c, E, Q, k
        return c, E, Q, None

    chains = []                                   # (line, chain, pre, post, drawn coords, edge per point, centreline)
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
            chains.append((l, seq, pre, post, c, E, Q))

    # a line ending at a terminus that curves in beside the bundle instead of along it (C at Terebelska Pętla from
    # Jana III Sobieskiego, curving north past the trunk's end): the station line stays square across the road and
    # is slid along it to where that line's curve runs exactly one lane outside the outermost line ending there on
    # its side (C next to G), so the capsule reaches out to it as one bar
    def crossings(c, b0, t_, across=False):
        """Positions across the road (along nrm) where the drawn polyline c crosses the line through b0 square to t_."""
        n_ = np.array([-t_[1], t_[0]])
        sv = (c - b0) @ t_
        out = []
        for j in range(1, len(c)):
            if sv[j - 1] == 0 or np.sign(sv[j - 1]) != np.sign(sv[j]):
                f = sv[j - 1] / (sv[j - 1] - sv[j]) if sv[j - 1] != sv[j] else 0
                q = c[j - 1] + (c[j] - c[j - 1]) * f
                if across:                                     # only where the line runs across the road
                    d_ = c[j] - c[j - 1]
                    if abs(d_ @ t_) > 0.7 * np.hypot(*d_):
                        continue
                out.append((float((q - b0) @ n_), q))
        return out

    for k, (i, base, t, _) in spots.items():
        extra = stations[k]["term"] - edges[i]["lines"]
        ending = stations[k]["term"] & edges[i]["lines"]
        if not extra or not ending:
            continue
        # the road's direction over its last 150 m into the junction (not the last few metres, which may bend)
        g_ = edges[i]["geom"]
        node = min((edges[i]["u"], edges[i]["v"]), key=lambda n: np.hypot(*(np.array(n) - base)))
        back = np.array(g_.interpolate(150 if edges[i]["u"] == node else g_.length - 150).coords[0])
        t = np.array(node) - back
        t /= np.hypot(*t)
        base = np.array(node, float)
        # where the lines ending there cross the road 40 m before the junction (they run straight on from there: past
        # the junction they are carried on along the road to wherever the capsule goes)
        hit0 = {l: x for l, _, _, _, c, _, _ in chains if l in ending
                for x, _ in crossings(c, base - t * 40, t) if abs(x) < 4 * SPACING + 60}
        if set(hit0) != ending:
            continue
        ends_x = list(hit0.values())
        up = t                                                     # along the road, past the junction (north)
        best = None
        for s_ in np.arange(-40, 301, 2.0):                          # as far along the road as C still fits by G
            b0 = base + up * s_
            for l, _, _, _, c, _, _ in chains:
                if l not in extra:
                    continue
                for x, q in crossings(c, b0, t, across=True):     # C coming in across the road, not its leg north
                    gap = (x - max(ends_x)) if x > max(ends_x) else (min(ends_x) - x) if x < min(ends_x) else -1
                    if gap < 0:
                        continue
                    sc = abs(gap - SPACING) - 0.15 * s_            # C one lane beyond G, as far north as possible
                    if abs(gap - SPACING) < 0.7 * SPACING and (best is None or sc < best[0]):
                        best = (sc, b0, q)
        if best is not None:
            shifted[k] = (best[1], best[2], t)

    tilted = {}
    for _ in [0]:
        pass                                   # station -> (end a, end b, unit direction a -> b)
    for k, (i, base, t, _) in spots.items():
        extra = stations[k]["term"] - edges[i]["lines"]
        if not extra or k in shifted:
            continue
        nrm = np.array([-t[1], t[0]])
        offs = [lane(orders, i, l) for l in edges[i]["lines"]]
        pa, pb = base + nrm * min(offs), base + nrm * max(offs)
        mid = (min(offs) + max(offs)) / 2
        for l in extra:
            best = None
            for j, e in enumerate(edges):
                if l in e["lines"] and j in _SECTIONS and e["geom"].distance(Point(base)) < 200:
                    g = LineString(_SECTIONS[j]).offset_curve(lane(orders, j, l))
                    if not g.is_empty and (best is None or g.distance(Point(base)) < best[0]):
                        best = (g.distance(Point(base)), np.array(g.interpolate(g.project(Point(base))).coords[0]))
            if best is None or best[0] > 3 * SPACING + 60:
                continue
            q = best[1]
            if float((q - base) @ nrm) >= mid:
                pb = q
            else:
                pa = q
        u = (pb - pa) / max(np.hypot(*(pb - pa)), 1e-9)
        tilted[k] = (pa, pb, u)


    for l, seq, pre, post, c, E, Q in chains:
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
            # solid / pale node to node: one value per piece between junctions and the line's termini (the station
            # lines), taken by length from the point-by-point mask, so the pale part starts exactly at the capsule
            # (Sidorska Pętla) and not wherever the shapes happen to thin out
            cuts = {j for j in range(1, len(c)) if E[j] != E[j - 1]}
            for k_, s_ in stations.items():
                if l not in s_["term"] or k_ not in spots:
                    continue
                _, b0, t0, _ = spots[k_]
                if k_ in shifted:
                    b0, _, t0 = shifted[k_]
                n0 = np.array([-t0[1], t0[0]])
                sv = (c - b0) @ t0
                for j in range(1, len(c)):
                    if np.sign(sv[j]) != np.sign(sv[j - 1]) and abs((c[j] - b0) @ n0) < 4 * SPACING + 60                             and np.hypot(*(c[j] - b0)) < 200:
                        cuts.add(j)
            cuts = [0] + sorted(cuts) + [len(c)]
            seg_l = np.append(np.hypot(*np.diff(c, axis=0).T), 0.0)
            for a_, b_ in zip(cuts, cuts[1:]):
                w_ = seg_l[a_:b_]
                if w_.sum() > 0:
                    ext[a_:b_] = (w_ * ext[a_:b_]).sum() / w_.sum() >= 0.5
            for k_, (b2, apex, t2) in shifted.items():           # a line meeting a moved capsule's tip (C at
                if l in stations[k_]["term"]:                     # Terebelska Pętla) is solid up to the capsule
                    near_ = (np.hypot(*(c - apex).T) < 110) & ((c - b2) @ t2 < 0)
                    ext[near_] = False
            cut = [0] + [k for k in range(1, len(c)) if ext[k] != ext[k - 1]] + [len(c)]
            for a_, b_ in zip(cut, cut[1:]):
                part = c[max(a_ - 1, 0):b_]
                if len(part) < 2:
                    continue
                if l in CASED_LINES:              # a light colour on the paper: a fine black outline
                    ax.plot(part[:, 0], part[:, 1], color=INK, lw=lw + 1.3, solid_capstyle="butt",
                            solid_joinstyle="round", zorder=(2.9 if ext[a_] else 3) - 0.01)
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

    # Dworzec PKP drawn as a box over the bottom of the U (box) or as a platform under it (under)
    station_shape = None
    if U_INFO and STATION_STYLE in ("box", "under"):
        N2_, dv_, pp_ = np.array(U_INFO["node"]), U_INFO["dv"], U_INFO["perp"]
        legs_ = (U_INFO["left"], U_INFO["right"])
        depth_ = float((U_INFO["Sx"] - N2_) @ dv_)
        cut_ = U_TOP_M + 30 if STATION_STYLE == "box" else depth_ - U_HALF_M * 1.1
        pts_ = [p for _, _, _, _, c, E_, _ in chains for p, e_ in zip(c, E_)
                if e_ in legs_ and (p - N2_) @ dv_ > cut_]
        if pts_:
            loc = np.array([[(p - N2_) @ pp_, (p - N2_) @ dv_] for p in pts_])
            pad = SPACING * 0.9
            (a0, b0), (a1, b1) = loc.min(axis=0) - pad, loc.max(axis=0) + pad
            to_map = lambda u_, w_: N2_ + pp_ * u_ + dv_ * w_  # noqa: E731
            if STATION_STYLE == "box":
                r_ = 60.0
                shp = box(a0 + r_, b0 + r_, a1 - r_, b1 - r_).buffer(r_, quad_segs=12)
            else:
                h_ = (b1 - b0) / 2
                shp = LineString([(a0 + h_, (b0 + b1) / 2), (a1 - h_, (b0 + b1) / 2)]).buffer(h_, quad_segs=16)
            ring = np.array([to_map(u_, w_) for u_, w_ in shp.exterior.coords])
            z_ = 5.5 if STATION_STYLE == "box" else 2.5
            ax.add_patch(plt.Polygon(ring, closed=True, facecolor="white", edgecolor=INK, linewidth=lw * 0.45,
                                     zorder=z_))
            station_shape = Polygon(ring)
            drawn.append(station_shape)
            if STATION_STYLE == "box":
                cen = np.array(station_shape.centroid.coords[0])
                ax.text(cen[0], cen[1], "Dworzec" + chr(10) + "PKP", fontproperties=font(9.5), color=INK,
                        ha="center", va="center", zorder=5.6, linespacing=1.05)

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
            if pale_end(l, cut_at) or l not in stations[cut_at]["term"]:
                var_end[cut_at].add(l)            # an extension's end: diamond bullet in the name only
            continue
        cand = [k for k, s_ in stations.items() if l in s_["lines"] and np.hypot(*(s_["xy"] - xy)) < 250]
        k = min(cand, key=lambda k: np.hypot(*(stations[k]["xy"] - xy))) if cand else None
        if k is not None:
            ends_at[k].add(l)
        # no end mark where the line turns at Dworzec PKP (it gets its own capsule)
        if any(s_["name"] in KEY_STATIONS and np.hypot(*(s_["xy"] - xy)) < 80 for s_ in stations.values()):
            continue
        nrm = np.array([-d[1], d[0]])               # every line end: a bar (diamonds only in name bullets)
        capsule(xy - nrm * SPACING * 0.5, xy + nrm * SPACING * 0.5)
        if not (k is not None and l in stations[k]["term"] and not pale_end(l, k)):
            if k is not None:
                var_end[k].add(l)

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
        if key_st and U_INFO.get("spur"):          # the spur's bottom: every lane crossing the line down from it
            O_c = U_INFO["Sx"] - U_INFO["dv"] * SPACING
            pts_ = [q for l, _, _, _, c, _, _ in chains if l in edges[i]["lines"]
                    for _, q in crossings(c, O_c, t) if 0 < (q - O_c) @ U_INFO["dv"] < 300]
        if not pts_:
            continue
        side = [float(np.dot(p - base, nrm)) for p in pts_]
        lo, hi = min(side) - SPACING * 0.55, max(side) + SPACING * 0.55
        a_, b_ = base + nrm * lo, base + nrm * hi
        if sid in tilted:                             # capsule along the tilted station line
            pa, pb, u_ = tilted[sid]
            a_, b_ = pa - u_ * SPACING * 0.55, pb + u_ * SPACING * 0.55
        if sid in shifted:                            # moved station line: from the lines cut there to the arc's peak
            b2, apex, t2 = shifted[sid]
            nrm = np.array([-t2[1], t2[0]])
            sd2 = [float((p - b2) @ nrm) for p in term_cut.get(sid, [])] or [0.0]
            for l, _, _, _, c, _, _ in chains:                     # lines passing through it (H, D): spanned too
                if l in edges[i]["lines"] and l not in ending_cut.get(sid, set()):
                    sd2 += [x for x, _ in crossings(c, b2, t2) if abs(x) < 4 * SPACING + 60]
            sa = float((apex - b2) @ nrm)
            far = max(sd2, key=lambda v: abs(v - sa))        # the lane farthest from the peak, padded outwards
            far += SPACING * 0.55 * (1 if far > sa else -1)
            a_, b_ = b2 + nrm * far, b2 + nrm * sa
            base = b2
        if key_st and station_shape is not None:     # Dworzec PKP as a box / platform (drawn above)
            if STATION_STYLE == "box":
                continue
            requests.append((0, -len(st_["lines"]), sid, base, station_shape.buffer(5), []))
            continue
        far_stop = np.hypot(*(st_["xy"] - base)) > 90           # placed by the fallback: the diamond marks it
        if not far_stop and not any(n_ == base_name(st_["name"]) and np.hypot(*(q_ - base)) < 150
                                    for n_, q_ in capsuled):
            capsuled.add((base_name(st_["name"]), tuple(base)))
            capsule(a_, b_)
        zone = LineString([a_, b_]).buffer(cap_r + 5)
        main = st_["term"]                           # real termini only (B only visits Dworzec PKP, in and out)
        main = main - var_end.get(sid, set())        # a line ending here on a pale extension: diamond
        bullets = [(l, "o") for l in sorted(main)] + [(l, "D") for l in sorted(var_end.get(sid, set()))]
        if key_st:                                   # lines calling without ending there: squares
            bullets += [(l, "s") for l in sorted(st_["lines"] - main - var_end.get(sid, set()))
                        if st_["trips"].get(l, 0) >= EXT_SHARE * line_trips[l]]
        requests.append((0 if key_st else 1, -len(st_["lines"]), sid, base, zone, bullets))

    def bullet_w(bullets, trips):
        """Width of a row of bullets in points (a diamond carries its trips a day after it)."""
        return sum(15 + (4 + 5.2 * len(str(trips.get(l, 0))) if mk == "D" else 0) for l, mk in bullets)

    # termini off the map (Grabanów ODR, Czosnówka, ...): named by the edge where their line leaves the map, next to
    # a thick black arrow beside the line(s), parallel to them, pointing out (on the side clear of other lines;
    # bullets as on the map: diamonds with the trips a day for ends of extensions)
    mapbox = box(x0, y0, x1, y1)
    exits = []                                    # (crossing point, outward unit vector, line, name, trips there)
    for (l, n), sts in END_STATIONS.items():
        if sum(sts.values()) < EXT_MIN_ENDS or (l, n) in SHORT_WORKINGS or not paths.get(l):
            continue
        k = max(sts, key=sts.get)
        if mapbox.contains(Point(stations[k]["xy"])):
            continue
        # where the trips to that terminus really leave the map (A to Styrzyniec along Warszawska, not via
        # Porosiuki, which is nearer as the crow flies): the last edge crossing of its busiest shape ending there
        pts_ = lambda h: [q for q in (h.geoms if hasattr(h, "geoms") else [h])  # noqa: E731
                          if not q.is_empty and q.geom_type == "Point"]
        ref, ref_n = Point(stations[k]["xy"]), 0
        for g_, n_ in shapes.get(l, []):
            for at_end in (False, True):
                if n_ > ref_n and Point(g_.coords[-1 if at_end else 0]).distance(Point(stations[k]["xy"])) < 150:
                    cr = pts_(g_.intersection(mapbox.exterior))
                    if cr:
                        ref, ref_n = (max if at_end else min)(cr, key=g_.project), n_
        best = None
        for g in paths[l]:
            for q in pts_(g.intersection(mapbox.exterior)):
                dq = q.distance(ref)
                if best is None or dq < best[0]:
                    best = (dq, g, q)
        if best is None:
            continue
        _, g, q = best
        d = g.project(q)
        a, b = np.array(g.interpolate(max(d - 15, 0)).coords[0]), np.array(g.interpolate(min(d + 15, g.length)).coords[0])
        v = (b - a) / max(np.hypot(*(b - a)), 1e-9)
        if mapbox.contains(Point(b)):
            v = -v                                # pointing out of the map
        inward = -1 if not mapbox.contains(Point(b)) else 1        # along g, into the map
        exits.append((np.array(q.coords[0]), v, g, d, inward, l, n,
                      sum(stations[s_]["trips"].get(l, 0) for s_ in sts)))
    groups = []                                   # one label per name and exit: [name, exits]
    for ex in sorted(exits, key=lambda x: (x[6], x[5])):
        for gr in groups:
            if gr[0] == ex[6] and np.hypot(*(gr[1][0][0] - ex[0])) < 300:
                gr[1].append(ex)
                break
        else:
            groups.append([ex[6], [ex]])
    off_requests = []                             # (name, base, zone, bullets, trips)
    from shapely.ops import substring
    for n, exs in groups:
        q, v, g, d, inward, l0, _, _ = exs[0]
        # the arrow: beside the line(s), parallel, from 380 m to 100 m inside the edge, pointing out
        a_d, b_d = d + inward * 380, d + inward * 100
        if not (0 <= min(a_d, b_d) and max(a_d, b_d) <= g.length):
            continue
        piece = substring(g, a_d, b_d)            # inside -> edge
        mid = np.array(piece.interpolate(0.5, normalized=True).coords[0])
        pa, pb = np.array(piece.coords[0]), np.array(piece.coords[-1])
        tv = (pb - pa) / max(np.hypot(*(pb - pa)), 1e-9)
        nv = np.array([-tv[1], tv[0]])
        rel = [float((np.array(gg.interpolate(gg.project(Point(mid))).coords[0]) - mid) @ nv)
               for _, _, gg, _, _, _, _, _ in exs]
        best = None
        for sd_ in (1, -1):
            off = sd_ * (max(sd_ * r_ for r_ in rel) + SPACING * 0.5 + 55)
            arr = piece.offset_curve(off)
            if arr.is_empty or arr.geom_type != "LineString":
                continue
            ac = np.array(arr.coords)
            if np.hypot(*(ac[-1] - pa)) < np.hypot(*(ac[0] - pa)):
                ac = ac[::-1]                     # inside -> edge
            hit = sum(gg.distance(LineString(ac)) < SPACING for gg in drawn)
            if best is None or hit < best[0]:
                best = (hit, ac, sd_, off)
        if best is None:
            continue
        _, ac, sd_, off = best
        head_l, head_w = 85.0, 42.0
        L_ = LineString(ac)
        shaft = substring(L_, 0, L_.length - head_l * 0.6)
        ax.plot(*np.array(shaft.coords).T, color=INK, lw=lw * 0.9, solid_capstyle="butt", zorder=6.4)
        tip_ = ac[-1]
        hb = np.array(L_.interpolate(L_.length - head_l).coords[0])
        hv = (tip_ - hb) / max(np.hypot(*(tip_ - hb)), 1e-9)
        hn = np.array([-hv[1], hv[0]])
        ax.add_patch(plt.Polygon([tip_, hb + hn * head_w, hb - hn * head_w], closed=True, color=INK, zorder=6.4,
                                 linewidth=0))
        drawn.append(L_.buffer(head_w * 0.6))
        bullets = [(l, "o" if tr_ >= EXT_SHARE * line_trips[l] else "D") for _, _, _, _, _, l, _, tr_ in exs]
        trips = {l: tr_ for _, _, _, _, _, l, _, tr_ in exs}
        amid = np.array(L_.interpolate(0.4, normalized=True).coords[0])
        off_requests.append((n, amid + nv * sd_ * 60, L_.buffer(head_w + 10), bullets, trips))

    # label placement: as close to the stop as possible (16 directions, sliding outwards from the capsule),
    # off other names at any cost, off the lines where that costs little distance
    occ = [g.buffer(lw * m_per_pt / 2 + 4) for g in drawn]
    occ_tree = STRtree(occ)
    frame = box(x0 + 60, y0 + 60, x1 - 60, y1 - 60)  # names stay on the map (stops outside it get none)
    placed = []
    labelled = set()
    todo = [(base_name(stations[sid]["name"]), 14 if prio == 0 else 12.5, base, zone, bullets, stations[sid]["trips"])
            for prio, _, sid, base, zone, bullets in sorted(requests, key=lambda r: (r[0], r[1]))
            if mapbox.contains(Point(base))]           # stops off the map: named at the edge instead (below)
    todo += [(n, 12.5, base, zone, bullets, trips) for n, base, zone, bullets, trips in off_requests]
    for name, size, base, zone, bullets, trips in todo:
        if name in labelled:
            continue
        labelled.add(name)
        w_m = len(name) * size * 0.52 * m_per_pt
        h_m = size * 1.25 * m_per_pt
        if bullets:
            w_m = max(w_m, bullet_w(bullets, trips) * m_per_pt)
            h_m += 17 * m_per_pt
        best = None
        spur_name = U_INFO.get("spur") and name in KEY_STATIONS     # under the spur, off its lines
        for ang in (range(200, 341, 10) if spur_name else range(0, 360, 22)):
            v = np.array([math.cos(math.radians(ang)), math.sin(math.radians(ang))])
            first = None
            for dist in np.arange(0, 320 if spur_name else 200, 5.0):
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
        row_w = bullet_w(bullets, trips) * m_per_pt
        bx0 = {"left": tx, "right": tx - row_w, "center": tx - row_w / 2}[ha]
        xo = 0.0
        for l, mk in bullets:  # lines ending here: bullets under the name (diamond: variant, with its trips a day)
            bx = bx0 + (xo + 7) * m_per_pt
            by = ty - (size * 0.62 + 9) * m_per_pt
            ax.scatter([bx], [by], marker=mk, s=(13 if mk == "o" else 11.5 if mk == "D" else 11) ** 2, color=colors.get(l, "#888"),
                       zorder=7, linewidths=0)
            ax.text(bx, by, l, fontproperties=font(8), color="white", ha="center", va="center", zorder=7.1)
            xo += 15
            if mk == "D":
                n_ = trips.get(l, 0)
                ax.text(bx + 8 * m_per_pt, by, f"{n_}", fontproperties=font(8.5), color="#666", ha="left",
                        va="center", zorder=7.1, path_effects=halo)
                xo += 4 + 5.2 * len(str(n_))

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
    row = (y1 - y0) * 0.047
    top = y1 - (y1 - y0) * 0.22
    ax.add_patch(plt.Rectangle((px0, top - row * (len(order) - 0.4)), pw, row * len(order), facecolor="white",
                               edgecolor=INK, linewidth=1.2, zorder=8))
    xd, dgap = px0 + pw * 0.40, pw * 0.02            # the legend's dash column
    for k, l in enumerate(sorted(order)):
        y = top - k * row
        ax.scatter([px0 + pw * 0.06], [y], s=24 ** 2, color=colors.get(l, "#888"), zorder=9, linewidths=0)
        ax.text(px0 + pw * 0.06, y, l, fontproperties=font(14), color="white", ha="center", va="center", zorder=10)
        t1, t2, e1, e2 = legend_name(l, names.get(l, l))
        sub = dict(fontproperties=font(9.5), color="#777", va="center", zorder=9)
        has_sub = bool(e1 or e2)
        ym, ys = (y + row * 0.13, y - row * 0.22) if has_sub else (y, None)
        if t2 is None:                               # circular line: its streets underneath
            ax.text(px0 + pw * 0.13, ym, t1, fontproperties=font(15), color=INK, va="center", zorder=9)
            if e1:
                ax.text(px0 + pw * 0.13, ys, e1, ha="left", **sub)
            continue
        # termini in two columns round an aligned dash; variant ends under the terminus they extend beyond
        ax.text(xd - dgap, ym, t1, fontproperties=font(15), color=INK, va="center", ha="right", zorder=9)
        ax.text(xd, ym, "–", fontproperties=font(15), color=INK, va="center", ha="center", zorder=9)
        ax.text(xd + dgap, ym, t2, fontproperties=font(15), color=INK, va="center", ha="left", zorder=9)
        if e1:
            ax.text(xd - dgap, ys, e1, ha="right", **sub)
        if e2:
            ax.text(xd + dgap, ys, e2, ha="left", **sub)
    # key: stop, one-way stop, extension
    ky = top - row * (len(order) - 0.4) - row * 0.8
    kx0, kx1 = px0 + pw * 0.02, px0 + pw * 0.10
    key_rows = [("dot", "przystanek"), ("arrow", "przystanek w jednym kierunku jazdy"),
                ("dash", "wydłużenia: mniej niż połowa kursów linii"),
                ("diamond", "koniec wydłużenia, obok: kursów linii na dobę"),
                ("square", "linia przejeżdża przez węzeł (nie kończy biegu)")]
    kl = lw * 2.2                                   # key drawn larger than on the map, to be legible
    for k, (kind, text) in enumerate(key_rows):
        y = ky - k * row * 0.68
        xm = (kx0 + kx1) / 2
        if kind == "square":
            ax.scatter([xm - (kx1 - kx0) * 0.18], [y], marker="s", s=19 ** 2, color="#8a8a8a", zorder=9.5,
                       linewidths=0)
            ax.text(xm - (kx1 - kx0) * 0.18, y, "A", fontproperties=font(12), color="white", ha="center",
                    va="center", zorder=9.6)
        elif kind == "diamond":
            ax.scatter([xm - (kx1 - kx0) * 0.18], [y], marker="D", s=21 ** 2, color="#8a8a8a", zorder=9.5,
                       linewidths=0)
            ax.text(xm - (kx1 - kx0) * 0.18, y, "E", fontproperties=font(12), color="white", ha="center",
                    va="center", zorder=9.6)
            ax.text(xm - (kx1 - kx0) * 0.18 + 18 / pt_per_m, y, "4", fontproperties=font(13), color="#666", ha="left",
                    va="center", zorder=9.6)
        else:
            draw_line(ax, np.array([[kx0, y], [kx1, y]]), "#8a8a8a", kl, kind == "dash", zorder=9)
        if kind == "dot":
            ax.scatter([xm], [y], s=(kl * 0.7) ** 2, color=INK, zorder=9.5, linewidths=0)
        elif kind == "arrow":
            ax.scatter([xm], [y], marker=(3, 0, -90), s=(kl * 1.35) ** 2, color=INK, zorder=9.5, linewidths=0)
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
    ap.add_argument("--config", default=str(CONFIG))
    ap.add_argument("--ext-style", choices=["dashed", "thin", "pale", "hollow"], default=None)
    ap.add_argument("--station", choices=["u", "untangle", "under", "spur", "box"], default=None)
    args = ap.parse_args()
    global EXT_STYLE, STATION_STYLE, STATION_CROSS_W, ABOVE_CROSS_W, U_TOP_M, U_S_R, U_HALF_M
    EXT_STYLE = args.ext_style or EXT_STYLE
    STATION_STYLE = args.station or STATION_STYLE
    if STATION_STYLE in ("untangle", "under", "spur"):   # no crossings at the station: they go to the junction
        STATION_CROSS_W, ABOVE_CROSS_W = 5.0, 1.0
    if STATION_STYLE == "box":                           # a compact U, its bottom hidden in the station box
        U_TOP_M, U_S_R, U_HALF_M = 200, 120, 90
    geoms, n_shape, line_of, stations, line_trips = load(args.date)
    pieces = corridors(geoms, n_shape, line_of, line_trips)
    sets = lines_on(pieces, geoms, line_of)
    edges = straighten_roundabouts(build_edges(pieces, sets),
                                   [s["xy"] for s in stations.values() if s["name"] in KEY_STATIONS])
    edges = fold_parallel(contract_short(repair_lines(edges, geoms, line_of)))
    edges = drop_turning_stubs(edges, stations)
    edges = edge_trips(edges, geoms, n_shape, line_of)
    shapes = defaultdict(list)                   # line -> [(shape, trips of the day)]
    for sid, g in geoms.items():
        shapes[line_of[sid]].append((g, int(n_shape[sid])))
    base_orders = None
    if STATION_STYLE != "u":
        # line orders away from the station as with the U (the layout the user approved): the order search is
        # run on the U first, then only the sections near the station are searched again for this style
        style, cw = STATION_STYLE, (STATION_CROSS_W, ABOVE_CROSS_W)
        STATION_STYLE, (STATION_CROSS_W, ABOVE_CROSS_W) = "u", (0.3, 3.0)
        e_u = station_access(copy.deepcopy(edges), copy.deepcopy(stations))
        drop_duplicate_sections(e_u)
        set_joins(e_u, shapes)
        base_orders = best_orders(sorted(line_trips.index), e_u, args.seed)
        for d_ in (PAIRS, NODE_W, U_INFO, LANE_SHIFT, JOINS):
            d_.clear()
        STATION_STYLE, (STATION_CROSS_W, ABOVE_CROSS_W) = style, cw
    edges = station_access(edges, stations)
    hairpin_tips(edges, stations)
    drop_duplicate_sections(edges)
    print(f"corridor pieces {len(pieces)}, edges {len(edges)}")
    set_joins(edges, shapes)
    if base_orders is not None and len(base_orders) == len(edges) and             all(set(o) == set(e["lines"]) for o, e in zip(base_orders, edges)):
        near = {i for i, e in enumerate(edges) if e["geom"].distance(Point(U_INFO["node"])) < STATION_REORDER_M}
        orders = refine_orders(edges, base_orders, near, args.seed)
    else:
        orders = best_orders(sorted(line_trips.index), edges, args.seed)
    for i, e in enumerate(edges):       # E beside B through the Narutowicza bundle (it leaves with B; user 2026-11-17)
        for l, nb in ORDER_NEXT_TO:
            if len(e["lines"]) >= 8 and {l, nb} <= e["lines"]:
                o_ = [x for x in orders[i] if x != l]
                o_.insert(o_.index(nb) + 1, l)
                orders[i] = o_
    through_lines_east(edges, orders)
    arm_shifts(edges, orders)
    relax_junction_corners(edges, orders)
    rebuild_hugs(edges, orders)
    for k, tip in HUG_TIPS.items():                 # the hairpin stop on its drawn tip
        for s_ in stations.values():
            if s_["name"] in HAIRPINS.values() and edges[k]["lines"] & s_["lines"]:
                s_["xy"] = tip
    render(edges, orders, stations, json.load(open(args.config, encoding="utf-8")), args.date, args.out, shapes, line_trips)


if __name__ == "__main__":
    main()
