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
import pandas as pd
from pyproj import Transformer
from shapely.geometry import LineString, Point, box
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
ROUND_M = 12        # corner rounding (densify step for Chaikin smoothing)
KEY_STATIONS = {"Dworzec PKP"}   # always drawn as interchange with a name
MIN_INTERCHANGE = 5               # stations served by this many lines get a capsule and a name
MAIN_TERMINUS = 0.5               # a station is a line's terminus when this share of its trips start or end there
CONNECT_M = 500     # a line broken into pieces is joined along the corridors if the gap is at most this long
DROP_M = 150        # pieces of a line shorter than this that cannot be joined are dropped (strays)
END_JOIN_M = 400    # a loose end of a line inside the network continues this far at most to reach the line again
PAIR_M = 80         # two loose ends of one line this close are bridged
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
    edges = {}
    for i, (p, s) in enumerate(zip(pieces, sets)):
        if s:
            edges[i] = dict(geom=p, lines=frozenset(s), u=key(p.coords[0]), v=key(p.coords[-1]))
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


# ---------------------------------------------------------------- 4: line order
def node_passes(edges):
    """For every node: incident edges (with direction away from the node) and, per line, the edges it uses."""
    nodes = defaultdict(list)
    for i, e in enumerate(edges):
        g = e["geom"]
        a = heading_at(g, min(15, g.length / 2))
        b = heading_at(g, max(0, g.length - min(15, g.length / 2))) + math.pi
        nodes[e["u"]].append((i, a, +1))   # +1: the edge leaves the node in its own direction
        nodes[e["v"]].append((i, b, -1))
    return nodes


def crossings(order, edges, nodes):
    rank = {l: k for k, l in enumerate(order)}

    def lat(i, l, sgn):  # lateral position (left positive) looking away from the node along edge i
        ls = sorted(edges[i]["lines"], key=rank.get)
        return sgn * ((len(ls) - 1) / 2 - ls.index(l))
    cost = 0
    for n, inc in nodes.items():
        if len(inc) < 2:
            continue
        per_line = defaultdict(list)
        for i, h, sgn in inc:
            for l in edges[i]["lines"]:
                per_line[l].append((i, h, sgn))
        ls = list(per_line)
        for x in range(len(ls)):
            for y in range(x + 1, len(ls)):
                a, b = ls[x], ls[y]
                ea = {i: (h, s) for i, h, s in per_line[a]}
                eb = {i: (h, s) for i, h, s in per_line[b]}
                for p in set(ea) & set(eb):
                    hp, sp = ea[p]
                    la_in, lb_in = -lat(p, a, sp), -lat(p, b, sp)   # looking into the node along p
                    for qa in ea:
                        for qb in eb:
                            if qa == p or qb == p:
                                continue
                            if qa == qb:     # continue together: lateral order must stay
                                hq, sq = ea[qa]
                                if (la_in > lb_in) != (lat(qa, a, sq) > lat(qa, b, sq)):
                                    cost += 1
                            else:            # split: the line turning more to the left should be on the left
                                ta = (ea[qa][0] - (hp + math.pi)) % (2 * math.pi)
                                tb = (eb[qb][0] - (hp + math.pi)) % (2 * math.pi)
                                if abs(ta - tb) > 1e-3 and (la_in > lb_in) != (ta > tb):
                                    cost += 1
    return cost


def best_order(lines, edges, seed, iters=5000, restarts=2):
    nodes = node_passes(edges)
    best, best_c = None, None
    for r in range(restarts):
        rnd = random.Random(seed * 100 + r)
        cur = sorted(lines)
        rnd.shuffle(cur)
        cur_c = crossings(cur, edges, nodes)
        T0 = 3.0
        for it in range(iters):
            cand = cur[:]
            if rnd.random() < 0.5:                       # swap two lines
                i, j = rnd.sample(range(len(cand)), 2)
                cand[i], cand[j] = cand[j], cand[i]
            else:                                         # move one line elsewhere
                i, j = rnd.sample(range(len(cand)), 2)
                cand.insert(j, cand.pop(i))
            c = crossings(cand, edges, nodes)
            T = T0 * (1 - it / iters) + 1e-3
            if c <= cur_c or rnd.random() < math.exp((cur_c - c) / T):
                cur, cur_c = cand, c
                if best_c is None or c < best_c:
                    best, best_c = cand[:], c
    return best, best_c


# ---------------------------------------------------------------- 5: drawing geometry
def lane_offset(e, l, rank):
    ls = sorted(e["lines"], key=rank.get)
    return ((len(ls) - 1) / 2 - ls.index(l)) * SPACING


def offset_coords(g, d):
    if abs(d) < 1e-6:
        return list(g.coords)
    o = g.offset_curve(d, join_style="round", quad_segs=4)
    ps = parts(o)
    if not ps:
        return list(g.coords)
    o = max(ps, key=lambda p: p.length)
    c = list(o.coords)
    # keep the direction of g
    if Point(c[0]).distance(Point(g.coords[0])) > Point(c[-1]).distance(Point(g.coords[0])):
        c = c[::-1]
    return c


def chaikin(coords, step=ROUND_M, rounds=3):
    g = LineString(coords)
    if g.length < 2 * step:
        return coords
    pts = np.array([g.interpolate(t).coords[0] for t in np.arange(0, g.length, step)] + [coords[-1]])
    for _ in range(rounds):
        q = 0.75 * pts[:-1] + 0.25 * pts[1:]
        r = 0.25 * pts[:-1] + 0.75 * pts[1:]
        mid = np.empty((2 * len(q), 2))
        mid[0::2], mid[1::2] = q, r
        pts = np.vstack([pts[:1], mid, pts[-1:]])
    return pts


def line_chains(l, edges, rank):
    """Continuous drawn paths of line l: its edges chained through nodes where it has exactly two edges."""
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
            i = nxt[0] if len(inc[node]) == 2 and nxt else None
        return seq
    starts = [n for n, es in inc.items() if len(es) != 2]
    for n in starts:
        for i in inc[n]:
            if i not in used:
                chains.append(walk(n, i))
    for i in mine:  # closed loops
        if i not in used:
            chains.append(walk(edges[i]["u"], i))
    out = []
    for seq in chains:
        coords = []
        for i, fwd in seq:
            c = offset_coords(edges[i]["geom"], lane_offset(edges[i], l, rank))
            coords += c if fwd else c[::-1]
        if len(coords) >= 2:
            out.append(chaikin(coords))
    return out


def station_spots(stations, edges, rank):
    """Station -> (edge index, point on the corridor, unit tangent, lines drawn there)."""
    geoms = [e["geom"] for e in edges]
    tree = STRtree(geoms)
    spots = {}
    for sid, s in stations.items():
        pt = Point(s["xy"])
        cand = [i for i in tree.query(pt.buffer(90)) if edges[i]["lines"] & s["lines"]]
        if not cand:
            continue
        i = min(cand, key=lambda i: geoms[i].distance(pt))
        g = geoms[i]
        d = g.project(pt)
        h = heading_at(g, d)
        base = np.array(g.interpolate(d).coords[0])
        spots[sid] = (i, base, np.array([math.cos(h), math.sin(h)]), s["lines"] & edges[i]["lines"])
    return spots


# ---------------------------------------------------------------- render
def render(edges, order, stations, line_cfg, date, out):
    rank = {l: k for k, l in enumerate(order)}
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
    drawn = []                                    # drawn lines / capsules (map metres), kept clear of labels
    for l in reversed(order):
        for c in line_chains(l, edges, rank):
            c = np.asarray(c)
            ax.plot(c[:, 0], c[:, 1], color=colors.get(l, "#888"), lw=lw, solid_capstyle="butt",
                    solid_joinstyle="round", zorder=3)
            drawn.append(LineString(c))

    # where each line ends: a terminus bar across its lane, and the stop there gets a name
    inc = defaultdict(list)
    for i, e in enumerate(edges):
        for l in e["lines"]:
            inc[(l, e["u"])].append(i)
            inc[(l, e["v"])].append(i)
    line_ends = [(l, np.array(n), es[0]) for (l, n), es in inc.items() if len(es) == 1]

    # stops, interchanges, labels (one name per stop name without the post number)
    spots = station_spots(stations, edges, rank)
    halo = [patheffects.withStroke(linewidth=4, foreground=LAND)]
    m_per_pt = 1 / pt_per_m
    base_name = lambda n: re.sub(r"\s+\d+[A-Za-z]?(\s*\(.*\))?$", "", n).strip()  # noqa: E731

    def capsule(p, q):
        ax.plot([p[0], q[0]], [p[1], q[1]], color=INK, lw=lw * 1.3, solid_capstyle="round", zorder=5)
        ax.plot([p[0], q[0]], [p[1], q[1]], color="white", lw=lw * 0.85, solid_capstyle="round", zorder=5.1)
        drawn.append(LineString([p, q]))

    # stations at line ends (within 200 m of the end, served by that line) count as termini of that line
    ends_at = defaultdict(set)
    for l, xy, i in line_ends:
        cand = [k for k, s_ in stations.items() if l in s_["lines"] and np.hypot(*(s_["xy"] - xy)) < 200]
        if cand:
            ends_at[min(cand, key=lambda k: np.hypot(*(stations[k]["xy"] - xy)))].add(l)
        g = edges[i]["geom"]
        at_start = np.hypot(*(np.array(g.coords[0]) - xy)) < 1
        h = heading_at(g, 0 if at_start else g.length)
        nrm = np.array([-math.sin(h), math.cos(h)])
        p = xy + nrm * lane_offset(edges[i], l, rank)
        capsule(p - nrm * SPACING * 0.5, p + nrm * SPACING * 0.5)

    requests = []                                 # labels to place
    for sid, (i, base, t, ls) in spots.items():
        s = stations[sid]
        nrm = np.array([-t[1], t[0]])
        term = s["term"] | ends_at.get(sid, set())
        big = s["name"] in KEY_STATIONS or len(s["lines"]) >= MIN_INTERCHANGE or term
        if not big:
            for o in [lane_offset(edges[i], l, rank) for l in ls]:
                p = base + nrm * o
                ax.scatter([p[0]], [p[1]], s=(lw * 0.7) ** 2, color=INK, zorder=4, linewidths=0)
            continue
        all_offs = [lane_offset(edges[i], l, rank) for l in edges[i]["lines"]]
        lo, hi = min(all_offs) - SPACING * 0.55, max(all_offs) + SPACING * 0.55
        capsule(base + nrm * lo, base + nrm * hi)
        prio = 0 if s["name"] in KEY_STATIONS else 1 if len(s["lines"]) >= MIN_INTERCHANGE else 2
        requests.append((prio, -len(s["lines"]), sid, base, max(abs(lo), abs(hi)), sorted(term, key=rank.get)))

    # label placement: 16 candidates round the stop, scored by overlap with lines, capsules and placed labels
    occ = [g.buffer(lw * m_per_pt / 2 + 4) for g in drawn]
    occ_tree = STRtree(occ)
    placed = []
    labelled = set()
    for prio, _, sid, base, half, bullets in sorted(requests, key=lambda r: (r[0], r[1])):
        name = base_name(stations[sid]["name"])
        if name in labelled:
            continue
        labelled.add(name)
        size = 14 if prio == 0 else 13 if prio == 1 else 11.5
        w_m = len(name) * size * 0.52 * m_per_pt
        h_m = size * 1.25 * m_per_pt
        if bullets:
            w_m = max(w_m, len(bullets) * 15 * m_per_pt)
            h_m += 17 * m_per_pt
        best = None
        for dist in (half + SPACING * 1.3, half + SPACING * 3.0):
            for ang in range(0, 360, 45):
                v = np.array([math.cos(math.radians(ang)), math.sin(math.radians(ang))])
                anchor = base + v * dist
                bx0 = anchor[0] if v[0] > 0.3 else anchor[0] - w_m if v[0] < -0.3 else anchor[0] - w_m / 2
                top = anchor[1] + (h_m if v[1] > 0.3 else h_m / 2 if abs(v[1]) <= 0.3 else 0)
                bb = box(bx0, top - h_m, bx0 + w_m, top)
                hit = sum(occ[k].intersection(bb).area for k in occ_tree.query(bb))
                hit += sum(q.intersection(bb).area for q in placed) * 10
                score = hit + dist * 0.01
                if best is None or score < best[0]:
                    best = (score, bb, top)
        _, bb, top = best
        placed.append(bb)
        ty = top - size * 0.62 * m_per_pt
        ax.text(bb.bounds[0], ty, name, fontproperties=font(size), color=INK, ha="left", va="center",
                path_effects=halo, zorder=7)
        for k, l in enumerate(bullets):     # lines ending here: bullets under the name
            bx = bb.bounds[0] + (k * 15 + 7) * m_per_pt
            by = ty - (size * 0.62 + 9) * m_per_pt
            ax.scatter([bx], [by], s=13 ** 2, color=colors.get(l, "#888"), zorder=7, linewidths=0)
            ax.text(bx, by, l, fontproperties=font(8), color="white", ha="center", va="center", zorder=7.1)

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
    args = ap.parse_args()
    geoms, n_shape, line_of, stations, line_trips = load(args.date)
    pieces = corridors(geoms, n_shape, line_of, line_trips)
    sets = lines_on(pieces, geoms, line_of)
    edges = repair_lines(build_edges(pieces, sets), geoms, line_of)
    print(f"corridor pieces {len(pieces)}, edges {len(edges)}")
    order, n_cross = best_order(sorted(line_trips.index), edges, args.seed)
    print(f"line order {''.join(order)}, crossings {n_cross}")
    render(edges, order, stations, json.load(open(CONFIG, encoding="utf-8")), args.date, args.out)


if __name__ == "__main__":
    main()
