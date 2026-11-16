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
STEP_M = 6          # drawn lines are sampled this often
SIMPLIFY_M = 30     # sections are simplified this much before their corners are rounded (street wiggles go)
CORNER_R = 150      # m: corners become circular arcs of this radius (smaller where the legs are short)
RAMP_M = 120        # a line changing lane at a junction eases over this length (no jag)
STATION_DROP_M = 160  # Dworzec PKP: the through lines continue straight on this far past the station (schematic)
CONTRACT_M = 40     # junctions joined by shorter edges are merged into one
SHIFT_COST = 0.15   # order search: cost per lane a line shifts through a junction (vs 1 per crossing)
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


def contract_short(edges, min_len=None):
    """Junction clusters joined by edges shorter than CONTRACT_M m become one junction (at their centre): lines
    then change order once, with room for the lane ramps, instead of curling over a few metres of corridor."""
    import networkx as nx
    min_len = min_len or CONTRACT_M
    G = nx.Graph()
    G.add_edges_from((e["u"], e["v"]) for e in edges if e["geom"].length < min_len and e["u"] != e["v"])
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
    print(f"junction clusters merged: {len(set(rep.values()))} ({len(edges) - len(out)} short edges removed)")
    return out


def station_access(edges, stations, geoms, line_of):
    """Dworzec PKP drawn schematically: buses run down Al. 1000-lecia to the station, some turn back there, the
    rest go on east to the roundabout. The corridor builder sees one loop through the station carrying every line
    to the roundabout; it is split where it passes the station into the access leg (all lines, drawn straight
    through the station) and the through leg (only the lines whose routes use it), which continues straight on
    STATION_DROP_M m before turning, so the station sits in a straight bundle all lines pass through."""
    key = lambda c: (round(c[0], 1), round(c[1], 1))  # noqa: E731
    by_line = defaultdict(list)
    for sid, g in geoms.items():
        by_line[line_of[sid]].append(g)
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
        legs = [LineString([g.interpolate(x) for x in np.linspace(0, d, 40)]),
                LineString([g.interpolate(x) for x in np.linspace(d, g.length, 80)])]
        on = []
        for leg in legs:
            probe = [leg.interpolate(f, normalized=True) for f in (0.3, 0.5, 0.7)]
            on.append(frozenset(l for l in e["lines"]
                                if any(all(sh.distance(p) < TOL for p in probe) for sh in by_line[l])))
        a = 0 if len(on[0]) >= len(on[1]) else 1          # access leg: the one every line uses
        acc, thr = legs[a], legs[1 - a]
        if a == 1:
            acc = LineString(list(acc.coords)[::-1])        # both legs start at the far ends, meet at the station
        else:
            thr = LineString(list(thr.coords)[::-1])
        acc_c, thr_c = np.array(acc.coords), np.array(thr.coords)  # acc: far end -> station; thr: far end -> station
        J = acc_c[-1]
        dv = J - np.array(acc.interpolate(acc.length - min(150, acc.length / 2)).coords[0])
        dv /= np.hypot(*dv)
        # through leg: straight on past the station, then its corner (farthest point from the chord) moved as well
        chord = LineString([thr_c[0], thr_c[-1]])
        k = int(np.argmax([chord.distance(Point(p)) for p in thr_c]))
        new_thr = np.vstack([thr_c[:k + 1], [thr_c[k] + dv * STATION_DROP_M, J + dv * STATION_DROP_M, J]])
        for st in stations.values():                      # stops on the through leg follow it
            if st["name"] != name and Point(st["xy"]).distance(thr) < 40:
                st["xy"] = st["xy"] + dv * STATION_DROP_M * 0.5
        for k_ in sts:                                    # the station: on the access leg, just before the split
            stations[k_]["xy"] = J - dv * 1.0
        nodes_ = (e["u"], e["v"]) if a == 0 else (e["v"], e["u"])
        acc_g = LineString([nodes_[0]] + [tuple(p) for p in acc_c[1:-1]] + [tuple(J)])
        thr_g = LineString([nodes_[1]] + [tuple(p) for p in new_thr[1:]])
        all_lines = e["lines"]
        thr_lines = on[1 - a] or all_lines
        out[i] = dict(geom=acc_g, lines=all_lines, u=nodes_[0], v=key(J))
        out.append(dict(geom=thr_g, lines=thr_lines, u=nodes_[1], v=key(J)))
        print(f"{name}: access leg {sorted(all_lines)} {acc_g.length:.0f} m, "
              f"through leg {sorted(thr_lines)} {thr_g.length:.0f} m")
    return out


# ---------------------------------------------------------------- 4: line order (per edge)
def node_passes(edges):
    """For every node: incident edges with their heading away from the node and +1 / -1 (edge leaves / enters)."""
    nodes = defaultdict(list)
    for i, e in enumerate(edges):
        g = e["geom"]
        k = min(25, g.length / 2)
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


def node_cost(inc, orders, edges):
    """Crossings at one node (pairs of lines whose left/right order is wrong) + SHIFT_COST per lane shifted."""
    per_line = defaultdict(dict)
    for i, h, sgn in inc:
        for l in orders[i]:
            per_line[l][i] = (h, sgn)
    cost = 0.0
    ls = list(per_line)
    for a in ls:                                         # lane shifts of a line passing straight through
        ea = per_line[a]
        if len(ea) == 2:
            (p, (hp, sp)), (q, (hq, sq)) = ea.items()
            cost += SHIFT_COST * abs(-left_of(orders, p, a, sp) - left_of(orders, q, a, sq))
    # pairs only among lines passing through once (two edges here): at terminal loops a line touches the node
    # three times and every order "crosses" somewhere, which would drown the real junctions
    ls = [l for l in ls if len(per_line[l]) == 2]
    for x in range(len(ls)):
        for y in range(x + 1, len(ls)):
            a, b = ls[x], ls[y]
            ea, eb = per_line[a], per_line[b]
            for p in set(ea) & set(eb):
                hp, sp = ea[p]
                entry = hp + math.pi                     # direction of travel into the node along p
                la, lb = -left_of(orders, p, a, sp), -left_of(orders, p, b, sp)
                for qa, (ha, sa) in ea.items():
                    for qb, (hb, sb) in eb.items():
                        if qa == p or qb == p:
                            continue
                        if qa == qb:                     # on together: left/right must stay
                            if (la > lb) != (left_of(orders, qa, a, sa) > left_of(orders, qa, b, sb)):
                                cost += 1
                        else:                            # split: the line turning further left keeps left
                            ta, tb = turn(entry, ha), turn(entry, hb)
                            if abs(ta - tb) > 0.05 and (la > lb) != (ta > tb):
                                cost += 1
    return cost


def global_start(lines, edges, nodes, seed, iters=3000):
    """A good common order to start from (annealing over one order for all edges)."""
    rnd = random.Random(seed)
    def orders_for(order):
        rank = {l: k for k, l in enumerate(order)}
        return [sorted(e["lines"], key=rank.get) for e in edges]
    def total(order):
        o = orders_for(order)
        return sum(node_cost(inc, o, edges) for inc in nodes.values())
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
    cost_at = {n: node_cost(inc, orders, edges) for n, inc in nodes.items()}
    multi = [i for i, o in enumerate(orders) if len(o) > 1]
    rnd = random.Random(seed + 7)
    cur_total = sum(cost_at.values())
    best, best_total = [o[:] for o in orders], cur_total
    for it in range(iters):
        i = rnd.choice(multi)
        a_, b_ = rnd.sample(range(len(orders[i])), 2)
        orders[i][a_], orders[i][b_] = orders[i][b_], orders[i][a_]
        touched = set(ends[i])
        new = {n: node_cost(nodes[n], orders, edges) for n in touched}
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
            i = nxt[0] if len(inc[node]) == 2 and nxt else None
        return seq
    for n in [n for n, es in inc.items() if len(es) != 2]:
        for i in inc[n]:
            if i not in used:
                chains.append(walk(n, i))
    for i in mine:  # closed loops
        if i not in used:
            chains.append(walk(edges[i]["u"], i))
    return chains


_SECTIONS = {}


def draw_path(seq, l, edges, orders):
    """Drawn coordinates of one chain (subway_like_map style): each section simplified, corners (inside sections
    and at junctions) replaced by circular arcs of radius CORNER_R, the lane offset eased over RAMP_M m wherever
    it changes, then applied along the normals, so a bundle bends as concentric arcs."""
    pts, seg_off = [], []
    for i, fwd in seq:
        if i not in _SECTIONS:
            _SECTIONS[i] = section_coords(edges[i], len(orders[i]))
        c = _SECTIONS[i] if fwd else _SECTIONS[i][::-1]
        o = lane(orders, i, l) * (1 if fwd else -1)
        if pts:
            c = c[1:]
        pts.extend(c)
        seg_off.extend([o] * len(c))
    P = np.array(pts)
    seg_off = np.array(seg_off[1:] if len(seg_off) > 1 else seg_off, dtype=float)  # offset of segment k: P[k]-P[k+1]
    if len(P) < 2:
        return P
    F, sid = fillet(P, CORNER_R)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(F, axis=0), axis=1))])
    keep = np.concatenate([[True], np.diff(s) > 1e-6])
    F, s, sid = F[keep], s[keep], np.array(sid)[keep]
    if s[-1] < 1:
        return F
    ds = np.append(np.arange(0, s[-1], STEP_M), s[-1])
    Q = np.column_stack([np.interp(ds, s, F[:, 0]), np.interp(ds, s, F[:, 1])])
    O = seg_off[np.clip(sid[np.clip(np.searchsorted(s, ds, side="right") - 1, 0, len(sid) - 1)], 0, len(seg_off) - 1)]
    O = smooth(O[:, None].astype(float), RAMP_M)[:, 0]
    t = np.gradient(Q, axis=0)
    t = smooth(t, STEP_M * 3)
    t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-9)[:, None]
    nrm = np.stack([-t[:, 1], t[:, 0]], axis=1)
    return Q + nrm * O[:, None]


def station_spots(stations, edges):
    """Station -> (edge index, point on the corridor, unit tangent, lines of the station on that edge)."""
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
        spots[sid] = (i, np.array(g.interpolate(d).coords[0]), np.array([math.cos(h), math.sin(h)]),
                      s["lines"] & edges[i]["lines"])
    return spots


# ---------------------------------------------------------------- render
def render(edges, orders, stations, line_cfg, date, out):
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
    drawn = []                                    # drawn lines / capsules (map metres), kept clear of labels
    paths = defaultdict(list)                     # line -> drawn LineStrings
    line_ends = []                                # (line, end point, direction pointing out of the line)
    for l in order:
        for seq in chain_seqs(l, edges):
            c = draw_path(seq, l, edges, orders)
            if len(c) < 2:
                continue
            ax.plot(c[:, 0], c[:, 1], color=colors.get(l, "#888"), lw=lw, solid_capstyle="butt",
                    solid_joinstyle="round", zorder=3)
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
                    line_ends.append((l, c[k], d / max(np.hypot(*d), 1e-9)))

    # stops, interchanges, labels (one name per stop name without the post number)
    spots = station_spots(stations, edges)
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

    # line ends: at a main terminus of the line an end bar, at the end of an extended variant a white diamond
    ends_at = defaultdict(set)                    # station -> lines whose drawing ends there
    for l, xy, d in line_ends:
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
            ax.scatter([xy[0]], [xy[1]], marker="D", s=(lw * 1.45) ** 2, facecolor="white",
                       edgecolor=colors.get(l, "#888"), linewidths=lw * 0.38, zorder=5.2)
            drawn.append(Point(xy).buffer(SPACING * 0.6))

    # names: termini only (main termini of a line, and Dworzec PKP); every other stop is a dot on its lines
    requests = []                                 # labels to place
    cap_r = lw * 1.3 / 2 * m_per_pt
    for sid, (i, base, t, ls) in spots.items():
        st_ = stations[sid]
        nrm = np.array([-t[1], t[0]])
        key_st = st_["name"] in KEY_STATIONS
        if not (st_["term"] or key_st):
            for l in ls:
                p = on_line(l, base)
                if p is not None:
                    ax.scatter([p[0]], [p[1]], s=(lw * 0.7) ** 2, color=INK, zorder=4, linewidths=0)
            continue
        pts_ = [p for p in (on_line(l, base) for l in edges[i]["lines"]) if p is not None]
        if not pts_:
            continue
        side = [float(np.dot(p - base, nrm)) for p in pts_]
        lo, hi = min(side) - SPACING * 0.55, max(side) + SPACING * 0.55
        capsule(base + nrm * lo, base + nrm * hi)
        zone = LineString([base + nrm * lo, base + nrm * hi]).buffer(cap_r + 5)
        bullets = sorted(st_["term"] | (ends_at.get(sid, set()) if key_st else set()))
        requests.append((0 if key_st else 1, -len(st_["lines"]), sid, base, zone, bullets))

    # label placement: as close to the stop as possible (16 directions, sliding outwards from the capsule),
    # off other names at any cost, off the lines where that costs little distance
    occ = [g.buffer(lw * m_per_pt / 2 + 4) for g in drawn]
    occ_tree = STRtree(occ)
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
                if bb.intersects(zone):
                    continue
                first = dist if first is None else first
                hit = sum(occ[k].intersection(bb).area for k in occ_tree.query(bb))
                hit += sum(q.buffer(10).intersection(bb).area for q in placed) * 100
                score = hit + dist * 25
                if best is None or score < best[0]:
                    best = (score, bb, top, "left" if v[0] > 0.3 else "right" if v[0] < -0.3 else "center")
        _, bb, top, ha = best
        placed.append(bb)
        ty = top - size * 0.62 * m_per_pt
        # text aligned towards the stop: a name left of its stop ends at the stop's side of the box
        tx = {"left": bb.bounds[0], "right": bb.bounds[2], "center": (bb.bounds[0] + bb.bounds[2]) / 2}[ha]
        ax.text(tx, ty, name, fontproperties=font(size), color=INK, ha=ha, va="center", path_effects=halo, zorder=7)
        row_w = len(bullets) * 15 * m_per_pt
        bx0 = {"left": tx, "right": tx - row_w, "center": tx - row_w / 2}[ha]
        for k, l in enumerate(bullets):     # lines ending here: bullets under the name
            bx = bx0 + (k * 15 + 7) * m_per_pt
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
    edges = contract_short(repair_lines(build_edges(pieces, sets), geoms, line_of))
    edges = station_access(edges, stations, geoms, line_of)
    print(f"corridor pieces {len(pieces)}, edges {len(edges)}")
    orders = best_orders(sorted(line_trips.index), edges, args.seed)
    render(edges, orders, stations, json.load(open(CONFIG, encoding="utf-8")), args.date, args.out)


if __name__ == "__main__":
    main()
