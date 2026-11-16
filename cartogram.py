"""Flow cartogram: every street piece drawn as thick as the number of buses passing it in one day.

All shapes are noded into one network; each piece counts the trips of every shape within 1 m of
its midpoint (both directions together). Dual carriageways stay two lines.

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
from shapely.geometry import LineString, box
from shapely.ops import linemerge

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

ROOT = Path(__file__).resolve().parent
GTFS = ROOT / "data" / "gtfs"
BASEMAP = ROOT / "data" / "basemap"
GPKG = ROOT / "data" / "gis" / "biala_podlaska_gtfs.gpkg"
CRS = "EPSG:2180"
FRAME_ALL = (773970, 466100, 790310, 476690)   # all stops + 300 m (as animate_buses.py)
LIMITS = ROOT / "data" / "gis" / "city_limits.gpkg"  # OSM relation 3016433 (Nominatim)
MAX_WIDTH_PT = 16
DAY_NAMES = ["poniedziałek", "wtorek", "środa", "czwartek", "piątek", "sobota", "niedziela"]


def piece_counts(date):
    """Pieces of the noded union of all shapes, each with the trips of every shape running on it.
    A piece belongs to a shape when its midpoint lies within 1 m of it (shapes cut at their own
    stop positions still line up)."""
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
    noded = unary_union(geoms)
    pieces = []
    for piece in (noded.geoms if hasattr(noded, "geoms") else [noded]):
        mid = piece.interpolate(0.5, normalized=True)
        hits = [i for i in tree.query(mid.buffer(1.0)) if geoms[i].distance(mid) < 1.0]
        if hits:
            pieces.append((piece, sum(int(n_shape[sids[i]]) for i in hits),
                           ",".join(sorted({sids[i].split("_")[0] for i in hits}))))
    return pieces, len(day)


def merged_segments(pieces):
    """Join pieces with the same count and line set into longer LineStrings (for the gpkg)."""
    groups = defaultdict(list)
    for g, n, ls in pieces:
        groups[(n, ls)].append(g)
    rows = []
    for (n, ls), segs in groups.items():
        m = linemerge(segs)
        for g in (m.geoms if hasattr(m, "geoms") else [m]):
            rows.append(dict(trips=n, lines=ls, n_lines=len(ls.split(",")), geometry=g))
    return gpd.GeoDataFrame(rows, crs=CRS)


def draw(pieces, n_trips, date, frame, out, limits=None):
    w, h = frame[2] - frame[0], frame[3] - frame[1]
    fig = plt.figure(figsize=(19.2, 19.2 * h / w), facecolor="#000000")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor("#000000")
    ax.set_xlim(frame[0], frame[2])
    ax.set_ylim(frame[1], frame[3])
    ax.set_aspect("equal")
    ax.set_axis_off()

    clip = gpd.GeoSeries([box(*frame).buffer(500)], crs=CRS).to_crs(4326).total_bounds.tolist()
    for name, style in [("water", dict(facecolor="#0b1a2b", edgecolor="none")),
                        ("forests", dict(facecolor="#07120a", edgecolor="none")),
                        ("buildings", dict(facecolor="#17171c", edgecolor="none")),
                        ("roads", dict(color="#26262c", linewidth=0.5)),
                        ("railways", dict(color="#3a3a44", linewidth=0.8))]:
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
    vmax = n.max()
    cmap = plt.get_cmap("plasma")
    widths = 0.6 + (MAX_WIDTH_PT - 0.6) * n / vmax
    colors = cmap(0.25 + 0.75 * np.sqrt(n / vmax))
    ax.add_collection(LineCollection([list(p[0].coords) for p in pieces], linewidths=widths, colors=colors,
                                     capstyle="round", joinstyle="round", zorder=3))

    d = pd.Timestamp(date)
    ax.text(0.015, 0.965, "BIAŁA PODLASKA", transform=ax.transAxes, color="#f0f0f5",
            fontsize=26, fontweight="bold", va="top")
    ax.text(0.015, 0.915, f"Ile autobusów MZK przejeżdża ulicą w ciągu doby · {DAY_NAMES[d.weekday()]} "
            f"{d:%d.%m.%Y} · {n_trips} kursów", transform=ax.transAxes, color="#a0a0aa", fontsize=14, va="top")
    # width legend
    for i, v in enumerate([v for v in (10, 50, 100, 200, 300, 400) if v <= vmax * 1.05] + [int(vmax)]):
        y = 0.30 - i * 0.04
        lw = 0.6 + (MAX_WIDTH_PT - 0.6) * v / vmax
        ax.plot([0.02, 0.07], [y, y], transform=ax.transAxes, lw=lw, solid_capstyle="round",
                color=cmap(0.25 + 0.75 * np.sqrt(v / vmax)), zorder=5)
        ax.text(0.08, y, f"{v}" + (" (max)" if v == int(vmax) else ""), transform=ax.transAxes,
                color="#c0c0c8", fontsize=12, va="center", zorder=5)
    ax.text(0.02, 0.30 + 0.035, "autobusy / dobę (oba kierunki)", transform=ax.transAxes,
            color="#c0c0c8", fontsize=12, va="bottom", zorder=5)
    ax.text(0.985, 0.02, "Dane: rozkład MZK Biała Podlaska (GTFS nieoficjalny) · mapa © OpenStreetMap",
            transform=ax.transAxes, color="#55555f", fontsize=9, ha="right")
    fig.savefig(out, dpi=150, facecolor="#000000")
    plt.close(fig)
    print(f"written {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="20261118")
    ap.add_argument("--city", action="store_true", help="frame on the city limits (+400 m) instead of all stops")
    args = ap.parse_args()

    pieces, n_trips = piece_counts(args.date)
    seg = merged_segments(pieces)
    seg.to_file(GPKG, layer="segments", driver="GPKG")
    print(f"segments: {len(seg)} lines, busiest {seg.trips.max()} buses/day "
          f"({seg.sort_values('trips').iloc[-1].lines}) -> {GPKG.name}:segments")
    out = ROOT / "png" / f"cartogram_{args.date}{'_city' if args.city else ''}.png"
    if args.city:
        limits = gpd.read_file(LIMITS).to_crs(CRS)
        x0, y0, x1, y1 = limits.total_bounds
        draw(pieces, n_trips, args.date, (x0 - 400, y0 - 400, x1 + 400, y1 + 400), out, limits)
    else:
        draw(pieces, n_trips, args.date, FRAME_ALL, out)


if __name__ == "__main__":
    main()
