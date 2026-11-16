"""Animate every MZK bus of one day, moving along its shape (from data/gtfs/).

Positions are interpolated between stop times on shape_dist_traveled, then along the shape.
Style follows gtfs_schedules_city/animate (black page, dark roads, outlined dots with trails).
Basemap layers come from data/basemap/ (OSM gpkg layers, git-ignored); missing layers are skipped.

Output: video/biala_podlaska_buses_<date>.mp4 + video/frame_<HHMM>.png (preview stills)

Usage:  python animate_buses.py [--date 20261118] [--start 04:30] [--end 23:30] [--duration 120]
                                [--fps 30] [--still 07:30]   (--still: one PNG, no video)
"""
import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
from pyproj import Transformer

matplotlib.use("Agg")
import matplotlib.animation as animation  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

ROOT = Path(__file__).resolve().parent
GTFS = ROOT / "data" / "gtfs"
OUT = ROOT / "video"
BASEMAP = ROOT / "data" / "basemap"  # OSM gpkg layers (git-ignored)
CRS = "EPSG:2180"
FRAME = (773970, 466100, 790310, 476690)  # same as TransitFrequency's render_frame (stop extent + 300 m)
TRAIL_FRAMES = 6
DAY_NAMES = ["poniedziałek", "wtorek", "środa", "czwartek", "piątek", "sobota", "niedziela"]

BG = {  # layer: (style, kwargs)
    "water": dict(facecolor="#0b1a2b", edgecolor="none"),
    "forests": dict(facecolor="#07120a", edgecolor="none"),
    "parks": dict(facecolor="#07120a", edgecolor="none"),
    "buildings": dict(facecolor="#141418", edgecolor="none"),
    "railways": dict(color="#3a3a44", linewidth=0.8),
    "roads": dict(color="#2a2a30", linewidth=0.6),
}


def secs(t):
    h, m, s = t.str.split(":", expand=True).astype(int).T.values
    return h * 3600 + m * 60 + s


def load(date):
    st = pd.read_csv(GTFS / "stop_times.txt", dtype={"trip_id": str, "stop_id": str})
    trips = pd.read_csv(GTFS / "trips.txt", dtype=str)
    cal = pd.read_csv(GTFS / "calendar_dates.txt", dtype=str)
    routes = pd.read_csv(GTFS / "routes.txt", dtype=str)
    shapes = pd.read_csv(GTFS / "shapes.txt", dtype={"shape_id": str})

    services = set(cal[(cal.date == date) & (cal.exception_type == "1")].service_id)
    trips = trips[trips.service_id.isin(services)]
    if trips.empty:
        raise SystemExit(f"no service on {date}")
    st = st[st.trip_id.isin(trips.trip_id)].sort_values(["trip_id", "stop_sequence"])
    st["t"] = secs(st.departure_time)

    tr = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    shp = {}
    for sid, g in shapes.sort_values(["shape_id", "shape_pt_sequence"]).groupby("shape_id"):
        x, y = tr.transform(g.shape_pt_lon.values, g.shape_pt_lat.values)
        shp[sid] = (g.shape_dist_traveled.values, np.asarray(x), np.asarray(y))

    color = dict(zip(routes.route_id, "#" + routes.route_color))
    tinfo = trips.set_index("trip_id")
    runs = []  # (t_array, dist_array, shape_id, route_id)
    for tid, g in st.groupby("trip_id"):
        runs.append((g.t.values, g.shape_dist_traveled.values, tinfo.shape_id[tid], tinfo.route_id[tid]))
    return runs, shp, color, routes


def positions(runs, shp, t):
    """Positions of all buses on the road at second t: (xy array, route ids)."""
    xs, ys, rids = [], [], []
    for ts, ds, sid, rid in runs:
        if ts[0] <= t <= ts[-1]:
            d = np.interp(t, ts, ds)
            sd, sx, sy = shp[sid]
            xs.append(np.interp(d, sd, sx))
            ys.append(np.interp(d, sd, sy))
            rids.append(rid)
    return np.c_[xs, ys] if xs else np.empty((0, 2)), rids


def draw_base(ax, shp, routes):
    from shapely.geometry import box
    clip = box(*FRAME).buffer(500)
    for name, style in BG.items():
        f = BASEMAP / f"{name}.gpkg"
        if not f.exists():
            continue
        g = gpd.read_file(f, bbox=gpd.GeoSeries([clip], crs=CRS).to_crs(4326).total_bounds.tolist())
        if g.empty:
            continue
        g = g.to_crs(CRS)
        g = g[g.geometry.notna()]
        g.plot(ax=ax, zorder=1, aspect=None, **style)
    # network: faint line colours
    rc = dict(zip(routes.route_id, "#" + routes.route_color))
    segs, cols = [], []
    for sid, (_, x, y) in shp.items():
        segs.append(np.c_[x, y])
        cols.append(rc.get(sid.split("_")[0], "#777777"))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=1.2, alpha=0.22, zorder=2))


def setup(shp, routes, date):
    w, h = FRAME[2] - FRAME[0], FRAME[3] - FRAME[1]
    fig = plt.figure(figsize=(19.2, 19.2 * h / w), facecolor="#000000")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor("#000000")
    ax.set_xlim(FRAME[0], FRAME[2])
    ax.set_ylim(FRAME[1], FRAME[3])
    ax.set_aspect("equal")
    ax.set_axis_off()
    draw_base(ax, shp, routes)

    trail = LineCollection([], linewidths=3, alpha=0.35, capstyle="round", zorder=3)
    ax.add_collection(trail)
    dots = ax.scatter([], [], s=70, edgecolors="#ffffff", linewidths=1.2, zorder=4)
    d = pd.Timestamp(date)
    ax.text(0.015, 0.965, "BIAŁA PODLASKA", transform=ax.transAxes, color="#f0f0f5",
            fontsize=26, fontweight="bold", va="top", zorder=5)
    ax.text(0.015, 0.915, f"autobusy MZK · {DAY_NAMES[d.weekday()]} {d:%d.%m.%Y}",
            transform=ax.transAxes, color="#a0a0aa", fontsize=14, va="top", zorder=5)
    clock = ax.text(0.015, 0.87, "", transform=ax.transAxes, color="#ffffff", fontsize=40,
                    family="monospace", va="top", zorder=5)
    count = ax.text(0.015, 0.785, "", transform=ax.transAxes, color="#a0a0aa", fontsize=14,
                    va="top", zorder=5)
    # legend: line letters in their colours
    for k, (_, r) in enumerate(routes.sort_values("route_id").iterrows()):
        ax.text(0.985 - (len(routes) - 1 - k) * 0.026, 0.035, r.route_short_name, transform=ax.transAxes,
                color="#000000" if r.route_text_color == "000000" else "#ffffff", fontsize=12,
                fontweight="bold", ha="center", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.3", fc="#" + r.route_color, ec="none"))
    ax.text(0.985, 0.07, "Dane: rozkład MZK Biała Podlaska (GTFS nieoficjalny) · mapa © OpenStreetMap",
            transform=ax.transAxes, color="#55555f", fontsize=9, ha="right", zorder=5)
    return fig, dots, trail, clock, count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="20261118")
    ap.add_argument("--start", default="04:30")
    ap.add_argument("--end", default="23:30")
    ap.add_argument("--duration", type=float, default=120, help="video length, s")
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--still", help="HH:MM: render one PNG instead of the video")
    args = ap.parse_args()

    runs, shp, color, routes = load(args.date)
    fig, dots, trail, clock, count = setup(shp, routes, args.date)
    OUT.mkdir(exist_ok=True)
    hm = lambda s: int(s[:2]) * 3600 + int(s[3:5]) * 60  # noqa: E731

    t0, t1 = hm(args.start), hm(args.end)
    n = int(args.duration * args.fps)
    step = (t1 - t0) / n
    history = []

    def frame(i, t=None):
        t = t0 + i * step if t is None else t
        xy, rids = positions(runs, shp, t)
        # trails: positions of the same trips over the last frames
        tr_segs, tr_cols = [], []
        for k, (ts, ds, sid, rid) in enumerate(runs):
            if ts[0] <= t <= ts[-1]:
                tt = np.linspace(max(ts[0], t - TRAIL_FRAMES * step), t, TRAIL_FRAMES + 1)
                dd = np.interp(tt, ts, ds)
                sd, sx, sy = shp[sid]
                tr_segs.append(np.c_[np.interp(dd, sd, sx), np.interp(dd, sd, sy)])
                tr_cols.append(color[rid])
        trail.set_segments(tr_segs)
        trail.set_colors(tr_cols)
        dots.set_offsets(xy)
        dots.set_facecolors([color[r] for r in rids] if rids else [])
        clock.set_text(f"{int(t // 3600):02d}:{int(t % 3600 // 60):02d}")
        count.set_text(f"autobusów w trasie: {len(rids)}")
        return dots, trail, clock, count

    if args.still:
        frame(0, hm(args.still))
        out = OUT / f"frame_{args.still.replace(':', '')}.png"
        fig.savefig(out, dpi=100, facecolor="#000000")
        print(f"written {out}")
        return

    out = OUT / f"biala_podlaska_buses_{args.date}.mp4"
    writer = animation.FFMpegWriter(fps=args.fps, codec="libx264", bitrate=6000,
                                    extra_args=["-pix_fmt", "yuv420p"])
    with writer.saving(fig, str(out), dpi=100):
        for i in range(n):
            frame(i)
            writer.grab_frame(facecolor="#000000")
            if i % (args.fps * 10) == 0:
                print(f"  {i}/{n}  {clock.get_text()}", flush=True)
    print(f"written {out}")


if __name__ == "__main__":
    main()
