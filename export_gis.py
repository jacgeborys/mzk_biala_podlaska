"""Export the GTFS feed in data/gtfs/ to a GeoPackage and render a network diagram.

Output:
  data/gis/biala_podlaska_gtfs.gpkg (EPSG:2180)
    stops     - stop posts, with lines served and departures per day type
    stations  - parent stations (NR_PRZYST groups)
    variants  - one connected LineString per shape (line variant), with trip counts
    lines     - one merged geometry per line (all its variants)
  png/network_diagram.png

Usage:  python export_gis.py
"""
from pathlib import Path

import geopandas as gpd
import matplotlib
import pandas as pd
from shapely.geometry import LineString
from shapely.ops import linemerge, unary_union

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent
GTFS = ROOT / "data" / "gtfs"
GPKG = ROOT / "data" / "gis" / "biala_podlaska_gtfs.gpkg"
PNG = ROOT / "png" / "network_diagram.png"
CRS = "EPSG:2180"
DAY_TYPES = ["RR", "SS", "NN", "XX"]


def read(name):
    return pd.read_csv(GTFS / name, dtype=str, keep_default_na=False)


def main():
    stops, routes, trips = read("stops.txt"), read("routes.txt"), read("trips.txt")
    st, shapes = read("stop_times.txt"), read("shapes.txt")
    for c in ("shape_pt_lat", "shape_pt_lon"):
        shapes[c] = shapes[c].astype(float)
    shapes["shape_pt_sequence"] = shapes.shape_pt_sequence.astype(int)

    # trip counts per shape and day type
    counts = (trips.groupby(["shape_id", "service_id"]).size().unstack(fill_value=0)
              .reindex(columns=DAY_TYPES, fill_value=0).add_prefix("trips_"))

    # variants: one connected line per shape
    rows = []
    for sid, g in shapes.sort_values(["shape_id", "shape_pt_sequence"]).groupby("shape_id"):
        t = trips[trips.shape_id == sid].iloc[0]
        rows.append(dict(shape_id=sid, route_id=t.route_id, variant=sid.split("_", 1)[1],
                         headsign=t.trip_headsign, direction_id=int(t.direction_id),
                         geometry=LineString(zip(g.shape_pt_lon, g.shape_pt_lat))))
    variants = gpd.GeoDataFrame(rows, crs="EPSG:4326").to_crs(CRS)
    variants = variants.join(counts, on="shape_id")
    variants = variants.merge(routes[["route_id", "route_long_name", "route_color"]], on="route_id")
    variants["route_color"] = "#" + variants.route_color
    variants["length_m"] = variants.length.round(1)

    # lines: all variants of a line merged into one (multi)linestring
    lines = []
    for rid, g in variants.groupby("route_id"):
        lines.append(dict(route_id=rid, route_long_name=g.route_long_name.iloc[0],
                          route_color=g.route_color.iloc[0], n_variants=len(g),
                          **{f"trips_{d}": int(g[f"trips_{d}"].sum()) for d in DAY_TYPES},
                          geometry=linemerge(unary_union(list(g.geometry)))))
    lines = gpd.GeoDataFrame(lines, crs=CRS)

    # stops with lines served and departures per day type
    stm = st.merge(trips[["trip_id", "route_id", "service_id"]], on="trip_id")
    dep = (stm.groupby(["stop_id", "service_id"]).size().unstack(fill_value=0)
           .reindex(columns=DAY_TYPES, fill_value=0).add_prefix("dep_"))
    served = stm.groupby("stop_id").route_id.agg(lambda s: ",".join(sorted(set(s)))).rename("lines")
    pts = gpd.GeoDataFrame(stops, geometry=gpd.points_from_xy(stops.stop_lon.astype(float),
                                                             stops.stop_lat.astype(float)),
                           crs="EPSG:4326").to_crs(CRS)
    posts = pts[pts.location_type == "0"].join(served, on="stop_id").join(dep, on="stop_id")
    stations = pts[pts.location_type == "1"].drop(columns=["location_type", "parent_station"])
    child = posts.groupby("parent_station")
    stations = stations.join(child.lines.agg(lambda s: ",".join(sorted(set(",".join(s).split(","))))),
                             on="stop_id")
    stations = stations.join(child[[f"dep_{d}" for d in DAY_TYPES]].sum(), on="stop_id")

    GPKG.parent.mkdir(parents=True, exist_ok=True)
    # layers are overwritten one by one (the file may be open in QGIS; keeps cartogram.py's `segments`)
    for name, gdf in [("stops", posts), ("stations", stations), ("variants", variants), ("lines", lines)]:
        gdf.to_file(GPKG, layer=name, driver="GPKG")
        print(f"  {name:9s} {len(gdf):4d} features")
    print(f"written {GPKG}")

    render(lines, variants, stations)


def render(lines, variants, stations):
    fig, ax = plt.subplots(figsize=(14, 12))
    # busier lines underneath, so rare lines stay visible on shared streets
    for _, r in lines.sort_values("trips_RR", ascending=False).iterrows():
        gpd.GeoSeries([r.geometry]).plot(ax=ax, color=r.route_color, linewidth=3, alpha=0.85,
                                         capstyle="round", joinstyle="round")
    stations.plot(ax=ax, color="white", edgecolor="black", markersize=14, linewidth=0.7, zorder=5)
    # label termini of the weekday network: stations at line ends
    ends = set()
    for g in variants.geometry:
        ends.update([g.coords[0], g.coords[-1]])
    from shapely.geometry import Point
    endpts = gpd.GeoSeries([Point(p) for p in ends], crs=CRS)
    near = stations[stations.geometry.apply(lambda s: endpts.distance(s).min() < 80)]
    for _, s in near.iterrows():
        ax.annotate(s.stop_name, (s.geometry.x, s.geometry.y), xytext=(4, 4), textcoords="offset points",
                    fontsize=7, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    handles = [plt.Line2D([], [], color=r.route_color, lw=4,
                          label=f"{r.route_id}  {r.route_long_name}\n     kursy: rob. {r.trips_RR}"
                                f" · sob. {r.trips_SS} · niedz. {r.trips_NN} · ferie {r.trips_XX}")
               for _, r in lines.sort_values("route_id").iterrows()]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, frameon=False,
              title="Linie MZK", labelspacing=0.9)
    ax.set_title("Biała Podlaska — sieć komunikacji miejskiej (MZK), GTFS 10.11–31.12.2026",
                 fontsize=13, loc="left")
    ax.set_axis_off()
    PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"written {PNG}")


if __name__ == "__main__":
    main()
