# GTFS for Biała Podlaska (MZK)

Unofficial GTFS feed built from the city bus timetable export of MZK Biała Podlaska
(https://mzkbp.pl/), which does not publish GTFS itself. The export is public and was obtained legally.

## Layout

| Path | What |
|---|---|
| `data/raw/1488.7z`, `data/raw/1488/` | Export as received (do not edit) |
| `convert_to_gtfs.py` | Raw export -> GTFS |
| `data/gtfs/`, `data/gtfs.zip` | Generated feed |
| `export_gis.py` | GTFS -> GeoPackage + network diagram |
| `data/gis/biala_podlaska_gtfs.gpkg` | EPSG:2180 layers: `stops` (posts, lines, departures per day type), `stations`, `variants` (one connected LineString per shape, trips per day type), `lines` (variants merged per line) |
| `png/network_diagram.png` | Quick network diagram (line colours, termini labels, trip counts in legend) |
| `data/qa/shape_legs.csv` | Per stop-to-stop leg: straight-line vs routed distance (shape QA) |
| `data/cache/` | OSM road graph + Overpass cache (git-ignored, re-downloaded if missing) |
| `data/validation/` | MobilityData validator report (git-ignored) |
| `tools/gtfs-validator-cli.jar` | MobilityData GTFS validator v8.0.1 (git-ignored) |

## Run

```bash
python convert_to_gtfs.py            # rebuild data/gtfs + data/gtfs.zip
python convert_to_gtfs.py --no-shapes   # skip OSM routing
python export_gis.py                 # GeoPackage + diagram (after the conversion)
java -jar tools/gtfs-validator-cli.jar -i data/gtfs.zip -o data/validation -c pl
```
Needs pandas, numpy, shapely, pyproj, osmnx (1.9), networkx; Java 17+ for the validator.

## Raw format (semicolon CSV, CP1250, CRLF)

- `przystanki.csv` — stop posts: `ID_PRZYST` (post id), `NR_PRZYST` (stop group), `NR_SLUPKA` (post no.),
  `NAZWA`, `POZ_X`/`POZ_Y` = WGS84 lon/lat × 100 000.
- `rozklad.csv` — one row per departure: stop, day type, line, variant, trip id (`ID_KRS`; `ID_KURSU` is a
  second unique id), `ODJAZD` HH:MM, trip note letters (`UWAGI_KURS`). **No stop sequence, no arrival times.**
- `warianty_tras.csv` — variants: direction `T`/`P`, `PODSTAWOWY` (basic variant), descriptions / headsign (`OPIS_TABL`).
- `typy_dni.csv` — day types RR weekday, SS Saturday, NN Sunday/holiday, XX school break.
- `kalendarz.csv` — day type per date (10.11.2026–31.12.2026).
- `uwagi.csv` — meaning of the trip note letters per line.

## Conversion choices

- **Stops**: each post is a stop (`stop_id` = `ID_PRZYST`, `stop_code` = `NR_PRZYST` + 2-digit post no.);
  each `NR_PRZYST` is a parent station `P<NR_PRZYST>` at the posts' centroid.
- **Service**: `service_id` = day type; `calendar_dates.txt` only (no `calendar.txt`).
- **Stop order**: all trips of a variant have identical running times, so stops are sorted by offset from the
  trip start. 33 groups of stops sharing a minute are ordered by the shortest path between neighbours.
- **Times**: arrival = departure = timetable minute, `timepoint=1`.
- **Shapes**: one per variant (`<line>_<variant>`), map-matched on the OSM `drive_service` network
  (Viterbi over the 6 nearest edges within 60 m of each stop); one-ways respected. Cost = length × road class
  weight (`ROAD_WEIGHTS`): primary/secondary/trunk 1.0, tertiary 1.05, residential/unclassified 1.1, service 1.5;
  `living_street` excluded. Mild on purpose: only removes side-street shortcuts (19 of 102 shapes changed, max +124 m).
  `shape_dist_traveled` in metres in shapes and stop_times. All stops lie within ~31 m of their shape.
- **Routes**: `route_type=3`, long name from the basic variant, colours assigned arbitrarily (`LINE_COLORS`).
- Trip notes (`UWAGI_KURS`) are not exported; headsigns come from the variant. Seasonal notes do not affect this
  period: "G" (to Grzybowa działki 1.IV–31.X) trips end at Grzybowa 1 in the export, "W" (not 1.XI) is before it.

## Status (2026-11-16)

- GeoPackage + `png/network_diagram.png` (`export_gis.py`); routing now uses road class weights, living streets excluded.

- Feed built: 13 lines (A–M, T weekends only), 318 posts / 172 stations, 1249 trips, 25 613 stop times, 102 shapes;
  valid 10.11–31.12.2026.
- MobilityData validator v8.0.1: **0 errors, 0 warnings**, 6 info `trip_headsign_matches_intermediate_stop`
  (loop trips, expected).
- `agency_url` / `feed_contact_url` = https://mzkbp.pl/ (operator site).
- Next: plug the feed into TransitFrequency (frequency / accessibility, new `--city`), subway_like_map /
  transit_scheme (network map) and gtfs_schedules_city (animation). Feed expires 31.12.2026 — analyses must
  use a date inside that window; a new export will be needed for 2027.
