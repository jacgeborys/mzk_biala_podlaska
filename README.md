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
| `animate_buses.py` | Animation of every bus of one day moving along its shape (dots + trails, clock, bus count) |
| `video/biala_podlaska_buses_<date>.mp4`, `video/frame_<HHMM>.png` | Animation (mp4 git-ignored) and preview stills |
| `data/basemap/` | OSM basemap layers (gpkg, EPSG:4326; roads, water, buildings, forests, railways, ...), git-ignored; fetched 2026-11-16 from the Geofabrik Lubelskie extract, bbox 51.990-52.083 N, 22.986-23.245 E (no script here regenerates it yet) |
| `data/qa/shape_legs.csv` | Per stop-to-stop leg: straight-line vs routed distance (shape QA) |
| `data/cache/` | OSM road graph + Overpass cache (git-ignored, re-downloaded if missing) |
| `data/validation/` | MobilityData validator report (git-ignored) |
| `tools/gtfs-validator-cli.jar` | MobilityData GTFS validator v8.0.1 (git-ignored) |

## Run

```bash
python convert_to_gtfs.py            # rebuild data/gtfs + data/gtfs.zip
python convert_to_gtfs.py --no-shapes   # skip OSM routing
python export_gis.py                 # GeoPackage + diagram (after the conversion)
python animate_buses.py --still 07:30   # one preview frame
python animate_buses.py              # video: 18.11.2026 04:30-23:30 in 120 s, 30 fps (~3600 frames)
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
  (Viterbi over the 12 nearest edges within 60 m of each stop; 6 was too few: at Białka 2 the right carriageway
  was 8th, behind a frontage road and both directions of side streets); one-ways respected. Cost = length × road class
  weight (`ROAD_WEIGHTS`): primary/secondary/trunk 1.0, tertiary 1.05, unclassified 1.1, residential 1.6,
  service 2.0; `living_street` excluded. Keeps buses on tertiary-or-above roads unless a stop is on a residential
  street (residential 1.1 -> 1.6 changed 69 of 102 shapes, typically +85-145 m; checked: C_j no longer U-turns on
  Jana III Sobieskiego, G_h Brzeska -> Narutowicza now via Al. 1000-lecia instead of a residential street).
  **Turn penalties**: routing runs on an edge-based graph (node = road edge, arc = edge-to-edge move), each move
  costs its turn angle: < 20° free, else 60 m × angle/90° (~8-10 s of bus time), >= 150° (U-turn) 200 m. Edge
  heading is measured over its first/last 15 m so OSM micro-jogs do not count as turns. Fixes zigzags through
  equal-class side streets: A/F/G/H Sidorska 5 -> Orzechowa 1 now straight up Włoska and left into Orzechowa
  (was Piaskowa -> Pływacka, 64 m shorter but 3 turns; 20 m per turn was not enough). Also removed U-turn loops
  at Dworzec PKP, Terebelska Pętla, Terebela, Al. Solidarności (most legs 30-150 m shorter).
  **Known bus ways** (`BUS_WAYS`, weight 1.0 whatever their class): OSM has almost no bus data for Biała Podlaska
  (no usable route relations, Overpass not relied on), so ways buses use against their class are listed by hand.
  So far: way 123907910, the unnamed frontage road north of the Białka trunk; E / K / T leave Białka at
  Powstańców and serve Białka 2 / Białka 3 there in both directions (user, Street View 2026-11-16).
  `shape_dist_traveled` in metres in shapes and stop_times. All stops lie within ~31 m of their shape.
- **Routes**: `route_type=3`, long name from the basic variant, colours assigned arbitrarily (`LINE_COLORS`).
- Trip notes (`UWAGI_KURS`) are not exported; headsigns come from the variant. Seasonal notes do not affect this
  period: "G" (to Grzybowa działki 1.IV–31.X) trips end at Grzybowa 1 in the export, "W" (not 1.XI) is before it.

## Status (2026-11-16)

- GeoPackage + `png/network_diagram.png` (`export_gis.py`); routing uses road class weights (residential 1.6),
  living streets excluded, turn penalties (Włoska / Orzechowa zigzag fixed), Białka frontage road via `BUS_WAYS`
  (both verified by the user).
- Feed built: 13 lines (A–M, T weekends only), 318 posts / 172 stations, 1249 trips, 25 613 stop times, 102 shapes;
  valid 10.11–31.12.2026.
- MobilityData validator v8.0.1: **0 errors, 0 warnings**, 6 info `trip_headsign_matches_intermediate_stop`
  (loop trips, expected).
- `agency_url` / `feed_contact_url` = https://mzkbp.pl/ (operator site).
- Animation: `animate_buses.py` (own implementation; positions interpolated on `shape_dist_traveled`, style after
  gtfs_schedules_city: black page, dark basemap, outlined dots with trails). Wednesday 18.11.2026: 474 trips, 23 buses
  on the road at 07:30.
- **Scope rule**: everything stays in this repo. Sibling projects (TransitFrequency, subway_like_map,
  gtfs_schedules_city, osm_basemap) are references only; an attempt to add the city there was reverted (2026-11-16).
- Next: frequency / accessibility analysis and a network map, implemented here. Feed expires 31.12.2026 — analyses
  must use a date inside that window; a new export will be needed for 2027.
