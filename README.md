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
| `data/gis/biala_podlaska_gtfs.gpkg` | EPSG:2180 layers: `stops` (posts, lines, departures per day type), `stations`, `variants` (one connected LineString per shape, trips per day type), `lines` (variants merged per line), `segments` + `segment_nodes` (simplified bus network from `cartogram.py`: edges with trips / lines / from_node / to_node, nodes with degree and busiest edge; reusable by other diagrams) |
| `png/network_diagram.png` | Quick preview, only with `python export_gis.py --png` |
| `animate_buses.py` | Animation of every bus of one day moving along its shape (dots + trails, clock, bus count) |
| `video/biala_podlaska_buses_<date>.mp4`, `video/frame_<HHMM>.png` | Animation (mp4 git-ignored) and preview stills |
| `data/basemap/` | OSM basemap layers (gpkg, EPSG:4326; roads, water, buildings, forests, railways, ...), git-ignored; fetched 2026-11-16 from the Geofabrik Lubelskie extract, bbox 51.990-52.083 N, 22.986-23.245 E (no script here regenerates it yet) |
| `cartogram.py` | Flow cartogram on a simplified node/edge network: shapes folded into one centreline per street (35 m: directions, carriageways, variants; a loop / out-and-back shape keeps one copy of its street), passes counted per 10 m bin along each piece (50 m median filter), pieces cut where the count steps, edges between nodes with one value each (length-weighted median); staggered crossings (< 40 m) and close loose ends contracted to one node, OSM roundabouts (junction=roundabout, 20 on the network) cut 10 m outside the ring and joined by straight spokes to one node at the trip-weighted mean of the cut points (busiest line through stays straight, e.g. Janowska at Okopowa), turning loops (< 250 m) dropped. Counts are flow-consistent (at most junctions the other edges add up exactly to the busiest). Drawing: no width cap; where a trunk splits / merges (other edges' trips <= 1.15 x trunk) the branches start side by side inside the trunk (slot by share of trips, leftmost turn on the left, round cap) and ease to their own centreline over 180 m; elsewhere a joint dot. Dworzec PKP: spur 469, last 116 m to the station square 633 (out-and-back shapes cut at their turning point, not 60 m before it) (everyone in and out, D/F/H also come from the east street, 164). TransitFrequency poster layout |
| `png/cartogram_<date>.png`, `png/cartogram_<date>_city.png` | Cartogram, all stops / city limits (dashed); `_city_nocap.png` is from the old version |
| `data/gis/city_limits.gpkg` | City limits, OSM relation 3016433 (via Nominatim, 49.4 km2) |
| `png/biala_podlaska_5min_departures_20261118.png` | 5-min walk frequency map, 300 dpi, made by TransitFrequency (see below) |
| `diagram.py` | Bus network poster built for a small town: corridors from the day's shapes (directions / carriageways / variants within 35 m folded; cut ends re-attached where the route continues), each route sample assigned to the nearest parallel corridor only; **OSM roundabouts straightened** to one node with straight spokes (as in the cartogram; not the one by Dworzec PKP, whose schematic needs it) and **parallel edges folded** (same two junctions, lengths within 1.5x + 50 m: two carriageways, the Białka frontage road beside the trunk, Sidorska / Al. Solidarności; repeated until stable, and short edges < 80 m closing a triangle contracted first: a split that rejoins at once, where a line on both sides counts as a terminal loop and the order search lets E / K swap for free, as on the Białka straight; edges with exactly the same lines fold up to 2.5x (one link, its two directions on different roads: C on the DK2 / Janowska ramps); **path fold**: an edge beside a path of up to 4 edges carrying all its lines, similar length, within 80 m (Hausdorff) is dropped (E's far-east variant: out along the north road and up the spur, back down it and along DK2; A's loop at Sławacinek is farther apart and stays)); lines repaired to stay continuous (corridor ends < 10 m apart linked first, so a gap closes in place: before, I was routed 286 m along Wolska to bridge a 3 m gap; gaps joined, close loose ends bridged, strays / stubs dropped), junction clusters joined by edges < 40 m merged into one (then parallel edges folded again); **per-edge line orders** (start from the best common order, then annealing per edge on junction crossings of lines passing through + a small cost per lane shifted; terminal loops excluded from the pair test), **junction corners relaxed** before drawing (a corner at a junction too tight for the bundles turning there drops the edge's next vertex, so the arc gets room for the innermost lane; decided per edge so lanes stay concentric: e.g. the 9-line bundle turning onto Al. 1000-lecia north of the station); lines drawn subway_like_map-style: each section simplified (30 m), corners replaced by circular arcs (radius up to 150 m, bends too tight for the bundle removed) so bundles bend as concentric arcs, lane changes eased over 120 m; **Dworzec PKP schematic**: a U hangs 690 m below the Sidorska / Al. 1000-lecia crossing; the legs leave the crossing straight down for 300 m side by side, as one bundle continuing Al. 1000-lecia (each pushed outwards by half the other's width), then move apart on an S-bend of two 160 m arcs (a smoothstep ease bent at ~37 m radius and crumpled the outer lanes) and meet in a half circle (before: both legs left diagonally from the crossing point, overlapped for ~80 m and the arms turned ~127° on ~50 m arcs, folding over each other); the arms (Sidorska west / east, Al. 1000-lecia north) approach the crossing shifted, over their last 300 m, to where their lines go on the legs (lane offsets kept through the turn: concentric arcs); the station is a vertical bar at its lowest point; A, C, I go round the U and split at the crossing north / west / east, D, F, H leave along Kościuszki straight east from the bottom, B (in and out along Sidorska west) ends at the bar; line order round the U: A / C / I pair their four sections at the crossing explicitly, by side (the more westerly arm with the west leg, the other with the east leg, so no line crosses itself: C, I Sidorska west + west leg, Al. 1000-lecia north + east leg; A Al. 1000-lecia north + west leg, Sidorska east + east leg; user 2026-11-17), crossings there cost x3, at the station (under the bar) x0.3; the crossing count also catches lines that pass a junction without sharing a section (interleaving ends); **no gaps at junctions**: a line piece that ends where the line goes on along other sections is drawn with a lead-in of the section its buses really continue on (`continuations()`: each route running along the section (near points 40 m and 150 m out, so routes merely sweeping past a merged junction do not count) is followed to the node and 120 m past it, the nearest other section of the line gets its trips; straightest only as fallback: E from Daleka north into Wolska, L north to Francuska Pętla, H from Al. Solidarności west along Sidorska, K Białka <-> cemetery at the Białka roundabout), so it eases out of that lane; audited 2026-11-17: at every junction where a line has 3+ sections the drawn join is the one most of its trips take; **extensions pale** (full width, colour washed 55 % towards the paper, drawn under the full-service lines; `--ext-style dashed|thin|pale|hollow`, chosen by the user from a comparison, 2026-11-17) point by point: where fewer than 25 % of the line's trips of the day pass (shapes within 45 m; uncovered schematic bits take their neighbours' value, runs under 150 m merged), so a line is solid right up to its main terminus (C to Robotnicza) and dashed beyond; black stop dots on the drawn lines, **arrowheads at stops a line passes one way only** (all its shapes of the day within 45 m run the same way: terminal loops such as E Wolska westbound -> Twarda 1 -> Twarda eastbound, I round Chłodna / Kołychawa; 12 on 18.11.2026), key for dot / arrow / dashed under the legend; **names for termini only** (main termini of a line + Dworzec PKP), white capsule + line bullets, placed as close to the stop as possible and aligned towards it; line ends: bar at a main terminus, white diamond at the end of an extended variant, whose last stop is named with diamond line bullets (round bullets: main termini); names kept inside the map frame |
| `png/bus_network_city.png`, `.pdf` | The poster (zoomed to the city limits, lines not trimmed); PDF = vector copy for print (sent to the transport authority) |
| `config/diagram_lines.json` | Diagram colours and legend names per line |
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
python cartogram.py [--city]             # flow cartogram (+ gpkg layer segments, one row per edge)
python diagram.py                    # network poster -> png/bus_network_city.png
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
  **No-bus streets** (`NO_BUS_STREETS`, by OSM name, excluded like living streets): Lotnicza and Torowa (user,
  2026-11-16): after the turn penalties D / F / H cut through them; they go round Rondo Błogosławionego Ks. and
  straight down Al. Jana Pawła II to Kościuszki. Residential streets that do carry buses (Orzechowa, Włoska,
  Kościuszki, Akademicka, Grzybowa, ...) are listed by a check in the session log: review that list when a route
  looks off.
  **Known bus ways** (`BUS_WAYS`, weight 1.0 whatever their class): OSM has almost no bus data for Biała Podlaska
  (no usable route relations, Overpass not relied on), so ways buses use against their class are listed by hand.
  So far: way 123907910, the unnamed frontage road north of the Białka trunk; E / K / T leave Białka at
  Powstańców and serve Białka 2 / Białka 3 there in both directions (user, Street View 2026-11-16).
  **Stops pinned to ways** (`STOP_WAYS`): Białka 2 / Białka 3 posts (170, 298, 353, 354) match only the frontage
  road (user, 2026-11-17: all buses use it); the Al. Solidarności 3 variants (E_a, E_b, E_k, K_E, K_F) had been
  matched to the trunk 11-22 m away. **Termini snapped hard** (`SNAP_END` 30 per metre vs 3): shapes start / end at
  the terminus post (were up to 52 m short: E inbound ended on Przyjaźni before Twarda 1; 75 shapes 5-360 m longer).
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
- Cartogram (`cartogram.py`): busiest street Narutowicza / Al. Jana Pawła II, 363 buses/day both directions
  (lines A B C D E F H I L). Counting per noded piece of the shape union (vertex-pair keys gave beads near stops,
  where each variant's shape is cut at its own stop post).

## Maps made in the sibling projects (edited in place with the user's OK, 2026-11-16)

Other projects are only changed for a reason and with the user's OK; their outputs are copied / written here.
- **TransitFrequency** (`D:\QGIS\TransitFrequency`, commits 187d9b5, ab3a8e4): city `bialapodlaska` in
  `script/cities.py` (reads `data/gtfs.zip` here + PolRegio/IC trains), `python script/run_city.py --city bialapodlaska
  --date 20261118` (6 min), poster `python script/render_poster.py --city bialapodlaska --date 18.11.2026 --dpi 300
  --out D:/QGIS/gtfs_biala_podlaska/png/biala_podlaska_5min_departures_20261118.png`. Two Biała-specific settings:
  `implicit_front_gate_m: 30` (detached-house blocks are mapped as one perimeter fence without gates, e.g. Orzechowa:
  1,817 houses get the walk area within 30 m, as through their front gate) and own colour classes: 25 linear classes
  1..205 (205 = the actual maximum, 61 buildings), `--min-lightness 50` so the lowest blues are less dark. Stop check: Brzegowa - Baza (depot) and Ekologiczna gated, Czosnówka 4 ~100 m off
  the walk network; all minor. Rerun after a new feed: `run_city.py ... --from gtfs`.
- **subway_like_map**: used for the first network posters, then **reverted completely** (2026-11-16, user: Biała is far
  smaller than the tram cities that renderer is tuned for; local revert commits 4b151d0..46c2cc0 there, its tree equals
  30ae893 again). The diagram is now `diagram.py` in this repo.
- Feed files and the GeoPackage are overwritten in place (QGIS keeps them open; deleting fails on Windows).
- The basemap in `data/basemap/` is a copy of `D:\QGIS\osm_basemapialapodlaska` (made by TransitFrequency).

- Diagram (2026-11-17): Chłodna / Twarda untangled (I no longer drawn on Wolska), Białka / Cmentarz Komunalny and
  Sidorska / Al. Solidarności as one corridor each (roundabouts straightened, parallel edges folded), A round the
  Dworzec PKP U the other way (no self-crossing), longer dashes on extensions; one-way stops as arrowheads (E's
  Wolska / Twarda terminal loop checked against all 8 E shapes: in along Wolska, out along Twarda).
  Crossing north of Dworzec PKP smoothed: U legs leave side by side straight down, arms shifted in (no folding).
  Loop ends join the section the buses really use (E Daleka -> Wolska curved from the north, L north to Francuska).
  E / K no longer swap on the Białka straight (triangles / repeated parallel fold); no kinks on Al. 1000-lecia north
  of the station (junction corners relaxed, U legs part on an S-bend of 160 m arcs).
  All loop-end joins audited against the routes (H / K at Białka and Al. Solidarności fixed); extensions drawn pale
  instead of dashed; I recoloured #009A8F -> #00838F (blue-teal, was too close to B's green on Sidorska west).
- Feed rebuilt 2026-11-17 (validator: 0 errors, 0 warnings): Białka stops pinned to the frontage road, termini
  reached by every shape. Node-by-node check of the diagram (Twarda 1, Daleka / Wolska, Daleka / Twarda, Białka
  roundabout, Cmentarz Komunalny, Sidorska / Al. Solidarności, Al. 1000-lecia / Sidorska, Francuska, DK2 /
  Janowska): every drawn join is the line's main relation there. GeoPackage and both cartograms regenerated from the
  rebuilt feed (busiest edge unchanged: 633 buses/day on the Dworzec PKP spur).
- Next: polish the subway-style map (legend texts, stations), maybe 10-min walk variant; feed expires
  31.12.2026 — a new export will be needed for 2027.
