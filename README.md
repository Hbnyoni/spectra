# SPECTRA

**Spatiotemporal Personal Exposure Characterization and TRAjectory Analyzer**

SPECTRA turns a raw wearable record, one position and one exposure reading per
minute, into a dwell-weighted, microenvironment-resolved exposure estimate for
each person, each day and each 100 m cell. It is the step between a monitor and
an epidemiological variable.

A personal monitor produces hundreds of thousands of readings per participant.
Averaging them naively fails in four specific ways, and each stage of SPECTRA
answers one:

| Failure | Stage that answers it |
|---|---|
| Readings are unevenly spaced, so a plain mean over-weights dense periods | **Weight**: dwell time per point, capped per observation |
| GPS drift makes a stationary person look mobile, and the reverse | **Quality**: haversine speed, implausible movement removed |
| A reading is carried across an arbitrarily long gap | **Weight**: per-record dwell cap |
| Speed alone cannot tell indoors from outdoors | **Classify**: speed plus building-footprint intersection |

## Installation

SPECTRA requires R 4.2 or later and the system libraries GDAL, GEOS and PROJ.

```r
install.packages(c("data.table", "sf", "terra", "shiny", "leaflet",
                   "ggplot2", "DBI", "RSQLite", "rmarkdown", "knitr", "digest"))
```

On Debian or Ubuntu the system libraries are:

```sh
sudo apt install libgdal-dev libgeos-dev libproj-dev libudunits2-dev
```

Clone the repository and source the engine:

```r
source("src/trajectory_analysis_core.R")
```

## Use

### Command line

The command-line driver runs the full pipeline on one input file and writes a
report and the segment tables:

```sh
Rscript src/run_analysis_cli.R \
    --input      examples/trajectory_example.csv \
    --buildings  buildings_cache.sqlite \
    --out        results/ \
    --dwell-cap  20 \
    --speed      100 \
    --buffer     20 \
    --grid       100
```

### Shiny application

```r
shiny::runApp("src")
```

The application provides upload, parameter control, interactive mapping and
report download for users who do not write R.

### Building footprint cache

The indoor test is a point-in-buffered-polygon query repeated once per record,
so footprints are held in an indexed SQLite cache built once and reused:

```sh
Rscript src/import_buildings_to_cache.R \
    --input  open_buildings_region.csv \
    --cache  buildings_cache.sqlite
```

Footprints may come from [OpenStreetMap](https://www.openstreetmap.org) or from
[Google Open Buildings](https://sites.research.google/open-buildings/). Which
one you use matters: see below.

## Input formats

CSV, GPX, RDS, GeoJSON and remote URLs, including Google Sheets share links. The
minimum a CSV needs is a participant identifier, a timestamp, longitude,
latitude and one measurement column.

## Parameters

| Parameter | Default | Purpose |
|---|---:|---|
| Teleport threshold | 100 km/h | Remove implausible coordinate jumps |
| Maximum dwell per record | 20 min | Limit contribution across long gaps |
| Outdoor speed threshold | 100 m/min | Speed-based classification |
| Building-footprint buffer | 20 m | Absorb GPS and footprint error |
| Grid-cell size | 100 m | Dwell-weighted spatial aggregation |
| Metric projection | EPSG:6933 | Distance, buffering and gridding |

## Outputs

- participant-day dwell-weighted means, with coverage and cells visited
- dwell-weighted means per 100 m cell
- segment summaries by hour, weekday, time of day, season, climate and mobility
- indoor and outdoor summaries
- a rendered R Markdown report

## A caution about the building basemap

The indoor/outdoor label is governed by building-map completeness, not only by
behaviour. In a published application to 343 participants in The Gambia, Kenya
and Mozambique, running the identical rule under OpenStreetMap and under Google
Open Buildings moved the outdoor share by up to 34.4 percentage points, and in
opposite directions in different countries.

The dwell-weighted exposure estimate was almost unaffected, differing by under
2.1% in every country. What changed was how much data survived a minimum
outdoor-time filter: a 120-minute floor retained 732 of 1,569 Gambian
participant-days under OpenStreetMap against 1,513 under Open Buildings.

Treat the indoor/outdoor fraction as a property of your basemap until you have
validated it. Run the pipeline under both maps; it is a one-argument change.

## Citing SPECTRA

If you use SPECTRA, please cite the software paper and the software itself.

## License

GNU General Public License v3.0. See `LICENSE.txt`.

Copyright (C) 2026 Place Alert Labs, Midlands State University.
