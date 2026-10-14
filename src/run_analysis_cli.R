#!/usr/bin/env Rscript
# ===========================================================================
# SPECTRA: Spatiotemporal Personal Exposure Characterization and TRAjectory
# Analyzer.
# Command line driver for the analysis engine.
#
# Copyright (c) 2026 Place Alert Labs, Midlands State University.
# All rights reserved. See LICENSE and NOTICE in the repository root.
#
# Author and lead developer: Handsome Bongani Nyoni.
# Disclosed to the Research and Innovation Hub, Midlands State University,
# on 17 September 2026 (SPECTRA IP Disclosure).
#
# Proprietary and confidential. No part of this file may be used, copied,
# modified or distributed without the written permission of the copyright
# holder. Academic use by request.
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ncdf4)
  library(RSQLite)
  library(lubridate)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

parse_args <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- list()
  for (arg in raw) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

as_num <- function(x, default) {
  if (is.null(x)) return(default)
  val <- suppressWarnings(as.numeric(x))
  if (is.na(val)) default else val
}

usage <- function() {
  cat(
"Standalone SPECTRA analysis runner

Required:
  --input=/path/to/trajectory.csv
  Input may be .csv, .rds, or .RDS.

Common options:
  --outdir=outputs_cli
  --lat=lat --lon=lon --pid=PID --date=date
  --pm25=pm25 --no2=no2
  --temp=temp --humidity=rh --country=country
  --csv-epsg=4326
  --cell-size=100
  --max-dwell=20
  --speed-limit=100
  --speed-outdoor-threshold=100
  --building-buffer=5
  --use-buildings=true
  --fetch-buildings=true
  --building-cache=buildings_cache.sqlite
  --preloaded-buildings=buildings_preloaded.csv
  --boundary=/path/to/boundary.geojson
  --exclude-indoor=false
  --building-raster=<path to a building presence GeoTIFF; overrides OSM footprints>

Example:
  Rscript run_analysis_cli.R --input=data.csv --outdir=quick_run --lat=Latitude --lon=Longitude --pid=participant --date=datetime --pm25=PM25 --no2=NO2 --cell-size=100 --fetch-buildings=true

RDS example using your segmentation files:
  Rscript run_analysis_cli.R --input=data/trajectory.rds --outdir=results/ --cell-size=100 --fetch-buildings=true
\n")
}

load_input_trajectory <- function(file_path, args) {
  ext <- tolower(tools::file_ext(file_path))

  if (ext == "csv") {
    return(load_csv_trajectory(
      file_path,
      lat_field = args$lat %||% "lat",
      lon_field = args$lon %||% "lon",
      pid_field = args$pid %||% "PID",
      date_field = args$date %||% "date",
      pm25_field = args$pm25 %||% "pm25",
      no2_field = args$no2 %||% "no2",
      csv_epsg = as.integer(as_num(args[["csv-epsg"]], 4326))
    ))
  }

  if (ext %in% c("rds", "rda")) {
    dt <- readRDS(file_path)
    dt <- as.data.table(dt)

    lat_field <- args$lat %||% "lat"
    lon_field <- args$lon %||% "lon"
    pid_field <- args$pid %||% "PID"
    date_field <- args$date %||% "date"
    pm25_field <- args$pm25 %||% "pm25"
    no2_field <- args$no2 %||% "no2"
    temp_field <- args$temp %||% "temp"
    humidity_field <- args$humidity %||% "rh"
    country_field <- args$country %||% "country"

    missing_required <- setdiff(c(lat_field, lon_field, pid_field, date_field), names(dt))
    if (length(missing_required) > 0) {
      stop("Missing required RDS columns: ", paste(missing_required, collapse = ", "),
           "\nAvailable columns: ", paste(names(dt), collapse = ", "))
    }

    if (lat_field != "lat") setnames(dt, lat_field, "lat")
    if (lon_field != "lon") setnames(dt, lon_field, "lon")
    if (pid_field != "PID") setnames(dt, pid_field, "PID")
    if (date_field != "date") setnames(dt, date_field, "date")
    if (!is.null(pm25_field) && pm25_field != "None" && pm25_field %in% names(dt) && pm25_field != "pm25") {
      setnames(dt, pm25_field, "pm25")
    }
    if (!is.null(no2_field) && no2_field != "None" && no2_field %in% names(dt) && no2_field != "no2") {
      setnames(dt, no2_field, "no2")
    }
    if (!is.null(temp_field) && temp_field != "None" && temp_field %in% names(dt) && !"temperature_c" %in% names(dt)) {
      setnames(dt, temp_field, "temperature_c")
    }
    if (!is.null(humidity_field) && humidity_field != "None" && humidity_field %in% names(dt) && !"humidity_pct" %in% names(dt)) {
      setnames(dt, humidity_field, "humidity_pct")
    }
    if (!is.null(country_field) && country_field != "None" && country_field %in% names(dt) && country_field != "country") {
      setnames(dt, country_field, "country")
    }

    if (!"pm25" %in% names(dt)) dt[, pm25 := NA_real_]
    if (!"no2" %in% names(dt)) dt[, no2 := NA_real_]

    dt[, lat := as.numeric(lat)]
    dt[, lon := as.numeric(lon)]
    dt[, PID := as.character(PID)]
    dt[, date := as.POSIXct(date, tz = "UTC")]
    dt[, pm25 := as.numeric(pm25)]
    dt[, no2 := as.numeric(no2)]
    if ("temperature_c" %in% names(dt)) dt[, temperature_c := as.numeric(temperature_c)]
    if ("humidity_pct" %in% names(dt)) dt[, humidity_pct := as.numeric(humidity_pct)]

    dt <- dt[!is.na(lat) & !is.na(lon) & !is.na(PID) & !is.na(date)]
    return(dt)
  }

  stop("Unsupported input extension: ", ext, ". Use .csv or .rds.")
}

write_exposure_netcdf <- function(result, file, exp1_label = "PM2.5", exp2_label = "NO2") {
  dg <- as.data.table(result$daily_grid)
  fg <- as.data.table(result$full_grid)
  gp <- result$grid_params

  if (nrow(dg) == 0 || nrow(fg) == 0) {
    stop("Cannot write NetCDF because daily_grid or full_grid is empty.")
  }

  mv <- -999.0
  mvi <- -1L

  lon_vals <- sort(unique(round(fg$lon, 6)))
  lat_vals <- sort(unique(round(fg$lat, 6)))
  nlon <- length(lon_vals)
  nlat <- length(lat_vals)

  dg[, date := as.Date(date)]
  dates <- sort(unique(dg$date))
  ref_date <- as.Date("1970-01-01")
  time_vals <- as.integer(dates - ref_date)
  ntime <- length(time_vals)

  use_compression <- tryCatch({
    test_file <- tempfile(fileext = ".nc")
    nc <- ncdf4::nc_create(
      test_file,
      list(ncdf4::ncvar_def("t", "", ncdf4::ncdim_def("x", "", 1), compression = 6)),
      force_v4 = TRUE
    )
    ncdf4::nc_close(nc)
    unlink(test_file)
    TRUE
  }, error = function(e) FALSE)
  cmp <- if (use_compression) 6L else NA_integer_

  lon_dim <- ncdf4::ncdim_def("lon", "degrees_east", lon_vals,
                              longname = "longitude of grid cell centre")
  lat_dim <- ncdf4::ncdim_def("lat", "degrees_north", lat_vals,
                              longname = "latitude of grid cell centre")
  time_dim <- ncdf4::ncdim_def("time", paste0("days since ", format(ref_date)),
                               time_vals, unlim = TRUE,
                               longname = "date of observation")

  var_grid_id <- ncdf4::ncvar_def("grid_id", "1", list(lon_dim, lat_dim),
                                  missval = mvi, prec = "integer", compression = cmp,
                                  longname = "Unique grid cell identifier")
  var_labeled <- ncdf4::ncvar_def("is_labeled", "1", list(lon_dim, lat_dim),
                                  missval = mvi, prec = "short", compression = cmp,
                                  longname = "1=cell has GPS observations, 0=never observed")
  var_crs <- ncdf4::ncvar_def("crs", "1", list(), prec = "integer",
                              longname = "WGS 84 coordinate reference system")
  var_dwell <- ncdf4::ncvar_def("dwell_time", "minutes",
                                list(lon_dim, lat_dim, time_dim), missval = mv,
                                compression = cmp,
                                longname = "Total dwell time in cell on this date")
  var_pm25 <- ncdf4::ncvar_def("pm25", "",
                               list(lon_dim, lat_dim, time_dim), missval = mv,
                               compression = cmp,
                               longname = paste0("Time-weighted mean ", exp1_label))
  var_no2 <- ncdf4::ncvar_def("no2", "",
                              list(lon_dim, lat_dim, time_dim), missval = mv,
                              compression = cmp,
                              longname = paste0("Time-weighted mean ", exp2_label))
  var_indoor <- ncdf4::ncvar_def("indoor_pct", "percent",
                                 list(lon_dim, lat_dim, time_dim), missval = mv,
                                 compression = cmp,
                                 longname = "Percentage of dwell classified as indoor")
  var_tstart <- ncdf4::ncvar_def("time_start", "seconds since 1970-01-01T00:00:00Z",
                                 list(lon_dim, lat_dim, time_dim), missval = mv,
                                 compression = cmp,
                                 longname = "Unix timestamp of first GPS observation")
  var_tend <- ncdf4::ncvar_def("time_end", "seconds since 1970-01-01T00:00:00Z",
                               list(lon_dim, lat_dim, time_dim), missval = mv,
                               compression = cmp,
                               longname = "Unix timestamp of last GPS observation")
  var_date <- ncdf4::ncvar_def("date", "", list(time_dim), prec = "char",
                               longname = "Calendar date of observation (YYYY-MM-DD)")

  nc <- ncdf4::nc_create(
    file,
    list(var_grid_id, var_labeled, var_crs, var_dwell, var_pm25, var_no2,
         var_indoor, var_tstart, var_tend, var_date),
    force_v4 = use_compression
  )
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)

  ncdf4::ncatt_put(nc, "time", "calendar", "standard")
  ncdf4::ncatt_put(nc, "time", "axis", "T")
  ncdf4::ncvar_put(nc, var_date, format(dates, "%Y-%m-%d"))
  ncdf4::ncatt_put(nc, "lon", "axis", "X")
  ncdf4::ncatt_put(nc, "lon", "standard_name", "longitude")
  ncdf4::ncatt_put(nc, "lat", "axis", "Y")
  ncdf4::ncatt_put(nc, "lat", "standard_name", "latitude")

  ncdf4::ncvar_put(nc, var_crs, 0L)
  ncdf4::ncatt_put(nc, "crs", "grid_mapping_name", "latitude_longitude")
  ncdf4::ncatt_put(nc, "crs", "epsg_code", "EPSG:4326")
  ncdf4::ncatt_put(nc, "crs", "spatial_ref", sf::st_crs(4326)$wkt)
  ncdf4::ncatt_put(nc, "crs", "crs_wkt", sf::st_crs(4326)$wkt)
  ncdf4::ncatt_put(nc, "crs", "semi_major_axis", 6378137.0)
  ncdf4::ncatt_put(nc, "crs", "inverse_flattening", 298.257223563)
  ncdf4::ncatt_put(nc, "crs", "longitude_of_prime_meridian", 0.0)

  for (vname in c("grid_id", "is_labeled", "dwell_time", "pm25", "no2",
                  "indoor_pct", "time_start", "time_end")) {
    ncdf4::ncatt_put(nc, vname, "grid_mapping", "crs")
    ncdf4::ncatt_put(nc, vname, "coordinates", "lon lat")
  }

  fg_lon_idx <- match(round(fg$lon, 6), lon_vals)
  fg_lat_idx <- match(round(fg$lat, 6), lat_vals)
  fg_valid <- which(!is.na(fg_lon_idx) & !is.na(fg_lat_idx))
  fg_idx_by_cell <- data.table(
    cell_id = fg$cell_id,
    lon_idx = fg_lon_idx,
    lat_idx = fg_lat_idx
  )

  gid_mat <- matrix(mvi, nrow = nlon, ncol = nlat)
  gid_mat[fg_lon_idx[fg_valid] + (fg_lat_idx[fg_valid] - 1L) * nlon] <- fg$cell_id[fg_valid]
  ncdf4::ncvar_put(nc, var_grid_id, gid_mat)

  labeled_mat <- matrix(0L, nrow = nlon, ncol = nlat)
  obs_cells <- fg_valid[fg$is_labeled[fg_valid] == 1L]
  labeled_mat[fg_lon_idx[obs_cells] + (fg_lat_idx[obs_cells] - 1L) * nlon] <- 1L
  ncdf4::ncvar_put(nc, var_labeled, labeled_mat)

  dg_idx <- merge(
    dg[, .(row_id = .I, cell_id)],
    fg_idx_by_cell,
    by = "cell_id",
    all.x = TRUE,
    sort = FALSE
  )[order(row_id)]
  dg_lon_idx <- dg_idx$lon_idx
  dg_lat_idx <- dg_idx$lat_idx

  for (ti in seq_len(ntime)) {
    rows <- which(dg$date == dates[ti] & !is.na(dg_lon_idx) & !is.na(dg_lat_idx))
    mat_idx <- dg_lon_idx[rows] + (dg_lat_idx[rows] - 1L) * nlon

    dwell_sl <- matrix(mv, nrow = nlon, ncol = nlat)
    pm25_sl <- matrix(mv, nrow = nlon, ncol = nlat)
    no2_sl <- matrix(mv, nrow = nlon, ncol = nlat)
    ip_sl <- matrix(mv, nrow = nlon, ncol = nlat)
    tstart_sl <- matrix(mv, nrow = nlon, ncol = nlat)
    tend_sl <- matrix(mv, nrow = nlon, ncol = nlat)

    if (length(rows)) {
      dwell_sl[mat_idx] <- dg$dwell[rows]
      ok_pm <- rows[!is.na(dg$pm25[rows])]
      if (length(ok_pm)) pm25_sl[dg_lon_idx[ok_pm] + (dg_lat_idx[ok_pm] - 1L) * nlon] <- dg$pm25[ok_pm]
      ok_no2 <- rows[!is.na(dg$no2[rows])]
      if (length(ok_no2)) no2_sl[dg_lon_idx[ok_no2] + (dg_lat_idx[ok_no2] - 1L) * nlon] <- dg$no2[ok_no2]
      ok_ip <- rows[!is.na(dg$indoor_pct[rows])]
      if (length(ok_ip)) ip_sl[dg_lon_idx[ok_ip] + (dg_lat_idx[ok_ip] - 1L) * nlon] <- dg$indoor_pct[ok_ip]

      if ("datetime_start" %in% names(dg)) {
        ts_vals <- suppressWarnings(as.numeric(as.POSIXct(
          dg$datetime_start[rows], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
        te_vals <- suppressWarnings(as.numeric(as.POSIXct(
          dg$datetime_end[rows], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
        ok_ts <- rows[!is.na(ts_vals)]
        if (length(ok_ts)) tstart_sl[dg_lon_idx[ok_ts] + (dg_lat_idx[ok_ts] - 1L) * nlon] <- ts_vals[!is.na(ts_vals)]
        ok_te <- rows[!is.na(te_vals)]
        if (length(ok_te)) tend_sl[dg_lon_idx[ok_te] + (dg_lat_idx[ok_te] - 1L) * nlon] <- te_vals[!is.na(te_vals)]
      }
    }

    ncdf4::ncvar_put(nc, var_dwell, dwell_sl, start = c(1, 1, ti), count = c(nlon, nlat, 1))
    ncdf4::ncvar_put(nc, var_pm25, pm25_sl, start = c(1, 1, ti), count = c(nlon, nlat, 1))
    ncdf4::ncvar_put(nc, var_no2, no2_sl, start = c(1, 1, ti), count = c(nlon, nlat, 1))
    ncdf4::ncvar_put(nc, var_indoor, ip_sl, start = c(1, 1, ti), count = c(nlon, nlat, 1))
    ncdf4::ncvar_put(nc, var_tstart, tstart_sl, start = c(1, 1, ti), count = c(nlon, nlat, 1))
    ncdf4::ncvar_put(nc, var_tend, tend_sl, start = c(1, 1, ti), count = c(nlon, nlat, 1))
  }

  ncdf4::ncatt_put(nc, 0, "Conventions", "CF-1.8")
  ncdf4::ncatt_put(nc, 0, "title", "SPECTRA Exposure Grid")
  ncdf4::ncatt_put(nc, 0, "institution", "CHEAQI-MNCH")
  ncdf4::ncatt_put(nc, 0, "source", "SPECTRA command-line analysis")
  ncdf4::ncatt_put(nc, 0, "history", paste("Created", as.character(Sys.time())))
  ncdf4::ncatt_put(nc, 0, "grid_cell_size_m", gp$cell_size)
  ncdf4::ncatt_put(nc, 0, "date_start", as.character(min(dates)))
  ncdf4::ncatt_put(nc, 0, "date_end", as.character(max(dates)))
  ncdf4::ncatt_put(nc, 0, "n_dates", ntime)
  ncdf4::ncatt_put(nc, 0, "n_cells_total", nrow(fg))
  ncdf4::ncatt_put(nc, 0, "n_cells_labeled", sum(fg$is_labeled == 1L, na.rm = TRUE))
  ncdf4::ncatt_put(nc, 0, "n_cells_unlabeled", sum(fg$is_labeled == 0L, na.rm = TRUE))
  ncdf4::ncatt_put(nc, 0, "missing_value", mv)

  ncdf4::nc_close(nc)
  on.exit(NULL, add = FALSE)
  invisible(file)
}

load_buildings_for_data <- function(dt, cache_db, preloaded_csv, fetch_missing = TRUE, buffer_deg = 0.02) {
  lat_vals <- dt$lat[is.finite(dt$lat)]
  lon_vals <- dt$lon[is.finite(dt$lon)]
  if (!length(lat_vals) || !length(lon_vals)) return(NULL)

  minlat <- min(lat_vals) - buffer_deg
  maxlat <- max(lat_vals) + buffer_deg
  minlon <- min(lon_vals) - buffer_deg
  maxlon <- max(lon_vals) + buffer_deg
  message(sprintf("Building bbox: lat %.5f..%.5f, lon %.5f..%.5f", minlat, maxlat, minlon, maxlon))

  if (!is.null(cache_db) && nzchar(cache_db)) {
    init_building_cache(cache_db)
    bldgs <- tryCatch(
      get_buildings_from_cache(cache_db, minlat, maxlat, minlon, maxlon, fetch_missing = fetch_missing),
      error = function(e) {
        message("SQLite/OSM building load failed: ", e$message)
        NULL
      }
    )
    if (!is.null(bldgs) && nrow(bldgs) > 0) {
      message("Loaded ", nrow(bldgs), " buildings from SQLite/OSM cache")
      return(bldgs)
    }
  }

  if (!is.null(preloaded_csv) && file.exists(preloaded_csv)) {
    bdt <- tryCatch(fread(preloaded_csv, showProgress = FALSE), error = function(e) NULL)
    if (!is.null(bdt) && "geometry_wkt" %in% names(bdt) && nrow(bdt) > 0) {
      bldgs <- tryCatch({
        x <- st_sf(geometry = st_as_sfc(bdt$geometry_wkt, crs = 4326L))
        x <- suppressWarnings(st_make_valid(x))
        x <- x[!st_is_empty(x), , drop = FALSE]
        bbox_poly <- st_as_sfc(st_bbox(c(xmin = minlon, ymin = minlat, xmax = maxlon, ymax = maxlat), crs = 4326L))
        hits <- st_intersects(x, bbox_poly, sparse = FALSE)[, 1]
        x[hits, , drop = FALSE]
      }, error = function(e) {
        message("Preloaded building CSV failed: ", e$message)
        NULL
      })
      if (!is.null(bldgs) && nrow(bldgs) > 0) {
        message("Loaded ", nrow(bldgs), " buildings from preloaded CSV")
        return(bldgs)
      }
    }
  }

  message("No buildings found for this dataset extent")
  NULL
}

main <- function() {
  args <- parse_args()
  if (isTRUE(as_bool(args$help, FALSE))) {
    usage()
    quit(status = 0)
  }
  if (is.null(args$input)) {
    usage()
    quit(status = 1)
  }

  source("trajectory_analysis_core.R")

  input <- args$input
  outdir <- args$outdir %||% "outputs_cli"
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message("Reading trajectory input: ", input)
  dt <- load_input_trajectory(input, args)
  message("Loaded ", format(nrow(dt), big.mark = ","), " valid points")

  boundary_sf <- NULL
  if (!is.null(args$boundary) && file.exists(args$boundary)) {
    message("Loading boundary: ", args$boundary)
    boundary_sf <- load_boundary_layer(args$boundary)
  }

  buildings_sf <- NULL
  if (as_bool(args[["use-buildings"]], TRUE)) {
    buildings_sf <- load_buildings_for_data(
      dt,
      cache_db = args[["building-cache"]] %||% "buildings_cache.sqlite",
      preloaded_csv = args[["preloaded-buildings"]] %||% "buildings_preloaded.csv",
      fetch_missing = as_bool(args[["fetch-buildings"]], TRUE)
    )
  }

  message("Running analysis...")
  result <- run_analysis(
    dt,
    boundary_sf = boundary_sf,
    buildings_sf = buildings_sf,
    fetch_osm = FALSE,
    cache_db_path = args[["building-cache"]] %||% "buildings_cache.sqlite",
    cell_size_m = as_num(args[["cell-size"]], 100),
    max_dwell_min = as_num(args[["max-dwell"]], 20),
    speed_threshold_kmh = as_num(args[["speed-limit"]], 100),
    remove_teleports = as_bool(args[["remove-teleports"]], TRUE),
    building_buffer_m = as_num(args[["building-buffer"]], 5),
    exclude_indoor = as_bool(args[["exclude-indoor"]], FALSE),
    speed_outdoor_threshold = as_num(args[["speed-outdoor-threshold"]], 100),
    use_building_validation = as_bool(args[["use-buildings"]], TRUE),
    building_raster = args[["building-raster"]]
  )

  message("Writing outputs to: ", normalizePath(outdir, mustWork = FALSE))
  fwrite(result$points_classified, file.path(outdir, "points_classified.csv"))
  fwrite(result$person_day_cell, file.path(outdir, "person_day_cell.csv"))
  fwrite(result$person_day, file.path(outdir, "person_day.csv"))
  fwrite(result$participant_exposure, file.path(outdir, "participant_exposure.csv"))
  fwrite(result$daily_grid, file.path(outdir, "daily_grid.csv"))
  fwrite(result$full_grid, file.path(outdir, "full_grid.csv"))
  fwrite(result$validation_stats, file.path(outdir, "validation_stats.csv"))
  fwrite(as.data.table(result$stats), file.path(outdir, "summary_stats.csv"))
  saveRDS(result, file.path(outdir, "analysis_result.rds"), compress = "gzip")

  nc_path <- file.path(outdir, "exposure_grid.nc")
  message("Writing NetCDF: ", nc_path)
  write_exposure_netcdf(
    result,
    nc_path,
    exp1_label = args[["exp1-label"]] %||% "PM2.5",
    exp2_label = args[["exp2-label"]] %||% "NO2"
  )

  message("Done.")
  message("Buildings used: ", if (is.null(result$buildings_sf)) 0L else nrow(result$buildings_sf))
  message("NetCDF: ", nc_path)
}

if (sys.nframe() == 0) {
  main()
}
