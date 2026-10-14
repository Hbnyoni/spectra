# ===========================================================================
# SPECTRA: Spatiotemporal Personal Exposure Characterization and TRAjectory
# Analyzer.
# Core analysis engine: trajectory loading, dwell and transit segmentation,
# indoor/outdoor classification, exposure joining and dose integration.
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
# ============================================================================
# SPECTRA — Spatiotemporal Personal Exposure Characterization and TRAjectory Analyzer
# Core Analysis Engine  |  CHEAQI-MNCH
# ============================================================================

library(data.table)

EPSG_WGS84  <- 4326
EPSG_METRIC <- 6933

# Safe CRS enforcement: sets `epsg` if CRS is absent OR has no valid WKT.
# NOTE: never use return() inside the tryCatch expr — it exits the whole function.
.ensure_crs <- function(x, epsg = 4326L) {
  needs_crs <- tryCatch({
    crs <- sf::st_crs(x)
    if (is.na(crs)) {
      TRUE
    } else {
      wkt <- tryCatch(crs$wkt, error = function(e) NA_character_)
      is.null(wkt) || is.na(wkt) || !nzchar(trimws(wkt))
    }
  }, error = function(e) TRUE)
  if (isTRUE(needs_crs)) sf::st_crs(x) <- epsg
  x
}

# Central safe wrapper for ALL sf::st_transform calls.
# 1. Ensures source has a valid CRS (defaults to assume_from if missing).
# 2. Attempts the transform.
# 3. If it still fails, force-sets the target CRS without reprojecting
#    (geometry stays in place but at least the app doesn't crash).
.safe_transform <- function(x, to_crs, assume_from = 4326L) {
  x <- .ensure_crs(x, assume_from)
  tryCatch(
    sf::st_transform(x, to_crs),
    error = function(e) {
      warning("[SPECTRA] st_transform failed (", conditionMessage(e),
              ") — forcing CRS to ", to_crs, " without reprojection.")
      suppressWarnings(sf::st_crs(x) <- to_crs)
      x
    }
  )
}

# Reproject to WGS84 using the safe wrapper.
.to_wgs84 <- function(x) {
  x <- .ensure_crs(x, 4326L)
  tryCatch({
    is_wgs84 <- tryCatch(isTRUE(sf::st_crs(x) == sf::st_crs(4326L)), error = function(e) FALSE)
    if (!is_wgs84) x <- .safe_transform(x, 4326L, assume_from = 4326L)
    x
  }, error = function(e) { suppressWarnings(sf::st_crs(x) <- 4326L); x })
}

# ============================================================================
# Data Loading & Standardisation
# ============================================================================

#' Load and standardise a CSV trajectory file
#'
#' @param file_path Path to CSV file
#' @param lat_field Column name for latitude
#' @param lon_field Column name for longitude
#' @param pid_field Column name for participant ID
#' @param date_field Column name for date/time
#' @param pm25_field Column name for PM2.5
#' @param no2_field Column name for NO2
#' @param csv_epsg EPSG code of the input CRS (default 4326)
#' @return data.table with standardised columns
#' Load a GPX file from a wearable / GPS device
#'
#' Reads track_points (or waypoints) and returns a standardised data.table.
#' Wearable extensions (heart_rate, cadence, speed, temperature) are preserved
#' as additional columns. pm25 / no2 are set to NA — map them from your own
#' merged sensor data using the field-mapping UI.
#'
#' @param file_path  Path to the .gpx file
#' @param pid_label  Default participant ID when the GPX has no track_fid/name
#' @return data.table with lat, lon, PID, date, pm25, no2 (+ any wearable cols)
load_gpx_trajectory <- function(file_path, pid_label = "participant_1") {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("sf package required to read GPX files.")

  pts <- tryCatch(
    sf::st_read(file_path, layer = "track_points", quiet = TRUE),
    error = function(e) NULL)
  if (is.null(pts) || nrow(pts) == 0)
    pts <- tryCatch(
      sf::st_read(file_path, layer = "waypoints", quiet = TRUE),
      error = function(e) stop("No track_points or waypoints found in GPX file."))
  if (nrow(pts) == 0) stop("GPX file contains no points.")

  coords <- sf::st_coordinates(pts)
  dt <- data.table(lon = coords[, 1], lat = coords[, 2])

  # Time
  if ("time" %in% names(pts)) {
    dt[, date := as.POSIXct(pts$time, tz = "UTC")]
  } else {
    dt[, date := seq(Sys.time(), by = "1 min", length.out = nrow(dt))]
    warning("No time field in GPX — synthetic timestamps assigned.")
  }

  # Participant ID
  if ("track_fid" %in% names(pts)) {
    dt[, PID := paste0("track_", as.character(pts$track_fid))]
  } else if ("name" %in% names(pts) && !all(is.na(pts$name))) {
    dt[, PID := as.character(pts$name)]
  } else {
    dt[, PID := pid_label]
  }

  # Elevation
  if ("ele" %in% names(pts)) dt[, elevation_m := as.numeric(pts$ele)]

  # Wearable extensions
  for (col in c("heart_rate", "cadence", "speed", "power", "temperature")) {
    if (col %in% names(pts)) dt[, (col) := as.numeric(pts[[col]])]
  }

  # Placeholder exposure columns — user maps their sensor data via field mapping
  dt[, pm25 := NA_real_]
  dt[, no2  := NA_real_]

  dt[!is.na(lat) & !is.na(lon) & !is.na(date) & !is.na(PID)]
}

#' Load trajectory data from an RDS file (data.frame or data.table)
#'
#' @param file_path Path to .rds file containing a data.frame or data.table
#' @return data.table
load_rds_trajectory <- function(file_path) {
  obj <- readRDS(file_path)
  # Reject if it looks like a SPECTRA analysis result rather than raw data
  if (is.list(obj) && "points_classified" %in% names(obj))
    stop("This RDS contains a previous analysis result, not raw trajectory data. ",
         "Use 'Restore Previous Session' below to reload it.")
  if (!is.data.frame(obj) && !is.data.table(obj))
    stop("RDS file must contain a data.frame or data.table of trajectory points.")
  as.data.table(obj)
}

#' Load trajectory data from a GeoJSON file (Point or LineString features)
#'
#' Supports FeatureCollections of Point, MultiPoint, or LineString geometries.
#' Coordinates are extracted as lon/lat columns; all other properties are preserved.
#'
#' @param file_path Path to .geojson / .json file
#' @return data.table with lon, lat columns plus all feature properties
load_geojson_trajectory <- function(file_path) {
  if (!requireNamespace("sf", quietly = TRUE))
    stop("sf package required to read GeoJSON files.")

  gdf <- sf::st_read(file_path, quiet = TRUE)
  if (nrow(gdf) == 0) stop("GeoJSON file contains no features.")

  gdf <- .ensure_crs(gdf, 4326L)
  gdf <- .to_wgs84(gdf)

  geom_type <- as.character(sf::st_geometry_type(gdf[1, ]))

  if (geom_type %in% c("POINT", "MULTIPOINT")) {
    # Explode multipoints to individual rows, then extract
    gdf <- sf::st_cast(gdf, "POINT", warn = FALSE)
    coords <- sf::st_coordinates(gdf)
    dt <- as.data.table(sf::st_drop_geometry(gdf))
    dt[, lon := coords[, 1]]
    dt[, lat := coords[, 2]]
  } else if (geom_type %in% c("LINESTRING", "MULTILINESTRING")) {
    # Explode lines to vertices
    gdf_pts <- sf::st_cast(gdf, "POINT", warn = FALSE)
    coords   <- sf::st_coordinates(gdf_pts)
    dt <- as.data.table(sf::st_drop_geometry(gdf_pts))
    dt[, lon := coords[, 1]]
    dt[, lat := coords[, 2]]
  } else {
    stop("GeoJSON must contain Point, MultiPoint, LineString, or MultiLineString features.")
  }

  dt[!is.na(lon) & !is.na(lat)]
}

#' Load trajectory data from a remote URL (CSV or Google Sheets)
#'
#' @param url  HTTP/HTTPS URL to a CSV file or Google Sheets share link
#' @return data.table
load_url_trajectory <- function(url) {
  # Convert Google Sheets view/edit links to CSV export
  if (grepl("docs.google.com/spreadsheets", url, fixed = TRUE)) {
    url <- sub("/edit.*$",    "/export?format=csv", url)
    url <- sub("/pub.*$",     "/export?format=csv", url)
    url <- sub("/view.*$",    "/export?format=csv", url)
    url <- sub("/htmlview.*$","/export?format=csv", url)
  }
  tryCatch(
    fread(url, showProgress = FALSE),
    error = function(e) stop(paste("Failed to load URL:", e$message))
  )
}

load_csv_trajectory <- function(file_path,
                                lat_field = "lat",
                                lon_field = "lon",
                                pid_field = "PID",
                                date_field = "date",
                                pm25_field = "pm25",
                                no2_field = "no2",
                                csv_epsg = 4326) {
  
  dt <- fread(file_path, showProgress = FALSE)
  
  # Validate required fields exist before renaming
  coord_fields <- c(lat_field, lon_field, pid_field, date_field)
  missing_required <- setdiff(coord_fields, names(dt))
  if (length(missing_required) > 0) {
    stop(paste("Missing required fields:", paste(missing_required, collapse = ", "),
               "\nAvailable columns:", paste(names(dt), collapse = ", ")))
  }

  # Rename to standard names (only if source != target to avoid data.table warning)
  if (lat_field  != "lat"  && lat_field  %in% names(dt)) setnames(dt, lat_field,  "lat")
  if (lon_field  != "lon"  && lon_field  %in% names(dt)) setnames(dt, lon_field,  "lon")
  if (pid_field  != "PID"  && pid_field  %in% names(dt)) setnames(dt, pid_field,  "PID")
  if (date_field != "date" && date_field %in% names(dt)) setnames(dt, date_field, "date")
  # Exposure fields are optional
  if (!is.null(pm25_field) && pm25_field != "None" && pm25_field != "pm25" &&
      pm25_field %in% names(dt)) setnames(dt, pm25_field, "pm25")
  if (!is.null(no2_field)  && no2_field  != "None" && no2_field  != "no2"  &&
      no2_field  %in% names(dt)) setnames(dt, no2_field,  "no2")
  if (!"pm25" %in% names(dt)) dt[, pm25 := NA_real_]
  if (!"no2"  %in% names(dt)) dt[, no2  := NA_real_]

  # Coerce types — use dt[["col"]] to avoid data.table scope resolving to base R functions
  dt[, lat := as.numeric(dt[["lat"]])]
  dt[, lon := as.numeric(dt[["lon"]])]
  dt[, PID := as.character(dt[["PID"]])]
  dt[, date := {
    d <- as.character(dt[["date"]])
    parsed <- suppressWarnings(
      lubridate::parse_date_time(d, orders = c(
        "ymd HMS", "ymd HM", "ymd",
        "dmy HMS", "dmy HM", "dmy",
        "mdy HMS", "mdy HM", "mdy",
        "ymd_HMS", "dmy_HMS"
      ), tz = "UTC", quiet = TRUE)
    )
    # Fallback for any remaining unparsed
    still_na <- is.na(parsed)
    if (any(still_na)) {
      parsed[still_na] <- suppressWarnings(as.POSIXct(d[still_na], tz = "UTC"))
    }
    parsed
  }]
  dt[, pm25 := as.numeric(pm25)]
  dt[, no2 := as.numeric(no2)]
  
  # Remove rows with NA in required navigation columns only (exposure is optional)
  dt <- dt[!is.na(lat) & !is.na(lon) & !is.na(PID) & !is.na(date)]

  # Reproject if needed
  csv_epsg <- if (is.null(csv_epsg) || is.na(csv_epsg)) 4326L else as.integer(csv_epsg)
  if (csv_epsg != 4326L) {
    if (requireNamespace("sf", quietly = TRUE)) {
      pts <- sf::st_as_sf(dt, coords = c("lon", "lat"), crs = csv_epsg)
      pts <- .safe_transform(pts, 4326L, csv_epsg)
      coords <- sf::st_coordinates(pts)
      dt[, lon := coords[, 1]]
      dt[, lat := coords[, 2]]
    } else {
      warning("sf package required for reprojection. Assuming WGS84.")
    }
  }
  
  return(dt)
}

# ============================================================================
# Boundary & Building Loading
# ============================================================================

#' Load a boundary polygon layer
#' @param file_path Path to GeoJSON/GPKG/SHP file
#' @return sf object or NULL
load_boundary_layer <- function(file_path) {
  if (is.null(file_path) || file_path == "") return(NULL)
  
  tryCatch({
    gdf <- sf::st_read(file_path, quiet = TRUE)
    if (nrow(gdf) == 0) return(NULL)
    
    gdf <- .safe_transform(gdf, 4326L)

    # Keep only polygons
    geom_types <- sf::st_geometry_type(gdf)
    gdf <- gdf[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
    # Subsetting can lose CRS — re-enforce
    if (is.na(sf::st_crs(gdf))) sf::st_crs(gdf) <- 4326L
    
    if (nrow(gdf) == 0) return(NULL)
    return(gdf)
  }, error = function(e) {
    warning(paste("Failed to load boundary:", e$message))
    return(NULL)
  })
}

#' Load a building footprint layer
#' @param file_path Path to GeoJSON/GPKG/SHP file
#' @return sf object or NULL
load_building_layer <- function(file_path) {
  if (is.null(file_path) || file_path == "") return(NULL)
  
  tryCatch({
    gdf <- sf::st_read(file_path, quiet = TRUE)
    if (nrow(gdf) == 0) return(NULL)
    
    gdf <- .safe_transform(gdf, 4326L)

    geom_types <- sf::st_geometry_type(gdf)
    gdf <- gdf[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
    if (is.na(sf::st_crs(gdf))) sf::st_crs(gdf) <- 4326L
    
    if (nrow(gdf) == 0) return(NULL)
    return(gdf)
  }, error = function(e) {
    warning(paste("Failed to load buildings:", e$message))
    return(NULL)
  })
}

# ============================================================================
# Building Cache — tile-based SQLite storage
# Each 0.25° × 0.25° tile (~28 km) is fetched once from OSM and stored
# permanently. Subsequent analyses for the same area read from local DB.
# ============================================================================

BUILDING_TILE_SIZE <- 0.25   # degrees

# Internal: build a tile ID from its lower-left corner
.tile_id <- function(lat_min, lon_min) {
  sprintf("%.2f_%.2f", lat_min, lon_min)
}

# Internal: parse a tile ID back to lat_min / lon_min
.parse_tile <- function(tile_id) {
  parts <- strsplit(tile_id, "_", fixed = TRUE)[[1]]
  list(lat_min = as.numeric(parts[1]), lon_min = as.numeric(parts[2]))
}

# Compute all tile IDs that cover a bounding box
get_needed_tiles <- function(lat_min, lat_max, lon_min, lon_max) {
  ts <- BUILDING_TILE_SIZE
  lat_starts <- seq(floor(lat_min / ts) * ts, floor(lat_max / ts) * ts, by = ts)
  lon_starts <- seq(floor(lon_min / ts) * ts, floor(lon_max / ts) * ts, by = ts)
  as.vector(outer(lat_starts, lon_starts,
                  FUN = function(la, lo) .tile_id(la, lo)))
}

#' Create the building cache database (idempotent)
.open_bldg_con <- function(db_path, flags = RSQLite::SQLITE_RW) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db_path, flags = flags)
  is_cloud <- identical(Sys.getenv("R_CONFIG_ACTIVE"), "shinyapps") ||
              grepl("^/srv/connect", getwd())
  if (!is_cloud) {
    RSQLite::dbExecute(con, "PRAGMA cache_size = -102400")
    RSQLite::dbExecute(con, "PRAGMA mmap_size = 4294967296")
  } else {
    RSQLite::dbExecute(con, "PRAGMA cache_size = -8192")   # 32 MB on cloud
    RSQLite::dbExecute(con, "PRAGMA mmap_size = 0")
  }
  RSQLite::dbExecute(con, "PRAGMA temp_store = MEMORY")
  con
}

init_building_cache <- function(db_path) {
  con <- .open_bldg_con(db_path, flags = RSQLite::SQLITE_RWC)
  on.exit(RSQLite::dbDisconnect(con))

  RSQLite::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS cache_tiles (
      tile_id     TEXT PRIMARY KEY,
      lat_min     REAL,
      lat_max     REAL,
      lon_min     REAL,
      lon_max     REAL,
      fetched_at  TEXT,
      n_buildings INTEGER,
      status      TEXT DEFAULT 'ok'
    )")

  RSQLite::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS buildings (
      id           INTEGER PRIMARY KEY AUTOINCREMENT,
      tile_id      TEXT,
      latitude     REAL,
      longitude    REAL,
      geometry_wkt TEXT,
      area_m2      REAL,
      confidence   REAL,
      osm_id       TEXT
    )")

  RSQLite::dbExecute(con,
    "CREATE INDEX IF NOT EXISTS idx_bldg_tile ON buildings(tile_id)")
  RSQLite::dbExecute(con,
    "CREATE INDEX IF NOT EXISTS idx_bldg_tile_lat_lon ON buildings(tile_id, latitude, longitude)")

  invisible(db_path)
}

#' Return cache statistics (tile count, building count)
get_cache_stats <- function(db_path) {
  if (!file.exists(db_path))
    return(list(n_tiles = 0, n_buildings = 0, n_failed = 0))
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(RSQLite::dbDisconnect(con))
  list(
    n_tiles    = RSQLite::dbGetQuery(con, "SELECT COUNT(*) AS n FROM cache_tiles WHERE status = 'ok'")$n,
    n_buildings= RSQLite::dbGetQuery(con, "SELECT COUNT(*) AS n FROM buildings")$n,
    n_failed   = RSQLite::dbGetQuery(con, "SELECT COUNT(*) AS n FROM cache_tiles WHERE status = 'failed'")$n
  )
}

#' Drop all cached buildings (admin use)
clear_building_cache <- function(db_path) {
  if (!file.exists(db_path)) return(invisible(NULL))
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(RSQLite::dbDisconnect(con))
  RSQLite::dbExecute(con, "DELETE FROM buildings")
  RSQLite::dbExecute(con, "DELETE FROM cache_tiles")
  invisible(NULL)
}

#' Import buildings from a local file (GeoJSON / GPKG / SHP) into the cache.
#'
#' Tiles the geometries into 0.25° grid cells and inserts them permanently.
#' Any tile already in the cache is replaced (re-import is safe).
#'
#' @param db_path     Path to buildings_cache.sqlite
#' @param file_path   Path to the input building file
#' @param progress_cb Optional function(n_done, n_total) for progress updates
#' @return Number of buildings imported (integer)
import_buildings_to_cache <- function(db_path, file_path, progress_cb = NULL) {
  message("Reading building file: ", file_path)
  gdf <- sf::st_read(file_path, quiet = TRUE)
  if (nrow(gdf) == 0) return(0L)

  gdf <- .safe_transform(gdf, EPSG_WGS84)
  gdf <- sf::st_make_valid(gdf)
  if (is.na(sf::st_crs(gdf))) sf::st_crs(gdf) <- EPSG_WGS84

  geom_types <- as.character(sf::st_geometry_type(gdf))
  gdf <- gdf[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
  # Subsetting can lose CRS — re-enforce
  if (is.na(sf::st_crs(gdf))) sf::st_crs(gdf) <- EPSG_WGS84
  if (nrow(gdf) == 0) return(0L)

  # Use centroids to assign each building to a tile
  # Note: st_centroid(st_geometry(gdf)) might lose CRS, so use gdf directly
  cents <- tryCatch(
    sf::st_coordinates(sf::st_centroid(gdf)),
    error = function(e) {
      # Fallback: extract coordinates manually
      sf::st_coordinates(sf::st_centroid(sf::st_geometry(gdf)))
    }
  )
  ts       <- BUILDING_TILE_SIZE
  tile_ids <- mapply(function(lat, lon)
    .tile_id(floor(lat / ts) * ts, floor(lon / ts) * ts),
    cents[, 2], cents[, 1])

  wkt_vec <- sf::st_as_text(sf::st_geometry(gdf))

  bldg_df <- data.frame(
    tile_id      = tile_ids,
    osm_id       = NA_character_,
    geometry_wkt = wkt_vec,
    stringsAsFactors = FALSE
  )

  if (!file.exists(db_path)) init_building_cache(db_path)
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(RSQLite::dbDisconnect(con))

  # Insert in batches of 2000
  n_total    <- nrow(bldg_df)
  batch_size <- 2000L
  n_done     <- 0L

  for (i in seq(1, n_total, by = batch_size)) {
    batch  <- bldg_df[i:min(i + batch_size - 1L, n_total), ]
    RSQLite::dbWriteTable(con, "buildings", batch, append = TRUE)
    n_done <- n_done + nrow(batch)
    if (!is.null(progress_cb)) progress_cb(n_done, n_total)
  }

  # Update cache_tiles for every tile that now has buildings
  tile_counts   <- table(bldg_df$tile_id)
  now_str       <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  for (tile in names(tile_counts)) {
    parts <- .parse_tile(tile)
    RSQLite::dbExecute(con,
      "INSERT OR REPLACE INTO cache_tiles
         (tile_id, lat_min, lat_max, lon_min, lon_max, fetched_at, n_buildings, status)
         VALUES (?, ?, ?, ?, ?, ?, ?, 'ok')",
      params = list(tile,
                    parts$lat_min, parts$lat_min + BUILDING_TILE_SIZE,
                    parts$lon_min, parts$lon_min + BUILDING_TILE_SIZE,
                    now_str, as.integer(tile_counts[[tile]])))
  }

  message(sprintf("Import complete: %d buildings across %d tiles",
                  n_done, length(tile_counts)))
  n_done
}

# Internal: fetch one tile from OSM and insert into an open DB connection
.fetch_one_tile <- function(con, tile_id, lat_min, lon_min) {
  if (!requireNamespace("osmdata", quietly = TRUE)) {
    warning("osmdata package required for OSM fetch.")
    return(FALSE)
  }
  lat_max  <- lat_min + BUILDING_TILE_SIZE
  lon_max  <- lon_min + BUILDING_TILE_SIZE
  now_str  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  tryCatch({
    q <- osmdata::opq(bbox = c(lon_min, lat_min, lon_max, lat_max),
                      timeout = 60) |>
      osmdata::add_osm_feature(key = "building")
    result    <- osmdata::osmdata_sf(q)
    buildings <- result$osm_polygons
    if (is.null(buildings) || nrow(buildings) == 0)
      buildings <- result$osm_multipolygons

    n_bldg <- 0L
    if (!is.null(buildings) && nrow(buildings) > 0) {
      buildings <- .safe_transform(buildings, EPSG_WGS84)
      buildings <- sf::st_make_valid(buildings)
      if (is.na(sf::st_crs(buildings))) sf::st_crs(buildings) <- EPSG_WGS84
      wkt_vec  <- sf::st_as_text(sf::st_geometry(buildings))
      osm_ids  <- if ("osm_id" %in% names(buildings))
                    as.character(buildings$osm_id)
                  else rep(NA_character_, nrow(buildings))
      bldg_df  <- data.frame(tile_id = tile_id, osm_id = osm_ids,
                              geometry_wkt = wkt_vec,
                              stringsAsFactors = FALSE)
      RSQLite::dbWriteTable(con, "buildings", bldg_df, append = TRUE)
      n_bldg <- nrow(bldg_df)
    }

    RSQLite::dbExecute(con,
      "INSERT OR REPLACE INTO cache_tiles
         (tile_id, lat_min, lat_max, lon_min, lon_max, fetched_at, n_buildings, status)
         VALUES (?, ?, ?, ?, ?, ?, ?, 'ok')",
      params = list(tile_id, lat_min, lat_max, lon_min, lon_max, now_str, n_bldg))

    message(sprintf("  Cached tile %s: %d buildings", tile_id, n_bldg))
    TRUE

  }, error = function(e) {
    warning(sprintf("OSM fetch failed for tile %s: %s", tile_id, e$message))
    RSQLite::dbExecute(con,
      "INSERT OR REPLACE INTO cache_tiles
         (tile_id, lat_min, lat_max, lon_min, lon_max, fetched_at, n_buildings, status)
         VALUES (?, ?, ?, ?, ?, ?, 0, 'failed')",
      params = list(tile_id, lat_min, lat_max, lon_min, lon_max, now_str))
    FALSE
  })
}

#' Get buildings for a bounding box — using the local cache.
#'
#' Missing tiles are fetched from OSM (if fetch_missing = TRUE) and stored
#' permanently so subsequent calls for the same area are instant.
#'
#' @param db_path       Path to buildings_cache.sqlite
#' @param lat_min/max   Bounding box latitudes
#' @param lon_min/max   Bounding box longitudes
#' @param fetch_missing Fetch tiles not yet in the cache from OSM?
#' @param retry_hours   Re-attempt previously failed tiles after this many hours
#' @param progress_cb   Optional function(i, total, tile_id) for progress updates
#' @return sf object with building polygons, or NULL
get_buildings_from_cache <- function(db_path, lat_min, lat_max, lon_min, lon_max,
                                     fetch_missing = TRUE,
                                     retry_hours   = 24,
                                     progress_cb   = NULL) {
  if (!file.exists(db_path)) init_building_cache(db_path)

  needed <- get_needed_tiles(lat_min, lat_max, lon_min, lon_max)
  if (length(needed) == 0) return(NULL)

  ph  <- paste(rep("?", length(needed)), collapse = ",")
  con <- .open_bldg_con(
    db_path,
    flags = if (isTRUE(fetch_missing)) RSQLite::SQLITE_RWC else RSQLite::SQLITE_RO
  )
  on.exit(RSQLite::dbDisconnect(con))

  cached_df <- RSQLite::dbGetQuery(con,
    sprintf("SELECT tile_id, status, fetched_at FROM cache_tiles WHERE tile_id IN (%s)", ph),
    params = as.list(needed))

  cached_ok <- cached_df$tile_id[cached_df$status == "ok"]

  # Re-try tiles that failed more than retry_hours ago
  retry_cutoff <- format(Sys.time() - retry_hours * 3600, "%Y-%m-%d %H:%M:%S")
  retry_tiles  <- cached_df$tile_id[
    cached_df$status == "failed" & cached_df$fetched_at < retry_cutoff]

  missing_tiles <- c(setdiff(needed, cached_df$tile_id), retry_tiles)

  if (fetch_missing && length(missing_tiles) > 0 &&
      requireNamespace("osmdata", quietly = TRUE)) {
    message(sprintf("Building cache: fetching %d new tile(s) from OSM...",
                    length(missing_tiles)))
    for (i in seq_along(missing_tiles)) {
      tile  <- missing_tiles[[i]]
      parts <- .parse_tile(tile)
      if (!is.null(progress_cb)) progress_cb(i, length(missing_tiles), tile)
      .fetch_one_tile(con, tile, parts$lat_min, parts$lon_min)
    }
  } else if (length(missing_tiles) > 0) {
    message(sprintf("Building cache: %d tile(s) not cached (fetch_missing = FALSE).",
                    length(missing_tiles)))
  }

  # Query all buildings in the requested tiles/extent. Do not cap this for the
  # CLI/batch workflow, because under-sampling buildings biases validation
  # toward Outdoor.
  rows <- RSQLite::dbGetQuery(con,
    sprintf(
      "SELECT geometry_wkt FROM buildings
       WHERE tile_id IN (%s)
         AND (
           latitude IS NULL OR longitude IS NULL OR
           (latitude BETWEEN ? AND ? AND longitude BETWEEN ? AND ?)
         )",
      ph),
    params = c(as.list(needed),
               list(lat_min, lat_max, lon_min, lon_max)))

  if (nrow(rows) == 0) return(NULL)

  geoms <- sf::st_as_sfc(rows$geometry_wkt, crs = EPSG_WGS84)
  # Ensure CRS is properly set on sfc before creating sf
  if (!is.null(geoms) && length(geoms) > 0) {
    geoms <- .ensure_crs(geoms, EPSG_WGS84)
  }
  sf::st_sf(geometry = geoms)
}

# ============================================================================
# OSM Building Fetching (legacy — now wraps the cache when db_path provided)
# ============================================================================

#' Fetch buildings for a data.table of GPS points (uses cache when available)
#' @param dt        data.table with lat/lon columns
#' @param buffer_km Buffer in kilometres around GPS points
#' @param cache_db  Path to buildings_cache.sqlite (kept for signature compat, not used)
#' @param preloaded_csv  Path to pre-loaded buildings WKT CSV (buildings_preloaded.csv)
#' @return sf object with building polygons clipped to data extent, or NULL
fetch_osm_buildings <- function(dt, buffer_km = 0.5, cache_db = NULL,
                                preloaded_csv = "buildings_preloaded.csv") {
  if (nrow(dt) == 0) return(NULL)

  lat_buffer <- buffer_km / 111
  mid_lat    <- (min(dt$lat) + max(dt$lat)) / 2
  lon_buffer <- buffer_km / max(1e-6, 111 * cos(mid_lat * pi / 180))

  lat_min <- min(dt$lat) - lat_buffer
  lat_max <- max(dt$lat) + lat_buffer
  lon_min <- min(dt$lon) - lon_buffer
  lon_max <- max(dt$lon) + lon_buffer

  # ── Primary source: pre-loaded WKT CSV ──────────────────────────────────────
  # Try the csv in the working directory or next to the app.R file
  csv_candidates <- c(
    preloaded_csv,
    file.path(dirname(sys.frame(1)$ofile %||% "."), "buildings_preloaded.csv"),
    "buildings_preloaded.csv"
  )
  csv_path <- Filter(file.exists, csv_candidates)[1]

  if (!is.na(csv_path) && !is.null(csv_path) && file.exists(csv_path)) {
    tryCatch({
      bdt <- data.table::fread(csv_path, showProgress = FALSE)
      if (!"geometry_wkt" %in% names(bdt) || nrow(bdt) == 0) stop("empty or wrong columns")

      # Filter to bounding box before parsing WKT (fast row filter on tile_id via lat/lon prefix)
      geoms <- sf::st_as_sfc(bdt$geometry_wkt, crs = EPSG_WGS84)
      geoms <- .ensure_crs(geoms, EPSG_WGS84)
      bldg  <- sf::st_sf(id = bdt$id, geometry = geoms)
      bldg  <- suppressWarnings(sf::st_make_valid(bldg))
      bldg  <- bldg[!sf::st_is_empty(bldg), , drop = FALSE]

      # Clip to data extent
      bbox_vec <- c(xmin = lon_min, ymin = lat_min, xmax = lon_max, ymax = lat_max)
      bbox_bbox <- sf::st_bbox(bbox_vec, crs = EPSG_WGS84)
      bbox_poly <- sf::st_as_sfc(bbox_bbox)
      # Explicitly set CRS on the sfc object after creation
      sf::st_crs(bbox_poly) <- EPSG_WGS84
      clipped <- tryCatch(
        suppressWarnings(sf::st_filter(bldg, bbox_poly)),
        error = function(e) {
          hits <- sf::st_intersects(bldg, bbox_poly, sparse = FALSE)[, 1]
          bldg[hits, , drop = FALSE]
        }
      )

      if (!is.null(clipped) && nrow(clipped) > 0) {
        message(sprintf("[SPECTRA] Loaded %d buildings from pre-built database", nrow(clipped)))
        return(clipped)
      }
    }, error = function(e) {
      message("[SPECTRA] Pre-loaded buildings CSV error: ", e$message, " — trying SQLite cache")
    })
  }

  # ── Fallback: SQLite tile cache ──────────────────────────────────────────────
  if (!is.null(cache_db) && file.exists(cache_db)) {
    tryCatch({
      result <- get_buildings_from_cache(cache_db,
                                      lat_min = lat_min, lat_max = lat_max,
                                      lon_min = lon_min, lon_max = lon_max,
                                      fetch_missing = FALSE)
      if (!is.null(result)) result <- .ensure_crs(result, EPSG_WGS84)
      return(result)
    }, error = function(e) {
      message("[SPECTRA] SQLite cache error: ", e$message)
    })
  }

  message("[SPECTRA] No buildings found — returning NULL")
  NULL
}

# ============================================================================
# Boundary Filtering
# ============================================================================

#' Filter trajectory points to within a boundary polygon
#' @param dt data.table with lat/lon
#' @param boundary_sf sf object with boundary polygons
#' @return filtered data.table
apply_boundary_filter <- function(dt, boundary_sf) {
  if (is.null(boundary_sf)) return(dt)
  
  pts <- sf::st_as_sf(dt, coords = c("lon", "lat"), crs = 4326L, remove = FALSE)
  boundary_sf <- .safe_transform(boundary_sf, 4326L)
  if (is.na(sf::st_crs(pts))) sf::st_crs(pts) <- 4326L
  
  joined <- sf::st_join(pts, boundary_sf, join = sf::st_within, left = FALSE)
  
  result <- as.data.table(sf::st_drop_geometry(joined))
  return(result)
}

# ============================================================================
# Haversine Distance
# ============================================================================

#' Calculate haversine distance between consecutive points
#' @param lat1,lon1,lat2,lon2 Numeric vectors
#' @return Distance in kilometres
haversine_distance <- function(lat1, lon1, lat2, lon2) {
  R <- 6371.0  # Earth's radius in km
  
  lat1_r <- lat1 * pi / 180
  lon1_r <- lon1 * pi / 180
  lat2_r <- lat2 * pi / 180
  lon2_r <- lon2 * pi / 180
  
  dlat <- lat2_r - lat1_r
  dlon <- lon2_r - lon1_r
  
  a <- sin(dlat / 2)^2 + cos(lat1_r) * cos(lat2_r) * sin(dlon / 2)^2
  c_val <- 2 * asin(sqrt(a))
  
  return(R * c_val)
}

# ============================================================================
# Indoor/Outdoor Classification
# ============================================================================

#' Classify trajectory points as indoor or outdoor
#' @param dt data.table with lat/lon (and optionally speed_kmh)
#' @param buildings_sf sf object with building polygons
#' @param buffer_m Buffer in metres around buildings
#' @param displacement_threshold_kmh Displacement threshold (km/h). 
#'        Primary classifier: speed >= this => OUTDOOR. Default 6 km/h (100 m/min)
#' @return data.table with is_indoor and indoor_confidence columns added
classify_indoor <- function(dt, buildings_sf, buffer_m = 0, 
                            displacement_threshold_kmh = 6) {
  result <- copy(dt)
  
  # Step 1: PRIMARY — Displacement-based classification
  if ("speed_kmh" %in% names(result)) {
    # Speed >= threshold => definitely outdoor
    result[, is_indoor := !(!is.na(speed_kmh) & speed_kmh >= displacement_threshold_kmh)]
  } else {
    # No speed data — default to indoor (will be corrected by buildings below)
    result[, is_indoor := TRUE]
  }
  
  # Step 2: SECONDARY — Building detection validates/corrects
  if (!is.null(buildings_sf) && nrow(buildings_sf) > 0) {
    # Ensure buildings have valid CRS before any spatial operations
    buildings_sf <- .ensure_crs(buildings_sf)
    buildings_sf <- .to_wgs84(buildings_sf)

    bldg_wgs84 <- buildings_sf
    bldg       <- .safe_transform(bldg_wgs84, EPSG_METRIC, EPSG_WGS84)
    if (is.na(sf::st_crs(bldg))) sf::st_crs(bldg) <- EPSG_METRIC

    if (buffer_m > 0) {
      bldg <- sf::st_buffer(bldg, dist = buffer_m)
      if (is.na(sf::st_crs(bldg))) sf::st_crs(bldg) <- EPSG_METRIC
    }

    # Pre-filter: only convert points inside the buildings bounding box to sf
    bb <- sf::st_bbox(bldg_wgs84)
    in_bbox <- result$lon >= bb["xmin"] & result$lon <= bb["xmax"] &
               result$lat >= bb["ymin"] & result$lat <= bb["ymax"]

    in_building <- rep(FALSE, nrow(result))
    if (any(in_bbox)) {
      pts_sub <- sf::st_as_sf(result[in_bbox], coords = c("lon", "lat"),
                               crs = EPSG_WGS84, remove = FALSE)
      pts_sub <- .safe_transform(pts_sub, EPSG_METRIC, EPSG_WGS84)

      # Chunk intersection to avoid peak-memory blowup on large point sets
      CHUNK_SIZE <- 50000L
      n_sub      <- nrow(pts_sub)
      in_bldg_sub <- logical(n_sub)
      for (cs in seq(1L, n_sub, by = CHUNK_SIZE)) {
        ce <- min(cs + CHUNK_SIZE - 1L, n_sub)
        ci <- sf::st_intersects(pts_sub[cs:ce, ], bldg, sparse = TRUE)
        in_bldg_sub[cs:ce] <- lengths(ci) > 0
      }
      in_building[in_bbox] <- in_bldg_sub
    }

    # Validation logic:
    # - If classified as INDOOR but NOT in building => correct to OUTDOOR (false positive)
    # - If classified as OUTDOOR and IS in building => keep OUTDOOR (may be moving indoors)
    result[is_indoor == TRUE & !in_building, is_indoor := FALSE]
  }
  
  # Confidence: high if displacement-based, lower if only tentative
  result[, indoor_confidence := as.numeric(is_indoor)]
  if ("speed_kmh" %in% names(result)) {
    # Speed-based classification is high confidence
    result[!is.na(speed_kmh) & speed_kmh >= displacement_threshold_kmh, indoor_confidence := 0.0]
  }
  
  return(result)
}

# ============================================================================
# Dwell Time Calculation
# ============================================================================

#' Calculate dwell time, speed, and filter teleports
#' @param dt data.table with lat/lon/PID/date
#' @param max_dwell_min Maximum dwell time in minutes
#' @param speed_threshold_kmh Speed threshold for teleport removal
#' @param remove_teleports Whether to remove teleporting points
#' @return data.table with dwell_min, speed_kmh, date_only columns
calculate_dwell_time <- function(dt, max_dwell_min = 10,
                                 speed_threshold_kmh = 120,
                                 remove_teleports = TRUE) {
  result <- copy(dt)
  result[, date_only := as.Date(date)]
  
  # Sort by PID and date
  setorder(result, PID, date)
  
  # Next point coordinates within same PID
  result[, lat_next := shift(lat, n = 1, type = "lead"), by = PID]
  result[, lon_next := shift(lon, n = 1, type = "lead"), by = PID]
  result[, date_next := shift(date, n = 1, type = "lead"), by = PID]
  result[, date_only_next := shift(date_only, n = 1, type = "lead"), by = PID]
  
  # Distance in km
  result[, distance_km := haversine_distance(
    lat, lon,
    ifelse(is.na(lat_next), lat, lat_next),
    ifelse(is.na(lon_next), lon, lon_next)
  )]
  
  # Time difference in hours
  result[, time_diff_hours := as.numeric(difftime(date_next, date, units = "hours"))]
  result[time_diff_hours <= 0, time_diff_hours := NA_real_]
  
  # Speed in km/h
  result[, speed_next_kmh := distance_km / time_diff_hours]
  result[, speed_prev_kmh := shift(speed_next_kmh, n = 1, type = "lag"), by = PID]
  result[, speed_kmh := pmax(speed_next_kmh, speed_prev_kmh, na.rm = TRUE)]
  result[!is.finite(speed_kmh), speed_kmh := NA_real_]
  
  # Remove teleports
  if (remove_teleports) {
    result <- result[is.na(speed_next_kmh) | speed_next_kmh <= speed_threshold_kmh]
  }
  
  # Dwell time in minutes
  result[, dwell_min := as.numeric(difftime(date_next, date, units = "mins"))]
  
  # Invalidate: NA, negative, or cross-day
  result[is.na(dwell_min) | dwell_min <= 0 | date_only != date_only_next,
         dwell_min := 0]
  
  # Clip to maximum
  result[dwell_min > max_dwell_min, dwell_min := max_dwell_min]

  # Remove temporary navigation columns — not needed downstream, frees memory
  result[, c("lat_next", "lon_next", "date_next", "date_only_next",
             "distance_km", "time_diff_hours", "speed_next_kmh",
             "speed_prev_kmh") := NULL]

  return(result)
}

# ============================================================================
# Grid Generation & Aggregation
# ============================================================================

#' Create a non-overlapping fishnet grid covering study area
#' @param dt data.table with lat/lon (in WGS84)
#' @param cell_size_m Grid cell size in metres
#' @return sf object with fishnet polygons (one row per cell), no overlaps
create_fishnet_grid <- function(dt, cell_size_m = 100) {
  if (nrow(dt) == 0) stop("No data to create grid")
  
  # Project to metric — ensure valid CRS
  dt_clean <- dt[!is.na(lon) & !is.na(lat)]
  if (nrow(dt_clean) == 0) stop("No valid coordinates for grid")
  
  pts <- sf::st_as_sf(dt_clean, coords = c("lon", "lat"), crs = EPSG_WGS84, remove = FALSE)
  pts <- .safe_transform(pts, EPSG_METRIC, EPSG_WGS84)
  bbox <- sf::st_bbox(pts)
  
  # Snap bounds to grid multiples
  xmin <- floor(bbox["xmin"] / cell_size_m) * cell_size_m
  ymin <- floor(bbox["ymin"] / cell_size_m) * cell_size_m
  xmax <- ceiling(bbox["xmax"] / cell_size_m) * cell_size_m
  ymax <- ceiling(bbox["ymax"] / cell_size_m) * cell_size_m
  
  # Create grid sequences
  x_seq <- seq(xmin, xmax, by = cell_size_m)
  y_seq <- seq(ymin, ymax, by = cell_size_m)
  
  # Create cells using expand.grid for simpler handling
  nx <- length(x_seq) - 1
  ny <- length(y_seq) - 1
  
  cell_list <- list()
  cell_id <- 1
  
  for (i in 1:nx) {
    for (j in 1:ny) {
      # Create polygon corners
      corners <- matrix(c(
        x_seq[i],     y_seq[j],
        x_seq[i+1],   y_seq[j],
        x_seq[i+1],   y_seq[j+1],
        x_seq[i],     y_seq[j+1],
        x_seq[i],     y_seq[j]
      ), ncol = 2, byrow = TRUE)
      
      cell_list[[cell_id]] <- sf::st_polygon(list(corners))
      cell_id <- cell_id + 1
    }
  }
  
  # Create sfc with proper CRS
  grid_sfc <- sf::st_sfc(cell_list, crs = EPSG_METRIC)
  
  # Create sf object directly
  grid_sf <- sf::st_sf(
    cell_id = 1:length(grid_sfc),
    geometry = grid_sfc
  )
  # Verify CRS is set
  if (is.na(sf::st_crs(grid_sf))) sf::st_crs(grid_sf) <- EPSG_METRIC
  
  # Get centroids in metric — st_centroid can lose CRS, so re-enforce
  centroids <- sf::st_centroid(grid_sf)
  if (is.na(sf::st_crs(centroids))) sf::st_crs(centroids) <- EPSG_METRIC
  coords_metric <- sf::st_coordinates(centroids)
  
  centroids_wgs84 <- .safe_transform(centroids, EPSG_WGS84, EPSG_METRIC)
  coords_wgs <- sf::st_coordinates(centroids_wgs84)
  
  # Add lat/lon to the grid
  grid_sf$lon <- coords_wgs[, 1]
  grid_sf$lat <- coords_wgs[, 2]
  
  # Ensure geometry retains metric CRS for spatial operations
  if (is.na(sf::st_crs(grid_sf))) sf::st_crs(grid_sf) <- EPSG_METRIC
  
  return(grid_sf)
}

#' Project to metric CRS and aggregate to grid cells
#' @param dt data.table with lat/lon/PID/date_only/pm25/no2/dwell_min/is_indoor
#' @param cell_size_m Grid cell size in metres
#' @return list with person_day_cell (data.table) and grid_params (list)
aggregate_to_grid <- function(dt, cell_size_m = 1000) {
  
  # Filter valid coordinates first
  dt_valid <- dt[!is.na(lon) & !is.na(lat) & !is.na(pm25)]
  if (nrow(dt_valid) == 0) stop("No valid coordinate/exposure data for gridding")
  
  # Project to metric — sf_project avoids allocating a full sf geometry object
  coords <- tryCatch(
    sf::sf_project(from = sf::st_crs(EPSG_WGS84), to = sf::st_crs(EPSG_METRIC),
                   pts  = cbind(dt_valid$lon, dt_valid$lat)),
    error = function(e) {
      pts <- sf::st_as_sf(dt_valid, coords = c("lon", "lat"), crs = EPSG_WGS84, remove = FALSE)
      pts <- .safe_transform(pts, EPSG_METRIC, EPSG_WGS84)
      sf::st_coordinates(pts)
    }
  )

  result <- copy(dt_valid)
  result[, x := coords[, 1]]
  result[, y := coords[, 2]]
  
  padding <- 5000
  xmin <- min(result$x) - padding
  xmax <- max(result$x) + padding
  ymin <- min(result$y) - padding
  ymax <- max(result$y) + padding
  
  ncols <- ceiling((xmax - xmin) / cell_size_m)
  nrows <- ceiling((ymax - ymin) / cell_size_m)
  
  result[, col_idx := as.integer(floor((x - xmin) / cell_size_m))]
  result[, row_idx := as.integer(floor((y - ymin) / cell_size_m))]
  result[, cell_id := row_idx * ncols + col_idx]
  
  # Filter in-bounds
  result <- result[col_idx >= 0 & col_idx < ncols & row_idx >= 0 & row_idx < nrows]
  
  # Weighted values
  result[, pm25_w := pm25 * dwell_min]
  result[, no2_w := no2 * dwell_min]
  
  # Aggregate per person-day-cell
  # country is optional — include it only when it exists
  by_cols <- c("PID", "date_only", "cell_id")
  if ("country" %in% names(result)) by_cols <- c(by_cols, "country")

  pdc <- result[, .(
    dwell_sum_min = sum(dwell_min, na.rm = TRUE),
    pm25_w_sum    = sum(pm25_w,   na.rm = TRUE),
    no2_w_sum     = sum(no2_w,    na.rm = TRUE),
    n_points      = .N,
    n_indoor      = sum(is_indoor, na.rm = TRUE),
    x_mean        = mean(x,       na.rm = TRUE),
    y_mean        = mean(y,       na.rm = TRUE),
    row_idx       = first(row_idx),
    col_idx       = first(col_idx),
    time_start    = min(date,     na.rm = TRUE),
    time_end      = max(date,     na.rm = TRUE)
  ), by = by_cols]
  
  # Time-weighted means
  pdc[, pm25_mean_tw := fifelse(dwell_sum_min > 0, pm25_w_sum / dwell_sum_min, NA_real_)]
  pdc[, no2_mean_tw := fifelse(dwell_sum_min > 0, no2_w_sum / dwell_sum_min, NA_real_)]
  pdc[, pct_indoor := (n_indoor / n_points) * 100]
  
  # Back-project centroids to WGS84
  coords_wgs <- tryCatch(
    sf::sf_project(from = sf::st_crs(EPSG_METRIC), to = sf::st_crs(EPSG_WGS84),
                   pts  = cbind(pdc$x_mean, pdc$y_mean)),
    error = function(e) {
      csf <- sf::st_as_sf(data.frame(x = pdc$x_mean, y = pdc$y_mean),
                          coords = c("x", "y"), crs = EPSG_METRIC)
      if (is.na(sf::st_crs(csf))) sf::st_crs(csf) <- EPSG_METRIC
      sf::st_coordinates(.safe_transform(csf, EPSG_WGS84, EPSG_METRIC))
    }
  )
  pdc[, lon := coords_wgs[, 1]]
  pdc[, lat := coords_wgs[, 2]]

  # Release x/y and weighted intermediate columns from the working copy
  result[, c("x", "y", "pm25_w", "no2_w", "col_idx", "row_idx", "cell_id") := NULL]
  gc(verbose = FALSE)
  
  grid_params <- list(
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    ncols = ncols, nrows = nrows,
    cell_size = cell_size_m
  )

  # ── Full continuous grid — ALL cells (observed + empty) ─────────────────
  # Aggregate pdc across persons and dates to get one row per cell
  has_country <- "country" %in% names(pdc)

  obs_agg <- pdc[, .(
    dwell           = sum(dwell_sum_min,                          na.rm = TRUE),
    pm25_w_sum      = sum(pm25_mean_tw * dwell_sum_min, na.rm = TRUE),
    no2_w_sum       = sum(no2_mean_tw  * dwell_sum_min, na.rm = TRUE),
    total_dwell_sum = sum(dwell_sum_min,                          na.rm = TRUE),
    indoor_pct      = mean(pct_indoor,                            na.rm = TRUE),
    country         = if (has_country) first(country) else NA_character_
  ), by = .(cell_id, col_idx, row_idx)]
  if (!has_country) obs_agg[, country := NULL]
  
  # Calculate final weighted means
  obs_agg[, pm25 := fifelse(total_dwell_sum > 0, pm25_w_sum / total_dwell_sum, NA_real_)]
  obs_agg[, no2  := fifelse(total_dwell_sum > 0, no2_w_sum / total_dwell_sum, NA_real_)]
  obs_agg[, c("pm25_w_sum", "no2_w_sum", "total_dwell_sum") := NULL]

  # Build every cell in the extent
  all_ci  <- rep(0L:(ncols - 1L), times = nrows)
  all_ri  <- rep(0L:(nrows - 1L), each  = ncols)
  full_grid <- data.table::data.table(
    col_idx = all_ci,
    row_idx = all_ri,
    cell_id = all_ri * ncols + all_ci
  )

  # Merge observed values (NA where no data) — include country when available
  obs_merge_cols <- c("cell_id", "dwell", "pm25", "no2", "indoor_pct")
  if (has_country) obs_merge_cols <- c(obs_merge_cols, "country")
  full_grid <- merge(full_grid, obs_agg[, obs_merge_cols, with = FALSE],
                     by = "cell_id", all.x = TRUE)

  # is_labeled: 1 = cell has GPS observations (labeled for ML), 0 = empty (unlabeled)
  full_grid[, is_labeled := as.integer(!is.na(dwell))]

  # Project cell centres to WGS84 — sf_project avoids full sf object allocation
  full_grid[, x_ctr := xmin + (col_idx + 0.5) * cell_size_m]
  full_grid[, y_ctr := ymin + (row_idx + 0.5) * cell_size_m]
  fg_wgs <- tryCatch(
    sf::sf_project(from = sf::st_crs(EPSG_METRIC), to = sf::st_crs(EPSG_WGS84),
                   pts  = cbind(full_grid$x_ctr, full_grid$y_ctr)),
    error = function(e) {
      fg_sf <- sf::st_as_sf(data.frame(x = full_grid$x_ctr, y = full_grid$y_ctr),
                             coords = c("x", "y"), crs = EPSG_METRIC)
      if (is.na(sf::st_crs(fg_sf))) sf::st_crs(fg_sf) <- EPSG_METRIC
      sf::st_coordinates(.safe_transform(fg_sf, EPSG_WGS84, EPSG_METRIC))
    }
  )
  full_grid[, lon := fg_wgs[, 1]]
  full_grid[, lat := fg_wgs[, 2]]
  full_grid[, c("x_ctr", "y_ctr") := NULL]

  # ── Date-resolved grid: one row per (date × cell), all PIDs merged ─────────
  # Used for 3-D NetCDF export (lon × lat × time) and GEE date-based extraction
  daily_grid <- pdc[, {
    tot_d <- sum(dwell_sum_min, na.rm = TRUE)
    list(
      dwell          = tot_d,
      pm25           = if (tot_d > 0) sum(pm25_mean_tw * dwell_sum_min, na.rm = TRUE) / tot_d else NA_real_,
      no2            = if (tot_d > 0) sum(no2_mean_tw  * dwell_sum_min, na.rm = TRUE) / tot_d else NA_real_,
      indoor_pct     = if (tot_d > 0) sum(pct_indoor   * dwell_sum_min, na.rm = TRUE) / tot_d else NA_real_,
      lon            = first(lon),
      lat            = first(lat),
      col_idx        = first(col_idx),
      row_idx        = first(row_idx),
      datetime_start = min(time_start, na.rm = TRUE),
      datetime_end   = max(time_end,   na.rm = TRUE)
    )
  }, by = .(date = date_only, cell_id)]
  daily_grid[, grid_id    := cell_id]
  daily_grid[, is_labeled := 1L]
  # Use true grid-cell centres for date-resolved exports. The aggregation above
  # keeps mean observed point coordinates; NetCDF and GEE extraction need the
  # stable cell-centre coordinates used by the full grid.
  centre_lookup <- full_grid[, .(
    cell_id,
    grid_lon     = lon,
    grid_lat     = lat,
    grid_col_idx = col_idx,
    grid_row_idx = row_idx
  )]
  daily_grid <- merge(daily_grid, centre_lookup, by = "cell_id", all.x = TRUE, sort = FALSE)
  daily_grid[!is.na(grid_lon), `:=`(
    lon     = grid_lon,
    lat     = grid_lat,
    col_idx = grid_col_idx,
    row_idx = grid_row_idx
  )]
  daily_grid[, c("grid_lon", "grid_lat", "grid_col_idx", "grid_row_idx") := NULL]
  # ISO 8601 strings for export (preserves timezone info)
  daily_grid[, datetime_start := format(datetime_start, "%Y-%m-%dT%H:%M:%SZ")]
  daily_grid[, datetime_end   := format(datetime_end,   "%Y-%m-%dT%H:%M:%SZ")]
  setcolorder(daily_grid, c("date", "datetime_start", "datetime_end",
                             "grid_id", "cell_id", "lon", "lat",
                             "col_idx", "row_idx", "dwell", "pm25", "no2",
                             "indoor_pct", "is_labeled"))

  return(list(person_day_cell = pdc, grid_params = grid_params,
              full_grid = full_grid, daily_grid = daily_grid))
}

# ============================================================================
# Daily Summaries
# ============================================================================

#' Create daily summaries per person
#' @param pdc data.table from aggregate_to_grid
#' @return data.table with daily summaries
create_daily_summaries <- function(pdc) {
  person_day <- pdc[, .(
    dwell_day_min = sum(dwell_sum_min, na.rm = TRUE),
    pm25_w_day_sum = sum(pm25_w_sum, na.rm = TRUE),
    no2_w_day_sum = sum(no2_w_sum, na.rm = TRUE),
    n_cells_visited = uniqueN(cell_id)
  ), by = .(PID, date_only)]
  
  person_day[, pm25_mean_tw_day := fifelse(
    dwell_day_min > 0, pm25_w_day_sum / dwell_day_min, NA_real_
  )]
  person_day[, no2_mean_tw_day := fifelse(
    dwell_day_min > 0, no2_w_day_sum / dwell_day_min, NA_real_
  )]
  person_day[, coverage_day := dwell_day_min / 1440]
  
  return(person_day)
}

# ============================================================================
# Per-Participant Average Daily Exposure (WHO comparison)
# ============================================================================

#' Compute per-participant average daily exposure from grid-cell daily averages.
#' Averages each participant's daily time-weighted exposures across all days,
#' enabling direct comparison with WHO Air Quality Guidelines (2021).
#'
#' WHO 2021 targets used:
#'   PM2.5: annual mean 5 µg/m³,  24-h mean 15 µg/m³
#'   NO2:   annual mean 10 µg/m³, 24-h mean 25 µg/m³
#'
#' @param person_day data.table from create_daily_summaries
#' @return data.table with one row per PID: averages, maxima, WHO ratios, and exceedance flags
create_participant_summaries <- function(person_day) {
  WHO_PM25_ANNUAL <- 5.0
  WHO_PM25_24H    <- 15.0
  WHO_NO2_ANNUAL  <- 10.0
  WHO_NO2_24H     <- 25.0

  ps <- person_day[, .(
    n_days             = .N,
    pm25_avg_daily     = mean(pm25_mean_tw_day, na.rm = TRUE),
    pm25_max_daily     = max(pm25_mean_tw_day,  na.rm = TRUE),
    no2_avg_daily      = mean(no2_mean_tw_day,  na.rm = TRUE),
    no2_max_daily      = max(no2_mean_tw_day,   na.rm = TRUE),
    dwell_avg_day_min  = mean(dwell_day_min,    na.rm = TRUE),
    n_cells_avg        = mean(n_cells_visited,  na.rm = TRUE)
  ), by = PID]

  ps[, pm25_who_annual_target := WHO_PM25_ANNUAL]
  ps[, pm25_who_24h_target    := WHO_PM25_24H]
  ps[, no2_who_annual_target  := WHO_NO2_ANNUAL]
  ps[, no2_who_24h_target     := WHO_NO2_24H]

  ps[, pm25_pct_of_who_annual := round(pm25_avg_daily / WHO_PM25_ANNUAL * 100, 1)]
  ps[, pm25_pct_of_who_24h    := round(pm25_avg_daily / WHO_PM25_24H    * 100, 1)]
  ps[, no2_pct_of_who_annual  := round(no2_avg_daily  / WHO_NO2_ANNUAL  * 100, 1)]
  ps[, no2_pct_of_who_24h     := round(no2_avg_daily  / WHO_NO2_24H     * 100, 1)]

  ps[, pm25_exceeds_who_annual := fifelse(!is.na(pm25_avg_daily), pm25_avg_daily > WHO_PM25_ANNUAL, NA)]
  ps[, pm25_exceeds_who_24h    := fifelse(!is.na(pm25_max_daily), pm25_max_daily > WHO_PM25_24H,   NA)]
  ps[, no2_exceeds_who_annual  := fifelse(!is.na(no2_avg_daily),  no2_avg_daily  > WHO_NO2_ANNUAL, NA)]
  ps[, no2_exceeds_who_24h     := fifelse(!is.na(no2_max_daily),  no2_max_daily  > WHO_NO2_24H,    NA)]

  ps
}

# ============================================================================
# Confusion Matrix / Quality Assessment
# ============================================================================

#' Calculate confusion matrix metrics
#' @param predicted Logical vector of predicted indoor
#' @param actual Logical vector of actual indoor
#' @return list with tp, fp, tn, fn, accuracy, precision, recall, specificity, f1
confusion_metrics <- function(predicted, actual) {
  tp <- sum(predicted & actual, na.rm = TRUE)
  fp <- sum(predicted & !actual, na.rm = TRUE)
  tn <- sum(!predicted & !actual, na.rm = TRUE)
  fn <- sum(!predicted & actual, na.rm = TRUE)
  
  n <- tp + fp + tn + fn
  accuracy <- if (n > 0) (tp + tn) / n else NA
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA
  }
  
  list(
    tp = tp, fp = fp, tn = tn, fn = fn,
    n_eval = n,
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    specificity = specificity,
    f1 = f1
  )
}

# ============================================================================
# Speed-Based + OSM Validated Indoor/Outdoor Classification
# (from indoor_outdoor_validated.R methodology)
# ============================================================================

#' Validated indoor/outdoor classification pipeline
#' Step 1: Speed-based classification (threshold in m/min)
#' Step 2: OSM building spatial validation
#' Step 3: Identify misclassifications between speed and OSM
#' Step 4: Correct misclassifications using OSM as ground truth
#'
#' @param dt data.table with lat/lon/PID/date (must already have speed_kmh from calculate_dwell_time)
#' @param buildings_sf sf object with building polygons (optional)
#' @param speed_threshold Speed threshold in m/min for indoor classification (default 100)
#' @param building_buffer_m Buffer in metres around buildings
#' @return data.table with speed_class, osm_label, misclass_type, corrected_class, is_indoor
validated_indoor_outdoor <- function(dt, buildings_sf, speed_threshold = 100,
                                     building_buffer_m = 5,
                                     require_building_for_indoor = TRUE,
                                     building_raster = NULL) {
  result <- copy(dt)

  # Convert speed from km/h to m/min for indoor/outdoor classification
  if ("speed_kmh" %in% names(result)) {
    result[, speed_m_min := speed_kmh * 1000 / 60]
  } else {
    result[, speed_m_min := NA_real_]
  }

  # Step 1: Speed-based classification
  result[, speed_class := fifelse(
    is.na(speed_m_min) | speed_m_min < speed_threshold,
    "Indoor", "Outdoor"
  )]

  # Step 2: building-based spatial validation.
  #
  # A raster presence mask is preferred when one is supplied.  OSM footprint
  # coverage is spatially incomplete across these study areas -- 10 of 12 Gambian
  # tiles and the main Kenyan participant tile returned nothing -- and because
  # Indoor requires BOTH low speed and a building, a missing footprint can only
  # ever create a false OUTDOOR.  Google Open Buildings v3 supplies 412k/630k/179k
  # footprints for the same extents against 18.6k/89.6k/33.6k from OSM, rasterised
  # at 5 m so that small dwellings survive (>=1.65 painted pixels per building).
  if (!is.null(building_raster) && file.exists(building_raster)) {
    r <- terra::rast(building_raster)
    pv <- terra::extract(r, cbind(result$lon, result$lat))[[1]]
    in_building <- !is.na(pv) & pv > 0
    result[, in_building := in_building]
    result[, osm_label := fifelse(in_building, "Indoor", "Outdoor")]

    # Step 3: disagreement between displacement and the building layer
    result[, misclass_type := "Correct"]
    result[speed_class == "Outdoor" & osm_label == "Indoor",
           misclass_type := "False Outdoor"]
    result[speed_class == "Indoor" & osm_label == "Outdoor",
           misclass_type := "False Indoor"]

    # Step 4: Indoor requires BOTH a stationary/slow fix and a building beneath it
    result[, corrected_class := fifelse(speed_class == "Indoor" & in_building,
                                        "Indoor", "Outdoor")]
    result[, is_indoor := corrected_class == "Indoor"]
    result[, indoor_confidence := fifelse(misclass_type == "Correct", 1.0, 0.7)]

  } else if (!is.null(buildings_sf) && nrow(buildings_sf) > 0) {
    buildings_sf <- .ensure_crs(buildings_sf)
    buildings_sf <- .to_wgs84(buildings_sf)
    bldg <- sf::st_make_valid(buildings_sf)
    # st_make_valid sometimes loses CRS — re-enforce before transform
    bldg <- .ensure_crs(bldg, EPSG_WGS84)

    if (building_buffer_m > 0) {
      bldg_metric <- .safe_transform(bldg, EPSG_METRIC, EPSG_WGS84)
      bldg_metric <- sf::st_buffer(bldg_metric, dist = building_buffer_m)
      if (is.na(sf::st_crs(bldg_metric))) sf::st_crs(bldg_metric) <- EPSG_METRIC
      bldg <- .safe_transform(bldg_metric, 4326L, EPSG_METRIC)
    }

    # Pre-filter: only convert points inside the buildings bounding box to sf
    bb      <- sf::st_bbox(bldg)
    in_bbox <- result$lon >= bb["xmin"] & result$lon <= bb["xmax"] &
               result$lat >= bb["ymin"] & result$lat <= bb["ymax"]

    in_building <- rep(FALSE, nrow(result))
    if (any(in_bbox)) {
      pts_sub <- sf::st_as_sf(result[in_bbox], coords = c("lon", "lat"),
                               crs = EPSG_WGS84, remove = FALSE)
      # Chunk intersection to avoid peak-memory blowup on large point sets
      CHUNK_SIZE  <- 50000L
      n_sub       <- nrow(pts_sub)
      in_bldg_sub <- logical(n_sub)
      for (cs in seq(1L, n_sub, by = CHUNK_SIZE)) {
        ce <- min(cs + CHUNK_SIZE - 1L, n_sub)
        ci <- sf::st_intersects(pts_sub[cs:ce, ], bldg, sparse = TRUE)
        in_bldg_sub[cs:ce] <- lengths(ci) > 0
      }
      in_building[in_bbox] <- in_bldg_sub
      rm(pts_sub, in_bldg_sub)
    }
    gc(verbose = FALSE)

    result[, in_building := in_building]
    result[, osm_label := fifelse(in_building, "Indoor", "Outdoor")]

    # Step 3: Identify misclassifications between speed and OSM (kept for diagnostics)
    result[, misclass_type := "Correct"]
    result[speed_class == "Outdoor" & osm_label == "Indoor",
           misclass_type := "False Outdoor"]
    result[speed_class == "Indoor" & osm_label == "Outdoor",
           misclass_type := "False Indoor"]

    # Step 4: GIS validation override.
    # Indoor requires BOTH low/stationary speed and a building footprint match.
    # False Indoor points (speed says Indoor, OSM says Outdoor) are corrected
    # to Outdoor for the final exposure split.
    result[, corrected_class := fifelse(speed_class == "Indoor" & in_building,
                                        "Indoor", "Outdoor")]
    result[, is_indoor := corrected_class == "Indoor"]
    result[, indoor_confidence := fifelse(misclass_type == "Correct", 1.0, 0.7)]
  } else {
    # No buildings: fall back to displacement only and make the diagnostics show
    # that no spatial validation happened.
    result[, in_building := NA]
    result[, osm_label := NA_character_]
    result[, misclass_type := NA_character_]
    result[, corrected_class := speed_class]
    result[, is_indoor := corrected_class == "Indoor"]
    result[, indoor_confidence := 0.5]
  }

  return(result)
}

# ============================================================================
# Segmentation Theme Computation
# (from Segmentation.R methodology)
# ============================================================================

#' Compute segmentation themes for trajectory data
#' Adds exposure_type, mobility_theme, temporal_theme, duration_theme,
#' and outdoor-specific themes (movement, context, stress)
#'
#' @param dt data.table with corrected_class, speed_m_min, date, dwell_min
#' @return data.table with additional segmentation theme columns
compute_segmentation_themes <- function(dt) {
  result <- copy(dt)

  if (!"corrected_class" %in% names(result)) {
    result[, corrected_class := fifelse(is_indoor, "Indoor", "Outdoor")]
  }

  result[, exposure_type := factor(corrected_class, levels = c("Indoor", "Outdoor"))]

  # --- Mobility theme based on speed quantiles ---
  if ("speed_m_min" %in% names(result)) {
    sv <- result$speed_m_min[is.finite(result$speed_m_min) & result$speed_m_min > 0]
    if (length(sv) > 10) {
      q25 <- quantile(sv, 0.25, na.rm = TRUE)
      q75 <- quantile(sv, 0.75, na.rm = TRUE)
    } else {
      q25 <- 10; q75 <- 50
    }
    result[, mobility_theme := fcase(
      is.na(speed_m_min) | speed_m_min <= q25, "Stationary",
      speed_m_min <= q75, "Perimeter Move",
      default = "Travel"
    )]
    result[, mobility_theme := factor(mobility_theme,
                                       levels = c("Stationary", "Perimeter Move", "Travel"))]
  }

  # --- Temporal theme ---
  result[, hour_local := hour(date)]
  result[, temporal_theme := fcase(
    hour_local < 5, "Night (0-4h)",
    hour_local < 12, "Morning (5-11h)",
    hour_local < 17, "Afternoon (12-16h)",
    default = "Evening (17-23h)"
  )]
  result[, temporal_theme := factor(temporal_theme,
                                     levels = c("Night (0-4h)", "Morning (5-11h)",
                                               "Afternoon (12-16h)", "Evening (17-23h)"))]
  result[, weekday_name := wday(date, label = TRUE, abbr = FALSE)]

  # --- Duration theme ---
  if ("dwell_min" %in% names(result)) {
    result[, duration_theme := fcase(
      dwell_min < 5, "Short (<5m)",
      dwell_min < 20, "Medium (5-20m)",
      default = "Long (>20m)"
    )]
    result[, duration_theme := factor(duration_theme,
                                       levels = c("Short (<5m)", "Medium (5-20m)", "Long (>20m)"))]
  }

  # --- Climate theme (requires temperature_c and humidity_pct columns) ---
  has_temp <- any(c("temperature_c", "temperature") %in% names(result))
  has_humid <- any(c("humidity_pct", "humidity") %in% names(result))

  if (has_temp && has_humid) {
    tcol <- intersect(c("temperature_c", "temperature"), names(result))[1]
    hcol <- intersect(c("humidity_pct", "humidity"), names(result))[1]
    tvals <- suppressWarnings(as.numeric(result[[tcol]]))
    hvals <- suppressWarnings(as.numeric(result[[hcol]]))
    tmed <- median(tvals, na.rm = TRUE)
    hmed <- median(hvals, na.rm = TRUE)
    if (!is.na(tmed) && !is.na(hmed)) {
      result[, temp_class := fifelse(tvals < tmed, "Cold", "Hot")]
      result[, humid_class := fifelse(hvals >= hmed, "Humid", "Dry")]
      result[, climate_theme := fifelse(
        is.na(temp_class) | is.na(humid_class), "Unknown",
        paste(temp_class, "and", humid_class)
      )]
    } else {
      result[, climate_theme := "Unknown"]
    }
  } else {
    result[, climate_theme := "Unknown"]
  }
  result[, climate_theme := factor(climate_theme,
                                    levels = c("Cold and Humid", "Cold and Dry",
                                              "Hot and Humid", "Hot and Dry", "Unknown"))]

  # --- Outdoor-specific themes ---
  outdoor_mask <- result$corrected_class == "Outdoor"

  if (sum(outdoor_mask) > 0 && "speed_m_min" %in% names(result)) {
    os <- result$speed_m_min[outdoor_mask]
    os_valid <- os[is.finite(os) & os > 0]
    if (length(os_valid) > 10) {
      oq33 <- quantile(os_valid, 0.33, na.rm = TRUE)
      oq67 <- quantile(os_valid, 0.67, na.rm = TRUE)
    } else {
      oq33 <- 30; oq67 <- 80
    }

    result[, outdoor_movement_type := NA_character_]
    result[outdoor_mask, outdoor_movement_type := fcase(
      speed_m_min <= oq33, "Slow Movement",
      speed_m_min <= oq67, "Moderate Movement",
      default = "Fast Movement"
    )]
    result[, outdoor_movement_type := factor(outdoor_movement_type,
                                              levels = c("Slow Movement", "Moderate Movement",
                                                        "Fast Movement"))]

    # Outdoor context (combines time-of-day with duration)
    if ("dwell_min" %in% names(result)) {
      result[, outdoor_context := NA_character_]
      result[outdoor_mask, outdoor_context := fcase(
        temporal_theme %in% c("Morning (5-11h)", "Afternoon (12-16h)") & dwell_min > 20,
          "Extended Daytime",
        temporal_theme %in% c("Morning (5-11h)", "Afternoon (12-16h)"),
          "Brief Daytime",
        temporal_theme == "Evening (17-23h)", "Evening",
        default = "Night"
      )]
      result[, outdoor_context := factor(outdoor_context,
                                          levels = c("Extended Daytime", "Brief Daytime",
                                                    "Evening", "Night"))]
    }
  }

  # --- Outdoor environmental stress (if temp + humidity available) ---
  if (has_temp && has_humid && sum(outdoor_mask) > 0) {
    tcol <- intersect(c("temperature_c", "temperature"), names(result))[1]
    hcol <- intersect(c("humidity_pct", "humidity"), names(result))[1]
    result[, heat_index := NA_real_]
    result[outdoor_mask, heat_index :=
      suppressWarnings(as.numeric(get(tcol))) +
      (0.5 * suppressWarnings(as.numeric(get(hcol))) / 10)]

    hi_vals <- result$heat_index[outdoor_mask]
    hi_valid <- hi_vals[is.finite(hi_vals)]
    if (length(hi_valid) > 10) {
      hi_q33 <- quantile(hi_valid, 0.33, na.rm = TRUE)
      hi_q67 <- quantile(hi_valid, 0.67, na.rm = TRUE)
    } else {
      hi_q33 <- 25; hi_q67 <- 35
    }

    result[, outdoor_stress_level := NA_character_]
    result[outdoor_mask, outdoor_stress_level := fcase(
      is.na(heat_index), "Unknown Stress",
      heat_index < hi_q33, "Low Stress",
      heat_index < hi_q67, "Moderate Stress",
      default = "High Stress"
    )]
    result[, outdoor_stress_level := factor(outdoor_stress_level,
                                             levels = c("Low Stress", "Moderate Stress",
                                                       "High Stress", "Unknown Stress"))]
  }

  return(result)
}

# ============================================================================
# Validation Pipeline Statistics
# ============================================================================

#' Compute validation pipeline statistics table
#' @param dt data.table with speed_class, osm_label, misclass_type, corrected_class
#' @param buildings_sf sf object with building polygons (for count)
#' @return data.frame with pipeline statistics
compute_validation_stats <- function(dt, buildings_sf = NULL) {
  n_total <- nrow(dt)
  has_speed <- "speed_class" %in% names(dt)
  has_osm <- "osm_label" %in% names(dt) && !all(is.na(dt$osm_label))
  has_misclass <- "misclass_type" %in% names(dt) && !all(is.na(dt$misclass_type))
  has_corrected <- "corrected_class" %in% names(dt)
  n_buildings <- if (!is.null(buildings_sf)) nrow(buildings_sf) else 0L
  n_speed_indoor <- if (has_speed) sum(dt$speed_class == "Indoor", na.rm = TRUE) else NA_integer_
  n_speed_outdoor <- if (has_speed) sum(dt$speed_class == "Outdoor", na.rm = TRUE) else NA_integer_
  n_osm_indoor <- if (has_osm) sum(dt$osm_label == "Indoor", na.rm = TRUE) else NA_integer_
  n_osm_outdoor <- if (has_osm) sum(dt$osm_label == "Outdoor", na.rm = TRUE) else NA_integer_
  n_correct <- if (has_misclass) sum(dt$misclass_type == "Correct", na.rm = TRUE) else NA_integer_
  n_false_indoor <- if (has_misclass) sum(dt$misclass_type == "False Indoor", na.rm = TRUE) else NA_integer_
  n_false_outdoor <- if (has_misclass) sum(dt$misclass_type == "False Outdoor", na.rm = TRUE) else NA_integer_
  n_final_indoor <- if (has_corrected) sum(dt$corrected_class == "Indoor", na.rm = TRUE) else NA_integer_
  n_final_outdoor <- if (has_corrected) sum(dt$corrected_class == "Outdoor", na.rm = TRUE) else NA_integer_
  osm_coverage_pct <- if (n_total > 0 && has_osm) round(100 * mean(!is.na(dt$osm_label)), 1) else NA_real_
  agreement_pct <- if (has_misclass) round(100 * mean(dt$misclass_type == "Correct", na.rm = TRUE), 1) else NA_real_
  validation_status <- if (!has_osm) {
    "Speed-only: no usable building labels"
  } else if (n_buildings == 0L) {
    "Speed-only: no buildings loaded"
  } else if (isTRUE(n_osm_indoor == 0L)) {
    "Building layer loaded, but no GPS points intersect buffered buildings"
  } else {
    "Building validation active"
  }

  data.frame(
    Step = c(
      "Total Points",
      "Step 1: Speed Indoor", "Step 1: Speed Indoor %",
      "Step 1: Speed Outdoor", "Step 1: Speed Outdoor %",
      "Step 2: OSM Buildings Downloaded",
      "Step 2: Building Validation Status",
      "Step 2: OSM Label Coverage %",
      "Step 2: Points In Buildings", "Step 2: In Buildings %",
      "Step 2: Points Outside Buildings",
      "Step 3: Correctly Classified", "Step 3: False Indoor",
      "Step 3: False Indoor %",
      "Step 3: False Outdoor", "Step 3: Misclassification Rate %",
      "Step 3: Speed/OSM Agreement %",
      "Step 4: Final Indoor", "Step 4: Final Indoor %",
      "Step 4: Final Outdoor", "Step 4: Final Outdoor %",
      "Net Change Indoor", "Net Change Outdoor"
    ),
    Value = c(
      n_total,
      n_speed_indoor,
      if (has_speed) round(100 * mean(dt$speed_class == "Indoor"), 1) else NA,
      n_speed_outdoor,
      if (has_speed) round(100 * mean(dt$speed_class == "Outdoor"), 1) else NA,
      n_buildings,
      validation_status,
      osm_coverage_pct,
      n_osm_indoor,
      if (has_osm) round(100 * mean(dt$osm_label == "Indoor", na.rm = TRUE), 1) else NA,
      n_osm_outdoor,
      n_correct,
      n_false_indoor,
      if (has_misclass) round(100 * mean(dt$misclass_type == "False Indoor", na.rm = TRUE), 1) else NA,
      n_false_outdoor,
      if (has_misclass) round(100 * mean(dt$misclass_type != "Correct", na.rm = TRUE), 1) else NA,
      agreement_pct,
      n_final_indoor,
      if (has_corrected) round(100 * mean(dt$corrected_class == "Indoor"), 1) else NA,
      n_final_outdoor,
      if (has_corrected) round(100 * mean(dt$corrected_class == "Outdoor"), 1) else NA,
      if (has_speed && has_corrected) n_final_indoor - n_speed_indoor else NA,
      if (has_speed && has_corrected) n_final_outdoor - n_speed_outdoor else NA
    ),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# Main Analysis Pipeline
# ============================================================================

#' Run the full trajectory analysis pipeline
#'
#' @param dt data.table with standardised trajectory data
#' @param boundary_sf sf object with boundary polygons (optional)
#' @param buildings_sf sf object with building polygons (optional)
#' @param fetch_osm Whether to fetch buildings from OSM
#' @param osm_buffer_km Buffer for OSM fetching
#' @param cell_size_m Grid cell size in metres
#' @param max_dwell_min Maximum dwell time per point
#' @param speed_threshold_kmh Teleport threshold
#' @param remove_teleports Remove teleporting points
#' @param building_buffer_m Buffer around buildings
#' @param exclude_indoor Exclude indoor-classified points
#' @return list with points, classified, person_day_cell, person_day, grid_params, stats
run_analysis <- function(dt,
                         boundary_sf = NULL,
                         buildings_sf = NULL,
                         fetch_osm = FALSE,
                         osm_buffer_km = 0.5,
                         cache_db_path = NULL,
                         cell_size_m = 1000,
                         max_dwell_min = 10,
                         speed_threshold_kmh = 120,
                         remove_teleports = TRUE,
                         building_buffer_m = 5,
                         exclude_indoor = FALSE,
                         speed_outdoor_threshold = 100,
                         use_building_validation = TRUE,
                         building_raster = NULL) {
  
  if (is.null(dt) || nrow(dt) == 0) {
    stop("No valid trajectory data loaded.")
  }
  
  working <- dt[!is.na(date) & !is.na(lat) & !is.na(lon) & !is.na(PID)]
  if (nrow(working) == 0) stop("No valid points after cleaning.")
  
  # Boundary filter
  if (!is.null(boundary_sf)) {
    working <- apply_boundary_filter(working, boundary_sf)
    if (nrow(working) == 0) stop("No points remain after boundary filtering.")
  }
  
  # Dwell time & speed
  working <- calculate_dwell_time(working, max_dwell_min, speed_threshold_kmh, remove_teleports)
  
  # Fetch OSM buildings (routed through cache when cache_db_path is set)
  if (is.null(buildings_sf) && fetch_osm) {
    buildings_sf <- fetch_osm_buildings(working,
                                         buffer_km = osm_buffer_km,
                                         cache_db  = cache_db_path)
  }

  # Ensure buildings always have a CRS before any spatial operations
  if (!is.null(buildings_sf) && inherits(buildings_sf, "sf")) {
    buildings_sf <- .ensure_crs(buildings_sf)
  }

  # Indoor classification — Step 2 (building check) is optional
  bldg_for_classification <- if (isTRUE(use_building_validation)) buildings_sf else NULL
  classified <- validated_indoor_outdoor(working, bldg_for_classification,
                                          speed_threshold = speed_outdoor_threshold,
                                          building_buffer_m = building_buffer_m,
                                          require_building_for_indoor = isTRUE(use_building_validation),
                                          building_raster = building_raster)

  # Add segmentation themes (mobility, temporal, outdoor-specific)
  classified <- compute_segmentation_themes(classified)

  # Compute validation pipeline statistics
  validation_stats <- compute_validation_stats(classified, buildings_sf)
  
  # Exclude indoor if requested
  if (exclude_indoor) {
    working <- classified[is_indoor == FALSE]
    if (nrow(working) == 0) {
      warning("All points were classified as indoor; proceeding with full dataset.")
      working <- copy(classified)
    }
  } else {
    working <- copy(classified)
  }

  if (nrow(working) == 0) {
    stop("No points remain after indoor filtering.")
  }
  
  # Grid aggregation
  grid_result <- aggregate_to_grid(working, cell_size_m)
  person_day_cell <- grid_result$person_day_cell
  grid_params     <- grid_result$grid_params
  full_grid       <- grid_result$full_grid
  daily_grid      <- grid_result$daily_grid
  
  # Daily summaries
  person_day <- create_daily_summaries(person_day_cell)

  # Per-participant average daily exposure with WHO comparison
  participant_exposure <- create_participant_summaries(person_day)

  # Statistics
  stats <- list(
    n_points_classified = nrow(classified),
    n_points_for_analysis = nrow(working),
    n_people = uniqueN(working$PID),
    n_days = uniqueN(as.Date(working$date)),
    total_dwell_min = sum(working$dwell_min, na.rm = TRUE),
    n_indoor_pred = sum(classified$is_indoor, na.rm = TRUE),
    n_outdoor_pred = sum(!classified$is_indoor, na.rm = TRUE)
  )
  
  # Release the large working copy — points_classified is the canonical output
  rm(working)
  gc(verbose = FALSE)

  list(
    points_classified   = classified,
    person_day_cell     = person_day_cell,
    person_day          = person_day,
    participant_exposure = participant_exposure,
    grid_params         = grid_params,
    stats               = stats,
    validation_stats    = validation_stats,
    buildings_sf        = buildings_sf,
    full_grid           = full_grid,
    daily_grid          = daily_grid
  )
}
