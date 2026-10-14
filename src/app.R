# ===========================================================================
# SPECTRA: Spatiotemporal Personal Exposure Characterization and TRAjectory
# Analyzer.
# Shiny web application: the interactive front end to the analysis engine.
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
# SPECTRA by CHEAQI-MNCH — Spatiotemporal Personal Exposure Characterization and TRAjectory Analyzer
# R Shiny Application
# ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(shinyjs)
  library(DT)
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(plotly)
  library(leaflet)
  library(ggplot2)
  library(viridis)
  library(sf)
  library(scales)
  library(forcats)
  library(RSQLite)
  library(digest)
  library(ncdf4)
  library(raster)
  library(emayili)
})

# Source the analysis core
# Try multiple possible locations for trajectory_analysis_core.R
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
core_file_options <- c(
  "trajectory_analysis_core.R",  # Same directory as app.R
  file.path(getwd(), "trajectory_analysis_core.R"),  # Current working directory
  file.path(getwd(), "shiny_app", "trajectory_analysis_core.R")  # If running from parent
)
if (!is.null(script_dir)) {
  core_file_options <- c(core_file_options,
    file.path(script_dir, "trajectory_analysis_core.R"))
}

core_file <- NULL
for (path in core_file_options) {
  if (!is.null(path) && file.exists(path)) {
    core_file <- path
    cat("Found trajectory_analysis_core.R at:", path, "\n")
    break
  }
}

if (is.null(core_file)) {
  stop("Cannot find trajectory_analysis_core.R. Please ensure you are running from the shiny_app directory.")
}

source(core_file, local = TRUE)

# ============================================================================
# Batched analysis wrapper (used on cloud to stay within 1 GB RAM)
# Splits by PID so each person's trajectory stays intact, classifies in
# batches of ~BATCH_PTS points, then aggregates the combined result.
# ============================================================================
run_analysis_batched <- function(dt, batch_pts = 100000L, progress_detail = NULL, ...) {
  pids      <- unique(dt$PID)
  pid_sizes <- dt[, .N, by = PID]
  setkey(pid_sizes, PID)

  batches      <- list()
  current_pids <- character(0)
  current_n    <- 0L

  for (pid in pids) {
    n <- pid_sizes[pid, N]
    if (current_n + n > batch_pts && length(current_pids) > 0) {
      batches <- c(batches, list(current_pids))
      current_pids <- character(0)
      current_n    <- 0L
    }
    current_pids <- c(current_pids, pid)
    current_n    <- current_n + n
  }
  if (length(current_pids) > 0) batches <- c(batches, list(current_pids))

  n_batches       <- length(batches)
  classified_list <- vector("list", n_batches)

  args <- list(...)
  buildings_sf             <- args$buildings_sf
  boundary_sf              <- args$boundary_sf
  fetch_osm                <- args$fetch_osm              %||% FALSE
  osm_buffer_km            <- args$osm_buffer_km          %||% 0.5
  cache_db_path            <- args$cache_db_path
  max_dwell_min            <- args$max_dwell_min          %||% 10
  speed_threshold_kmh      <- args$speed_threshold_kmh    %||% 120
  remove_teleports         <- args$remove_teleports       %||% TRUE
  building_buffer_m        <- args$building_buffer_m      %||% 5
  exclude_indoor           <- args$exclude_indoor         %||% FALSE
  speed_outdoor_threshold  <- args$speed_outdoor_threshold %||% 100
  use_building_validation  <- args$use_building_validation %||% TRUE
  cell_size_m              <- args$cell_size_m            %||% 1000

  bldg_for_class <- if (isTRUE(use_building_validation)) buildings_sf else NULL

  for (i in seq_len(n_batches)) {
    batch_dt <- dt[PID %in% batches[[i]]]

    if (!is.null(boundary_sf))
      batch_dt <- apply_boundary_filter(batch_dt, boundary_sf)
    if (nrow(batch_dt) == 0) next

    batch_dt  <- calculate_dwell_time(batch_dt, max_dwell_min,
                                      speed_threshold_kmh, remove_teleports)
    batch_cls <- validated_indoor_outdoor(batch_dt, bldg_for_class,
                                          speed_threshold    = speed_outdoor_threshold,
                                          building_buffer_m  = building_buffer_m,
                                          require_building_for_indoor = isTRUE(use_building_validation))
    batch_cls <- compute_segmentation_themes(batch_cls)
    classified_list[[i]] <- batch_cls

    rm(batch_dt, batch_cls)
    gc(verbose = FALSE)
  }

  classified <- rbindlist(Filter(Negate(is.null), classified_list), fill = TRUE)
  if (nrow(classified) == 0) stop("No points remain after batched classification.")

  validation_stats <- compute_validation_stats(classified, buildings_sf)

  working <- if (exclude_indoor) {
    w <- classified[is_indoor == FALSE]
    if (nrow(w) == 0) copy(classified) else w
  } else {
    copy(classified)
  }

  grid_result         <- aggregate_to_grid(working, cell_size_m)
  person_day_cell     <- grid_result$person_day_cell
  grid_params         <- grid_result$grid_params
  full_grid           <- grid_result$full_grid
  daily_grid          <- grid_result$daily_grid
  person_day          <- create_daily_summaries(person_day_cell)
  participant_exposure <- create_participant_summaries(person_day)

  stats <- list(
    n_points_classified   = nrow(classified),
    n_points_for_analysis = nrow(working),
    n_people              = uniqueN(working$PID),
    n_days                = uniqueN(as.Date(working$date)),
    total_dwell_min       = sum(working$dwell_min, na.rm = TRUE),
    n_indoor_pred         = sum(classified$is_indoor, na.rm = TRUE),
    n_outdoor_pred        = sum(!classified$is_indoor, na.rm = TRUE)
  )

  list(
    points_classified    = classified,
    person_day_cell      = person_day_cell,
    person_day           = person_day,
    participant_exposure = participant_exposure,
    grid_params          = grid_params,
    stats                = stats,
    validation_stats     = validation_stats,
    buildings_sf         = buildings_sf,
    full_grid            = full_grid,
    daily_grid           = daily_grid
  )
}

# ============================================================================
# Constants
# ============================================================================

TRUE_LABELS  <- c("1", "true", "t", "yes", "y", "indoor", "inside", "in")
FALSE_LABELS <- c("0", "false", "f", "no", "n", "outdoor", "outside", "out")

MAX_UPLOAD_MB <- 1024  # 1GB
BUILDING_CACHE_DB <- if (grepl("^/srv/connect", getwd())) {
  file.path(tempdir(), "buildings_cache.sqlite")
} else {
  "buildings_cache.sqlite"
}

# ============================================================================
# Server-level buildings tile cache (shared across all sessions)
# Tiles are loaded from SQLite once and held in RAM; max 60 tiles (~30min LRU)
# ============================================================================
.bldg_tile_cache     <- new.env(hash = TRUE, parent = emptyenv())
.bldg_tile_cache_ts  <- new.env(hash = TRUE, parent = emptyenv())  # load timestamps
BLDG_CACHE_MAX_TILES <- 60L   # evict oldest when over this
.preloaded_buildings_cache <- new.env(hash = TRUE, parent = emptyenv())

.load_preloaded_buildings_all <- function(path = "buildings_preloaded.csv") {
  if (exists("buildings", envir = .preloaded_buildings_cache, inherits = FALSE)) {
    return(.preloaded_buildings_cache$buildings)
  }
  if (!file.exists(path)) return(NULL)

  bdt <- tryCatch(data.table::fread(path, showProgress = FALSE),
                  error = function(e) NULL)
  if (is.null(bdt) || nrow(bdt) == 0 || !"geometry_wkt" %in% names(bdt)) return(NULL)

  bldgs <- tryCatch({
    geoms <- sf::st_as_sfc(bdt$geometry_wkt, crs = 4326L)
    x <- sf::st_sf(
      id = if ("id" %in% names(bdt)) bdt$id else seq_len(nrow(bdt)),
      geometry = geoms
    )
    x <- suppressWarnings(sf::st_make_valid(x))
    x <- x[!sf::st_is_empty(x), , drop = FALSE]
    sf::st_crs(x) <- 4326L
    x
  }, error = function(e) {
    message("Preloaded buildings parse error: ", e$message)
    NULL
  })
  if (!is.null(bldgs)) .preloaded_buildings_cache$buildings <- bldgs
  bldgs
}

.get_preloaded_buildings_for_bbox <- function(lat_min, lat_max, lon_min, lon_max,
                                              path = "buildings_preloaded.csv") {
  bldgs <- .load_preloaded_buildings_all(path)
  if (is.null(bldgs) || nrow(bldgs) == 0) return(NULL)

  bbox_poly <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = lon_min, ymin = lat_min, xmax = lon_max, ymax = lat_max),
    crs = 4326L
  ))
  sf::st_crs(bbox_poly) <- 4326L

  hits <- tryCatch(
    sf::st_intersects(bldgs, bbox_poly, sparse = FALSE)[, 1],
    error = function(e) {
      message("Preloaded buildings bbox intersection error: ", e$message)
      rep(FALSE, nrow(bldgs))
    }
  )
  out <- bldgs[hits, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  out
}

.cache_get_buildings <- function(lat_min, lat_max, lon_min, lon_max) {
  if (!file.exists(BUILDING_CACHE_DB)) {
    message("Buildings cache: ", BUILDING_CACHE_DB, " not found")
    return(NULL)
  }
  needed_tiles <- get_needed_tiles(lat_min, lat_max, lon_min, lon_max)
  if (length(needed_tiles) == 0) return(NULL)

  con <- .open_bldg_con(BUILDING_CACHE_DB, flags = RSQLite::SQLITE_RO)
  on.exit(RSQLite::dbDisconnect(con))

  all_sf <- list()

  for (tile in needed_tiles) {
    # Serve from cache if already loaded
    if (exists(tile, envir = .bldg_tile_cache, inherits = FALSE)) {
      .bldg_tile_cache_ts[[tile]] <- Sys.time()
      all_sf[[tile]] <- .bldg_tile_cache[[tile]]
      next
    }

    # Load from SQLite (tile index hit — fast). Prefer centroid columns when
    # present, but fall back to tile-only reads for legacy caches that only
    # stored geometry_wkt.
    rows <- tryCatch(
      RSQLite::dbGetQuery(con,
        "SELECT geometry_wkt FROM buildings
         WHERE tile_id = ?
           AND (
             latitude IS NULL OR longitude IS NULL OR
             (latitude BETWEEN ? AND ? AND longitude BETWEEN ? AND ?)
           )
         LIMIT 8000",
        params = list(tile, lat_min, lat_max, lon_min, lon_max)),
      error = function(e) {
        tryCatch(
          RSQLite::dbGetQuery(con,
            "SELECT geometry_wkt FROM buildings WHERE tile_id = ? LIMIT 8000",
            params = list(tile)),
          error = function(e2) NULL
        )
      }
    )
    if (is.null(rows) || nrow(rows) == 0) next

    geoms   <- tryCatch(sf::st_as_sfc(rows$geometry_wkt, crs = 4326L), error = function(e) NULL)
    if (is.null(geoms)) next
    tile_sf <- sf::st_sf(geometry = geoms)

    # Evict oldest tile if cache is full
    cache_keys <- ls(.bldg_tile_cache)
    if (length(cache_keys) >= BLDG_CACHE_MAX_TILES) {
      ts_vals  <- sapply(cache_keys, function(k) as.numeric(.bldg_tile_cache_ts[[k]]))
      oldest   <- cache_keys[which.min(ts_vals)]
      rm(list = oldest, envir = .bldg_tile_cache)
      rm(list = oldest, envir = .bldg_tile_cache_ts)
    }

    assign(tile, tile_sf, envir = .bldg_tile_cache)
    assign(tile, Sys.time(), envir = .bldg_tile_cache_ts)
    all_sf[[tile]] <- tile_sf
    message("Buildings cache: loaded tile ", tile, " (", nrow(tile_sf), " polygons)")
  }

  if (length(all_sf) == 0) return(NULL)

  # Merge tiles, enforce CRS, clip to bbox
  combined <- do.call(rbind, all_sf)
  sf::st_crs(combined) <- 4326L   # always set — rbind can drop CRS in some sf versions

  bbox_vec <- c(xmin = lon_min, ymin = lat_min, xmax = lon_max, ymax = lat_max)
  bbox_bbox <- sf::st_bbox(bbox_vec, crs = 4326L)
  bbox_poly <- sf::st_as_sfc(bbox_bbox)
  # Explicitly set CRS on the sfc object after creation
  sf::st_crs(bbox_poly) <- 4326L
  hits <- tryCatch(
    sf::st_intersects(combined, bbox_poly, sparse = FALSE)[, 1],
    error = function(e) rep(TRUE, nrow(combined))   # fallback: keep all
  )
  combined <- combined[hits, , drop = FALSE]
  if (nrow(combined) == 0) return(NULL)
  combined
}

.add_building_polygons <- function(map, buildings_sf, group = "Building Polygons") {
  if (is.null(buildings_sf) || !inherits(buildings_sf, "sf") || nrow(buildings_sf) == 0) {
    return(map)
  }

  bldgs <- buildings_sf
  if (is.na(sf::st_crs(bldgs))) sf::st_crs(bldgs) <- 4326L
  bldgs <- tryCatch(sf::st_transform(bldgs, 4326L), error = function(e) bldgs)
  bldgs <- suppressWarnings(sf::st_make_valid(bldgs))
  bldgs <- bldgs[!sf::st_is_empty(bldgs), , drop = FALSE]
  if (nrow(bldgs) == 0) return(map)

  leaflet::addPolygons(
    map,
    data = bldgs,
    group = group,
    color = "#0F6FAD",
    weight = 1,
    opacity = 0.75,
    fillColor = "#4A90E2",
    fillOpacity = 0.18,
    smoothFactor = 0.5,
    popup = "Building footprint"
  )
}

# ============================================================================
# Custom CSS
# ============================================================================

app_css <- "
/* ============================================================
   SPECTRA — Dark Navy, Sky Blue & Orange Theme
   CHEAQI-MNCH  ·  2025
   ============================================================
   Palette:
     --navy-deep: #061120   body background (deepest navy)
     --navy:      #0A1A30   sidebar, surfaces
     --navy-mid:  #0E2240   box headers
     --sky:       #0077B6   sky blue accent (borders, highlights)
     --sky-bright:#00B4D8   bright sky blue (icons)
     --orange:    #F97316   orange (nav active, primary buttons)
     --orange-dk: #C2410C   dark orange (hover)
     --white:     #FFFFFF   white content cards
     --text:      #E8F0FE   light text on dark backgrounds
   ============================================================ */

/* ── Base ─────────────────────────────────────────────────── */
body {
  background-color: #061120;
  background-image: url('honeycomb_bg.svg');
  background-size: 156px 135px;
  background-repeat: repeat;
  font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
  font-size: 13px;
  color: #E8F0FE;
}

/* ── Header — hidden entirely ─────────────────────────────── */
.main-header { display: none !important; }

/* Prevent header text overflow (kept for safety) */
.main-header .logo > span {
  display: flex !important;
  align-items: center !important;
  overflow: hidden !important;
  white-space: nowrap !important;
  max-width: 100% !important;
}
.main-header .logo > span span {
  overflow: hidden !important;
  text-overflow: ellipsis !important;
  white-space: nowrap !important;
}

/* ── Sidebar ─────────────────────────────────────────────── */
.skin-blue .main-sidebar {
  background: linear-gradient(180deg, #061120 0%, #0A1A30 100%) !important;
  width: 220px !important;
}
.skin-blue .main-sidebar .sidebar-menu > li > a {
  font-size: 12.5px !important;
  padding: 9px 12px 9px 15px !important;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  color: rgba(255,255,255,0.82) !important;
  border-left: 3px solid transparent;
  transition: all 0.2s;
}
.skin-blue .main-sidebar .sidebar-menu > li.active > a,
.skin-blue .main-sidebar .sidebar-menu > li > a:hover {
  background: rgba(255,255,255,0.1) !important;
  border-left-color: #F97316 !important;
  color: #fff !important;
}
.sidebar-menu > li > a > .fa,
.sidebar-menu > li > a > .glyphicon,
.sidebar-menu > li > a > .ion {
  font-size: 13px !important;
  width: 18px !important;
  min-width: 18px !important;
  margin-right: 6px !important;
}
.main-sidebar .sidebar-menu .treeview-menu > li > a {
  font-size: 12px !important;
  padding-left: 30px !important;
}
/* Sidebar badge chips */
.sidebar-menu > li > a .badge {
  font-size: 9px !important;
  padding: 2px 5px !important;
  float: right;
  margin-top: 2px;
}
/* Sidebar footer */
.skin-blue .main-sidebar .sidebar-footer,
.main-sidebar > div[style*='padding: 10px'] {
  font-size: 11px !important;
  color: rgba(255,255,255,0.55) !important;
  line-height: 1.5 !important;
}

/* ── Content area — honeycomb background ─────────────────── */
.content-wrapper {
  background: transparent !important;
  margin-left: 0 !important;
}
.content-wrapper { position: relative; }
.content-wrapper > * { position: relative; z-index: 1; }
.content {
  padding: 10px 14px !important;
}

/* ── Section banners ─────────────────────────────────────── */
.section-banner {
  background: linear-gradient(135deg, #061120 0%, #0A1A30 50%, #0E2240 100%);
  color: white;
  padding: 10px 18px !important;
  border-radius: 10px;
  margin-bottom: 14px;
  box-shadow: 0 4px 20px rgba(0,119,182,0.25);
  border: 1px solid rgba(0,119,182,0.20);
  display: flex;
  align-items: center;
  gap: 12px;
  overflow: hidden;
}
.section-banner h2 {
  font-size: 1.05em !important;
  font-weight: 700 !important;
  margin: 0 !important;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  letter-spacing: 0.5px;
}
.section-banner h2 .fa {
  font-size: 0.95em !important;
  margin-right: 6px;
}
.section-banner p {
  font-size: 0.78em !important;
  margin: 0 !important;
  opacity: 0.82;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}
/* Orange accent rule on the right */
.section-banner::after {
  content: '';
  position: absolute;
  right: 0; top: 0; bottom: 0;
  width: 4px;
  background: linear-gradient(180deg, #F97316, #00B4D8);
  border-radius: 0 10px 10px 0;
}
.section-banner { position: relative; }

/* ── Box cards — white content cards ─────────────────────── */
.box {
  border-radius: 10px !important;
  box-shadow: 0 4px 20px rgba(0,0,0,0.25) !important;
  border-top: 3px solid #0077B6 !important;
  background: #ffffff !important;
  margin-bottom: 14px !important;
  color: #0A1A30 !important;
}
.box-header {
  padding: 8px 12px !important;
  border-radius: 10px 10px 0 0 !important;
  background: #0E2240 !important;
  border-bottom: 1px solid rgba(0,119,182,0.25) !important;
}
.box-header .box-title {
  font-size: 13px !important;
  font-weight: 600 !important;
  color: #E8F0FE !important;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  display: block;
  max-width: 100%;
}
.box-header .box-title .fa {
  font-size: 12px !important;
  margin-right: 5px;
  color: #00B4D8 !important;
}
.box-body {
  padding: 10px 12px !important;
  font-size: 12.5px !important;
  color: #0A1A30 !important;
  background: #ffffff !important;
}
.box-body p, .box-body label, .box-body small {
  font-size: 12px !important;
  line-height: 1.4 !important;
  color: #1a2e40 !important;
}
.box-body h3, .box-body h4, .box-body h5 {
  color: #0A1A30 !important;
}
/* Wrap long labels in inputs */
.form-group label {
  white-space: normal !important;
  word-break: break-word !important;
  font-size: 12px !important;
  font-weight: 600 !important;
  margin-bottom: 3px !important;
  color: #0A1A30 !important;
}
/* Checkboxes and radio buttons inside box */
.box-body .checkbox label, .box-body .radio label {
  color: #1a2e40 !important;
}
/* Select inputs and text inputs */
.box-body .form-control {
  color: #0A1A30 !important;
  background: #f0f4f8 !important;
  border: 1px solid #b0c4d8 !important;
}
/* hrule inside box */
.box-body hr { border-color: #dde3ea !important; }
.form-control {
  font-size: 12px !important;
  height: 30px !important;
  padding: 4px 8px !important;
  border-radius: 5px !important;
}
textarea.form-control { height: auto !important; }
.shiny-input-container { margin-bottom: 8px !important; }

/* ── Info boxes (dashboard KPIs) ─────────────────────────── */
.info-box {
  min-height: 58px !important;
  border-radius: 8px !important;
  box-shadow: 0 2px 8px rgba(0,0,0,0.08) !important;
  margin-bottom: 10px !important;
}
.info-box .info-box-icon {
  width: 56px !important;
  height: 58px !important;
  line-height: 58px !important;
  font-size: 24px !important;
  border-radius: 8px 0 0 8px !important;
}
.info-box-content {
  padding: 4px 8px !important;
  margin-left: 56px !important;
}
.info-box-text {
  font-size: 10.5px !important;
  white-space: normal !important;
  line-height: 1.2 !important;
}
.info-box-number {
  font-size: 15px !important;
  font-weight: 700 !important;
}

/* ── Tabs ────────────────────────────────────────────────── */
.nav-tabs-custom > .nav-tabs > li > a {
  font-size: 12px !important;
  padding: 6px 12px !important;
}
.nav-tabs-custom > .nav-tabs > li.active > a {
  border-top-color: #0077B6 !important;
  color: #061120 !important;
  font-weight: 700 !important;
}

/* ── Buttons ─────────────────────────────────────────────── */
.btn {
  font-size: 12px !important;
  padding: 5px 12px !important;
  border-radius: 6px !important;
}
.btn-primary {
  background: linear-gradient(135deg, #0077B6, #0E2240) !important;
  border: 1px solid rgba(0,180,216,0.35) !important;
  font-weight: 600 !important;
  color: #fff !important;
}
.btn-primary:hover {
  background: linear-gradient(135deg, #0096D6, #1A3560) !important;
  border-color: #00B4D8 !important;
}
.btn-success {
  background: linear-gradient(135deg, #16863F, #22C55E) !important;
  border: none !important;
  color: #fff !important;
  font-weight: 600 !important;
}
.btn-success:hover {
  background: linear-gradient(135deg, #22C55E, #16863F) !important;
}
.btn-danger  { border-radius: 6px !important; }
.btn-warning { border-radius: 6px !important; }
.btn-success { border-radius: 6px !important; }
.btn-sm { font-size: 11px !important; padding: 4px 9px !important; }

/* ── Alert boxes ─────────────────────────────────────────── */
.alert {
  font-size: 11.5px !important;
  padding: 7px 10px !important;
  border-radius: 6px !important;
  margin-bottom: 8px !important;
  line-height: 1.4 !important;
}

/* ── Metric cards ────────────────────────────────────────── */
.metric-card {
  background: white;
  border-radius: 8px;
  padding: 10px 12px;
  text-align: center;
  box-shadow: 0 2px 6px rgba(0,0,0,0.10);
  border-left: 3px solid #0077B6;
  margin-bottom: 10px;
}
.metric-value {
  font-size: 1.7em;
  font-weight: 700;
  color: #061120;
  line-height: 1.1;
}
.metric-label {
  font-size: 0.78em;
  color: #555;
  margin-top: 3px;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

/* ── DT tables ───────────────────────────────────────────── */
.dataTables_wrapper .dataTables_filter input {
  color: #333 !important;
  background: white !important;
  border: 1px solid #ccc !important;
  border-radius: 4px !important;
  padding: 3px 7px !important;
  font-size: 12px !important;
}
table.dataTable td, table.dataTable th {
  font-size: 12px !important;
  padding: 4px 8px !important;
}

/* ── radioGroupButtons compact ───────────────────────────── */
.btn-group-justified .btn {
  font-size: 11.5px !important;
  padding: 5px 4px !important;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

/* ── Horizontal rules ────────────────────────────────────── */
hr { margin: 8px 0 !important; border-color: #dee3ea !important; }

/* ── Shiny progress + notifications — always on top ─────── */
.shiny-progress-container { z-index: 9999 !important; }
.shiny-progress .progress  { background: rgba(6,17,32,0.92) !important; border: 1px solid rgba(0,119,182,0.35) !important; border-radius: 10px !important; }
.shiny-progress .progress-text { color: #E8F0FE !important; }
.shiny-progress .progress-bar  { background: linear-gradient(90deg, #0077B6, #F97316) !important; }
.shiny-notification-panel  { z-index: 9999 !important; }
.shiny-notification         { background: #0E2240 !important; border-left: 4px solid #0077B6 !important; color: #E8F0FE !important; border-radius: 8px !important; }
.shiny-notification-error   { border-left-color: #E84040 !important; }
.shiny-notification-warning { border-left-color: #F97316 !important; }
.shiny-notification-message { border-left-color: #22C55E !important; }

/* ── Login page ──────────────────────────────────────────── */
.login-page {
  width: 100vw;
  height: 100vh;
  background-color: #061120;
  background-image: url('honeycomb_bg.svg');
  background-size: 156px 135px;
  background-repeat: repeat;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  position: relative;
}
.login-page::before {
  content: '';
  position: absolute;
  inset: 0;
  background: rgba(6, 17, 32, 0.55);
  z-index: 0;
}
.login-branding, .login-box, .login-footer {
  position: relative;
  z-index: 1;
}
.login-branding {
  text-align: center;
  margin-bottom: 20px;
  color: white;
}
.login-branding h1 {
  font-family: 'Segoe UI', Arial, sans-serif !important;
  font-size: 7em;
  font-weight: 900;
  color: #F97316;
  margin-bottom: 4px;
  letter-spacing: 14px;
  text-shadow: 0 4px 24px rgba(0,0,0,0.65);
  line-height: 1;
}
.login-branding h3 {
  font-family: 'Segoe UI', Arial, sans-serif !important;
  color: rgba(255,255,255,0.90);
  font-size: 1.15em;
  font-weight: 300;
  line-height: 1.6;
  max-width: 580px;
  margin: 10px auto 6px auto;
  letter-spacing: 0.6px;
}
.login-branding h3 b { color: #F97316; font-weight: 700; }
.login-branding p { color: rgba(255,255,255,0.6); font-size: 1em; margin: 0; letter-spacing: 1.5px; }
.login-box {
  background: rgba(10,26,48,0.78);
  backdrop-filter: blur(20px);
  -webkit-backdrop-filter: blur(20px);
  border-radius: 18px;
  padding: 36px 34px;
  width: 440px;
  max-width: 94vw;
  box-shadow: 0 10px 48px rgba(0,0,0,0.60);
  border: 1px solid rgba(0,119,182,0.32);
  color: white;
}
.login-box h2 {
  text-align: center;
  font-size: 1.55em;
  font-weight: 700;
  margin-bottom: 22px;
  letter-spacing: 0.5px;
}
.login-box .form-control {
  background: rgba(255,255,255,0.18) !important;
  border: 1px solid rgba(255,255,255,0.30) !important;
  color: white !important;
  border-radius: 9px !important;
  height: 44px !important;
  font-size: 15px !important;
}
.login-box .form-control::placeholder { color: rgba(255,255,255,0.55); font-size: 14px; }
.login-box .btn-login, .login-box .btn-signup {
  width: 100%;
  padding: 13px;
  font-size: 16px;
  font-weight: 700;
  border-radius: 9px;
  margin-top: 12px;
  cursor: pointer;
  transition: all 0.25s;
  letter-spacing: 0.5px;
}
.login-box .btn-login {
  background: linear-gradient(135deg, #F97316, #C2410C);
  border: none;
  color: #ffffff;
}
.login-box .btn-login:hover { background: linear-gradient(135deg, #FB923C, #EA580C); }
.login-box .btn-signup {
  background: transparent;
  border: 1.5px solid rgba(255,255,255,0.45);
  color: white;
}
.login-box .btn-signup:hover { background: rgba(255,255,255,0.12); }
.login-toggle { text-align: center; margin-top: 14px; color: rgba(255,255,255,0.70); font-size: 14px; }
.login-toggle a { color: #F97316; cursor: pointer; text-decoration: none; font-weight: 700; }
.login-footer {
  margin-top: 20px;
  text-align: center;
  color: rgba(255,255,255,0.50);
  font-size: 0.95em;
}
.login-footer a { color: rgba(255,255,255,0.70); text-decoration: none; }
.login-footer a:hover { color: #F97316; }

/* ── Horizontal navigation bar ─────────────────────────────────────────── */
.spectra-nav-wrap {
  position: fixed;
  top: 50px;
  left: 0;
  width: 100%;
  z-index: 900;
  background: linear-gradient(135deg, #061120 0%, #0A1A30 60%, #0E2240 100%);
  box-shadow: 0 3px 12px rgba(0,0,0,0.45);
  border-bottom: 3px solid #F97316;
  padding: 0;
  display: flex;
  align-items: stretch;
}
.spectra-nav-brand {
  display: flex;
  flex-direction: column;
  justify-content: center;
  padding: 8px 18px;
  min-width: 140px;
  border-right: 1px solid rgba(255,255,255,0.12);
  background: rgba(0,0,0,0.2);
}
.spectra-nav-logo {
  font-family: 'Segoe UI', Arial, sans-serif;
  font-size: 1.4em;
  font-weight: 900;
  color: #F97316;
  letter-spacing: 3px;
  line-height: 1;
}
.spectra-nav-byline {
  font-size: 0.62em;
  color: rgba(255,255,255,0.45);
  font-style: italic;
  margin-top: 2px;
  line-height: 1.2;
}
.spectra-nav-row {
  display: flex;
  flex: 1;
  gap: 0;
}

/* ── Tab buttons — large card style ─────────────────────────── */
.snav-btn {
  display: inline-flex !important;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 4px;
  padding: 10px 14px !important;
  flex: 1;
  min-height: 70px;
  font-size: 10.5px !important;
  font-weight: 700 !important;
  color: rgba(255,255,255,0.75) !important;
  background: transparent !important;
  border: none !important;
  border-right: 1px solid rgba(255,255,255,0.08) !important;
  text-decoration: none !important;
  cursor: pointer;
  transition: all 0.18s ease;
  text-align: center;
  letter-spacing: 0.4px;
  text-transform: uppercase;
  line-height: 1.2;
  border-radius: 0 !important;
  position: relative;
}
.snav-btn::after {
  content: '';
  position: absolute;
  bottom: 0; left: 0; right: 0;
  height: 3px;
  background: #F97316;
  opacity: 0;
  transition: opacity 0.18s;
}
.snav-btn:hover {
  background: rgba(255,255,255,0.09) !important;
  color: #ffffff !important;
  text-decoration: none !important;
}
.snav-btn:hover::after { opacity: 0.6; }
.snav-btn.snav-active {
  background: rgba(249,115,22,0.15) !important;
  color: #F97316 !important;
}
.snav-btn.snav-active::after { opacity: 1; }
.snav-btn .fa {
  font-size: 20px !important;
  display: block;
  margin-bottom: 3px;
  opacity: 0.9;
}
.snav-btn:hover .fa, .snav-btn.snav-active .fa { opacity: 1; }

/* About tab — visually distinct, pushed to far right */
.snav-about {
  min-width: 70px !important;
  flex: 0 0 auto !important;
  background: rgba(0,0,0,0.2) !important;
  border-left: 1px solid rgba(255,255,255,0.12) !important;
  border-right: none !important;
}

/* No sidebar, no header */
.content-wrapper {
  margin-left: 0 !important;
  padding-top: 0 !important;
}
.content {
  padding-top: 0 !important;
}
.sidebar-toggle { display: none !important; }
.main-footer { margin-left: 0 !important; }

/* Home / Menu button in each section banner */
.section-banner { position: relative; }
.home-btn {
  position: absolute !important;
  top: 50% !important;
  right: 18px !important;
  transform: translateY(-50%) !important;
  background: rgba(255,255,255,0.15) !important;
  border: 1px solid rgba(255,255,255,0.4) !important;
  color: #fff !important;
  border-radius: 6px !important;
  padding: 5px 14px !important;
  font-size: 12px !important;
  font-weight: 600 !important;
  cursor: pointer !important;
  transition: background 0.15s !important;
  box-shadow: none !important;
  line-height: 1.4 !important;
}
.home-btn:hover,
.home-btn:focus,
.home-btn:active {
  background: rgba(255,255,255,0.3) !important;
  color: #fff !important;
  outline: none !important;
  box-shadow: none !important;
}

/* ── Landing / Home page ─────────────────────────────────── */
.spectra-home-wrap {
  position: relative;
  z-index: 2;
  min-height: 100vh;
  margin: -10px -14px;
  padding: 60px 40px 80px;
  background-color: #061120;
  background-image: url('honeycomb_bg.svg');
  background-size: 156px 135px;
  background-repeat: repeat;
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  overflow: hidden;
}
.spectra-home-wrap::before {
  content: '';
  position: absolute;
  inset: 0;
  background: rgba(6, 17, 32, 0.55);
  z-index: 0;
}
/* Watermark logo — top right, 20% opacity */
.spectra-home-watermark {
  position: absolute;
  top: -10px;
  right: -10px;
  width: 320px;
  opacity: 0.20;
  pointer-events: none;
  z-index: 1;
}
/* All content above the overlay */
.spectra-home-wrap > *:not(.spectra-home-wrap::before) {
  position: relative;
  z-index: 2;
}
/* SPECTRA title — orange style matching login */
.spectra-home-logo {
  font-family: 'Segoe UI', Arial, sans-serif;
  font-size: 4.8em;
  font-weight: 900;
  color: #F97316;
  letter-spacing: 10px;
  text-shadow: 0 3px 18px rgba(0,0,0,0.55);
  margin-bottom: 4px;
  line-height: 1;
  position: relative;
  z-index: 2;
}
/* Tagline — same white light style as login */
.spectra-home-tagline {
  font-family: 'Segoe UI', Arial, sans-serif;
  font-size: 0.92em;
  font-weight: 300;
  color: rgba(255,255,255,0.85);
  line-height: 1.6;
  max-width: 520px;
  margin: 6px auto 0;
  letter-spacing: 0.5px;
  position: relative;
  z-index: 2;
}
.spectra-home-tagline b { color: #F97316; font-weight: 700; }

/* 4-card row */
.spectra-home-cards {
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: 22px;
  width: 100%;
  max-width: 980px;
  margin: 70px auto 0;    /* pushed well down from title */
  position: relative;
  z-index: 2;
}
/* 3-card row */
.spectra-home-cards-row2 {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 22px;
  width: 100%;
  max-width: 735px;
  margin: 22px auto 0;
  position: relative;
  z-index: 2;
}
/* Home module card buttons */
.home-card-btn {
  background: rgba(14,34,64,0.82) !important;
  border: 1px solid rgba(255,255,255,0.18) !important;
  border-bottom: 4px solid rgba(249,115,22,0.55) !important;
  border-radius: 14px !important;
  color: #ffffff !important;
  width: 100% !important;
  min-height: 190px !important;
  padding: 36px 16px 28px !important;
  display: inline-flex !important;
  flex-direction: column !important;
  align-items: center !important;
  justify-content: center !important;
  cursor: pointer !important;
  transition: all 0.22s ease !important;
  box-shadow: 0 6px 24px rgba(0,0,0,0.28) !important;
  white-space: normal !important;
  text-align: center !important;
  line-height: 1.4 !important;
}
.home-card-btn:hover,
.home-card-btn:focus,
.home-card-btn:active {
  background: rgba(249,115,22,0.18) !important;
  border-color: rgba(249,115,22,0.70) !important;
  border-bottom-color: #F97316 !important;
  box-shadow: 0 12px 36px rgba(0,0,0,0.40) !important;
  transform: translateY(-6px) !important;
  color: #F97316 !important;
  outline: none !important;
}
.home-card-btn .fa {
  font-size: 2.8em !important;
  display: block !important;
  margin-bottom: 16px !important;
  opacity: 0.9;
}
.home-card-btn:hover .fa { opacity: 1; color: #F97316 !important; }
.home-card-label {
  font-size: 12px;
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: 1px;
  display: block;
  line-height: 1.45;
}
.home-card-btn-muted {
  background: rgba(255,255,255,0.06) !important;
  border-color: rgba(255,255,255,0.12) !important;
}

/* ── Live map panel (sidebar style in config) ────────────── */
.live-map-sidebar {
  position: sticky;
  top: 140px;
}

/* ── Config section boxes ────────────────────────────────── */
.config-section-num {
  display: inline-block;
  background: #0E2240;
  color: #00B4D8;
  border-radius: 50%;
  width: 22px; height: 22px;
  text-align: center;
  line-height: 22px;
  font-size: 11px;
  font-weight: 700;
  margin-right: 6px;
}
"

# ============================================================================
# SMTP Email Configuration
# ============================================================================
# Configure your SMTP server details here to enable email verification.
# Supports Gmail, Outlook, Yahoo, or any SMTP server.
#
# For Gmail: host="smtp.gmail.com", port=465, use_ssl=TRUE
#   - Use an App Password (not your regular password):
#     https://myaccount.google.com/apppasswords
#
# For Outlook/Hotmail: host="smtp-mail.outlook.com", port=587, use_ssl=FALSE (STARTTLS)
# For Yahoo: host="smtp.mail.yahoo.com", port=465, use_ssl=TRUE

SMTP_CONFIG <- list(
  host     = Sys.getenv("SPECTRA_SMTP_HOST", "smtp.gmail.com"),
  port     = as.integer(Sys.getenv("SPECTRA_SMTP_PORT", "465")),
  username = Sys.getenv("SPECTRA_SMTP_USER", ""),     # e.g. yourname@gmail.com
  password = Sys.getenv("SPECTRA_SMTP_PASS", ""),     # App password
  from     = Sys.getenv("SPECTRA_SMTP_FROM", ""),     # sender email (if blank, uses username)
  use_ssl  = as.logical(Sys.getenv("SPECTRA_SMTP_SSL", "TRUE"))
)

# Load SMTP config from local file if it exists (overrides env vars)
if (file.exists("smtp_config.R")) {
  source("smtp_config.R", local = TRUE)
}

smtp_configured <- function() {
  SMTP_CONFIG$username != "" && SMTP_CONFIG$password != ""
}

send_verification_email <- function(to_email, to_name, code) {
  if (!smtp_configured()) {
    message("[SPECTRA] SMTP not configured. Verification code for ", to_email, ": ", code)
    return(list(success = FALSE, msg = "SMTP not configured. Code shown on screen."))
  }

  from_addr <- if (SMTP_CONFIG$from != "") SMTP_CONFIG$from else SMTP_CONFIG$username

  email_body <- paste0(
    "<html><body style='font-family: Arial, sans-serif; background: #f0f4f8; padding: 20px;'>",
    "<div style='max-width: 500px; margin: auto; background: white; border-radius: 15px; padding: 30px; box-shadow: 0 4px 15px rgba(0,0,0,0.1);'>",
    "<div style='text-align: center; margin-bottom: 20px;'>",
    "<h1 style='color: #1e3c72; margin: 0;'>SPECTRA</h1>",
    "<p style='color: #666; margin: 5px 0 0 0;'>Spatiotemporal Personal Exposure Characterization and TRAjectory Analyzer</p>",
    "</div>",
    "<hr style='border: none; border-top: 2px solid #4a90e2; margin: 20px 0;'>",
    "<p>Hello <strong>", htmltools::htmlEscape(to_name), "</strong>,</p>",
    "<p>Thank you for registering with SPECTRA. Please use the verification code below to activate your account:</p>",
    "<div style='text-align: center; margin: 25px 0;'>",
    "<div style='display: inline-block; background: linear-gradient(135deg, #1e3c72, #2a5298); color: white; font-size: 2em; font-weight: bold; letter-spacing: 8px; padding: 15px 30px; border-radius: 10px;'>",
    htmltools::htmlEscape(code),
    "</div>",
    "</div>",
    "<p style='color: #666; font-size: 0.9em;'>If you did not create an account, you can safely ignore this email.</p>",
    "<hr style='border: none; border-top: 1px solid #eee; margin: 20px 0;'>",
    "<p style='color: #999; font-size: 0.8em; text-align: center;'>SPECTRA by CHEAQI-MNCH</p>",
    "</div></body></html>"
  )

  tryCatch({
    smtp <- if (isTRUE(SMTP_CONFIG$use_ssl)) {
      server(host = SMTP_CONFIG$host, port = SMTP_CONFIG$port,
             username = SMTP_CONFIG$username, password = SMTP_CONFIG$password)
    } else {
      server(host = SMTP_CONFIG$host, port = SMTP_CONFIG$port,
             username = SMTP_CONFIG$username, password = SMTP_CONFIG$password)
    }

    msg <- envelope() |>
      from(from_addr) |>
      to(to_email) |>
      subject("SPECTRA - Email Verification Code") |>
      html(email_body)

    smtp(msg)
    message("[SPECTRA] Verification email sent to: ", to_email)
    return(list(success = TRUE, msg = "Verification code sent to your email."))
  }, error = function(e) {
    message("[SPECTRA] Email send failed: ", e$message)
    return(list(success = FALSE, msg = paste("Could not send email:", e$message)))
  })
}

# ============================================================================
# SQLite User Database Setup
# ============================================================================

init_user_db <- function(db_path = "spectra_users.sqlite") {
  con <- dbConnect(SQLite(), db_path)
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS users (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      email TEXT UNIQUE NOT NULL,
      password_hash TEXT NOT NULL,
      full_name TEXT,
      organization TEXT,
      role TEXT DEFAULT 'user',
      verified INTEGER DEFAULT 0,
      verification_code TEXT,
      created_at TEXT DEFAULT (datetime('now')),
      last_login TEXT
    )
  ")
  # Add verified column if upgrading from old schema
  tryCatch({
    dbExecute(con, "ALTER TABLE users ADD COLUMN verified INTEGER DEFAULT 0")
  }, error = function(e) NULL)
  tryCatch({
    dbExecute(con, "ALTER TABLE users ADD COLUMN verification_code TEXT")
  }, error = function(e) NULL)

  # Create default admin if not exists
  admin_exists <- dbGetQuery(con, "SELECT COUNT(*) as n FROM users WHERE role = 'admin'")$n
  if (admin_exists == 0) {
    # The seed password is read from the environment. No default is shipped: if
    # SPECTRA_ADMIN_PASSWORD is unset, no administrator account is created.
    admin_pw <- Sys.getenv("SPECTRA_ADMIN_PASSWORD", "")
    if (nzchar(admin_pw)) {
      admin_hash <- digest(admin_pw, algo = "sha256", serialize = FALSE)
      dbExecute(con, "INSERT OR IGNORE INTO users (email, password_hash, full_name, role, verified) VALUES (?, ?, ?, ?, 1)",
                params = list(Sys.getenv("SPECTRA_ADMIN_EMAIL", "admin@spectra.cheaqi.org"),
                              admin_hash, "SPECTRA Admin", "admin"))
    } else {
      message("No administrator account created. Set SPECTRA_ADMIN_PASSWORD to seed one.")
    }
  }
  dbDisconnect(con)
}

verify_user <- function(email, password, db_path = "spectra_users.sqlite") {
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))
  pw_hash <- digest(password, algo = "sha256", serialize = FALSE)
  user <- dbGetQuery(con, "SELECT * FROM users WHERE email = ? AND password_hash = ?",
                     params = list(email, pw_hash))
  if (nrow(user) == 1) {
    if (!isTRUE(user$verified == 1)) {
      return(list(error = "email_not_verified", code = user$verification_code))
    }
    dbExecute(con, "UPDATE users SET last_login = datetime('now') WHERE id = ?",
              params = list(user$id))
    return(user)
  }
  return(NULL)
}

generate_verification_code <- function() {
  paste0(sample(0:9, 6, replace = TRUE), collapse = "")
}

register_user <- function(email, password, full_name, organization = "", db_path = "spectra_users.sqlite") {
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))
  existing <- dbGetQuery(con, "SELECT * FROM users WHERE email = ?", params = list(email))
  if (nrow(existing) > 0) {
    # If unverified, regenerate code and let them verify
    if (!isTRUE(existing$verified[1] == 1)) {
      vcode <- generate_verification_code()
      dbExecute(con, "UPDATE users SET verification_code = ? WHERE email = ?",
                params = list(vcode, email))
      return(list(success = TRUE, msg = "Verification code regenerated.", code = vcode))
    }
    return(list(success = FALSE, msg = "Email already registered"))
  }
  pw_hash <- digest(password, algo = "sha256", serialize = FALSE)
  vcode <- generate_verification_code()
  tryCatch({
    dbExecute(con, "INSERT INTO users (email, password_hash, full_name, organization, verified, verification_code) VALUES (?, ?, ?, ?, 0, ?)",
              params = list(email, pw_hash, full_name, organization, vcode))
    return(list(success = TRUE, msg = "Account created!", code = vcode))
  }, error = function(e) {
    return(list(success = FALSE, msg = paste("Registration failed:", e$message)))
  })
}

verify_email_code <- function(email, code, db_path = "spectra_users.sqlite") {
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))
  user <- dbGetQuery(con, "SELECT * FROM users WHERE email = ? AND verification_code = ?",
                     params = list(email, code))
  if (nrow(user) == 1) {
    dbExecute(con, "UPDATE users SET verified = 1, verification_code = NULL WHERE email = ?",
              params = list(email))
    return(TRUE)
  }
  return(FALSE)
}

# Initialize the database on app start
init_user_db()

# ============================================================================
# UI
# ============================================================================

login_page <- div(class = "login-page",
  div(class = "login-branding",
      img(src = "cheaqi_logo.png", height = "60px", style = "margin-bottom: 12px; opacity:0.9;"),
      h1("SPECTRA"),
      h3(HTML("<b>S</b>patiotemporal <b>P</b>ersonal <b>E</b>xposure <b>C</b>haracterization And <b>TRA</b>jectory Analyzer")),
      p("By CHEAQI-MNCH")
  ),

  div(class = "login-box", id = "login_form_box",
      h2("Welcome Back"),
      div(id = "login_form",
          textInput("login_email", NULL, placeholder = "Email Address"),
          passwordInput("login_password", NULL, placeholder = "Password"),
          div(style = "margin-top: 5px;",
              uiOutput("login_msg")),
          actionButton("btn_login", "Log In", class = "btn-login"),
          div(class = "login-toggle",
              p(HTML("New To SPECTRA? "), tags$a(id = "show_signup", onclick = "Shiny.setInputValue('toggle_form', 'signup', {priority: 'event'})", "Create Account"))
          )
      ),
      div(id = "signup_form", style = "display: none;",
          h2("Create Account"),
          textInput("signup_name", NULL, placeholder = "Full Name"),
          textInput("signup_email", NULL, placeholder = "Email Address"),
          textInput("signup_org", NULL, placeholder = "Organization (Optional)"),
          passwordInput("signup_password", NULL, placeholder = "Create Password"),
          passwordInput("signup_password2", NULL, placeholder = "Confirm Password"),
          div(style = "margin-top: 5px;",
              uiOutput("signup_msg")),
          actionButton("btn_signup", "Sign Up", class = "btn-signup"),
          div(class = "login-toggle",
              p(HTML("Already Have An Account? "), tags$a(id = "show_login", onclick = "Shiny.setInputValue('toggle_form', 'login', {priority: 'event'})", "Log In"))
          )
      ),
      div(id = "verify_form", style = "display: none;",
          h2("Email Verification"),
          div(id = "verify_info_text",
              p(style = "color: rgba(255,255,255,0.7); font-size: 0.85em;",
                "A 6-digit verification code has been sent to your email. Check your inbox and spam folder.")),
          div(id = "verify_fallback_code", style = "display: none;",
              p(style = "color: #F97316; font-size: 0.82em;",
                icon("exclamation-triangle"),
                "SMTP Not Configured. Your verification code is shown below:"),
              div(id = "verify_code_box",
                  style = "text-align: center; background: rgba(255,255,255,0.15); padding: 12px; border-radius: 8px; font-size: 1.5em; font-weight: bold; color: #F97316; letter-spacing: 6px; margin: 10px 0;",
                  "------")),
          textInput("verify_code_input", NULL, placeholder = "Enter 6-Digit Code"),
          div(style = "margin-top: 5px;",
              uiOutput("verify_msg")),
          actionButton("btn_verify", "Verify Email", class = "btn-login"),
          div(class = "login-toggle",
              p(HTML("Back To "), tags$a(id = "back_to_login", onclick = "Shiny.setInputValue('toggle_form', 'login', {priority: 'event'})", "Log In"))
          )
      )
  ),

  div(class = "login-footer",
      HTML("Developed By <a href='https://dsi-africa.org/project/23' target='_blank'>CHEAQI-MNCH</a> &nbsp;|&nbsp; ",
           "<a href='https://www.linkedin.com/company/place-alert-labs/' target='_blank'>Follow Us On LinkedIn</a>")
  )
)

main_app <- dashboardPage(
  title = "SPECTRA",
  skin = "blue",
  
  # --- Header ---
  dashboardHeader(
    title = span(
      img(src = "cheaqi_logo.png", height = "30px", style = "margin-right: 6px; vertical-align: middle;"),
      span("SPECTRA", style = "color: #F97316; font-weight: 800; font-size: 1.1em; letter-spacing: 1px;"),
      span(" · Personal Exposure Trajectory Analyzer", style = "font-size: 0.68em; font-weight: 400; opacity: 0.85;")
    ),
    titleWidth = 420
  ),
  
  # --- Sidebar (hidden — nav is horizontal in body) ---
  dashboardSidebar(
    collapsed = TRUE,
    tags$style(".main-sidebar { display:none !important; } .content-wrapper { margin-left:0 !important; }"),
    sidebarMenu(id = "tabs",
      menuItem("Home",                  tabName = "home"),
      menuItem("Input & Configuration", tabName = "config"),
      menuItem("Dashboard",             tabName = "dashboard"),
      menuItem("Spatial",               tabName = "spatial"),
      menuItem("Quality",               tabName = "quality"),
      menuItem("Segmentation",          tabName = "segmentation"),
      menuItem("Data & Export",         tabName = "dataexport"),
      menuItem("About",                 tabName = "about")
    )
  ),

  # --- Body ---
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML(app_css))),

    tabItems(

      # ====================================================================
      # TAB 0: HOME (landing page — shown on startup)
      # ====================================================================
      tabItem(tabName = "home",
        div(class = "spectra-home-wrap",

          # Watermark logo — top right, 20% opacity
          tags$img(src = "cheaqi_logo.png", class = "spectra-home-watermark"),

          # Title block — gold/white exactly like login page
          div(class = "spectra-home-logo", "SPECTRA"),
          p(class = "spectra-home-tagline",
            tags$b("S"), "patiotemporal ",
            tags$b("P"), "ersonal ",
            tags$b("E"), "xposure ",
            tags$b("C"), "haracterization and ",
            tags$b("TRA"), "jectory ",
            tags$b("A"), "nalyzer"
          ),
          p(style = "color:rgba(255,255,255,0.45); font-size:11px; margin-top:8px; letter-spacing:1px; position:relative; z-index:2;",
            "by CHEAQI-MNCH \u00b7 Supported by DSI Africa"),

          # Row 1 — 4 cards
          div(class = "spectra-home-cards",
            actionButton("card_config",
              tagList(icon("upload"), tags$span(class="home-card-label", "Input &", tags$br(), "Configuration")),
              class = "home-card-btn"),
            actionButton("card_dashboard",
              tagList(icon("chart-bar"), tags$span(class="home-card-label", "Dashboard")),
              class = "home-card-btn"),
            actionButton("card_spatial",
              tagList(icon("globe"), tags$span(class="home-card-label", "Geospatial", tags$br(), "& Tracks")),
              class = "home-card-btn"),
            actionButton("card_quality",
              tagList(icon("check-circle"), tags$span(class="home-card-label", "Quality", tags$br(), "Assessment")),
              class = "home-card-btn")
          ),

          # Row 2 — 3 cards
          div(class = "spectra-home-cards-row2",
            actionButton("card_segmentation",
              tagList(icon("layer-group"), tags$span(class="home-card-label", "Segmentation")),
              class = "home-card-btn"),
            actionButton("card_dataexport",
              tagList(icon("download"), tags$span(class="home-card-label", "Data &", tags$br(), "Export")),
              class = "home-card-btn"),
            actionButton("card_about",
              tagList(icon("info-circle"), tags$span(class="home-card-label", "About")),
              class = "home-card-btn home-card-btn-muted")
          ),

          div(style = "margin-top: 40px; color: #aaa; font-size: 11px;",
            "Click a module to begin \u00b7 Start by uploading your GPS trajectory data under Input & Configuration"
          )
        )
      ),

      # ====================================================================
      # TAB 1: INPUT & CONFIGURATION
      # ====================================================================
      tabItem(tabName = "config",

        div(class = "section-banner",
            actionButton("home_from_config", tagList(icon("home"), " Menu"),
                         class = "home-btn"),
            h2(icon("upload"), "INPUT & CONFIGURATION"),
            p("Upload your dataset, map fields, configure parameters and run analysis — preview your data on the live map")
        ),

        # Two-column layout: left = config sections, right = live map
        fluidRow(

          # ── LEFT: all config sections ─────────────────────────────────
          column(8,

            # Section 1: Data Input
            box(title = tagList(
                  tags$span(class = "config-section-num", "1"),
                  icon("upload"), " Data Input"),
                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,

                p("Upload trajectory data for initial analysis"),

                radioGroupButtons("data_source",
                  label = "File Format:",
                  choices = c("CSV" = "csv",
                              "RDS" = "rds",
                              "GeoJSON" = "geojson",
                              "GPX" = "gpx",
                              "URL" = "url"),
                  selected = "csv", status = "primary", justified = TRUE),

                conditionalPanel(
                  condition = "input.data_source == 'csv'",
                  fileInput("traj_file", "Upload CSV File",
                            accept = c(".csv"),
                            placeholder = "Drag and drop CSV here")
                ),
                conditionalPanel(
                  condition = "input.data_source == 'rds'",
                  fileInput("rds_traj_file", "Upload RDS File",
                            accept = c(".rds"),
                            placeholder = "Drag and drop .rds here"),
                  div(class = "alert alert-info", style = "font-size: 0.82em;",
                      icon("info-circle"),
                      HTML(" RDS file must contain a <b>data.frame or data.table</b> of GPS
                            trajectory points (not a session file). Column mapping is applied
                            in the next section."))
                ),
                conditionalPanel(
                  condition = "input.data_source == 'geojson'",
                  fileInput("geojson_file", "Upload GeoJSON File",
                            accept = c(".geojson", ".json"),
                            placeholder = "Drag and drop .geojson here"),
                  div(class = "alert alert-info", style = "font-size: 0.82em;",
                      icon("info-circle"),
                      HTML(" Supports Point, MultiPoint, LineString, and MultiLineString
                            features. Coordinates are extracted automatically; all properties
                            are available for field mapping below."))
                ),
                conditionalPanel(
                  condition = "input.data_source == 'gpx'",
                  fileInput("gpx_file", "Upload GPX File (wearable / GPS device)",
                            accept = c(".gpx"),
                            placeholder = "Drag and drop GPX here"),
                  div(class = "alert alert-info", style = "font-size: 0.82em;",
                      icon("info-circle"),
                      HTML(" Heart rate, cadence, speed and temperature are extracted
                            automatically. Map your sensor columns to
                            <b>Exposure 1 / Exposure 2</b> below."))
                ),
                conditionalPanel(
                  condition = "input.data_source == 'url'",
                  textInput("data_url", "Data URL",
                            placeholder = "https://... or Google Sheets share link"),
                  actionButton("load_url_btn", icon("cloud-download-alt"),
                               " Load from URL",
                               class = "btn-primary btn-sm",
                               style = "width: 100%; margin-bottom: 8px;"),
                  div(class = "alert alert-info", style = "font-size: 0.82em;",
                      icon("info-circle"),
                      HTML(" Supports any direct CSV URL or a Google Sheets
                            <em>Share \u2192 Publish to web \u2192 CSV</em> link."))
                ),
                conditionalPanel(
                  condition = "output.file_loaded",
                  uiOutput("file_info_ui"),
                  hr(),
                  verbatimTextOutput("file_preview_text")
                )
            ),

            # Session restore \u2014 separate from initial data load
            box(title = tagList(icon("history"), " Restore Previous Session"),
                width = NULL, status = "info", collapsible = TRUE, collapsed = TRUE,
                p(style = "font-size: 12px;",
                  "Upload a previously downloaded SPECTRA session (.rds) to restore all
                   analysis results instantly without re-running the pipeline."),
                fileInput("rds_session_file", NULL,
                          accept = c(".rds"),
                          placeholder = "Drag and drop session .rds here"),
                div(class = "alert alert-warning", style = "font-size: 0.82em; margin-bottom:0;",
                    icon("info-circle"),
                    HTML(" Session files are downloaded from <b>Data &amp; Export &rarr;
                          Downloads &rarr; Session RDS</b>."))
            ),

            # Section 2: Field Mapping
            box(title = tagList(
                  tags$span(class = "config-section-num", "2"),
                  icon("map-signs"), " Field Mapping"),
                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,

                p("Map your data columns to analysis fields"),

                h4(icon("map-marker-alt", style = "color: red;"), "Coordinates"),
                fluidRow(
                  column(6, uiOutput("lat_field_ui")),
                  column(6, uiOutput("lon_field_ui"))
                ),
                h4(icon("id-card"), "Identification"),
                fluidRow(
                  column(6, uiOutput("pid_field_ui")),
                  column(6, uiOutput("date_field_ui"))
                ),
                h4(icon("flask"), "Exposure Variables"),
                uiOutput("dynamic_exposure_ui"),
                actionButton("add_exposure_btn",
                  tagList(icon("plus-circle"), " Add Exposure Variable"),
                  class = "btn-sm btn-default",
                  style = "margin-bottom:10px;"),
                h4(icon("cog"), "Optional"),
                fluidRow(
                  column(6, uiOutput("truth_field_ui")),
                  column(6, numericInput("csv_epsg", "EPSG Code", value = 4326, min = 1000, max = 100000))
                ),
                fluidRow(
                  column(4, uiOutput("country_field_ui")),
                  column(4, uiOutput("temp_field_ui")),
                  column(4, uiOutput("humid_field_ui"))
                ),
                h4(icon("heartbeat"), "Wearable / Biometric (Optional)"),
                fluidRow(
                  column(4, uiOutput("heart_rate_field_ui")),
                  column(4, uiOutput("steps_field_ui")),
                  column(4, uiOutput("activity_field_ui"))
                )
            ),

            # Section 3: Analysis Parameters
            box(title = tagList(
                  tags$span(class = "config-section-num", "3"),
                  icon("sliders-h"), " Analysis Parameters"),
                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,

                p("Spatial and temporal analysis settings"),

                fluidRow(
                  column(6,
                    numericInput("cell_size", "Grid Size (m)", value = 100, min = 10, max = 10000, step = 100),
                    numericInput("max_dwell", "Max Dwell (min)", value = 20, min = 1, max = 120, step = 1),
                    numericInput("building_buffer", "Buffer (m)", value = 5, min = 0, max = 50, step = 1),
                    numericInput("speed_outdoor_threshold", "Outdoor Speed (m/min)", value = 100,
                                 min = 10, max = 500, step = 10)
                  ),
                  column(6,
                    numericInput("speed_limit", "Speed Limit (km/h)", value = 100, min = 10, max = 500, step = 10),
                    checkboxInput("remove_teleports", "Remove Teleports", value = TRUE),
                    checkboxInput("exclude_indoor", "Exclude Indoor", value = FALSE)
                  )
                ),
                hr(),
                div(class = "alert alert-info", style = "margin-bottom: 0; padding: 8px 12px;",
                  checkboxInput("use_building_validation",
                                HTML("<strong>Use Building Validation (Step 2)</strong>"),
                                value = TRUE),
                  conditionalPanel(
                    condition = "input.use_building_validation",
                    div(style = "font-size: 11px; color: #0A1A30; margin-top: 2px;",
                        icon("info-circle"),
                        " Indoor requires speed below threshold",
                        strong("AND"), "point inside a building footprint.")
                  ),
                  conditionalPanel(
                    condition = "!input.use_building_validation",
                    div(style = "font-size: 11px; color: #856404; margin-top: 2px;",
                        icon("exclamation-triangle"),
                        " Speed-only: Indoor when speed < threshold, regardless of building location.")
                  )
                )
            ),

            # Section 4: Quality Assurance
            box(title = tagList(
                  tags$span(class = "config-section-num", "4"),
                  icon("check-double"), " Quality Assurance"),
                width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,

                p("Ground-truth labels and validation"),

                checkboxInput("use_label_wizard", "Use label wizard", value = TRUE),
                conditionalPanel(
                  condition = "input.use_label_wizard",
                  h5(icon("home", style = "color: orange;"), "Indoor labels"),
                  textInput("indoor_labels", NULL,
                            value = "1,true,t,yes,y,indoor,inside,in"),
                  h5(icon("tree", style = "color: green;"), "Outdoor labels"),
                  textInput("outdoor_labels", NULL,
                            value = "0,false,f,no,n,outdoor,outside,out")
                ),
                hr(),
                h5(icon("crosshairs", style = "color: red;"), "Accuracy Parameters"),
                sliderInput("indoor_sensitivity", "Indoor Sensitivity",
                            min = 0, max = 1, value = 0.85, step = 0.01),
                sliderInput("outdoor_specificity", "Outdoor Specificity",
                            min = 0, max = 1, value = 0.99, step = 0.01)
            ),

            # Run Analysis
            box(width = NULL, status = "success", solidHeader = FALSE,
                div(style = "text-align: center; padding: 16px;",
                    actionButton("run_analysis",
                                 label = tagList(icon("rocket"), " RUN ANALYSIS"),
                                 class = "btn-success btn-lg",
                                 style = "width: 80%; font-size: 1.25em; padding: 14px; font-weight: 700;"),
                    br(), br(),
                    uiOutput("run_status_ui")
                )
            )

          ), # end col-8

          # ── RIGHT: Live Data Preview Map ──────────────────────────────
          column(4,
            div(class = "live-map-sidebar",
              box(title = tagList(icon("map-marker-alt"), " Live Data Preview"),
                  width = NULL, status = "info", solidHeader = TRUE,
                  selectInput("live_map_style", "Map Style",
                              choices = c("CartoDB Positron" = "CartoDB.Positron",
                                          "OpenStreetMap"   = "OpenStreetMap",
                                          "CartoDB Dark"    = "CartoDB.DarkMatter",
                                          "Esri Satellite"  = "Esri.WorldImagery"),
                              selected = "CartoDB.Positron"),
                  checkboxInput("show_buildings_live", "Show Building Polygons", FALSE),
                  shinycssloaders::withSpinner(
                    leafletOutput("live_map", height = "480px"),
                    type = 4, color = "#4a90e2"
                  ),
                  hr(),
                  uiOutput("live_map_info_ui")
              )
            )
          ) # end col-4

        ) # end fluidRow
      ), # end tabItem config

      # (livemap tab removed — content is now in the config tab right column)

      # ====================================================================
      # TAB: RUN ANALYSIS (hidden — run button is in config tab)
      # ====================================================================
      tabItem(tabName = "run",
        div(class = "section-banner",
            actionButton("home_from_run", tagList(icon("home"), " Menu"), class = "home-btn"),
            h2(icon("play-circle"), "ANALYSIS EXECUTION"),
            p("Monitor analysis progress")
        ),
        
        fluidRow(
          column(8, offset = 2,
            box(width = NULL, status = "primary",
                div(style = "text-align: center; padding: 30px;",
                    h3("Analysis Pipeline"),
                    p("1. Load & validate data"),
                    p("2. Apply boundary filter"),
                    p("3. Calculate dwell time & speed"),
                    p("4. Fetch/load building footprints"),
                    p("5. Classify indoor/outdoor"),
                    p("6. Aggregate to grid cells"),
                    p("7. Generate statistics"),
                    hr(),
                    actionButton("run_analysis2",
                                 label = tagList(icon("rocket"), "RUN ANALYSIS"),
                                 class = "btn-primary btn-lg",
                                 style = "width: 60%; font-size: 1.3em; padding: 15px;"),
                    br(), br(),
                    uiOutput("run_progress_ui")
                )
            )
          )
        )
      ),
      
      # ====================================================================
      # TAB 4: DASHBOARD
      # ====================================================================
      tabItem(tabName = "dashboard",

        div(class = "section-banner",
            actionButton("home_from_dashboard", tagList(icon("home"), " Menu"), class = "home-btn"),
            h2(icon("chart-bar"), "ANALYSIS OUTPUTS DASHBOARD"),
            p("Comprehensive Results & Key Performance Indicators")
        ),
        
        # Key Metrics Row
        fluidRow(
          infoBoxOutput("metric_classified", width = 2),
          infoBoxOutput("metric_analysis", width = 2),
          infoBoxOutput("metric_people", width = 2),
          infoBoxOutput("metric_indoor", width = 3),
          infoBoxOutput("metric_dwell", width = 3)
        ),
        
        # Filters
        fluidRow(
          column(8,
            box(title = "Interactive Filters", width = NULL, status = "info",
                         collapsible = TRUE, collapsed = FALSE,
                fluidRow(
                  column(6, uiOutput("filter_pid_ui")),
                  column(6, uiOutput("filter_date_ui"))
                )
            )
          ),
          column(4,
            box(title = "Filtered Summary", width = NULL, status = "success",
                uiOutput("filtered_summary_ui")
            )
          )
        ),
        
        # Charts Row
        fluidRow(
          column(6,
            box(title = "Dwell Time Distribution", width = NULL, status = "primary",
                shinycssloaders::withSpinner(
                  plotlyOutput("plot_dwell_hist", height = "400px"),
                  type = 4, color = "#4a90e2"
                )
            )
          ),
          column(6,
            box(title = "Indoor/Outdoor Split", width = NULL, status = "primary",
                shinycssloaders::withSpinner(
                  plotlyOutput("plot_indoor_pie", height = "400px"),
                  type = 4, color = "#4a90e2"
                )
            )
          )
        ),
        
        # Dwell per person
        fluidRow(
          column(5,
            box(title = "Dwell Time Per Participant", width = NULL, status = "primary",
                DT::dataTableOutput("dwell_per_person_table")
            )
          ),
          column(7,
            box(title = "Top Participants by Dwell Time", width = NULL, status = "primary",
                shinycssloaders::withSpinner(
                  plotlyOutput("plot_dwell_bar", height = "400px"),
                  type = 4, color = "#4a90e2"
                )
            )
          )
        )
      ),
      
      # ====================================================================
      # TAB 5: SPATIAL VIEW
      # ====================================================================
      tabItem(tabName = "spatial",

        div(class = "section-banner",
            actionButton("home_from_spatial", tagList(icon("home"), " Menu"), class = "home-btn"),
            h2(icon("globe"), "GEOSPATIAL ANALYSIS"),
            p("Interactive map with classification, exposure levels, raster grids, and time animation")
        ),
        
        tabsetPanel(
          tabPanel("Point Map",
            br(),
            fluidRow(
              column(9,
                box(title = "Trajectory Map", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      leafletOutput("spatial_map", height = "650px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(3,
                box(title = "Map Controls", width = NULL, status = "info",
                    selectInput("spatial_map_style", "Map Style",
                                choices = c("CartoDB Positron" = "CartoDB.Positron",
                                            "OpenStreetMap" = "OpenStreetMap",
                                            "CartoDB Dark" = "CartoDB.DarkMatter",
                                            "Esri Imagery" = "Esri.WorldImagery"),
                                selected = "CartoDB.Positron"),
                    
                    radioButtons("spatial_color_by", "Color Points By:",
                                 choices = c("Classification", "Exposure 1", "Exposure 2", "Speed")),
                    
                    sliderInput("spatial_sample", "Max Points",
                                min = 500, max = 20000, value = 5000, step = 500),
                    
                    hr(),
                    checkboxInput("show_heatmap",   "Show Density Heatmap", FALSE),
                    checkboxInput("show_buildings_spatial", "Show Building Polygons", FALSE)
                )
              )
            )
          ),

          tabPanel("Raster Grid (Dwell & Exposure)",
            br(),
            fluidRow(
              column(9,
                box(title = "Grid Cell Raster Map", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      leafletOutput("raster_grid_map", height = "650px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(3,
                box(title = "Raster Controls", width = NULL, status = "info",
                    radioButtons("raster_variable", "Grid Variable:",
                                 choices = c("Dwell Time (min)" = "dwell",
                                             "Exposure 1" = "pm25",
                                             "Exposure 2" = "no2",
                                             "Indoor %" = "indoor_pct")),
                    uiOutput("raster_date_ui")
                )
              )
            )
          ),

          tabPanel("Trajectory Lines",
            br(),
            fluidRow(
              column(9,
                box(title = "Movement Lines", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      leafletOutput("trajectory_map", height = "600px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(3,
                box(title = "Trajectory Controls", width = NULL, status = "info",
                    uiOutput("traj_pid_select_ui"),
                    radioButtons("traj_color_by", "Color Lines By:",
                                 choices = c("Participant", "Day", "Indoor/Outdoor")),
                    checkboxInput("traj_show_points", "Show GPS Points", TRUE),
                    selectInput("traj_map_style", "Map Style",
                                choices = c("CartoDB Positron" = "CartoDB.Positron",
                                            "OpenStreetMap" = "OpenStreetMap",
                                            "CartoDB Dark" = "CartoDB.DarkMatter"),
                                selected = "CartoDB.Positron"),
                    hr(),
                    uiOutput("traj_info_ui")
                )
              )
            )
          ),

          tabPanel("Time Animation",
            br(),
            fluidRow(
              column(9,
                box(title = "Exposure Animation", width = NULL, status = "primary",
                    uiOutput("animation_video_ui"),
                    div(style = "text-align:center; margin-top:10px;",
                      actionButton("btn_render_animation", tagList(icon("film"), " Render Animation"),
                                   class = "btn-primary", style = "margin-right:8px;"),
                      uiOutput("anim_render_status_ui")
                    )
                )
              ),
              column(3,
                box(title = "Animation Settings", width = NULL, status = "info",
                    radioButtons("anim_color", "Color By:",
                                 choices = c("Exposure 1" = "pm25",
                                             "Exposure 2" = "no2",
                                             "Classification" = "class")),
                    sliderInput("anim_duration", "Max Duration (sec)", min = 10, max = 60, value = 30, step = 5),
                    sliderInput("anim_fps", "Frames Per Second", min = 2, max = 10, value = 5, step = 1),
                    hr(),
                    p(style = "font-size:11px; color:#666;",
                      icon("info-circle"),
                      " Click Render to build the animation. Once ready, press Play in the video player above.")
                )
              )
            )
          )
        )
      ),
      
      
      # ====================================================================
      # TAB 7: QUALITY ASSESSMENT
      # ====================================================================
      tabItem(tabName = "quality",

        div(class = "section-banner",
            actionButton("home_from_quality", tagList(icon("home"), " Menu"), class = "home-btn"),
            h2(icon("check-circle"), "QUALITY ASSESSMENT"),
            p("Model performance, confusion matrix and per-participant metrics")
        ),
        
        fluidRow(
          column(6,
            box(title = "Overall Performance Metrics", width = NULL, status = "primary",
                uiOutput("qa_overall_ui")
            )
          ),
          column(6,
            box(title = "Confusion Matrix Heatmap", width = NULL, status = "primary",
                shinycssloaders::withSpinner(
                  plotlyOutput("plot_confusion_matrix", height = "400px"),
                  type = 4, color = "#4a90e2"
                )
            )
          )
        ),
        
        fluidRow(
          column(6,
            box(title = "Per-Participant Metrics", width = NULL, status = "primary",
                DT::dataTableOutput("qa_per_person_table")
            )
          ),
          column(6,
            box(title = "Per-Participant Accuracy Chart", width = NULL, status = "primary",
                shinycssloaders::withSpinner(
                  plotlyOutput("plot_per_person_accuracy", height = "400px"),
                  type = 4, color = "#4a90e2"
                )
            )
          )
        )
      ),
      
      # ====================================================================
      # TAB 8: SEGMENTATION & VALIDATION
      # ====================================================================
      tabItem(tabName = "segmentation",

        div(class = "section-banner",
            actionButton("home_from_segmentation", tagList(icon("home"), " Menu"), class = "home-btn"),
            h2(icon("layer-group"), "SEGMENTATION & VALIDATION"),
            p("Indoor/outdoor classification with temporal, mobility and climate segmentation")
        ),

        tabsetPanel(

          # --- Sub-tab 1: Validation Pipeline ---
          tabPanel("Validation Pipeline",
            br(),
            fluidRow(
              column(6,
                box(title = "Validation Pipeline Statistics", width = NULL, status = "primary",
                    DT::dataTableOutput("seg_validation_stats_table")
                )
              ),
              column(6,
                box(title = "Speed vs OSM Confusion Matrix", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_confusion_heatmap", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            ),
            fluidRow(
              column(6,
                box(title = "Misclassification Breakdown", width = NULL, status = "info",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_misclass_pie", height = "350px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(6,
                box(title = "Speed-Only vs Corrected Comparison", width = NULL, status = "info",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_method_comparison", height = "350px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            )
          ),

          # --- Sub-tab 2: Indoor/Outdoor Distribution ---
          tabPanel("Indoor/Outdoor Distribution",
            br(),
            fluidRow(
              column(6,
                box(title = "Corrected Indoor vs Outdoor", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_indoor_outdoor_bar", height = "450px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(6,
                box(title = "Exposure Type Over Time", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_timeseries", height = "450px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            )
          ),

          # --- Sub-tab 3: Outdoor Analysis ---
          tabPanel("Outdoor Analysis",
            br(),
            fluidRow(
              column(6,
                box(title = "Outdoor Movement Intensity", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_outdoor_movement", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(6,
                box(title = "Outdoor Exposure Context", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_outdoor_context", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            ),
            fluidRow(
              column(6,
                box(title = "Outdoor Environmental Stress", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_outdoor_stress", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(6,
                box(title = "Outdoor Corrected Data Table", width = NULL, status = "info",
                    DT::dataTableOutput("seg_outdoor_table")
                )
              )
            ),
            fluidRow(
              column(12,
                box(title = "Download Corrected Datasets", width = NULL, status = "success",
                    fluidRow(
                      column(4,
                        downloadButton("download_outdoor_corrected", "Download Outdoor Corrected (CSV)",
                                       class = "btn-primary", style = "width: 100%;")
                      ),
                      column(4,
                        downloadButton("download_indoor_corrected", "Download Indoor Corrected (CSV)",
                                       class = "btn-primary", style = "width: 100%;")
                      ),
                      column(4,
                        downloadButton("download_full_segmented", "Download Full Segmented (CSV)",
                                       class = "btn-primary", style = "width: 100%;")
                      )
                    )
                )
              )
            )
          ),

          # --- Sub-tab 4: Temporal & Mobility ---
          tabPanel("Temporal & Mobility",
            br(),
            fluidRow(
              column(12,
                box(title = "Exposure Variable", width = NULL, status = "info",
                    collapsible = TRUE, collapsed = FALSE,
                    fluidRow(
                      column(4,
                        selectInput("seg_exposure", "Select Exposure Variable",
                                    choices = c("Exposure 1" = "pm25", "Exposure 2" = "no2"),
                                    selected = "pm25")
                      ),
                      column(8,
                        p(style = "margin-top: 28px; color: #666;",
                          icon("info-circle"),
                          "Charts below show mean exposure levels across different themes to identify peak exposure periods and patterns.")
                      )
                    )
                )
              )
            ),
            fluidRow(
              column(6,
                box(title = "Exposure by Time of Day", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_temporal_pattern", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(6,
                box(title = "Exposure by Mobility Theme", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_mobility_theme", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            ),
            fluidRow(
              column(6,
                box(title = "Weekday x Time-of-Day Exposure Heatmap", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_weekday_heatmap", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(6,
                box(title = "Exposure by Duration Theme", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_duration_theme", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            ),
            fluidRow(
              column(6,
                box(title = "Hourly Exposure Profile", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_hourly_profile", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(6,
                box(title = "Peak Exposure Hours Summary", width = NULL, status = "warning",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_peak_hours", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            )
          ),

          # --- Sub-tab 5: Climate Themes ---
          tabPanel("Climate Themes",
            br(),
            fluidRow(
              column(6,
                box(title = "Climate Conditions by Exposure Type", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_climate_exposure", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              ),
              column(6,
                box(title = "Climate Theme Distribution", width = NULL, status = "primary",
                    shinycssloaders::withSpinner(
                      plotlyOutput("seg_climate_pie", height = "400px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            )
          ),

          # --- Sub-tab 6: Misclassification Map ---
          tabPanel("Misclassification Map",
            br(),
            fluidRow(
              column(12,
                box(title = "Spatial Misclassification Map", width = NULL, status = "primary",
                    fluidRow(
                      column(6, p("Green = Correctly classified, Orange = False Indoor, Purple = False Outdoor")),
                      column(6, checkboxInput("show_buildings_misclass",
                                              "Show Building Polygons", FALSE))
                    ),
                    shinycssloaders::withSpinner(
                      leafletOutput("seg_misclass_map", height = "600px"),
                      type = 4, color = "#4a90e2"
                    )
                )
              )
            )
          )
        )
      ),

      # ====================================================================
      # TAB 9: DATA & EXPORT
      # ====================================================================
      tabItem(tabName = "dataexport",

        div(class = "section-banner",
            actionButton("home_from_dataexport", tagList(icon("home"), " Menu"), class = "home-btn"),
            h2(icon("download"), "DATA & EXPORT"),
            p("Browse, review and download all analysis results")
        ),

        tabsetPanel(
          tabPanel("Data Tables",
            br(),
            tabsetPanel(
              tabPanel("Classified Points",    br(), DT::dataTableOutput("table_classified")),
              tabPanel("Grid Cells",           br(), DT::dataTableOutput("table_grid")),
              tabPanel("Daily Summaries",      br(), DT::dataTableOutput("table_daily")),
              tabPanel("Participant Exposure", br(),
                p(style = "color:#0A1A30; font-size:12px; margin-bottom:8px;",
                  icon("info-circle"),
                  "Average daily exposure per participant (mean of daily grid-averaged values) with WHO Air Quality Guideline (2021) comparison."),
                DT::dataTableOutput("table_participant_exposure")
              )
            )
          ),
          tabPanel("Downloads",
            br(),

            # ── Primary Downloads: RDS + NC (fast, recommended) ───────────────
            fluidRow(
              column(12,
                box(title = tagList(icon("star"), " Primary Downloads — Recommended"),
                    width = NULL, status = "primary",
                    p(style = "color:#0A1A30; font-size:12px;",
                      "These two files are all you need. The RDS restores your session; the
                       NetCDF grid is the spatial output for modelling and GIS workflows."),
                    fluidRow(
                      column(6,
                        div(style = paste0("border:2px solid #0077B6; border-radius:8px;",
                                           " padding:14px; background:#f0f6ff;"),
                          h4(style = "margin-top:0; color:#0A1A30;",
                             icon("save"), " Grid RDS"),
                          tags$ul(style = "font-size:11.5px; color:#333; padding-left:16px;",
                            tags$li("Full spatial grid as R data.table"),
                            tags$li("Convertible to NC or CSV later"),
                            tags$li("Columns: cell_id, lon, lat, dwell, pm25, no2, is_labeled, country"),
                            tags$li("Compact — gzip compressed")
                          ),
                          downloadButton("download_result_rds",
                                         tagList(icon("download"), " Download Grid (.rds)"),
                                         class = "btn-primary",
                                         style = "width:100%; font-weight:bold;")
                        )
                      ),
                      column(6,
                        div(style = paste0("border:2px solid #00B4D8; border-radius:8px;",
                                           " padding:14px; background:#f0fbff;"),
                          h4(style = "margin-top:0; color:#0A1A30;",
                             icon("cube"), " Exposure Grid NetCDF"),
                          tags$ul(style = "font-size:11.5px; color:#333; padding-left:16px;",
                            tags$li("Full spatial extent — ALL grid cells"),
                            tags$li(HTML("<b>is_labeled=1</b>: cells with GPS data (labeled)")),
                            tags$li(HTML("<b>is_labeled=0</b>: empty cells (unlabeled dataset)")),
                            tags$li("Variables: dwell_time, exposure_1/2, indoor_pct"),
                            tags$li("CF-1.8 compliant, compatible with Python/QGIS/R")
                          ),
                          downloadButton("download_raster_nc",
                                         tagList(icon("download"), " Download Grid (.nc)"),
                                         class = "btn-info",
                                         style = "width:100%; font-weight:bold;")
                        )
                      )
                    ),
                    fluidRow(
                      column(12,
                        div(style = paste0("border:2px solid #2d8a4e; border-radius:8px;",
                                           " padding:14px; background:#f0fff4; margin-top:10px;"),
                          h4(style = "margin-top:0; color:#0A1A30;",
                             icon("images"), " All Figures Bundle"),
                          tags$ul(style = "font-size:11.5px; color:#333; padding-left:16px;",
                            tags$li("All key charts exported as high-resolution PNG"),
                            tags$li("Indoor/outdoor, exposure profiles, heatmaps, QA plots"),
                            tags$li("ZIP archive — ready to use in reports or presentations")
                          ),
                          downloadButton("download_all_figures",
                                         tagList(icon("download"), " Download All Figures (.zip)"),
                                         class = "btn-success",
                                         style = "width:100%; font-weight:bold;")
                        )
                      )
                    )
                )
              )
            ),

            # ── Advanced Exports (collapsible) ────────────────────────────────
            fluidRow(
              column(12,
                box(title = tagList(icon("archive"), " Advanced Exports (CSV / GeoJSON / Excel / HTML)"),
                    width = NULL, status = "primary",
                    collapsible = TRUE, collapsed = TRUE,
                    fluidRow(
                      column(4,
                        strong("Data Tables (CSV)"),
                        br(), br(),
                        downloadButton("download_classified_csv", "Classified Points",
                                       class = "btn-primary btn-sm",
                                       style = "width:100%; margin-bottom:6px;"),
                        downloadButton("download_grid_csv", "Grid Cells (Observed)",
                                       class = "btn-primary btn-sm",
                                       style = "width:100%; margin-bottom:6px;"),
                        downloadButton("download_full_grid_csv", "Full Grid (All Cells)",
                                       class = "btn-success btn-sm",
                                       style = "width:100%; margin-bottom:6px;"),
                        downloadButton("download_daily_csv", "Daily Summaries",
                                       class = "btn-primary btn-sm",
                                       style = "width:100%; margin-bottom:6px;"),
                        downloadButton("download_participant_csv", "Participant Exposure + WHO",
                                       class = "btn-success btn-sm",
                                       style = "width:100%; margin-bottom:6px;"),
                        downloadButton("download_raster_csv", "Full Grid (CSV)",
                                       class = "btn-info btn-sm",
                                       style = "width:100%;")
                      ),
                      column(4,
                        strong("Spatial / Summary"),
                        br(), br(),
                        downloadButton("download_points_geojson", "Points (GeoJSON)",
                                       class = "btn-primary btn-sm",
                                       style = "width:100%; margin-bottom:6px;"),
                        downloadButton("download_stats_json", "Statistics (JSON)",
                                       class = "btn-primary btn-sm",
                                       style = "width:100%; margin-bottom:6px;"),
                        downloadButton("download_summary", "Summary Report (CSV)",
                                       class = "btn-primary btn-sm",
                                       style = "width:100%;")
                      ),
                      column(4,
                        strong("Excel / HTML Report"),
                        br(), br(),
                        downloadButton("download_master_excel",
                                       tagList(icon("file-excel"), " Master Excel (.xlsx)"),
                                       class = "btn-success btn-sm",
                                       style = "width:100%; margin-bottom:6px;"),
                        downloadButton("download_full_report",
                                       tagList(icon("file-alt"), " Full Report (.html)"),
                                       class = "btn-warning btn-sm",
                                       style = "width:100%;")
                      )
                    )
                )
              )
            )
          ) # end tabPanel("Downloads")
        ) # end outer tabsetPanel
      ), # end tabItem dataexport

      # ── ABOUT ─────────────────────────────────────────────────────────────
      tabItem(tabName = "about",

        div(class = "section-banner",
            actionButton("home_from_about", tagList(icon("home"), " Menu"), class = "home-btn"),
            h2(icon("info-circle"), "ABOUT SPECTRA"),
            p("Spatiotemporal Personal Exposure Characterization and TRAjectory Analyzer")
        ),

        fluidRow(
          column(8,
            box(title = tagList(icon("microscope"), " About the Platform"),
                width = NULL, status = "primary", solidHeader = TRUE,
                tags$p(style = "font-size:14px; line-height:1.75;",
                  "SPECTRA (Spatiotemporal Personal Exposure Characterization and TRAjectory Analyzer) is an
                  open-source software platform for personal environmental exposure assessment, developed by
                  the CHEAQI-MNCH research group with support from DSI Africa."
                ),
                tags$p(style = "font-size:14px; line-height:1.75;",
                  "SPECTRA addresses a critical methodological gap in environmental health research: the need
                  for reproducible, fine-grained tools that link individual mobility trajectories to real-time
                  environmental measurements across dynamic microenvironments. At its core, SPECTRA processes
                  GPS trajectory data fused with co-located environmental sensor measurements to generate
                  spatially and temporally resolved exposure profiles at the person-day level."
                ),
                tags$p(style = "font-size:14px; line-height:1.75;",
                  "The platform\u2019s indoor/outdoor classification engine applies a validated four-stage pipeline
                  \u2014 integrating instantaneous speed analysis, displacement thresholds, and OpenStreetMap building
                  footprint intersections \u2014 to assign environmental context to each trajectory point.
                  Exposure is then characterised across five segmentation dimensions: temporal patterns,
                  mobility mode, activity duration, meteorological conditions, and outdoor microenvironment type."
                ),
                tags$p(style = "font-size:14px; line-height:1.75;",
                  "SPECTRA is deployed as an interactive R Shiny web application with a user authentication system.
                  Outputs include time-weighted mean exposure statistics, interactive spatial visualisations,
                  raster grids, and exportable datasets in CSV, GeoJSON, and NetCDF formats. The platform is
                  designed to be domain-extensible, accommodating a wide range of environmental stressors
                  including chemical, physical, and biological agents."
                ),
                tags$p(style = "font-size:14px; line-height:1.75;",
                  "It provides a standardised, citable methodology for exposure assessment in urban epidemiology,
                  environmental justice research, active travel studies, and child and maternal health
                  investigations in low- and middle-income settings."
                )
            ),

            box(title = tagList(icon("cogs"), " Key Features"),
                width = NULL, status = "info", solidHeader = TRUE,
                fluidRow(
                  column(6,
                    tags$ul(style = "font-size:13px; line-height:1.9;",
                      tags$li("GPS trajectory ingestion & preprocessing"),
                      tags$li("Indoor/outdoor classification (4-stage pipeline)"),
                      tags$li("OpenStreetMap building footprint integration"),
                      tags$li("Five-dimension exposure segmentation")
                    )
                  ),
                  column(6,
                    tags$ul(style = "font-size:13px; line-height:1.9;",
                      tags$li("Interactive spatial visualisations & live map"),
                      tags$li("CSV, GeoJSON, NetCDF, and Excel exports"),
                      tags$li("Raster grid & daily summary generation"),
                      tags$li("User authentication & multi-session support")
                    )
                  )
                )
            )
          ),

          column(4,
            box(title = tagList(icon("users"), " Development Team"),
                width = NULL, status = "warning", solidHeader = TRUE,
                div(style = "text-align:center; padding: 10px 0;",
                  tags$img(src = "cheaqi_logo.png", height = "60px",
                           style = "margin-bottom:12px; display:block; margin-left:auto; margin-right:auto;"),
                  tags$h4("CHEAQI-MNCH", style = "color:#1a3a6e; font-weight:900; margin:4px 0;"),
                  tags$p(style = "font-size:11px; color:#666; margin:0; line-height:1.5;",
                    "Characterizing Effects of Air Quality in Maternal,",
                    tags$br(), "Newborn and Child Health (CHEAQI\u2013MNCH)",
                    tags$br(), "in sub-Saharan Africa")
                ),
                hr(),
                tags$p(tags$strong("Lead Developer:"),
                       tags$br(), "Handsome Bongani Nyoni"),
                tags$p(tags$strong("Supported by:"),
                       tags$br(), "DSI Africa"),
                hr(),
                tags$h5(icon("envelope"), " Contact", style = "font-weight:700;"),
                tags$p(
                  tags$a(href = "mailto:nyonih@staff.msu.ac.zw",
                         icon("envelope"), " nyonih@staff.msu.ac.zw",
                         style = "color:#1a3a6e;"),
                  tags$br(),
                  tags$a(href = "mailto:nynhan001@myuct.ac.za",
                         icon("envelope"), " nynhan001@myuct.ac.za",
                         style = "color:#1a3a6e;")
                ),
                hr(),
                tags$h5(icon("share-alt"), " Connect", style = "font-weight:700;"),
                tags$p(
                  tags$a(href = "https://www.linkedin.com/in/bongani-handsome-nyoni-62a132119/",
                         target = "_blank",
                         tags$span(icon("linkedin"), " Handsome Bongani Nyoni",
                                   style = "color:#0077b5; font-weight:600;"))
                )
            ),

            box(title = tagList(icon("balance-scale"), " License & Citation"),
                width = NULL, status = "info",
                tags$p(style = "font-size:12px;",
                  tags$strong("License:"), " Open Source",
                  tags$br(),
                  tags$strong("Version:"), " SPECTRA 1.0",
                  tags$br(),
                  tags$strong("Year:"), " 2025"
                ),
                tags$p(style = "font-size:12px; color:#555;",
                  "When using SPECTRA in publications, please cite the platform and the CHEAQI-MNCH research group."
                )
            )
          )
        )
      ) # end tabItem about

    ) # end tabItems
  ) # end dashboardBody
)

ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML(app_css))),
  uiOutput("page_ui")
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {

  # Increase max upload size to 1GB
  options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)

  # Initialise building cache DB on startup
  tryCatch(init_building_cache(BUILDING_CACHE_DB), error = function(e) NULL)

  # Buildings are now pre-loaded from buildings_preloaded.csv — no UI management needed

  # ====================================================================
  # Authentication State
  # ====================================================================
  auth <- reactiveValues(
    logged_in = FALSE,
    user = NULL
  )

  # Pending verification state
  pending_verify <- reactiveValues(email = NULL, code = NULL)

  # Toggle between login, signup, and verify forms
  observeEvent(input$toggle_form, {
    if (input$toggle_form == "signup") {
      shinyjs::hide("login_form")
      shinyjs::hide("verify_form")
      shinyjs::show("signup_form")
    } else if (input$toggle_form == "verify") {
      shinyjs::hide("login_form")
      shinyjs::hide("signup_form")
      shinyjs::show("verify_form")
    } else {
      shinyjs::show("login_form")
      shinyjs::hide("signup_form")
      shinyjs::hide("verify_form")
    }
  })

  # Login handler
  observeEvent(input$btn_login, {
    email <- trimws(input$login_email)
    password <- input$login_password
    if (email == "" || password == "") {
      output$login_msg <- renderUI(div(style = "color: #ff6b6b;", "Please enter email and password"))
      return()
    }
    user <- verify_user(email, password)
    if (is.list(user) && !is.null(user$error) && user$error == "email_not_verified") {
      # Resend verification email and redirect to verification form
      pending_verify$email <- email
      pending_verify$code <- user$code
      # Look up user name for the email
      con_tmp <- dbConnect(SQLite(), "spectra_users.sqlite")
      uname <- dbGetQuery(con_tmp, "SELECT full_name FROM users WHERE email = ?", params = list(email))$full_name
      dbDisconnect(con_tmp)
      email_result <- send_verification_email(email, uname %||% "User", user$code)
      if (email_result$success) {
        shinyjs::hide("verify_fallback_code")
        shinyjs::show("verify_info_text")
      } else {
        shinyjs::html("verify_code_box", as.character(user$code))
        shinyjs::show("verify_fallback_code")
        shinyjs::hide("verify_info_text")
      }
      shinyjs::hide("login_form")
      shinyjs::show("verify_form")
      return()
    }
    if (!is.null(user) && is.data.frame(user)) {
      auth$logged_in <- TRUE
      auth$user <- user
    } else {
      output$login_msg <- renderUI(div(style = "color: #ff6b6b;", "Invalid email or password"))
    }
  })

  # Signup handler
  observeEvent(input$btn_signup, {
    name <- trimws(input$signup_name)
    email <- trimws(input$signup_email)
    org <- trimws(input$signup_org)
    pw1 <- input$signup_password
    pw2 <- input$signup_password2
    if (name == "" || email == "" || pw1 == "") {
      output$signup_msg <- renderUI(div(style = "color: #ff6b6b;", "Please fill in all required fields"))
      return()
    }
    if (!grepl("@", email)) {
      output$signup_msg <- renderUI(div(style = "color: #ff6b6b;", "Please enter a valid email"))
      return()
    }
    if (nchar(pw1) < 6) {
      output$signup_msg <- renderUI(div(style = "color: #ff6b6b;", "Password must be at least 6 characters"))
      return()
    }
    if (pw1 != pw2) {
      output$signup_msg <- renderUI(div(style = "color: #ff6b6b;", "Passwords do not match"))
      return()
    }
    result <- register_user(email, pw1, name, org)
    if (result$success) {
      # Try sending the verification code by email
      pending_verify$email <- email
      pending_verify$code <- result$code
      email_result <- send_verification_email(email, name, result$code)
      if (email_result$success) {
        # Email sent — hide code from screen
        shinyjs::hide("verify_fallback_code")
        shinyjs::show("verify_info_text")
      } else {
        # Email failed — show code on screen as fallback
        shinyjs::html("verify_code_box", as.character(result$code))
        shinyjs::show("verify_fallback_code")
        shinyjs::hide("verify_info_text")
      }
      shinyjs::hide("signup_form")
      shinyjs::show("verify_form")
    } else {
      output$signup_msg <- renderUI(div(style = "color: #ff6b6b;", result$msg))
    }
  })

  # Email verification handler
  observeEvent(input$btn_verify, {
    code_input <- trimws(input$verify_code_input)
    if (code_input == "" || is.null(pending_verify$email)) {
      output$verify_msg <- renderUI(div(style = "color: #ff6b6b;", "Please enter the verification code"))
      return()
    }
    if (verify_email_code(pending_verify$email, code_input)) {
      output$verify_msg <- renderUI(div(style = "color: #2ecc71;", icon("check-circle"), " Email verified successfully! You can now log in."))
      shinyjs::delay(1500, {
        shinyjs::hide("verify_form")
        shinyjs::show("login_form")
        pending_verify$email <- NULL
        pending_verify$code <- NULL
      })
    } else {
      output$verify_msg <- renderUI(div(style = "color: #ff6b6b;", "Invalid verification code. Please try again."))
    }
  })

  # Page routing
  output$page_ui <- renderUI({
    if (auth$logged_in) {
      main_app
    } else {
      login_page
    }
  })
  
  # ── Horizontal nav: switch active tab when nav buttons clicked ─────────────
  # ── Tab switching via direct DOM manipulation ─────────────────────────────
  # updateTabItems() fails when main_app is rendered inside renderUI because
  # shinydashboard's JS binding is not re-initialised by Shiny.bindAll().
  # Instead we directly activate the tab pane div via shinyjs::runjs.
  go_tab <- function(tabId) {
    shinyjs::runjs(paste0(
      "$('.tab-pane').removeClass('active');",
      "$('#shiny-tab-", tabId, "').addClass('active');"
    ))
  }

  # ── Home landing card navigation ──────────────────────────────────────────
  observeEvent(input$card_config,       go_tab("config"))
  observeEvent(input$card_dashboard,    go_tab("dashboard"))
  observeEvent(input$card_spatial,      go_tab("spatial"))
  observeEvent(input$card_quality,      go_tab("quality"))
  observeEvent(input$card_segmentation, go_tab("segmentation"))
  observeEvent(input$card_dataexport,   go_tab("dataexport"))
  observeEvent(input$card_about,        go_tab("about"))

  # ── Back-to-home buttons (one per tab) ────────────────────────────────────
  observeEvent(input$home_from_config,       go_tab("home"))
  observeEvent(input$home_from_run,          go_tab("home"))
  observeEvent(input$home_from_dashboard,    go_tab("home"))
  observeEvent(input$home_from_spatial,      go_tab("home"))
  observeEvent(input$home_from_quality,      go_tab("home"))
  observeEvent(input$home_from_segmentation, go_tab("home"))
  observeEvent(input$home_from_dataexport,   go_tab("home"))
  observeEvent(input$home_from_about,        go_tab("home"))

  # Start on home landing tab (after a brief delay to let the DOM render)
  shinyjs::delay(200, go_tab("home"))

  # ====================================================================
  # Reactive Values
  # ====================================================================
  rv <- reactiveValues(
    raw_data        = NULL,   # Loaded data.table
    columns         = NULL,   # Column names
    analysis_result = NULL,   # Full analysis result list
    truth_field     = NULL,   # Ground truth column name
    status_msg      = "",     # Status message
    buildings_sf    = NULL,   # Buildings from buildings.csv (filtered to data bbox)
    data_source     = "csv",  # "csv" | "gpx" | "url"
    exp1_label      = "",
    exp2_label      = "",
    # Dynamic exposure variables: list of list(label, col)
    n_exposures     = 2L,
    exp_labels      = c("", ""),
    exp_cols        = c(NA_character_, NA_character_)
  )

  # ── Dynamic Exposure Variable UI ─────────────────────────────────────────
  observeEvent(input$add_exposure_btn, {
    rv$n_exposures <- rv$n_exposures + 1L
    rv$exp_labels  <- c(rv$exp_labels, paste0("Exposure ", rv$n_exposures))
    rv$exp_cols    <- c(rv$exp_cols,   NA_character_)
  })

  output$dynamic_exposure_ui <- renderUI({
    cols <- c("(none)", rv$columns %||% character(0))
    n    <- rv$n_exposures
    rows <- lapply(seq_len(n), function(i) {
      fluidRow(
        column(5,
          textInput(paste0("exp_label_", i),
                    paste0("Exposure ", i, " — Name"),
                    value = rv$exp_labels[i])
        ),
        column(5,
          selectInput(paste0("exp_col_", i),
                      "Column",
                      choices = cols,
                      selected = if (!is.na(rv$exp_cols[i])) rv$exp_cols[i] else "(none)")
        ),
        if (i > 2) column(2,
          div(style = "margin-top:25px;",
            actionButton(paste0("rm_exp_", i), icon("times"),
                         class = "btn-xs btn-danger"))
        )
      )
    })
    tagList(rows)
  })

  # ====================================================================
  # File / URL Loading
  # ====================================================================

  # Helper: load buildings using the server-level RAM cache (tiles loaded once, reused across sessions)
  .load_buildings_for_data <- function(dt) {
    lat_col <- if ("lat" %in% names(dt)) "lat" else if ("latitude" %in% names(dt)) "latitude" else return(NULL)
    lon_col <- if ("lon" %in% names(dt)) "lon" else if ("longitude" %in% names(dt)) "longitude" else return(NULL)

    lat_vals <- dt[[lat_col]]; lat_vals <- lat_vals[is.finite(lat_vals)]
    lon_vals <- dt[[lon_col]]; lon_vals <- lon_vals[is.finite(lon_vals)]
    if (length(lat_vals) == 0) return(NULL)

    buf    <- 0.02
    minlat <- min(lat_vals) - buf;  maxlat <- max(lat_vals) + buf
    minlon <- min(lon_vals) - buf;  maxlon <- max(lon_vals) + buf
    message(sprintf(
      "Buildings lookup bbox: lat %.5f..%.5f, lon %.5f..%.5f",
      minlat, maxlat, minlon, maxlon
    ))

    bldgs <- tryCatch(.cache_get_buildings(minlat, maxlat, minlon, maxlon),
                      error = function(e) { message("Buildings cache error: ", e$message); NULL })
    if (!is.null(bldgs) && nrow(bldgs) > 0) {
      message("Buildings lookup: using SQLite cache (", nrow(bldgs), " polygons)")
      return(bldgs)
    }

    bldgs <- tryCatch(
      .get_preloaded_buildings_for_bbox(minlat, maxlat, minlon, maxlon),
      error = function(e) { message("Buildings preloaded fallback error: ", e$message); NULL }
    )
    if (!is.null(bldgs) && nrow(bldgs) > 0) {
      message("Buildings lookup: using preloaded CSV (", nrow(bldgs), " polygons)")
      return(bldgs)
    }

    if (requireNamespace("osmdata", quietly = TRUE)) {
      message("Buildings lookup: local sources empty; trying live OSM fetch/cache")
      bldgs <- tryCatch(
        get_buildings_from_cache(BUILDING_CACHE_DB,
                                 lat_min = minlat, lat_max = maxlat,
                                 lon_min = minlon, lon_max = maxlon,
                                 fetch_missing = TRUE),
        error = function(e) {
          message("Buildings live OSM fetch error: ", e$message)
          NULL
        }
      )
      if (!is.null(bldgs) && nrow(bldgs) > 0) {
        message("Buildings lookup: using live OSM/cache (", nrow(bldgs), " polygons)")
        return(bldgs)
      }
    } else {
      message("Buildings lookup: osmdata is not installed; live OSM fetch unavailable")
    }

    message("Buildings lookup: no building polygons found for dataset extent")
    NULL
  }

  # Helper: store loaded data and trigger building load
  .set_data <- function(dt, source_label) {
    rv$raw_data    <- dt
    rv$columns     <- names(dt)
    rv$data_source <- source_label
    rv$status_msg  <- paste0("Loaded ", format(nrow(dt), big.mark = ","),
                              " rows, ", ncol(dt), " columns [", source_label, "]")
    # Load buildings in background after data is set
    shinyjs::delay(100, {
      bldgs <- .load_buildings_for_data(dt)
      rv$buildings_sf <- bldgs
      if (!is.null(bldgs)) {
        showNotification(
          paste0(format(nrow(bldgs), big.mark=","), " building polygons loaded"),
          type = "message", duration = 4)
      } else {
        showNotification(
          "0 building polygons found for this dataset extent. Classification will use displacement only.",
          type = "warning", duration = 8)
      }
    })
  }

  # CSV
  observeEvent(input$traj_file, {
    req(input$traj_file)
    tryCatch({
      dt <- fread(input$traj_file$datapath, showProgress = FALSE, nrows = 100000)
      .set_data(dt, "csv")
    }, error = function(e) {
      showNotification(paste("Failed to load CSV:", e$message), type = "error")
    })
  })

  # RDS — raw trajectory data.frame / data.table
  observeEvent(input$rds_traj_file, {
    req(input$rds_traj_file)
    tryCatch({
      dt <- load_rds_trajectory(input$rds_traj_file$datapath)
      .set_data(dt, "rds")
      showNotification(paste0("RDS loaded: ", format(nrow(dt), big.mark = ","),
                               " rows, ", ncol(dt), " columns."), type = "message")
    }, error = function(e) {
      showNotification(paste("RDS load failed:", e$message), type = "error", duration = 8)
    })
  })

  # GeoJSON — spatial trajectory features (Point / LineString)
  observeEvent(input$geojson_file, {
    req(input$geojson_file)
    tryCatch({
      dt <- load_geojson_trajectory(input$geojson_file$datapath)
      .set_data(dt, "geojson")
      showNotification(paste0("GeoJSON loaded: ", format(nrow(dt), big.mark = ","),
                               " points."), type = "message")
    }, error = function(e) {
      showNotification(paste("GeoJSON load failed:", e$message), type = "error", duration = 8)
    })
  })

  # GPX / Wearable
  observeEvent(input$gpx_file, {
    req(input$gpx_file)
    tryCatch({
      dt <- load_gpx_trajectory(input$gpx_file$datapath,
                                 pid_label = tools::file_path_sans_ext(
                                   basename(input$gpx_file$name)))
      .set_data(dt, "gpx")
      showNotification(paste0("GPX loaded: ", format(nrow(dt), big.mark = ","),
                               " track points."), type = "message")
    }, error = function(e) {
      showNotification(paste("GPX load failed:", e$message), type = "error")
    })
  })

  # Remote URL
  observeEvent(input$load_url_btn, {
    req(nchar(trimws(input$data_url %||% "")) > 0)
    tryCatch({
      dt <- load_url_trajectory(trimws(input$data_url))
      .set_data(dt, "url")
      showNotification(paste0("URL loaded: ", format(nrow(dt), big.mark = ","),
                               " rows."), type = "message")
    }, error = function(e) {
      showNotification(paste("URL load failed:", e$message), type = "error")
    })
  })

  # Session restore — reloads a previous SPECTRA analysis result RDS
  observeEvent(input$rds_session_file, {
    req(input$rds_session_file)
    tryCatch({
      result <- readRDS(input$rds_session_file$datapath)
      if (!is.list(result) || !"points_classified" %in% names(result))
        stop("This file does not look like a SPECTRA session. Use the format buttons above for raw trajectory data.")
      rv$analysis_result <- result
      rv$buildings_sf    <- result$buildings_sf
      rv$status_msg      <- paste0("Session restored: ",
                                    format(nrow(result$points_classified), big.mark = ","),
                                    " classified points")
      showNotification(tagList(icon("check-circle"), " Session restored!"),
                       type = "message", duration = 4)
      go_tab("dashboard")
    }, error = function(e) {
      showNotification(paste("Session restore failed:", e$message), type = "error", duration = 8)
    })
  })

  output$file_loaded <- reactive({ !is.null(rv$raw_data) })
  outputOptions(output, "file_loaded", suspendWhenHidden = FALSE)
  
  output$file_info_ui <- renderUI({
    req(rv$raw_data)
    dt <- rv$raw_data
    size_mb <- file.info(input$traj_file$datapath)$size / 1024^2
    
    div(
      div(class = "alert alert-success",
          icon("check-circle"),
          paste0("File loaded: ", nrow(dt), " rows")),
      fluidRow(
        column(6, div(class = "metric-card",
                      div(class = "metric-value", paste0(nrow(dt), "+")),
                      div(class = "metric-label", "Rows"))),
        column(6, div(class = "metric-card",
                      div(class = "metric-value", sprintf("%.2f MB", size_mb)),
                      div(class = "metric-label", "Size")))
      )
    )
  })
  
  output$file_preview_text <- renderPrint({
    req(rv$raw_data)
    cat("First 5 rows:\n")
    print(head(rv$raw_data, 5))
  })
  
  # ====================================================================
  # Dynamic Field Mapping UI
  # ====================================================================
  
  # Helper to find best matching column
  find_match <- function(cols, patterns) {
    for (p in patterns) {
      m <- grep(p, cols, ignore.case = TRUE, value = TRUE)
      if (length(m) > 0) return(m[1])
    }
    return(cols[1])
  }
  
  output$lat_field_ui <- renderUI({
    req(rv$columns)
    sel <- find_match(rv$columns, c("^lat$", "lat", "latitude"))
    selectInput("lat_field", "Latitude", choices = rv$columns, selected = sel)
  })
  
  output$lon_field_ui <- renderUI({
    req(rv$columns)
    sel <- find_match(rv$columns, c("^lon$", "^lng$", "lon", "longitude"))
    selectInput("lon_field", "Longitude", choices = rv$columns, selected = sel)
  })
  
  output$pid_field_ui <- renderUI({
    req(rv$columns)
    sel <- find_match(rv$columns, c("^PID$", "^pid$", "person", "id", "participant"))
    selectInput("pid_field", "Person ID", choices = rv$columns, selected = sel)
  })
  
  output$date_field_ui <- renderUI({
    req(rv$columns)
    sel <- find_match(rv$columns, c("^date$", "datetime", "time", "timestamp"))
    selectInput("date_field", "Date/Time", choices = rv$columns, selected = sel)
  })
  
  output$pm25_field_ui <- renderUI({
    req(rv$columns)
    sel <- find_match(rv$columns, c("^pm25$", "pm2.5", "pm_25", "exposure_var1"))
    selectInput("pm25_field", "Exposure 1", choices = rv$columns, selected = sel)
  })
  
  output$no2_field_ui <- renderUI({
    req(rv$columns)
    sel <- find_match(rv$columns, c("^no2$", "nitrogen", "exposure_var2"))
    selectInput("no2_field", "Exposure 2", choices = rv$columns, selected = sel)
  })
  
  output$truth_field_ui <- renderUI({
    req(rv$columns)
    choices <- c("None", rv$columns)
    sel <- "None"
    for (p in c("truth", "ground", "actual", "label", "indoor_truth")) {
      m <- grep(p, rv$columns, ignore.case = TRUE, value = TRUE)
      if (length(m) > 0) { sel <- m[1]; break }
    }
    selectInput("truth_field", "Ground-truth", choices = choices, selected = sel)
  })
  
  # ====================================================================
  # Live Map
  # ====================================================================
  
  output$live_map <- renderLeaflet({
    leaflet() |>
      addProviderTiles("CartoDB.Positron") |>
      setView(lng = 35, lat = -15, zoom = 5)
  })
  
  # Reactive: sampled GPS points (only recomputes when data or field selection changes)
  live_map_pts <- reactive({
    req(rv$raw_data, input$lat_field, input$lon_field)
    dt <- rv$raw_data
    if (!(input$lat_field %in% names(dt)) || !(input$lon_field %in% names(dt))) return(NULL)
    lat <- as.numeric(dt[[input$lat_field]])
    lon <- as.numeric(dt[[input$lon_field]])
    valid <- !is.na(lat) & !is.na(lon)
    if (sum(valid) == 0) return(NULL)
    idx <- which(valid)
    if (length(idx) > 3000) idx <- sample(idx, 3000)
    list(lat = lat[idx], lon = lon[idx])
  })

  # Draw points and optional building footprints whenever data changes
  observeEvent(list(live_map_pts(), rv$buildings_sf, input$show_buildings_live), {
    pts <- live_map_pts()
    if (is.null(pts)) return()
    tile_provider <- isolate(input$live_map_style) %||% "CartoDB.Positron"
    m <- leafletProxy("live_map") |>
      clearMarkers() |>
      clearShapes() |>
      clearTiles() |>
      addProviderTiles(tile_provider) |>
      addCircleMarkers(
        lng = pts$lon, lat = pts$lat,
        radius = 3, color = "#FF6B6B",
        fillOpacity = 0.7, stroke = FALSE,
        popup = paste0("Lat: ", round(pts$lat, 5), "<br>Lon: ", round(pts$lon, 5)),
        group = "GPS Points"
      ) |>
      fitBounds(
        lng1 = min(pts$lon), lat1 = min(pts$lat),
        lng2 = max(pts$lon), lat2 = max(pts$lat)
      )

    if (isTRUE(input$show_buildings_live)) {
      m <- .add_building_polygons(m, rv$buildings_sf)
    }
  }, ignoreNULL = FALSE)

  # Swap tile layer only (no data reprocessing)
  observeEvent(input$live_map_style, {
    req(input$live_map_style)
    leafletProxy("live_map") |>
      clearTiles() |>
      addProviderTiles(input$live_map_style)
  }, ignoreInit = TRUE)

  output$live_map_info_ui <- renderUI({
    req(rv$raw_data, input$lat_field, input$lon_field)
    dt <- rv$raw_data
    lat <- as.numeric(dt[[input$lat_field]])
    lon <- as.numeric(dt[[input$lon_field]])
    valid <- sum(!is.na(lat) & !is.na(lon))
    
    div(
      h4("Data Summary"),
      p(icon("map-marker-alt"), paste0(valid, " valid GPS points")),
      p(icon("arrows-alt"), sprintf("Lat: %.4f - %.4f", min(lat, na.rm = TRUE), max(lat, na.rm = TRUE))),
      p(icon("arrows-alt-h"), sprintf("Lon: %.4f - %.4f", min(lon, na.rm = TRUE), max(lon, na.rm = TRUE)))
    )
  })
  
  # ====================================================================
  # Run Analysis
  # ====================================================================
  
  run_the_analysis <- function() {
    # Only require the core navigation fields; exposure fields are optional
    if (is.null(rv$raw_data)) {
      showNotification("Please upload data first.", type = "warning", duration = 4)
      return()
    }
    if (is.null(input$lat_field) || is.null(input$lon_field) ||
        is.null(input$pid_field) || is.null(input$date_field)) {
      showNotification("Please map Lat, Lon, Person ID, and Date fields before running.", type = "warning", duration = 5)
      return()
    }

    showNotification(
      tagList(icon("spinner", class = "fa-spin"), " Analysis running — please wait..."),
      type = "message", duration = NULL, id = "analysis_running_notif"
    )

    withProgress(message = "SPECTRA — Analysis Progress", value = 0, {

      tryCatch({
        # ── Step 1 (10%): Load & standardise data ─────────────────────────────
        incProgress(0.10, detail = "10% — Reading and standardising data...")

        src <- rv$data_source %||% "csv"

        if (src == "csv") {
          req(input$traj_file)
          dt <- load_csv_trajectory(
            input$traj_file$datapath,
            lat_field  = input$lat_field,
            lon_field  = input$lon_field,
            pid_field  = input$pid_field,
            date_field = input$date_field,
            pm25_field = input$pm25_field,
            no2_field  = input$no2_field,
            csv_epsg   = input$csv_epsg
          )
        } else {
          # GPX / URL — already standardised in rv$raw_data; re-apply field mapping
          req(rv$raw_data)
          dt <- copy(rv$raw_data)
          if (!("lat"  %in% names(dt))) setnames(dt, input$lat_field,  "lat",  skip_absent = TRUE)
          if (!("lon"  %in% names(dt))) setnames(dt, input$lon_field,  "lon",  skip_absent = TRUE)
          if (!("PID"  %in% names(dt))) setnames(dt, input$pid_field,  "PID",  skip_absent = TRUE)
          if (!("date" %in% names(dt))) setnames(dt, input$date_field, "date", skip_absent = TRUE)
          if (!("pm25" %in% names(dt)) && !is.null(input$pm25_field) &&
              input$pm25_field != "None" && input$pm25_field %in% names(dt))
            setnames(dt, input$pm25_field, "pm25", skip_absent = TRUE)
          if (!("no2"  %in% names(dt)) && !is.null(input$no2_field) &&
              input$no2_field != "None" && input$no2_field %in% names(dt))
            setnames(dt, input$no2_field, "no2", skip_absent = TRUE)
          if (!"pm25" %in% names(dt)) dt[, pm25 := NA_real_]
          if (!"no2"  %in% names(dt)) dt[, no2  := NA_real_]
          dt[, lat  := as.numeric(lat)]
          dt[, lon  := as.numeric(lon)]
          dt[, PID  := as.character(PID)]
          dt[, pm25 := as.numeric(pm25)]
          dt[, no2  := as.numeric(no2)]
          dt <- dt[!is.na(lat) & !is.na(lon) & !is.na(PID) & !is.na(date)]
        }

        # Optional standard columns
        if (!is.null(input$country_field) && input$country_field != "None" &&
            input$country_field %in% names(dt))
          setnames(dt, input$country_field, "country", skip_absent = TRUE)
        if (!is.null(input$temp_field) && input$temp_field != "None" &&
            input$temp_field %in% names(dt)) {
          setnames(dt, input$temp_field, "temperature_c", skip_absent = TRUE)
          dt[, temperature_c := as.numeric(temperature_c)]
        }
        if (!is.null(input$humid_field) && input$humid_field != "None" &&
            input$humid_field %in% names(dt)) {
          setnames(dt, input$humid_field, "humidity_pct", skip_absent = TRUE)
          dt[, humidity_pct := as.numeric(humidity_pct)]
        }

        # Wearable / biometric columns
        if (!is.null(input$heart_rate_field) && input$heart_rate_field != "None" &&
            input$heart_rate_field %in% names(dt)) {
          setnames(dt, input$heart_rate_field, "heart_rate_bpm", skip_absent = TRUE)
          dt[, heart_rate_bpm := as.numeric(heart_rate_bpm)]
        }
        if (!is.null(input$steps_field) && input$steps_field != "None" &&
            input$steps_field %in% names(dt)) {
          setnames(dt, input$steps_field, "steps", skip_absent = TRUE)
          dt[, steps := as.numeric(steps)]
        }
        if (!is.null(input$activity_field) && input$activity_field != "None" &&
            input$activity_field %in% names(dt)) {
          setnames(dt, input$activity_field, "activity_type", skip_absent = TRUE)
          dt[, activity_type := as.character(activity_type)]
        }

        # Store exposure variable labels for use in outputs
        rv$exp1_label <- trimws(input$exp1_label %||% "Exposure 1")
        rv$exp2_label <- trimws(input$exp2_label %||% "Exposure 2")

        # ── Step 2 (20%): Load boundary shapefile ─────────────────────────────
        incProgress(0.10, detail = "20% — Loading boundary layer...")
        boundary_sf <- NULL
        if (!is.null(input$boundary_file)) {
          boundary_sf <- load_boundary_layer(input$boundary_file$datapath)
        }

        # ── Step 3 (35%): Buildings from RAM cache ────────────────────────────
        incProgress(0.15, detail = "35% — Fetching building footprints...")
        buildings_sf <- if (!is.null(rv$buildings_sf)) rv$buildings_sf else .load_buildings_for_data(dt)
        if (!is.null(buildings_sf) && is.null(rv$buildings_sf)) rv$buildings_sf <- buildings_sf
        n_bldg <- if (!is.null(buildings_sf)) nrow(buildings_sf) else 0
        fetch_osm_flag <- FALSE
        osm_buffer     <- 0.5

        # ── Step 4 (45%): Ground truth mapping ───────────────────────────────
        incProgress(0.10, detail = "45% — Mapping ground truth labels...")
        truth_field_name <- NULL
        if (!is.null(input$truth_field) && input$truth_field != "None") {
          truth_field_name <- input$truth_field
        }

        # ── Step 5 (50–90%): Core trajectory analysis ────────────────────────
        is_cloud <- grepl("^/srv/connect", getwd())

        incProgress(0.05, detail = paste0("50% — Running analysis on ",
                                          format(nrow(dt), big.mark = ","),
                                          " points & ", format(n_bldg, big.mark = ","),
                                          " buildings..."))

        analysis_args <- list(
          boundary_sf              = boundary_sf,
          buildings_sf             = buildings_sf,
          fetch_osm                = fetch_osm_flag,
          osm_buffer_km            = osm_buffer,
          cache_db_path            = BUILDING_CACHE_DB,
          cell_size_m              = input$cell_size,
          max_dwell_min            = input$max_dwell,
          speed_threshold_kmh      = input$speed_limit,
          remove_teleports         = input$remove_teleports,
          building_buffer_m        = input$building_buffer,
          exclude_indoor           = input$exclude_indoor,
          speed_outdoor_threshold  = input$speed_outdoor_threshold,
          use_building_validation  = isTRUE(input$use_building_validation)
        )

        result <- if (is_cloud) {
          do.call(run_analysis_batched, c(list(dt = dt, batch_pts = 100000L), analysis_args))
        } else {
          do.call(run_analysis, c(list(dt), analysis_args))
        }

        # ── Step 6 (92%): Store results ───────────────────────────────────────
        incProgress(0.42, detail = "92% — Storing results...")

        if (!is.null(truth_field_name) && truth_field_name %in% names(result$points_classified)) {
          indoor_labels <- unlist(strsplit(input$indoor_labels, ","))
          vals <- tolower(trimws(as.character(result$points_classified[[truth_field_name]])))
          result$points_classified[, truth_is_indoor := vals %in% tolower(trimws(indoor_labels))]
          rv$truth_field <- "truth_is_indoor"
        } else {
          rv$truth_field <- NULL
        }

        rv$buildings_sf    <- result$buildings_sf %||% rv$buildings_sf
        result$buildings_sf <- NULL   # already saved in rv$buildings_sf; free the duplicate
        rv$analysis_result <- result
        rv$status_msg      <- "Analysis completed successfully!"
        gc(verbose = FALSE)

        # ── Step 7 (100%): Done ───────────────────────────────────────────────
        incProgress(0.08, detail = "100% — Complete!")
        removeNotification("analysis_running_notif")
        showNotification(
          tagList(icon("check-circle"), " Analysis complete!"),
          type = "message", duration = 6
        )
        go_tab("dashboard")

      }, error = function(e) {
        removeNotification("analysis_running_notif")
        showNotification(
          tagList(icon("exclamation-triangle"), " Analysis failed: ", e$message),
          type = "error", duration = 12
        )
        rv$status_msg <- paste("Error:", e$message)
      })
    })
  }
  
  observeEvent(input$run_analysis, { run_the_analysis() })
  observeEvent(input$run_analysis2, { run_the_analysis() })
  
  output$run_status_ui <- renderUI({
    if (rv$status_msg == "") return(NULL)
    
    if (grepl("Error", rv$status_msg)) {
      div(class = "alert alert-danger", icon("exclamation-circle"), rv$status_msg)
    } else if (grepl("success", rv$status_msg, ignore.case = TRUE)) {
      div(class = "alert alert-success", icon("check-circle"), rv$status_msg)
    } else {
      div(class = "alert alert-info", icon("info-circle"), rv$status_msg)
    }
  })
  
  output$run_progress_ui <- renderUI({
    output$run_status_ui()
  })
  
  # ====================================================================
  # Helper: Filtered Data
  # ====================================================================
  
  filtered_data <- reactive({
    req(rv$analysis_result)
    dt <- as.data.table(rv$analysis_result$points_classified)
    
    if (!is.null(input$filter_pids) && length(input$filter_pids) > 0) {
      dt <- dt[PID %in% input$filter_pids]
    }
    
    if (!is.null(input$filter_dates)) {
      dates <- input$filter_dates
      if (length(dates) == 2) {
        dt <- dt[as.Date(date) >= dates[1] & as.Date(date) <= dates[2]]
      }
    }
    
    dt
  })
  
  # ====================================================================
  # Dashboard: Info Boxes
  # ====================================================================
  
  output$metric_classified <- renderInfoBox({
    stats <- rv$analysis_result$stats %||% list(n_points_classified = 0)
    infoBox("Classified", format(stats$n_points_classified, big.mark = ","),
            icon = icon("map-marker-alt"), color = "blue")
  })
  
  output$metric_analysis <- renderInfoBox({
    stats <- rv$analysis_result$stats %||% list(n_points_for_analysis = 0)
    infoBox("Analysis Pts", format(stats$n_points_for_analysis, big.mark = ","),
            icon = icon("chart-line"), color = "aqua")
  })
  
  output$metric_people <- renderInfoBox({
    stats <- rv$analysis_result$stats %||% list(n_people = 0)
    infoBox("People", stats$n_people,
            icon = icon("users"), color = "green")
  })
  
  output$metric_indoor <- renderInfoBox({
    stats <- rv$analysis_result$stats %||% list(n_indoor_pred = 0)
    infoBox("Indoor Predicted", format(stats$n_indoor_pred, big.mark = ","),
            icon = icon("home"), color = "purple")
  })
  
  output$metric_dwell <- renderInfoBox({
    stats <- rv$analysis_result$stats %||% list(total_dwell_min = 0)
    infoBox("Total Dwell (min)", format(round(stats$total_dwell_min), big.mark = ","),
            icon = icon("clock"), color = "orange")
  })
  
  # ====================================================================
  # Dashboard: Filters
  # ====================================================================
  
  output$filter_pid_ui <- renderUI({
    req(rv$analysis_result)
    pids <- sort(unique(rv$analysis_result$points_classified$PID))
    pickerInput("filter_pids", "Select Participants",
                choices = pids, selected = pids,
                options = list(`actions-box` = TRUE, `live-search` = TRUE,
                               `selected-text-format` = "count > 3"),
                multiple = TRUE)
  })
  
  output$filter_date_ui <- renderUI({
    req(rv$analysis_result)
    dates <- as.Date(rv$analysis_result$points_classified$date)
    dateRangeInput("filter_dates", "Date Range",
                   start = min(dates, na.rm = TRUE),
                   end = max(dates, na.rm = TRUE),
                   min = min(dates, na.rm = TRUE),
                   max = max(dates, na.rm = TRUE))
  })
  
  output$filtered_summary_ui <- renderUI({
    dt <- filtered_data()
    if (is.null(dt) || nrow(dt) == 0) return(p("No data matches filters"))
    
    indoor_n <- sum(dt$is_indoor, na.rm = TRUE)
    outdoor_n <- sum(!dt$is_indoor, na.rm = TRUE)
    total <- nrow(dt)
    pct <- if (total > 0) round(indoor_n / total * 100, 1) else 0
    
    div(
      fluidRow(
        column(6, div(class = "metric-card",
                      div(class = "metric-value", format(indoor_n, big.mark = ",")),
                      div(class = "metric-label", paste0("Indoor (", pct, "%)")))),
        column(6, div(class = "metric-card",
                      div(class = "metric-value", format(outdoor_n, big.mark = ",")),
                      div(class = "metric-label", paste0("Outdoor (", round(100 - pct, 1), "%)"))))
      )
    )
  })
  
  # ====================================================================
  # Dashboard: Charts
  # ====================================================================
  
  output$plot_dwell_hist <- renderPlotly({
    dt <- filtered_data()
    req(nrow(dt) > 0)
    
    dt[, label := ifelse(is_indoor, "Indoor", "Outdoor")]
    
    p <- plot_ly(dt, x = ~dwell_min, color = ~label,
                 colors = c("Indoor" = "#2a5298", "Outdoor" = "#4a90e2"),
                 type = "histogram", nbinsx = 40) |>
      layout(title = "Dwell Time Distribution",
             xaxis = list(title = "Dwell Time (minutes)"),
             yaxis = list(title = "Frequency"),
             barmode = "overlay",
             plot_bgcolor = "rgba(0,0,0,0)",
             paper_bgcolor = "rgba(0,0,0,0)")
    p
  })
  
  output$plot_indoor_pie <- renderPlotly({
    dt <- filtered_data()
    req(nrow(dt) > 0)
    
    counts <- data.frame(
      Category = c("Indoor", "Outdoor"),
      Count = c(sum(dt$is_indoor), sum(!dt$is_indoor))
    )
    
    plot_ly(counts, labels = ~Category, values = ~Count, type = "pie",
            marker = list(colors = c("#2a5298", "#4a90e2")),
            hole = 0.4) |>
      layout(title = "Indoor/Outdoor Distribution",
             plot_bgcolor = "rgba(0,0,0,0)",
             paper_bgcolor = "rgba(0,0,0,0)")
  })
  
  output$dwell_per_person_table <- DT::renderDataTable({
    dt <- filtered_data()
    req(nrow(dt) > 0)
    
    dwell_pp <- dt[, .(
      Total_Dwell_Min = sum(dwell_min, na.rm = TRUE),
      Indoor_Dwell_Min = sum(dwell_min[is_indoor == TRUE], na.rm = TRUE)
    ), by = PID]
    
    dwell_pp[, Outdoor_Dwell_Min := Total_Dwell_Min - Indoor_Dwell_Min]
    dwell_pp[, Total_Dwell_Hours := round(Total_Dwell_Min / 60, 2)]
    dwell_pp[, Indoor_Pct := round(Indoor_Dwell_Min / Total_Dwell_Min * 100, 1)]
    
    DT::datatable(dwell_pp, options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE) |>
      DT::formatRound(c("Total_Dwell_Min", "Indoor_Dwell_Min", "Outdoor_Dwell_Min"), 1)
  })
  
  output$plot_dwell_bar <- renderPlotly({
    dt <- filtered_data()
    req(nrow(dt) > 0)
    
    dwell_pp <- dt[, .(
      Total_Hours = sum(dwell_min, na.rm = TRUE) / 60,
      Indoor_Pct = sum(dwell_min[is_indoor == TRUE], na.rm = TRUE) / sum(dwell_min, na.rm = TRUE) * 100
    ), by = PID]
    
    top20 <- dwell_pp[order(-Total_Hours)][1:min(20, nrow(dwell_pp))]
    
    plot_ly(top20, x = ~reorder(PID, -Total_Hours), y = ~Total_Hours,
            color = ~Indoor_Pct, type = "bar",
            colors = viridis::viridis(10)) |>
      layout(title = "Top 20 Participants by Dwell Time",
             xaxis = list(title = "Participant ID"),
             yaxis = list(title = "Total Dwell Time (hours)"),
             plot_bgcolor = "rgba(0,0,0,0)",
             paper_bgcolor = "rgba(0,0,0,0)")
  })
  
  # ====================================================================
  # Spatial View
  # ====================================================================
  
  output$spatial_map <- renderLeaflet({
    dt <- filtered_data()
    req(nrow(dt) > 0)
    
    # Sample for performance
    n <- min(nrow(dt), input$spatial_sample %||% 5000)
    if (nrow(dt) > n) {
      dt <- dt[sample(.N, n)]
    }
    
    # Colour palette
    color_by <- input$spatial_color_by %||% "Classification"
    
    if (color_by == "Classification") {
      dt[, color := ifelse(is_indoor, "#2a5298", "#4a90e2")]
      dt[, popup_label := ifelse(is_indoor, "Indoor", "Outdoor")]
    } else if (grepl("Exposure 1", color_by)) {
      pal <- colorNumeric("YlOrRd", dt$pm25, na.color = "#ccc")
      dt[, color := pal(pm25)]
      dt[, popup_label := paste0(rv$exp1_label %||% "Exposure 1", ": ", round(pm25, 2))]
    } else if (grepl("Exposure 2", color_by)) {
      pal <- colorNumeric("YlGnBu", dt$no2, na.color = "#ccc")
      dt[, color := pal(no2)]
      dt[, popup_label := paste0(rv$exp2_label %||% "Exposure 2", ": ", round(no2, 2))]
    } else {
      pal <- colorNumeric("RdYlGn", dt$speed_kmh, na.color = "#ccc", reverse = TRUE)
      dt[, color := pal(speed_kmh)]
      dt[, popup_label := paste0("Speed: ", round(speed_kmh, 1), " km/h")]
    }
    
    m <- leaflet() |>
      addProviderTiles(input$spatial_map_style %||% "CartoDB.Positron")

    overlay_groups <- character(0)

    if (isTRUE(input$show_buildings_spatial)) {
      m <- .add_building_polygons(m, rv$buildings_sf)
      overlay_groups <- c(overlay_groups, "Building Polygons")
    }

    m <- m |>
      addCircleMarkers(
        data = dt,
        lng = ~lon, lat = ~lat,
        radius = 3, color = ~color,
        fillOpacity = 0.7, stroke = FALSE,
        group = "Points",
        popup = ~paste0(
          "<b>PID:</b> ", PID, "<br>",
          "<b>Date:</b> ", date, "<br>",
          "<b>", popup_label, "</b><br>",
          "<b>Dwell:</b> ", round(dwell_min, 1), " min"
        )
      )

    # Density heatmap overlay
    if (isTRUE(input$show_heatmap)) {
      pal_heat <- colorNumeric("YlOrRd", domain = dt$dwell_min, na.color = "transparent")
      m <- m |>
        addCircleMarkers(
          lng = ~lon, lat = ~lat,
          radius = 12, color = ~pal_heat(dwell_min),
          fillOpacity = 0.3, stroke = FALSE,
          group = "Density"
        )
      overlay_groups <- c(overlay_groups, "Density")
    }

    if (length(overlay_groups) > 0) {
      m <- m |>
        addLayersControl(overlayGroups = overlay_groups,
                         options = layersControlOptions(collapsed = FALSE))
    }

    m
  })

  # ====================================================================
  # Raster Grid Map (Dwell Time & Exposure)
  # ====================================================================

  output$raster_date_ui <- renderUI({
    req(rv$analysis_result)
    dates <- sort(unique(rv$analysis_result$person_day_cell$date_only))
    selectInput("raster_date", "Select Date:", choices = c("All Days" = "all", as.character(dates)))
  })

  raster_grid_data <- reactive({
    req(rv$analysis_result)
    pdc <- as.data.table(rv$analysis_result$person_day_cell)
    if (!is.null(input$raster_date) && input$raster_date != "all") {
      pdc <- pdc[date_only == as.Date(input$raster_date)]
    }
    # Aggregate across PIDs per cell
    grid <- pdc[, .(
      dwell = sum(dwell_sum_min, na.rm = TRUE),
      pm25 = weighted.mean(pm25_mean_tw, dwell_sum_min, na.rm = TRUE),
      no2 = weighted.mean(no2_mean_tw, dwell_sum_min, na.rm = TRUE),
      indoor_pct = mean(pct_indoor, na.rm = TRUE),
      lon = mean(lon, na.rm = TRUE),
      lat = mean(lat, na.rm = TRUE),
      row_idx = first(row_idx),
      col_idx = first(col_idx)
    ), by = cell_id]
    grid[is.nan(pm25), pm25 := 0]
    grid[is.nan(no2), no2 := 0]
    grid
  })

  output$raster_grid_map <- renderLeaflet({
    grid <- raster_grid_data()
    req(nrow(grid) > 0)

    var_col <- input$raster_variable %||% "dwell"
    vals <- grid[[var_col]]
    var_label <- switch(var_col,
      dwell = "Dwell (min)", pm25 = "Exposure 1", no2 = "Exposure 2", indoor_pct = "Indoor %")

    pal <- colorNumeric("YlOrRd", domain = vals, na.color = "#ccc")
    cell_size_deg <- (rv$analysis_result$grid_params$cell_size %||% 100) / 111000

    half <- cell_size_deg / 2
    cell_colors  <- pal(grid[[var_col]])
    cell_popups  <- paste0(
      "<b>Cell:</b> ", grid$cell_id, "<br>",
      "<b>Dwell:</b> ", round(grid$dwell, 1), " min<br>",
      "<b>Exposure 1:</b> ", round(grid$pm25, 2), "<br>",
      "<b>Exposure 2:</b> ", round(grid$no2, 2), "<br>",
      "<b>Indoor %:</b> ", round(grid$indoor_pct, 1)
    )

    leaflet(grid) |>
      addProviderTiles("CartoDB.Positron") |>
      addRectangles(
        lng1 = grid$lon - half, lat1 = grid$lat - half,
        lng2 = grid$lon + half, lat2 = grid$lat + half,
        color = cell_colors, weight = 0.5,
        fillColor = cell_colors, fillOpacity = 0.6,
        popup = cell_popups
      ) |>
      addLegend(pal = pal, values = vals, title = var_label, position = "bottomright") |>
      fitBounds(
        lng1 = min(grid$lon), lat1 = min(grid$lat),
        lng2 = max(grid$lon), lat2 = max(grid$lat)
      )
  })

  # NetCDF Export — full continuous grid (labeled + unlabeled cells)
  output$download_raster_nc <- downloadHandler(
    filename = function() paste0("exposure_grid_", Sys.Date(), ".nc"),
    content = function(file) {
      req(rv$analysis_result, rv$analysis_result$daily_grid, rv$analysis_result$full_grid)
      dg <- as.data.table(rv$analysis_result$daily_grid)   # date × cell
      fg <- as.data.table(rv$analysis_result$full_grid)    # all cells (labeled + unlabeled)
      gp <- rv$analysis_result$grid_params

      mv  <- -999.0
      mvi <- -1L

      # ── Spatial axes from full_grid (covers all cells incl. unlabeled) ───────
      lon_vals <- sort(unique(round(fg$lon, 6)))
      lat_vals <- sort(unique(round(fg$lat, 6)))
      nlon     <- length(lon_vals)
      nlat     <- length(lat_vals)

      # ── Time axis — CF-compliant "days since" reference ──────────────────────
      dates     <- sort(unique(dg$date))
      ref_date  <- as.Date("1970-01-01")
      time_vals <- as.integer(dates - ref_date)   # days since 1970-01-01
      ntime     <- length(time_vals)

      # ── Compression probe ────────────────────────────────────────────────────
      use_compression <- tryCatch({
        ncdf4::nc_create(tempfile(),
          list(ncdf4::ncvar_def("t", "", ncdf4::ncdim_def("x", "", 1), compression = 6)),
          force_v4 = TRUE)
        TRUE
      }, error = function(e) FALSE)
      cmp <- if (use_compression) 6L else NA_integer_

      # ── Dimensions ───────────────────────────────────────────────────────────
      lon_dim  <- ncdf4::ncdim_def("lon",  "degrees_east",  lon_vals,
                                   longname = "longitude of grid cell centre")
      lat_dim  <- ncdf4::ncdim_def("lat",  "degrees_north", lat_vals,
                                   longname = "latitude of grid cell centre")
      time_dim <- ncdf4::ncdim_def("time", paste0("days since ", format(ref_date)),
                                   time_vals, unlim = TRUE,
                                   longname = "date of observation")

      # ── 2-D static variables (grid_id and is_labeled — one value per cell) ───
      var_grid_id <- ncdf4::ncvar_def("grid_id", "1", list(lon_dim, lat_dim),
                                      missval = mvi, prec = "integer", compression = cmp,
                                      longname = "Unique grid cell identifier (cell_id)")
      var_labeled <- ncdf4::ncvar_def("is_labeled", "1", list(lon_dim, lat_dim),
                                      missval = mvi, prec = "short", compression = cmp,
                                      longname = "1=cell has GPS observations across any date, 0=never observed")
      var_crs     <- ncdf4::ncvar_def("crs", "1", list(),
                                      prec = "integer",
                                      longname = "WGS 84 coordinate reference system")

      # ── 3-D time-resolved variables (lon × lat × time) ───────────────────────
      var_dwell  <- ncdf4::ncvar_def("dwell_time", "minutes",
                                     list(lon_dim, lat_dim, time_dim), missval = mv,
                                     compression = cmp,
                                     longname = "Total dwell time in cell on this date (minutes)")
      var_pm25   <- ncdf4::ncvar_def("pm25", "",
                                     list(lon_dim, lat_dim, time_dim), missval = mv,
                                     compression = cmp,
                                     longname = paste0("Time-weighted mean ",
                                                       rv$exp1_label %||% "PM2.5", " on this date"))
      var_no2    <- ncdf4::ncvar_def("no2", "",
                                     list(lon_dim, lat_dim, time_dim), missval = mv,
                                     compression = cmp,
                                     longname = paste0("Time-weighted mean ",
                                                       rv$exp2_label %||% "NO2", " on this date"))
      var_indoor <- ncdf4::ncvar_def("indoor_pct", "percent",
                                     list(lon_dim, lat_dim, time_dim), missval = mv,
                                     compression = cmp,
                                     longname = "Percentage of dwell classified as indoor on this date")
      # time_start / time_end: seconds since 1970-01-01 (Unix epoch) for each
      # cell's first and last GPS observation on that date — enables sub-daily GEE matching
      var_tstart <- ncdf4::ncvar_def("time_start", "seconds since 1970-01-01T00:00:00Z",
                                     list(lon_dim, lat_dim, time_dim), missval = mv,
                                     compression = cmp,
                                     longname = "Unix timestamp of first GPS observation in cell on this date")
      var_tend   <- ncdf4::ncvar_def("time_end",   "seconds since 1970-01-01T00:00:00Z",
                                     list(lon_dim, lat_dim, time_dim), missval = mv,
                                     compression = cmp,
                                     longname = "Unix timestamp of last GPS observation in cell on this date")
      # Human-readable date string (YYYY-MM-DD) — 1-D coordinate on time dimension
      var_date   <- ncdf4::ncvar_def("date", "",
                                     list(time_dim), prec = "char",
                                     longname = "Calendar date of observation (YYYY-MM-DD)")

      nc <- ncdf4::nc_create(file,
                             list(var_grid_id, var_labeled, var_crs,
                                  var_dwell, var_pm25, var_no2, var_indoor,
                                  var_tstart, var_tend, var_date),
                             force_v4 = use_compression)

      # ── CF calendar attribute on time ─────────────────────────────────────────
      ncdf4::ncatt_put(nc, "time", "calendar",      "standard")
      ncdf4::ncatt_put(nc, "time", "axis",          "T")

      # Write human-readable date strings (one per time step)
      ncdf4::ncvar_put(nc, var_date, format(dates, "%Y-%m-%d"))
      ncdf4::ncatt_put(nc, "lon",  "axis",          "X")
      ncdf4::ncatt_put(nc, "lon",  "standard_name", "longitude")
      ncdf4::ncatt_put(nc, "lat",  "axis",          "Y")
      ncdf4::ncatt_put(nc, "lat",  "standard_name", "latitude")

      # ── WGS84 georeferencing / CF grid mapping ───────────────────────────────
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

      # ── Lookup vectors for full_grid cells ───────────────────────────────────
      fg_lon_idx <- match(round(fg$lon, 6), lon_vals)
      fg_lat_idx <- match(round(fg$lat, 6), lat_vals)
      fg_valid   <- which(!is.na(fg_lon_idx) & !is.na(fg_lat_idx))
      fg_idx_by_cell <- data.table::data.table(
        cell_id = fg$cell_id,
        lon_idx = fg_lon_idx,
        lat_idx = fg_lat_idx
      )

      # ── 2-D: grid_id ─────────────────────────────────────────────────────────
      gid_mat <- matrix(mvi, nrow = nlon, ncol = nlat)
      gid_mat[fg_lon_idx[fg_valid] + (fg_lat_idx[fg_valid] - 1L) * nlon] <- fg$cell_id[fg_valid]
      ncdf4::ncvar_put(nc, var_grid_id, gid_mat)
      rm(gid_mat)

      # ── 2-D: is_labeled (any-date aggregation from full_grid) ────────────────
      labeled_mat <- matrix(0L, nrow = nlon, ncol = nlat)
      obs_cells   <- fg_valid[fg$is_labeled[fg_valid] == 1L]
      labeled_mat[fg_lon_idx[obs_cells] + (fg_lat_idx[obs_cells] - 1L) * nlon] <- 1L
      ncdf4::ncvar_put(nc, var_labeled, labeled_mat)
      rm(labeled_mat)

      # ── 3-D: fill one time slice per date ────────────────────────────────────
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
        rows     <- which(dg$date == dates[ti] &
                          !is.na(dg_lon_idx) & !is.na(dg_lat_idx))
        li       <- dg_lon_idx[rows]
        ai       <- dg_lat_idx[rows]
        mat_idx  <- li + (ai - 1L) * nlon

        dwell_sl  <- matrix(mv, nrow = nlon, ncol = nlat)
        pm25_sl   <- matrix(mv, nrow = nlon, ncol = nlat)
        no2_sl    <- matrix(mv, nrow = nlon, ncol = nlat)
        ip_sl     <- matrix(mv, nrow = nlon, ncol = nlat)
        tstart_sl <- matrix(mv, nrow = nlon, ncol = nlat)
        tend_sl   <- matrix(mv, nrow = nlon, ncol = nlat)

        if (length(rows)) {
          dwell_sl[mat_idx] <- dg$dwell[rows]
          ok_pm  <- rows[!is.na(dg$pm25[rows])];   if (length(ok_pm))
            pm25_sl[dg_lon_idx[ok_pm]  + (dg_lat_idx[ok_pm]  - 1L)*nlon] <- dg$pm25[ok_pm]
          ok_no2 <- rows[!is.na(dg$no2[rows])];    if (length(ok_no2))
            no2_sl[dg_lon_idx[ok_no2]  + (dg_lat_idx[ok_no2] - 1L)*nlon] <- dg$no2[ok_no2]
          ok_ip  <- rows[!is.na(dg$indoor_pct[rows])]; if (length(ok_ip))
            ip_sl[dg_lon_idx[ok_ip]    + (dg_lat_idx[ok_ip]  - 1L)*nlon] <- dg$indoor_pct[ok_ip]
          # datetime_start / datetime_end stored as Unix seconds (numeric)
          if ("datetime_start" %in% names(dg)) {
            ts_vals <- suppressWarnings(
              as.numeric(as.POSIXct(dg$datetime_start[rows], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
            te_vals <- suppressWarnings(
              as.numeric(as.POSIXct(dg$datetime_end[rows],   format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
            ok_ts <- rows[!is.na(ts_vals)]; if (length(ok_ts))
              tstart_sl[dg_lon_idx[ok_ts] + (dg_lat_idx[ok_ts] - 1L)*nlon] <- ts_vals[!is.na(ts_vals)]
            ok_te <- rows[!is.na(te_vals)]; if (length(ok_te))
              tend_sl[dg_lon_idx[ok_te]   + (dg_lat_idx[ok_te] - 1L)*nlon] <- te_vals[!is.na(te_vals)]
          }
        }

        ncdf4::ncvar_put(nc, var_dwell,  dwell_sl,  start = c(1,1,ti), count = c(nlon,nlat,1))
        ncdf4::ncvar_put(nc, var_pm25,   pm25_sl,   start = c(1,1,ti), count = c(nlon,nlat,1))
        ncdf4::ncvar_put(nc, var_no2,    no2_sl,    start = c(1,1,ti), count = c(nlon,nlat,1))
        ncdf4::ncvar_put(nc, var_indoor, ip_sl,     start = c(1,1,ti), count = c(nlon,nlat,1))
        ncdf4::ncvar_put(nc, var_tstart, tstart_sl, start = c(1,1,ti), count = c(nlon,nlat,1))
        ncdf4::ncvar_put(nc, var_tend,   tend_sl,   start = c(1,1,ti), count = c(nlon,nlat,1))
        rm(dwell_sl, pm25_sl, no2_sl, ip_sl, tstart_sl, tend_sl)
      }

      # ── CF-1.8 global attributes ──────────────────────────────────────────────
      ncdf4::ncatt_put(nc, 0, "Conventions",       "CF-1.8")
      ncdf4::ncatt_put(nc, 0, "title",             "SPECTRA Exposure Grid")
      ncdf4::ncatt_put(nc, 0, "institution",       "CHEAQI-MNCH")
      ncdf4::ncatt_put(nc, 0, "source",            "SPECTRA — Spatiotemporal Personal Exposure Characterization and TRAjectory Analyzer")
      ncdf4::ncatt_put(nc, 0, "history",           paste("Created", as.character(Sys.time())))
      ncdf4::ncatt_put(nc, 0, "grid_cell_size_m",  gp$cell_size)
      ncdf4::ncatt_put(nc, 0, "date_start",        as.character(min(dates)))
      ncdf4::ncatt_put(nc, 0, "date_end",          as.character(max(dates)))
      ncdf4::ncatt_put(nc, 0, "n_dates",           ntime)
      ncdf4::ncatt_put(nc, 0, "n_cells_total",     nrow(fg))
      ncdf4::ncatt_put(nc, 0, "n_cells_labeled",   sum(fg$is_labeled == 1L, na.rm = TRUE))
      ncdf4::ncatt_put(nc, 0, "n_cells_unlabeled", sum(fg$is_labeled == 0L, na.rm = TRUE))
      ncdf4::ncatt_put(nc, 0, "missing_value",     mv)
      ncdf4::ncatt_put(nc, 0, "comment", paste0(
        "3-D grid (lon x lat x time). grid_id and is_labeled are 2-D (static across dates). ",
        "dwell_time/pm25/no2/indoor_pct are 3-D — one slice per date. ",
        "Unlabeled cells (is_labeled=0) have missing_value in all 3-D variables. ",
        "time axis uses CF 'days since 1970-01-01' for GEE compatibility."))

      ncdf4::nc_close(nc)
    }
  )

  output$download_raster_csv <- downloadHandler(
    filename = function() paste0("exposure_grid_daily_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$analysis_result, rv$analysis_result$daily_grid)
      data.table::fwrite(rv$analysis_result$daily_grid, file)
    }
  )

  output$download_full_grid_csv <- downloadHandler(
    filename = function() paste0("exposure_grid_daily_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$analysis_result, rv$analysis_result$daily_grid)
      data.table::fwrite(rv$analysis_result$daily_grid, file)
    }
  )

  # ── Animation: render frames and serve as auto-playing HTML ────────────────
  # ── Animation: render PNG frames → HTML slideshow (no system deps needed) ──
  anim_render_status <- reactiveVal("")   # "" | "rendering" | "ready" | "error: ..."
  anim_frame_count   <- reactiveVal(0L)
  anim_frame_ms      <- reactiveVal(200L)

  output$animation_video_ui <- renderUI({
    status <- anim_render_status()
    n      <- anim_frame_count()
    ms     <- anim_frame_ms()
    if (status == "ready" && n > 0) {
      # Build JS-driven slideshow from pre-rendered PNG frames served from /anim_frames/
      frames_js <- paste0('["', paste(sprintf("anim_frames/frame_%04d.png", seq_len(n)), collapse='","'), '"]')
      tags$div(
        tags$div(id = "spectra_anim_wrap",
          style = "position:relative; background:#061120; border-radius:10px; overflow:hidden;",
          tags$img(id = "spectra_anim_img",
                   src  = sprintf("anim_frames/frame_%04d.png", 1L),
                   style= "width:100%; display:block; border-radius:10px;"),
          tags$div(id = "spectra_anim_label",
                   style= paste0("position:absolute; top:8px; left:12px; color:#F97316;",
                                 " font-family:'Segoe UI',sans-serif; font-size:13px;",
                                 " font-weight:700; background:rgba(0,0,0,0.5);",
                                 " padding:2px 8px; border-radius:4px;"),
                   "Frame 1")
        ),
        tags$div(style = "text-align:center; margin-top:8px; display:flex; gap:8px; justify-content:center;",
          tags$button(id = "anim_play_btn", onclick = "spectraAnimPlay()",
                      style = paste0("background:linear-gradient(135deg,#0077B6,#0E2240);",
                                     " color:#F97316; border:1px solid #F97316; border-radius:6px;",
                                     " padding:6px 20px; font-weight:700; cursor:pointer;"),
                      HTML('<i class="fa fa-play"></i> Play')),
          tags$button(id = "anim_pause_btn", onclick = "spectraAnimPause()",
                      style = paste0("background:rgba(255,255,255,0.1); color:white;",
                                     " border:1px solid rgba(255,255,255,0.3); border-radius:6px;",
                                     " padding:6px 20px; font-weight:700; cursor:pointer;"),
                      HTML('<i class="fa fa-pause"></i> Pause')),
          tags$button(id = "anim_restart_btn", onclick = "spectraAnimRestart()",
                      style = paste0("background:rgba(255,255,255,0.1); color:white;",
                                     " border:1px solid rgba(255,255,255,0.3); border-radius:6px;",
                                     " padding:6px 20px; font-weight:700; cursor:pointer;"),
                      HTML('<i class="fa fa-redo"></i> Restart'))
        ),
        tags$script(HTML(sprintf("
          (function() {
            var frames = %s;
            var ms     = %d;
            var cur    = 0;
            var timer  = null;
            function showFrame(i) {
              cur = i %% frames.length;
              var img = document.getElementById('spectra_anim_img');
              var lbl = document.getElementById('spectra_anim_label');
              if (img) img.src = frames[cur];
              if (lbl) lbl.innerText = 'Frame ' + (cur+1) + ' / ' + frames.length;
            }
            window.spectraAnimPlay = function() {
              if (timer) return;
              timer = setInterval(function() { showFrame(cur + 1); }, ms);
            };
            window.spectraAnimPause = function() {
              clearInterval(timer); timer = null;
            };
            window.spectraAnimRestart = function() {
              clearInterval(timer); timer = null; showFrame(0);
            };
            // auto-play on load
            setTimeout(window.spectraAnimPlay, 300);
          })();
        ", frames_js, ms)))
      )
    } else {
      tags$div(class = "alert alert-info", style = "margin:20px; font-size:13px;",
        icon("film"),
        HTML(" Configure settings on the right, then click <b>Render Animation</b>.
              Frames are rendered and played directly in the browser — no extra software needed.")
      )
    }
  })

  output$anim_render_status_ui <- renderUI({
    s <- anim_render_status()
    n <- anim_frame_count()
    if (s == "rendering") {
      span(icon("spinner", class = "fa-spin"), " Rendering frames\u2026", style = "color:#F97316;")
    } else if (s == "ready") {
      span(icon("check-circle"), sprintf(" %d frames ready \u2014 auto-playing above", n), style = "color:#27ae60;")
    } else if (nchar(s) > 0) {
      span(icon("exclamation-triangle"), s, style = "color:#e74c3c;")
    } else NULL
  })

  observeEvent(input$btn_render_animation, {
    req(rv$analysis_result)
    anim_render_status("rendering")
    anim_frame_count(0L)

    dt <- as.data.table(rv$analysis_result$points_classified)
    req(nrow(dt) > 0)

    tryCatch({
      color_col <- input$anim_color     %||% "pm25"
      fps       <- as.integer(input$anim_fps %||% 5L)
      dur_sec   <- as.integer(input$anim_duration %||% 30L)
      n_frames  <- fps * dur_sec

      dt[, time_hour := lubridate::floor_date(date, "hour")]
      time_steps <- sort(unique(dt$time_hour))
      if (length(time_steps) > n_frames) {
        idx <- round(seq(1, length(time_steps), length.out = n_frames))
        time_steps <- time_steps[idx]
      }
      n_frames <- length(time_steps)

      out_dir <- file.path(getwd(), "www", "anim_frames")
      unlink(out_dir, recursive = TRUE)
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

      exp_label <- rv$exp1_label %||% "Exposure 1"

      withProgress(message = "Rendering animation frames\u2026", value = 0, {
        for (i in seq_along(time_steps)) {
          t   <- time_steps[i]
          sub <- dt[time_hour == t]
          if (nrow(sub) == 0) next

          if (color_col == "class") {
            p <- ggplot2::ggplot(sub, ggplot2::aes(x = lon, y = lat, color = corrected_class)) +
              ggplot2::scale_color_manual(values = c(indoor = "#F97316", outdoor = "#00B4D8"),
                                          name = "Class", na.value = "#888")
          } else {
            p <- ggplot2::ggplot(sub, ggplot2::aes(x = lon, y = lat,
                   color = if (color_col %in% names(sub)) .data[[color_col]] else lat)) +
              ggplot2::scale_color_viridis_c(option = "plasma", name = exp_label, na.value = "#888")
          }

          p <- p +
            ggplot2::geom_point(size = 1.8, alpha = 0.8) +
            ggplot2::labs(
              title    = format(t, "%Y-%m-%d  %H:%M"),
              subtitle = sprintf("Frame %d / %d", i, n_frames),
              x = "Longitude", y = "Latitude"
            ) +
            ggplot2::theme_minimal(base_size = 10) +
            ggplot2::theme(
              plot.background  = ggplot2::element_rect(fill = "#061120", color = NA),
              panel.background = ggplot2::element_rect(fill = "#061120", color = NA),
              panel.grid       = ggplot2::element_line(color = "rgba(255,255,255,0.06)"),
              text             = ggplot2::element_text(color = "white", family = "sans"),
              axis.text        = ggplot2::element_text(color = "#aaaaaa"),
              plot.title       = ggplot2::element_text(color = "#F97316", face = "bold"),
              plot.subtitle    = ggplot2::element_text(color = "#aaaaaa", size = 8),
              legend.background = ggplot2::element_rect(fill = "transparent"),
              legend.text      = ggplot2::element_text(color = "white")
            )

          fpath <- file.path(out_dir, sprintf("frame_%04d.png", i))
          ggplot2::ggsave(fpath, p, width = 8, height = 5, dpi = 80, bg = "#061120")
          incProgress(1 / n_frames)
        }
      })

      actual_frames <- length(list.files(out_dir, pattern = "\\.png$"))
      anim_frame_count(actual_frames)
      anim_frame_ms(as.integer(1000L / fps))
      anim_render_status("ready")

    }, error = function(e) {
      anim_render_status(paste("Render failed:", conditionMessage(e)))
    })
  })

  # ====================================================================
  # Trajectory Lines
  # ====================================================================
  
  output$traj_pid_select_ui <- renderUI({
    req(rv$analysis_result)
    pids <- sort(unique(rv$analysis_result$points_classified$PID))
    default <- head(pids, 3)
    pickerInput("traj_pids", "Select Participants",
                choices = pids, selected = default,
                options = list(`actions-box` = TRUE, `live-search` = TRUE),
                multiple = TRUE)
  })
  
  output$trajectory_map <- renderLeaflet({
    req(rv$analysis_result, input$traj_pids)
    
    dt <- as.data.table(rv$analysis_result$points_classified)
    dt <- dt[PID %in% input$traj_pids]
    req(nrow(dt) > 0)
    
    setorder(dt, PID, date)
    
    tile <- input$traj_map_style %||% "CartoDB.Positron"
    
    m <- leaflet() |>
      addProviderTiles(tile)
    
    # Generate colour palette for participants
    pids <- unique(dt$PID)
    pal_colors <- rainbow(length(pids))
    names(pal_colors) <- pids
    
    color_by <- input$traj_color_by %||% "Participant"
    
    for (pid in pids) {
      pid_dt <- dt[PID == pid]
      if (nrow(pid_dt) < 2) next
      
      # Determine line colour
      if (color_by == "Participant") {
        line_color <- pal_colors[pid]
      } else if (color_by == "Day") {
        line_color <- pal_colors[pid]  # Different days handled below
      } else {
        line_color <- pal_colors[pid]
      }
      
      if (color_by == "Day") {
        days       <- unique(pid_dt$date_only)
        day_colors <- setNames(rainbow(length(days)), as.character(days))
        pid_dt[, .day_col := day_colors[as.character(date_only)]]
        # One addPolylines call per unique colour (day), not per row
        for (d in days) {
          day_dt <- pid_dt[date_only == d]
          if (nrow(day_dt) < 2) next
          m <- m |> addPolylines(
            lng = day_dt$lon, lat = day_dt$lat,
            color = day_colors[as.character(d)], weight = 3, opacity = 0.8,
            group = "Trajectory Lines",
            popup = paste0("<b>", pid, "</b> - ", d)
          )
        }
      } else if (color_by == "Indoor/Outdoor") {
        # Two-segment approach: indoor vs outdoor
        m <- m |>
          addPolylines(
            lng = pid_dt$lon, lat = pid_dt$lat,
            color = line_color, weight = 3, opacity = 0.7,
            group = "Trajectory Lines",
            popup = paste0("<b>", pid, "</b>")
          )
      } else {
        m <- m |>
          addPolylines(
            lng = pid_dt$lon, lat = pid_dt$lat,
            color = line_color, weight = 3, opacity = 0.8,
            group = "Trajectory Lines",
            popup = paste0("<b>", pid, "</b>")
          )
      }
      
      # Add points if requested
      if (isTRUE(input$traj_show_points)) {
        point_color <- if (color_by == "Indoor/Outdoor") {
          ifelse(pid_dt$is_indoor, "#e74c3c", "#2ecc71")
        } else {
          line_color
        }
        
        m <- m |>
          addCircleMarkers(
            lng = pid_dt$lon, lat = pid_dt$lat,
            radius = 2, color = point_color,
            fillOpacity = 0.6, stroke = FALSE,
            group = "GPS Points",
            popup = paste0(
              "<b>PID:</b> ", pid_dt$PID, "<br>",
              "<b>Time:</b> ", pid_dt$date, "<br>",
              "<b>Class:</b> ", ifelse(pid_dt$is_indoor, "Indoor", "Outdoor"), "<br>",
              "<b>Speed:</b> ", round(pid_dt$speed_kmh, 1), " km/h"
            )
          )
      }
    }
    
    # Fit bounds
    m <- m |>
      fitBounds(
        lng1 = min(dt$lon), lat1 = min(dt$lat),
        lng2 = max(dt$lon), lat2 = max(dt$lat)
      ) |>
      addLayersControl(
        overlayGroups = c("Trajectory Lines", "GPS Points"),
        options = layersControlOptions(collapsed = FALSE)
      )

    m
  })

  output$traj_info_ui <- renderUI({
    req(rv$analysis_result, input$traj_pids)
    dt <- as.data.table(rv$analysis_result$points_classified)
    dt <- dt[PID %in% input$traj_pids]
    
    div(
      h5("Trajectory Info"),
      p(icon("users"), paste0(length(input$traj_pids), " participant(s)")),
      p(icon("map-marker-alt"), paste0(nrow(dt), " total points")),
      p(icon("calendar"), paste0(uniqueN(dt$date_only), " day(s)")),
      hr(),
      p(style = "font-size: 0.85em; color: #888;",
        "Lines connect GPS points chronologically. ",
        "Use to identify commuting routes and activity spaces.")
    )
  })
  
  # ====================================================================
  # Quality Assessment
  # ====================================================================
  
  output$qa_overall_ui <- renderUI({
    req(rv$analysis_result)
    dt <- as.data.table(rv$analysis_result$points_classified)
    
    if (is.null(rv$truth_field) || !(rv$truth_field %in% names(dt))) {
      return(div(class = "alert alert-info",
                 icon("info-circle"),
                 "No ground-truth data available. Upload data with ground-truth labels in Section 5."))
    }
    
    metrics <- confusion_metrics(dt$is_indoor, dt[[rv$truth_field]])
    
    div(
      fluidRow(
        column(6, div(class = "metric-card",
                      div(class = "metric-value", percent(metrics$accuracy, accuracy = 0.1)),
                      div(class = "metric-label", "Accuracy"))),
        column(6, div(class = "metric-card",
                      div(class = "metric-value", percent(metrics$precision, accuracy = 0.1)),
                      div(class = "metric-label", "Precision")))
      ),
      fluidRow(
        column(6, div(class = "metric-card",
                      div(class = "metric-value", percent(metrics$recall, accuracy = 0.1)),
                      div(class = "metric-label", "Recall"))),
        column(6, div(class = "metric-card",
                      div(class = "metric-value", percent(metrics$f1, accuracy = 0.1)),
                      div(class = "metric-label", "F1-Score")))
      ),
      fluidRow(
        column(6, div(class = "metric-card",
                      div(class = "metric-value", percent(metrics$specificity, accuracy = 0.1)),
                      div(class = "metric-label", "Specificity"))),
        column(6, div(class = "metric-card",
                      div(class = "metric-value", format(metrics$n_eval, big.mark = ",")),
                      div(class = "metric-label", "Eval Points")))
      ),
      hr(),
      h5("Confusion Matrix Values"),
      tableOutput("confusion_table")
    )
  })
  
  output$confusion_table <- renderTable({
    req(rv$analysis_result, rv$truth_field)
    dt <- as.data.table(rv$analysis_result$points_classified)
    if (!(rv$truth_field %in% names(dt))) return(NULL)
    
    metrics <- confusion_metrics(dt$is_indoor, dt[[rv$truth_field]])
    
    data.frame(
      ` ` = c("Actual Indoor", "Actual Outdoor"),
      `Predicted Indoor` = c(metrics$tp, metrics$fp),
      `Predicted Outdoor` = c(metrics$fn, metrics$tn),
      check.names = FALSE
    )
  }, bordered = TRUE, striped = TRUE)
  
  output$plot_confusion_matrix <- renderPlotly({
    req(rv$analysis_result, rv$truth_field)
    dt <- as.data.table(rv$analysis_result$points_classified)
    if (!(rv$truth_field %in% names(dt))) {
      return(plot_ly() |> layout(title = "No ground-truth data available"))
    }
    
    metrics <- confusion_metrics(dt$is_indoor, dt[[rv$truth_field]])
    
    m <- matrix(c(metrics$tp, metrics$fn, metrics$fp, metrics$tn),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("Indoor", "Outdoor"), c("Indoor", "Outdoor")))
    
    plot_ly(z = m, x = colnames(m), y = rownames(m), type = "heatmap",
            colors = "Blues", text = m, texttemplate = "%{text}",
            showscale = TRUE) |>
      layout(title = "Confusion Matrix",
             xaxis = list(title = "Predicted"),
             yaxis = list(title = "Actual"))
  })
  
  output$qa_per_person_table <- DT::renderDataTable({
    req(rv$analysis_result, rv$truth_field)
    dt <- as.data.table(rv$analysis_result$points_classified)
    if (!(rv$truth_field %in% names(dt))) return(NULL)
    
    per_person <- dt[, {
      m <- confusion_metrics(is_indoor, get(rv$truth_field))
      .(accuracy = m$accuracy, precision = m$precision,
        recall = m$recall, specificity = m$specificity,
        f1 = m$f1, n_eval = m$n_eval)
    }, by = PID]
    
    DT::datatable(per_person, options = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE) |>
      DT::formatPercentage(c("accuracy", "precision", "recall", "specificity", "f1"), 1)
  })
  
  output$plot_per_person_accuracy <- renderPlotly({
    req(rv$analysis_result, rv$truth_field)
    dt <- as.data.table(rv$analysis_result$points_classified)
    if (!(rv$truth_field %in% names(dt))) {
      return(plot_ly() |> layout(title = "No ground-truth data"))
    }
    
    per_person <- dt[, {
      m <- confusion_metrics(is_indoor, get(rv$truth_field))
      .(Accuracy = m$accuracy, Precision = m$precision,
        Recall = m$recall, F1 = m$f1)
    }, by = PID]
    
    top20 <- per_person[order(-Accuracy)][1:min(20, nrow(per_person))]
    
    plot_ly(top20, x = ~PID) |>
      add_bars(y = ~Accuracy, name = "Accuracy", marker = list(color = "#2a5298")) |>
      add_bars(y = ~Precision, name = "Precision", marker = list(color = "#4a90e2")) |>
      add_bars(y = ~Recall, name = "Recall", marker = list(color = "#7ab5e8")) |>
      add_bars(y = ~F1, name = "F1", marker = list(color = "#a8d1f0")) |>
      layout(title = "Per-Participant Classification Metrics",
             barmode = "group",
             xaxis = list(title = "Participant ID"),
             yaxis = list(title = "Score"),
             plot_bgcolor = "rgba(0,0,0,0)",
             paper_bgcolor = "rgba(0,0,0,0)")
  })
  
  # ====================================================================
  # Data Tables
  # ====================================================================
  
  output$table_classified <- DT::renderDataTable({
    req(rv$analysis_result)
    dt <- rv$analysis_result$points_classified
    display_cols <- intersect(c("PID", "date", "lat", "lon", "pm25", "no2",
                                "speed_kmh", "dwell_min", "is_indoor",
                                "indoor_confidence", "date_only"),
                              names(dt))
    DT::datatable(dt[, ..display_cols], 
                  filter = "top",
                  options = list(pageLength = 20, scrollX = TRUE,
                                 searching = TRUE, search = list(regex = TRUE)),
                  rownames = FALSE)
  })
  
  output$table_grid <- DT::renderDataTable({
    req(rv$analysis_result)
    DT::datatable(rv$analysis_result$person_day_cell,
                  filter = "top",
                  options = list(pageLength = 20, scrollX = TRUE,
                                 searching = TRUE, search = list(regex = TRUE)),
                  rownames = FALSE)
  })
  
  output$table_daily <- DT::renderDataTable({
    req(rv$analysis_result)
    DT::datatable(rv$analysis_result$person_day,
                  filter = "top",
                  options = list(pageLength = 20, scrollX = TRUE,
                                 searching = TRUE, search = list(regex = TRUE)),
                  rownames = FALSE)
  })

  output$table_participant_exposure <- DT::renderDataTable({
    req(rv$analysis_result)
    pe <- rv$analysis_result$participant_exposure
    if (is.null(pe) || nrow(pe) == 0) return(NULL)
    dt_out <- DT::datatable(pe,
                  filter = "top",
                  options = list(pageLength = 20, scrollX = TRUE,
                                 searching = TRUE, search = list(regex = TRUE)),
                  rownames = FALSE) |>
      DT::formatRound(c("pm25_avg_daily","pm25_max_daily","no2_avg_daily","no2_max_daily",
                         "dwell_avg_day_min","n_cells_avg"), 2) |>
      DT::formatStyle("pm25_exceeds_who_annual",
                      backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#ffd6d6", "#d4f5e2"))) |>
      DT::formatStyle("no2_exceeds_who_annual",
                      backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#ffd6d6", "#d4f5e2")))
    dt_out
  })
  
  # ====================================================================
  # Export / Downloads
  # ====================================================================
  
  # ── Grid RDS download — full_grid data.table only ─────────────────────────
  output$download_result_rds <- downloadHandler(
    filename = function() paste0("SPECTRA_grid_", Sys.Date(), ".rds"),
    content = function(file) {
      req(rv$analysis_result, rv$analysis_result$full_grid)
      saveRDS(rv$analysis_result$full_grid, file, compress = "gzip")
    }
  )

  output$download_classified_csv <- downloadHandler(
    filename = function() paste0("classified_points_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$analysis_result)
      fwrite(rv$analysis_result$points_classified, file)
    }
  )
  
  output$download_grid_csv <- downloadHandler(
    filename = function() paste0("grid_cells_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$analysis_result)
      fwrite(rv$analysis_result$person_day_cell, file)
    }
  )
  
  output$download_daily_csv <- downloadHandler(
    filename = function() paste0("daily_summaries_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$analysis_result)
      fwrite(rv$analysis_result$person_day, file)
    }
  )

  output$download_participant_csv <- downloadHandler(
    filename = function() paste0("participant_exposure_who_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$analysis_result)
      fwrite(rv$analysis_result$participant_exposure, file)
    }
  )
  
  output$download_points_geojson <- downloadHandler(
    filename = function() paste0("classified_points_", Sys.Date(), ".geojson"),
    content = function(file) {
      req(rv$analysis_result)
      dt <- rv$analysis_result$points_classified
      pts <- sf::st_as_sf(dt, coords = c("lon", "lat"), crs = EPSG_WGS84, remove = FALSE)
      sf::st_write(pts, file, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
    }
  )
  
  output$download_stats_json <- downloadHandler(
    filename = function() paste0("analysis_stats_", Sys.Date(), ".json"),
    content = function(file) {
      req(rv$analysis_result)
      jsonlite::write_json(rv$analysis_result$stats, file, pretty = TRUE, auto_unbox = TRUE)
    }
  )
  
  output$download_summary <- downloadHandler(
    filename = function() paste0("analysis_summary_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$analysis_result)
      stats <- rv$analysis_result$stats
      summary_df <- data.frame(
        Metric = names(stats),
        Value = as.character(unlist(stats))
      )
      fwrite(summary_df, file)
    }
  )

  # ====================================================================
  # Master Excel Export
  # ====================================================================

  output$download_master_excel <- downloadHandler(
    filename = function() paste0("SPECTRA_master_export_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(rv$analysis_result)

      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        showNotification("openxlsx package required. Run: install.packages('openxlsx')",
                         type = "error", duration = 10)
        return(NULL)
      }

      result <- rv$analysis_result
      dt_seg <- as.data.table(result$points_classified)

      wb <- openxlsx::createWorkbook()

      add_sheet <- function(name, data) {
        openxlsx::addWorksheet(wb, name)
        openxlsx::writeData(wb, name, as.data.frame(data))
        openxlsx::freezePane(wb, name, firstRow = TRUE)
        hs <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#1e3c72",
                                     halign = "CENTER", textDecoration = "Bold",
                                     border = "Bottom")
        openxlsx::addStyle(wb, name, style = hs, rows = 1,
                           cols = seq_len(ncol(data)), gridExpand = TRUE)
      }

      add_sheet("Classified_Points",      result$points_classified)
      add_sheet("Grid_Cells",             result$person_day_cell)
      add_sheet("Daily_Summaries",        result$person_day)
      if (!is.null(result$participant_exposure) && nrow(result$participant_exposure) > 0)
        add_sheet("Participant_WHO_Exposure", result$participant_exposure)

      add_sheet("Full_Segmented",     dt_seg)
      add_sheet("Outdoor_Corrected",  dt_seg[corrected_class == "Outdoor"])
      add_sheet("Indoor_Corrected",   dt_seg[corrected_class == "Indoor"])

      stats_df <- data.frame(Metric = names(result$stats),
                              Value  = as.character(unlist(result$stats)))
      add_sheet("Summary_Statistics", stats_df)

      if (!is.null(result$validation_stats) && nrow(result$validation_stats) > 0)
        add_sheet("Validation_Pipeline", result$validation_stats)

      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

  # ====================================================================
  # Comprehensive HTML Report
  # ====================================================================

  output$download_full_report <- downloadHandler(
    filename = function() paste0("SPECTRA_report_", Sys.Date(), ".html"),
    content = function(file) {
      req(rv$analysis_result)

      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        showNotification("rmarkdown package required. Run: install.packages('rmarkdown')",
                         type = "error", duration = 10)
        return(NULL)
      }

      tmp_data <- tempfile(fileext = ".rds")
      saveRDS(rv$analysis_result, tmp_data)
      on.exit(unlink(tmp_data), add = TRUE)

      template_path <- file.path(dirname(sys.frame(1)$ofile %||% getwd()),
                                  "report_template.Rmd")
      if (!file.exists(template_path))
        template_path <- file.path(getwd(), "report_template.Rmd")
      if (!file.exists(template_path))
        template_path <- file.path(getwd(), "shiny_app", "report_template.Rmd")

      if (!file.exists(template_path)) {
        showNotification("Cannot find report_template.Rmd",
                         type = "error", duration = 10)
        return(NULL)
      }

      tmp_rmd <- tempfile(fileext = ".Rmd")
      file.copy(template_path, tmp_rmd)
      on.exit(unlink(tmp_rmd), add = TRUE)

      withProgress(message = "Generating report...", value = 0.5, {
        rmarkdown::render(
          input       = tmp_rmd,
          output_file = file,
          params      = list(data_path = tmp_data),
          envir       = new.env(parent = globalenv()),
          quiet       = TRUE
        )
      })
    }
  )

  # ====================================================================
  # Optional Field Mapping UI (Country, Temperature, Humidity)
  # ====================================================================

  output$country_field_ui <- renderUI({
    req(rv$columns)
    choices <- c("None", rv$columns)
    sel <- "None"
    for (p in c("^country$", "country", "Country", "region")) {
      m <- grep(p, rv$columns, ignore.case = TRUE, value = TRUE)
      if (length(m) > 0) { sel <- m[1]; break }
    }
    selectInput("country_field", "Country (Opt)", choices = choices, selected = sel)
  })

  output$temp_field_ui <- renderUI({
    req(rv$columns)
    choices <- c("None", rv$columns)
    sel <- "None"
    for (p in c("temperature", "temp_c", "temp")) {
      m <- grep(p, rv$columns, ignore.case = TRUE, value = TRUE)
      if (length(m) > 0) { sel <- m[1]; break }
    }
    selectInput("temp_field", "Temperature (Opt)", choices = choices, selected = sel)
  })

  output$humid_field_ui <- renderUI({
    req(rv$columns)
    choices <- c("None", rv$columns)
    sel <- "None"
    for (p in c("humidity", "humid", "rh")) {
      m <- grep(p, rv$columns, ignore.case = TRUE, value = TRUE)
      if (length(m) > 0) { sel <- m[1]; break }
    }
    selectInput("humid_field", "Humidity (Opt)", choices = choices, selected = sel)
  })

  output$heart_rate_field_ui <- renderUI({
    req(rv$columns)
    choices <- c("None", rv$columns)
    sel <- "None"
    for (p in c("heart_rate", "heartrate", "hr", "bpm", "pulse")) {
      m <- grep(p, rv$columns, ignore.case = TRUE, value = TRUE)
      if (length(m) > 0) { sel <- m[1]; break }
    }
    selectInput("heart_rate_field", "Heart Rate (Opt)", choices = choices, selected = sel)
  })

  output$steps_field_ui <- renderUI({
    req(rv$columns)
    choices <- c("None", rv$columns)
    sel <- "None"
    for (p in c("steps", "step_count", "cadence")) {
      m <- grep(p, rv$columns, ignore.case = TRUE, value = TRUE)
      if (length(m) > 0) { sel <- m[1]; break }
    }
    selectInput("steps_field", "Steps / Cadence (Opt)", choices = choices, selected = sel)
  })

  output$activity_field_ui <- renderUI({
    req(rv$columns)
    choices <- c("None", rv$columns)
    sel <- "None"
    for (p in c("activity", "activity_type", "sport", "exercise")) {
      m <- grep(p, rv$columns, ignore.case = TRUE, value = TRUE)
      if (length(m) > 0) { sel <- m[1]; break }
    }
    selectInput("activity_field", "Activity Type (Opt)", choices = choices, selected = sel)
  })

  # ====================================================================
  # Segmentation & Validation Tab - Server Logic
  # ====================================================================

  seg_data <- reactive({
    req(rv$analysis_result)
    as.data.table(rv$analysis_result$points_classified)
  })

  # --- Validation Pipeline Stats Table ---
  output$seg_validation_stats_table <- DT::renderDataTable({
    req(rv$analysis_result)
    stats <- rv$analysis_result$validation_stats
    DT::datatable(stats, options = list(dom = "t", pageLength = 25),
                  rownames = FALSE) |>
      DT::formatStyle("Step", fontWeight = "bold")
  })

  # --- Confusion Heatmap: Speed vs OSM ---
  output$seg_confusion_heatmap <- renderPlotly({
    dt <- seg_data()
    req("speed_class" %in% names(dt), "osm_label" %in% names(dt))
    dt_valid <- dt[!is.na(osm_label)]
    req(nrow(dt_valid) > 0)

    cm <- table(
      Speed = factor(dt_valid$speed_class, levels = c("Indoor", "Outdoor")),
      OSM = factor(dt_valid$osm_label, levels = c("Indoor", "Outdoor"))
    )
    cm_pct <- prop.table(cm, margin = 1) * 100

    labels_mat <- matrix(
      paste0(as.numeric(cm), "\n(", round(as.numeric(cm_pct), 1), "%)"),
      nrow = 2
    )

    plot_ly(
      x = colnames(cm), y = rownames(cm),
      z = matrix(as.numeric(cm_pct), nrow = 2),
      type = "heatmap", colors = "Blues",
      text = labels_mat, texttemplate = "%{text}",
      showscale = TRUE
    ) |>
      layout(
        title = "Speed vs OSM Classification\n(Row-wise %)",
        xaxis = list(title = "OSM Ground Truth"),
        yaxis = list(title = "Speed-Based Prediction")
      )
  })

  # --- Misclassification Pie ---
  output$seg_misclass_pie <- renderPlotly({
    dt <- seg_data()
    req("misclass_type" %in% names(dt))
    dt_valid <- dt[!is.na(misclass_type)]
    req(nrow(dt_valid) > 0)

    counts <- dt_valid[, .N, by = misclass_type]

    plot_ly(counts, labels = ~misclass_type, values = ~N, type = "pie",
            marker = list(colors = c("#2ecc71", "#e67e22", "#9b59b6")),
            hole = 0.4) |>
      layout(title = "Misclassification Breakdown")
  })

  # --- Methods Comparison Bar ---
  output$seg_method_comparison <- renderPlotly({
    dt <- seg_data()
    req("speed_class" %in% names(dt), "corrected_class" %in% names(dt))

    speed_in <- sum(dt$speed_class == "Indoor")
    speed_out <- sum(dt$speed_class == "Outdoor")
    corr_in <- sum(dt$corrected_class == "Indoor")
    corr_out <- sum(dt$corrected_class == "Outdoor")

    comparison <- data.frame(
      Method = rep(c("Speed-Only", "Corrected (Speed+OSM)"), each = 2),
      Class = rep(c("Indoor", "Outdoor"), 2),
      Count = c(speed_in, speed_out, corr_in, corr_out)
    )

    plot_ly(comparison, x = ~Method, y = ~Count, color = ~Class, type = "bar",
            colors = c("Indoor" = "#2a5298", "Outdoor" = "#e67e22")) |>
      layout(title = "Speed-Only vs Corrected Classification",
             barmode = "group",
             xaxis = list(title = ""),
             yaxis = list(title = "Number of Points"))
  })

  # --- Indoor/Outdoor Bar Chart ---
  output$seg_indoor_outdoor_bar <- renderPlotly({
    dt <- seg_data()
    req("corrected_class" %in% names(dt))

    counts <- dt[, .(N = .N), by = corrected_class]
    counts[, pct := round(N / sum(N) * 100, 1)]

    colors <- c("Indoor" = "#1b9e77", "Outdoor" = "#d95f02")

    plot_ly(counts, x = ~corrected_class, y = ~N, type = "bar",
            text = ~paste0(format(N, big.mark = ","), " (", pct, "%)"),
            textposition = "auto",
            marker = list(color = colors[counts$corrected_class])) |>
      layout(title = "Indoor vs Outdoor - Corrected Dataset",
             xaxis = list(title = "Exposure Type"),
             yaxis = list(title = "Number of Points"))
  })

  # --- Time Series ---
  output$seg_timeseries <- renderPlotly({
    dt <- seg_data()
    req("corrected_class" %in% names(dt), "date" %in% names(dt))

    daily <- dt[, .(
      Indoor = sum(corrected_class == "Indoor"),
      Outdoor = sum(corrected_class == "Outdoor")
    ), by = .(day = as.Date(date))]
    daily[, Indoor_pct := round(Indoor / (Indoor + Outdoor) * 100, 1)]
    daily[, Outdoor_pct := 100 - Indoor_pct]

    plot_ly(daily, x = ~day) |>
      add_lines(y = ~Indoor_pct, name = "Indoor %",
                line = list(color = "#1b9e77")) |>
      add_lines(y = ~Outdoor_pct, name = "Outdoor %",
                line = list(color = "#d95f02")) |>
      layout(title = "Daily Indoor/Outdoor Proportion",
             xaxis = list(title = "Date"),
             yaxis = list(title = "Percentage (%)", range = c(0, 100)))
  })

  # --- Outdoor Movement Intensity ---
  output$seg_outdoor_movement <- renderPlotly({
    dt <- seg_data()
    req("outdoor_movement_type" %in% names(dt))

    outdoor <- dt[corrected_class == "Outdoor" & !is.na(outdoor_movement_type)]
    req(nrow(outdoor) > 0)

    counts <- outdoor[, .(N = .N), by = outdoor_movement_type]
    counts[, pct := round(N / sum(N) * 100, 1)]

    plot_ly(counts, x = ~outdoor_movement_type, y = ~N, type = "bar",
            text = ~paste0(pct, "%"), textposition = "auto",
            marker = list(color = c("#fee5d9", "#fcae91", "#fb6a4a"))) |>
      layout(title = "Outdoor Movement Intensity",
             xaxis = list(title = ""),
             yaxis = list(title = "Number of Points"))
  })

  # --- Outdoor Context ---
  output$seg_outdoor_context <- renderPlotly({
    dt <- seg_data()
    req("outdoor_context" %in% names(dt))

    outdoor <- dt[corrected_class == "Outdoor" & !is.na(outdoor_context)]
    req(nrow(outdoor) > 0)

    counts <- outdoor[, .(N = .N), by = outdoor_context]
    counts[, pct := round(N / sum(N) * 100, 1)]

    plot_ly(counts, y = ~outdoor_context, x = ~N, type = "bar",
            orientation = "h",
            text = ~paste0(pct, "%"), textposition = "auto",
            marker = list(color = viridis::viridis(nrow(counts)))) |>
      layout(title = "Outdoor Exposure Context",
             xaxis = list(title = "Number of Points"),
             yaxis = list(title = ""))
  })

  # --- Outdoor Environmental Stress ---
  output$seg_outdoor_stress <- renderPlotly({
    dt <- seg_data()
    if (!("outdoor_stress_level" %in% names(dt))) {
      return(plot_ly() |>
               layout(title = "No temperature/humidity data for stress analysis",
                      annotations = list(
                        text = "Map Temperature and Humidity fields to enable this chart",
                        x = 0.5, y = 0.5, showarrow = FALSE, xref = "paper", yref = "paper"
                      )))
    }

    outdoor <- dt[corrected_class == "Outdoor" &
                    !is.na(outdoor_stress_level) &
                    outdoor_stress_level != "Unknown Stress"]
    req(nrow(outdoor) > 0)

    counts <- outdoor[, .(N = .N), by = outdoor_stress_level]
    counts[, pct := round(N / sum(N) * 100, 1)]

    stress_colors <- c("Low Stress" = "#91bfdb",
                        "Moderate Stress" = "#ffffbf",
                        "High Stress" = "#fc8d59")

    plot_ly(counts, x = ~outdoor_stress_level, y = ~N, type = "bar",
            text = ~paste0(pct, "%"), textposition = "auto",
            marker = list(color = stress_colors[counts$outdoor_stress_level])) |>
      layout(title = "Outdoor Environmental Stress Levels",
             xaxis = list(title = ""),
             yaxis = list(title = "Number of Points"))
  })

  # --- Outdoor Data Table ---
  output$seg_outdoor_table <- DT::renderDataTable({
    dt <- seg_data()
    req("corrected_class" %in% names(dt))

    outdoor <- dt[corrected_class == "Outdoor"]
    display_cols <- intersect(
      c("PID", "date", "lat", "lon", "pm25", "no2", "speed_kmh",
        "dwell_min", "corrected_class", "speed_class", "osm_label",
        "misclass_type", "outdoor_movement_type", "outdoor_context",
        "mobility_theme", "temporal_theme"),
      names(outdoor)
    )
    DT::datatable(outdoor[, ..display_cols],
                  options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE)
  })

  # --- Helper: get selected exposure variable column and label ---
  seg_exposure_col <- reactive({
    pol <- input$seg_exposure %||% "pm25"
    label <- if (pol == "pm25") "PM2.5" else "NO2"
    unit <- if (pol == "pm25") "\u00b5g/m\u00b3" else "ppb"
    list(col = pol, label = label, unit = unit)
  })

  # --- Temporal Activity Pattern (Exposure by Time of Day) ---
  output$seg_temporal_pattern <- renderPlotly({
    dt <- seg_data()
    pol <- seg_exposure_col()
    req("temporal_theme" %in% names(dt), pol$col %in% names(dt))

    agg <- dt[!is.na(get(pol$col)), .(
      mean_exp = round(mean(get(pol$col), na.rm = TRUE), 2),
      median_exp = round(median(get(pol$col), na.rm = TRUE), 2),
      n = .N
    ), by = .(temporal_theme, exposure_type = corrected_class)]

    plot_ly(agg, x = ~temporal_theme, y = ~mean_exp, color = ~exposure_type,
            type = "bar",
            text = ~paste0("Mean: ", mean_exp, " ", pol$unit, "\nN=", n),
            textposition = "auto",
            colors = c("Indoor" = "#1b9e77", "Outdoor" = "#d95f02")) |>
      layout(title = paste(pol$label, "Exposure by Time of Day"),
             barmode = "group",
             xaxis = list(title = ""),
             yaxis = list(title = paste("Mean", pol$label, paste0("(", pol$unit, ")"))))
  })

  # --- Mobility Theme (Exposure by Mobility) ---
  output$seg_mobility_theme <- renderPlotly({
    dt <- seg_data()
    pol <- seg_exposure_col()
    req("mobility_theme" %in% names(dt), pol$col %in% names(dt))

    agg <- dt[!is.na(mobility_theme) & !is.na(get(pol$col)), .(
      mean_exp = round(mean(get(pol$col), na.rm = TRUE), 2),
      sd_exp = round(sd(get(pol$col), na.rm = TRUE), 2),
      n = .N
    ), by = .(mobility_theme, exposure_type = corrected_class)]

    plot_ly(agg, x = ~mobility_theme, y = ~mean_exp, color = ~exposure_type,
            type = "bar",
            error_y = list(type = "data", array = ~sd_exp, visible = TRUE),
            text = ~paste0("Mean: ", mean_exp, "\nSD: ", sd_exp, "\nN=", n),
            colors = c("Indoor" = "#1b9e77", "Outdoor" = "#d95f02")) |>
      layout(title = paste(pol$label, "Exposure by Mobility Theme"),
             barmode = "group",
             xaxis = list(title = ""),
             yaxis = list(title = paste("Mean", pol$label, paste0("(", pol$unit, ")"))))
  })

  # --- Weekday x Time-of-Day Exposure Heatmap ---
  output$seg_weekday_heatmap <- renderPlotly({
    dt <- seg_data()
    pol <- seg_exposure_col()
    req("weekday_name" %in% names(dt), "temporal_theme" %in% names(dt), pol$col %in% names(dt))

    agg <- dt[!is.na(get(pol$col)), .(
      mean_exp = round(mean(get(pol$col), na.rm = TRUE), 2)
    ), by = .(weekday_name, temporal_theme)]

    weekdays_order <- c("Monday", "Tuesday", "Wednesday", "Thursday",
                         "Friday", "Saturday", "Sunday")
    time_order <- c("Night (0-4h)", "Morning (5-11h)",
                     "Afternoon (12-16h)", "Evening (17-23h)")

    wide <- dcast(agg, weekday_name ~ temporal_theme, value.var = "mean_exp", fill = 0)
    wd_present <- intersect(weekdays_order, wide$weekday_name)
    wide <- wide[match(wd_present, wide$weekday_name), ]

    mat_cols <- intersect(time_order, names(wide))
    mat <- as.matrix(wide[, ..mat_cols])

    plot_ly(x = mat_cols, y = wide$weekday_name, z = mat,
            type = "heatmap", colors = "YlOrRd",
            hovertemplate = paste("Day: %{y}<br>Period: %{x}<br>Mean",
                                  pol$label, ": %{z}", pol$unit,
                                  "<extra></extra>")) |>
      layout(title = paste("Weekday x Time-of-Day Mean", pol$label, "Exposure"),
             xaxis = list(title = ""),
             yaxis = list(title = "", autorange = "reversed"))
  })

  # --- Duration Theme (Exposure) ---
  output$seg_duration_theme <- renderPlotly({
    dt <- seg_data()
    pol <- seg_exposure_col()
    req("duration_theme" %in% names(dt), pol$col %in% names(dt))

    agg <- dt[!is.na(duration_theme) & !is.na(get(pol$col)), .(
      mean_exp = round(mean(get(pol$col), na.rm = TRUE), 2),
      n = .N
    ), by = .(duration_theme, exposure_type = corrected_class)]

    plot_ly(agg, x = ~duration_theme, y = ~mean_exp, color = ~exposure_type,
            type = "bar",
            text = ~paste0("Mean: ", mean_exp, " ", pol$unit, "\nN=", n),
            textposition = "auto",
            colors = c("Indoor" = "#1b9e77", "Outdoor" = "#d95f02")) |>
      layout(title = paste(pol$label, "Exposure by Duration Theme"),
             barmode = "group",
             xaxis = list(title = ""),
             yaxis = list(title = paste("Mean", pol$label, paste0("(", pol$unit, ")"))))
  })

  # --- Hourly Exposure Profile ---
  output$seg_hourly_profile <- renderPlotly({
    dt <- seg_data()
    pol <- seg_exposure_col()
    req("hour_local" %in% names(dt), pol$col %in% names(dt))

    hourly <- dt[!is.na(get(pol$col)), .(
      mean_exp = round(mean(get(pol$col), na.rm = TRUE), 2),
      p25 = round(quantile(get(pol$col), 0.25, na.rm = TRUE), 2),
      p75 = round(quantile(get(pol$col), 0.75, na.rm = TRUE), 2)
    ), by = .(hour = hour_local, exposure_type = corrected_class)]

    plot_ly() |>
      add_trace(data = hourly[exposure_type == "Indoor"],
                x = ~hour, y = ~mean_exp, type = "scatter", mode = "lines+markers",
                name = "Indoor", line = list(color = "#1b9e77"),
                marker = list(color = "#1b9e77")) |>
      add_trace(data = hourly[exposure_type == "Outdoor"],
                x = ~hour, y = ~mean_exp, type = "scatter", mode = "lines+markers",
                name = "Outdoor", line = list(color = "#d95f02"),
                marker = list(color = "#d95f02")) |>
      layout(title = paste("Hourly", pol$label, "Exposure Profile (24h)"),
             xaxis = list(title = "Hour of Day", dtick = 1, range = c(-0.5, 23.5)),
             yaxis = list(title = paste("Mean", pol$label, paste0("(", pol$unit, ")"))))
  })

  # --- Peak Exposure Hours Summary ---
  output$seg_peak_hours <- renderPlotly({
    dt <- seg_data()
    pol <- seg_exposure_col()
    req("hour_local" %in% names(dt), pol$col %in% names(dt))

    hourly <- dt[!is.na(get(pol$col)), .(
      mean_exp = round(mean(get(pol$col), na.rm = TRUE), 2)
    ), by = .(hour = hour_local)]
    setorder(hourly, -mean_exp)

    # Color: top 5 hours = red, rest = blue
    hourly[, rank := seq_len(.N)]
    hourly[, bar_color := fifelse(rank <= 5, "#e74c3c", "#3498db")]
    setorder(hourly, hour)

    plot_ly(hourly, x = ~factor(hour), y = ~mean_exp, type = "bar",
            marker = list(color = ~bar_color),
            text = ~paste0(mean_exp, " ", pol$unit),
            textposition = "outside") |>
      layout(title = paste("Peak", pol$label, "Exposure Hours (Top 5 in Red)"),
             xaxis = list(title = "Hour of Day"),
             yaxis = list(title = paste("Mean", pol$label, paste0("(", pol$unit, ")"))))
  })

  # --- Climate by Exposure Type (Exposure Levels) ---
  output$seg_climate_exposure <- renderPlotly({
    dt <- seg_data()
    pol <- seg_exposure_col()
    req("climate_theme" %in% names(dt), pol$col %in% names(dt))

    dt_valid <- dt[climate_theme != "Unknown" & !is.na(get(pol$col))]
    if (nrow(dt_valid) == 0) {
      return(plot_ly() |>
               layout(title = "No climate data available",
                      annotations = list(
                        text = "Map Temperature and Humidity fields to enable this chart",
                        x = 0.5, y = 0.5, showarrow = FALSE, xref = "paper", yref = "paper"
                      )))
    }

    agg <- dt_valid[, .(
      mean_exp = round(mean(get(pol$col), na.rm = TRUE), 2),
      n = .N
    ), by = .(climate_theme, exposure_type = corrected_class)]

    plot_ly(agg, x = ~climate_theme, y = ~mean_exp, color = ~exposure_type,
            type = "bar",
            text = ~paste0("Mean: ", mean_exp, " ", pol$unit, "\nN=", n),
            colors = c("Indoor" = "#1b9e77", "Outdoor" = "#d95f02")) |>
      layout(title = paste(pol$label, "Exposure by Climate Condition"),
             barmode = "group",
             xaxis = list(title = "", tickangle = -30),
             yaxis = list(title = paste("Mean", pol$label, paste0("(", pol$unit, ")"))))
  })

  # --- Climate Pie (Mean Exposure by Climate) ---
  output$seg_climate_pie <- renderPlotly({
    dt <- seg_data()
    pol <- seg_exposure_col()
    req("climate_theme" %in% names(dt), pol$col %in% names(dt))

    agg <- dt[climate_theme != "Unknown" & !is.na(get(pol$col)),
              .(mean_exp = round(mean(get(pol$col), na.rm = TRUE), 2)),
              by = climate_theme]
    if (nrow(agg) == 0) {
      return(plot_ly() |>
               layout(title = "No climate data available"))
    }

    plot_ly(agg, labels = ~climate_theme, values = ~mean_exp, type = "pie",
            marker = list(colors = c("#91bfdb", "#fee090", "#5ab4ac", "#d8b365")),
            hole = 0.3,
            textinfo = "label+value",
            hovertemplate = paste("%{label}<br>Mean", pol$label, ": %{value}", pol$unit, "<extra></extra>")) |>
      layout(title = paste("Mean", pol$label, "by Climate Theme"))
  })

  # --- Misclassification Map ---
  output$seg_misclass_map <- renderLeaflet({
    dt <- seg_data()
    req("misclass_type" %in% names(dt))

    dt_valid <- dt[!is.na(misclass_type)]
    req(nrow(dt_valid) > 0)

    n <- min(nrow(dt_valid), 10000)
    if (nrow(dt_valid) > n) dt_valid <- dt_valid[sample(.N, n)]

    color_map <- c("Correct" = "#2ecc71",
                    "False Indoor" = "#e67e22",
                    "False Outdoor" = "#9b59b6")
    dt_valid[, point_color := color_map[misclass_type]]

    m <- leaflet() |>
      addProviderTiles("CartoDB.Positron")

    overlay_groups <- character(0)
    if (isTRUE(input$show_buildings_misclass)) {
      m <- .add_building_polygons(m, rv$buildings_sf)
      overlay_groups <- c(overlay_groups, "Building Polygons")
    }

    m <- m |>
      addCircleMarkers(
        data = dt_valid,
        lng = ~lon, lat = ~lat,
        radius = 3, color = ~point_color,
        fillOpacity = 0.6, stroke = FALSE,
        group = "Points",
        popup = ~paste0(
          "<b>PID:</b> ", PID, "<br>",
          "<b>Speed Class:</b> ", speed_class, "<br>",
          "<b>OSM Label:</b> ", osm_label, "<br>",
          "<b>Status:</b> ", misclass_type, "<br>",
          "<b>Corrected:</b> ", corrected_class
        )
      ) |>
      addLegend(position = "bottomright",
                colors = color_map, labels = names(color_map),
                title = "Classification") |>
      fitBounds(
        lng1 = min(dt_valid$lon), lat1 = min(dt_valid$lat),
        lng2 = max(dt_valid$lon), lat2 = max(dt_valid$lat)
      )

    if (length(overlay_groups) > 0) {
      m <- m |>
        addLayersControl(overlayGroups = overlay_groups,
                         options = layersControlOptions(collapsed = FALSE))
    }

    m
  })

  # --- Download Handlers for Corrected Datasets ---
  output$download_outdoor_corrected <- downloadHandler(
    filename = function() paste0("outdoor_corrected_", Sys.Date(), ".csv"),
    content = function(file) {
      dt <- seg_data()
      outdoor <- dt[corrected_class == "Outdoor"]
      fwrite(outdoor, file)
    }
  )

  output$download_indoor_corrected <- downloadHandler(
    filename = function() paste0("indoor_corrected_", Sys.Date(), ".csv"),
    content = function(file) {
      dt <- seg_data()
      indoor <- dt[corrected_class == "Indoor"]
      fwrite(indoor, file)
    }
  )

  output$download_full_segmented <- downloadHandler(
    filename = function() paste0("full_segmented_", Sys.Date(), ".csv"),
    content = function(file) {
      dt <- seg_data()
      fwrite(dt, file)
    }
  )

  # ── Download All Figures (ZIP of PNGs) ────────────────────────────────────
  output$download_all_figures <- downloadHandler(
    filename = function() paste0("SPECTRA_figures_", Sys.Date(), ".zip"),
    content = function(zipfile) {
      req(rv$analysis_result)
      dt  <- seg_data()
      pol <- seg_exposure_col()
      exp_col <- pol$col
      exp_label <- pol$label
      exp_unit  <- pol$unit

      tmp <- file.path(tempdir(), paste0("spectra_figs_", as.integer(Sys.time())))
      dir.create(tmp, showWarnings = FALSE)

      save_fig <- function(p, name) {
        tryCatch(
          ggplot2::ggsave(file.path(tmp, name), plot = p,
                          width = 10, height = 6, dpi = 150, bg = "white"),
          error = function(e) NULL
        )
      }

      clrs_io <- c("Indoor" = "#2196F3", "Outdoor" = "#FF9800")

      # 1. Indoor/Outdoor distribution bar
      if ("corrected_class" %in% names(dt)) {
        agg1 <- dt[, .N, by = corrected_class]
        agg1[, pct := round(N / sum(N) * 100, 1)]
        p1 <- ggplot2::ggplot(agg1, ggplot2::aes(corrected_class, N, fill = corrected_class)) +
          ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
          ggplot2::geom_text(ggplot2::aes(label = paste0(N, "\n(", pct, "%)")),
                             vjust = -0.3, size = 3.5) +
          ggplot2::scale_fill_manual(values = clrs_io) +
          ggplot2::labs(title = "Indoor / Outdoor Distribution",
                        x = NULL, y = "Count") +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p1, "01_indoor_outdoor_distribution.png")
      }

      # 2. Daily indoor/outdoor % time series
      if (all(c("corrected_class", "date") %in% names(dt))) {
        agg2 <- dt[, .(n = .N), by = .(date = as.Date(date), corrected_class)]
        agg2[, pct := n / sum(n) * 100, by = date]
        p2 <- ggplot2::ggplot(agg2, ggplot2::aes(date, pct, color = corrected_class, group = corrected_class)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::scale_color_manual(values = clrs_io) +
          ggplot2::scale_y_continuous(limits = c(0, 100)) +
          ggplot2::labs(title = "Daily Indoor / Outdoor Percentage",
                        x = "Date", y = "% of Points", color = NULL) +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p2, "02_daily_io_timeseries.png")
      }

      # 3. Exposure by temporal theme
      if (all(c("temporal_theme", "corrected_class", exp_col) %in% names(dt))) {
        agg3 <- dt[!is.na(get(exp_col)),
                   .(mean_exp = mean(get(exp_col), na.rm = TRUE)),
                   by = .(temporal_theme, corrected_class)]
        lvls3 <- c("Night", "Morning", "Afternoon", "Evening")
        agg3[, temporal_theme := factor(temporal_theme, levels = intersect(lvls3, temporal_theme))]
        p3 <- ggplot2::ggplot(agg3, ggplot2::aes(temporal_theme, mean_exp, fill = corrected_class)) +
          ggplot2::geom_col(position = "dodge", width = 0.7) +
          ggplot2::scale_fill_manual(values = clrs_io) +
          ggplot2::labs(title = paste("Mean", exp_label, "by Time of Day"),
                        x = "Period", y = paste0("Mean ", exp_label, " (", exp_unit, ")"),
                        fill = NULL) +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p3, "03_exposure_temporal_theme.png")
      }

      # 4. Exposure by mobility theme
      if (all(c("mobility_theme", "corrected_class", exp_col) %in% names(dt))) {
        agg4 <- dt[!is.na(get(exp_col)),
                   .(mean_exp = mean(get(exp_col), na.rm = TRUE),
                     sd_exp   = sd(get(exp_col),   na.rm = TRUE)),
                   by = .(mobility_theme, corrected_class)]
        p4 <- ggplot2::ggplot(agg4, ggplot2::aes(mobility_theme, mean_exp, fill = corrected_class)) +
          ggplot2::geom_col(position = "dodge", width = 0.7) +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(mean_exp - sd_exp, 0),
                                              ymax = mean_exp + sd_exp),
                                 position = ggplot2::position_dodge(0.7), width = 0.2) +
          ggplot2::scale_fill_manual(values = clrs_io) +
          ggplot2::labs(title = paste("Mean", exp_label, "by Mobility Theme"),
                        x = "Mobility", y = paste0("Mean ", exp_label, " (", exp_unit, ")"),
                        fill = NULL) +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p4, "04_exposure_mobility_theme.png")
      }

      # 5. Weekday × time-of-day heatmap
      if (all(c("weekday_name", "temporal_theme", exp_col) %in% names(dt))) {
        agg5 <- dt[!is.na(get(exp_col)),
                   .(mean_exp = mean(get(exp_col), na.rm = TRUE)),
                   by = .(weekday_name, temporal_theme)]
        day_order <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
        tod_order <- c("Night", "Morning", "Afternoon", "Evening")
        agg5[, weekday_name  := factor(weekday_name,  levels = intersect(day_order, weekday_name))]
        agg5[, temporal_theme := factor(temporal_theme, levels = intersect(tod_order, temporal_theme))]
        p5 <- ggplot2::ggplot(agg5, ggplot2::aes(temporal_theme, weekday_name, fill = mean_exp)) +
          ggplot2::geom_tile(color = "white") +
          ggplot2::scale_fill_gradientn(colors = c("#ffffb2","#fecc5c","#fd8d3c","#e31a1c"),
                                        name = exp_unit) +
          ggplot2::labs(title = paste("Mean", exp_label, "— Weekday × Time of Day"),
                        x = "Time of Day", y = "Weekday") +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p5, "05_weekday_tod_heatmap.png")
      }

      # 6. Exposure by duration theme
      if (all(c("duration_theme", "corrected_class", exp_col) %in% names(dt))) {
        agg6 <- dt[!is.na(get(exp_col)),
                   .(mean_exp = mean(get(exp_col), na.rm = TRUE)),
                   by = .(duration_theme, corrected_class)]
        p6 <- ggplot2::ggplot(agg6, ggplot2::aes(duration_theme, mean_exp, fill = corrected_class)) +
          ggplot2::geom_col(position = "dodge", width = 0.7) +
          ggplot2::scale_fill_manual(values = clrs_io) +
          ggplot2::labs(title = paste("Mean", exp_label, "by Duration Theme"),
                        x = "Duration", y = paste0("Mean ", exp_label, " (", exp_unit, ")"),
                        fill = NULL) +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p6, "06_exposure_duration_theme.png")
      }

      # 7. 24-hour exposure profile
      if (all(c("hour_local", "corrected_class", exp_col) %in% names(dt))) {
        agg7 <- dt[!is.na(get(exp_col)),
                   .(mean_exp = mean(get(exp_col), na.rm = TRUE)),
                   by = .(hour_local, corrected_class)]
        p7 <- ggplot2::ggplot(agg7, ggplot2::aes(hour_local, mean_exp, color = corrected_class, group = corrected_class)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::geom_point(size = 2) +
          ggplot2::scale_color_manual(values = clrs_io) +
          ggplot2::scale_x_continuous(breaks = seq(0, 23, 3)) +
          ggplot2::labs(title = paste("24-Hour", exp_label, "Exposure Profile"),
                        x = "Hour of Day", y = paste0("Mean ", exp_label, " (", exp_unit, ")"),
                        color = NULL) +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p7, "07_hourly_exposure_profile.png")
      }

      # 8. Peak exposure hours
      if (all(c("hour_local", exp_col) %in% names(dt))) {
        agg8 <- dt[!is.na(get(exp_col)),
                   .(mean_exp = mean(get(exp_col), na.rm = TRUE)),
                   by = hour_local][order(-mean_exp)]
        top5 <- head(agg8$hour_local, 5)
        agg8[, is_peak := hour_local %in% top5]
        agg8 <- agg8[order(hour_local)]
        p8 <- ggplot2::ggplot(agg8, ggplot2::aes(factor(hour_local), mean_exp, fill = is_peak)) +
          ggplot2::geom_col(width = 0.8, show.legend = FALSE) +
          ggplot2::scale_fill_manual(values = c("FALSE" = "#90caf9", "TRUE" = "#e53935")) +
          ggplot2::labs(title = paste("Peak", exp_label, "Exposure Hours (Top 5 in Red)"),
                        x = "Hour of Day", y = paste0("Mean ", exp_label, " (", exp_unit, ")")) +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p8, "08_peak_exposure_hours.png")
      }

      # 9. Climate exposure (only if data present)
      if (all(c("climate_theme", "corrected_class", exp_col) %in% names(dt))) {
        agg9 <- dt[climate_theme != "Unknown" & !is.na(get(exp_col)),
                   .(mean_exp = mean(get(exp_col), na.rm = TRUE)),
                   by = .(climate_theme, corrected_class)]
        if (nrow(agg9) > 0) {
          p9 <- ggplot2::ggplot(agg9, ggplot2::aes(climate_theme, mean_exp, fill = corrected_class)) +
            ggplot2::geom_col(position = "dodge", width = 0.7) +
            ggplot2::scale_fill_manual(values = clrs_io) +
            ggplot2::labs(title = paste("Mean", exp_label, "by Climate Condition"),
                          x = "Climate Theme", y = paste0("Mean ", exp_label, " (", exp_unit, ")"),
                          fill = NULL) +
            ggplot2::theme_minimal(base_size = 13) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -30, hjust = 0))
          save_fig(p9, "09_climate_exposure.png")
        }
      }

      # 10. Speed vs OSM confusion matrix heatmap
      if (all(c("speed_class", "osm_label") %in% names(dt))) {
        dt_cm <- dt[!is.na(osm_label)]
        if (nrow(dt_cm) > 0) {
          cm10 <- as.data.frame(table(Speed = dt_cm$speed_class, OSM = dt_cm$osm_label))
          p10 <- ggplot2::ggplot(cm10, ggplot2::aes(OSM, Speed, fill = Freq)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::geom_text(ggplot2::aes(label = Freq), size = 5) +
            ggplot2::scale_fill_gradient(low = "#e3f2fd", high = "#1565c0", name = "Count") +
            ggplot2::labs(title = "Speed vs OSM Classification Confusion Matrix",
                          x = "OSM Building Label", y = "Speed-Based Classification") +
            ggplot2::theme_minimal(base_size = 13)
          save_fig(p10, "10_confusion_matrix_heatmap.png")
        }
      }

      # 11. Dwell time distribution histogram
      if ("dwell_time_min" %in% names(dt)) {
        dt_dwell <- dt[dwell_time_min > 0 & dwell_time_min < quantile(dwell_time_min, 0.99, na.rm = TRUE)]
        p11 <- ggplot2::ggplot(dt_dwell, ggplot2::aes(dwell_time_min)) +
          ggplot2::geom_histogram(bins = 40, fill = "#42a5f5", color = "white") +
          ggplot2::labs(title = "Dwell Time Distribution",
                        x = "Dwell Time (min)", y = "Count") +
          ggplot2::theme_minimal(base_size = 13)
        save_fig(p11, "11_dwell_time_histogram.png")
      }

      # 12. Indoor % pie chart
      if ("corrected_class" %in% names(dt)) {
        agg12 <- dt[, .N, by = corrected_class]
        agg12[, pct := round(N / sum(N) * 100, 1)]
        agg12[, label := paste0(corrected_class, "\n", pct, "%")]
        p12 <- ggplot2::ggplot(agg12, ggplot2::aes(x = "", y = N, fill = corrected_class)) +
          ggplot2::geom_bar(stat = "identity", width = 1) +
          ggplot2::coord_polar("y") +
          ggplot2::scale_fill_manual(values = clrs_io) +
          ggplot2::geom_text(ggplot2::aes(label = label),
                             position = ggplot2::position_stack(vjust = 0.5), size = 4) +
          ggplot2::labs(title = "Indoor / Outdoor Split (%)") +
          ggplot2::theme_void(base_size = 13) +
          ggplot2::theme(legend.position = "none")
        save_fig(p12, "12_indoor_outdoor_pie.png")
      }

      # Zip all PNGs
      png_files <- list.files(tmp, pattern = "\\.png$", full.names = TRUE)
      zip(zipfile, files = png_files, flags = "-j")
    }
  )

  # ── Pre-compute all analytical outputs (no lazy rendering) ─────────────────
  for (.oid in c(
    "plot_dwell_hist", "plot_indoor_pie", "plot_dwell_bar", "dwell_per_person_table",
    "confusion_table", "plot_confusion_matrix", "plot_per_person_accuracy",
    "qa_per_person_table", "qa_overall_ui",
    "seg_validation_stats_table", "seg_confusion_heatmap", "seg_misclass_pie",
    "seg_method_comparison", "seg_indoor_outdoor_bar", "seg_timeseries",
    "seg_outdoor_movement", "seg_outdoor_context", "seg_outdoor_stress", "seg_outdoor_table",
    "seg_temporal_pattern", "seg_mobility_theme", "seg_weekday_heatmap",
    "seg_duration_theme", "seg_hourly_profile", "seg_peak_hours",
    "seg_climate_exposure", "seg_climate_pie",
    "table_classified", "table_grid", "table_daily", "table_participant_exposure"
  )) {
    outputOptions(output, .oid, suspendWhenHidden = FALSE)
  }
  rm(.oid)

}

# ============================================================================
# Launch
# ============================================================================

shinyApp(ui = ui, server = server)
