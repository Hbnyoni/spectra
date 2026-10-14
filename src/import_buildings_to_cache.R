#!/usr/bin/env Rscript
# ===========================================================================
# SPECTRA: Spatiotemporal Personal Exposure Characterization and TRAjectory
# Analyzer.
# Building footprint ingestion into the local spatial cache.
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
# Import buildings.csv (102GB) into buildings_cache.sqlite (one-time operation)
# ============================================================================
# Optimized for SPEED: 250K batches, WAL mode, mmap, direct fread

library(data.table)
library(RSQLite)
library(sf)

source("trajectory_analysis_core.R")

BUILDINGS_CSV_PATH <- "buildings.csv"
CACHE_DB_PATH      <- "buildings_cache.sqlite"
BATCH_SIZE         <- 250000  # 5x larger (fewer DB round-trips)
BUILDING_TILE_SIZE <- 0.25

cat("\n")
cat(strrep("=", 80), "\n")
cat("BUILDINGS CACHE IMPORT: CSV → SQLite Database (OPTIMIZED)\n")
cat(strrep("=", 80), "\n")

# Step 1: Initialize cache database
cat("\n[1/4] Initializing database...\n")
init_building_cache(CACHE_DB_PATH)
cat("     Database initialized at:", CACHE_DB_PATH, "\n")

# Step 2: Connect and prepare
cat("\n[2/4] Connecting to database...\n")
con <- RSQLite::dbConnect(RSQLite::SQLite(), CACHE_DB_PATH)
on.exit(RSQLite::dbDisconnect(con))

# Ultra-fast SQLite settings for bulk import
RSQLite::dbExecute(con, "PRAGMA synchronous = OFF")
RSQLite::dbExecute(con, "PRAGMA cache_size = 100000")
RSQLite::dbExecute(con, "PRAGMA temp_store = MEMORY")
RSQLite::dbExecute(con, "PRAGMA journal_mode = WAL")
RSQLite::dbExecute(con, "PRAGMA mmap_size = 30000000000")

chunk_num <- 0
total_imported <- 0
tiles_to_insert <- data.table()

cat("\n[3/4] Reading and importing buildings.csv (437M rows)...\n")

tryCatch({
  # Stream CSV instead of loading entire file
  cat("     Streaming CSV (NOT loading full file into memory)...\n\n")
  
  # Open file connection
  con_file <- file(BUILDINGS_CSV_PATH, "r")
  on.exit(close(con_file), add = TRUE)
  
  # Skip header
  header <- readLines(con_file, n = 1)
  
  batch_num <- 0
  repeat {
    # Read exactly BATCH_SIZE lines
    lines <- readLines(con_file, n = BATCH_SIZE)
    if (length(lines) == 0) break
    
    batch_num <- batch_num + 1
    
    # Parse this batch
    chunk <- data.table::fread(text = c(header, lines), sep = ",", 
                               select = c("latitude", "longitude", "geometry", 
                                        "area_in_meters", "confidence"))
    
    if (nrow(chunk) == 0) break
    
    # Compute tile IDs (vectorized)
    chunk$tile_id <- sprintf("%.2f_%.2f",
      floor(chunk$latitude / BUILDING_TILE_SIZE) * BUILDING_TILE_SIZE,
      floor(chunk$longitude / BUILDING_TILE_SIZE) * BUILDING_TILE_SIZE)
    
    # Prepare data frame
    chunk_df <- data.frame(
      tile_id = chunk$tile_id,
      latitude = chunk$latitude,
      longitude = chunk$longitude,
      geometry_wkt = chunk$geometry,
      area_m2 = chunk$area_in_meters,
      confidence = chunk$confidence,
      osm_id = NA_character_
    )
    
    # Bulk insert
    RSQLite::dbAppendTable(con, "buildings", chunk_df)
    
    # Collect tiles (deduped later)
    tiles_to_insert <- rbind(tiles_to_insert,
      unique(data.table(tile_id = chunk$tile_id, 
                        lat = chunk$latitude, 
                        lon = chunk$longitude)))
    
    total_imported <- total_imported + nrow(chunk_df)
    
    # Progress every 10 batches (250K rows)
    if (batch_num %% 10 == 0) {
      cat(sprintf("     [Batch %d] Imported %s rows\n", batch_num, 
                  format(total_imported, big.mark = ",")))
    }
  }
  
  cat("\n     Total imported:", format(total_imported, big.mark = ","), "buildings\n")
  
  # Batch insert all unique tiles
  cat("\n     Updating tile metadata...\n")
  if (nrow(tiles_to_insert) > 0) {
    tiles_to_insert <- unique(tiles_to_insert)
    for (j in 1:nrow(tiles_to_insert)) {
      RSQLite::dbExecute(con,
        "INSERT OR IGNORE INTO cache_tiles 
         (tile_id, lat_min, lat_max, lon_min, lon_max, fetched_at, n_buildings, status)
         VALUES (?, ?, ?, ?, ?, datetime('now'),
                 (SELECT COUNT(*) FROM buildings WHERE tile_id = ?), 'cached')",
        list(tiles_to_insert$tile_id[j],
             tiles_to_insert$lat[j],
             tiles_to_insert$lat[j] + BUILDING_TILE_SIZE,
             tiles_to_insert$lon[j],
             tiles_to_insert$lon[j] + BUILDING_TILE_SIZE,
             tiles_to_insert$tile_id[j]))
    }
  }
  
}, error = function(e) {
  cat("\n[ERROR]", e$message, "\n")
  quit(status = 1)
})

# Re-enable normal synchronous mode
RSQLite::dbExecute(con, "PRAGMA synchronous = NORMAL")

# Step 4: Create indices
cat("\n[4/4] Creating indices for fast queries...\n")

tryCatch({
  indices <- c(
    "CREATE INDEX IF NOT EXISTS idx_bldg_tile ON buildings(tile_id)",
    "CREATE INDEX IF NOT EXISTS idx_bldg_latlon ON buildings(latitude, longitude)",
    "CREATE INDEX IF NOT EXISTS idx_cache_tile ON cache_tiles(tile_id)"
  )
  
  RSQLite::dbBegin(con)
  for (idx in indices) {
    RSQLite::dbExecute(con, idx)
  }
  RSQLite::dbCommit(con)
  
  # Vacuum and optimize
  RSQLite::dbExecute(con, "VACUUM")
  RSQLite::dbExecute(con, "ANALYZE")
  
  # Get final stats
  stats <- RSQLite::dbGetQuery(con, 
    "SELECT 
      (SELECT COUNT(*) FROM buildings) as total_buildings,
      (SELECT COUNT(DISTINCT tile_id) FROM cache_tiles) as total_tiles,
      (SELECT COUNT(*) FROM cache_tiles WHERE status = 'cached') as cached_tiles
    ")
  
  cat(strrep("=", 80), "\n")
  cat("✓ IMPORT COMPLETE!\n")
  cat("\n")
  cat("Database Summary:\n")
  cat("  Total buildings:", format(stats$total_buildings, big.mark = ","), "\n")
  cat("  Total tiles:", format(stats$total_tiles, big.mark = ","), "\n")
  cat("  Cached tiles:", format(stats$cached_tiles, big.mark = ","), "\n")
  cat("  Database size:",
      format(file.info(CACHE_DB_PATH)$size, units = "GB"), "\n")
  cat("\n")
  cat("You can now DELETE buildings.csv (102GB) — it's no longer needed!\n")
  cat("Update app.R to use get_buildings_from_cache() instead of reading CSV.\n")
  cat(strrep("=", 80), "\n")
  cat("\n")
  
}, error = function(e) {
  cat("\n[ERROR during indexing]", e$message, "\n")
  quit(status = 1)
})

cat("✓ Script completed successfully.\n")
