#!/usr/bin/env Rscript
# ===========================================================================
# SPECTRA: Spatiotemporal Personal Exposure Characterization and TRAjectory
# Analyzer.
# Batch segmentation driver.
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

base_dir <- Sys.getenv("SPECTRA_DATA_DIR", "data/Segmentation")

runs <- data.frame(
  country = c("Gambia", "Kenya", "Mozambique"),
  input = file.path(
    base_dir,
    c("Gambia", "Kenya", "Mozambique"),
    c("Gambia_raw.rds", "Kenya_raw.rds", "Mozambique_raw.rds")
  ),
  outdir = file.path(
    base_dir,
    "cli_outputs",
    c("Gambia", "Kenya", "Mozambique")
  ),
  stringsAsFactors = FALSE
)

cell_size <- "100"
fetch_buildings <- "true"
use_buildings <- "true"
building_cache <- file.path(getwd(), "buildings_cache.sqlite")
preloaded_buildings <- file.path(getwd(), "buildings_preloaded.csv")

for (i in seq_len(nrow(runs))) {
  cat("\n============================================================\n")
  cat("Running", runs$country[[i]], "\n")
  cat("Input: ", runs$input[[i]], "\n", sep = "")
  cat("Output:", runs$outdir[[i]], "\n")
  cat("============================================================\n")

  if (!file.exists(runs$input[[i]])) {
    stop("Missing input file: ", runs$input[[i]])
  }

  cmd <- c(
    "run_analysis_cli.R",
    paste0("--input=", runs$input[[i]]),
    paste0("--outdir=", runs$outdir[[i]]),
    "--lat=lat",
    "--lon=lon",
    "--pid=PID",
    "--date=date",
    "--pm25=pm25",
    "--no2=no2",
    "--temp=temp",
    "--humidity=rh",
    "--country=country",
    paste0("--cell-size=", cell_size),
    "--max-dwell=20",
    "--speed-limit=100",
    "--speed-outdoor-threshold=100",
    "--building-buffer=20",
    paste0("--use-buildings=", use_buildings),
    paste0("--fetch-buildings=", fetch_buildings),
    paste0("--building-cache=", building_cache),
    paste0("--preloaded-buildings=", preloaded_buildings)
  )

  status <- system2("Rscript", cmd)
  if (!identical(status, 0L)) {
    stop("Run failed for ", runs$country[[i]], " with status ", status)
  }
}

cat("\nAll batch runs completed.\n")
