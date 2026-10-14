# ===========================================================================
# SPECTRA: Spatiotemporal Personal Exposure Characterization and TRAjectory
# Analyzer.
# Dependency installation helper.
#
# Copyright (c) 2026 Place Alert Labs, Midlands State University.
# All rights reserved. See LICENSE and NOTICE in the repository root.
#
# Author and lead developer: Handsome Bongani Nyoni.
#
# Proprietary and confidential. No part of this file may be used, copied,
# modified or distributed without the written permission of the copyright
# holder. Academic use by request.
# ===========================================================================
# ============================================================================
# CHEAQI-MNCH Trajectory Exposure Analysis - R Shiny App
# Package Installation Script
# ============================================================================

cat("Installing required packages for CHEAQI-MNCH Trajectory Exposure R Shiny App...\n\n")

# Core packages
required_packages <- c(
  # Shiny & UI
  "shiny",
  "shinydashboard",
  "shinycssloaders",
  "shinyWidgets",
  "shinyjs",
  "DT",
  

  # Spatial
  "sf",
  "osmdata",
  
  # Mapping
  "leaflet",
  "leaflet.extras",
  
  # Plotting
  "plotly",
  "ggplot2",
  "viridis",
  
  # Data manipulation
  "data.table",
  "dplyr",
  "tidyr",
  "lubridate",
  
  # Utils
  "geosphere",
  "htmltools",
  "htmlwidgets",
  "jsonlite",
  "readr",
  "scales"
)

# Install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("  Installing: ", pkg, "\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(paste0("  Already installed: ", pkg, "\n"))
  }
}

cat("Checking and installing packages...\n")
for (pkg in required_packages) {
  tryCatch(
    install_if_missing(pkg),
    error = function(e) {
      cat(paste0("  WARNING: Failed to install ", pkg, ": ", e$message, "\n"))
    }
  )
}

cat("\n============================================\n")
cat("Installation complete!\n")
cat("Run the app with: shiny::runApp('shiny_app')\n")
cat("Or open app.R in RStudio and click 'Run App'\n")
cat("============================================\n")
