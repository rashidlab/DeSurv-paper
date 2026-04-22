#!/usr/bin/env Rscript
# code/00_helpers_gridcv.R — Helpers for Grid CV pipeline
#
# Sources the standard 00_helpers.R, then overrides RESULTS_DIR
# to use a "gridcv" subfolder.

source("code/00_helpers.R")

# Override results directory to gridcv subfolder
RESULTS_DIR  <- file.path(RESULTS_BASE, "gridcv")
FIGURE_DIR   <- file.path("paper/figures", "gridcv")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# Update precomputed_path to use new RESULTS_DIR
precomputed_path <- function(name) {
  file.path(RESULTS_DIR, paste0(name, ".rds"))
}

# load_precomputed: check gridcv subfolder first, fall back to base for shared files
load_precomputed <- function(name) {
  path <- precomputed_path(name)
  if (file.exists(path)) {
    return(readRDS(path))
  }
  # Allow fallback only for shared step-2 outputs
  shared_base_files <- c("tar_data_tcgacptac", "sim_figs_by_scenario")
  if (name %in% shared_base_files) {
    base_path <- file.path(RESULTS_BASE, paste0(name, ".rds"))
    if (file.exists(base_path)) {
      message("  Loading from base directory: ", name)
      return(readRDS(base_path))
    }
  }
  stop("Required result not found: ", path,
       "\nRun earlier pipeline steps first.")
}

message("  Grid CV results dir: ", RESULTS_DIR)
