#!/usr/bin/env Rscript
# code/00_helpers_loco_tol.R — Helpers for LOCO + tolerance-rule pipeline
#
# Writes new results to results/precomputed/loco_tol/.
# Falls back to results/precomputed/loco/ for BO results and other
# upstream files that are unchanged by the tolerance rule.

source("code/00_helpers.R")

# Override results directory to loco_tol subfolder
RESULTS_DIR  <- file.path(RESULTS_BASE, "loco_tol")
FIGURE_DIR   <- file.path("paper/figures", "loco_tol")
LOCO_DIR     <- file.path(RESULTS_BASE, "loco")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# Update precomputed_path to use new RESULTS_DIR
precomputed_path <- function(name) {
  file.path(RESULTS_DIR, paste0(name, ".rds"))
}

# Files that can be loaded from the loco/ directory (BO and upstream outputs)
LOCO_SHARED_FILES <- c(
  "tar_data_tcgacptac",
  "desurv_bo_results_tcgacptac",
  "desurv_bo_results_alpha0_tcgacptac",
  "desurv_bo_results_elbowk_tcgacptac",
  "fit_std_tcgacptac",
  "std_nmf_selected_k_tcgacptac"
)

# load_precomputed: check loco_tol first, fall back to loco for shared files
load_precomputed <- function(name) {
  path <- precomputed_path(name)
  if (file.exists(path)) {
    return(readRDS(path))
  }
  # Fall back to loco/ for BO results and upstream data
  if (name %in% LOCO_SHARED_FILES) {
    loco_path <- file.path(LOCO_DIR, paste0(name, ".rds"))
    if (file.exists(loco_path)) {
      message("  Loading from loco/: ", name)
      return(readRDS(loco_path))
    }
  }
  # Fall back to base for shared step-2 outputs
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

message("  LOCO-tol results dir: ", RESULTS_DIR)
message("  Falling back to loco/ for: ", paste(LOCO_SHARED_FILES, collapse = ", "))
