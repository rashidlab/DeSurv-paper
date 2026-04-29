#!/usr/bin/env Rscript
# code/06_sensitivity_analysis.R — K-sensitivity grid analysis
#
# Evaluates DeSurv across a grid of k x alpha combinations in external
# validation, with and without adjustment for PurIST/DeCAF classifiers.
# Also computes k=3 vs k=7 factor nesting analysis.
#
# This analysis was originally run via standalone scripts in inst/:
#   inst/cv_grid_validation_analysis_amber_optimal.R → adj_p_270, hcor_270, master_rows_270
#   inst/cv_grid_validation_analysis_allgenes.R      → adj_p_all, hcor_all, master_rows_all
#   inst/k3_k7_nesting_heatmap.R                     → k3_k7_summary
#   inst/create_ksens_rds.R                          → lam300_summary, production_summary
#
# Inputs:  model fits, validation data, cv_grid infrastructure
# Outputs: results/cv_grid/*.rds (10 files)
#
# Runtime: ~2-4 hours on HPC, not supported in quick mode
#
# NOTE: The pre-computed results in results/cv_grid/ are sufficient for
# paper rendering. Full re-computation requires the cv_grid infrastructure
# from the inst/ scripts.

message("=== Step 6: Sensitivity Analysis ===")
source("code/00_helpers.R")

cv_grid_files <- c(
  "adj_p_270_matrix.rds", "adj_p_all_matrix.rds", "unadj_p_270_matrix.rds",
  "hcor_270_matrix.rds", "hcor_all_matrix.rds",
  "master_rows_270.rds", "master_rows_all.rds",
  "k3_k7_summary.rds", "lam300_summary.rds", "production_summary.rds"
)

existing <- file.exists(file.path(CV_GRID_DIR, cv_grid_files))
message(sprintf("  %d/%d cv_grid result files found", sum(existing), length(cv_grid_files)))

if (all(existing)) {
  message("  All sensitivity analysis results present. Skipping computation.")
} else {
  missing <- cv_grid_files[!existing]
  message("  Missing files: ", paste(missing, collapse = ", "))
  if (CONFIG$quick) {
    message("  Sensitivity analysis not supported in quick mode. Skipping.")
  } else {
    message("  To regenerate, run the inst/ scripts:")
    message("    Rscript inst/cv_grid_validation_analysis_amber_optimal.R")
    message("    Rscript inst/cv_grid_validation_analysis_allgenes.R")
    message("    Rscript inst/k3_k7_nesting_heatmap.R")
    message("    Rscript inst/create_ksens_rds.R")
  }
}

message("=== Step 6 complete ===")
