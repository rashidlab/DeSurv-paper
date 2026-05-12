#!/usr/bin/env Rscript
# code/07b_merge_simulations.R — Merge partial simulation results
#
# Reads all results/precomputed/sim_partial_*.rds files produced by SLURM
# array jobs, combines them, and produces the final sim_figs_by_scenario.rds
# that downstream scripts (paper/) expect.
#
# Supports both per-replicate files (sim_partial_<scenario>_<analysis>_rep*.rds)
# and legacy per-(scenario,analysis) files (sim_partial_<scenario>_<analysis>.rds).
#
# Usage:
#   Rscript code/07b_merge_simulations.R
#
# Inputs:  results/precomputed/sim_partial_*.rds (600 files from array jobs)
# Outputs: results/precomputed/sim_figs_by_scenario.rds

message("=== Step 7b: Merge Simulation Results ===")
source("code/00_helpers.R")
library(DeSurv)
library(dplyr)
library(purrr)
library(tibble)
library(survival)

# Source simulation helper functions
sim_files <- list.files("R/simulation_functions", pattern = "[.]R$", full.names = TRUE)
purrr::walk(sim_files, source)
source("R/get_top_genes.R")
source("R/bo_helpers.R")
source("R/sim_figs.R")
source("R/sim_helpers.R")

# ── Find partial result files ─────────────────────────────────────────────
partial_files <- sort(list.files(
  RESULTS_DIR,
  pattern = "^sim_partial_.*\\.rds$",
  full.names = TRUE
))

if (!length(partial_files)) {
  stop("No sim_partial_*.rds files found in ", RESULTS_DIR,
       "\nRun the simulation array jobs first.")
}

message(sprintf("  Found %d partial result files", length(partial_files)))

# ── Check for missing replicates ──────────────────────────────────────────
expected_scenarios <- c("R0_easy", "R00_null", "R_mixed")
expected_analyses  <- c("bo_tune_ntop", "bo_tune_ntop_alpha0")

per_rep_pattern <- "^sim_partial_(.+)_(.+)_rep(\\d+)\\.rds$"
is_per_rep <- grepl(per_rep_pattern, basename(partial_files))

if (any(is_per_rep)) {
  message("  Per-replicate mode detected")
  found <- basename(partial_files)[is_per_rep]
  missing <- c()
  for (s in expected_scenarios) {
    for (a in expected_analyses) {
      for (r in seq_len(SIM_DATASETS_PER_SCENARIO)) {
        fname <- sprintf("sim_partial_%s_%s_rep%03d.rds", s, a, r)
        if (!fname %in% found) missing <- c(missing, fname)
      }
    }
  }
  if (length(missing)) {
    warning(sprintf("Missing %d partial file(s) (results will be incomplete). First few: %s",
                    length(missing), paste(head(missing, 3), collapse = ", ")))
  }
} else {
  message("  Legacy per-(scenario,analysis) mode detected")
  expected_names <- as.vector(outer(
    expected_scenarios, expected_analyses,
    function(s, a) paste0("sim_partial_", s, "_", a, ".rds")
  ))
  missing <- setdiff(expected_names, basename(partial_files))
  if (length(missing)) {
    warning("Missing partial files: ", paste(missing, collapse = ", "))
  }
}

# ── Load and combine ──────────────────────────────────────────────────────
message("  Loading and combining partial results...")
partial_tables <- lapply(partial_files, readRDS)
combined_table <- dplyr::bind_rows(partial_tables)
message(sprintf("  Combined: %d rows across %d files", nrow(combined_table), length(partial_files)))

# ── Build figures by scenario ─────────────────────────────────────────────
message("  Building simulation figures...")
sim_figs_by_scenario <- build_sim_figs_by_scenario(combined_table)

# ── Save final result ─────────────────────────────────────────────────────
out_path <- precomputed_path("sim_figs_by_scenario")
saveRDS(sim_figs_by_scenario, out_path)
message(sprintf("  Saved: %s", out_path))

message("=== Step 7b complete ===")
