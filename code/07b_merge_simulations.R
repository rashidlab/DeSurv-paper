#!/usr/bin/env Rscript
# code/07b_merge_simulations.R — Merge partial simulation results
#
# Reads all results/precomputed/sim_partial_*.rds files produced by SLURM
# array jobs, combines them, and produces the final sim_figs_by_scenario.rds
# that downstream scripts (08_figures.R, paper) expect.
#
# Usage:
#   Rscript code/07b_merge_simulations.R
#
# Inputs:  results/precomputed/sim_partial_*.rds (6 files from array jobs)
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
source("sim_figs.R")
source("code/sim_helpers.R")

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

message(sprintf("  Found %d partial result files:", length(partial_files)))
for (f in partial_files) message("    ", basename(f))

# ── Expected partials ─────────────────────────────────────────────────────
expected_scenarios <- c("R0_easy", "R00_null", "R_mixed")
expected_analyses <- c("bo_tune_ntop", "bo_tune_ntop_alpha0")
expected_names <- outer(expected_scenarios, expected_analyses,
                        function(s, a) paste0("sim_partial_", s, "_", a, ".rds"))
missing <- setdiff(as.vector(expected_names), basename(partial_files))
if (length(missing)) {
  warning("Missing partial files (results will be incomplete): ",
          paste(missing, collapse = ", "))
}

# ── Load and combine ──────────────────────────────────────────────────────
message("  Loading and combining partial results...")
partial_tables <- lapply(partial_files, readRDS)
combined_table <- dplyr::bind_rows(partial_tables)
message(sprintf("  Combined: %d rows", nrow(combined_table)))

# ── Build figures by scenario ─────────────────────────────────────────────
message("  Building simulation figures...")
sim_figs_by_scenario <- build_sim_figs_by_scenario(combined_table)

# ── Save final result ─────────────────────────────────────────────────────
out_path <- precomputed_path("sim_figs_by_scenario")
saveRDS(sim_figs_by_scenario, out_path)
message(sprintf("  Saved: %s", out_path))

message("=== Step 7b complete ===")
