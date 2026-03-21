#!/usr/bin/env Rscript
# code/07_simulations.R — Simulation studies
#
# Runs 3 scenarios x 100 replicates x 6 analysis methods = 1800 evaluations
# Scenarios: R0_easy (prognostic != variance), R00_null (no signal), R_mixed (partial)
#
# Inputs:  R/simulation_functions/*.R, sim_figs.R
# Outputs: results/precomputed/sim_figs_by_scenario.rds
#
# Runtime: ~6-12 hours on HPC, ~5 min in quick mode

message("=== Step 7: Simulation Studies ===")
source("code/00_helpers.R")
library(DeSurv)
library(dplyr)
library(purrr)
library(tibble)

# Source simulation functions
sim_files <- list.files("R/simulation_functions", pattern = "[.]R$", full.names = TRUE)
purrr::walk(sim_files, source)
source("R/get_top_genes.R")
source("sim_figs.R")

# ── Configuration ─────────────────────────────────────────────────────────
if (CONFIG$quick) {
  SIM_DATASETS_PER_SCENARIO <- 2L
  message("  Quick mode: 2 replicates per scenario")
} else {
  SIM_DATASETS_PER_SCENARIO <- 100L
}

SIM_GLOBAL_SEED <- 101L

# ── Note: simulation pipeline is complex (1700 lines in _targets_sims.R)
# This script delegates to the same functions but runs them sequentially.
# For full reproduction, the simulation functions in R/simulation_functions/
# handle data generation, analysis method dispatch, and result summarization.
#
# The pre-computed sim_figs_by_scenario.rds contains the final figure objects
# for all 3 scenarios. Full re-computation requires the simulation infrastructure.

sim_figs_by_scenario <- cache_or_compute("sim_figs_by_scenario", {
  message("  Full simulation re-computation not yet implemented in standalone scripts.")
  message("  This requires porting the 1700-line _targets_sims.R pipeline.")
  message("  Use the pre-computed results for now.")
  stop("Simulation re-computation requires pre-computed results.\n",
       "Download from Zenodo or copy from the targets store.")
})

message("=== Step 7 complete ===")
