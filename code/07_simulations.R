#!/usr/bin/env Rscript
# code/07_simulations.R — Simulation studies
#
# Two execution modes:
#
#   Sequential (local/quick):
#     Rscript code/07_simulations.R
#     Runs 3 scenarios x 100 replicates x 2 methods = 600 evaluations.
#     Outputs: results/precomputed/sim_figs_by_scenario.rds
#
#   SLURM array (HPC, one task per scenario×method pair):
#     sbatch slurm/run_simulations.sh
#     Env vars: DESURV_SIM_SCENARIO, DESURV_SIM_ANALYSIS
#     Outputs: results/precomputed/sim_partial_<scenario>_<analysis>.rds
#     Then run: Rscript code/07b_merge_simulations.R to combine.
#
# Scenarios:
#   R0_easy  — prognostic programs explain low variance (primary)
#   R00_null — no survival signal (null control)
#   R_mixed  — partial overlap between variance and prognosis
#
# Analysis methods:
#   bo_tune_ntop        — DeSurv with BO-tuned alpha + ntop
#   bo_tune_ntop_alpha0 — NMF with BO-tuned ntop
#
# Inputs:  R/simulation_functions/*.R, code/sim_helpers.R, sim_figs.R
# Runtime: ~6-12 hours on HPC (30 cores per analysis), ~10-15 min in quick mode

message("=== Step 7: Simulation Studies ===")
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

# ── SLURM array mode detection ────────────────────────────────────────────
SIM_SCENARIO_FILTER <- Sys.getenv("DESURV_SIM_SCENARIO", "")
SIM_ANALYSIS_FILTER <- Sys.getenv("DESURV_SIM_ANALYSIS", "")
SIM_ARRAY_MODE <- nzchar(SIM_SCENARIO_FILTER) && nzchar(SIM_ANALYSIS_FILTER)

if (SIM_ARRAY_MODE) {
  message(sprintf("  SLURM array mode: scenario=%s, analysis=%s",
                  SIM_SCENARIO_FILTER, SIM_ANALYSIS_FILTER))
}

# ── Quick mode overrides ──────────────────────────────────────────────────
if (CONFIG$quick) {
  SIM_DATASETS_PER_SCENARIO <- 2L
  SIM_CV_NSTARTS <<- 2L
  SIM_BO_N_INIT  <<- 5L
  SIM_BO_N_ITER  <<- 5L
  SIM_BO_CANDIDATE_POOL <<- 200L
  message("  Quick mode: 2 replicates, 5 BO iterations, 2 init starts")
}

# ── Scenarios ─────────────────────────────────────────────────────────────
SIMULATION_SCENARIOS <- list(
  list(scenario_id = "R0_easy", scenario = "R0",
       description = "Default easy/sanity scenario",
       replicates = SIM_DATASETS_PER_SCENARIO,
       seed_offset = SIM_GLOBAL_SEED, overrides = list()),
  list(scenario_id = "R00_null", scenario = "R00",
       description = "No survival associated programs",
       replicates = SIM_DATASETS_PER_SCENARIO,
       seed_offset = SIM_GLOBAL_SEED, overrides = list()),
  list(scenario_id = "R_mixed", scenario = "R_mixed",
       description = "Mixed marker/background survival genes",
       replicates = SIM_DATASETS_PER_SCENARIO,
       seed_offset = SIM_GLOBAL_SEED + 3000L, overrides = list())
)

if (SIM_ARRAY_MODE) {
  SIMULATION_SCENARIOS <- Filter(
    function(s) s$scenario_id == SIM_SCENARIO_FILTER,
    SIMULATION_SCENARIOS
  )
  if (!length(SIMULATION_SCENARIOS))
    stop("Unknown DESURV_SIM_SCENARIO=", SIM_SCENARIO_FILTER,
         "\nValid: R0_easy, R00_null, R_mixed")
}

# ── Analysis specs ────────────────────────────────────────────────────────
SIM_ANALYSIS_SPECS <- list(
  list(analysis_id = "bo_tune_ntop", label = "DeSurv BO tuning ntop",
       mode = "bayesopt",
       bounds = modifyList(SIM_DESURV_BO_BOUNDS,
                           list(ntop = list(lower = 50, upper = 250, type = "integer"))),
       bo_fixed = list(ngene = NULL, lambdaW_grid = 0, lambdaH_grid = 0)),
  list(analysis_id = "bo_tune_ntop_alpha0", label = "NMF BO tuning ntop",
       mode = "bayesopt",
       bounds = modifyList(SIM_DESURV_BO_BOUNDS,
                           list(alpha_grid = list(lower = 0, upper = 0, type = "continuous"),
                                ntop = list(lower = 50, upper = 250, type = "integer"))),
       bo_fixed = list(ngene = NULL, alpha_grid = 0, lambdaW_grid = 0, lambdaH_grid = 0),
       final_overrides = list(alpha = 0))
)

if (SIM_ARRAY_MODE) {
  SIM_ANALYSIS_SPECS <- Filter(
    function(a) a$analysis_id == SIM_ANALYSIS_FILTER,
    SIM_ANALYSIS_SPECS
  )
  if (!length(SIM_ANALYSIS_SPECS))
    stop("Unknown DESURV_SIM_ANALYSIS=", SIM_ANALYSIS_FILTER,
         "\nValid: bo_tune_ntop, bo_tune_ntop_alpha0")
}

# ── Core simulation runner (used by both modes) ───────────────────────────
run_sim_pipeline <- function(scenarios, analysis_specs) {
  dataset_specs <- build_simulation_dataset_specs(
    scenarios = scenarios,
    replicates_per_scenario = SIM_DATASETS_PER_SCENARIO,
    base_seed = SIM_GLOBAL_SEED
  )
  message(sprintf("  %d dataset specs across %d scenario(s)",
                  length(dataset_specs), length(scenarios)))

  datasets <- lapply(dataset_specs, generate_simulation_dataset)

  n_total <- length(datasets) * length(analysis_specs)
  message(sprintf("  Running %d analyses (%d datasets x %d methods)...",
                  n_total, length(datasets), length(analysis_specs)))

  results <- vector("list", n_total)
  idx <- 0L
  for (i in seq_along(datasets)) {
    for (j in seq_along(analysis_specs)) {
      idx <- idx + 1L
      ds <- datasets[[i]]
      spec <- analysis_specs[[j]]
      message(sprintf("  [%d/%d] %s rep %d - %s",
                      idx, n_total,
                      ds$metadata$scenario_id,
                      ds$metadata$replicate,
                      spec$analysis_id))
      results[[idx]] <- tryCatch(
        run_simulation_analysis(ds, spec),
        error = function(e) {
          warning(sprintf("  FAILED: %s rep %d %s: %s",
                          ds$metadata$scenario_id, ds$metadata$replicate,
                          spec$analysis_id, conditionMessage(e)))
          NULL
        }
      )
    }
  }

  message("  Summarizing simulation results...")
  summarize_simulation_results(results)
}

# ── Execute ───────────────────────────────────────────────────────────────
if (SIM_ARRAY_MODE) {
  partial_name <- sprintf("sim_partial_%s_%s", SIM_SCENARIO_FILTER, SIM_ANALYSIS_FILTER)
  cache_or_compute(partial_name, {
    run_sim_pipeline(SIMULATION_SCENARIOS, SIM_ANALYSIS_SPECS)
  })
} else {
  sim_figs_by_scenario <- cache_or_compute("sim_figs_by_scenario", {
    sim_results_table <- run_sim_pipeline(SIMULATION_SCENARIOS, SIM_ANALYSIS_SPECS)
    message("  Building simulation figures...")
    build_sim_figs_by_scenario(sim_results_table)
  })
}

message("=== Step 7 complete ===")
