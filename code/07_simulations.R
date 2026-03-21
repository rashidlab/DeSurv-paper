#!/usr/bin/env Rscript
# code/07_simulations.R — Simulation studies
#
# Runs 3 scenarios x 100 replicates x 4 analysis methods = 1200 evaluations
# (quick mode: 2 replicates x 4 methods = 24 evaluations)
#
# Scenarios:
#   R0_easy  — prognostic programs explain low variance (primary)
#   R00_null — no survival signal (null control)
#   R_mixed  — partial overlap between variance and prognosis
#
# Analysis methods per replicate:
#   bo            — DeSurv with BO-tuned alpha
#   bo_alpha0     — Standard NMF (alpha=0) with BO
#   bo_tune_ntop  — DeSurv with BO-tuned alpha + ntop
#   bo_tune_ntop_alpha0 — NMF with BO-tuned ntop
#
# Inputs:  R/simulation_functions/*.R, code/sim_helpers.R, sim_figs.R
# Outputs: results/precomputed/sim_figs_by_scenario.rds
#
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

# Simulation infrastructure (data generation, analysis dispatch, metrics)
source("code/sim_helpers.R")

# ── Configuration ─────────────────────────────────────────────────────────
if (CONFIG$quick) {
  SIM_DATASETS_PER_SCENARIO <- 2L
  SIM_CV_NSTARTS <<- 2L
  SIM_BO_N_INIT  <<- 5L
  SIM_BO_N_ITER  <<- 5L
  SIM_BO_CANDIDATE_POOL <<- 200L
  message("  Quick mode: 2 replicates, 5 BO iterations, 2 init starts")
} else {
  SIM_DATASETS_PER_SCENARIO <- 100L
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

# ── Analysis specs ────────────────────────────────────────────────────────
SIM_ANALYSIS_SPECS <- list(
  list(analysis_id = "bo", label = "Bayesian optimization",
       mode = "bayesopt", bounds = SIM_DESURV_BO_BOUNDS,
       bo_fixed = list(ngene = NULL, ntop = SIM_DEFAULT_NTOP, lambdaW_grid = 0, lambdaH_grid = 0)),
  list(analysis_id = "bo_alpha0", label = "BO with NMF (alpha=0)",
       mode = "bayesopt",
       bounds = modifyList(SIM_DESURV_BO_BOUNDS,
                           list(alpha_grid = list(lower = 0, upper = 0, type = "continuous"))),
       bo_fixed = list(ngene = NULL, ntop = SIM_DEFAULT_NTOP, alpha_grid = 0, lambdaW_grid = 0, lambdaH_grid = 0),
       final_overrides = list(alpha = 0)),
  list(analysis_id = "bo_tune_ntop", label = "BO tuning ntop",
       mode = "bayesopt",
       bounds = modifyList(SIM_DESURV_BO_BOUNDS,
                           list(ntop = list(lower = 50, upper = 250, type = "integer"))),
       bo_fixed = list(ngene = NULL, lambdaW_grid = 0, lambdaH_grid = 0)),
  list(analysis_id = "bo_tune_ntop_alpha0", label = "BO NMF tuning ntop",
       mode = "bayesopt",
       bounds = modifyList(SIM_DESURV_BO_BOUNDS,
                           list(alpha_grid = list(lower = 0, upper = 0, type = "continuous"),
                                ntop = list(lower = 50, upper = 250, type = "integer"))),
       bo_fixed = list(ngene = NULL, alpha_grid = 0, lambdaW_grid = 0, lambdaH_grid = 0),
       final_overrides = list(alpha = 0))
)

# ── Run simulations ───────────────────────────────────────────────────────
sim_figs_by_scenario <- cache_or_compute("sim_figs_by_scenario", {

  # 1. Build dataset specs
  message("  Building simulation dataset specs...")
  dataset_specs <- build_simulation_dataset_specs(
    scenarios = SIMULATION_SCENARIOS,
    replicates_per_scenario = SIM_DATASETS_PER_SCENARIO,
    base_seed = SIM_GLOBAL_SEED
  )
  message(sprintf("  %d dataset specs across %d scenarios",
                  length(dataset_specs), length(SIMULATION_SCENARIOS)))

  # 2. Generate datasets
  message("  Generating simulation datasets...")
  datasets <- lapply(dataset_specs, function(spec) {
    generate_simulation_dataset(spec[[1]])
  })

  # 3. Run analyses (cross product: datasets x analysis_specs)
  n_total <- length(datasets) * length(SIM_ANALYSIS_SPECS)
  message(sprintf("  Running %d analyses (%d datasets x %d methods)...",
                  n_total, length(datasets), length(SIM_ANALYSIS_SPECS)))

  results <- vector("list", n_total)
  idx <- 0L
  for (i in seq_along(datasets)) {
    for (j in seq_along(SIM_ANALYSIS_SPECS)) {
      idx <- idx + 1L
      ds <- datasets[[i]]
      spec <- SIM_ANALYSIS_SPECS[[j]]
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

  # 4. Summarize results
  message("  Summarizing simulation results...")
  sim_results_table <- summarize_simulation_results(results)

  # 5. Build figures by scenario
  message("  Building simulation figures...")
  build_sim_figs_by_scenario(sim_results_table)
})

message("=== Step 7 complete ===")
