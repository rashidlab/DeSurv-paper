#!/usr/bin/env Rscript
# code/03_bayesian_optimization.R — Bayesian optimization for hyperparameter selection
#
# Runs BO for:
#   1. DeSurv (supervised, alpha tuned)
#   2. Standard NMF baseline (alpha fixed at 0)
#   3. Standard NMF at elbow-selected k
#
# Inputs:  results/precomputed/tar_data_tcgacptac.rds
# Outputs: results/precomputed/desurv_bo_results_tcgacptac.rds
#          results/precomputed/tar_k_selection_tcgacptac.rds
#          results/precomputed/tar_params_best_tcgacptac.rds
#          results/precomputed/tar_data_filtered_tcgacptac.rds
#          results/desurv_bo_finalists_tcgacptac.rds (+ _alpha0_): finalist
#            re-evaluations and the paired 1-SE k selection
#          (+ alpha0 and elbow variants)
#
# Runtime: ~4-8 hours on HPC (30 cores), ~5 min in quick mode

message("=== Step 3: Bayesian Optimization ===")
source("code/00_helpers.R")
library(DeSurv)
library(NMF)
source("R/bo_helpers.R")
source("R/pick_k_elbow.R")
source("R/reported_settings.R")

# ── Configuration ─────────────────────────────────────────────────────────
# Rerun plan rev 3 (docs/plans/2026-09-30-rerun-plan-rev3.md), section 2:
# every BO evaluation scores the reported estimator (the code/04 recipe: seed
# fits -> consensus seed -> final fit, DeSurv::desurv_reported_fit) in each
# training fold, with code/04's seed-fit count and tolerances, on 3 repeated CV
# partitions with the same partitions and per-fold seed banks for every
# configuration; the incumbent is the GP posterior-mean argmax; k is chosen
# among re-evaluated finalists by the 1-SE rule on the empirical paired SE.
BO_ESTIMATOR    <- "reported"           # "mean_starts" = pre-rev-2 evaluator
BO_INCUMBENT    <- "posterior_mean"     # "observed" = pre-rev-2 argmax
REPORTED        <- desurv_reported_settings(CONFIG$quick)    # = final-fit settings
BO_N_SEED_FITS  <- REPORTED$n_seed_fits # code/04 NINIT_FULL: 100 full, 3 quick
CV_SEED_BANK_BASE <- 4000000L           # plan section 9: per-fold bank 4,000,001 onward
if (CONFIG$quick) {
  NINIT          <- 2L                  # used only if BO_ESTIMATOR == "mean_starts"
  BO_N_INIT      <- 4L
  BO_N_ITER      <- 4L
  NCORES_GRID    <- 1L
  PARALLEL       <- FALSE
  STD_NMF_NRUN   <- 2L
  STD_NMF_K_GRID <- 2:6
  BO_K_UPPER     <- 6L
  BO_NGENE       <- 500L
  BO_CV_SEEDS    <- 20261101L           # one partition in quick mode
  FINALISTS_PER_K <- 2L
} else {
  NINIT          <- 30L                 # used only if BO_ESTIMATOR == "mean_starts"
  BO_N_INIT      <- 50L
  BO_N_ITER      <- 100L
  NCORES_GRID    <- CONFIG$ncores
  PARALLEL       <- CONFIG$ncores > 1
  STD_NMF_NRUN   <- 30L
  STD_NMF_K_GRID <- 2:12
  BO_K_UPPER     <- 12L
  BO_NGENE       <- 3000L
  BO_CV_SEEDS    <- c(20261101L, 20261102L, 20261103L)   # plan section 9
  FINALISTS_PER_K <- 3L
}
# Seed-fit and final-fit tolerances inside CV = the final fit's (one source).
BO_TOL   <- REPORTED$tol
BO_MAXIT <- REPORTED$maxit

BO_BOUNDS <- list(
  k_grid     = list(lower = 2L, upper = BO_K_UPPER, type = "integer"),
  alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
  nu_grid    = list(lower = 0, upper = 1, type = "continuous")
)

BO_FIXED <- list(
  ngene        = BO_NGENE,
  lambdaW_grid = 0,
  lambdaH_grid = 0
)
if (identical(BO_ESTIMATOR, "mean_starts")) BO_FIXED$n_starts <- NINIT

# ── ntop handling ──────────────────────────────────────────────────────────
BO_BOUNDS$ntop <- list(lower = CONFIG$ntop_lower, upper = CONFIG$ntop_upper,
                       type = "integer")
message("  ntop BO bounds: [", CONFIG$ntop_lower, ", ", CONFIG$ntop_upper, "]")

BO_COMMON <- list(
  preprocess        = TRUE,
  method_trans_train = "rank",
  engine            = "warmstart",
  nfolds            = 5L,
  tol               = BO_TOL,
  maxit             = BO_MAXIT,
  max_refinements   = 0L,
  tol_gain          = 0.002,
  plateau           = 1L,
  top_k             = 10L,
  shrink_base       = 0.3,
  importance_gain   = 0.1,
  coarse_control    = list(n_init = BO_N_INIT, n_iter = BO_N_ITER, candidate_pool = 4000L,
                           exploration_weight = 0.01, seed = 123),
  refine_control    = NULL,
  verbose           = TRUE,
  parallel_grid     = PARALLEL,
  ncores_grid       = NCORES_GRID,
  estimator         = BO_ESTIMATOR,
  incumbent         = BO_INCUMBENT,
  cv_seeds          = BO_CV_SEEDS,
  n_seed_fits       = BO_N_SEED_FITS,
  seed_bank_base    = CV_SEED_BANK_BASE,
  reported_control  = list(tol_init = REPORTED$tol_init, imaxit = REPORTED$imaxit)
)

# Finalists: top FINALISTS_PER_K configurations per k by GP posterior mean,
# re-evaluated on the identical repeated CV (same partitions and seed banks).
run_bo_finalists <- function(bo_results) {
  DeSurv::desurv_bo_finalists(
    bo_results,
    X = tar_data$ex,
    y = tar_data$sampInfo$time,
    d = tar_data$sampInfo$event,
    dataset = tar_data$sampInfo$dataset,
    samp_keeps = tar_data$samp_keeps,
    n_per_k = FINALISTS_PER_K,
    cv_overrides = list(parallel_grid = PARALLEL, ncores_grid = NCORES_GRID)
  )
}

# ngene is fixed in BO_FIXED (not tuned), so the filtered data use that value.
bo_ngene <- function(params) {
  if (!is.null(params$ngene) && !is.na(params$ngene)) {
    as.integer(round(params$ngene))
  } else {
    as.integer(BO_FIXED$ngene)
  }
}

# ── Load training data ────────────────────────────────────────────────────
tar_data <- load_precomputed("tar_data_tcgacptac")

# ── 1. DeSurv BO (supervised) ────────────────────────────────────────────
desurv_bo_results <- cache_or_compute("desurv_bo_results_tcgacptac", {
  do.call(DeSurv::desurv_cv_bayesopt_refine, c(
    list(
      X = tar_data$ex,
      y = tar_data$sampInfo$time,
      d = tar_data$sampInfo$event,
      dataset = tar_data$sampInfo$dataset,
      samp_keeps = tar_data$samp_keeps,
      coarse_bounds = BO_BOUNDS,
      bo_fixed = BO_FIXED
    ),
    BO_COMMON
  ))
})

# Finalists and k selection (1-SE rule on the empirical paired CV SE)
desurv_bo_finalists <- cache_or_compute("desurv_bo_finalists_tcgacptac", {
  run_bo_finalists(desurv_bo_results)
})

tar_k_selection <- cache_or_compute("tar_k_selection_tcgacptac", {
  k_selection_from_finalists(desurv_bo_finalists)
})

# All hyperparameters from the selected k's own (finalist) row.
tar_params_best <- cache_or_compute("tar_params_best_tcgacptac", {
  select_bo_params_finalists(desurv_bo_finalists)
})
stopifnot(identical(tar_params_best$k, tar_k_selection$k_selected))

# Preprocess training data with BO-selected ngene
source("R/preprocess_helpers.R")
tar_data_filtered <- cache_or_compute("tar_data_filtered_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = bo_ngene(tar_params_best),
                           method_trans_train = "rank")
})

message("  DeSurv BO: k=", tar_params_best$k, " (runner-up k=", tar_k_selection$runner_up_k,
        "), alpha=", round(tar_params_best$alpha, 3))

# ── 2. Standard NMF BO (alpha = 0) ──────────────────────────────────────
desurv_bo_results_alpha0 <- cache_or_compute("desurv_bo_results_alpha0_tcgacptac", {
  bounds_a0 <- BO_BOUNDS
  bounds_a0$alpha_grid <- NULL
  bo_fixed_a0 <- c(BO_FIXED, list(alpha_grid = 0))
  do.call(DeSurv::desurv_cv_bayesopt_refine, c(
    list(
      X = tar_data$ex,
      y = tar_data$sampInfo$time,
      d = tar_data$sampInfo$event,
      dataset = tar_data$sampInfo$dataset,
      samp_keeps = tar_data$samp_keeps,
      coarse_bounds = bounds_a0,
      bo_fixed = bo_fixed_a0
    ),
    BO_COMMON
  ))
})

# Same selection as the main BO (finalists + paired 1-SE), so the control's
# hyperparameters also come from its selected k's own row.
desurv_bo_finalists_alpha0 <- cache_or_compute("desurv_bo_finalists_alpha0_tcgacptac", {
  run_bo_finalists(desurv_bo_results_alpha0)
})
tar_k_selection_alpha0 <- k_selection_from_finalists(desurv_bo_finalists_alpha0)

tar_params_best_alpha0 <- cache_or_compute("tar_params_best_alpha0_tcgacptac", {
  params <- select_bo_params_finalists(desurv_bo_finalists_alpha0)
  params$alpha <- 0
  params
})

ngene_alpha0 <- bo_ngene(tar_params_best_alpha0)

tar_data_filtered_alpha0 <- cache_or_compute("tar_data_filtered_alpha0_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = ngene_alpha0,
                           method_trans_train = "rank")
})

message("  Alpha=0 BO: k=", tar_params_best_alpha0$k)

# ── 3. Standard NMF for elbow k selection ────────────────────────────────
STD_NMF_SEED <- 20261201L   # rev 3: fixed seed for the standard-NMF rank survey (was unseeded)
fit_std <- cache_or_compute("fit_std_tcgacptac", {
  # Numeric seed: NMF derives L'Ecuyer RNG streams from it for the nrun runs, so the
  # result is reproducible in parallel mode too (?NMF::nmf, argument `seed`).
  NMF::nmf(tar_data_filtered$ex, STD_NMF_K_GRID,
           nrun = STD_NMF_NRUN, method = "lee",
           seed = STD_NMF_SEED,
           .options = paste0("p", STD_NMF_NRUN))
})

std_nmf_selected_k <- cache_or_compute("std_nmf_selected_k_tcgacptac", {
  ranks <- fit_std$measures$rank
  resids <- fit_std$measures$residuals
  pick_k_elbow(ranks, resids)
})

message("  Standard NMF elbow k=", std_nmf_selected_k)

# ── 4. BO at elbow k (for comparison) ────────────────────────────────────
desurv_bo_results_elbowk <- cache_or_compute("desurv_bo_results_elbowk_tcgacptac", {
  bounds_ek <- BO_BOUNDS
  bounds_ek$k_grid <- NULL
  bo_fixed_ek <- c(BO_FIXED, list(k_grid = std_nmf_selected_k))
  do.call(DeSurv::desurv_cv_bayesopt_refine, c(
    list(
      X = tar_data$ex,
      y = tar_data$sampInfo$time,
      d = tar_data$sampInfo$event,
      dataset = tar_data$sampInfo$dataset,
      samp_keeps = tar_data$samp_keeps,
      coarse_bounds = bounds_ek,
      bo_fixed = bo_fixed_ek
    ),
    BO_COMMON
  ))
})

tar_params_best_elbowk <- cache_or_compute("tar_params_best_elbowk_tcgacptac", {
  standardize_bo_params(desurv_bo_results_elbowk$overall_best$params)
})

ngene_elbowk <- bo_ngene(tar_params_best_elbowk)

tar_data_filtered_elbowk <- cache_or_compute("tar_data_filtered_elbowk_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = ngene_elbowk,
                           method_trans_train = "rank")
})

# ── 5. Fair matched-rank comparator: alpha = 0 at DeSurv's own k ─────────
# The matched-rank comparison in the Results uses standard NMF at k = 3, which
# receives neither DeSurv's consensus initialization nor its tuning budget, so
# that contrast cannot separate supervision from tuning. This block removes both
# confounds by composing the two pinning idioms already used above: alpha is
# dropped from the search as in section 2 (:141-143) and k is dropped as in
# section 4 (:199-201). The comparator therefore gets the same optimizer, the
# same BO budget, the same ngene and the same folds as DeSurv. Only supervision
# is switched off, and rank is held at DeSurv's selected k.
desurv_bo_results_a0k3 <- cache_or_compute("desurv_bo_results_a0k3_tcgacptac", {
  bounds_a0k3 <- BO_BOUNDS
  bounds_a0k3$alpha_grid <- NULL
  bounds_a0k3$k_grid     <- NULL
  bo_fixed_a0k3 <- c(BO_FIXED, list(alpha_grid = 0, k_grid = tar_params_best$k))
  do.call(DeSurv::desurv_cv_bayesopt_refine, c(
    list(
      X = tar_data$ex,
      y = tar_data$sampInfo$time,
      d = tar_data$sampInfo$event,
      dataset = tar_data$sampInfo$dataset,
      samp_keeps = tar_data$samp_keeps,
      coarse_bounds = bounds_a0k3,
      bo_fixed = bo_fixed_a0k3
    ),
    BO_COMMON
  ))
})

# k and alpha were both fixed, so BO does not return them; set them explicitly.
# select_bo_k_by_cv_se() is deliberately NOT called here: with k out of
# coarse_bounds the history carries no k_grid column and the helper returns
# reason = "missing_k" (R/bo_helpers.R). The elbow-k block above skips it for
# the same reason.
tar_params_best_a0k3 <- cache_or_compute("tar_params_best_a0k3_tcgacptac", {
  params <- standardize_bo_params(desurv_bo_results_a0k3$overall_best$params)
  params$k     <- tar_params_best$k
  params$alpha <- 0
  params
})

# ngene is fixed in BO_FIXED for every variant, so the comparator shares
# tar_data_filtered_tcgacptac rather than getting its own filtered object.
# (ngene_value above is scoped inside a cache_or_compute expr and does not
# exist here on a cache hit, so recompute the expectation from BO_FIXED.)
local({
  .want <- as.integer(BO_FIXED$ngene)
  .got  <- tar_params_best_a0k3$ngene
  if (!is.null(.got) && !is.na(.got)) stopifnot(as.integer(round(.got)) == .want)
})

message("  Fair matched-rank BO (alpha=0, k fixed): k=", tar_params_best_a0k3$k,
        " lambda=", signif(tar_params_best_a0k3$lambda, 3),
        " nu=", signif(tar_params_best_a0k3$nu, 3),
        " ntop=", tar_params_best_a0k3$ntop)

message("=== Step 3 complete ===")
