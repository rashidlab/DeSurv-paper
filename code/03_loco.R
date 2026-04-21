#!/usr/bin/env Rscript
# code/03_loco.R — BO with leave-one-cohort-out (LOCO) CV
#
# Same as code/03_bayesian_optimization.R but uses 2-fold LOCO CV
# (train on TCGA → validate on CPTAC, and vice versa) instead of
# random 5-fold stratified CV.
#
# K selection uses best mean LOCO c-index (no 1-SE rule — meaningless
# with 2 folds).
#
# Outputs to: results/precomputed/loco/
#
# Usage:
#   DESURV_NCORES=30 Rscript code/03_loco.R
#   DESURV_QUICK=TRUE Rscript code/03_loco.R   # smoke test

message("=== Step 3 (LOCO CV): Bayesian Optimization ===")

# Force subfolder to "loco" via env var, then source helpers
Sys.setenv(DESURV_CV_SUBFOLDER = "loco")
source("code/00_helpers_loco.R")
library(DeSurv)
library(NMF)
source("R/bo_helpers.R")
source("R/bo_helpers_loco.R")
source("R/pick_k_elbow.R")

# ── Configuration ─────────────────────────────────────────────────────────
if (CONFIG$quick) {
  NINIT        <- 2L
  BO_N_INIT    <- 4L
  BO_N_ITER    <- 4L
  NCORES_GRID  <- 1L
  PARALLEL     <- FALSE
  STD_NMF_NRUN <- 2L
  STD_NMF_K_GRID <- 2:6
  BO_K_UPPER   <- 6L
} else {
  NINIT        <- 30L
  BO_N_INIT    <- 50L
  BO_N_ITER    <- 150L       # more iterations to compensate for noisier 2-fold
  NCORES_GRID  <- CONFIG$ncores
  PARALLEL     <- CONFIG$ncores > 1
  STD_NMF_NRUN <- 30L
  STD_NMF_K_GRID <- 2:12
  BO_K_UPPER   <- 12L
}

BO_BOUNDS <- list(
  k_grid     = list(lower = 2L, upper = BO_K_UPPER, type = "integer"),
  alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
  nu_grid    = list(lower = 0, upper = 1, type = "continuous")
)

BO_FIXED <- list(
  n_starts     = NINIT,
  ngene        = 3000L,
  lambdaW_grid = 0,
  lambdaH_grid = 0
)

# ── Build LOCO folds ──────────────────────────────────────────────────────
tar_data <- load_precomputed("tar_data_tcgacptac")

# samp_keeps gives the indices of samples used in CV
kept_datasets <- tar_data$sampInfo$dataset[tar_data$samp_keeps]
cohorts <- unique(kept_datasets)
if (length(cohorts) != 2) {
  stop("LOCO CV requires exactly 2 cohorts, found: ",
       paste(cohorts, collapse = ", "))
}

# Build folds vector of length = number of kept samples (not all samples).
# desurv_cv applies samp_keeps first, then checks length(folds) == ncol(X_kept).
loco_folds <- as.integer(factor(kept_datasets))

message("  LOCO folds: ", paste(cohorts, "->", seq_along(cohorts), collapse = ", "))
message("  Fold sizes: ", paste(table(loco_folds[tar_data$samp_keeps]), collapse = ", "))

BO_COMMON <- list(
  preprocess        = TRUE,
  method_trans_train = "rank",
  engine            = "warmstart",
  nfolds            = 2L,
  folds             = loco_folds,
  tol               = 1e-5,
  maxit             = 4000L,
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
  ncores_grid       = NCORES_GRID
)

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

# K selection: best mean LOCO c-index, then tolerance rule for alpha
tar_k_selection <- cache_or_compute("tar_k_selection_tcgacptac", {
  select_bo_k_best_mean(desurv_bo_results, tol = 0.005)
})
message("  LOCO k selection: k=", tar_k_selection$k_selected,
        " (best mean=", round(tar_k_selection$best_mean, 4),
        ", selected alpha=", round(tar_k_selection$selected_alpha, 4),
        ", selected cindex=", round(tar_k_selection$selected_cindex, 4), ")")

# Extract params from the tolerance-selected evaluation
tar_params_best <- cache_or_compute("tar_params_best_tcgacptac", {
  params <- tar_k_selection$params
  if (is.null(params)) {
    # Fallback to overall_best if selection failed
    params <- standardize_bo_params(desurv_bo_results$overall_best$params)
  }
  params$k <- tar_k_selection$k_selected
  params
})

# Preprocess training data
source("R/targets_config.R")
tar_data_filtered <- cache_or_compute("tar_data_filtered_tcgacptac", {
  ngene_value <- if (!is.null(tar_params_best$ngene) && !is.na(tar_params_best$ngene)) {
    as.integer(round(tar_params_best$ngene))
  } else {
    3000L
  }
  preprocess_training_data(data = tar_data, ngene = ngene_value,
                           method_trans_train = "rank")
})

message("  DeSurv LOCO BO: k=", tar_params_best$k, ", alpha=", round(tar_params_best$alpha, 3))

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

tar_k_selection_alpha0 <- select_bo_k_best_mean(desurv_bo_results_alpha0)

tar_params_best_alpha0 <- cache_or_compute("tar_params_best_alpha0_tcgacptac", {
  params <- standardize_bo_params(desurv_bo_results_alpha0$overall_best$params)
  if (!is.null(tar_k_selection_alpha0$k_selected)) {
    params$k <- tar_k_selection_alpha0$k_selected
  }
  params$alpha <- 0
  params
})

ngene_alpha0 <- if (!is.null(tar_params_best_alpha0$ngene) && !is.na(tar_params_best_alpha0$ngene)) {
  as.integer(round(tar_params_best_alpha0$ngene))
} else {
  3000L
}

tar_data_filtered_alpha0 <- cache_or_compute("tar_data_filtered_alpha0_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = ngene_alpha0,
                           method_trans_train = "rank")
})

message("  Alpha=0 LOCO BO: k=", tar_params_best_alpha0$k)

# ── 3. Standard NMF for elbow k selection ────────────────────────────────
fit_std <- cache_or_compute("fit_std_tcgacptac", {
  NMF::nmf(tar_data_filtered$ex, STD_NMF_K_GRID,
           nrun = STD_NMF_NRUN, method = "lee",
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

ngene_elbowk <- if (!is.null(tar_params_best_elbowk$ngene) && !is.na(tar_params_best_elbowk$ngene)) {
  as.integer(round(tar_params_best_elbowk$ngene))
} else {
  3000L
}

tar_data_filtered_elbowk <- cache_or_compute("tar_data_filtered_elbowk_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = ngene_elbowk,
                           method_trans_train = "rank")
})

message("=== Step 3 (LOCO) complete ===")
