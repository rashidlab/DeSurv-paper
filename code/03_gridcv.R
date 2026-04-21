#!/usr/bin/env Rscript
# code/03_gridcv.R — Grid CV with standard 5-fold stratified CV
#
# Instead of Bayesian optimization, evaluates a full grid of
# (k, alpha, lambda, nu) via standard 5-fold CV.
# ntop is NOT tuned; consensus initialization uses ntop=100.
#
# Outputs to: results/precomputed/gridcv/
#
# Usage:
#   DESURV_NCORES=30 Rscript code/03_gridcv.R
#   DESURV_QUICK=TRUE Rscript code/03_gridcv.R   # smoke test

message("=== Step 3 (Grid CV): Grid Search ===")

source("code/00_helpers_gridcv.R")

library(DeSurv)
library(NMF)
source("R/bo_helpers.R")
source("R/pick_k_elbow.R")

# ── Configuration ─────────────────────────────────────────────────────────
if (CONFIG$quick) {
  NINIT        <- 2L
  NCORES_GRID  <- 1L
  PARALLEL     <- FALSE
  STD_NMF_NRUN <- 2L
  STD_NMF_K_GRID <- 2:6
  K_GRID       <- c(2L, 4L)
  ALPHA_GRID   <- c(0.0, 0.5, 1.0)
  LAMBDA_GRID  <- c(0.01, 1)
  NU_GRID      <- c(0.0, 0.5)
} else {
  NINIT        <- 30L
  NCORES_GRID  <- CONFIG$ncores
  PARALLEL     <- CONFIG$ncores > 1
  STD_NMF_NRUN <- 30L
  STD_NMF_K_GRID <- 2:12
  K_GRID       <- 2:12
  ALPHA_GRID   <- seq(0, 1, by = 0.1)
  LAMBDA_GRID  <- c(0.05, 0.1, 0.5, 1, 1.5)
  NU_GRID      <- c(0.05, 0.1, 0.15)
}

# Allow env var overrides for grid
if (nzchar(Sys.getenv("GRIDCV_K")))
  K_GRID <- eval(parse(text = Sys.getenv("GRIDCV_K")))
if (nzchar(Sys.getenv("GRIDCV_ALPHA")))
  ALPHA_GRID <- eval(parse(text = Sys.getenv("GRIDCV_ALPHA")))
if (nzchar(Sys.getenv("GRIDCV_LAMBDA")))
  LAMBDA_GRID <- eval(parse(text = Sys.getenv("GRIDCV_LAMBDA")))
if (nzchar(Sys.getenv("GRIDCV_NU")))
  NU_GRID <- eval(parse(text = Sys.getenv("GRIDCV_NU")))

FIXED_PARAMS <- list(
  n_starts     = NINIT,
  ngene        = 3000L,
  lambdaW_grid = 0,
  lambdaH_grid = 0
)

grid_full <- expand.grid(
  k      = K_GRID,
  alpha  = ALPHA_GRID,
  lambda = LAMBDA_GRID,
  nu     = NU_GRID,
  KEEP.OUT.ATTRS = FALSE
)

message("  Grid size: ", nrow(grid_full),
        " (k=", length(K_GRID),
        " x alpha=", length(ALPHA_GRID),
        " x lambda=", length(LAMBDA_GRID),
        " x nu=", length(NU_GRID), ")")

# ── Load training data ────────────────────────────────────────────────────
tar_data <- load_precomputed("tar_data_tcgacptac")

# ── Grid CV helper ────────────────────────────────────────────────────────
eval_one_grid_point <- function(i, grid, tar_data, fixed) {
  row <- grid[i, ]
  args <- list(
    X = tar_data$ex,
    y = tar_data$sampInfo$time,
    d = tar_data$sampInfo$event,
    dataset = tar_data$sampInfo$dataset,
    samp_keeps = tar_data$samp_keeps,
    k_grid     = row$k,
    alpha_grid = row$alpha,
    lambda_grid = row$lambda,
    nu_grid    = row$nu,
    lambdaW_grid = fixed$lambdaW_grid,
    lambdaH_grid = fixed$lambdaH_grid,
    n_starts   = fixed$n_starts,
    ngene      = fixed$ngene,
    preprocess = TRUE,
    method_trans_train = "rank",
    engine     = "warmstart",
    nfolds     = 5L,
    tol        = 1e-5,
    maxit      = 4000L,
    cv_only    = TRUE,
    verbose    = FALSE
  )
  result <- tryCatch(
    do.call(DeSurv::desurv_cv, args),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    data.frame(
      k = row$k, alpha = row$alpha, lambda = row$lambda, nu = row$nu,
      mean_cindex = NA_real_, se_cindex = NA_real_, status = "error",
      message = conditionMessage(result),
      stringsAsFactors = FALSE
    )
  } else {
    se_val <- if (!is.null(result$summary$se_cindex)) result$summary$se_cindex[1L] else NA_real_
    data.frame(
      k = row$k, alpha = row$alpha, lambda = row$lambda, nu = row$nu,
      mean_cindex = result$summary$mean_cindex[1L],
      se_cindex = se_val,
      status = "ok", message = "",
      stringsAsFactors = FALSE
    )
  }
}

run_grid_cv <- function(grid, tar_data, fixed, ncores, parallel) {
  n <- nrow(grid)
  message("  Evaluating ", n, " grid points ...")
  t0 <- Sys.time()

  if (parallel && ncores > 1) {
    results <- parallel::mclapply(
      seq_len(n),
      eval_one_grid_point,
      grid = grid, tar_data = tar_data,
      fixed = fixed,
      mc.cores = ncores
    )
  } else {
    results <- lapply(
      seq_len(n),
      function(i) {
        if (i %% 50 == 0) message("    ", i, "/", n)
        eval_one_grid_point(i, grid, tar_data, fixed)
      }
    )
  }

  history <- do.call(rbind, results)
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  n_ok <- sum(history$status == "ok", na.rm = TRUE)
  message(sprintf("  Grid CV done: %d/%d ok (%.1f min)", n_ok, n, as.numeric(elapsed)))
  history
}

# Package results into the same format as desurv_cv_bo_refine
# so downstream code (k selection, param extraction) works unchanged
package_grid_results <- function(history) {
  ok <- history$status == "ok" & !is.na(history$mean_cindex)
  if (!any(ok)) stop("Grid CV failed: no successful evaluations.")

  best_idx <- which.max(history$mean_cindex[ok])
  best_row <- history[ok, , drop = FALSE][best_idx, , drop = FALSE]
  best_params <- c(
    k_grid      = best_row$k,
    alpha_grid  = best_row$alpha,
    lambda_grid = best_row$lambda,
    nu_grid     = best_row$nu
  )

  # Rename columns to match BO history format
  history_bo <- history
  names(history_bo)[names(history_bo) == "k"]      <- "k_grid"
  names(history_bo)[names(history_bo) == "alpha"]   <- "alpha_grid"
  names(history_bo)[names(history_bo) == "lambda"]  <- "lambda_grid"
  names(history_bo)[names(history_bo) == "nu"]      <- "nu_grid"

  structure(list(
    runs         = list(),
    bounds       = list(),
    history      = history_bo,
    overall_best = list(params = best_params, mean_cindex = best_row$mean_cindex),
    stop_reason  = "grid_complete"
  ), class = "desurv_cv_bo_refine")
}

# ── 1. DeSurv grid CV (supervised, alpha tuned) ─────────────────────────
desurv_grid_history <- cache_or_compute("desurv_grid_history_tcgacptac", {
  run_grid_cv(grid_full, tar_data, FIXED_PARAMS,
              NCORES_GRID, PARALLEL)
})

desurv_bo_results <- cache_or_compute("desurv_bo_results_tcgacptac", {
  package_grid_results(desurv_grid_history)
})

# K selection (1-SE rule)
tar_k_selection <- cache_or_compute("tar_k_selection_tcgacptac", {
  select_bo_k_by_cv_se(desurv_bo_results)
})

tar_params_best <- cache_or_compute("tar_params_best_tcgacptac", {
  params <- standardize_bo_params(desurv_bo_results$overall_best$params)
  if (!is.null(tar_k_selection$k_selected)) {
    params$k <- tar_k_selection$k_selected
  }
  params
})

source("R/targets_config.R")
tar_data_filtered <- cache_or_compute("tar_data_filtered_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = 3000L,
                           method_trans_train = "rank")
})

message("  DeSurv Grid CV: k=", tar_params_best$k,
        ", alpha=", round(tar_params_best$alpha, 3))

# ── 2. Alpha=0 grid CV ──────────────────────────────────────────────────
grid_alpha0 <- expand.grid(
  k      = K_GRID,
  alpha  = 0,
  lambda = LAMBDA_GRID,
  nu     = NU_GRID,
  KEEP.OUT.ATTRS = FALSE
)
message("  Alpha=0 grid size: ", nrow(grid_alpha0))

desurv_grid_history_alpha0 <- cache_or_compute("desurv_grid_history_alpha0_tcgacptac", {
  run_grid_cv(grid_alpha0, tar_data, FIXED_PARAMS,
              NCORES_GRID, PARALLEL)
})

desurv_bo_results_alpha0 <- cache_or_compute("desurv_bo_results_alpha0_tcgacptac", {
  package_grid_results(desurv_grid_history_alpha0)
})

tar_k_selection_alpha0 <- select_bo_k_by_cv_se(desurv_bo_results_alpha0)

tar_params_best_alpha0 <- cache_or_compute("tar_params_best_alpha0_tcgacptac", {
  params <- standardize_bo_params(desurv_bo_results_alpha0$overall_best$params)
  if (!is.null(tar_k_selection_alpha0$k_selected)) {
    params$k <- tar_k_selection_alpha0$k_selected
  }
  params$alpha <- 0
  params
})

tar_data_filtered_alpha0 <- cache_or_compute("tar_data_filtered_alpha0_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = 3000L,
                           method_trans_train = "rank")
})

message("  Alpha=0 Grid CV: k=", tar_params_best_alpha0$k)

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

# ── 4. Grid CV at elbow k (for comparison) ───────────────────────────────
grid_elbowk <- expand.grid(
  k      = std_nmf_selected_k,
  alpha  = ALPHA_GRID,
  lambda = LAMBDA_GRID,
  nu     = NU_GRID,
  KEEP.OUT.ATTRS = FALSE
)
message("  Elbow k grid size: ", nrow(grid_elbowk))

desurv_grid_history_elbowk <- cache_or_compute("desurv_grid_history_elbowk_tcgacptac", {
  run_grid_cv(grid_elbowk, tar_data, FIXED_PARAMS,
              NCORES_GRID, PARALLEL)
})

desurv_bo_results_elbowk <- cache_or_compute("desurv_bo_results_elbowk_tcgacptac", {
  package_grid_results(desurv_grid_history_elbowk)
})

tar_params_best_elbowk <- cache_or_compute("tar_params_best_elbowk_tcgacptac", {
  standardize_bo_params(desurv_bo_results_elbowk$overall_best$params)
})

tar_data_filtered_elbowk <- cache_or_compute("tar_data_filtered_elbowk_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = 3000L,
                           method_trans_train = "rank")
})

message("=== Step 3 (Grid CV) complete ===")
