#!/usr/bin/env Rscript
# code/04_gridcv.R — Fit final models using Grid CV-selected hyperparameters
#
# Identical to code/04_loco.R but reads/writes from results/precomputed/gridcv/

message("=== Step 4 (Grid CV): Fit Models ===")
source("code/00_helpers_gridcv.R")
library(DeSurv)
library(NMF)
library(glmnet)
library(survival)
library(dplyr)
source("R/fit_cox_model.R")

# ── Configuration ─────────────────────────────────────────────────────────
if (CONFIG$quick) {
  NINIT_FULL <- 3L
  RUN_TOL    <- 1e-4
  RUN_MAXIT  <- 500L
} else {
  NINIT_FULL <- 100L
  RUN_TOL    <- 1e-5
  RUN_MAXIT  <- 4000L
}

# ── Load prerequisites ────────────────────────────────────────────────────
tar_data_filtered  <- load_precomputed("tar_data_filtered_tcgacptac")
tar_params_best    <- load_precomputed("tar_params_best_tcgacptac")
fit_std            <- load_precomputed("fit_std_tcgacptac")
std_nmf_selected_k <- load_precomputed("std_nmf_selected_k_tcgacptac")
tar_data_filtered_elbowk <- load_precomputed("tar_data_filtered_elbowk_tcgacptac")
tar_params_best_elbowk   <- load_precomputed("tar_params_best_elbowk_tcgacptac")

extract_param <- function(params, name, default) {
  val <- params[[name]]
  if (is.null(val) || is.na(val)) default else val
}

lambdaW_value <- extract_param(tar_params_best, "lambdaW", 0)
lambdaH_value <- extract_param(tar_params_best, "lambdaH", 0)
# ntop not tuned; use 100 for consensus initialization
ntop_value <- 100L

# ── Helper: multi-start fitting ──────────────────────────────────────────
run_seed_fits <- function(data, params, lambdaW, lambdaH, ninit, label) {
  seeds <- seq_len(ninit)
  fit_one_seed <- function(seed_i) {
    fit_i <- try(
      desurv_fit(
        X = data$ex, y = data$sampInfo$time, d = data$sampInfo$event,
        k = params$k, alpha = params$alpha, lambda = params$lambda,
        nu = params$nu, lambdaW = lambdaW, lambdaH = lambdaH,
        seed = seed_i, tol = RUN_TOL / 100, tol_init = RUN_TOL,
        maxit = RUN_MAXIT, imaxit = RUN_MAXIT, ninit = 1,
        parallel_init = FALSE, verbose = FALSE
      ),
      silent = TRUE
    )
    if (!inherits(fit_i, "try-error") && inherits(fit_i, "desurv_fit")) {
      list(fit = fit_i, cindex = if (!is.null(fit_i$cindex)) fit_i$cindex else NA_real_)
    } else {
      list(fit = NULL, cindex = NA_real_)
    }
  }
  if (CONFIG$ncores > 1) {
    results <- parallel::mclapply(seeds, fit_one_seed, mc.cores = CONFIG$ncores)
  } else {
    results <- lapply(seeds, fit_one_seed)
  }
  fits <- lapply(results, `[[`, "fit")
  scores <- vapply(results, function(r) r$cindex, numeric(1))
  keep <- !vapply(fits, is.null, logical(1))
  if (!any(keep)) stop("No successful fits for ", label)
  message(sprintf("  %s: %d/%d successful seed fits", label, sum(keep), ninit))
  list(fits = fits[keep], seeds = seeds[keep], cindex = scores[keep])
}

# ── Helper: consensus → final fit ────────────────────────────────────────
consensus_and_fit <- function(seed_fits, data, params, lambdaW, lambdaH, ntop, label) {
  init_vals <- DeSurv::desurv_consensus_seed(
    fits = seed_fits$fits, X = data$ex, ntop = ntop,
    k = params$k, min_frequency = max(1L, as.integer(ceiling(0.3 * length(seed_fits$fits))))
  )
  fit <- desurv_fit(
    X = data$ex, y = data$sampInfo$time, d = data$sampInfo$event,
    k = params$k, alpha = params$alpha, lambda = params$lambda,
    nu = params$nu, lambdaW = lambdaW, lambdaH = lambdaH,
    W0 = init_vals$W0, H0 = init_vals$H0, beta0 = init_vals$beta0,
    seed = NULL, tol = RUN_TOL / 100, maxit = RUN_MAXIT, verbose = FALSE
  )
  message(sprintf("  %s: final fit cindex=%.3f", label, fit$cindex))
  fit
}

# ── 1. DeSurv (supervised) ─────────────────────────────────────────
desurv_seed_fits <- cache_or_compute("desurv_seed_fits_tcgacptac", {
  run_seed_fits(tar_data_filtered, tar_params_best,
                lambdaW_value, lambdaH_value, NINIT_FULL, "DeSurv")
})

tar_fit_desurv <- cache_or_compute("tar_fit_desurv_tcgacptac", {
  consensus_and_fit(desurv_seed_fits, tar_data_filtered, tar_params_best,
                    lambdaW_value, lambdaH_value, ntop_value, "DeSurv")
})

# ── 1b. Alpha=0 model ─────────────────────────────────────────────
tar_params_best_alpha0    <- load_precomputed("tar_params_best_alpha0_tcgacptac")
tar_data_filtered_alpha0  <- load_precomputed("tar_data_filtered_alpha0_tcgacptac")

lambdaW_alpha0 <- extract_param(tar_params_best_alpha0, "lambdaW", 0)
lambdaH_alpha0 <- extract_param(tar_params_best_alpha0, "lambdaH", 0)
ntop_alpha0 <- 100L

desurv_seed_fits_alpha0 <- cache_or_compute("desurv_seed_fits_alpha0_tcgacptac", {
  run_seed_fits(tar_data_filtered_alpha0, tar_params_best_alpha0,
                lambdaW_alpha0, lambdaH_alpha0, NINIT_FULL, "Alpha=0")
})

tar_fit_desurv_alpha0 <- cache_or_compute("tar_fit_desurv_alpha0_tcgacptac", {
  consensus_and_fit(desurv_seed_fits_alpha0, tar_data_filtered_alpha0,
                    tar_params_best_alpha0, lambdaW_alpha0, lambdaH_alpha0,
                    ntop_alpha0, "Alpha=0")
})

# ── 2. Standard NMF at DeSurv k ───────────────────────────────────
fit_std_desurvk <- cache_or_compute("fit_std_desurvk_tcgacptac", {
  selected_k <- tar_params_best$k
  fit_nmf <- fit_std$fit[[as.character(selected_k)]]
  W <- fit_nmf@fit@W
  X <- tar_data_filtered$ex
  XtW <- t(X) %*% W
  df <- data.frame(
    time = tar_data_filtered$sampInfo$time,
    event = tar_data_filtered$sampInfo$event
  )
  beta <- fit_cox_model(XtW, df, nfold = 5)
  list(W = W, H = fit_nmf@fit@H, beta = beta)
})

# ── 3. Standard NMF at elbow k ────────────────────────────────────
fit_std_elbowk <- cache_or_compute("fit_std_elbowk_tcgacptac", {
  fit_nmf <- fit_std$fit[[as.character(std_nmf_selected_k)]]
  W <- fit_nmf@fit@W
  X <- tar_data_filtered_elbowk$ex
  XtW <- t(X) %*% W
  df <- data.frame(
    time = tar_data_filtered_elbowk$sampInfo$time,
    event = tar_data_filtered_elbowk$sampInfo$event
  )
  beta <- fit_cox_model(XtW, df, nfold = 5)
  list(W = W, H = fit_nmf@fit@H, beta = beta)
})

message("=== Step 4 (Grid CV) complete ===")
