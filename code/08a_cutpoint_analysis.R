#!/usr/bin/env Rscript
# code/08a_cutpoint_analysis.R — CV-based cutpoint selection for dichotomized risk
#
# Computes the optimal z-score cutpoint that converts the continuous DeSurv
# (and standard NMF at DeSurv-k) linear predictor into a binary risk
# stratification. Outputs the per-cutpoint summary and the LP statistics
# (training mean / sd plus optimal z) needed by both the main-text median
# survival KM curves (08_figures.R) and the SI appendix figures
# (08b_si_figures.R).
#
# Inputs:  model fits + filtered training data from steps 02-04
# Outputs: results/precomputed/<sub>/desurv_cutpoint_summary_tcgacptac.rds
#          results/precomputed/<sub>/desurv_lp_stats_tcgacptac.rds
#          results/precomputed/<sub>/std_desurvk_cutpoint_summary_tcgacptac.rds
#          results/precomputed/<sub>/std_desurvk_lp_stats_tcgacptac.rds
#
# Runtime: ~5-15 min (depends on n_starts)

message("=== Step 8a: Cutpoint analysis ===")
source("code/00_helpers.R")
library(dplyr)
library(survival)
library(NMF)
library(glmnet)

source("R/cv_grid_helpers.R")
source("R/fit_cox_model.R")

# ── Load prerequisites ────────────────────────────────────────────────────
tar_fit_desurv    <- load_precomputed("tar_fit_desurv_tcgacptac")
fit_std_desurvk   <- load_precomputed("fit_std_desurvk_tcgacptac")
tar_data_filtered <- load_precomputed("tar_data_filtered_tcgacptac")
tar_params_best   <- load_precomputed("tar_params_best_tcgacptac")

# ntop for LP computation (matches main pipeline)
ntop_for_lp <- if (!is.null(CONFIG$ntop_value)) CONFIG$ntop_value else tar_params_best$ntop

# ── Helpers ──────────────────────────────────────────────────────────────
summarize_cutpoint <- function(cv_result, z_grid = seq(-2.0, 2.0, by = 0.2)) {
  evaluate_cutpoint_zscores(cv_result, z_grid) |>
    dplyr::group_by(z_cutpoint) |>
    dplyr::summarise(
      mean_cindex_dichot = mean(cindex_dichot, na.rm = TRUE),
      se_cindex_dichot   = sd(cindex_dichot, na.rm = TRUE) /
        sqrt(sum(!is.na(cindex_dichot))),
      mean_abs_logrank_z = mean(abs(logrank_z), na.rm = TRUE),
      se_abs_logrank_z   = sd(abs(logrank_z), na.rm = TRUE) /
        sqrt(sum(!is.na(logrank_z))),
      .groups = "drop"
    )
}

pick_optimal <- function(summary_tbl, metric = c("logrank", "cindex")) {
  metric <- match.arg(metric)
  col <- if (metric == "logrank") "mean_abs_logrank_z" else "mean_cindex_dichot"
  summary_tbl |>
    dplyr::slice_max(.data[[col]], n = 1, with_ties = FALSE) |>
    dplyr::pull(z_cutpoint)
}

build_lp_stats <- function(fit, data_filtered, ntop, summary_tbl) {
  lp <- compute_lp(fit$W, fit$beta, data_filtered$ex, ntop)
  lp_mean <- mean(lp, na.rm = TRUE)
  lp_sd   <- sd(lp, na.rm = TRUE)
  z_lr <- pick_optimal(summary_tbl, "logrank")
  z_ci <- pick_optimal(summary_tbl, "cindex")
  list(
    lp_mean                   = lp_mean,
    lp_sd                     = lp_sd,
    optimal_z_cutpoint        = z_lr,
    cutpoint_abs              = z_lr * lp_sd + lp_mean,
    optimal_z_cutpoint_cindex = z_ci,
    cutpoint_abs_cindex       = z_ci * lp_sd + lp_mean
  )
}

# ── DeSurv cutpoint (5-fold CV at the BO-selected (k, alpha)) ────────────
desurv_cv_cutpoint_result <- cache_or_compute(
  "desurv_cv_cutpoint_result_tcgacptac",
  run_cv_grid_point(
    data         = tar_data_filtered,
    k            = tar_params_best$k,
    alpha        = tar_params_best$alpha,
    fixed_params = list(
      lambda  = tar_params_best$lambda,
      nu      = tar_params_best$nu,
      lambdaW = 0,
      lambdaH = 0,
      ntop    = ntop_for_lp
    ),
    nfolds   = 5,
    n_starts = 30,
    seed     = 125,
    verbose  = FALSE
  )
)

desurv_cutpoint_summary <- cache_or_compute(
  "desurv_cutpoint_summary_tcgacptac",
  summarize_cutpoint(desurv_cv_cutpoint_result)
)

desurv_lp_stats <- cache_or_compute(
  "desurv_lp_stats_tcgacptac",
  build_lp_stats(tar_fit_desurv, tar_data_filtered,
                 ntop_for_lp, desurv_cutpoint_summary)
)

# ── Standard NMF (at DeSurv-selected k) cutpoint ─────────────────────────
std_desurvk_cv_cutpoint_result <- cache_or_compute(
  "std_desurvk_cv_cutpoint_result_tcgacptac",
  run_cv_grid_point_std_nmf(
    data   = tar_data_filtered,
    k      = tar_params_best$k,
    nrun   = 30,
    nfolds = 5,
    seed   = 123
  )
)

std_desurvk_cutpoint_summary <- cache_or_compute(
  "std_desurvk_cutpoint_summary_tcgacptac",
  summarize_cutpoint(std_desurvk_cv_cutpoint_result)
)

std_desurvk_lp_stats <- cache_or_compute(
  "std_desurvk_lp_stats_tcgacptac",
  build_lp_stats(fit_std_desurvk, tar_data_filtered,
                 ntop = NULL, std_desurvk_cutpoint_summary)
)

message(sprintf("  DeSurv optimal z (logrank): %.2f   (cindex): %.2f",
                desurv_lp_stats$optimal_z_cutpoint,
                desurv_lp_stats$optimal_z_cutpoint_cindex))
message(sprintf("  Std-NMF optimal z (logrank): %.2f   (cindex): %.2f",
                std_desurvk_lp_stats$optimal_z_cutpoint,
                std_desurvk_lp_stats$optimal_z_cutpoint_cindex))

message("=== Step 8a complete ===")
