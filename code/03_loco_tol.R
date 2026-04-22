#!/usr/bin/env Rscript
# code/03_loco_tol.R — Re-select hyperparameters from existing LOCO BO results
#                      using the tolerance rule for alpha
#
# Reads BO results and standard NMF fits from results/precomputed/loco/.
# Writes new k-selection, params, and preprocessed data to
# results/precomputed/loco_tol/.
#
# No expensive BO or NMF fitting is performed here.
#
# Usage:
#   Rscript code/03_loco_tol.R

message("=== Step 3 (LOCO-tol): Re-select from existing BO results ===")

source("code/00_helpers_loco_tol.R")
library(DeSurv)
library(NMF)
source("R/bo_helpers.R")
source("R/bo_helpers_loco.R")
source("R/pick_k_elbow.R")
source("R/targets_config.R")

# ── 1. DeSurv (supervised) ──────────────────────────────────────────────
desurv_bo_results <- load_precomputed("desurv_bo_results_tcgacptac")

tar_k_selection <- cache_or_compute("tar_k_selection_tcgacptac", {
  select_bo_k_best_mean(desurv_bo_results, tol = 0.005)
})
message("  LOCO k selection: k=", tar_k_selection$k_selected,
        " (best mean=", round(tar_k_selection$best_mean, 4),
        ", selected alpha=", round(tar_k_selection$selected_alpha, 4),
        ", selected cindex=", round(tar_k_selection$selected_cindex, 4), ")")

tar_params_best <- cache_or_compute("tar_params_best_tcgacptac", {
  params <- tar_k_selection$params
  if (is.null(params)) {
    params <- standardize_bo_params(desurv_bo_results$overall_best$params)
  }
  params$k <- tar_k_selection$k_selected
  params
})

tar_data <- load_precomputed("tar_data_tcgacptac")

tar_data_filtered <- cache_or_compute("tar_data_filtered_tcgacptac", {
  ngene_value <- if (!is.null(tar_params_best$ngene) && !is.na(tar_params_best$ngene)) {
    as.integer(round(tar_params_best$ngene))
  } else {
    3000L
  }
  preprocess_training_data(data = tar_data, ngene = ngene_value,
                           method_trans_train = "rank")
})

message("  DeSurv: k=", tar_params_best$k,
        ", alpha=", round(tar_params_best$alpha, 4),
        ", lambda=", round(tar_params_best$lambda, 4),
        ", nu=", round(tar_params_best$nu, 4))

# ── 2. Alpha=0 (unsupervised NMF via BO for lambda/nu) ──────────────────
desurv_bo_results_alpha0 <- load_precomputed("desurv_bo_results_alpha0_tcgacptac")

# Alpha=0 has no alpha to select — just pick best k, then best lambda/nu at that k
tar_k_selection_alpha0 <- select_bo_k_best_mean(desurv_bo_results_alpha0, tol = 0.005)

tar_params_best_alpha0 <- cache_or_compute("tar_params_best_alpha0_tcgacptac", {
  params <- tar_k_selection_alpha0$params
  if (is.null(params)) {
    params <- standardize_bo_params(desurv_bo_results_alpha0$overall_best$params)
  }
  if (!is.null(tar_k_selection_alpha0$k_selected)) {
    params$k <- tar_k_selection_alpha0$k_selected
  }
  params$alpha <- 0
  params
})

tar_data_filtered_alpha0 <- cache_or_compute("tar_data_filtered_alpha0_tcgacptac", {
  ngene_alpha0 <- if (!is.null(tar_params_best_alpha0$ngene) && !is.na(tar_params_best_alpha0$ngene)) {
    as.integer(round(tar_params_best_alpha0$ngene))
  } else {
    3000L
  }
  preprocess_training_data(data = tar_data, ngene = ngene_alpha0,
                           method_trans_train = "rank")
})

message("  Alpha=0: k=", tar_params_best_alpha0$k)

# ── 3. Standard NMF elbow k (reuse from loco/) ──────────────────────────
fit_std <- load_precomputed("fit_std_tcgacptac")
std_nmf_selected_k <- load_precomputed("std_nmf_selected_k_tcgacptac")
message("  Standard NMF elbow k=", std_nmf_selected_k)

# ── 4. BO at elbow k ────────────────────────────────────────────────────
desurv_bo_results_elbowk <- load_precomputed("desurv_bo_results_elbowk_tcgacptac")

# For fixed-k BO, just pick the best lambda/nu/alpha combination
tar_params_best_elbowk <- cache_or_compute("tar_params_best_elbowk_tcgacptac", {
  standardize_bo_params(desurv_bo_results_elbowk$overall_best$params)
})

tar_data_filtered_elbowk <- cache_or_compute("tar_data_filtered_elbowk_tcgacptac", {
  ngene_elbowk <- if (!is.null(tar_params_best_elbowk$ngene) && !is.na(tar_params_best_elbowk$ngene)) {
    as.integer(round(tar_params_best_elbowk$ngene))
  } else {
    3000L
  }
  preprocess_training_data(data = tar_data, ngene = ngene_elbowk,
                           method_trans_train = "rank")
})

message("=== Step 3 (LOCO-tol) complete ===")
