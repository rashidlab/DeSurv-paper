#!/usr/bin/env Rscript
# code/03b_select_k.R — K-selection and parameter extraction from cached BO results
#
# Applies the 1-SE rule to select k, extracts best parameters, and
# preprocesses training data.  Requires BO results from step 03 to already
# exist in RESULTS_DIR.
#
# Inputs:  BO results from step 03 (desurv_bo_results_*.rds, fit_std_*.rds)
# Outputs: tar_k_selection_tcgacptac.rds
#          tar_params_best_tcgacptac.rds
#          tar_data_filtered_tcgacptac.rds
#          (+ alpha0 and elbow variants)

message("=== Step 3b: K-selection & params ===")
source("code/00_helpers.R")
library(DeSurv)
source("R/bo_helpers.R")
source("R/targets_config.R")

# ── Load training data ────────────────────────────────────────────────────
tar_data <- load_precomputed("tar_data_tcgacptac")

# ── 1. DeSurv k-selection ────────────────────────────────────────────────
desurv_bo_results <- load_precomputed("desurv_bo_results_tcgacptac")

tar_k_selection <- cache_or_compute("tar_k_selection_tcgacptac", {
  select_bo_k_by_cv_se(desurv_bo_results)
})
message("  1-SE rule: k_selected=", tar_k_selection$k_selected,
        ", k_best=", tar_k_selection$k_best,
        ", reason=", tar_k_selection$reason)

tar_params_best <- cache_or_compute("tar_params_best_tcgacptac", {
  if (CONFIG$param_rule == "1se") {
    params <- select_bo_params_1se(desurv_bo_results)
    message("  Param rule '1se': took full 1-SE row from BO history.")
  } else {
    params <- standardize_bo_params(desurv_bo_results$overall_best$params)
    if (!is.null(tar_k_selection$k_selected)) {
      params$k <- tar_k_selection$k_selected
    }
  }
  params
})

tar_data_filtered <- cache_or_compute("tar_data_filtered_tcgacptac", {
  ngene_value <- if (!is.null(tar_params_best$ngene) && !is.na(tar_params_best$ngene)) {
    as.integer(round(tar_params_best$ngene))
  } else {
    3000L
  }
  preprocess_training_data(data = tar_data, ngene = ngene_value,
                           method_trans_train = "rank")
})

message("  DeSurv: k=", tar_params_best$k, ", alpha=", round(tar_params_best$alpha, 3))

# ── 2. Alpha=0 k-selection ──────────────────────────────────────────────
desurv_bo_results_alpha0 <- load_precomputed("desurv_bo_results_alpha0_tcgacptac")

tar_k_selection_alpha0 <- select_bo_k_by_cv_se(desurv_bo_results_alpha0)

tar_params_best_alpha0 <- cache_or_compute("tar_params_best_alpha0_tcgacptac", {
  if (CONFIG$param_rule == "1se") {
    params <- select_bo_params_1se(desurv_bo_results_alpha0)
    message("  Alpha=0 param rule '1se': took full 1-SE row from BO history.")
  } else {
    params <- standardize_bo_params(desurv_bo_results_alpha0$overall_best$params)
    if (!is.null(tar_k_selection_alpha0$k_selected)) {
      params$k <- tar_k_selection_alpha0$k_selected
    }
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

message("  Alpha=0: k=", tar_params_best_alpha0$k)

# ── 3. Elbow-k params ───────────────────────────────────────────────────
desurv_bo_results_elbowk <- load_precomputed("desurv_bo_results_elbowk_tcgacptac")

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

message("=== Step 3b complete ===")
