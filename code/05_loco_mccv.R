#!/usr/bin/env Rscript
# code/05_loco_mccv.R — External validation using LOCO-MCCV-fitted models
#
# Identical to code/05_loco.R but reads/writes from results/precomputed/loco_mccv/

message("=== Step 5 (LOCO MCCV): External Validation ===")

# Source base loco helpers then override to loco_mccv subfolder
source("code/00_helpers_loco.R")
RESULTS_DIR  <- file.path(RESULTS_BASE, "loco_mccv")
FIGURE_DIR   <- file.path("paper/figures", "loco_mccv")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
precomputed_path <- function(name) {
  file.path(RESULTS_DIR, paste0(name, ".rds"))
}
load_precomputed <- function(name) {
  path <- precomputed_path(name)
  if (file.exists(path)) return(readRDS(path))
  shared_base_files <- c("tar_data_tcgacptac", "sim_figs_by_scenario")
  if (name %in% shared_base_files) {
    base_path <- file.path(RESULTS_BASE, paste0(name, ".rds"))
    if (file.exists(base_path)) {
      message("  Loading from base directory: ", name)
      return(readRDS(base_path))
    }
  }
  stop("Required result not found: ", path,
       "\nRun earlier pipeline steps first.")
}

library(DeSurv)
library(survival)
source("R/load_data.R")
source("R/load_data_internal.R")
source("R/targets_config.R")
source("R/predict_validation_scores.R")
source("R/cv_grid_helpers.R")

# ── Validation cohorts ────────────────────────────────────────────────────
VAL_DATASETS <- c("Dijk", "Moffitt_GEO_array", "PACA_AU_array",
                   "PACA_AU_seq", "Puleo_array")

# ── Load prerequisites ────────────────────────────────────────────────────
tar_data_filtered        <- load_precomputed("tar_data_filtered_tcgacptac")
tar_data_filtered_elbowk <- load_precomputed("tar_data_filtered_elbowk_tcgacptac")
tar_fit_desurv           <- load_precomputed("tar_fit_desurv_tcgacptac")
fit_std_desurvk          <- load_precomputed("fit_std_desurvk_tcgacptac")
fit_std_elbowk           <- load_precomputed("fit_std_elbowk_tcgacptac")
tar_params_best          <- load_precomputed("tar_params_best_tcgacptac")

# ntop for validation: NULL means all genes
ntop_value <- if (!is.null(tar_params_best$ntop) && !is.na(tar_params_best$ntop)) {
  as.integer(round(tar_params_best$ntop))
} else {
  NULL
}

# ── Helper: preprocess and filter validation data ────────────────────────
load_val_datasets <- function(val_names, training_data, method_trans_train = "rank") {
  lapply(val_names, function(ds) {
    raw <- load_data_internal(dataname = ds)
    preprocess_validation_data(
      dataset = raw,
      genes = rownames(training_data$ex),
      method_trans_train = method_trans_train,
      dataname = ds,
      zero_fill_missing = TRUE
    )
  })
}

# ── Helper: compute C-index for a model across validation cohorts ────────
compute_val_cindex <- function(fit, data_val_list, ntop = NULL) {
  results <- lapply(data_val_list, function(dv) {
    if (!is.null(ntop)) {
      top_res <- DeSurv::desurv_get_top_genes(fit$W, ntop)
      top_genes <- unique(unlist(top_res$top_genes))
      W_sub <- fit$W[top_genes, , drop = FALSE]
      XtW <- t(dv$ex[top_genes, , drop = FALSE]) %*% W_sub
      lp <- as.numeric(XtW %*% fit$beta)
    } else {
      lp <- as.numeric(t(dv$ex) %*% fit$W %*% fit$beta)
    }
    keep <- dv$sampInfo$keep == 1
    ci <- tryCatch(
      survival::concordance(
        survival::Surv(dv$sampInfo$time[keep], dv$sampInfo$event[keep]) ~ lp[keep],
        reverse = TRUE
      )$concordance,
      error = function(e) NA_real_
    )
    data.frame(dataset = dv$dataname, cindex = ci, stringsAsFactors = FALSE)
  })
  do.call(rbind, results)
}

# ── Helper: extract latent scores ────────────────────────────────────────
extract_val_latent <- function(fit, data_val_list, ntop = NULL) {
  lapply(data_val_list, function(dv) {
    keep <- dv$sampInfo$keep == 1
    X_keep <- dv$ex[, keep, drop = FALSE]
    Z <- t(X_keep) %*% fit$W
    lp <- as.numeric(Z %*% fit$beta)
    list(
      dataset = dv$dataname,
      latent = Z,
      risk_score = lp,
      survival = dv$sampInfo[keep, c("time", "event"), drop = FALSE]
    )
  })
}

# ── Preprocess validation data ────────────────────────────────────────────
data_val_filtered <- cache_or_compute("data_val_filtered_tcgacptac", {
  load_val_datasets(VAL_DATASETS, tar_data_filtered)
})

data_val_filtered_elbowk <- data_val_filtered

# ── DeSurv validation ────────────────────────────────────────────────────
val_cindex_desurv <- cache_or_compute("val_cindex_desurv_tcgacptac", {
  compute_val_cindex(tar_fit_desurv, data_val_filtered, ntop = ntop_value)
})

val_latent_desurv <- cache_or_compute("val_latent_desurv_tcgacptac", {
  extract_val_latent(tar_fit_desurv, data_val_filtered, ntop = ntop_value)
})

# ── Standard NMF at DeSurv k validation ──────────────────────────────────
val_cindex_std_desurvk <- cache_or_compute("val_cindex_std_desurvk_tcgacptac", {
  compute_val_cindex(fit_std_desurvk, data_val_filtered, ntop = ntop_value)
})

val_latent_std_desurvk <- cache_or_compute("val_latent_std_desurvk_tcgacptac", {
  extract_val_latent(fit_std_desurvk, data_val_filtered)
})

# ── Standard NMF at elbow k validation ───────────────────────────────────
val_cindex_std_elbowk <- cache_or_compute("val_cindex_std_elbowk_tcgacptac", {
  compute_val_cindex(fit_std_elbowk, data_val_filtered_elbowk, ntop = ntop_value)
})

# ── Alpha=0 validation ──────────────────────────────────────────────────
tar_fit_desurv_alpha0    <- load_precomputed("tar_fit_desurv_alpha0_tcgacptac")
tar_data_filtered_alpha0 <- load_precomputed("tar_data_filtered_alpha0_tcgacptac")
data_val_filtered_alpha0 <- data_val_filtered

val_cindex_desurv_alpha0 <- cache_or_compute("val_cindex_desurv_alpha0_tcgacptac", {
  compute_val_cindex(tar_fit_desurv_alpha0, data_val_filtered_alpha0, ntop = ntop_value)
})

val_latent_desurv_alpha0 <- cache_or_compute("val_latent_desurv_alpha0_tcgacptac", {
  extract_val_latent(tar_fit_desurv_alpha0, data_val_filtered_alpha0, ntop = ntop_value)
})

message("=== Step 5 (LOCO MCCV) complete ===")
