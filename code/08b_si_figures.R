#!/usr/bin/env Rscript
# code/08b_si_figures.R — Generate SI appendix-specific figure outputs
#
# Builds the PDF figures referenced directly by paper/si_appendix.Rmd:
#
#   figures/cutpoint_curve_logrank_<bo_label>.pdf
#   figures/cutpoint_curve_cindex_<bo_label>.pdf
#   figures/km_dichot/km_val_<dataset>_logrank_<bo_label>.pdf
#   figures/km_dichot/km_val_pooled_logrank_<bo_label>.pdf
#   figures/km_dichot/subtype_overlap_pooled_logrank_<bo_label>.pdf
#
# The cutpoint summary, optimal z, and LP statistics are loaded from
# code/08a_cutpoint_analysis.R outputs so the SI figures and the main-text
# median KM curve agree on the dichotomization rule.
#
# Inputs:  cutpoint summary / lp stats from step 08a, validation data
#          and DeSurv fit from earlier steps.
# Outputs: PDF files under figures/ and figures/km_dichot/.
#
# Runtime: ~2 min

message("=== Step 8b: SI appendix figures ===")
source("code/00_helpers.R")
library(ggplot2)
library(cowplot)
library(survival)
library(survminer)
library(dplyr)
library(glmnet)

source("R/cv_grid_helpers.R")
source("R/figure_targets.R")
source("R/targets_config.R")

# ── Identifier used in output filenames ──────────────────────────────────
# Mirrors the bo_label convention from the original DeSurv-paper repo. The
# rl pipeline only ever uses tcgacptac, but expose it as a variable so the
# script can be reused if a second cohort is added.
bo_label <- "tcgacptac"

# ── Output directories ───────────────────────────────────────────────────
fig_dir <- "figures"
km_dir  <- file.path(fig_dir, "km_dichot")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(km_dir,  recursive = TRUE, showWarnings = FALSE)

# ── Load prerequisites ───────────────────────────────────────────────────
tar_fit_desurv      <- load_precomputed("tar_fit_desurv_tcgacptac")
tar_data_filtered   <- load_precomputed("tar_data_filtered_tcgacptac")
data_val_filtered   <- load_precomputed("data_val_filtered_tcgacptac")
tar_params_best     <- load_precomputed("tar_params_best_tcgacptac")

desurv_cutpoint_summary <- load_precomputed("desurv_cutpoint_summary_tcgacptac")
desurv_lp_stats         <- load_precomputed("desurv_lp_stats_tcgacptac")

ntop_for_lp <- if (!is.null(CONFIG$ntop_value)) CONFIG$ntop_value else tar_params_best$ntop
ntop_label  <- if (is.null(ntop_for_lp) || is.na(ntop_for_lp)) NA_real_ else as.numeric(ntop_for_lp)

# ── Survival validation list (with PACA_AU array+seq merged) ─────────────
# data_val_filtered from step 05 is an unnamed lapply result. The dichotomized
# KM and subtype-overlap helpers loop over names, so we restore them from
# each element's $dataname before merging the two PACA_AU platforms.
val_named <- data_val_filtered
nms <- vapply(val_named, function(x) {
  if (!is.null(x$dataname) && nzchar(x$dataname)) x$dataname
  else infer_dataset_name(x, fallback = "unknown")
}, character(1))
names(val_named) <- nms
data_val_filtered_surv <- merge_paca_au_datasets(val_named)

# ── Cutpoint selection curves (logrank + cindex) ─────────────────────────
p_logrank <- plot_cutpoint_curve_logrank(
  desurv_cutpoint_summary,
  k = tar_params_best$k,
  alpha = tar_params_best$alpha,
  ntop = ntop_label,
  optimal_z = desurv_lp_stats$optimal_z_cutpoint
)
ggplot2::ggsave(
  file.path(fig_dir, sprintf("cutpoint_curve_logrank_%s.pdf", bo_label)),
  p_logrank, width = 5, height = 4
)

p_cindex <- plot_cutpoint_curve_cindex(
  desurv_cutpoint_summary,
  k = tar_params_best$k,
  alpha = tar_params_best$alpha,
  ntop = ntop_label,
  optimal_z = desurv_lp_stats$optimal_z_cutpoint_cindex
)
ggplot2::ggsave(
  file.path(fig_dir, sprintf("cutpoint_curve_cindex_%s.pdf", bo_label)),
  p_cindex, width = 5, height = 4
)

# ── KM curves on dichotomized validation cohorts ─────────────────────────
# fit_entry mirrors the structure cv_grid_helpers expects: a fit list (W, beta)
# with metadata describing how to compute and z-standardize the LP.
fit_entry <- list(
  fit         = list(W = tar_fit_desurv$W, beta = tar_fit_desurv$beta),
  ntop        = ntop_for_lp,
  z_cutpoint  = desurv_lp_stats$optimal_z_cutpoint,
  lp_mean     = desurv_lp_stats$lp_mean,
  lp_sd       = desurv_lp_stats$lp_sd,
  k           = tar_params_best$k,
  alpha       = tar_params_best$alpha %||% NA_real_
)

for (ds_name in names(data_val_filtered_surv)) {
  p <- plot_km_validation(fit_entry, data_val_filtered_surv[[ds_name]], ds_name)
  if (is.null(p)) next
  fpath <- file.path(km_dir,
                     sprintf("km_val_%s_logrank_%s.pdf", ds_name, bo_label))
  save_ggsurvplot(p, fpath, width = 7, height = 5)
}

p_pooled <- plot_km_validation_pooled(fit_entry, data_val_filtered_surv)
if (!is.null(p_pooled)) {
  fpath <- file.path(km_dir, sprintf("km_val_pooled_logrank_%s.pdf", bo_label))
  save_ggsurvplot(p_pooled, fpath, width = 7, height = 5)
}

# ── Subtype-overlap stacked bar (pooled validation, log-rank cutpoint) ───
W    <- tar_fit_desurv$W
beta <- tar_fit_desurv$beta
z_cut   <- desurv_lp_stats$optimal_z_cutpoint
lp_mean <- desurv_lp_stats$lp_mean
lp_sd   <- desurv_lp_stats$lp_sd

pooled_rows <- list()
for (ds_name in names(data_val_filtered_surv)) {
  val_ds <- data_val_filtered_surv[[ds_name]]
  common_genes <- intersect(rownames(W), rownames(val_ds$ex))
  if (length(common_genes) < 2) next

  W_common <- W[common_genes, , drop = FALSE]
  X_val <- val_ds$ex[common_genes, , drop = FALSE]

  time_val  <- val_ds$sampInfo$time
  event_val <- val_ds$sampInfo$event
  valid_idx <- which(is.finite(time_val) & !is.na(event_val) & time_val > 0)
  if (length(valid_idx) < 2) next

  val_group <- compute_risk_group(
    W_common, beta, X_val[, valid_idx, drop = FALSE],
    ntop_for_lp, lp_mean, lp_sd, z_cut
  )
  val_si <- val_ds$sampInfo[valid_idx, ]

  if (!all(c("PurIST", "DeCAF") %in% names(val_si))) next
  ds_df <- data.frame(
    group  = val_group,
    PurIST = val_si$PurIST,
    DeCAF  = val_si$DeCAF,
    stringsAsFactors = FALSE
  )
  ds_df <- ds_df[complete.cases(ds_df), ]
  if (nrow(ds_df) > 0) pooled_rows[[ds_name]] <- ds_df
}

if (length(pooled_rows) > 0) {
  pooled_df <- do.call(rbind, pooled_rows)
  if (nrow(pooled_df) > 0) {
    p <- plot_subtype_overlap(
      pooled_df,
      dataset_label = "Pooled validation (logrank cutpoint)"
    )
    fpath <- file.path(km_dir,
                       sprintf("subtype_overlap_pooled_logrank_%s.pdf", bo_label))
    ggplot2::ggsave(fpath, p, width = 8, height = 5)
  }
}

message("=== Step 8b complete ===")
