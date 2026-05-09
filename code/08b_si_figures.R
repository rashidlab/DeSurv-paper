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
fit_std_desurvk      <- load_precomputed("fit_std_desurvk_tcgacptac")
tar_data_filtered   <- load_precomputed("tar_data_filtered_tcgacptac")
data_val_filtered   <- load_precomputed("data_val_filtered_tcgacptac")
tar_params_best     <- load_precomputed("tar_params_best_tcgacptac")

desurv_cutpoint_summary <- load_precomputed("desurv_cutpoint_summary_tcgacptac")
desurv_lp_stats         <- load_precomputed("desurv_lp_stats_tcgacptac")
std_desurvk_cutpoint_summary <- load_precomputed("std_desurvk_cutpoint_summary_tcgacptac")
std_desurvk_lp_stats         <- load_precomputed("std_desurvk_lp_stats_tcgacptac")

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


# ── Subtype-overlap stacked bar std nmf k=3 ───
W    <- fit_std_desurvk$W
beta <- fit_std_desurvk$beta
z_cut   <- std_desurvk_lp_stats$optimal_z_cutpoint
lp_mean <- std_desurvk_lp_stats$lp_mean
lp_sd   <- std_desurvk_lp_stats$lp_sd

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
    NULL, lp_mean, lp_sd, z_cut
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
      dataset_label = "Standard NMF k=3, pooled validation (logrank cutpoint)"
    )
    fpath <- file.path(km_dir,
                       sprintf("subtype_overlap_pooled_logrank_std_%s.pdf", bo_label))
    ggplot2::ggsave(fpath, p, width = 8, height = 5)
  }
}
# ── SI complete figure PDFs ─────────────────────────────────────────────────
library(NMF)
library(magick)
library(ggrepel)

set_fig_font <- function(plot_obj, size = 10) {
  if (inherits(plot_obj, "ggplot")) plot_obj + ggplot2::theme(text = ggplot2::element_text(size = size)) else plot_obj
}

# ── SI S4: NMF rank selection diagnostics ─────────────────────────────────
fig_res  <- load_precomputed("fig_residuals_tcgacptac")
fig_coph <- load_precomputed("fig_cophenetic_tcgacptac")
fig_sil  <- load_precomputed("fig_silhouette_tcgacptac")
fit_std_loaded <- load_precomputed("fit_std_tcgacptac")

legend_plot_s4 <- plot(fit_std_loaded) + ggplot2::theme(legend.position = "bottom")
legend_grob_s4 <- ggplotGrob(legend_plot_s4)
legend_s4 <- legend_grob_s4$grobs[
  sapply(legend_grob_s4$grobs, function(x) x$name) == "guide-box"][[1]]

ggsave(
  file.path(fig_dir, "si_fig_nmf_diagnostics_tcgacptac.pdf"),
  cowplot::plot_grid(
    cowplot::plot_grid(
      set_fig_font(fig_res,  10) + theme(legend.position = "none"),
      set_fig_font(fig_coph, 10) + theme(legend.position = "none"),
      set_fig_font(fig_sil,  10) + theme(legend.position = "none"),
      ncol = 3, labels = c("A", "B", "C")
    ),
    legend_s4, nrow = 2, rel_heights = c(1, 0.15)
  ),
  width = 6, height = 3.5
)
message("Saved si_fig_nmf_diagnostics_tcgacptac.pdf")

# ── SI S5: Cutpoint selection and KM validation (composite) ──────────────────
read_pdf_grob <- function(path) {
  cowplot::ggdraw() + cowplot::draw_image(magick::image_read(path, density = 300))
}
ggsave(
  file.path(fig_dir, sprintf("si_fig_cutpoint_km_%s.pdf", bo_label)),
  cowplot::plot_grid(
    cowplot::plot_grid(
      read_pdf_grob(file.path(fig_dir, sprintf("cutpoint_curve_logrank_%s.pdf", bo_label))),
      read_pdf_grob(file.path(km_dir,  sprintf("km_val_pooled_logrank_%s.pdf", bo_label))),
      ncol = 2, labels = c("A", "B"), label_size = 14),
    cowplot::plot_grid(
      read_pdf_grob(file.path(km_dir, sprintf("km_val_Dijk_logrank_%s.pdf", bo_label))),
      read_pdf_grob(file.path(km_dir, sprintf("km_val_Moffitt_GEO_array_logrank_%s.pdf", bo_label))),
      ncol = 2, labels = c("C", "D"), label_size = 14),
    cowplot::plot_grid(
      read_pdf_grob(file.path(km_dir, sprintf("km_val_PACA_AU_logrank_%s.pdf", bo_label))),
      read_pdf_grob(file.path(km_dir, sprintf("km_val_Puleo_array_logrank_%s.pdf", bo_label))),
      ncol = 2, labels = c("E", "F"), label_size = 14),
    nrow = 3, rel_heights = c(1, 1, 1)
  ),
  width = 10, height = 11
)
message(sprintf("Saved si_fig_cutpoint_km_%s.pdf", bo_label))

# ── SI S6: Subtype overlap composite ─────────────────────────────────────────
img_so_desurv <- magick::image_read(
  file.path(km_dir, sprintf("subtype_overlap_pooled_logrank_%s.pdf", bo_label)), density = 300)
img_so_nmf <- magick::image_read(
  file.path(km_dir, sprintf("subtype_overlap_pooled_logrank_std_%s.pdf", bo_label)), density = 300)
ggsave(
  file.path(fig_dir, sprintf("si_fig_subtype_overlap_%s.pdf", bo_label)),
  cowplot::ggdraw() + cowplot::draw_image(
    magick::image_trim(magick::image_append(c(img_so_desurv, img_so_nmf), stack = TRUE))),
  width = 8, height = 10
)
message(sprintf("Saved si_fig_subtype_overlap_%s.pdf", bo_label))

# ── SI S7: DeSurv gene overlap heatmap at k=5 ────────────────────────────────
fig_desurv_elbowk <- load_precomputed("fig_gene_overlap_heatmap_desurv_elbowk_tcgacptac")
ggsave(
  file.path(fig_dir, sprintf("si_fig_desurv_k5_heatmap_%s.pdf", bo_label)),
  cowplot::plot_grid(fig_desurv_elbowk$plot, cowplot::ggdraw(fig_desurv_elbowk$legend),
                     ncol = 2, rel_widths = c(4, 1)),
  width = 6, height = 5
)
message(sprintf("Saved si_fig_desurv_k5_heatmap_%s.pdf", bo_label))

# ── SI S8: NMF gene overlap heatmap at k=5 ───────────────────────────────────
fig_nmf_elbowk <- load_precomputed("fig_gene_overlap_heatmap_std_elbowk_tcgacptac")
ggsave(
  file.path(fig_dir, sprintf("si_fig_nmf_k5_heatmap_%s.pdf", bo_label)),
  cowplot::plot_grid(fig_nmf_elbowk$plot, cowplot::ggdraw(fig_nmf_elbowk$legend),
                     ncol = 2, rel_widths = c(4, 1)),
  width = 6, height = 5
)
message(sprintf("Saved si_fig_nmf_k5_heatmap_%s.pdf", bo_label))

# ── SI S9: NMF gene overlap heatmap at k=7 ───────────────────────────────────
fig_nmf_alpha0 <- load_precomputed("fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac")
ggsave(
  file.path(fig_dir, sprintf("si_fig_nmf_k7_heatmap_%s.pdf", bo_label)),
  cowplot::plot_grid(fig_nmf_alpha0$plot, cowplot::ggdraw(fig_nmf_alpha0$legend),
                     ncol = 2, rel_widths = c(4, 1)),
  width = 6, height = 5
)
message(sprintf("Saved si_fig_nmf_k7_heatmap_%s.pdf", bo_label))

# ── SI S10: Variance vs survival at k=5 ──────────────────────────────────────
fit_std_elbowk  <- load_precomputed("fit_std_elbowk_tcgacptac")
tar_data_elbowk <- load_precomputed("tar_data_filtered_elbowk_tcgacptac")

compute_var_expl <- function(X, W, H) {
  X <- as.matrix(X); W <- as.matrix(W); H <- as.matrix(H)
  total_ss <- sum(X^2, na.rm = TRUE)
  rss_full <- sum((X - W %*% H)^2, na.rm = TRUE)
  tibble::tibble(
    factor = seq_len(ncol(W)),
    variance_explained = purrr::map_dbl(seq_len(ncol(W)), function(j) {
      rss_mj <- sum((X - W[, -j, drop = FALSE] %*% H[-j, , drop = FALSE])^2, na.rm = TRUE)
      (rss_mj - rss_full) / total_ss
    })
  )
}
compute_surv_expl <- function(X, W, time, event) {
  keep    <- intersect(rownames(X), rownames(W))
  XtW     <- t(X[keep, ]) %*% W[keep, ]
  ll_full <- coxph(Surv(time, event) ~ XtW)$loglik[2]
  tibble::tibble(
    factor = seq_len(ncol(XtW)),
    delta_loglik = purrr::map_dbl(seq_len(ncol(XtW)), function(j) {
      XtW_mj  <- XtW[, -j, drop = FALSE]
      red_fit <- if (ncol(XtW_mj) == 0) coxph(Surv(time, event) ~ 1) else coxph(Surv(time, event) ~ XtW_mj)
      ll_full - red_fit$loglik[2]
    })
  )
}

df_varsurvk5 <- dplyr::left_join(
  compute_var_expl(tar_data_elbowk$ex, fit_std_elbowk$W, fit_std_elbowk$H),
  compute_surv_expl(tar_data_elbowk$ex, fit_std_elbowk$W,
                    tar_data_elbowk$sampInfo$time, tar_data_elbowk$sampInfo$event),
  by = "factor"
) |> dplyr::mutate(factor_label = paste0("N", factor))

ggsave(
  file.path(fig_dir, sprintf("si_fig_varsurvival_k5_%s.pdf", bo_label)),
  ggplot2::ggplot(df_varsurvk5, ggplot2::aes(x = variance_explained, y = delta_loglik,
                                              label = factor_label)) +
    ggplot2::geom_point(size = 4, color = "red") +
    ggrepel::geom_text_repel(size = 4, max.overlaps = Inf, box.padding = 0.4,
                              point.padding = 0.3, segment.size = 0.3) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = "Conditional variance explained\n(semi-partial R²)",
      y = expression(atop(Delta ~ "partial log-likelihood", "(full vs. k-1 factor model)"))
    ) +
    ggplot2::theme_classic(base_size = 10),
  width = 4.5, height = 3.5
)
message(sprintf("Saved si_fig_varsurvival_k5_%s.pdf", bo_label))

message("=== Step 8b complete ===")
