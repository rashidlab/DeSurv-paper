#!/usr/bin/env Rscript
# code/09b_si_figures.R — Generate SI appendix-specific figure outputs
#
# Builds the PDF figures referenced directly by paper/si_appendix.Rmd:
#
#   figures/cutpoint_curve_logrank_<bo_label>.pdf
#   figures/cutpoint_curve_cindex_<bo_label>.pdf
#   figures/km_val_<dataset>_logrank_<bo_label>.pdf
#   figures/km_val_pooled_logrank_<bo_label>.pdf
#   figures/subtype_overlap_pooled_logrank_<bo_label>.pdf
#
# The cutpoint summary, optimal z, and LP statistics are loaded from
# code/08_cutpoint_analysis.R outputs so the SI figures and the main-text
# median KM curve agree on the dichotomization rule.
#
# Inputs:  cutpoint summary / lp stats from step 08, validation data
#          and DeSurv fit from earlier steps.
# Outputs: PDF files under figures/.
#
# Runtime: ~2 min

message("=== Step 9b: SI appendix figures ===")
source("code/00_helpers.R")
library(ggplot2)
library(cowplot)
library(survival)
library(survminer)
library(dplyr)
library(glmnet)

source("R/cv_grid_helpers.R")
source("R/preprocess_helpers.R")

# ── Identifier used in output filenames ──────────────────────────────────
# Mirrors the bo_label convention from the original DeSurv-paper repo. The
# rl pipeline only ever uses tcgacptac, but expose it as a variable so the
# script can be reused if a second cohort is added.
bo_label <- "tcgacptac"

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

ntop_for_lp <- tar_params_best$ntop
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
  file.path(FIGURE_DIR, sprintf("cutpoint_curve_logrank_%s.pdf", bo_label)),
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
  file.path(FIGURE_DIR, sprintf("cutpoint_curve_cindex_%s.pdf", bo_label)),
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
  fpath <- file.path(FIGURE_DIR,
                     sprintf("km_val_%s_logrank_%s.pdf", ds_name, bo_label))
  save_ggsurvplot(p, fpath, width = 7, height = 5)
}

p_pooled <- plot_km_validation_pooled(fit_entry, data_val_filtered_surv)
if (!is.null(p_pooled)) {
  fpath <- file.path(FIGURE_DIR, sprintf("km_val_pooled_logrank_%s.pdf", bo_label))
  save_ggsurvplot(p_pooled, fpath, width = 7, height = 5)
}

# ── Subtype-overlap stacked bar ───────────────────────────────────────────
collect_subtype_overlap <- function(W, beta, ntop, lp_mean, lp_sd, z_cut, data_val) {
  pooled_rows <- list()
  for (ds_name in names(data_val)) {
    val_ds <- data_val[[ds_name]]
    common_genes <- intersect(rownames(W), rownames(val_ds$ex))
    if (length(common_genes) < 2) next
    W_common <- W[common_genes, , drop = FALSE]
    X_val    <- val_ds$ex[common_genes, , drop = FALSE]
    time_val  <- val_ds$sampInfo$time
    event_val <- val_ds$sampInfo$event
    valid_idx <- which(is.finite(time_val) & !is.na(event_val) & time_val > 0)
    if (length(valid_idx) < 2) next
    val_group <- compute_risk_group(
      W_common, beta, X_val[, valid_idx, drop = FALSE],
      ntop, lp_mean, lp_sd, z_cut
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
  if (length(pooled_rows) > 0) do.call(rbind, pooled_rows) else NULL
}

# DeSurv subtype overlap
pooled_df <- collect_subtype_overlap(
  tar_fit_desurv$W, tar_fit_desurv$beta, ntop_for_lp,
  desurv_lp_stats$lp_mean, desurv_lp_stats$lp_sd,
  desurv_lp_stats$optimal_z_cutpoint,
  data_val_filtered_surv
)
if (!is.null(pooled_df) && nrow(pooled_df) > 0) {
  p <- plot_subtype_overlap(pooled_df,
    dataset_label = "Pooled validation (logrank cutpoint)")
  ggplot2::ggsave(
    file.path(FIGURE_DIR, sprintf("subtype_overlap_pooled_logrank_%s.pdf", bo_label)),
    p, width = 8, height = 5)
}

# Standard NMF subtype overlap
pooled_df <- collect_subtype_overlap(
  fit_std_desurvk$W, fit_std_desurvk$beta, NULL,
  std_desurvk_lp_stats$lp_mean, std_desurvk_lp_stats$lp_sd,
  std_desurvk_lp_stats$optimal_z_cutpoint,
  data_val_filtered_surv
)
if (!is.null(pooled_df) && nrow(pooled_df) > 0) {
  p <- plot_subtype_overlap(pooled_df,
    dataset_label = "Standard NMF k=3, pooled validation (logrank cutpoint)")
  ggplot2::ggsave(
    file.path(FIGURE_DIR, sprintf("subtype_overlap_pooled_logrank_std_%s.pdf", bo_label)),
    p, width = 8, height = 5)
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

legend_s4 <- cowplot::get_legend(
  set_fig_font(fig_res, 10) + ggplot2::theme(legend.position = "bottom")
)

ggsave(
  file.path(FIGURE_DIR, "si_fig_nmf_diagnostics_tcgacptac.pdf"),
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
  file.path(FIGURE_DIR, sprintf("si_fig_cutpoint_km_%s.pdf", bo_label)),
  cowplot::plot_grid(
    cowplot::plot_grid(
      read_pdf_grob(file.path(FIGURE_DIR, sprintf("cutpoint_curve_logrank_%s.pdf", bo_label))),
      read_pdf_grob(file.path(FIGURE_DIR,  sprintf("km_val_pooled_logrank_%s.pdf", bo_label))),
      ncol = 2, labels = c("A", "B"), label_size = 14),
    cowplot::plot_grid(
      read_pdf_grob(file.path(FIGURE_DIR, sprintf("km_val_Dijk_logrank_%s.pdf", bo_label))),
      read_pdf_grob(file.path(FIGURE_DIR, sprintf("km_val_Moffitt_GEO_array_logrank_%s.pdf", bo_label))),
      ncol = 2, labels = c("C", "D"), label_size = 14),
    cowplot::plot_grid(
      read_pdf_grob(file.path(FIGURE_DIR, sprintf("km_val_PACA_AU_logrank_%s.pdf", bo_label))),
      read_pdf_grob(file.path(FIGURE_DIR, sprintf("km_val_Puleo_array_logrank_%s.pdf", bo_label))),
      ncol = 2, labels = c("E", "F"), label_size = 14),
    nrow = 3, rel_heights = c(1, 1, 1)
  ),
  width = 10, height = 11
)
message(sprintf("Saved si_fig_cutpoint_km_%s.pdf", bo_label))

# ── SI S6: Subtype overlap composite ─────────────────────────────────────────
img_so_desurv <- magick::image_read(
  file.path(FIGURE_DIR, sprintf("subtype_overlap_pooled_logrank_%s.pdf", bo_label)), density = 300)
img_so_nmf <- magick::image_read(
  file.path(FIGURE_DIR, sprintf("subtype_overlap_pooled_logrank_std_%s.pdf", bo_label)), density = 300)
ggsave(
  file.path(FIGURE_DIR, sprintf("si_fig_subtype_overlap_%s.pdf", bo_label)),
  cowplot::ggdraw() + cowplot::draw_image(
    magick::image_trim(magick::image_append(c(img_so_desurv, img_so_nmf), stack = TRUE))),
  width = 8, height = 10
)
message(sprintf("Saved si_fig_subtype_overlap_%s.pdf", bo_label))

# ── SI S7: DeSurv gene overlap heatmap at k=5 ────────────────────────────────
fig_desurv_elbowk <- load_precomputed("fig_gene_overlap_heatmap_desurv_elbowk_tcgacptac")
ggsave(
  file.path(FIGURE_DIR, sprintf("si_fig_desurv_k5_heatmap_%s.pdf", bo_label)),
  cowplot::plot_grid(fig_desurv_elbowk$plot, cowplot::ggdraw(fig_desurv_elbowk$legend),
                     ncol = 2, rel_widths = c(4, 1)),
  width = 6, height = 5
)
message(sprintf("Saved si_fig_desurv_k5_heatmap_%s.pdf", bo_label))

# ── SI S8: NMF gene overlap heatmap at k=5 ───────────────────────────────────
fig_nmf_elbowk <- load_precomputed("fig_gene_overlap_heatmap_std_elbowk_tcgacptac")
ggsave(
  file.path(FIGURE_DIR, sprintf("si_fig_nmf_k5_heatmap_%s.pdf", bo_label)),
  cowplot::plot_grid(fig_nmf_elbowk$plot, cowplot::ggdraw(fig_nmf_elbowk$legend),
                     ncol = 2, rel_widths = c(4, 1)),
  width = 6, height = 5
)
message(sprintf("Saved si_fig_nmf_k5_heatmap_%s.pdf", bo_label))

# ── SI S9: NMF gene overlap heatmap at k=7 ───────────────────────────────────
fig_nmf_alpha0 <- load_precomputed("fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac")
ggsave(
  file.path(FIGURE_DIR, sprintf("si_fig_nmf_k7_heatmap_%s.pdf", bo_label)),
  cowplot::plot_grid(fig_nmf_alpha0$plot, cowplot::ggdraw(fig_nmf_alpha0$legend),
                     ncol = 2, rel_widths = c(4, 1)),
  width = 6, height = 5
)
message(sprintf("Saved si_fig_nmf_k7_heatmap_%s.pdf", bo_label))

# ── SI S10: Variance vs survival at k=5 ──────────────────────────────────────
fit_std_elbowk  <- load_precomputed("fit_std_elbowk_tcgacptac")
tar_data_elbowk <- load_precomputed("tar_data_filtered_elbowk_tcgacptac")

df_varsurvk5 <- build_variance_survival_df(
  X        = tar_data_elbowk$ex,
  scores   = fit_std_elbowk$W,
  loadings = fit_std_elbowk$H,
  time     = tar_data_elbowk$sampInfo$time,
  event    = tar_data_elbowk$sampInfo$event,
  method   = "NMF"
) |> dplyr::mutate(factor_label = paste0("N", factor))

ggsave(
  file.path(FIGURE_DIR, sprintf("si_fig_varsurvival_k5_%s.pdf", bo_label)),
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

# ── CV grid C-index figure (step 06 output) ──────────────────────────────────
cv_summary_path <- file.path(CV_GRID_DIR, "cv_grid_summary.csv")
cv_val_path     <- file.path(CV_GRID_DIR, "cv_grid_val_summary.csv")
cv_alpha_path   <- file.path(CV_GRID_DIR, "cv_grid_best_alpha.csv")

if (all(file.exists(cv_summary_path, cv_val_path, cv_alpha_path))) {
  message("Generating cv_cindex_by_k_primary.pdf...")
  cv_summary  <- read.csv(cv_summary_path,  stringsAsFactors = FALSE)
  val_summary <- read.csv(cv_val_path,      stringsAsFactors = FALSE)
  best_alpha  <- read.csv(cv_alpha_path,    stringsAsFactors = FALSE)

  CV_GRID_LAMBDA <- 0.349
  CV_GRID_NU     <- 0.056
  configs <- data.frame(
    ntop   = c(270L, NA_integer_),
    lambda = c(CV_GRID_LAMBDA, CV_GRID_LAMBDA),
    nu     = c(CV_GRID_NU, CV_GRID_NU),
    label  = c(
      sprintf("ntop270_lam%.3f_nu%.3f", CV_GRID_LAMBDA, CV_GRID_NU),
      sprintf("ntopALL_lam%.3f_nu%.3f", CV_GRID_LAMBDA, CV_GRID_NU)
    ),
    stringsAsFactors = FALSE
  )

  plots <- plot_cindex_by_k(
    cv_grid_summary     = cv_summary,
    cv_grid_best_alpha  = best_alpha,
    cv_grid_val_summary = val_summary,
    configs             = configs
  )

  for (i in seq_along(plots)) {
    ggplot2::ggsave(
      file.path(FIGURE_DIR, sprintf("cv_cindex_by_k_%s.pdf", configs$label[i])),
      plots[[i]], width = 10, height = 5
    )
  }

  gt           <- ggplot2::ggplotGrob(plots[[1]])
  legend_grob  <- gt$grobs[[which(vapply(gt$grobs, `[[`, character(1), "name") == "guide-box")]]
  plots_no_leg <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
  panel_grid   <- cowplot::plot_grid(plotlist = plots_no_leg, nrow = 2,
                                     labels = c("A", "B"), label_size = 12)
  combined     <- cowplot::plot_grid(panel_grid, legend_grob,
                                     ncol = 1, rel_heights = c(1, 0.05))
  ggplot2::ggsave(file.path(FIGURE_DIR, "cv_cindex_by_k_primary.pdf"),
                  combined, width = 10, height = 10)
  message("Saved cv_cindex_by_k_primary.pdf")
} else {
  message("Skipping cv_cindex_by_k_primary.pdf: cv_grid results not found in ", CV_GRID_DIR)
  message("  Run 'make cv-grid' (or sbatch slurm/run_cv_grid.sh on HPC) to generate them.")
}

message("=== Step 9b complete ===")
