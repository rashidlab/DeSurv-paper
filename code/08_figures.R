#!/usr/bin/env Rscript
# code/08_figures.R — Generate all figure objects from pre-computed results
#
# Rebuilds figure objects stored in results/precomputed/ from model fits
# and validation results produced by earlier pipeline steps.
#
# The CV-based cutpoint analysis lives in code/08a_cutpoint_analysis.R so
# that both this script and code/08b_si_figures.R can reuse the same
# (k, alpha)-conditional cutpoint selection without recomputing it.
#
# Inputs:  model fits, validation results from steps 04-05
#          cutpoint summaries / LP stats from step 08a
# Outputs: figure objects in results/precomputed/fig_*.rds
#
# Runtime: ~2 min (cutpoint CV is now in step 08a)

message("=== Step 8: Figures ===")
source("code/00_helpers.R")
library(ggplot2)
library(ggrepel)
library(cowplot)
library(survival)
library(survminer)
library(NMF)
library(glmnet)
library(dplyr)

# Source figure-building functions
source("R/figure_targets.R")
source("R/cluster_alignment.R")
source("R/compare_models.R")
source("R/cv_grid_helpers.R")
source("R/get_top_genes.R")
source("R/plot_survival.R")
source("R/fit_cox_model.R")

# ── Load prerequisites ────────────────────────────────────────────────────
tar_fit_desurv        <- load_precomputed("tar_fit_desurv_tcgacptac")
tar_fit_desurv_elbowk <- load_precomputed("tar_fit_desurv_elbowk_tcgacptac")
fit_std_desurvk       <- load_precomputed("fit_std_desurvk_tcgacptac")
fit_std_elbowk        <- load_precomputed("fit_std_elbowk_tcgacptac")
fit_std              <- load_precomputed("fit_std_tcgacptac")
tar_data_filtered   <- load_precomputed("tar_data_filtered_tcgacptac")
data_val_filtered   <- load_precomputed("data_val_filtered_tcgacptac")
desurv_bo_results   <- load_precomputed("desurv_bo_results_tcgacptac")
tar_params_best     <- load_precomputed("tar_params_best_tcgacptac")

# ── NMF diagnostic plots (residuals, cophenetic, silhouette) ─────────────
fig_residuals <- make_nmf_metric_plot(fit_std, "residuals")
saveRDS(fig_residuals, file.path(RESULTS_DIR, "fig_residuals_tcgacptac.rds"))

fig_cophenetic <- make_nmf_metric_plot(fit_std, "cophenetic")
saveRDS(fig_cophenetic, file.path(RESULTS_DIR, "fig_cophenetic_tcgacptac.rds"))

fig_silhouette <- make_nmf_metric_plot(fit_std, "silhouette")
saveRDS(fig_silhouette, file.path(RESULTS_DIR, "fig_silhouette_tcgacptac.rds"))

# ── BO heatmap (GP-predicted C-index over k x alpha grid) ────────────────
curve <- extract_gp_curve(desurv_bo_results, tar_params_best)
fig_bo_heat <- ggplot(curve, aes(x = k, y = alpha, fill = mean)) +
    geom_tile(color = NA) +
    scale_x_continuous(breaks = seq(2, 12, by = 1)) +
    scale_fill_viridis_c(
      name = "CV C-index",
      option = "D",
      guide = guide_colorbar(barheight = unit(3, "cm"), barwidth = unit(0.4, "cm"))
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      labels = function(x) ifelse(x == 0, "0 (NMF)", as.character(x))
    ) +
    labs(x = "Factorization rank (k)", y = "Supervision strength") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.text  = element_text(color = "black")
    )
saveRDS(fig_bo_heat, file.path(RESULTS_DIR, "fig_bo_heat_tcgacptac.rds"))

# maxed bo heat
curve = extract_gp_curve_maxed(desurv_bo_results)
fig_bo_heat_maxed <- ggplot(curve, aes(x = k, y = alpha, fill = mean)) +
  geom_tile(color = NA) +
  scale_x_continuous(breaks = seq(2, 12, by = 1)) +
  scale_fill_viridis_c(
    name = "CV C-index",
    option = "D",
    guide = guide_colorbar(barheight = unit(3, "cm"), barwidth = unit(0.4, "cm"))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) ifelse(x == 0, "0 (NMF)", as.character(x))
  ) +
  labs(x = "Factorization rank (k)", y = "Supervision strength") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(color = "black")
  )
saveRDS(fig_bo_heat_maxed, file.path(RESULTS_DIR, "fig_bo_heat_maxed_tcgacptac.rds"))

# ── Variance explained vs survival contribution scatter ───────────────────
df_nmf <- build_variance_survival_df(
  X = tar_data_filtered$ex,
  scores = fit_std_desurvk$W,
  loadings = fit_std_desurvk$H,
  time = tar_data_filtered$sampInfo$time,
  event = tar_data_filtered$sampInfo$event,
  method = "NMF"
)
df_desurv <- build_variance_survival_df(
  X = tar_data_filtered$ex,
  scores = tar_fit_desurv$W,
  loadings = tar_fit_desurv$H,
  time = tar_data_filtered$sampInfo$time,
  event = tar_data_filtered$sampInfo$event,
  method = "DeSurv"
)
df_plot <- dplyr::bind_rows(df_nmf, df_desurv) |>
  dplyr::mutate(
    factor_label = dplyr::case_when(
      method == "NMF" ~ paste0("N", factor),
      method == "DeSurv" ~ paste0("D", factor),
      TRUE ~ paste0("F", factor)
    )
  )
fig_variation_explained <- ggplot(df_plot,
         aes(x = variance_explained, y = delta_loglik,
             label = factor_label, color = method)) +
    geom_point(size = 4) +
    geom_text_repel(
      size = 4, max.overlaps = Inf, box.padding = 0.6,
      point.padding = 0.4, segment.size = 0.3, force = 2
    ) +
    scale_color_manual(values = c("NMF" = "red", "DeSurv" = "blue")) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Conditional variance explained\n(semi-partial R\u00b2)",
      y = expression(atop(Delta ~ "partial log-likelihood",
                          "(full vs. k-1 factor model)")),
      color = "Method"
    ) +
    theme_classic(base_size = 10)
saveRDS(fig_variation_explained, file.path(RESULTS_DIR, "fig_variation_explained_tcgacptac.rds"))

# ── Gene overlap heatmaps ─────────────────────────────────────────────────
# Reference gene lists from combined subtype data
top_genes_path <- "data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(top_genes_path)) {
  # Fall back to the original DeSurv-paper repo
  top_genes_path <- "../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData"
}
if (!file.exists(top_genes_path)) {
  stop("Reference gene list not found. Copy data/derv/cmbSubtypes_formatted.RData ",
       "from DeSurv-paper or set the path manually.")
}
load(top_genes_path)  # loads: top_genes, colors, subtypeList, etc.

# Factor labels — only apply when k matches the original run (k=3)
desurv_k <- ncol(tar_fit_desurv$W)
std_k    <- ncol(fit_std_desurvk$W)
heatmap_factor_labels <- if (desurv_k == 3) {
  c("D1 Classical/restCAF", "D2 proCAF", "D3 Basal-like")
} else {
  paste0("D", seq_len(desurv_k))
}
heatmap_factor_labels_std <- if (std_k == 3) {
  c("N1 Classical", "N2 Exocrine", "N3 Microenviron.")
} else {
  paste0("N", seq_len(std_k))
}

# ntop for top gene extraction: fixed mode uses CONFIG$ntop_value; BO mode uses
# the BO-selected ntop from tar_params_best
ntop_value <- if (!is.null(CONFIG$ntop_value)) CONFIG$ntop_value else tar_params_best$ntop

# Compute top genes per model
tar_tops_desurv        <- get_top_genes(W = tar_fit_desurv$W, ntop = ntop_value)
tar_tops_desurv_elbowk <- get_top_genes(W = tar_fit_desurv_elbowk$W, ntop = ntop_value)
tar_tops_std_desurvk   <- get_top_genes(W = fit_std_desurvk$W, ntop = ntop_value)
tar_tops_std_elbowk    <- get_top_genes(W = fit_std_elbowk$W, ntop = ntop_value)

tar_fit_desurv_alpha0 <- load_precomputed("tar_fit_desurv_alpha0_tcgacptac")
tar_tops_desurv_alpha0 <- get_top_genes(W = tar_fit_desurv_alpha0$W, ntop = ntop_value)

# DeSurv heatmap
fig_gene_overlap_heatmap_desurv <- make_gene_overlap_heatmap(
      tar_fit_desurv, tar_tops_desurv$top_genes, top_genes,
      factor_labels = heatmap_factor_labels, title = "DeSurv", fontsize_row = 7
    )
saveRDS(fig_gene_overlap_heatmap_desurv,
        file.path(RESULTS_DIR, "fig_gene_overlap_heatmap_desurv_tcgacptac.rds"))

# Standard NMF at DeSurv k heatmap
fig_gene_overlap_heatmap_std_desurvk <- make_gene_overlap_heatmap(
      fit_std_desurvk, tar_tops_std_desurvk$top_genes, top_genes,
      factor_labels = heatmap_factor_labels_std, title = "NMF", fontsize_row = 7
    )
saveRDS(fig_gene_overlap_heatmap_std_desurvk,
        file.path(RESULTS_DIR, "fig_gene_overlap_heatmap_std_desurvk_tcgacptac.rds"))

# Standard NMF at elbow k heatmap
fig_gene_overlap_heatmap_std_elbowk <- make_gene_overlap_heatmap(
      fit_std_elbowk, tar_tops_std_elbowk$top_genes, top_genes
    )
saveRDS(fig_gene_overlap_heatmap_std_elbowk,
        file.path(RESULTS_DIR, "fig_gene_overlap_heatmap_std_elbowk_tcgacptac.rds"))

# DeSurv at elbow k heatmap
fig_gene_overlap_heatmap_desurv_elbowk <- make_gene_overlap_heatmap(
      tar_fit_desurv_elbowk, tar_tops_desurv_elbowk$top_genes, top_genes
    )
saveRDS(fig_gene_overlap_heatmap_desurv_elbowk,
        file.path(RESULTS_DIR, "fig_gene_overlap_heatmap_desurv_elbowk_tcgacptac.rds"))

# DeSurv alpha=0 heatmap
fig_gene_overlap_heatmap_desurv_alpha0 <- make_gene_overlap_heatmap(
      tar_fit_desurv_alpha0, tar_tops_desurv_alpha0$top_genes, top_genes
    )
saveRDS(fig_gene_overlap_heatmap_desurv_alpha0,
        file.path(RESULTS_DIR, "fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac.rds"))

# ── HR forest plot ────────────────────────────────────────────────────────
desurv_df <- compute_hrs(data_val_filtered, tar_fit_desurv, "DeSurv")
nmf_df    <- compute_hrs(data_val_filtered, fit_std_desurvk, "NMF")

df <- rbind(desurv_df, nmf_df)
pd <- position_dodge(width = 0.6)

df$label <- sprintf("%.2f (%.2f\u2013%.2f)", df$HR, df$lower, df$upper)
df$factor_name <- dplyr::case_when(
  df$method == "DeSurv" ~ paste0("D", df$factor),
  df$method == "NMF"    ~ paste0("N", df$factor),
  TRUE                  ~ as.character(df$factor)
)
fig_hr_forest <- ggplot(df, aes(x = HR, y = factor_name, color = dataset, group = dataset)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               linewidth = 0.5, color = "grey60") +
    geom_errorbarh(
      aes(xmin = lower, xmax = upper),
      height = 0.25, linewidth = 0.8, position = pd
    ) +
    geom_point(size = 2.8, position = pd) +
    scale_x_log10(
      breaks = c(0.5, 1, 2, 4),
      labels = c("0.5", "1", "2", "4")
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_cartesian(xlim = c(min(df$lower), max(df$upper) * 1.4)) +
    labs(x = "Hazard ratio (95% CI)", y = NULL) +
    facet_wrap(~method)
saveRDS(fig_hr_forest, file.path(RESULTS_DIR, "fig_hr_forest_tcgacptac.rds"))

# ── NMF vs DeSurv Spearman correlation heatmap ────────────────────────────
c_mat <- cor(fit_std_desurvk$W, tar_fit_desurv$W, method = "spearman")
rownames(c_mat) <- heatmap_factor_labels_std
colnames(c_mat) <- heatmap_factor_labels

ph_args <- list(
  mat = c_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 8,
  fontsize_number = 8,
  number_color = "black",
  breaks = seq(-0.5, 1, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  silent = TRUE
)

ph <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = FALSE)))
ph_grob <- ph$gtable
mat_idx  <- which(ph_grob$layout$name == "matrix")
mat_l <- ph_grob$layout$l[mat_idx]; mat_r <- ph_grob$layout$r[mat_idx]
mat_t <- ph_grob$layout$t[mat_idx]; mat_b <- ph_grob$layout$b[mat_idx]
w_pt     <- as.numeric(grid::convertUnit(ph_grob$widths,  "pt"))
h_pt     <- as.numeric(grid::convertUnit(ph_grob$heights, "pt"))
x_center <- (sum(w_pt[seq_len(mat_l - 1)]) + sum(w_pt[mat_l:mat_r]) / 2) / sum(w_pt)
y_center <- 1 - (sum(h_pt[seq_len(mat_t - 1)]) + sum(h_pt[mat_t:mat_b]) / 2) / sum(h_pt)

hm <- cowplot::ggdraw(ph_grob)
pheat <- cowplot::plot_grid(
  hm,
  cowplot::ggdraw() + cowplot::draw_label("NMF", angle = 90, size = 9, y = y_center, fontface = "bold"),
  cowplot::ggdraw() + cowplot::draw_label("DeSurv", size = 9, x = x_center, fontface = "bold"),
  NULL,
  ncol = 2, rel_widths = c(1, 0.1), rel_heights = c(1, 0.1)
)

ph_leg <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = TRUE)))
leg_idx <- which(ph_leg$gtable$layout$name == "legend")
legend_grob <- if (length(leg_idx) > 0) ph_leg$gtable$grobs[[leg_idx[1]]] else grid::nullGrob()

fig_desurv_std_correlation <- list(plot = pheat, legend = legend_grob)
saveRDS(fig_desurv_std_correlation,
        file.path(RESULTS_DIR, "fig_desurv_std_correlation_tcgacptac.rds"))


# ── ntop NMF vs DeSurv Spearman correlation heatmap ────────────────────────────
tops_de = tar_tops_desurv$top_genes
tops_std = tar_tops_std_desurvk$top_genes

c_mat <- cor(fit_std_desurvk$W[unlist(tops_std),], tar_fit_desurv$W[unlist(tops_de),], method = "spearman")
rownames(c_mat) <- heatmap_factor_labels_std
colnames(c_mat) <- heatmap_factor_labels

ph_args <- list(
  mat = c_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 8,
  fontsize_number = 8,
  number_color = "black",
  breaks = seq(-0.5, 1, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  silent = TRUE
)

ph <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = FALSE)))
ph_grob <- ph$gtable
mat_idx  <- which(ph_grob$layout$name == "matrix")
mat_l <- ph_grob$layout$l[mat_idx]; mat_r <- ph_grob$layout$r[mat_idx]
mat_t <- ph_grob$layout$t[mat_idx]; mat_b <- ph_grob$layout$b[mat_idx]
w_pt     <- as.numeric(grid::convertUnit(ph_grob$widths,  "pt"))
h_pt     <- as.numeric(grid::convertUnit(ph_grob$heights, "pt"))
x_center <- (sum(w_pt[seq_len(mat_l - 1)]) + sum(w_pt[mat_l:mat_r]) / 2) / sum(w_pt)
y_center <- 1 - (sum(h_pt[seq_len(mat_t - 1)]) + sum(h_pt[mat_t:mat_b]) / 2) / sum(h_pt)

hm <- cowplot::ggdraw(ph_grob)
pheat <- cowplot::plot_grid(
  hm,
  cowplot::ggdraw() + cowplot::draw_label("NMF", angle = 90, size = 9, y = y_center, fontface = "bold"),
  cowplot::ggdraw() + cowplot::draw_label("DeSurv", size = 9, x = x_center, fontface = "bold"),
  NULL,
  ncol = 2, rel_widths = c(1, 0.1), rel_heights = c(1, 0.1)
)

ph_leg <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = TRUE)))
leg_idx <- which(ph_leg$gtable$layout$name == "legend")
legend_grob <- if (length(leg_idx) > 0) ph_leg$gtable$grobs[[leg_idx[1]]] else grid::nullGrob()

fig_desurv_std_correlation_ntop <- list(plot = pheat, legend = legend_grob)
saveRDS(fig_desurv_std_correlation_ntop,
        file.path(RESULTS_DIR, "fig_desurv_std_correlation_tcgacptac_ntop.rds"))


# ── top 50 NMF vs DeSurv Spearman correlation heatmap ────────────────────────────
tops_de = tar_tops_desurv$top_genes[1:50,]
tops_std = tar_tops_std_desurvk$top_genes[1:50,]

c_mat <- cor(fit_std_desurvk$W[unlist(tops_std),], tar_fit_desurv$W[unlist(tops_de),], method = "spearman")
rownames(c_mat) <- heatmap_factor_labels_std
colnames(c_mat) <- heatmap_factor_labels

ph_args <- list(
  mat = c_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 8,
  fontsize_number = 8,
  number_color = "black",
  breaks = seq(-0.5, 1, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  silent = TRUE
)

ph <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = FALSE)))
ph_grob <- ph$gtable
mat_idx  <- which(ph_grob$layout$name == "matrix")
mat_l <- ph_grob$layout$l[mat_idx]; mat_r <- ph_grob$layout$r[mat_idx]
mat_t <- ph_grob$layout$t[mat_idx]; mat_b <- ph_grob$layout$b[mat_idx]
w_pt     <- as.numeric(grid::convertUnit(ph_grob$widths,  "pt"))
h_pt     <- as.numeric(grid::convertUnit(ph_grob$heights, "pt"))
x_center <- (sum(w_pt[seq_len(mat_l - 1)]) + sum(w_pt[mat_l:mat_r]) / 2) / sum(w_pt)
y_center <- 1 - (sum(h_pt[seq_len(mat_t - 1)]) + sum(h_pt[mat_t:mat_b]) / 2) / sum(h_pt)

hm <- cowplot::ggdraw(ph_grob)
pheat <- cowplot::plot_grid(
  hm,
  cowplot::ggdraw() + cowplot::draw_label("NMF", angle = 90, size = 9, y = y_center, fontface = "bold"),
  cowplot::ggdraw() + cowplot::draw_label("DeSurv", size = 9, x = x_center, fontface = "bold"),
  NULL,
  ncol = 2, rel_widths = c(1, 0.1), rel_heights = c(1, 0.1)
)

ph_leg <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = TRUE)))
leg_idx <- which(ph_leg$gtable$layout$name == "legend")
legend_grob <- if (length(leg_idx) > 0) ph_leg$gtable$grobs[[leg_idx[1]]] else grid::nullGrob()

fig_desurv_std_correlation_top50 <- list(plot = pheat, legend = legend_grob)
saveRDS(fig_desurv_std_correlation_top50,
        file.path(RESULTS_DIR, "fig_desurv_std_correlation_tcgacptac_top50.rds"))



# ── Median survival KM curves (pooled validation, log-rank optimal cutpoint) ──
# Cutpoint selection (cv_res, z-cutpoint) and LP standardisation stats
# (mean, sd) are produced by code/08a_cutpoint_analysis.R. Run that step
# before this one, or `run_pipeline.R` will do so automatically.
ntop_for_lp <- if (!is.null(CONFIG$ntop_value)) CONFIG$ntop_value else tar_params_best$ntop

lp_stats     <- load_precomputed("desurv_lp_stats_tcgacptac")
lp_stats_std <- load_precomputed("std_desurvk_lp_stats_tcgacptac")

# DeSurv KM figure
fig_median_survival_desurv <- splot_cutpoint(
  data_val_filtered, tar_fit_desurv, lp_stats, ntop = ntop_for_lp
)
saveRDS(fig_median_survival_desurv,
        file.path(RESULTS_DIR, "fig_median_survival_desurv_tcgacptac.rds"))

# Standard NMF at DeSurv k KM figure
fig_median_survival_std_desurvk <- splot_cutpoint(
  data_val_filtered, fit_std_desurvk, lp_stats_std, ntop = NULL
)
saveRDS(fig_median_survival_std_desurvk,
        file.path(RESULTS_DIR, "fig_median_survival_std_desurvk_tcgacptac.rds"))

message("=== Step 8 complete ===")
