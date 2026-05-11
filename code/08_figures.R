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


# Standard NMF at DeSurv k heatmap
fig_gene_overlap_heatmap_std_desurvk <- make_gene_overlap_heatmap(
      fit_std_desurvk, tar_tops_std_desurvk$top_genes, top_genes,
      factor_labels = heatmap_factor_labels_std, title = "NMF", fontsize_row = 7
    )


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
mat_idx <- which(ph_grob$layout$name == "matrix")
mat_l   <- ph_grob$layout$l[mat_idx]
mat_r   <- ph_grob$layout$r[mat_idx]
mat_t   <- ph_grob$layout$t[mat_idx]
mat_b   <- ph_grob$layout$b[mat_idx]

ph_grob <- gtable::gtable_add_cols(ph_grob, grid::unit(12, "pt"), pos = 0)
ph_grob <- gtable::gtable_add_grob(ph_grob,
  grid::textGrob("NMF", rot = 90, gp = grid::gpar(fontface = "bold", fontsize = 9)),
  t = mat_t, b = mat_b, l = 1, r = 1, name = "nmf-label")
ph_grob <- gtable::gtable_add_rows(ph_grob, grid::unit(12, "pt"), pos = 0)
ph_grob <- gtable::gtable_add_grob(ph_grob,
  grid::textGrob("DeSurv", gp = grid::gpar(fontface = "bold", fontsize = 9)),
  t = 1, b = 1, l = mat_l + 1, r = mat_r + 1, name = "desurv-label")

pheat <- cowplot::ggdraw(ph_grob)

ph_leg <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = TRUE)))
leg_idx <- which(ph_leg$gtable$layout$name == "legend")
legend_grob <- if (length(leg_idx) > 0) ph_leg$gtable$grobs[[leg_idx[1]]] else grid::nullGrob()

fig_desurv_std_correlation <- list(plot = pheat, legend = legend_grob)



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
mat_idx <- which(ph_grob$layout$name == "matrix")
mat_l   <- ph_grob$layout$l[mat_idx]
mat_r   <- ph_grob$layout$r[mat_idx]
mat_t   <- ph_grob$layout$t[mat_idx]
mat_b   <- ph_grob$layout$b[mat_idx]

ph_grob <- gtable::gtable_add_cols(ph_grob, grid::unit(12, "pt"), pos = 0)
ph_grob <- gtable::gtable_add_grob(ph_grob,
  grid::textGrob("NMF", rot = 90, gp = grid::gpar(fontface = "bold", fontsize = 9)),
  t = mat_t, b = mat_b, l = 1, r = 1, name = "nmf-label")
ph_grob <- gtable::gtable_add_rows(ph_grob, grid::unit(12, "pt"), pos = 0)
ph_grob <- gtable::gtable_add_grob(ph_grob,
  grid::textGrob("DeSurv", gp = grid::gpar(fontface = "bold", fontsize = 9)),
  t = 1, b = 1, l = mat_l + 1, r = mat_r + 1, name = "desurv-label")

pheat <- cowplot::ggdraw(ph_grob)

ph_leg <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = TRUE)))
leg_idx <- which(ph_leg$gtable$layout$name == "legend")
legend_grob <- if (length(leg_idx) > 0) ph_leg$gtable$grobs[[leg_idx[1]]] else grid::nullGrob()

fig_desurv_std_correlation_ntop <- list(plot = pheat, legend = legend_grob)



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
mat_idx <- which(ph_grob$layout$name == "matrix")
mat_l   <- ph_grob$layout$l[mat_idx]
mat_r   <- ph_grob$layout$r[mat_idx]
mat_t   <- ph_grob$layout$t[mat_idx]
mat_b   <- ph_grob$layout$b[mat_idx]

ph_grob <- gtable::gtable_add_cols(ph_grob, grid::unit(12, "pt"), pos = 0)
ph_grob <- gtable::gtable_add_grob(ph_grob,
  grid::textGrob("NMF", rot = 90, gp = grid::gpar(fontface = "bold", fontsize = 9)),
  t = mat_t, b = mat_b, l = 1, r = 1, name = "nmf-label")
ph_grob <- gtable::gtable_add_rows(ph_grob, grid::unit(12, "pt"), pos = 0)
ph_grob <- gtable::gtable_add_grob(ph_grob,
  grid::textGrob("DeSurv", gp = grid::gpar(fontface = "bold", fontsize = 9)),
  t = 1, b = 1, l = mat_l + 1, r = mat_r + 1, name = "desurv-label")

pheat <- cowplot::ggdraw(ph_grob)

ph_leg <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = TRUE)))
leg_idx <- which(ph_leg$gtable$layout$name == "legend")
legend_grob <- if (length(leg_idx) > 0) ph_leg$gtable$grobs[[leg_idx[1]]] else grid::nullGrob()

fig_desurv_std_correlation_top50 <- list(plot = pheat, legend = legend_grob)




# ── Median survival KM curves (pooled validation, log-rank optimal cutpoint) ──
# Cutpoint selection (cv_res, z-cutpoint) and LP standardisation stats
# (mean, sd) are produced by code/08a_cutpoint_analysis.R. Run that step
# before this one, or `run_pipeline.R` will do so automatically.
ntop_for_lp <- if (!is.null(CONFIG$ntop_value)) CONFIG$ntop_value else tar_params_best$ntop

lp_stats     <- load_precomputed("desurv_lp_stats_tcgacptac")
lp_stats_std <- load_precomputed("std_desurvk_lp_stats_tcgacptac")

# DeSurv KM figure
fig_median_survival_desurv <- cache_or_compute(
  "fig_median_survival_desurv_tcgacptac",
  splot_cutpoint(data_val_filtered, tar_fit_desurv, lp_stats, ntop = ntop_for_lp)
)

# Standard NMF at DeSurv k KM figure
fig_median_survival_std_desurvk <- cache_or_compute(
  "fig_median_survival_std_desurvk_tcgacptac",
  splot_cutpoint(data_val_filtered, fit_std_desurvk, lp_stats_std, ntop = NULL)
)

# ── Complete figure PDFs ──────────────────────────────────────────────────
fig_dir <- "figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Fig 3: PDAC factor structure ──────────────────────────────────────────
plot_3a <- fig_gene_overlap_heatmap_desurv$plot +
  theme(plot.margin = margin(t = 14, r = 2, b = 2, l = 2))
plot_3b <- fig_gene_overlap_heatmap_std_desurvk$plot +
  theme(plot.margin = margin(t = 14, r = 2, b = 2, l = 2))
legend_ab <- gtable::gtable_add_padding(
  fig_gene_overlap_heatmap_desurv$legend, padding = unit(c(0, 10, 0, 0), "pt"))
top_row_3 <- plot_grid(
  plot_3a, plot_3b, cowplot::ggdraw(legend_ab),
  ncol = 3, labels = c("A", "B", ""), align = "hv",
  label_size = 12, rel_widths = c(3.5, 3.5, 0.3)
)

plot_3c <- fig_variation_explained +
  theme(legend.position    = c(1, 0.5), legend.justification = c(1, 0),
        legend.background  = element_rect(color = "black"),
        axis.title         = element_text(size = 8),
        plot.margin        = margin(2, 30, 2, 2))

legend_d_plot <- ggplot(
  data.frame(x = 0, y = seq(-0.5, 1, length.out = 100)), aes(x = x, y = y, fill = y)
) +
  geom_tile() +
  scale_fill_gradientn(
    colors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
    limits = c(-0.5, 1), breaks = c(-0.4, 0, 0.4, 0.8), name = "Spearman\ncorrelation"
  ) +
  guides(fill = guide_colorbar(barwidth = unit(0.3, "cm"), barheight = unit(2, "cm"),
                               title.position = "top", title.hjust = 0.5)) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 6), legend.text = element_text(size = 6))
legend_d_grob <- cowplot::get_legend(legend_d_plot)

plot_3d <- plot_grid(
  fig_desurv_std_correlation_top50$plot + theme(plot.margin = margin(2, 2, 2, 20)),
  plot_grid(NULL, cowplot::ggdraw(legend_d_grob), nrow = 2, rel_heights = c(0.08, 0.92)),
  ncol = 2, rel_widths = c(4, 1)
)
bottom_row_3 <- plot_grid(plot_3c, plot_3d, ncol = 2, labels = c("C", "D"),
                          label_size = 12, rel_widths = c(0.55, 0.45))
ggsave(
  file.path(fig_dir, "fig3_tcgacptac.pdf"),
  plot_grid(top_row_3, bottom_row_3, nrow = 2, rel_heights = c(1.3, 0.7)),
  width = 7, height = 7
)
message("Saved fig3_tcgacptac.pdf")

# ── Fig 4: External validation ─────────────────────────────────────────────
d <- fig_hr_forest$data
d$dataset <- dplyr::recode(d$dataset,
  "Puleo_array" = "Puleo", "Moffitt_GEO_array" = "Moffitt",
  "PACA_AU_seq" = "PACA seq", "PACA_AU_array" = "PACA array")
pooled <- d %>%
  dplyr::mutate(logHR = log(HR), se = (log(upper) - log(lower)) / (2 * 1.96), w = 1 / se^2) %>%
  dplyr::group_by(method, factor_name) %>%
  dplyr::summarise(pooled_logHR = sum(logHR * w) / sum(w),
                   pooled_se    = sqrt(1 / sum(w)), .groups = "drop") %>%
  dplyr::mutate(HR = exp(pooled_logHR),
                lower = exp(pooled_logHR - 1.96 * pooled_se),
                upper = exp(pooled_logHR + 1.96 * pooled_se),
                dataset = "Pooled")
keep_cols <- c("factor_name", "HR", "lower", "upper", "dataset", "method", "row_type")
d$row_type <- "cohort"; pooled$row_type <- "pooled"
all_data <- rbind(d[, keep_cols], pooled[, keep_cols])
all_data$method <- factor(all_data$method, levels = c("DeSurv", "NMF"),
                          labels = c("DeSurv (D1-D3)", "Standard NMF (N1-N3)"))
all_data$dataset <- factor(all_data$dataset,
  levels = c("Puleo", "PACA seq", "PACA array", "Moffitt", "Dijk", "Pooled"))
cohort_cols   <- c("Dijk" = "#E69F00", "Moffitt" = "#56B4E9", "PACA array" = "#009E73",
                   "PACA seq" = "#0072B2", "Puleo" = "#CC79A7", "Pooled" = "#000000")
cohort_shapes <- c("Dijk" = 16, "Moffitt" = 17, "PACA array" = 15,
                   "PACA seq" = 25, "Puleo" = 8, "Pooled" = 18)
cohort_sizes  <- c("Dijk" = 2.2, "Moffitt" = 2.2, "PACA array" = 2.2,
                   "PACA seq" = 2.2, "Puleo" = 2.2, "Pooled" = 4)
base_size <- 8; dodge_w <- 0.6

desurv_dat <- all_data[all_data$method == "DeSurv (D1-D3)", ]
desurv_dat$factor_name <- factor(desurv_dat$factor_name, levels = c("D1", "D2", "D3"))
nmf_dat    <- all_data[all_data$method == "Standard NMF (N1-N3)", ]
nmf_dat$factor_name <- factor(nmf_dat$factor_name, levels = c("N1", "N2", "N3"))

make_forest_panel <- function(dat, title, highlight_level) {
  hl <- data.frame(factor_name = factor(highlight_level, levels = levels(dat$factor_name)), x = 1)
  ggplot(dat, aes(x = HR, y = factor_name)) +
    geom_tile(data = hl, aes(x = x, y = factor_name),
              width = 100, height = 0.8, fill = "steelblue", alpha = 0.08, inherit.aes = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lower, xmax = upper, colour = dataset, linewidth = row_type),
                  orientation = "y", width = 0, position = position_dodge(width = dodge_w)) +
    geom_point(aes(colour = dataset, shape = dataset, size = dataset),
               position = position_dodge(width = dodge_w)) +
    scale_linewidth_manual(values = c("cohort" = 0.6, "pooled" = 1.2), guide = "none") +
    scale_colour_manual(values = cohort_cols, name = NULL, drop = FALSE) +
    scale_shape_manual(values = cohort_shapes, name = NULL, drop = FALSE) +
    scale_size_manual(values = cohort_sizes, name = NULL, drop = FALSE) +
    scale_x_log10(limits = c(0.35, 3)) +
    labs(x = "Hazard ratio (95% CI)", y = NULL, title = title) +
    theme_classic(base_size = base_size) +
    theme(legend.position = "none",
          plot.title    = element_text(face = "bold", size = 9, hjust = 0.5),
          axis.text.y   = element_text(size = 8, face = "bold"),
          axis.text.x   = element_text(size = 8, face = "bold"),
          axis.title.x  = element_text(size = 9),
          plot.margin   = margin(2, 6, 2, 4))
}

p_leg_4 <- ggplot(all_data, aes(x = HR, y = factor_name,
                                colour = dataset, shape = dataset, size = dataset)) +
  geom_point() +
  scale_colour_manual(values = cohort_cols, name = NULL) +
  scale_shape_manual(values = cohort_shapes, name = NULL) +
  scale_size_manual(values = cohort_sizes, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.5)),
         shape  = guide_legend(nrow = 1), size = guide_legend(nrow = 1)) +
  theme_void(base_size = base_size) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))
forest_legend <- gtable::gtable_filter(ggplotGrob(p_leg_4), "guide-box")

plot_forest_4 <- plot_grid(
  plot_grid(make_forest_panel(desurv_dat, "DeSurv (D1-D3)", "D1"),
            make_forest_panel(nmf_dat, "Standard NMF (N1-N3)", "N1"),
            ncol = 2, align = "hv", axis = "tb"),
  forest_legend, nrow = 2, rel_heights = c(15, 1)
)

theme_pnas <- theme_classic(base_size = base_size) +
  theme(plot.title = element_text(face = "bold"),
        plot.margin = margin(6, 10, 6, 15),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.spacing.x  = unit(6, "pt"), legend.key.width = unit(12, "pt"))
km_text_size <- base_size

stack_surv <- function(surv_obj, title) {
  p <- surv_obj$plot + theme_pnas +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.title.y = element_text(size = km_text_size),
          axis.text.x  = element_text(size = km_text_size),
          axis.text.y  = element_text(size = km_text_size),
          plot.margin  = margin(4, 10, 0, 10), plot.title = element_text(size = 9)) +
    labs(title = title)
  t <- surv_obj$table + theme_pnas +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.x  = element_text(size = km_text_size),
          axis.text.y  = element_text(size = km_text_size),
          text         = element_text(size = km_text_size),
          axis.title.x = element_text(size = km_text_size),
          plot.title   = element_text(size = km_text_size),
          plot.margin  = margin(4, 10, 2, 10)) +
    labs(x = "Time (months)")
  t$layers[[1]]$aes_params$size <- km_text_size / ggplot2::.pt
  plot_grid(p, t, ncol = 1, rel_heights = c(2.5, 1.5), align = "v", axis = "lr")
}

km_legend_plot <- ggplot(
  data.frame(x = 1:2, y = 1:2,
             group = factor(c("Low", "High"), levels = c("Low", "High"))),
  aes(x = x, y = y, colour = group)
) +
  geom_line() +
  scale_colour_manual(values = c("Low" = "violetred2", "High" = "turquoise4"),
                      name = "Risk group") +
  theme_pnas +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = km_text_size),
        legend.title = element_text(size = km_text_size)) +
  guides(colour = guide_legend(nrow = 1))
km_legend_gg   <- ggplotGrob(km_legend_plot)
km_legend_grob <- km_legend_gg$grobs[
  sapply(km_legend_gg$grobs, function(x) x$name) == "guide-box"][[1]]

km_block_4 <- plot_grid(
  stack_surv(fig_median_survival_desurv, "DeSurv"),
  stack_surv(fig_median_survival_std_desurvk, "NMF"),
  ggdraw(km_legend_grob),
  nrow = 3, labels = c("B", "C", ""), label_size = 12, rel_heights = c(5, 5, 0.6)
)
ggsave(
  file.path(fig_dir, "fig4_tcgacptac.pdf"),
  plot_grid(plot_forest_4, km_block_4, ncol = 2, rel_widths = c(1.4, 1),
            labels = c("A", ""), label_size = 12),
  width = 7, height = 4.8
)
message("Saved fig4_tcgacptac.pdf")

message("=== Step 8 complete ===")
