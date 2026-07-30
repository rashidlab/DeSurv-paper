#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Standalone builder for Fig 2 panels A-D (biological factor structure),
# reproducing code/09a_figures.R's fig3_tcgacptac.pdf from cached objects in
# the flat results/ store WITHOUT the DeSurv package.
#
#   A/B: gene-overlap heatmaps (factor top genes vs established PDAC programs)
#        for DeSurv (A) and matched-rank standard NMF (B). Reference programs
#        come from data/derv/cmbSubtypes_formatted.RData (fallback to the
#        sibling DeSurv-paper checkout, exactly as 09a does).
#   C:   fig_variation_explained (cached ggplot; Shapley reconstruction share
#        vs survival contribution).
#   D:   Spearman correspondence of NMF vs DeSurv W-loadings (pure cor()).
#
# get_top_genes (R/get_top_genes.R) and make_gene_overlap_heatmap
# (R/figure_plot_helpers.R) are DeSurv-free; make_spearman_heatmap is inlined
# verbatim from 09a. Panel E (D1 score vs GATA6 RNA-ISH in COMPASS) reads the
# raw per-sample points cached in treated_cohort_stats.rds$gata6$points by
# code/11 -- aggregate-level derived values, so no restricted-data dependency here.
#
# Output: figures/standalone_preview/fig2_full.pdf  (7x9, panels A-E)
# ---------------------------------------------------------------------------

suppressMessages({
  library(ggplot2); library(cowplot); library(pheatmap); library(gtable)
  library(grid); library(RColorBrewer); library(dplyr)
})
source("R/get_top_genes.R")
source("R/figure_plot_helpers.R")   # make_gene_overlap_heatmap (DeSurv-free)

out_dir <- "figures/standalone_preview"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- make_spearman_heatmap: inlined verbatim from code/09a_figures.R:253 -----
make_spearman_heatmap <- function(c_mat, row_labels, col_labels) {
  rownames(c_mat) <- row_labels; colnames(c_mat) <- col_labels
  ph_args <- list(mat = c_mat, cluster_rows = FALSE, cluster_cols = FALSE,
    show_colnames = TRUE, show_rownames = TRUE, fontsize = 8, fontsize_number = 8,
    number_color = "black", breaks = seq(-0.5, 1, length.out = 101),
    display_numbers = TRUE, number_format = "%.2f", silent = TRUE)
  ph <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = FALSE)))
  ph_grob <- ph$gtable
  mat_idx <- which(ph_grob$layout$name == "matrix")
  mat_l <- ph_grob$layout$l[mat_idx]; mat_r <- ph_grob$layout$r[mat_idx]
  mat_t <- ph_grob$layout$t[mat_idx]; mat_b <- ph_grob$layout$b[mat_idx]
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
  list(plot = pheat, legend = legend_grob)
}

# --- inputs -----------------------------------------------------------------
tar_fit_desurv          <- readRDS("results/tar_fit_desurv_tcgacptac.rds")
fit_std_desurvk         <- readRDS("results/fit_std_desurvk_tcgacptac.rds")
tar_params_best         <- readRDS("results/tar_params_best_tcgacptac.rds")
fig_variation_explained <- readRDS("results/fig_variation_explained_tcgacptac.rds")

top_genes_path <- "data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(top_genes_path)) top_genes_path <- "../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(top_genes_path)) stop("Reference program file cmbSubtypes_formatted.RData not found")
load(top_genes_path)   # -> top_genes (reference programs), colors, subtypeList

ntop_value <- tar_params_best$ntop
heatmap_factor_labels     <- c("D1 Classical/restCAF", "D2 proCAF", "D3 Basal-like")
heatmap_factor_labels_std <- c("N1 Tumor", "N2 Exocrine", "N3 Microenviron.")

tar_tops_desurv      <- get_top_genes(W = tar_fit_desurv$W,  ntop = ntop_value)
tar_tops_std_desurvk <- get_top_genes(W = fit_std_desurvk$W, ntop = ntop_value)

# --- Panels A/B: gene-overlap heatmaps --------------------------------------
# Fixed, biologically prespecified reference panel displayed in BOTH panels
# (identical rows and order), chosen independently of the observed rank-biserial
# values. Single source of truth in R/fig2_display_panel.R (shared with the
# production pipeline code/09a_figures.R and the SI source table).
source("R/fig2_display_panel.R")

fig_gene_overlap_heatmap_desurv <- make_gene_overlap_heatmap(
  tar_fit_desurv, tar_tops_desurv$top_genes, top_genes,
  factor_labels = heatmap_factor_labels, title = "DeSurv", fontsize_row = 7,
  display_sigs = fig2_display_sigs)
fig_gene_overlap_heatmap_std_desurvk <- make_gene_overlap_heatmap(
  fit_std_desurvk, tar_tops_std_desurvk$top_genes, top_genes,
  factor_labels = heatmap_factor_labels_std, title = "NMF", fontsize_row = 7,
  display_sigs = fig2_display_sigs)

# --- Panel D: NMF vs DeSurv W-loading correspondence ------------------------
fig_desurv_std_correlation <- make_spearman_heatmap(
  cor(fit_std_desurvk$W, tar_fit_desurv$W, method = "spearman"),
  heatmap_factor_labels_std, heatmap_factor_labels)

# --- Panel E: D1 score vs GATA6 RNA-ISH in COMPASS --------------------------
# Raw per-sample points are cached in treated_cohort_stats.rds$gata6$points by
# code/11 (aggregate-level derived values; no raw expression), so panel E builds
# here with no restricted-data dependency.
g6 <- readRDS("results/treated_cohort_stats.rds")$gata6
plot_3e <- ggplot(g6$points, aes(x = factor(gata6), y = D1)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey92", linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, size = 1.3, alpha = 0.75, colour = "#08519c") +
  annotate("text", x = 0.6, y = max(g6$points$D1), hjust = 0, vjust = 1, size = 2.9,
           label = sprintf("Spearman~italic(r)==%.2f", g6$rho), parse = TRUE) +
  annotate("text", x = 0.6, y = max(g6$points$D1) - 0.45, hjust = 0, vjust = 1, size = 2.9,
           label = sprintf("italic(P)<0.001*','~n==%d", g6$n), parse = TRUE) +
  labs(x = "GATA6 RNA-ISH level", y = "DeSurv D1 score (z)") +
  theme_classic(base_size = 9) +
  theme(axis.title = element_text(size = 8))

# --- assemble (verbatim code/09a_figures.R:336-379) -------------------------
plot_3a <- fig_gene_overlap_heatmap_desurv$plot +
  theme(plot.margin = margin(t = 14, r = 2, b = 2, l = 2))
plot_3b <- fig_gene_overlap_heatmap_std_desurvk$plot +
  theme(plot.margin = margin(t = 14, r = 2, b = 2, l = 2))
legend_ab <- gtable::gtable_add_padding(
  fig_gene_overlap_heatmap_desurv$legend, padding = unit(c(0, 6, 0, 6), "pt"))
# Wider legend column so the "Rank-biserial enrichment" title is not clipped on
# the right of panel B.
top_row_3 <- plot_grid(plot_3a, plot_3b, cowplot::ggdraw(legend_ab),
  ncol = 3, labels = c("A", "B", ""), align = "hv", label_size = 12,
  rel_widths = c(3.4, 3.4, 0.85))

plot_3c <- fig_variation_explained +
  theme(legend.position = c(1, 0.5), legend.justification = c(1, 0),
        legend.background = element_rect(color = "black"),
        axis.title = element_text(size = 8), plot.margin = margin(2, 30, 2, 2))

legend_d_plot <- ggplot(data.frame(x = 0, y = seq(-0.5, 1, length.out = 100)),
                        aes(x = x, y = y, fill = y)) + geom_tile() +
  scale_fill_gradientn(
    colors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100),
    limits = c(-0.5, 1), breaks = c(-0.4, 0, 0.4, 0.8), name = "Spearman\ncorrelation") +
  guides(fill = guide_colorbar(barwidth = unit(0.3, "cm"), barheight = unit(2, "cm"),
                               title.position = "top", title.hjust = 0.5)) +
  theme_void() + theme(legend.position = "right",
        legend.title = element_text(size = 6), legend.text = element_text(size = 6))
legend_d_grob <- cowplot::get_legend(legend_d_plot)

plot_3d <- plot_grid(
  fig_desurv_std_correlation$plot + theme(plot.margin = margin(2, 2, 2, 20)),
  plot_grid(NULL, cowplot::ggdraw(legend_d_grob), nrow = 2, rel_heights = c(0.08, 0.92)),
  ncol = 2, rel_widths = c(4, 1))
# Bottom row: C, D, E side by side -> keeps a ~square aspect so the figure renders
# full-width (out.width='\textwidth') without distortion or overflowing the page.
bottom_row_3 <- plot_grid(plot_3c, plot_3d, plot_3e, ncol = 3, labels = c("C", "D", "E"),
                          label_size = 12, rel_widths = c(0.37, 0.33, 0.30))

fig2 <- plot_grid(top_row_3, bottom_row_3, nrow = 2, rel_heights = c(1.3, 0.72))
ggsave(file.path(out_dir, "fig2_full.pdf"), fig2, width = 7, height = 7)
cat("Wrote", file.path(out_dir, "fig2_full.pdf"), "- Fig 2 panels A-E (7x7, full-width).\n")
