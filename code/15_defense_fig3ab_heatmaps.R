#!/usr/bin/env Rscript
# code/15_defense_fig3ab_heatmaps.R — Defense version of main-text Fig 3A/3B
#
# Slide 13 ("Supervision reorganizes PDAC factors") shows the DeSurv (A) and
# standard-NMF (B) factor-vs-signature Spearman correlation heatmaps. The deck
# previously cropped these out of the paper composite, leaving the row/column
# labels too small to read at slide scale. This rebuilds panels A+B directly
# from the cached fits with enlarged fonts and stitches them side-by-side with
# a shared colour bar.
#
# Run:  Rscript code/15_defense_fig3ab_heatmaps.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

source("code/00_helpers.R")
source("R/get_top_genes.R")
source("R/figure_plot_helpers.R")

OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read1 <- function(name) readRDS(file.path("results", paste0(name, ".rds")))

tar_fit_desurv  <- read1("tar_fit_desurv_tcgacptac")
fit_std_desurvk <- read1("fit_std_desurvk_tcgacptac")
ntop_value      <- read1("tar_params_best_tcgacptac")$ntop

# Reference gene lists from combined subtype data (same fallback as 09a).
top_genes_path <- "data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(top_genes_path))
  top_genes_path <- "../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData"
load(top_genes_path)  # loads: top_genes, ...

heatmap_factor_labels     <- c("D1 Classical/restCAF", "D2 proCAF", "D3 Basal-like")
heatmap_factor_labels_std <- c("N1 Tumor", "N2 Exocrine", "N3 Microenviron.")

tops_desurv      <- get_top_genes(W = tar_fit_desurv$W,  ntop = ntop_value)
tops_std_desurvk <- get_top_genes(W = fit_std_desurvk$W, ntop = ntop_value)

# Enlarged fonts for slide legibility (fontsize_row = row/col labels;
# fontsize = title; legend_fontsize = colour-bar labels).
FS_ROW <- 13; FS_TITLE <- 17; FS_LEG <- 13

hm_desurv <- make_gene_overlap_heatmap(
  tar_fit_desurv, tops_desurv$top_genes, top_genes,
  factor_labels = heatmap_factor_labels, title = "DeSurv",
  fontsize_row = FS_ROW, fontsize = FS_TITLE, legend_fontsize = FS_LEG)

hm_nmf <- make_gene_overlap_heatmap(
  fit_std_desurvk, tops_std_desurvk$top_genes, top_genes,
  factor_labels = heatmap_factor_labels_std, title = "NMF",
  fontsize_row = FS_ROW, fontsize = FS_TITLE, legend_fontsize = FS_LEG)

plot_a <- hm_desurv$plot + theme(plot.margin = margin(t = 14, r = 2, b = 2, l = 2))
plot_b <- hm_nmf$plot    + theme(plot.margin = margin(t = 14, r = 2, b = 2, l = 2))
legend_ab <- gtable::gtable_add_padding(
  hm_desurv$legend, padding = unit(c(0, 10, 0, 0), "pt"))

fig <- plot_grid(
  plot_a, plot_b, ggdraw(legend_ab),
  ncol = 3, labels = c("A", "B", ""), align = "hv",
  label_size = 22, rel_widths = c(3.5, 3.5, 0.45)
)

ggsave(file.path(OUT, "fig3_AB_heatmaps.pdf"), fig,
       width = 13, height = 6.2, device = cairo_pdf)
ggsave(file.path(OUT, "fig3_AB_heatmaps.png"), fig,
       width = 13, height = 6.2, dpi = 300, bg = "white")
message("Saved fig3_AB_heatmaps (13x6.2 in, enlarged fonts)")
