#!/usr/bin/env Rscript
# code/18_defense_backup_panels.R — Backup-slide figure assets for the defense deck
#
# The Q&A backup slides (B2-B14) mostly reuse the paper's SI figures. Reveal.js
# needs raster PNGs (the self-contained HTML can't embed PDFs), so this script:
#   1. exports the BO C-index heatmap (Fig 2D) from its cached ggplot, and
#   2. rasterizes the relevant SI figure PDFs to PNG via poppler's pdftoppm.
#
# Requires poppler (pdftoppm) on PATH. Greek glyphs in a couple of PDF titles
# (alpha/lambda) drop out in conversion because the PDFs reference the Symbol
# font; the affected values are restated in the slide text, so this is cosmetic.
#
# Run:  Rscript code/18_defense_backup_panels.R

message("=== Defense backup-slide panels ===")
suppressPackageStartupMessages(library(ggplot2))
source("code/00_helpers.R")

OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1. BO heatmap (Fig 2D / B6), enlarged for slide scale ───────────────────
p_bo <- load_precomputed("fig_bo_heat_tcgacptac") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_text(size = 17), legend.text = element_text(size = 15),
        axis.title = element_text(size = 19), axis.text = element_text(size = 15))
ggsave(file.path(OUT, "bk_bo_heat.png"), p_bo, width = 8, height = 5.6,
       dpi = 300, bg = "white")
message("Saved bk_bo_heat.png")

# ── 1b. C-index vs rank, panel B only (n_top = ALL) for B9 ──────────────────
# Rebuilt from the CV-grid CSVs (not the SI PDF) so the lambda/xi symbols render
# and the fonts are slide-scale. The deck shows only this n_top = ALL panel.
suppressPackageStartupMessages({library(dplyr); library(cowplot)})
source("R/cv_grid_helpers.R")
cvdir <- file.path("results", "cv_grid")
if (file.exists(file.path(cvdir, "cv_grid_summary.csv"))) {
  cfgB <- data.frame(ntop = NA_integer_, lambda = 0.349, nu = 0.056,
                     label = "ntopALL", stringsAsFactors = FALSE)
  pl <- plot_cindex_by_k(
    cv_grid_summary     = read.csv(file.path(cvdir, "cv_grid_summary.csv")),
    cv_grid_best_alpha  = read.csv(file.path(cvdir, "cv_grid_best_alpha.csv")),
    cv_grid_val_summary = read.csv(file.path(cvdir, "cv_grid_val_summary.csv")),
    configs             = cfgB)
  pB <- pl[[1]] + theme_classic(base_size = 20) +
    theme(legend.position = "bottom",
          legend.text  = element_text(size = 19),
          plot.subtitle = element_text(size = 19),
          strip.text   = element_text(size = 19, face = "bold"),
          axis.title   = element_text(size = 20),
          axis.text    = element_text(size = 15, color = "black")) +
    geom_line(linewidth = 1) + geom_point(size = 3)
  ggsave(file.path(OUT, "bk_cindex_by_k.png"), pB, width = 11, height = 5,
         dpi = 300, bg = "white")
  message("Saved bk_cindex_by_k.png (panel B / n_top=ALL, enlarged)")
} else {
  warning("cv_grid CSVs not found — skipping bk_cindex_by_k")
}

# ── 2. Rasterize SI figure PDFs (figures/*.pdf) -> figures/defense/bk_*.png ──
# Only the convergence figure is still sourced from its SI PDF; the other dense
# backups (B7/B8/B10/B11/B12) are rebuilt from source with large fonts in
# code/19_defense_backup_figs.R.
pdf_map <- c(
  si_fig_converge_tcgacptac = "bk_converge"  # B4
)
if (nzchar(Sys.which("pdftoppm"))) {
  for (src in names(pdf_map)) {
    pdf <- file.path("figures", paste0(src, ".pdf"))
    if (!file.exists(pdf)) { warning("missing ", pdf); next }
    system2("pdftoppm", c("-png", "-r", "200", "-singlefile",
                          shQuote(pdf), shQuote(file.path(OUT, pdf_map[[src]]))))
    message("Converted ", pdf_map[[src]])
  }
} else {
  warning("pdftoppm not found on PATH — skipping SI PDF conversions")
}
message("=== Defense backup-slide panels complete ===")
