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

# ── 2. Rasterize SI figure PDFs (figures/*.pdf) -> figures/defense/bk_*.png ──
pdf_map <- c(
  si_fig_converge_tcgacptac      = "bk_converge",      # B4
  si_fig_sim_null_mixed_tcgacptac = "bk_null_mixed",   # B7
  si_fig_nmf_diagnostics_tcgacptac = "bk_nmf_diag",    # B8
  cv_cindex_by_k_primary         = "bk_cindex_by_k",   # B9
  si_fig_nmf_k7_heatmap_tcgacptac = "bk_nmf_k7",       # B10
  si_fig_cutpoint_km_tcgacptac   = "bk_cutpoint_km",   # B11
  si_fig_subtype_overlap_tcgacptac = "bk_subtype_overlap" # B12
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
