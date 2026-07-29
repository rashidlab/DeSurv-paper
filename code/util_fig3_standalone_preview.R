#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Standalone preview builder for Fig 3 panels B (DeSurv KM) and C (supervised
# heatmap) — the two Fig-3 panels that are reconstructable in the clean
# checkout WITHOUT the DeSurv package or the (absent) precomputed model store.
#
# WHY THIS EXISTS: results/precomputed/ is empty in this checkout, so panel A
# (the frozen-LP HR forest, built by compute_hrs from data_val_filtered +
# tar_fit_desurv) cannot be rebuilt here — that data lives on Longleaf.
# Panels B and C, however, come from loose cached objects in results/:
#   - fig_median_survival_desurv_tcgacptac.rds  (ggsurvplot; panel B)
#   - desurv_vs_supervised_tuned.rds$axis_decomposition  (panel C values)
# Neither needs DeSurv, DiceKriging, or ggplot2<4.0, so both build with the
# default-library ggplot2 4.0 / survminer / cowplot.
#
# Panel C here is the *new* Fig-3C element (replaces the old NMF KM). The
# authoritative full-figure regeneration (A+B+C) still happens via
# code/09a_figures.R on Longleaf; see figures/FIGURE_REGENERATION.md.
#
# Output: figures/standalone_preview/{fig3C_supervised_heatmap,fig3B_desurv_km,
#         fig3_BC_preview}.pdf
# ---------------------------------------------------------------------------

suppressMessages({
  library(ggplot2)
  library(cowplot)
  library(survminer)
})

out_dir <- "figures/standalone_preview"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Panel C: tuned-supervised vs DeSurv-program correspondence heatmap ------
# |Spearman r| of each tuned supervised risk score with the three DeSurv
# programs. Message: supervised scores concentrate on D1 and are near-orthogonal
# to the separately-resolved D2 (proCAF) stromal program. No DeSurv-LP row
# (its trained score has beta(D2)=0, so it would also read ~0 on D2).
axd <- readRDS("results/desurv_vs_supervised_tuned.rds")$axis_decomposition

hm_df <- data.frame(
  method  = factor(rep(c("Supervised PCA", "Penalized Cox"), each = 3),
                   levels = c("Supervised PCA", "Penalized Cox")),
  program = factor(rep(c("D1", "D2", "D3"), 2), levels = c("D1", "D2", "D3")),
  r = c(abs(as.numeric(axd["Supervised PCA", c("D1", "D2", "D3")])),
        abs(as.numeric(axd["Sparse Cox",     c("D1", "D2", "D3")])))
)

fig3C <- ggplot(hm_df, aes(program, method, fill = r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.1) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1),
                      name = expression("|" * italic(r) * "|")) +
  labs(x = NULL, y = NULL, title = "Supervised score vs DeSurv program") +
  theme_minimal(base_size = 9) +
  theme(plot.title   = element_text(size = 9, hjust = 0.5),
        panel.grid   = element_blank(),
        axis.text    = element_text(color = "black"),
        legend.position = "right")

ggsave(file.path(out_dir, "fig3C_supervised_heatmap.pdf"), fig3C,
       width = 3.4, height = 2.2)

# --- Panel B: pooled DeSurv KM (load cached ggsurvplot, no recompute) --------
km <- readRDS("results/fig_median_survival_desurv_tcgacptac.rds")
fig3B <- km$plot +
  ggtitle("DeSurv") +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

ggsave(file.path(out_dir, "fig3B_desurv_km.pdf"), fig3B, width = 3.4, height = 3.0)

# --- B | C composite: the new right column of Fig 3 -------------------------
bc <- plot_grid(fig3B, fig3C, ncol = 1, labels = c("B", "C"),
                rel_heights = c(5, 4.2))
ggsave(file.path(out_dir, "fig3_BC_preview.pdf"), bc, width = 3.6, height = 6.4)

cat("Wrote:\n",
    file.path(out_dir, "fig3C_supervised_heatmap.pdf"), "\n",
    file.path(out_dir, "fig3B_desurv_km.pdf"), "\n",
    file.path(out_dir, "fig3_BC_preview.pdf"), "\n",
    "NOTE: panel A (HR forest) needs the precomputed store (empty here) -> Longleaf.\n")
