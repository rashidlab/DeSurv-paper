#!/usr/bin/env Rscript
# code/14_defense_fig3c_reconstruction.R — Defense version of main-text Fig 3C
#
# Deck-scale "reconstruction vs prognosis" scatter, matching the MANUSCRIPT
# framing (paper/04_results_REVISED.Rmd, Fig 3C). The x-axis is each factor's
# per-factor contribution to RECONSTRUCTION (normalized reconstruction Shapley
# share, which partitions the model's reconstruction and sums to 100%); the
# y-axis is its survival contribution (Type III / leave-one-out partial
# log-likelihood). This supersedes the earlier centered-variance reframe
# (former code/13_defense_fig3c_variance.R): the deck now uses the same metric
# and numbers as the paper.
#
# We reuse the cached manuscript ggplot object so the values are identical to
# Fig 3C; only the fonts/sizes are scaled up for slides.
#
# Run:  Rscript code/14_defense_fig3c_reconstruction.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Cached manuscript Fig 3C object (built in code/09a_figures.R).
# $data columns: variance_explained (reconstruction Shapley share, 0-1),
#                delta_loglik, factor_label, method
fig_src <- readRDS(file.path("results", "fig_variation_explained_tcgacptac.rds"))
df_plot <- fig_src$data

cat("\nPer-factor values (reconstruction Shapley x-axis):\n")
print(within(df_plot, {
  variance_explained <- round(variance_explained, 3)
  delta_loglik       <- round(delta_loglik, 2)
})[, c("method", "factor_label", "variance_explained", "delta_loglik")])

# ── Plot (deck-scale fonts; same visual grammar as the centered-variance deck
#    figure it replaces, so the slide look is consistent) ──────────────────────
fig <- ggplot(df_plot, aes(variance_explained, delta_loglik,
                           label = factor_label, color = method)) +
  geom_point(size = 5) +
  geom_text_repel(size = 5.2, fontface = "bold", max.overlaps = Inf,
                  box.padding = 0.7, point.padding = 0.5,
                  segment.size = 0.3, force = 2, show.legend = FALSE) +
  scale_color_manual(values = c(NMF = "#c0392b", DeSurv = "#2166ac"),
                     name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    x = "Contribution to reconstruction\n(Shapley share)",
    y = expression(atop(Delta ~ "partial log-likelihood",
                        "(full vs. " * italic(k) * "-1 factor model)"))
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.position = c(0.5, 0.93),
    legend.direction = "horizontal",
    legend.text  = element_text(size = 15, face = "bold"),
    plot.margin  = margin(12, 16, 10, 12)
  )

ggsave(file.path(OUT, "new_3c_reconstruction.pdf"), fig,
       width = 8, height = 6.2, device = cairo_pdf)
ggsave(file.path(OUT, "new_3c_reconstruction.png"), fig,
       width = 8, height = 6.2, dpi = 300, bg = "white", type = "cairo")
message("Saved new_3c_reconstruction (8x6.2 in)")
