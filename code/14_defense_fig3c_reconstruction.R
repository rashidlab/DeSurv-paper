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
  library(cowplot)
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

# ── Broken y-axis: D1's survival contribution (Delta-l ~ 66) dwarfs every other
#    factor (all < 1.2), collapsing them onto the axis. We split the panel so the
#    near-zero factors get full resolution (N2 ~ 1.2 and D3 ~ 1.0 separate clearly
#    from N1/N3/D2 ~ 0) while D1 still reads as dramatically higher. ─────────────
pal <- c(NMF = "#c0392b", DeSurv = "#2166ac")
df_plot$method <- factor(df_plot$method, levels = c("DeSurv", "NMF"))

# Split the data by panel so geom_text_repel only labels each panel's own points
# (clipping alone leaves stray labels for off-screen points).
df_hi <- df_plot[df_plot$delta_loglik > 10, ]    # D1 only
df_lo <- df_plot[df_plot$delta_loglik <= 10, ]    # every other factor

x_scale <- scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                              expand = expansion(mult = c(0.08, 0.12)))

base_layers <- list(
  geom_point(size = 6.5),
  geom_text_repel(size = 7, fontface = "bold", max.overlaps = Inf,
                  box.padding = 0.7, point.padding = 0.5,
                  min.segment.length = Inf, force = 2, show.legend = FALSE),
  scale_color_manual(values = pal, name = NULL, drop = FALSE),
  x_scale
)

# Top panel: only the D1 region (zoomed), no x-axis, legend lives here.
p_top <- ggplot(df_hi, aes(variance_explained, delta_loglik,
                           label = factor_label, color = method)) +
  base_layers +
  coord_cartesian(ylim = c(62, 70)) +
  scale_y_continuous(breaks = c(65)) +
  theme_classic(base_size = 21) +
  theme(
    axis.title   = element_blank(),
    axis.text.y  = element_text(color = "black"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    legend.position = "none",
    plot.margin  = margin(8, 16, 0, 70)
  )

# Bottom panel: every other factor, full resolution; carries the x-axis.
p_bot <- ggplot(df_lo, aes(variance_explained, delta_loglik,
                           label = factor_label, color = method)) +
  base_layers +
  coord_cartesian(ylim = c(-0.12, 1.5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0, 1.5)) +
  labs(x = "Contribution to reconstruction\n(Shapley share)") +
  theme_classic(base_size = 21) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    axis.text    = element_text(color = "black"),
    legend.position = "none",
    plot.margin  = margin(0, 16, 10, 70)
  )

stacked <- plot_grid(p_top, p_bot, ncol = 1, rel_heights = c(1, 1.9),
                     align = "v", axis = "lr")

# Shared rotated y-axis title + broken-axis slash marks straddling the y-axis line.
y_title <- expression(atop(Delta ~ "partial log-likelihood",
                           "(full vs. " * italic(k) * "-1 factor model)"))
boundary <- 1.9 / 2.9                     # y (NPC) of the top/bottom panel join
brk_x <- c(0.160, 0.196)                  # straddles the y-axis line (NPC x ~ 0.177)
body <- ggdraw(stacked) +
  draw_label(y_title, x = 0.045, y = 0.55, angle = 90,
             fontface = "bold", size = 18, hjust = 0.5) +
  draw_line(x = brk_x, y = boundary + c(-0.004, 0.016), linewidth = 1.1) +
  draw_line(x = brk_x, y = boundary + c(0.010, 0.030), linewidth = 1.1)

# Shared horizontal legend, placed below the plot.
legend_row <- get_legend(
  p_bot + theme(legend.position = "bottom", legend.direction = "horizontal",
                legend.text = element_text(size = 18, face = "bold"))
)
fig <- plot_grid(body, legend_row, ncol = 1, rel_heights = c(1, 0.08))

ggsave(file.path(OUT, "new_3c_reconstruction.pdf"), fig,
       width = 8, height = 6.2, device = cairo_pdf)
ggsave(file.path(OUT, "new_3c_reconstruction.png"), fig,
       width = 8, height = 6.2, dpi = 300, bg = "white", type = "cairo")
message("Saved new_3c_reconstruction (8x6.2 in, broken y-axis)")
