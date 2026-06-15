#!/usr/bin/env Rscript
# code/17_defense_fig2_panels.R — Defense version of main-text Fig 2 panels A-C
#
# Slide "DeSurv recovers the true prognostic program and rank" shows simulation
# panels A (test C-index), B (precision), C (selected-k distribution) in a row.
# The deck previously used a crop-and-stitch of the paper composite, which left
# a white seam through panel B's y-axis title. This rebuilds the row directly
# from the cached panel objects (results/sim_figs_by_scenario.rds), so there is
# no seam, and bumps the fonts for slide legibility.
#
# Run:  Rscript code/17_defense_fig2_panels.R

message("=== Defense Fig 2 panels (simulation recovery) ===")
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

source("code/00_helpers.R")
OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Primary ("alternative") scenario tuned with BO over n_top — same selection as
# code/09c_sim_figures.R uses for the paper Fig 2.
sim_figs     <- load_precomputed("sim_figs_by_scenario")
scenario_ids <- sapply(sim_figs, function(x) x$scenario_id)
analysis_ids <- sapply(sim_figs, function(x) x$analysis_id)
alt <- sim_figs[[which(scenario_ids == "R0_easy" & analysis_ids == "bo_tune_ntop")]]

# Enlarged fonts for slide scale; legend on top to match the cached panels.
slide_theme <- theme(
  plot.title   = element_text(size = 19, face = "bold", hjust = 0.5),
  axis.title   = element_text(size = 18),
  axis.text    = element_text(size = 15, color = "black"),
  legend.title = element_text(size = 16),
  legend.text  = element_text(size = 15),
  strip.text   = element_text(size = 16, face = "bold"),
  legend.position = "top"
)

pA <- alt$cindex_box    + labs(title = NULL) + scale_y_continuous(limits = c(.5, 1)) + slide_theme
pB <- alt$precision_box + labs(title = NULL) + slide_theme
pC <- alt$k_hist        + slide_theme

row <- plot_grid(pA, pB, pC, ncol = 3, labels = c("A", "B", "C"),
                 label_size = 24, rel_widths = c(1, 1, 1.12))

ggsave(file.path(OUT, "fig2_ABC_row.pdf"), row, width = 13.2, height = 4.4,
       device = cairo_pdf)
ggsave(file.path(OUT, "fig2_ABC_row.png"), row, width = 13.2, height = 4.4,
       dpi = 300, bg = "white")
message("Saved fig2_ABC_row (13.2x4.4 in, rebuilt from source — no stitch seam)")
