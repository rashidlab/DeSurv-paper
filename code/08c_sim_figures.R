#!/usr/bin/env Rscript
# code/08c_sim_figures.R — Simulation figure PDFs
#
# Generates complete multi-panel figure PDFs for all simulation-based figures:
#   figures/fig2_tcgacptac.pdf                  (main text Fig 2)
#   figures/si_fig_converge_tcgacptac.pdf       (SI S1)
#   figures/si_fig_sim_null_mixed_tcgacptac.pdf (SI S2)
#
# Inputs: sim_figs_by_scenario, desurv_seed_fits (step 07),
#         fig_bo_heat (step 08)
# Runtime: ~1 min

message("=== Step 8c: Simulation figures ===")
source("code/00_helpers.R")
library(ggplot2)
library(cowplot)
library(dplyr)

fig_dir <- "figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set_fig_font <- function(plot_obj, size = 10) {
  if (inherits(plot_obj, "ggplot")) {
    plot_obj + ggplot2::theme(text = ggplot2::element_text(size = size))
  } else {
    plot_obj
  }
}

# ── Load simulation results ──────────────────────────────────────────────
sim_figs_by_scenario <- load_precomputed("sim_figs_by_scenario")
scenario_ids <- sapply(sim_figs_by_scenario, function(x) x$scenario_id)
analysis_ids <- sapply(sim_figs_by_scenario, function(x) x$analysis_id)

alt_plots   <- sim_figs_by_scenario[[which(scenario_ids == "R0_easy"  & analysis_ids == "bo_tune_ntop")]]
null_plots  <- sim_figs_by_scenario[[which(scenario_ids == "R00_null" & analysis_ids == "bo_tune_ntop")]]
mixed_plots <- sim_figs_by_scenario[[which(scenario_ids == "R_mixed"  & analysis_ids == "bo_tune_ntop")]]

# ── Fig 2: Simulation panels A-C + BO heatmap D ──────────────────────────
fig_bo_heat <- load_precomputed("fig_bo_heat_tcgacptac")
fig_bo_heat <- fig_bo_heat +
  guides(fill = guide_colorbar(barheight = unit(2, "cm"), barwidth = unit(0.3, "cm"))) +
  theme(legend.title = element_text(size = 7), legend.text = element_text(size = 7))
leg_grob    <- ggplotGrob(fig_bo_heat)
heat_legend <- leg_grob$grobs[[which(leg_grob$layout$name == "guide-box-right")]]
fig_bo_heat <- fig_bo_heat + theme(legend.position = "none")

upper <- plot_grid(
  alt_plots$cindex_box    + labs(title = NULL) + scale_y_continuous(limits = c(.5, 1)),
  alt_plots$precision_box + labs(title = NULL),
  NULL,
  ncol = 3, labels = c("A", "B", ""), rel_widths = c(4, 4, .8)
)
lower <- plot_grid(
  alt_plots$k_hist, fig_bo_heat, heat_legend,
  ncol = 3, labels = c("C", "D", ""), rel_widths = c(4, 4, .8)
)
ggsave(file.path(fig_dir, "fig2_tcgacptac.pdf"),
       plot_grid(upper, lower, nrow = 2),
       width = 6.5, height = 6)
message("Saved fig2_tcgacptac.pdf")

# ── SI S1: Convergence trajectories ──────────────────────────────────────
desurv_seed_fits <- load_precomputed("desurv_seed_fits_tcgacptac")
set.seed(147)
seeds       <- sample(seq_along(desurv_seed_fits$fits), 5, replace = FALSE)
seed_labels <- setNames(paste0("Init ", seq_along(seeds)), as.character(seeds))

lossit <- dplyr::bind_rows(lapply(seeds, function(i) {
  loss_vec <- desurv_seed_fits$fits[[i]]$lossit
  data.frame(lossit  = loss_vec,
             slossit = desurv_seed_fits$fits[[i]]$slossit,
             nlossit = desurv_seed_fits$fits[[i]]$nlossit,
             iter    = seq_along(loss_vec), init = i)
})) |>
  dplyr::filter(iter < 5000) |>
  dplyr::group_by(init) |>
  dplyr::arrange(iter) |>
  dplyr::mutate(rel_delta_loss = lossit / dplyr::first(lossit),
                init_label = seed_labels[as.character(init)]) |>
  dplyr::filter(!is.na(rel_delta_loss))

ggsave(
  file.path(fig_dir, "si_fig_converge_tcgacptac.pdf"),
  ggplot(lossit, aes(y = rel_delta_loss, x = iter, color = init_label)) +
    geom_line(linewidth = 0.8) +
    theme_minimal(base_size = 10) +
    labs(x = "Iteration", y = "Loss normalized to initial value", color = "Initialization"),
  width = 5, height = 3
)
message("Saved si_fig_converge_tcgacptac.pdf")

# ── SI S2: Null and mixed simulation scenarios ────────────────────────────
null_f  <- lapply(null_plots,  set_fig_font, size = 10)
mixed_f <- lapply(mixed_plots, set_fig_font, size = 10)

null_cindex <- null_f$cindex_box + labs(title = NULL) + scale_y_continuous(limits = c(0.3, 1))
null_cindex$data$scenario_id <- "Null scenario (β = 0)"
null_khist  <- null_f$k_hist + labs(title = NULL)
null_khist$data$scenario_id  <- "Null scenario (β = 0)"
null_row <- plot_grid(null_cindex, null_khist, ncol = 2, labels = c("A", "B"), rel_widths = c(1, 1.2))

mixed_cindex <- mixed_f$cindex_box + labs(title = NULL) + scale_y_continuous(limits = c(0.3, 1)) +
  theme(plot.clip = "off", plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))
mixed_cindex$data$scenario_id <- "Mixed scenario"
mixed_prec <- if (!is.null(mixed_f$precision_breakdown)) {
  mixed_f$precision_breakdown + labs(title = NULL)
} else {
  mixed_f$precision_box + labs(title = NULL)
}
mixed_row1 <- plot_grid(mixed_cindex, mixed_prec, ncol = 2, labels = c("C", "D"), rel_widths = c(1, 2))

mixed_beta  <- if (!is.null(mixed_f$matched_beta_box)) mixed_f$matched_beta_box + labs(title = NULL) else NULL
mixed_khist <- mixed_f$k_hist + labs(title = NULL)
mixed_khist$data$scenario_id <- "Mixed scenario"

fig_null_mixed <- if (!is.null(mixed_beta)) {
  plot_grid(null_row, mixed_row1,
            plot_grid(mixed_beta, mixed_khist, ncol = 2, labels = c("E", "F"), rel_widths = c(1, 1.2)),
            nrow = 3, rel_heights = c(1, 1, 1))
} else {
  plot_grid(null_row, mixed_row1,
            plot_grid(mixed_khist, ncol = 1, labels = "E"),
            nrow = 3, rel_heights = c(1, 1, 0.7))
}
ggsave(file.path(fig_dir, "si_fig_sim_null_mixed_tcgacptac.pdf"), fig_null_mixed, width = 6.5, height = 7)
message("Saved si_fig_sim_null_mixed_tcgacptac.pdf")

message("=== Step 8c complete ===")
