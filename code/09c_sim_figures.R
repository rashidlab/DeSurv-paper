#!/usr/bin/env Rscript
# code/09c_sim_figures.R — Simulation figure PDFs
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
source("R/theme_nature.R")   # single source of truth: theme_nature(), desurv_method_cols

# Recolor a prebuilt sim panel (step 07) to the semantic DeSurv/NMF pair and put
# it on the shared theme. Legends are suppressed here; a single shared Method
# legend is composed once per figure (avoids the duplicate fill+colour legends).
recolor <- function(p) p +
  ggplot2::scale_fill_manual(values = desurv_method_cols, guide = "none") +
  ggplot2::scale_colour_manual(values = desurv_method_cols, guide = "none") +
  theme_nature()


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

# ── SI: BO tuning surface (relocated out of the main-text simulations figure) ───────
# GP-predicted CV C-index over factorization rank (k) x supervision strength;
# a model-selection sensitivity display, now an SI figure rather than main text.
fig_bo_heat <- load_precomputed("fig_bo_heat_tcgacptac") +
  guides(fill = guide_colorbar(barheight = unit(2.4, "cm"), barwidth = unit(0.35, "cm"))) +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
ggsave(file.path(FIGURE_DIR, "si_fig_bo_tuning_surface_tcgacptac.pdf"),
       fig_bo_heat, width = 5, height = 4)
message("Saved si_fig_bo_tuning_surface_tcgacptac.pdf")

# ── Fig 2 (main text): simulation panels A-D ──────────────────────────────
# 4B: recovery of true prognostic genes (relabeled away from "Precision")
panel_b <- recolor(alt_plots$precision_box) +
  labs(title = NULL, y = "Proportion of selected genes\nfrom true prognostic program")

# 4D: scenario gradient in gene recovery. R0_easy = prognostic program explains LOW
# variance (separated from dominant variance); R_mixed = partial overlap between
# variance and prognosis (per code/07_simulations.R). DeSurv recovers the true
# program across both regimes while standard NMF, which chases variance, does not.
# Precision is undefined in the null scenario (no true prognostic program), so it
# is excluded here (null is covered by 4A C-index = 0.5 and the SI).
prec_grad <- do.call(rbind, lapply(list(alt_plots, mixed_plots), function(b)
  b$cindex_box$data[, c("scenario_id", "method", "precision")]))
prec_grad <- prec_grad[!is.na(prec_grad$precision), ]
prec_grad$scenario <- factor(prec_grad$scenario_id, levels = c("R0_easy", "R_mixed"),
  labels = c("Prognostic\n(low variance)", "Partial\noverlap"))
prec_grad$method <- factor(prec_grad$method, levels = c("DeSurv", "NMF"))
panel_d <- ggplot(prec_grad, aes(x = scenario, y = precision, fill = method)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, position = position_dodge(0.8)) +
  scale_fill_manual(values = desurv_method_cols, name = "Method") +   # DeSurv #0072B2, NMF #D55E00
  labs(x = NULL, y = "Proportion from true\nprognostic program", title = NULL) +
  theme_nature() +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8))

# One shared Method legend (from panel_d) placed in the upper row's third column.
method_legend <- cowplot::get_legend(
  panel_d + theme_nature() + theme(legend.position = "right"))
upper <- plot_grid(
  recolor(alt_plots$cindex_box) + labs(title = NULL) + scale_y_continuous(limits = c(.5, 1)),
  panel_b, method_legend,
  ncol = 3, labels = c("a", "b", ""), rel_widths = c(4, 4, 1)
)
lower <- plot_grid(
  recolor(alt_plots$k_hist) + labs(title = NULL),
  panel_d + theme(legend.position = "none"),
  ncol = 2, labels = c("c", "d"), rel_widths = c(1, 1.3)
)
ggsave(file.path(FIGURE_DIR, "fig2_tcgacptac.pdf"),
       plot_grid(upper, lower, nrow = 2),
       width = 6.5, height = 6)
message("Saved fig2_tcgacptac.pdf")
# NOTE: the k x alpha BO tuning surface (fig_bo_heat) is no longer a
# main-text panel; relocate it to the SI (model-selection sensitivity) in a follow-up.

# ── SI S1: Convergence trajectories ──────────────────────────────────────
# Earlier versions plotted five randomly chosen initializations on a linear
# 10^-5 scale. That was misleading in two ways, both raised by the PI: the
# five curves were not representative (across all 100 runs the final relative
# decrease spans roughly three orders of magnitude, so a different five produce
# a completely different-looking panel and y-range), and a linear axis cannot
# display that spread. We therefore plot ALL runs on a log axis, which shows
# the monotone decrease and the true between-initialization variability at once.
desurv_seed_fits <- load_precomputed("desurv_seed_fits_tcgacptac")

lossit <- dplyr::bind_rows(lapply(seq_along(desurv_seed_fits$fits), function(i) {
  v <- desurv_seed_fits$fits[[i]]$lossit
  data.frame(iter = seq_along(v), init = i,
             rel_dec = (v[1] - v) / v[1])
})) |>
  dplyr::filter(iter < 5000, rel_dec > 0)

.med <- dplyr::group_by(lossit, iter) |>
  dplyr::summarise(rel_dec = median(rel_dec), .groups = "drop") |>
  dplyr::filter(iter <= median(sapply(desurv_seed_fits$fits, function(x) length(x$lossit))))

ggsave(
  file.path(FIGURE_DIR, "si_fig_converge_tcgacptac.pdf"),
  ggplot(lossit, aes(x = iter, y = rel_dec, group = init)) +
    geom_line(linewidth = 0.25, alpha = 0.30, colour = "grey45") +
    geom_line(data = .med, aes(group = 1), linewidth = 0.9, colour = "#2C6FBB") +
    scale_y_log10(labels = function(y) formatC(y, format = "e", digits = 0)) +
    annotation_logticks(sides = "l", size = 0.2) +
    labs(x = "Iteration",
         y = expression("Relative decrease in objective, " * (l[0] - l[t]) / l[0]),
         caption = paste0("All ", length(desurv_seed_fits$fits),
                          " initializations (grey); median trajectory (blue)")) +
    theme_nature(),
  width = 5, height = 3.2
)
message("Saved si_fig_converge_tcgacptac.pdf")

# ── SI S2: Null and mixed simulation scenarios ────────────────────────────
null_f  <- lapply(null_plots,  set_fig_font, size = 10)
mixed_f <- lapply(mixed_plots, set_fig_font, size = 10)

null_cindex <- null_f$cindex_box + labs(title = NULL) + scale_y_continuous(limits = c(0.3, 1))
null_cindex$data$scenario_id <- "Null scenario (β = 0)"
null_khist  <- null_f$k_hist + labs(title = NULL)
null_khist$data$scenario_id  <- "Null scenario (β = 0)"
null_row <- plot_grid(null_cindex, null_khist, ncol = 2, labels = c("a", "b"), rel_widths = c(1, 1.2))

mixed_cindex <- mixed_f$cindex_box + labs(title = NULL) + scale_y_continuous(limits = c(0.3, 1)) +
  theme(plot.clip = "off", plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))
mixed_cindex$data$scenario_id <- "Mixed scenario"
mixed_prec <- if (!is.null(mixed_f$precision_breakdown)) {
  mixed_f$precision_breakdown + labs(title = NULL)
} else {
  mixed_f$precision_box + labs(title = NULL)
}
mixed_row1 <- plot_grid(mixed_cindex, mixed_prec, ncol = 2, labels = c("c", "d"), rel_widths = c(1, 2))

mixed_beta  <- if (!is.null(mixed_f$matched_beta_box)) mixed_f$matched_beta_box + labs(title = NULL) else NULL
mixed_khist <- mixed_f$k_hist + labs(title = NULL)
mixed_khist$data$scenario_id <- "Mixed scenario"

fig_null_mixed <- if (!is.null(mixed_beta)) {
  plot_grid(null_row, mixed_row1,
            plot_grid(mixed_beta, mixed_khist, ncol = 2, labels = c("e", "f"), rel_widths = c(1, 1.2)),
            nrow = 3, rel_heights = c(1, 1, 1))
} else {
  plot_grid(null_row, mixed_row1,
            plot_grid(mixed_khist, ncol = 1, labels = "e"),
            nrow = 3, rel_heights = c(1, 1, 0.7))
}
ggsave(file.path(FIGURE_DIR, "si_fig_sim_null_mixed_tcgacptac.pdf"), fig_null_mixed, width = 6.5, height = 7)
message("Saved si_fig_sim_null_mixed_tcgacptac.pdf")

message("=== Step 8c complete ===")
