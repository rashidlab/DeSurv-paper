#!/usr/bin/env Rscript
# code/16_defense_fig4_panels.R — Defense version of main-text Fig 4 (external validation)
#
# Section 5 of the deck has two slides:
#   Slide "Generalization to 5 cohorts"     -> fig4_forest      (Fig 4A)
#   Slide "Risk stratification in the clinic" -> fig4_km_desurv / fig4_km_nmf (Fig 4B/4C)
#
# The paper composite (code/09a_figures.R, fig4_tcgacptac.pdf) packs the forest
# and both KM panels into a 7x4.8 in figure with base_size 8 — illegible at
# slide scale. This rebuilds each panel standalone from the same cached fits and
# validation data with enlarged fonts, matching the defense-deck conventions
# (no A/B labels; KM panels saved separately so NMF can be revealed on click).
#
# Run:  Rscript code/16_defense_fig4_panels.R

message("=== Defense Fig 4 panels (external validation) ===")
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(survival)
  library(survminer)
  library(dplyr)
})

source("code/00_helpers.R")
source("R/get_top_genes.R")
source("R/fit_cox_model.R")
source("R/cv_grid_helpers.R")
source("R/figure_plot_helpers.R")

OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── Load prerequisites (same cached objects as 09a) ─────────────────────────
tar_fit_desurv    <- load_precomputed("tar_fit_desurv_tcgacptac")
fit_std_desurvk   <- load_precomputed("fit_std_desurvk_tcgacptac")
data_val_filtered <- load_precomputed("data_val_filtered_tcgacptac")
tar_params_best   <- load_precomputed("tar_params_best_tcgacptac")

# ════════════════════════════════════════════════════════════════════════════
# Fig 4A — per-factor HR forest plot across the 5 validation cohorts (+ pooled)
# ════════════════════════════════════════════════════════════════════════════
desurv_df <- compute_hrs(data_val_filtered, tar_fit_desurv, "DeSurv")
nmf_df    <- compute_hrs(data_val_filtered, fit_std_desurvk, "NMF")
df <- rbind(desurv_df, nmf_df)
df$factor_name <- dplyr::case_when(
  df$method == "DeSurv" ~ paste0("D", df$factor),
  df$method == "NMF"    ~ paste0("N", df$factor),
  TRUE                  ~ as.character(df$factor)
)
df$dataset <- dplyr::recode(df$dataset,
  "Puleo_array" = "Puleo", "Moffitt_GEO_array" = "Moffitt",
  "PACA_AU_seq" = "PACA seq", "PACA_AU_array" = "PACA array")

# Inverse-variance pooled HR per factor (fixed-effect meta-analysis).
pooled <- df %>%
  dplyr::mutate(logHR = log(HR), se = (log(upper) - log(lower)) / (2 * 1.96),
                w = 1 / se^2) %>%
  dplyr::group_by(method, factor_name) %>%
  dplyr::summarise(pooled_logHR = sum(logHR * w) / sum(w),
                   pooled_se    = sqrt(1 / sum(w)), .groups = "drop") %>%
  dplyr::mutate(HR = exp(pooled_logHR),
                lower = exp(pooled_logHR - 1.96 * pooled_se),
                upper = exp(pooled_logHR + 1.96 * pooled_se),
                dataset = "Pooled")
keep_cols <- c("factor_name", "HR", "lower", "upper", "dataset", "method", "row_type")
df$row_type <- "cohort"; pooled$row_type <- "pooled"
all_data <- rbind(df[, keep_cols], pooled[, keep_cols])
all_data$method <- factor(all_data$method, levels = c("DeSurv", "NMF"),
                          labels = c("DeSurv (D1–D3)", "Standard NMF (N1–N3)"))
all_data$dataset <- factor(all_data$dataset,
  levels = c("Puleo", "PACA seq", "PACA array", "Moffitt", "Dijk", "Pooled"))

cohort_cols   <- c("Dijk" = "#E69F00", "Moffitt" = "#56B4E9", "PACA array" = "#009E73",
                   "PACA seq" = "#0072B2", "Puleo" = "#CC79A7", "Pooled" = "#000000")
cohort_shapes <- c("Dijk" = 16, "Moffitt" = 17, "PACA array" = 15,
                   "PACA seq" = 25, "Puleo" = 8, "Pooled" = 18)
cohort_sizes  <- c("Dijk" = 3.4, "Moffitt" = 3.4, "PACA array" = 3.4,
                   "PACA seq" = 3.4, "Puleo" = 3.4, "Pooled" = 5.8)
BASE <- 18; dodge_w <- 0.6

desurv_dat <- all_data[all_data$method == "DeSurv (D1–D3)", ]
desurv_dat$factor_name <- factor(desurv_dat$factor_name, levels = c("D1", "D2", "D3"))
nmf_dat <- all_data[all_data$method == "Standard NMF (N1–N3)", ]
nmf_dat$factor_name <- factor(nmf_dat$factor_name, levels = c("N1", "N2", "N3"))

# Standalone single-method panel (no legend — a shared legend strip sits below
# both on the slide). `box_level` optionally draws a red rectangle around one
# factor row (the deck reveals this on a later click to spotlight D1).
make_forest_panel <- function(dat, title, highlight_level, box_level = NULL) {
  hl <- data.frame(factor_name = factor(highlight_level, levels = levels(dat$factor_name)),
                   x = 1)
  p <- ggplot(dat, aes(x = HR, y = factor_name)) +
    geom_tile(data = hl, aes(x = x, y = factor_name), width = 100, height = 0.8,
              fill = "steelblue", alpha = 0.08, inherit.aes = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
    geom_errorbar(aes(xmin = lower, xmax = upper, colour = dataset, linewidth = row_type),
                  orientation = "y", width = 0, position = position_dodge(width = dodge_w)) +
    geom_point(aes(colour = dataset, shape = dataset, size = dataset),
               position = position_dodge(width = dodge_w)) +
    scale_linewidth_manual(values = c("cohort" = 0.9, "pooled" = 1.8), guide = "none") +
    scale_colour_manual(values = cohort_cols, name = NULL, drop = FALSE) +
    scale_shape_manual(values = cohort_shapes, name = NULL, drop = FALSE) +
    scale_size_manual(values = cohort_sizes, name = NULL, drop = FALSE) +
    scale_x_log10(limits = c(0.35, 3), breaks = c(0.5, 1, 2)) +
    labs(x = "Hazard ratio (95% CI)", y = NULL, title = title) +
    theme_classic(base_size = BASE) +
    theme(legend.position = "none",
          plot.title   = element_text(face = "bold", size = BASE + 2, hjust = 0.5),
          axis.text.y  = element_text(size = BASE, face = "bold"),
          axis.text.x  = element_text(size = BASE - 1, face = "bold"),
          axis.title.x = element_text(size = BASE),
          plot.margin  = margin(4, 10, 12, 6))
  if (!is.null(box_level)) {
    yi <- match(box_level, levels(dat$factor_name))
    p <- p + annotate("rect", xmin = 0.35, xmax = 3, ymin = yi - 0.45, ymax = yi + 0.45,
                      colour = "red", fill = NA, linewidth = 1.6)
  }
  p
}

# Shared horizontal cohort legend, saved as its own strip for the slide footer.
p_leg <- ggplot(all_data, aes(x = HR, y = factor_name,
                              colour = dataset, shape = dataset, size = dataset)) +
  geom_point() +
  scale_colour_manual(values = cohort_cols, name = NULL) +
  scale_shape_manual(values = cohort_shapes, name = NULL) +
  scale_size_manual(values = cohort_sizes, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4.5)),
         shape  = guide_legend(nrow = 1), size = guide_legend(nrow = 1)) +
  theme_void(base_size = BASE) +
  theme(legend.position = "bottom", legend.text = element_text(size = BASE))
forest_legend <- gtable::gtable_filter(ggplotGrob(p_leg), "guide-box")

save_forest <- function(p, stem, w = 5.2, h = 4.8) {
  ggsave(file.path(OUT, paste0(stem, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(stem, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}

save_forest(make_forest_panel(nmf_dat, "Standard NMF (N1–N3)", "N1"),
            "fig4_forest_nmf")
save_forest(make_forest_panel(desurv_dat, "DeSurv (D1–D3)", "D1"),
            "fig4_forest_desurv")
save_forest(make_forest_panel(desurv_dat, "DeSurv (D1–D3)", "D1", box_level = "D1"),
            "fig4_forest_desurv_box")
save_forest(cowplot::ggdraw(forest_legend), "fig4_forest_legend", w = 9, h = 0.55)
message("Saved fig4_forest_{nmf,desurv,desurv_box,legend}")

# ════════════════════════════════════════════════════════════════════════════
# Fig 4B/4C — pooled-validation KM curves at the CV-selected log-rank cutpoint
# ════════════════════════════════════════════════════════════════════════════
km_desurv <- load_precomputed("fig_median_survival_desurv_tcgacptac")
km_nmf    <- load_precomputed("fig_median_survival_std_desurvk_tcgacptac")

KM <- 20
theme_km <- theme_classic(base_size = KM) +
  theme(plot.title = element_text(face = "bold", size = KM + 2, hjust = 0.5),
        plot.margin = margin(6, 12, 4, 12))

# Matched left/right x-expansion on BOTH the curve and the risk table keeps the
# panels aligned and insets x=0 so the leftmost at-risk count no longer collides
# with the y-axis. Bumps every baked-in text layer (in-plot HR annotation; the
# at-risk numbers) so they read at slide scale.
x_common <- scale_x_continuous(breaks = seq(0, 150, 25),
                               expand = expansion(mult = c(0.10, 0.03)))
bump_text <- function(gg, sz) {
  for (i in seq_along(gg$layers))
    if (inherits(gg$layers[[i]]$geom, "GeomText"))
      gg$layers[[i]]$aes_params$size <- sz
  gg
}

stack_surv <- function(surv_obj, title) {
  # The baked-in HR annotation sits at y=0.85; on the shorter panel its top line
  # crowds the plot border, so drop it a little.
  pp <- surv_obj$plot
  for (i in seq_along(pp$layers))
    if (inherits(pp$layers[[i]]$geom, "GeomText") && !is.null(pp$layers[[i]]$data$y))
      pp$layers[[i]]$data$y <- 0.72
  p <- bump_text(pp, 6) + theme_km + x_common +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.title.y = element_text(size = KM),
          axis.text    = element_text(size = KM),
          plot.margin  = margin(4, 12, 0, 12)) +
    labs(title = title)
  t <- bump_text(surv_obj$table, 6) + theme_km + x_common +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text    = element_text(size = KM - 3),
          text         = element_text(size = KM),
          axis.title.x = element_text(size = KM),
          plot.title   = element_text(size = KM),
          plot.margin  = margin(4, 12, 4, 12)) +
    labs(x = "Time (months)")
  plot_grid(p, t, ncol = 1, rel_heights = c(2.1, 1.6), align = "v", axis = "lr")
}

km_legend_plot <- ggplot(
  data.frame(x = 1:2, y = 1:2,
             group = factor(c("Low", "High"), levels = c("Low", "High"))),
  aes(x = x, y = y, colour = group)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("Low" = "violetred2", "High" = "turquoise4"),
                      name = "Risk group") +
  theme_void(base_size = KM) +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = KM),
        legend.title = element_text(size = KM)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 2)))
km_legend_grob <- gtable::gtable_filter(ggplotGrob(km_legend_plot), "guide-box")

save_km <- function(surv_obj, title, stem) {
  # Trailing NULL row adds bottom whitespace so the legend's descenders aren't
  # clipped at the image edge.
  panel <- plot_grid(stack_surv(surv_obj, title), ggdraw(km_legend_grob), NULL,
                     ncol = 1, rel_heights = c(10, 1, 0.3))
  # Landscape-ish (aspect < 1) so each panel fills a half-column without its
  # bottom running off the slide.
  ggsave(file.path(OUT, paste0(stem, ".pdf")), panel, width = 5.8, height = 5.2,
         device = cairo_pdf)
  ggsave(file.path(OUT, paste0(stem, ".png")), panel, width = 5.8, height = 5.2,
         dpi = 300, bg = "white")
}

save_km(km_desurv, "DeSurv", "fig4_km_desurv")
save_km(km_nmf,    "Standard NMF", "fig4_km_nmf")
message("Saved fig4_km_desurv and fig4_km_nmf (5.2x5.6 in, enlarged fonts)")
message("=== Defense Fig 4 panels complete ===")
