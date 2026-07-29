#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Standalone builder for the COMPLETE Fig 3 (external validation), assembled
# from cached objects in the flat results/ store WITHOUT the DeSurv package,
# DiceKriging, or a ggplot2 downgrade.
#
# This works because every input Fig 3 needs is already cached in results/
# (RESULTS_DIR is flat "results/", per code/00_helpers.R:28 — NOT
# results/precomputed/, which the stale CLAUDE.md doc references):
#   - tar_fit_desurv_tcgacptac.rds   (desurv_fit; $W only — loads w/o DeSurv)
#   - fit_std_desurvk_tcgacptac.rds  (matched-rank standard NMF fit)
#   - data_val_filtered_tcgacptac.rds (5 external validation cohorts)
#   - fig_median_survival_desurv_tcgacptac.rds (ggsurvplot; panel B)
#   - desurv_vs_supervised_tuned.rds$axis_decomposition (panel C values)
# compute_hrs() (the forest engine) is pure base-R/survival — no DeSurv — so
# panel A rebuilds locally. The authoritative build still lives in
# code/09a_figures.R (which additionally recomputes cutpoint KMs via DeSurv);
# this util reproduces its Fig-3 output faithfully for local iteration.
#
# Layout mirrors 09a: plot_grid(forest[A], km_block[B/C], rel_widths=c(1.4,1)).
# Output: figures/standalone_preview/fig3_full.pdf
# ---------------------------------------------------------------------------

suppressMessages({
  library(ggplot2); library(cowplot); library(survminer)
  library(survival); library(dplyr); library(gtable)
})
# The cached KM risk table styles its strata labels with element_markdown, which
# needs ggtext loaded to render. ggtext is not a declared pipeline dependency, so
# load it defensively: if absent, we drop the number-at-risk table (panel B keeps
# the survival curve). Add "ggtext" to code/01_install.R to always include it.
has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
if (has_ggtext) suppressMessages({
  library(ggtext)
})

out_dir <- "figures/standalone_preview"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- inputs -----------------------------------------------------------------
tar_fit_desurv    <- readRDS("results/tar_fit_desurv_tcgacptac.rds")
fit_std_desurvk   <- readRDS("results/fit_std_desurvk_tcgacptac.rds")
data_val_filtered <- readRDS("results/data_val_filtered_tcgacptac.rds")
km                <- readRDS("results/fig_median_survival_desurv_tcgacptac.rds")
axd               <- readRDS("results/desurv_vs_supervised_tuned.rds")$axis_decomposition

# --- compute_hrs: frozen-LP per-factor Cox HR per cohort (DeSurv-free) -------
# (verbatim logic from R/figure_plot_helpers.R:303)
compute_hrs <- function(data_val_filtered, fit, method) {
  df <- list()
  for (i in seq_along(data_val_filtered)) {
    dat  <- data_val_filtered[[i]]
    keep <- intersect(rownames(dat$ex), rownames(fit$W))
    XtW  <- t(dat$ex[keep, ]) %*% fit$W[keep, ]
    hr <- lower <- upper <- numeric(ncol(XtW))
    for (j in seq_len(ncol(XtW))) {
      s <- summary(coxph(Surv(dat$sampInfo$time, dat$sampInfo$event) ~ scale(XtW[, j])))
      hr[j] <- exp(s$coefficients[1]); lower[j] <- s$conf.int[3]; upper[j] <- s$conf.int[4]
    }
    df[[i]] <- data.frame(factor = seq_len(ncol(XtW)), HR = hr, lower = lower,
                          upper = upper, dataset = dat$dataname)
  }
  out <- do.call(rbind, df); out$method <- method; out
}

df <- rbind(compute_hrs(data_val_filtered, tar_fit_desurv, "DeSurv"),
            compute_hrs(data_val_filtered, fit_std_desurvk, "NMF"))
df$factor_name <- ifelse(df$method == "DeSurv", paste0("D", df$factor), paste0("N", df$factor))

# --- Panel A: forest (verbatim styling from 09a:384-458) --------------------
d <- df
d$dataset <- dplyr::recode(d$dataset, "Puleo_array" = "Puleo", "Moffitt_GEO_array" = "Moffitt",
                           "PACA_AU_seq" = "PACA seq", "PACA_AU_array" = "PACA array")
pooled <- d %>%
  dplyr::mutate(logHR = log(HR), se = (log(upper) - log(lower)) / (2 * 1.96), w = 1 / se^2) %>%
  dplyr::group_by(method, factor_name) %>%
  dplyr::summarise(pooled_logHR = sum(logHR * w) / sum(w), pooled_se = sqrt(1 / sum(w)),
                   .groups = "drop") %>%
  dplyr::mutate(HR = exp(pooled_logHR), lower = exp(pooled_logHR - 1.96 * pooled_se),
                upper = exp(pooled_logHR + 1.96 * pooled_se), dataset = "Pooled")
keep_cols <- c("factor_name", "HR", "lower", "upper", "dataset", "method", "row_type")
d$row_type <- "cohort"; pooled$row_type <- "pooled"
all_data <- rbind(d[, keep_cols], pooled[, keep_cols])
all_data$dataset <- factor(all_data$dataset,
  levels = c("Puleo", "PACA seq", "PACA array", "Moffitt", "Dijk", "Pooled"))
cohort_cols   <- c("Dijk"="#E69F00","Moffitt"="#56B4E9","PACA array"="#009E73",
                   "PACA seq"="#0072B2","Puleo"="#CC79A7","Pooled"="#000000")
cohort_shapes <- c("Dijk"=16,"Moffitt"=17,"PACA array"=15,"PACA seq"=25,"Puleo"=8,"Pooled"=18)
cohort_sizes  <- c("Dijk"=2.2,"Moffitt"=2.2,"PACA array"=2.2,"PACA seq"=2.2,"Puleo"=2.2,"Pooled"=4)
base_size <- 8; dodge_w <- 0.6
desurv_dat <- all_data[all_data$method == "DeSurv", ]
desurv_dat$factor_name <- factor(desurv_dat$factor_name, levels = c("D1","D2","D3"))
nmf_dat <- all_data[all_data$method == "NMF", ]
nmf_dat$factor_name <- factor(nmf_dat$factor_name, levels = c("N1","N2","N3"))

make_forest_panel <- function(dat, title, highlight_level) {
  hl <- data.frame(factor_name = factor(highlight_level, levels = levels(dat$factor_name)), x = 1)
  ggplot(dat, aes(x = HR, y = factor_name)) +
    geom_tile(data = hl, aes(x = x, y = factor_name), width = 100, height = 0.8,
              fill = "steelblue", alpha = 0.08, inherit.aes = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lower, xmax = upper, colour = dataset, linewidth = row_type),
                  orientation = "y", width = 0, position = position_dodge(width = dodge_w)) +
    geom_point(aes(colour = dataset, shape = dataset, size = dataset),
               position = position_dodge(width = dodge_w)) +
    scale_linewidth_manual(values = c("cohort" = 0.6, "pooled" = 1.2), guide = "none") +
    scale_colour_manual(values = cohort_cols, name = NULL, drop = FALSE) +
    scale_shape_manual(values = cohort_shapes, name = NULL, drop = FALSE) +
    scale_size_manual(values = cohort_sizes, name = NULL, drop = FALSE) +
    scale_x_log10(limits = c(0.35, 3)) +
    labs(x = "Hazard ratio (95% CI)", y = NULL, title = title) +
    theme_classic(base_size = base_size) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
          axis.text = element_text(size = 8, face = "bold"),
          axis.title.x = element_text(size = 9), plot.margin = margin(2, 6, 2, 4))
}
p_leg <- ggplot(all_data, aes(HR, factor_name, colour = dataset, shape = dataset, size = dataset)) +
  geom_point() +
  scale_colour_manual(values = cohort_cols, name = NULL) +
  scale_shape_manual(values = cohort_shapes, name = NULL) +
  scale_size_manual(values = cohort_sizes, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.5)),
         shape = guide_legend(nrow = 1), size = guide_legend(nrow = 1)) +
  theme_void(base_size = base_size) + theme(legend.position = "bottom",
                                            legend.text = element_text(size = 8))
forest_legend <- gtable::gtable_filter(ggplotGrob(p_leg), "guide-box")
plot_forest <- plot_grid(
  plot_grid(make_forest_panel(desurv_dat, "DeSurv (D1-D3)", "D1"),
            make_forest_panel(nmf_dat, "Standard NMF (N1-N3)", "N1"),
            ncol = 2, align = "hv", axis = "tb"),
  forest_legend, nrow = 2, rel_heights = c(15, 1))

# --- Panel B: DeSurv KM, plot stacked over number-at-risk table --------------
km_plot <- km$plot + ggtitle("DeSurv") + theme(plot.title = element_text(size = 9, hjust = 0.5))
km_b <- if (has_ggtext) {
  plot_grid(km_plot, km$table + theme(plot.title = element_text(size = 8)),
            ncol = 1, rel_heights = c(4, 1.2))
} else {
  message("ggtext not available: dropping number-at-risk table from panel B")
  km_plot
}

# --- Panel C: tuned-supervised vs DeSurv-program correspondence heatmap ------
hm_df <- data.frame(
  method  = factor(rep(c("Supervised PCA", "Penalized Cox"), each = 3),
                   levels = c("Supervised PCA", "Penalized Cox")),
  program = factor(rep(c("D1", "D2", "D3"), 2), levels = c("D1", "D2", "D3")),
  r = c(abs(as.numeric(axd["Supervised PCA", c("D1","D2","D3")])),
        abs(as.numeric(axd["Sparse Cox",     c("D1","D2","D3")]))))
fig_c <- ggplot(hm_df, aes(program, method, fill = r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.1) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1),
                      name = expression("|" * italic(r) * "|")) +
  labs(x = NULL, y = NULL, title = "Supervised score vs DeSurv program") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0.5), panel.grid = element_blank(),
        axis.text = element_text(color = "black"), legend.position = "right")

km_block <- plot_grid(km_b, fig_c, ncol = 1, labels = c("B", "C"), rel_heights = c(5, 4.2))

# --- compose full Fig 3 -----------------------------------------------------
fig3 <- plot_grid(plot_forest, km_block, ncol = 2, rel_widths = c(1.4, 1),
                  labels = c("A", ""))
ggsave(file.path(out_dir, "fig3_full.pdf"), fig3, width = 9.2, height = 6.4)
cat("Wrote", file.path(out_dir, "fig3_full.pdf"),
    "- full Fig 3 (A forest + B KM + C supervised heatmap), built locally without DeSurv.\n")
