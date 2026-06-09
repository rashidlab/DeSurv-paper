#!/usr/bin/env Rscript
# code/13_defense_fig3c_variance.R — Defense version of main-text Fig 3C
#
# Rebuilds the "variance explained vs survival contribution" scatter with a
# CORRECTED x-axis. The published Fig 3C plots per-factor reconstruction ENERGY
# of an UNCENTERED, rank-normalized matrix ( ||W_j H_j||^2 / ||X||^2 ), which is
# dominated by the grand-mean/baseline level and therefore deflates factors that
# encode genuine cross-sample VARIATION (notably the exocrine factor N2: 6.9%).
#
# This version puts true cross-sample variance on the x-axis: each gene is
# centered across samples before the per-factor share is computed. On this
# basis the exocrine factor moves 6.9% -> ~39% and the prognostic DeSurv factor
# D1 sits at ~1% — the "prognostic signal is low-variance" thesis made explicit.
#
# y-axis (Type III / leave-one-out partial log-likelihood) matches the original.
#
# Run:  Rscript code/13_defense_fig3c_variance.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(survival)
})

OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read1 <- function(name) readRDS(file.path("results", paste0(name, ".rds")))

fit_std    <- read1("fit_std_desurvk_tcgacptac")
fit_desurv <- read1("tar_fit_desurv_tcgacptac")
dat        <- read1("tar_data_filtered_tcgacptac")

X     <- as.matrix(dat$ex)
Xc    <- X - rowMeans(X)              # center each gene across samples
time  <- dat$sampInfo$time
event <- dat$sampInfo$event
var_den <- sum(Xc^2)                   # total cross-sample variance (SS)

# ── Per-factor centered variance share + survival contribution ──────────────
build_df <- function(W, H, method) {
  W <- as.matrix(W); H <- as.matrix(H); k <- ncol(W)

  # x: centered variance explained by factor j's rank-1 component
  var_exp <- vapply(seq_len(k), function(j) {
    comp  <- W[, j] %o% H[j, ]
    compc <- comp - rowMeans(comp)
    sum(compc^2) / var_den
  }, numeric(1))

  # y: leave-one-out partial log-likelihood (full vs k-1 factor Cox model)
  XtW     <- t(X) %*% W                 # subjects x factors
  ll_full <- coxph(Surv(time, event) ~ XtW)$loglik[2]
  delta_ll <- vapply(seq_len(k), function(j) {
    XtW_mj  <- XtW[, -j, drop = FALSE]
    reduced <- if (ncol(XtW_mj) == 0L) coxph(Surv(time, event) ~ 1)
               else coxph(Surv(time, event) ~ XtW_mj)
    ll_full - reduced$loglik[2]
  }, numeric(1))

  data.frame(method = method, factor = seq_len(k),
             variance_explained = var_exp, delta_loglik = delta_ll)
}

df_nmf    <- build_df(fit_std$W,    fit_std$H,    "NMF")
df_desurv <- build_df(fit_desurv$W, fit_desurv$H, "DeSurv")

df_plot <- rbind(df_nmf, df_desurv)
df_plot$factor_label <- ifelse(df_plot$method == "NMF",
                               paste0("N", df_plot$factor),
                               paste0("D", df_plot$factor))

cat("\nPer-factor values (centered variance x-axis):\n")
print(within(df_plot, {
  variance_explained <- round(variance_explained, 3)
  delta_loglik       <- round(delta_loglik, 2)
}))

# ── Plot (deck-scale fonts; same visual grammar as published 3C) ────────────
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
    x = "Cross-sample variance explained\n(gene-centered)",
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

ggsave(file.path(OUT, "new_3c_variance_explained.pdf"), fig,
       width = 8, height = 6.2, device = cairo_pdf)
ggsave(file.path(OUT, "new_3c_variance_explained.png"), fig,
       width = 8, height = 6.2, dpi = 300, bg = "white", type = "cairo")
message("Saved new_3c_variance_explained (8x6.2 in)")
