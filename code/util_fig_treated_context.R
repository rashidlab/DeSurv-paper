#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Generates figures/fig_treated_context.pdf (manuscript Supplementary Fig. S13):
# a three-panel translational figure built entirely from the tracked cache
# results/treated_cohort_stats.rds (no restricted data needed here).
#   A: O'Kane/COMPASS prognostic transportability of D1 in treated metastatic
#      PDAC -- per ACTUAL treatment arm (FFX, GA, GA/experimental; never pooled)
#      plus the arm-stratified common estimate (diamond).
#   B: Linehan/Rash/Accept bulk treated-cohort program associations with OS
#      (D1 marginal, D2 adjusted for PurIST+DeCAF, D3 marginal).
#   C: Linehan paired pre/post biopsies -- D2 (proCAF) change after treatment.
# ---------------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(cowplot) })

t   <- readRDS("results/treated_cohort_stats.rds")
out <- "figures/standalone_preview"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
blue <- "#08519c"

## --- Panel A: O'Kane per-arm D1 forest + arm-stratified diamond ------------
oa <- t$okane$per_arm; st <- t$okane$stratified
dfA <- rbind(
  data.frame(label = as.character(oa$arm), hr = oa$hr, lo = oa$lo, hi = oa$hi, kind = "arm"),
  data.frame(label = "Arm-stratified", hr = st["hr"], lo = st["lo"], hi = st["hi"], kind = "summary"))
dfA$label <- factor(dfA$label, levels = rev(c("FFX", "GA", "GA/experimental", "Arm-stratified")))
pA <- ggplot(dfA, aes(hr, label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.18, linewidth = 0.4) +
  geom_point(aes(shape = kind, size = kind, colour = kind)) +
  scale_shape_manual(values = c(arm = 16, summary = 18), guide = "none") +
  scale_size_manual(values = c(arm = 2.4, summary = 4.2), guide = "none") +
  scale_colour_manual(values = c(arm = blue, summary = "black"), guide = "none") +
  scale_x_log10(limits = c(0.3, 1.9), breaks = c(0.4, 0.6, 1.0, 1.5)) +
  labs(x = "D1 hazard ratio per SD (95% CI)", y = NULL,
       title = "O'Kane/COMPASS: D1 transportability in\ntreated metastatic PDAC (interaction P = 0.40)") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 8.5, face = "bold"), axis.text.y = element_text(size = 8))

## --- Panel B: bulk cohorts, D1 marginal / D2 adjusted / D3 marginal --------
os <- t$os
selB <- rbind(os[os$program == "D1" & os$model == "marginal", ],
              os[os$program == "D2" & os$model == "adjusted", ],
              os[os$program == "D3" & os$model == "marginal", ])
selB$cohort  <- factor(selB$cohort, levels = rev(c("Linehan", "Rash", "Accept")))
selB$program <- factor(selB$program, levels = c("D1", "D2", "D3"),
                       labels = c("D1 (marginal)", "D2 (adjusted)", "D3 (marginal)"))
pB <- ggplot(selB, aes(hr, cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.2, linewidth = 0.4) +
  geom_point(size = 2, colour = blue) +
  facet_wrap(~program, nrow = 1) +
  scale_x_log10() +
  labs(x = "Hazard ratio per SD (95% CI)", y = NULL,
       title = "Bulk treated cohorts (Linehan, Rash, Accept): program associations with overall survival") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 8.5, face = "bold"), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA))

## --- Panel C: Linehan paired pre/post D2 ----------------------------------
pp <- t$paired$points
long <- data.frame(pair = rep(seq_len(nrow(pp)), 2),
                   time = factor(rep(c("Pre", "Post"), each = nrow(pp)), levels = c("Pre", "Post")),
                   D2 = c(pp$D2_pre, pp$D2_post))
pC <- ggplot(long, aes(time, D2, group = pair)) +
  geom_line(alpha = 0.35, colour = "grey40") +
  geom_point(alpha = 0.65, size = 1.4, colour = blue) +
  annotate("text", x = 0.7, y = max(long$D2), hjust = 0, vjust = 1, size = 2.9,
           label = sprintf("paired Wilcoxon\nP = %.3f, n = %d", t$paired$wilcox_p, nrow(pp))) +
  labs(x = NULL, y = "D2 (proCAF) score (z)",
       title = "Linehan paired biopsies:\nD2 change after treatment") +
  theme_classic(base_size = 9) + theme(plot.title = element_text(size = 8.5, face = "bold"))

## --- compose (reading order A/B/C): A O'Kane | B paired on top, C bulk bottom
top <- plot_grid(pA, pC, ncol = 2, labels = c("a", "b"), label_size = 12, rel_widths = c(1.2, 0.8))
fig <- plot_grid(top, pB, nrow = 2, labels = c("", "c"), label_size = 12, rel_heights = c(1, 0.62))
ggsave("figures/fig_treated_context.pdf", fig, width = 7.2, height = 6)
cat("Wrote figures/fig_treated_context.pdf (S13: A O'Kane per-arm, B paired D2, C bulk forest)\n")
