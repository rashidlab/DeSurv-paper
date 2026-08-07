#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Generates figures/fig_treated_context.pdf (manuscript Fig. 4, \label{fig:treated}):
# a two-panel translational figure built entirely from the tracked cache
# results/treated_cohort_stats.rds (no restricted data needed here).
#   a: O'Kane/COMPASS prognostic transportability of D1 in treated metastatic
#      PDAC -- per ACTUAL treatment arm (FFX, GA, GA/experimental; never pooled)
#      plus the arm-stratified common estimate (diamond).
#   b: Linehan paired pre/post biopsies -- D2 (proCAF) change after treatment.
# The former Panel B (bulk trial x arm survival forest) was removed; see the
# rationale block below.
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
## PANEL B (bulk trial x arm survival forest) REMOVED, 2026-08-06.
## Rationale. The treated section had three components of very unequal quality:
##   (i)  O'Kane D1 transportability -- arm-stratified WITHIN one cohort, clean;
##   (ii) paired pre/post D2 -- within-patient, clean, and a statement about
##        program LEVELS, not about survival;
##   (iii) the bulk trial x arm survival forest -- 7 cells of 14-52 patients,
##        nothing surviving BH correction.
## Only (iii) produced a D2 survival claim in treated disease, and that claim
## formally contradicts the untreated validation: D2 is a well-powered null
## untreated (HR 1.02, 0.92-1.12, n=570/388 events) but protective treated
## (HR 0.76, 0.63-0.91, strata = trial x arm), differing at P = 0.0065. The
## mundane alternative -- resected untreated versus biopsy treated tissue
## differing in tumour/stroma content -- cannot be excluded, because neither
## treated cohort carries a purity field. D2 also has an elastic-net coefficient
## of exactly zero in the trained model, so a D2 survival claim is the most
## exposed statement available. Dropping (iii) removes the contradiction at its
## source instead of arguing it away. (i) and (ii) are retained and neither
## conflicts with the untreated results: D1 attenuates (P = 0.073, n.s.), which
## is the stage attenuation PurIST/Moffitt/Bailey also show.
## The underlying analysis is retained in code/24_treated_arm_level.R and
## results/treated_arm_level_stats.rds as reviewer-response material.

## --- Panel C: Linehan paired pre/post D2 ----------------------------------
pp <- t$paired$points
long <- data.frame(pair = rep(seq_len(nrow(pp)), 2),
                   time = factor(rep(c("Pre", "Post"), each = nrow(pp)), levels = c("Pre", "Post")),
                   D2 = c(pp$D2_pre, pp$D2_post))
# PI asked for the mean trajectory to be drawn over the per-patient lines.
mn <- data.frame(time = factor(c("Pre", "Post"), levels = c("Pre", "Post")),
                 D2   = c(mean(pp$D2_pre), mean(pp$D2_post)))
pC <- ggplot(long, aes(time, D2, group = pair)) +
  geom_line(alpha = 0.30, colour = "grey55") +
  geom_point(alpha = 0.55, size = 1.3, colour = blue) +
  geom_line(data = mn, aes(group = 1), colour = "#B2182B", linewidth = 1.4) +
  geom_point(data = mn, aes(group = 1), colour = "#B2182B", size = 2.6) +
  annotate("text", x = 0.7, y = max(long$D2), hjust = 0, vjust = 1, size = 2.9,
           label = sprintf("paired Wilcoxon\nP = %.3f, n = %d", t$paired$wilcox_p, nrow(pp))) +
  labs(x = NULL, y = "D2 (proCAF) score (z)",
       title = "Linehan paired biopsies:\nD2 change after treatment") +
  theme_classic(base_size = 9) + theme(plot.title = element_text(size = 8.5, face = "bold"))

## --- compose: a O'Kane arm-level forest | b Linehan paired pre/post ---------
fig <- plot_grid(pA, pC, ncol = 2, labels = c("a", "b"), label_size = 12, rel_widths = c(1.2, 0.8))
ggsave("figures/fig_treated_context.pdf", fig, width = 7.2, height = 3.4)
cat("Wrote figures/fig_treated_context.pdf (Fig. 4: a OKane arm-stratified D1, b paired D2)\n")
