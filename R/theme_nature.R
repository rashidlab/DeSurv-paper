# ---------------------------------------------------------------------------
# Shared publication theme + palettes for the DeSurv display items, styled to
# Nature Cancer conventions: sans-serif type at small final sizes, white
# background, no panel box, minimal gridlines, thin strokes, and a colorblind-
# safe (Okabe-Ito) palette with fixed semantic colors reused across figures.
#
# Single source of truth: sourced by code/09a_figures.R, 09b_si_figures.R,
# 09c_sim_figures.R, code/util_fig2_standalone.R and R/figure_plot_helpers.R so
# every figure regenerates in one consistent style.
#
# REFINED for the clarity pass (2026-07): added desurv_risk_cols (KM Low/High),
# desurv_discrete() (small qualitative series, e.g. random-init traces) and
# desurv_sequential() (magnitude-only |r| heatmaps), plus nat_height() and a
# theme_nature_heatmap() base_size hook so pheatmap grobs match the ggplot type.
# Nothing here changes any data; colors and type only.
# ---------------------------------------------------------------------------
suppressMessages({ library(ggplot2) })

# --- Okabe-Ito colorblind-safe base palette --------------------------------
okabe_ito <- c(black = "#000000", orange = "#E69F00", skyblue = "#56B4E9",
               green = "#009E73", yellow = "#F0E442", blue = "#0072B2",
               vermillion = "#D55E00", purple = "#CC79A7", grey = "#999999")

# --- Fixed semantic colors (reused across every figure) --------------------
# Method contrast (DeSurv vs standard NMF): colorblind-safe blue vs vermillion.
# USE THESE EVERYWHERE a DeSurv/NMF contrast appears. Do not fall back to the
# ggplot defaults ("blue"/"red") or to brewer Paired ("#1f78b4"/"#e31a1c").
desurv_method_cols <- c("DeSurv" = "#0072B2", "NMF" = "#D55E00")
# Validation datasets (kept from the existing forest palette).
desurv_cohort_cols <- c("Dijk" = "#E69F00", "Moffitt" = "#56B4E9",
                        "PACA array" = "#009E73", "PACA seq" = "#0072B2",
                        "Puleo" = "#CC79A7", "Pooled" = "#000000")
# Kaplan-Meier risk dichotomy (higher DeSurv predictor = higher risk). Blue vs
# vermillion is the canonical colorblind-safe pair; it never co-occurs with the
# method contrast inside a KM panel, so there is no ambiguity. Replaces the
# ad-hoc violetred2/turquoise4 (main text) and survminer blue/red (SI) so Low
# and High read identically in every survival panel.
desurv_risk_cols <- c("Low" = "#0072B2", "High" = "#D55E00")
# Single-series accent (e.g., GATA6, treated-cohort points).
desurv_accent <- "#0072B2"

# Perceptually-ordered diverging ramp for enrichment/correlation heatmaps,
# centered at 0 (RdBu; robust in print and to common color-vision deficiencies).
# Use for signed statistics (rank-biserial, Spearman) and ALWAYS with limits
# symmetric about 0 so that 0 maps to the white midpoint.
desurv_diverging <- function(n = 100)
  grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdBu")))(n)

# Sequential single-hue ramp for magnitude-only heatmaps (|r| in [0, 1]);
# anchored on the DeSurv blue so it reads consistently with the method color.
desurv_sequential <- function(n = 100)
  grDevices::colorRampPalette(c("#f7fbff", "#0072B2"))(n)

# Qualitative scale for small categorical series (<= 6), e.g. the random-init
# convergence traces in SI S1. Okabe-Ito minus black/grey/yellow for maximum
# separability in print.
desurv_discrete <- function()
  unname(okabe_ito[c("vermillion", "orange", "green", "blue", "purple", "skyblue")])

# --- Theme -----------------------------------------------------------------
# base_size defaults to 7 pt (Nature body-text figure size). Drop-in for
# theme_classic()/theme_minimal().
theme_nature <- function(base_size = 7, base_family = "Helvetica",
                         grid = FALSE, legend = "right") {
  th <- theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text            = element_text(colour = "black"),
      axis.text       = element_text(colour = "black", size = rel(0.9)),
      axis.line       = element_line(colour = "black", linewidth = 0.3),
      axis.ticks      = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.6, "pt"),
      plot.title      = element_text(size = rel(1.0), face = "plain", hjust = 0),
      plot.subtitle   = element_text(size = rel(0.9)),
      strip.background = element_blank(),
      strip.text      = element_text(size = rel(0.95), face = "plain"),
      legend.position = legend,
      legend.key.size = unit(8, "pt"),
      legend.title    = element_text(size = rel(0.9)),
      legend.text     = element_text(size = rel(0.85)),
      legend.background = element_blank(),
      legend.key      = element_blank(),
      plot.tag        = element_text(size = rel(1.3), face = "bold"),
      plot.margin     = margin(2, 2, 2, 2)
    )
  if (grid) th <- th + theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
                             panel.grid.minor = element_blank())
  th
}

# pheatmap does not accept a ggplot theme; call this to get the matching type
# sizes to pass through fontsize / fontsize_row / fontsize_col so the Fig 2 A/B/D
# heatmap grobs sit at the same 6-7 pt as the ggplot panels around them.
theme_nature_heatmap <- function(base_size = 7)
  list(fontsize = base_size, fontsize_row = base_size, fontsize_col = base_size,
       fontsize_number = base_size, fontfamily = "Helvetica")

# Nature Cancer column widths (inches) for ggsave(width = ...).
nat_width <- function(cols = c("single", "onehalf", "double")) {
  cols <- match.arg(cols)
  mm <- c(single = 88, onehalf = 120, double = 180)[[cols]]
  mm / 25.4
}

# Height cap (inches) for ggsave(height = ...); Nature Cancer max page is 240 mm.
nat_height <- function(mm = 240) mm / 25.4
