# ---------------------------------------------------------------------------
# Shared publication theme + palettes for the DeSurv display items, styled to
# Nature Cancer conventions: sans-serif type at small final sizes, white
# background, no panel box, minimal gridlines, thin strokes, and a colorblind-
# safe (Okabe-Ito) palette with fixed semantic colors reused across figures.
#
# Single source of truth: sourced by code/09a_figures.R, 09b_si_figures.R,
# 09c_sim_figures.R, code/util_fig2_standalone.R and R/figure_plot_helpers.R so
# every figure regenerates in one consistent style.
# ---------------------------------------------------------------------------
suppressMessages({ library(ggplot2) })

# --- Okabe-Ito colorblind-safe base palette --------------------------------
okabe_ito <- c(black = "#000000", orange = "#E69F00", skyblue = "#56B4E9",
               green = "#009E73", yellow = "#F0E442", blue = "#0072B2",
               vermillion = "#D55E00", purple = "#CC79A7", grey = "#999999")

# --- Fixed semantic colors (reused across every figure) --------------------
# Method contrast (DeSurv vs standard NMF): colorblind-safe blue vs vermillion.
desurv_method_cols <- c("DeSurv" = "#0072B2", "NMF" = "#D55E00")
# Validation datasets (kept from the existing forest palette).
desurv_cohort_cols <- c("Dijk" = "#E69F00", "Moffitt" = "#56B4E9",
                        "PACA array" = "#009E73", "PACA seq" = "#0072B2",
                        "Puleo" = "#CC79A7", "Pooled" = "#000000")
# Single-series accent (e.g., GATA6, treated-cohort points).
desurv_accent <- "#0072B2"
# Perceptually-ordered diverging ramp for enrichment/correlation heatmaps,
# centered at 0 (RdBu; robust in print and to common color-vision deficiencies).
desurv_diverging <- function(n = 100)
  grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdBu")))(n)

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

# Nature Cancer column widths (inches) for ggsave(width = ...).
nat_width <- function(cols = c("single", "onehalf", "double")) {
  cols <- match.arg(cols)
  mm <- c(single = 88, onehalf = 120, double = 180)[[cols]]
  mm / 25.4
}
