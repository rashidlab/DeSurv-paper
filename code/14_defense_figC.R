#!/usr/bin/env Rscript
# code/14_defense_figC.R — Defense Slide 4 "reconstruction != prognosis" figure
#
# Concept (per defense_talk_plan.md, Slide 4 thesis hook):
#   left   : an expression matrix X (genes x patients)
#   middle : X branches into MULTIPLE factorizations, each with GOOD (and nearly
#            identical) reconstruction  -- distinct gene-program matrices W
#   right  : beside each solution, a Kaplan-Meier panel whose risk-group
#            separation is only MODERATE and VARIES across solutions
#
# Point: many factorizations reconstruct X equally well, yet imply different
# (and merely moderate) survival separation -> reconstruction does not pin down
# prognosis. Fully schematic; illustrative curves only, no data.
#
# Run:  Rscript code/14_defense_figC.R

suppressPackageStartupMessages({
  library(ggplot2)
})

OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PAL <- list(
  X       = "#3b3b3b",   # expression matrix (greyscale)
  W       = "#1b9e77",   # gene-program matrix W (green)
  arrow   = "grey45",
  surv_lo = "#2166ac",   # low-risk / good prognosis (blue)
  surv_hi = "#b2182b"    # high-risk / poor prognosis (red)
)

# ── tile helper: render a matrix as identity-colored tiles in a block ────────
ramp_hex <- function(m, hi) {
  pal <- grDevices::colorRampPalette(c("#ffffff", hi))(100)
  idx <- pmax(1, pmin(100, ceiling(m / max(m) * 100)))
  matrix(pal[idx], nrow(m), ncol(m))
}
tiles <- function(m, hi, x0, ytop, cellw, cellh) {
  hex <- ramp_hex(m, hi)
  d <- expand.grid(r = seq_len(nrow(m)), c = seq_len(ncol(m)))
  data.frame(x = x0 + (d$c - 0.5) * cellw,
             y = ytop - (d$r - 0.5) * cellh,
             fill = hex[cbind(d$r, d$c)],
             w = cellw, h = cellh)
}

# ── realistic KM step curve from simulated event times ──────────────────────
stepify <- function(t, S) {
  n  <- length(t)
  xs <- rep(t, each = 2)[-1]
  ys <- rep(S, each = 2)[-(2 * n)]
  data.frame(t = xs, S = ys)
}
# Simulate n exponential event times; censor at tmax; build the Kaplan-Meier
# step (flat at 1 until the first event, irregular drops of 1/n at each event,
# administrative-censoring plateau at the right). Mapped into a panel band.
km_sim <- function(rate, n, x0, xspan, ybase, yspan, grp, seed, tmax = 12) {
  set.seed(seed)
  evt <- sort(rexp(n, rate))
  evt <- evt[evt <= tmax]
  tt  <- c(0, evt, tmax)
  SS  <- c(1, 1 - seq_along(evt) / n, 1 - length(evt) / n)
  st  <- stepify(tt, SS)
  data.frame(x = x0 + (st$t / tmax) * xspan,
             y = ybase + st$S * yspan, grp = grp)
}

# ════════════════════════════════════════════════════════════════════════════
fig_C <- function() {
  set.seed(21)

  # ---- expression matrix X (genes x patients) ----
  Xm   <- matrix(rgamma(14 * 10, 2, 2), 14, 10)
  x_x0 <- 1.0; x_ytop <- 14.5; x_w <- 5.6; x_h <- 11.5
  X_tiles <- tiles(Xm, PAL$X, x_x0, x_ytop, x_w / ncol(Xm), x_h / nrow(Xm))
  x_right <- x_x0 + x_w
  x_mid_y <- x_ytop - x_h / 2                            # fan origin height

  # ---- three branch rows ----
  yc   <- c(13, 8, 3)                                    # band centers
  # illustrative rates (low-risk, high-risk); all MODERATE but varying
  rates <- list(c(lo = 0.055, hi = 0.110),               # moderate
                c(lo = 0.065, hi = 0.092),               # weaker
                c(lo = 0.058, hi = 0.098))               # less pronounced
  seeds   <- list(c(101, 102), c(103, 104), c(105, 106))
  sep_lab <- c("moderate", "weaker", "moderate")
  recon_lab <- c("R² ≈ 0.951", "R² ≈ 0.949", "R² ≈ 0.952")

  w_x0 <- 11.5; w_w <- 2.4; w_h <- 2.6                   # W glyph
  km_x0 <- 18.2; km_xspan <- 8.6; km_h <- 3.0           # KM panel
  n_km <- 45

  W_tiles <- list(); KM <- list()
  for (i in 1:3) {
    Wm <- matrix(rgamma(12 * 3, 1.4 + 0.3 * i, 2), 12, 3)  # distinct programs
    W_tiles[[i]] <- tiles(Wm, PAL$W, w_x0, yc[i] + w_h / 2, w_w / 3, w_h / 12)
    ybase <- yc[i] - km_h / 2
    KM[[i]] <- rbind(
      km_sim(rates[[i]]["lo"], n_km, km_x0, km_xspan, ybase, km_h, "low",  seeds[[i]][1]),
      km_sim(rates[[i]]["hi"], n_km, km_x0, km_xspan, ybase, km_h, "high", seeds[[i]][2])
    )
    KM[[i]]$row <- i
  }
  KMall <- do.call(rbind, KM)

  at  <- function(...) annotate("text", ...)
  seg <- function(...) annotate("segment", ...)
  frame <- function(x0, y0, w, h)
    annotate("rect", xmin = x0, xmax = x0 + w, ymin = y0, ymax = y0 + h,
             fill = NA, color = "grey75", linewidth = 0.5)
  ar <- arrow(length = unit(0.22, "cm"), type = "closed")

  p <- ggplot() +
    # ---- X matrix + labels (above the block, clear of tiles) ----
    geom_tile(data = X_tiles, aes(x, y, fill = fill, width = w, height = h),
              color = "white", linewidth = 0.3) +
    frame(x_x0, x_ytop - x_h, x_w, x_h) +
    at(x = x_x0 + x_w / 2, y = x_ytop + 0.95, label = "Expression matrix",
       fontface = "bold", size = 8) +
    at(x = x_x0 + x_w / 2, y = x_ytop + 0.35, label = "genes × patients",
       size = 6.4, color = "grey30") +
    # ---- fan caption, lifted above the whole fan, clear of arrows ----
    at(x = 8.8, y = 13.6, label = "many solutions,\nall fit X well",
       size = 5.5, fontface = "italic", color = "grey35", lineheight = 0.9)

  # ---- fan of branch arrows: X -> each W glyph ----
  for (i in 1:3)
    p <- p + seg(x = x_right + 0.15, xend = w_x0 - 0.3,
                 y = x_mid_y, yend = yc[i],
                 color = PAL$arrow, linewidth = 1.0, arrow = ar)

  # ---- per-branch W glyphs, reconstruction tags, KM panels ----
  for (i in 1:3) {
    ybase <- yc[i] - km_h / 2
    p <- p +
      geom_tile(data = W_tiles[[i]],
                aes(x, y, fill = fill, width = w, height = h),
                color = "white", linewidth = 0.2) +
      frame(w_x0, yc[i] - w_h / 2, w_w, w_h) +
      at(x = w_x0 + w_w / 2, y = yc[i] + w_h / 2 + 0.5,
         label = paste0("W", i), fontface = "bold", size = 6.7) +
      at(x = w_x0 + w_w / 2, y = yc[i] - w_h / 2 - 0.5,
         label = recon_lab[i], size = 5.4, color = "grey30") +
      # W glyph -> KM panel (stop short of the axis / survival label)
      seg(x = w_x0 + w_w + 0.2, xend = km_x0 - 1.2, y = yc[i], yend = yc[i],
          color = PAL$arrow, linewidth = 0.9, arrow = ar) +
      # KM axes (light L)
      seg(x = km_x0, xend = km_x0, y = ybase, yend = ybase + km_h,
          color = "grey55", linewidth = 0.5) +
      seg(x = km_x0, xend = km_x0 + km_xspan, y = ybase, yend = ybase,
          color = "grey55", linewidth = 0.5) +
      # separation descriptor (right of panel)
      at(x = km_x0 + km_xspan + 0.35, y = yc[i], label = sep_lab[i],
         hjust = 0, size = 5.6, fontface = "italic", color = "grey25")
  }

  # ---- KM curves ----
  p <- p +
    geom_path(data = KMall,
              aes(x, y, group = interaction(row, grp), color = grp),
              linewidth = 1.1) +
    scale_color_manual(values = c(low = PAL$surv_lo, high = PAL$surv_hi),
                       guide = "none")

  # ---- inline legend (top KM), axis titles ----
  topy <- yc[1] + km_h / 2 + 0.45
  p <- p +
    seg(x = km_x0 + 0.3, xend = km_x0 + 1.0, y = topy, yend = topy,
        color = PAL$surv_lo, linewidth = 1.1) +
    at(x = km_x0 + 1.15, y = topy, label = "low risk", hjust = 0,
       size = 5.0, color = PAL$surv_lo) +
    seg(x = km_x0 + 3.3, xend = km_x0 + 4.0, y = topy, yend = topy,
        color = PAL$surv_hi, linewidth = 1.1) +
    at(x = km_x0 + 4.15, y = topy, label = "high risk", hjust = 0,
       size = 5.0, color = PAL$surv_hi) +
    at(x = km_x0 - 0.6, y = yc[2], label = "survival", angle = 90,
       size = 5.4, color = "grey40") +
    at(x = km_x0 + km_xspan / 2, y = yc[3] - km_h / 2 - 0.75, label = "time",
       size = 5.4, color = "grey40")

  p +
    scale_fill_identity() +
    coord_cartesian(xlim = c(-0.5, 32), ylim = c(-0.3, 15.6), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(8, 14, 8, 10))
}

# ── render ───────────────────────────────────────────────────────────────────
fig <- fig_C()
ggsave(file.path(OUT, "new_C_reconstruction_vs_prognosis.pdf"), fig,
       width = 13.5, height = 7.6, device = cairo_pdf)
ggsave(file.path(OUT, "new_C_reconstruction_vs_prognosis.png"), fig,
       width = 13.5, height = 7.6, dpi = 300, bg = "white")
message("Saved new_C_reconstruction_vs_prognosis (13.5x7.6 in)")
