#!/usr/bin/env Rscript
# code/12_defense_new_figures.R — Background/methods figures for the defense deck
#
# Generates the six "new" schematic/conceptual figures described in
# docs/defense_talk_plan.R (Section 4):
#   A  bulk-tumor mixture cartoon
#   B  NMF X = WH decomposition visual
#   C  variance-vs-prognosis misalignment cartoon   (highest value)
#   D  supervision-on-W-not-H diagram               (highest value)
#   E  simulation ground-truth schematic
#   F  CONSORT-style cohort flow
#
# Self-contained: depends only on ggplot2 + cowplot (already used in 09a) and,
# for figure E, R/simulation_functions/simulate_W.R. Does NOT read pipeline
# results in results/. Outputs vector PDF + 300-dpi PNG to figures/defense/.
#
# Run:  Rscript code/12_defense_new_figures.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── Shared palette (reuse across slides for visual consistency) ─────────────
PAL <- list(
  W      = "#1b9e77",  # gene-program matrix W (green)
  H      = "#d95f02",  # sample-loading matrix H (orange)
  X      = "#3b3b3b",  # expression matrix X (greyscale)
  cox    = "#c0392b",  # survival / Cox gradient (red) — used on Slides 6,8,10
  recon  = "#2c3e50",  # reconstruction gradient (slate)
  surv_lo = "#2166ac", # long survival (blue)
  surv_hi = "#b2182b", # short survival (red)
  cell = c(Malignant = "#e41a1c", CAF = "#377eb8",
           Immune = "#4daf4a", Exocrine = "#984ea3")
)

# Save helper: vector PDF + matching 300-dpi PNG.
# Use the Cairo devices so UTF-8 math glyphs (β, ≈, ×, →, −, ²) render correctly;
# the base pdf() device drops them via mbcsToSbcs on Windows.
save_fig <- function(plot, name, width, height) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), plot, width = width, height = height,
         dpi = 300, bg = "white", type = "cairo")
  message("Saved ", name, " (", width, "x", height, " in)")
}

# Melt a matrix to long form without reshape2
melt_mat <- function(m, row = "row", col = "col", val = "value") {
  df <- expand.grid(row = seq_len(nrow(m)), col = seq_len(ncol(m)))
  df$value <- as.vector(m)
  names(df) <- c(row, col, val)
  df
}

# ════════════════════════════════════════════════════════════════════════════
# Figure A — Bulk-tumor mixture cartoon  (Slide 2)
# ════════════════════════════════════════════════════════════════════════════
fig_A <- function() {
  set.seed(11)
  # Genes (rows) with recognizable cell-type marker structure
  genes <- c("KRT17", "S100A2", "TFF1", "COL1A1", "ACTA2",
             "CD8A", "PTPRC", "PRSS1", "CPA1", "ACTB")
  g <- length(genes)
  types <- c("Malignant", "CAF", "Immune", "Exocrine")
  # Mixture fractions: an exocrine-dominant (low-purity) sample
  frac  <- c(Malignant = 0.12, CAF = 0.25, Immune = 0.18, Exocrine = 0.45)

  # Each cell type has its OWN expression profile: high on its marker genes,
  # low elsewhere; ACTB is a shared housekeeping gene high in all.
  prof <- matrix(0.07, g, 4, dimnames = list(genes, types))
  prof[1:3, "Malignant"] <- c(0.95, 0.82, 0.70)
  prof[4:5, "CAF"]       <- c(0.92, 0.78)
  prof[6:7, "Immune"]    <- c(0.88, 0.80)
  prof[8:9, "Exocrine"]  <- c(0.90, 0.83)
  prof["ACTB", ]         <- c(0.60, 0.55, 0.50, 0.58)
  prof[] <- pmin(1, prof + matrix(runif(g * 4, 0, 0.05), g, 4))  # keep dims
  bulk <- as.vector(prof %*% frac)            # observed bulk = weighted mixture
  bulk <- bulk / max(bulk)

  # Draw everything in ONE ggplot (manual coordinates) so there is no plot_grid
  # margin whitespace and the rotated "genes" label can sit exactly where we want.
  ramp_hex <- function(v, hi) {
    pal <- grDevices::colorRampPalette(c("#ffffff", hi))(100)
    pal[pmax(1, pmin(100, ceiling(v * 100)))]
  }
  colw <- 2.6
  col_tiles <- function(v, hi, bx)
    data.frame(x = bx, y = -seq_along(v), fill = ramp_hex(v, hi))
  frame <- function(bx) annotate("rect", xmin = bx - colw / 2, xmax = bx + colw / 2,
    ymin = -(g + 0.5), ymax = -0.5, fill = NA, color = "grey75", linewidth = 0.5)

  name_x  <- 2.2                              # gene names right-aligned here
  genes_x <- 0.15                             # rotated axis title, far left
  bx <- c(mal = 4.1, caf = 8.1, imm = 12.1, exo = 16.1)
  bx_bulk <- 20.7
  ymid <- -(g + 1) / 2

  dat <- rbind(
    col_tiles(prof[, "Malignant"], PAL$cell["Malignant"], bx["mal"]),
    col_tiles(prof[, "CAF"],       PAL$cell["CAF"],       bx["caf"]),
    col_tiles(prof[, "Immune"],    PAL$cell["Immune"],    bx["imm"]),
    col_tiles(prof[, "Exocrine"],  PAL$cell["Exocrine"],  bx["exo"]),
    col_tiles(bulk,                PAL$X,                 bx_bulk)
  )

  at <- function(...) annotate("text", ...)
  ty <- 1.5; sy <- 0.4                        # title / subtitle y
  cx <- (bx["mal"] + bx_bulk) / 2

  ggplot(dat, aes(x, y, fill = fill)) +
    geom_tile(width = colw, height = 1, color = "white", linewidth = 0.5) +
    scale_fill_identity() +
    frame(bx["mal"]) + frame(bx["caf"]) + frame(bx["imm"]) + frame(bx["exo"]) + frame(bx_bulk) +
    # gene names + rotated "genes" axis title (left of the names)
    at(x = name_x, y = -seq_len(g), label = genes, hjust = 1, size = 4.9, color = "grey15") +
    at(x = genes_x, y = ymid, label = "genes", angle = 90, size = 5.6,
       fontface = "italic", color = "grey40") +
    # column headers + mixture fractions
    at(x = bx["mal"], y = ty, label = "Malignant", fontface = "bold", size = 6.2) +
    at(x = bx["mal"], y = sy, label = "12%", size = 5.2, color = "grey30") +
    at(x = bx["caf"], y = ty, label = "CAF", fontface = "bold", size = 6.2) +
    at(x = bx["caf"], y = sy, label = "25%", size = 5.2, color = "grey30") +
    at(x = bx["imm"], y = ty, label = "Immune", fontface = "bold", size = 6.2) +
    at(x = bx["imm"], y = sy, label = "18%", size = 5.2, color = "grey30") +
    at(x = bx["exo"], y = ty, label = "Exocrine", fontface = "bold", size = 6.2) +
    at(x = bx["exo"], y = sy, label = "45%", size = 5.2, color = "grey30") +
    at(x = bx_bulk, y = ty, label = "Bulk", fontface = "bold", size = 6.2) +
    at(x = bx_bulk, y = sy, label = "observed", size = 5.2, color = "grey30") +
    # operators
    at(x = (bx["mal"] + bx["caf"]) / 2, y = ymid, label = "+", fontface = "bold", size = 9, color = "grey45") +
    at(x = (bx["caf"] + bx["imm"]) / 2, y = ymid, label = "+", fontface = "bold", size = 9, color = "grey45") +
    at(x = (bx["imm"] + bx["exo"]) / 2, y = ymid, label = "+", fontface = "bold", size = 9, color = "grey45") +
    at(x = (bx["exo"] + bx_bulk) / 2,   y = ymid, label = "=", fontface = "bold", size = 9, color = "grey45") +
    # title + caption
    at(x = cx, y = 3.1, label = "Bulk expression = weighted sum of cell-type expression profiles",
       fontface = "bold", size = 7) +
    at(x = cx, y = -(g + 1.7),
       label = "each cell type has its own profile over the same genes; the bulk sample observes only their mixture",
       size = 4.7, fontface = "italic", color = "grey30") +
    coord_equal(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(6, 12, 6, 8))
}

# ════════════════════════════════════════════════════════════════════════════
# Figure B — NMF X = WH decomposition  (Slide 3)
# ════════════════════════════════════════════════════════════════════════════
fig_B <- function() {
  set.seed(3)
  p <- 16; n <- 14; k <- 3   # p genes, n samples, k programs
  Wt <- matrix(rgamma(p * k, 2, 2), p, k)
  Ht <- matrix(rgamma(k * n, 2, 2), k, n)
  Xt <- Wt %*% Ht

  # Draw all three matrices in ONE ggplot on a shared coordinate grid (1 unit =
  # 1 cell) so cells stay square (coord_equal) with NO plot_grid whitespace.
  # Each block keeps its own white->hue ramp via scale_fill_identity. All three
  # blocks are vertically CENTERED on the equation midline (ymid) so the operators
  # and H's "k" row label line up with the matrices they connect.
  ramp_hex <- function(m, hi) {
    pal <- grDevices::colorRampPalette(c("#ffffff", hi))(100)
    idx <- pmax(1, pmin(100, ceiling(m / max(m) * 100)))
    matrix(pal[idx], nrow(m), ncol(m))
  }
  # ytop = y of row 0; row r is drawn at ytop - r. Centering a block of nrow rows
  # on the midline means ytop = ymid + (nrow + 1)/2.
  tiles <- function(m, hi, x0, ytop) {
    hex <- ramp_hex(m, hi)
    d <- expand.grid(r = seq_len(nrow(m)), c = seq_len(ncol(m)))
    data.frame(x = x0 + d$c, y = ytop - d$r, fill = hex[cbind(d$r, d$c)])
  }
  frame <- function(x0, ncol_m, nrow_m, ytop)
    annotate("rect", xmin = x0 + 0.5, xmax = x0 + ncol_m + 0.5,
             ymin = ytop - (nrow_m + 0.5), ymax = ytop - 0.5,
             fill = NA, color = "grey75", linewidth = 0.5)

  # Wider W->H gap so the x operator and H's "k" row label both fit without colliding
  GAP_XW <- 2.6; GAP_WH <- 4.4
  x0X <- 0; x0W <- n + GAP_XW; x0H <- x0W + k + GAP_WH
  ymid <- -(p + 1) / 2
  ytopX <- 0                            # X/W span full height (rows -1..-p)
  ytopH <- 0                            # H TOP-ALIGNED with X/W (rows -1..-k);
                                        # integer rows keep all cells square
  yH    <- ytopH - (k + 1) / 2          # vertical center of H (operator/label line)
  dat <- rbind(tiles(Xt, PAL$X, x0X, ytopX),
               tiles(Wt, PAL$W, x0W, ytopX),
               tiles(Ht, PAL$H, x0H, ytopH))

  cxX <- x0X + (n + 1) / 2; cxW <- x0W + (k + 1) / 2; cxH <- x0H + (n + 1) / 2
  sym1   <- (x0X + n + 0.5 + x0W + 0.5) / 2   # operator centered between X and W
  sym2   <- (x0W + k + 0.5) + 1.4             # operator just right of W
  kH_lab <- (x0H + 0.5) - 0.5                 # H row label just left of H
  axfs <- 6.6                                  # row / column label font size
  at <- function(...) annotate("text", ...)

  ggplot(dat, aes(x, y, fill = fill)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_identity() +
    frame(x0X, n, p, ytopX) + frame(x0W, k, p, ytopX) + frame(x0H, n, k, ytopH) +
    # titles + subtitles
    at(x = cxX, y = 2.4, label = "X", fontface = "bold", size = 9) +
    at(x = cxX, y = 0.9, label = "expression", size = 6, color = "grey30") +
    at(x = cxW, y = 2.4, label = "W", fontface = "bold", size = 9) +
    at(x = cxW, y = 0.9, label = "gene programs", size = 6, color = "grey30") +
    at(x = cxH, y = 2.4, label = "H", fontface = "bold", size = 9) +
    at(x = cxH, y = 0.9, label = "loadings", size = 6, color = "grey30") +
    # operators: ~ centered on the tall X/W blocks; x aligned to H's center
    at(x = sym1, y = ymid, label = "≈", fontface = "bold", size = 13, color = "grey40") +
    at(x = sym2, y = yH,   label = "×", fontface = "bold", size = 12, color = "grey40") +
    # row / column labels (enlarged for projection)
    at(x = x0X - 0.5, y = ymid, label = "genes (p)", angle = 90, size = axfs, color = "grey30") +
    at(x = cxX, y = -(p + 1.5), label = "samples (n)", size = axfs, color = "grey30") +
    at(x = cxW, y = -(p + 1.5), label = "k", size = axfs, color = "grey30") +
    at(x = kH_lab - 0.5, y = yH, label = "k", angle = 90, size = axfs, color = "grey30") +
    at(x = cxH, y = ytopH - (k + 1.4), label = "samples (n)", size = axfs, color = "grey30") +
    coord_equal(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(8, 12, 8, 14))
}

# ════════════════════════════════════════════════════════════════════════════
# Figure C — variance-vs-prognosis misalignment  (Slide 4; callback Slide 12)
# ════════════════════════════════════════════════════════════════════════════
fig_C <- function() {
  set.seed(42)
  n <- 170
  v <- rnorm(n, 0, 3.0)    # highest-variance latent coord (purity/exocrine)
  u <- rnorm(n, 0, 0.85)   # prognostic latent coord (low variance)
  risk <- u                 # survival risk driven by the LOW-variance axis

  # rotate so neither latent axis aligns with the plot axes
  th <- 27 * pi / 180
  R <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  XY <- cbind(v, u) %*% t(R)
  dat <- data.frame(X = XY[, 1], Y = XY[, 2], risk = risk)

  # axis arrows (unit directions scaled), anchored at centroid
  c0 <- c(mean(dat$X), mean(dat$Y))
  var_dir  <- as.vector(R %*% c(1, 0))   # high-variance direction
  prog_dir <- as.vector(R %*% c(0, 1))   # prognostic direction
  arr <- function(dir, len) data.frame(
    x = c0[1] - dir[1] * len, y = c0[2] - dir[2] * len,
    xend = c0[1] + dir[1] * len, yend = c0[2] + dir[2] * len)
  a_var  <- arr(var_dir, 6.0)
  a_prog <- arr(prog_dir, 2.1)

  ggplot(dat, aes(X, Y)) +
    geom_point(aes(color = risk), size = 2.6, alpha = 0.9) +
    scale_color_gradient(low = PAL$surv_lo, high = PAL$surv_hi,
                         name = "patient risk", breaks = c(min(risk), max(risk)),
                         labels = c("low\n(long surv.)", "high\n(short surv.)")) +
    # variance axis (what unsupervised NMF chases)
    geom_segment(data = a_var, aes(x, y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.3, "cm"), ends = "both", type = "closed"),
                 linewidth = 1.2, color = "grey25") +
    annotate("text", x = a_var$xend, y = a_var$yend + 0.5,
             label = "Highest-variance axis\n(tumor purity / exocrine)",
             fontface = "bold", size = 4.1, hjust = 0.7, color = "grey20") +
    # prognostic axis (what we actually care about)
    geom_segment(data = a_prog, aes(x, y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.3, "cm"), ends = "both", type = "closed"),
                 linewidth = 1.2, color = PAL$cox) +
    annotate("text", x = a_prog$xend + 1.4, y = a_prog$yend + 0.2,
             label = "True prognostic axis", fontface = "bold",
             size = 4.1, color = PAL$cox) +
    labs(title = "Variance ≠ prognosis",
         subtitle = "Unsupervised NMF projects onto the grey axis and misses the survival gradient") +
    coord_equal() +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 10.5, color = "grey30"),
          axis.title = element_blank(), axis.text = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right", legend.title = element_text(face = "bold"))
}

# ════════════════════════════════════════════════════════════════════════════
# Figure D — supervision enters W, not H  (Slide 7)
# ════════════════════════════════════════════════════════════════════════════
fig_D <- function() {
  box <- function(xc, yc, w, h, fill, label, sub = NULL, lab_size = 6) {
    list(
      geom_rect(data = data.frame(xmin = xc - w / 2, xmax = xc + w / 2,
                                  ymin = yc - h / 2, ymax = yc + h / 2),
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = fill, color = "grey20", linewidth = 0.5, alpha = 0.9,
                inherit.aes = FALSE),
      annotate("text", x = xc, y = yc + ifelse(is.null(sub), 0, 0.18),
               label = label, fontface = "bold", size = lab_size, color = "white",
               lineheight = 0.9),
      if (!is.null(sub)) annotate("text", x = xc, y = yc - 0.28, label = sub,
                                  size = 3.2, color = "white") else NULL
    )
  }
  seg <- function(x, xend, y, yend, color, lwd = 1.1, lty = 1) {
    annotate("segment", x = x, xend = xend, y = y, yend = yend, colour = color,
             linewidth = lwd, linetype = lty,
             arrow = arrow(length = unit(0.26, "cm"), type = "closed"))
  }

  ggplot() +
    # Top row: X ~ W H  (W bottom edge at y = 7.3, H bottom edge at y = 7.5)
    box(2.0, 8, 1.6, 1.4, PAL$X, "X", "expression") +
    annotate("text", x = 3.15, y = 8, label = "≈", size = 11) +
    box(4.3, 8, 1.0, 1.4, PAL$W, "W", "programs") +
    annotate("text", x = 5.2, y = 8, label = "×", size = 9) +
    box(6.5, 8, 1.7, 1.0, PAL$H, "H", "loadings") +
    # Reconstruction loss node (acts on both W and H)
    box(4.7, 5.2, 4.0, 1.1, PAL$recon,
        "Reconstruction loss\n||X − WH||²", lab_size = 4.4) +
    seg(4.3, 4.3, 5.75, 7.30, PAL$recon) +        # recon -> W (vertical)
    seg(5.7, 6.5, 5.75, 7.50, PAL$recon) +        # recon -> H
    annotate("text", x = 8.4, y = 5.2, label = "acts on\nboth W and H",
             size = 3.4, color = PAL$recon, fontface = "italic", lineheight = 0.9) +
    # Survival / Cox node (acts on W only)
    box(2.3, 2.4, 3.2, 1.1, PAL$cox,
        "Cox survival loss\nZ = W'X", lab_size = 4.4) +
    seg(2.9, 3.95, 2.95, 7.30, PAL$cox, lwd = 1.6) +   # cox -> W (bold red)
    annotate("text", x = 2.55, y = 5.2, label = "survival gradient\nacts on W only",
             size = 3.5, color = PAL$cox, fontface = "bold", angle = 70,
             lineheight = 0.9) +
    # crossed-out (suppressed) gradient to H
    seg(3.2, 5.9, 2.95, 7.50, "grey60", lwd = 0.9, lty = "dashed") +
    annotate("text", x = 4.7, y = 4.7, label = "X", size = 7, fontface = "bold",
             color = "grey55") +
    annotate("text", x = 7.3, y = 3.5,
             label = "H keeps its\nmixture-coefficient\ninterpretation",
             size = 3.4, color = "grey35", fontface = "italic", lineheight = 0.9) +
    # portability note
    annotate("text", x = 5.0, y = 0.9,
             label = "Programs W fixed at training → new samples scored by projection   Z_new = (W*)' X_new",
             size = 3.6, fontface = "italic", color = "grey20") +
    coord_equal(xlim = c(0.2, 10), ylim = c(0.3, 9)) +
    theme_void()
}

# ════════════════════════════════════════════════════════════════════════════
# Figure E — simulation ground-truth schematic  (Slide 9)
# ════════════════════════════════════════════════════════════════════════════
fig_E <- function() {
  src <- file.path("R", "simulation_functions", "simulate_W.R")
  set.seed(7)
  if (file.exists(src)) {
    source(src)
    sim <- simulate_W_marker_background(
      G = 240, K = 3, markers_per_factor = 24, B_size = 72,
      normalize_cols = TRUE, seed = 7)
    W <- sim$W
    grp <- rep("Noise genes", nrow(W))
    names(grp) <- rownames(W)
    for (kk in seq_along(sim$marker_sets))
      grp[sim$marker_sets[[kk]]] <- paste0("F", kk, " markers")
    grp[sim$background] <- "Background"
  } else {
    # Fallback construction if the sim helper is unavailable
    G <- 240; K <- 3; mk <- 24; B <- 72
    W <- matrix(rgamma(G * K, 1, 20), G, K)
    grp <- rep("Noise genes", G)
    idx <- 1
    for (kk in 1:K) { rng <- idx:(idx + mk - 1)
      W[rng, kk] <- rgamma(mk, 3, 0.8); grp[rng] <- paste0("F", kk, " markers"); idx <- idx + mk }
    brng <- idx:(idx + B - 1); W[brng, ] <- rgamma(B * K, 2, 1); grp[brng] <- "Background"
    W <- sweep(W, 2, sqrt(colSums(W^2)), "/")
  }

  # order genes by group for a block-structured heatmap
  ord <- order(factor(grp, levels = c("F1 markers", "F2 markers", "F3 markers",
                                       "Background", "Noise genes")))
  Wo <- W[ord, ]; grpo <- grp[ord]
  d <- melt_mat(Wo, row = "gene", col = "factor", val = "loading")

  # group boundaries for side annotation
  rle <- rle(as.character(grpo))
  ends <- cumsum(rle$lengths); starts <- ends - rle$lengths + 1
  segs <- data.frame(grp = rle$values, ymid = (starts + ends) / 2)

  hm <- ggplot(d, aes(factor, gene, fill = loading)) +
    geom_raster() +
    scale_fill_gradient(low = "#ffffff", high = PAL$W, name = "loading") +
    scale_x_continuous(breaks = 1:3, labels = c("F1", "F2", "F3"),
                       position = "top", expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    labs(title = "Ground-truth W (simulation)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_text(face = "bold", size = 12),
          panel.grid = element_blank(),
          legend.position = "left")

  # side labels for gene groups + survival arrow spanning the F1-marker block
  f1 <- segs[segs$grp == "F1 markers", ]
  f1_lo <- starts[rle$values == "F1 markers"]
  f1_hi <- ends[rle$values == "F1 markers"]
  side <- ggplot(segs) +
    geom_text(aes(x = 0.1, y = ymid, label = grp), hjust = 0, size = 3.7) +
    geom_segment(data = f1,
                 aes(x = 1.65, xend = 1.65, y = f1_lo, yend = f1_hi),
                 arrow = arrow(length = unit(0.18, "cm"), ends = "both"),
                 color = PAL$cox, linewidth = 1.1) +
    annotate("text", x = 1.8, y = f1$ymid,
             label = "drives\nsurvival\n(β1 = 2)", color = PAL$cox,
             fontface = "bold", size = 3.4, hjust = 0, lineheight = 0.9) +
    scale_y_reverse(limits = c(nrow(Wo), 1)) +
    coord_cartesian(xlim = c(0, 2.6), clip = "off") +
    theme_void()

  plot_grid(hm, side, nrow = 1, rel_widths = c(1, 0.62))
}

# ════════════════════════════════════════════════════════════════════════════
# Figure F — CONSORT-style cohort flow  (Slide 14 / backup B1)
# ════════════════════════════════════════════════════════════════════════════
fig_F <- function() {
  nodebox <- function(xc, yc, w, h, fill, label, size = 3.8) {
    list(
      geom_rect(data = data.frame(a = xc - w / 2, b = xc + w / 2,
                                  c = yc - h / 2, d = yc + h / 2),
                aes(xmin = a, xmax = b, ymin = c, ymax = d),
                fill = fill, color = "grey25", linewidth = 0.5, inherit.aes = FALSE),
      annotate("text", x = xc, y = yc, label = label, size = size, lineheight = 0.95)
    )
  }
  down <- function(x, y0, y1) annotate("segment", x = x, xend = x, y = y0, yend = y1,
    arrow = arrow(length = unit(0.22, "cm"), type = "closed"), linewidth = 0.8)

  train_fill <- "#dbe9f6"; val_fill <- "#e8f4e1"; excl_fill <- "#f3e3e3"; ana_fill <- "#fbf3d6"

  ggplot() +
    annotate("text", x = 2.5, y = 9.6, label = "Training", fontface = "bold", size = 5) +
    annotate("text", x = 7.5, y = 9.6, label = "External validation", fontface = "bold", size = 5) +
    # training column
    nodebox(2.5, 8.7, 4.2, 0.9, train_fill,
            "TCGA-PAAD (n=181) + CPTAC-3 (n=140)\n321 pooled samples") +
    down(2.5, 8.25, 7.55) +
    nodebox(4.7, 7.3, 2.0, 0.8, excl_fill, "48 excluded\n(missing surv./QC)", size = 3.2) +
    nodebox(2.5, 6.4, 4.2, 0.9, ana_fill, "273 analytic samples\n139 events  ·  EPV = 46") +
    # validation column
    nodebox(7.5, 8.7, 4.4, 0.9, val_fill,
            "5 cohorts: Dijk, Moffitt, PACA-AU\n(array+seq), Puleo · 979 pooled") +
    down(7.5, 8.25, 7.55) +
    nodebox(9.6, 7.3, 2.0, 0.8, excl_fill,
            "363 excluded\n(non-PDAC/QC)", size = 3.2) +
    nodebox(7.5, 6.4, 4.4, 0.9, ana_fill, "616 analytic samples\n414 events") +
    # exclusion connectors
    annotate("segment", x = 2.5, xend = 3.7, y = 7.3, yend = 7.3, linewidth = 0.5, color = "grey45") +
    annotate("segment", x = 7.5, xend = 8.6, y = 7.3, yend = 7.3, linewidth = 0.5, color = "grey45") +
    # per-cohort detail (offset right of the converging arrow)
    annotate("text", x = 7.7, y = 5.2,
             label = "Dijk 90 (81 ev) · Moffitt 123 (83)\nPACA array 63 (38) · PACA seq 52 (31)\nPuleo 288 (181)",
             size = 2.9, color = "grey30", lineheight = 0.95) +
    # arrows converge into the shared downstream box
    down(2.5, 5.95, 4.5) +
    annotate("segment", x = 7.5, xend = 6.0, y = 5.95, yend = 4.5,
             arrow = arrow(length = unit(0.22, "cm"), type = "closed"), linewidth = 0.8) +
    nodebox(4.5, 4.0, 6.8, 0.9, "#eceaf4",
            "1,970 shared genes  ·  within-sample rank transform  ·  n_top = 270") +
    coord_cartesian(xlim = c(-0.2, 10.8), ylim = c(3.3, 10)) +
    theme_void()
}

# ── Render all ──────────────────────────────────────────────────────────────
save_fig(fig_A(), "new_A_bulk_mixture",        width = 12,  height = 7.6)
save_fig(fig_B(), "new_B_nmf_decomposition",   width = 12,  height = 6.5)
save_fig(fig_C(), "new_C_variance_vs_prognosis", width = 8, height = 6)
save_fig(fig_D(), "new_D_supervise_W_not_H",   width = 8.5, height = 6.5)
save_fig(fig_E(), "new_E_sim_ground_truth",    width = 7.5, height = 6)
save_fig(fig_F(), "new_F_cohort_flow",         width = 9,   height = 6)

message("=== Defense figures written to ", OUT, " ===")
