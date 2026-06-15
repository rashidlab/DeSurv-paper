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
  surv_lo = "#2166ac", # long survival / good prognosis (blue)
  surv_hi = "#b2182b", # short survival / poor prognosis (red)
  surv_neutral = "#bdbdbd", # outcome-neutral (grey, hollow)
  # program identity — RESERVE red/blue for survival; avoid orange; subtypes by label, not color
  cell = c(Malignant = "#7B3FA0",  # tumor (purple)
           CAF       = "#7F5539",  # stroma/CAF (brown)
           Immune    = "#4DAF4A",  # immune (green)
           Exocrine  = "#969696")  # exocrine/purity (grey)
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
# Figure C — the same genes, factored two ways  (Slide 4; callback Slide 12)
#   Same heatmap idiom as Fig A (gene rows, same cell-type ORIGIN colors), but
#   the cell-type columns are MERGED/REGROUPED into each method's 3 factors.
#   Std NMF merges into composition (Tumor / Microenvironment / Exocrine).
#   DeSurv SPLITS the tumor rows (Classical genes -> D1, Basal genes -> D3) and
#   the CAF rows (restCAF -> D1, proCAF -> D2), and has NO exocrine factor
#   (the exocrine rows stay blank on the right). Each cell keeps its Fig A
#   origin color; intensity = loading. Survival bar: blue good / red poor /
#   grey-hollow neutral.
# ════════════════════════════════════════════════════════════════════════════
fig_C <- function() {
  genes  <- c("TFF1", "LGALS4", "KRT17", "S100A2", "COL1A1", "POSTN", "PTPRC", "PRSS1")
  origin <- c("Mal", "Mal", "Mal", "Mal", "CAF", "CAF", "Imm", "Exo")  # cell-type of origin
  g <- length(genes)
  ocol <- c(Mal = unname(PAL$cell["Malignant"]), CAF = unname(PAL$cell["CAF"]),
            Imm = unname(PAL$cell["Immune"]),     Exo = unname(PAL$cell["Exocrine"]))

  # per-method gene x 3-factor loadings (0 = blank/white)
  nmf <- matrix(0, g, 3)
  nmf[1:4, 1] <- c(.90, .82, .86, .80)         # Tumor  = Classical+Basal genes
  nmf[5:7, 2] <- c(.85, .82, .83)              # Microenv = CAF + Immune
  nmf[8,   3] <- .90                            # Exocrine
  des <- matrix(0, g, 3)
  des[c(1, 2, 5), 1] <- c(.90, .82, .80)        # D1 Classical + restCAF
  des[c(6, 7),    2] <- c(.82, .83)             # D2 proCAF + Immune
  des[c(3, 4),    3] <- c(.88, .83)             # D3 Basal       (exocrine -> none)

  ramp1 <- function(v, hi) {
    pal <- grDevices::colorRampPalette(c("#ffffff", hi))(100)
    pal[max(1, min(100, ceiling(v * 100)))]
  }
  cellfill <- function(vals) mapply(function(v, og) ramp1(v, ocol[[og]]), vals, origin)

  # layout mirrors Fig A: contiguous evenly-spaced columns, white tile borders,
  # grey75 frames, header + colored subtitle (prognosis, like A's "%"),
  # rotated "genes" axis title, presentation fonts, landscape. A thin dashed
  # divider separates the two factorizations.
  colw <- 2.2
  genes_x <- 0.15; name_x <- 2.1
  xcols <- c(4, 7, 10, 13, 16, 19)
  nmf_x <- xcols[1:3]; des_x <- xcols[4:6]
  divx <- (xcols[3] + xcols[4]) / 2
  ymid <- -(g + 1) / 2
  tile_df <- function(M, xs) do.call(rbind, lapply(1:3, function(j)
    data.frame(x = xs[j], y = -(1:g), fill = cellfill(M[, j]))))
  dat <- rbind(tile_df(nmf, nmf_x), tile_df(des, des_x))

  at <- function(...) annotate("text", ...)
  frame <- function(cx) annotate("rect", xmin = cx - colw / 2, xmax = cx + colw / 2,
    ymin = -(g + 0.5), ymax = -0.5, fill = NA, color = "grey75", linewidth = 0.5)
  Sgd <- PAL$surv_lo; Spr <- PAL$surv_hi; Sne <- PAL$surv_neutral
  # column header: factor name (bold) + colored prognosis subtitle (like A's %)
  hdr <- function(cx, name, word, wcol) list(
    at(x = cx, y = 1.7, label = name, fontface = "bold", size = 5.4, lineheight = 0.85),
    at(x = cx, y = 0.4, label = word, size = 4.6, color = wcol, fontface = "bold"))
  cx_all <- (genes_x + max(xcols) + colw / 2) / 2

  ggplot(dat, aes(x, y, fill = fill)) +
    geom_tile(width = colw, height = 1, color = "white", linewidth = 0.5) +
    scale_fill_identity() +
    lapply(xcols, frame) +
    # thin dashed divider between the two factorizations
    annotate("segment", x = divx, xend = divx, y = 2.6, yend = -(g + 0.6),
             color = "grey55", linewidth = 0.6, linetype = "22") +
    # gene names + rotated "genes" axis title (Fig A style)
    at(x = name_x, y = -(1:g), label = genes, hjust = 1, size = 4.9, color = "grey15") +
    at(x = genes_x, y = ymid, label = "genes", angle = 90, size = 5.6,
       fontface = "italic", color = "grey40") +
    # method labels over each triplet
    at(x = mean(nmf_x), y = 3.3, label = "Standard NMF", fontface = "bold", size = 6.0) +
    at(x = mean(des_x), y = 3.3, label = "DeSurv",       fontface = "bold", size = 6.0) +
    # column headers + colored prognosis subtitles
    hdr(nmf_x[1], "Tumor", "neutral", "grey45") +
    hdr(nmf_x[2], "Micro-\nenviron.", "neutral", "grey45") +
    hdr(nmf_x[3], "Exocrine", "neutral", "grey45") +
    hdr(des_x[1], "Classical\n+ restCAF", "good", Sgd) +
    hdr(des_x[2], "proCAF\n+ Immune", "poor", Spr) +
    hdr(des_x[3], "Basal", "poor", Spr) +
    # title + caption (Fig A style)
    at(x = cx_all, y = 4.6, label = "Same genes, regrouped into each method's 3 factors",
       fontface = "bold", size = 7) +
    at(x = cx_all, y = -(g + 1.8),
       label = paste("DeSurv splits the tumor rows (Classical → D1, Basal → D3) and the CAF rows (restCAF/proCAF);",
                     "it models no exocrine factor — the exocrine row stays blank on the right.", sep = "\n"),
       size = 4.7, fontface = "italic", color = "grey30", lineheight = 0.95) +
    coord_equal(xlim = c(-0.4, max(xcols) + colw / 2 + 0.3),
                ylim = c(-(g + 2.7), 5.2), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(6, 12, 6, 8))
}

# ════════════════════════════════════════════════════════════════════════════
# Figure D — supervision enters W, not H  (Slide 7)
# ════════════════════════════════════════════════════════════════════════════
fig_D <- function() {
  box <- function(xc, yc, w, h, fill, label, sub = NULL, lab_size = 6,
                  sub_size = 4.4) {
    list(
      geom_rect(data = data.frame(xmin = xc - w / 2, xmax = xc + w / 2,
                                  ymin = yc - h / 2, ymax = yc + h / 2),
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = fill, color = "grey20", linewidth = 0.5, alpha = 0.9,
                inherit.aes = FALSE),
      annotate("text", x = xc, y = yc + ifelse(is.null(sub), 0, 0.24),
               label = label, fontface = "bold", size = lab_size, color = "white",
               lineheight = 0.9),
      if (!is.null(sub)) annotate("text", x = xc, y = yc - 0.36, label = sub,
                                  size = sub_size, color = "white") else NULL
    )
  }
  seg <- function(x, xend, y, yend, color, lwd = 1.1, lty = 1) {
    annotate("segment", x = x, xend = xend, y = y, yend = yend, colour = color,
             linewidth = lwd, linetype = lty,
             arrow = arrow(length = unit(0.26, "cm"), type = "closed"))
  }

  ggplot() +
    # Top row: X ~ W H  (boxes widened so the larger subtitles fit)
    # Top row: all three boxes top-aligned at y = 8.9; X and W made taller (h=1.8),
    # H shorter (h=1.1) so its bottom hangs higher. × sits at H's vertical middle.
    box(2.0, 8.0,  2.0, 1.8, PAL$X, "X", "expression", lab_size = 9, sub_size = 5) +
    annotate("text", x = 3.25, y = 8.0,  label = "≈", size = 12) +
    box(4.4, 8.0,  1.4, 1.8, PAL$W, "W", "programs", lab_size = 9, sub_size = 5) +
    annotate("text", x = 5.35, y = 8.35, label = "×", size = 10) +
    box(6.7, 8.35, 2.0, 1.1, PAL$H, "H", "loadings", lab_size = 9, sub_size = 5) +
    # Reconstruction loss node (right) — in line with the Cox node
    box(6.6, 5.2, 2.8, 1.1, PAL$recon,
        "Reconstruction loss\n||X − WH||²", lab_size = 5.6) +
    seg(5.7, 4.7, 5.75, 7.10, PAL$recon) +        # recon -> W
    seg(6.9, 6.85, 5.75, 7.80, PAL$recon) +       # recon -> H
    annotate("text", x = 6.6, y = 4.15, label = "acts on\nboth W and H",
             size = 4.9, color = PAL$recon, fontface = "bold", lineheight = 0.9) +
    # Survival / Cox node (left, acts on W only) — same level as reconstruction box
    box(2.4, 5.2, 3.0, 1.1, PAL$cox,
        "Cox survival loss\nZ = W'X", lab_size = 5.6) +
    seg(3.4, 4.15, 5.75, 7.10, PAL$cox, lwd = 1.7) +   # cox -> W (bold red)
    annotate("text", x = 2.4, y = 4.15, label = "survival gradient\nacts on W only",
             size = 4.9, color = PAL$cox, fontface = "bold", lineheight = 0.9) +
    coord_equal(xlim = c(0.7, 8.2), ylim = c(3.6, 9.0), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# ════════════════════════════════════════════════════════════════════════════
# Figure E — simulation ground-truth schematic  (Slide 9)
# ════════════════════════════════════════════════════════════════════════════
fig_E <- function() {
  src <- file.path("R", "simulation_functions", "simulate_W.R")
  set.seed(7)
  # G chosen so the noise block is small (noise = G - K*markers - B_size = 36),
  # keeping the figure from being mostly empty white rows at the bottom.
  if (file.exists(src)) {
    source(src)
    sim <- simulate_W_marker_background(
      G = 180, K = 3, markers_per_factor = 24, B_size = 72,
      normalize_cols = TRUE, seed = 7)
    W <- sim$W
    grp <- rep("Noise genes", nrow(W))
    names(grp) <- rownames(W)
    for (kk in seq_along(sim$marker_sets))
      grp[sim$marker_sets[[kk]]] <- paste0("F", kk, " markers")
    grp[sim$background] <- "Background"
  } else {
    # Fallback construction if the sim helper is unavailable
    G <- 180; K <- 3; mk <- 24; B <- 72
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
    # sqrt-scaled fill lifts the very-low-loading noise genes into a visible tint;
    # legend shown qualitatively (low -> high) because a sqrt axis would otherwise
    # place numeric ticks at uneven positions, and exact loadings aren't the point
    scale_fill_gradient(low = "#ffffff", high = PAL$W, name = "Loading",
                        trans = "sqrt", breaks = range(d$loading),
                        labels = c("low", "high")) +
    scale_x_continuous(breaks = 1:3, labels = c("F1", "F2", "F3"),
                       position = "top", expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    labs(title = "Ground-truth W (simulation)") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_text(face = "bold", size = 16),
          panel.grid = element_blank(),
          legend.position = "left",
          legend.title = element_text(size = 14, margin = margin(b = 12)),
          legend.text = element_text(size = 12))

  # side panel: group labels flush against the heatmap; the F1-marker block gets a
  # survival span-arrow hugging the heatmap edge + a red "drives survival" label.
  f1_lo  <- starts[rle$values == "F1 markers"]
  f1_hi  <- ends[rle$values == "F1 markers"]
  f1_mid <- segs$ymid[segs$grp == "F1 markers"]
  non_f1 <- segs[segs$grp != "F1 markers", ]

  side <- ggplot() +
    # double-headed survival arrow spanning the F1-marker rows, right at the edge
    annotate("segment", x = 0.12, xend = 0.12, y = f1_lo - 0.5, yend = f1_hi + 0.5,
             arrow = arrow(length = unit(0.2, "cm"), ends = "both", type = "closed"),
             color = PAL$cox, linewidth = 1.5) +
    annotate("text", x = 0.3, y = f1_mid,
             label = "F1 markers\n(drive survival, β1 = 2)",
             color = PAL$cox, fontface = "bold", size = 4.6, hjust = 0,
             lineheight = 0.9) +
    geom_text(data = non_f1, aes(x = 0.3, y = ymid, label = grp),
              hjust = 0, size = 4.6, color = "grey15") +
    scale_y_reverse(limits = c(nrow(Wo) + 0.5, 0.5), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 2.7), clip = "off") +
    theme_void()

  plot_grid(hm, side, nrow = 1, rel_widths = c(1, 0.62), align = "h", axis = "tb")
}

# ════════════════════════════════════════════════════════════════════════════
# Figure F — CONSORT-style cohort flow  (Slide 14 / backup B1)
# ════════════════════════════════════════════════════════════════════════════
fig_F <- function() {
  nodebox <- function(xc, yc, w, h, fill, label, size = 5) {
    list(
      geom_rect(data = data.frame(a = xc - w / 2, b = xc + w / 2,
                                  c = yc - h / 2, d = yc + h / 2),
                aes(xmin = a, xmax = b, ymin = c, ymax = d),
                fill = fill, color = "grey25", linewidth = 0.5, inherit.aes = FALSE),
      annotate("text", x = xc, y = yc, label = label, size = size, lineheight = 0.95)
    )
  }
  # orthogonal arrow with a head that lands on the target box edge
  arr <- function(x0, y0, x1, y1, lwd = 0.9)
    annotate("segment", x = x0, xend = x1, y = y0, yend = y1,
             arrow = arrow(length = unit(0.26, "cm"), type = "closed"), linewidth = lwd)
  # plain (headless) connector segment
  ln <- function(x0, y0, x1, y1, lwd = 0.9)
    annotate("segment", x = x0, xend = x1, y = y0, yend = y1, linewidth = lwd, colour = "grey25")

  train_fill <- "#dbe9f6"; val_fill <- "#e8f4e1"; excl_fill <- "#f3e3e3"; ana_fill <- "#fbf3d6"
  xT <- 2.8; xV <- 8.5; xS <- (xT + xV) / 2   # training / validation / shared columns

  ggplot() +
    annotate("text", x = xT, y = 10.2, label = "Training", fontface = "bold", size = 6.5) +
    annotate("text", x = xV, y = 10.2, label = "External validation", fontface = "bold", size = 6.5) +
    # top pooled boxes (bottom edge at 8.7)
    nodebox(xT, 9.2, 5.2, 1.0, train_fill,
            "TCGA-PAAD (n=181) + CPTAC (n=140)\n321 pooled samples") +
    nodebox(xV, 9.2, 5.2, 1.0, val_fill,
            "5 cohorts: Dijk, Moffitt, PACA-AU\n(array+seq), Puleo · 979 pooled") +
    # analytic boxes (top edge 7.1, bottom edge 6.1)
    nodebox(xT, 6.6, 5.2, 1.0, ana_fill, "273 analytic samples\n139 events") +
    nodebox(xV, 6.6, 5.2, 1.0, ana_fill, "616 analytic samples\n414 events") +
    # exclusion boxes (to the right of each column)
    nodebox(5.6, 7.95, 2.6, 0.9, excl_fill, "48 excluded\n(missing surv./QC)", size = 4.3) +
    nodebox(11.3, 7.95, 2.6, 0.9, excl_fill, "363 excluded\n(non-PDAC/QC)", size = 4.3) +
    # pooled -> analytic: vertical arrows, heads touching analytic box tops (7.1)
    arr(xT, 8.7, xT, 7.1) +
    arr(xV, 8.7, xV, 7.1) +
    # exclusion branches: horizontal arrows, heads touching exclusion box left edges
    arr(xT, 7.95, 4.3, 7.95) +
    arr(xV, 7.95, 10.0, 7.95) +
    # merge: down lines from analytic bottoms (6.1) to a horizontal bus (5.2),
    # then a single vertical arrow into the shared box top (4.6) — no diagonals
    ln(xT, 6.1, xT, 5.2) +
    ln(xV, 6.1, xV, 5.2) +
    ln(xT, 5.2, xV, 5.2) +
    arr(xS, 5.2, xS, 4.6) +
    # shared downstream box (top edge 4.6)
    nodebox(xS, 4.1, 7.8, 1.0, "#eceaf4",
            "1,970 shared genes  ·  within-sample rank transform",
            size = 4.7) +
    # per-cohort detail, left-aligned to the right of the validation merge line,
    # anchored just below the analytic box so it stays clear of the shared box
    annotate("text", x = 8.9, y = 6.0, hjust = 0, vjust = 1, size = 4.6,
             color = "grey30", lineheight = 1.05,
             label = paste("Dijk 90 (81 events)", "Moffitt 123 (83)",
                           "PACA-AU array 63 (38)", "PACA-AU seq 52 (31)",
                           "Puleo 288 (181)", sep = "\n")) +
    coord_cartesian(xlim = c(0, 12.9), ylim = c(3.3, 10.6)) +
    theme_void()
}

# ── Render all ──────────────────────────────────────────────────────────────
save_fig(fig_A(), "new_A_bulk_mixture",        width = 12,  height = 7.6)
save_fig(fig_B(), "new_B_nmf_decomposition",   width = 12,  height = 6.5)
# NOTE: Figure C is now generated by code/14_defense_figC.R (matrix -> multiple
# good-reconstruction solutions -> varying-moderate KM curves). Do not re-render
# it here, or it will clobber that version.
save_fig(fig_D(), "new_D_supervise_W_not_H",   width = 8.5, height = 6.1)
save_fig(fig_E(), "new_E_sim_ground_truth",    width = 7.5, height = 6)
save_fig(fig_F(), "new_F_cohort_flow",         width = 11,  height = 7.5)

message("=== Defense figures written to ", OUT, " ===")
