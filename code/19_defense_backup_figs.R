#!/usr/bin/env Rscript
# code/19_defense_backup_figs.R — Large-font rebuilds of dense backup figures
#
# The SI-PDF rasterizations used on backup slides 27/28/30/31/32 had small
# baked-in fonts and (for the tall ones) overflowed the slide. This rebuilds
# them from source with enlarged base fonts and slide-fitting aspect ratios:
#   bk_nmf_diag      (B8)  NMF rank-selection diagnostics, 3-in-a-row
#   bk_null_mixed    (B7)  null + mixed simulation panels, 2x3 landscape
#   bk_nmf_k7        (B10) NMF (alpha=0) k=7 gene-overlap heatmap
#   bk_subtype_overlap (B12) risk-group subtype composition, DeSurv vs NMF
#   bk_cutpoint_curve (B11) log-rank cutpoint-selection curve
#
# Run:  Rscript code/19_defense_backup_figs.R

message("=== Defense backup figure rebuilds (large fonts) ===")
suppressPackageStartupMessages({
  library(ggplot2); library(cowplot); library(survival); library(survminer)
})
source("code/00_helpers.R")
source("R/get_top_genes.R")
source("R/figure_plot_helpers.R")
source("R/cv_grid_helpers.R")

OUT <- file.path("figures", "defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
big <- theme(text = element_text(size = 19))

# ── B8: NMF rank-selection diagnostics (3-in-a-row, enlarged) ───────────────
fig_res  <- load_precomputed("fig_residuals_tcgacptac")
fig_coph <- load_precomputed("fig_cophenetic_tcgacptac")
fig_sil  <- load_precomputed("fig_silhouette_tcgacptac")
strip <- function(p) p + theme_classic(base_size = 20) +
  theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold"))
leg <- get_legend(fig_res + theme_classic(base_size = 20) +
                    theme(legend.position = "bottom", legend.text = element_text(size = 18),
                          legend.title = element_text(size = 18)))
diag_row <- plot_grid(strip(fig_res), strip(fig_coph), strip(fig_sil),
                      ncol = 3, labels = c("A", "B", "C"), label_size = 24)
ggsave(file.path(OUT, "bk_nmf_diag.png"),
       plot_grid(diag_row, leg, ncol = 1, rel_heights = c(1, 0.12)),
       width = 12, height = 4.7, dpi = 150, bg = "white")  # 1800 px wide
message("Saved bk_nmf_diag.png")

# ── B7: Null + mixed simulation panels (2x3 landscape, enlarged) ────────────
sim_figs     <- load_precomputed("sim_figs_by_scenario")
scenario_ids <- sapply(sim_figs, function(x) x$scenario_id)
analysis_ids <- sapply(sim_figs, function(x) x$analysis_id)
nullp  <- sim_figs[[which(scenario_ids == "R00_null" & analysis_ids == "bo_tune_ntop")]]
mixedp <- sim_figs[[which(scenario_ids == "R_mixed"  & analysis_ids == "bo_tune_ntop")]]
# All panels drop their own legend; one shared Method legend sits below the grid.
ftheme <- function(p) if (inherits(p, "ggplot"))
  p + theme_classic(base_size = 18) +
    theme(plot.title = element_text(size = 18), legend.position = "none") else p
nz <- function(x, y) if (!is.null(x)) x else y
mixed_prec <- nz(mixedp$precision_breakdown, mixedp$precision_box)
panels <- list(
  ftheme(nullp$cindex_box)  + labs(title = "Null: C-index") + scale_y_continuous(limits = c(0.3, 1)),
  ftheme(nullp$k_hist)      + labs(title = "Null: selected k"),
  ftheme(mixedp$cindex_box) + labs(title = "Mixed: C-index") + scale_y_continuous(limits = c(0.3, 1)),
  ftheme(mixed_prec)        + labs(title = "Mixed: precision"),
  ftheme(nz(mixedp$matched_beta_box, mixedp$k_hist)) + labs(title = "Mixed: matched |β|"),
  ftheme(mixedp$k_hist)     + labs(title = "Mixed: selected k")
)
shared_leg <- get_legend(
  nullp$cindex_box + theme_classic(base_size = 18) +
    theme(legend.position = "bottom", legend.title = element_text(size = 18),
          legend.text = element_text(size = 18)))
grid6 <- plot_grid(plotlist = panels, ncol = 3, labels = LETTERS[1:6], label_size = 22)
# Taller 2x3 + one shared legend; intrinsic width < 1920 px (10 in x 180 dpi).
ggsave(file.path(OUT, "bk_null_mixed.png"),
       plot_grid(grid6, shared_leg, ncol = 1, rel_heights = c(1, 0.07)),
       width = 10, height = 6.8, dpi = 180, bg = "white")
message("Saved bk_null_mixed.png")

# ── B10: NMF (alpha=0) k=7 gene-overlap heatmap, enlarged fonts ─────────────
tg_path <- "data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(tg_path)) tg_path <- "../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData"
load(tg_path)  # top_genes
fit_a0  <- load_precomputed("tar_fit_desurv_alpha0_tcgacptac")
ntop_v  <- load_precomputed("tar_params_best_tcgacptac")$ntop
tops_a0 <- get_top_genes(W = fit_a0$W, ntop = ntop_v)
hm <- make_gene_overlap_heatmap(fit_a0, tops_a0$top_genes, top_genes,
                                title = "NMF (α = 0), k = 7",
                                fontsize_row = 12, fontsize = 16, legend_fontsize = 12)
ggsave(file.path(OUT, "bk_nmf_k7.png"),
       plot_grid(hm$plot + theme(plot.margin = margin(8, 2, 2, 2)),
                 ggdraw(hm$legend), ncol = 2, rel_widths = c(4, 0.6)),
       width = 7.5, height = 6.6, dpi = 240, bg = "white")  # 1800 px wide
message("Saved bk_nmf_k7.png")

# ── k=5 (elbow) gene-overlap heatmaps, NMF and DeSurv, enlarged fonts ────────
save_k5 <- function(fit, title, stem) {
  tops <- get_top_genes(W = fit$W, ntop = ntop_v)
  h <- make_gene_overlap_heatmap(fit, tops$top_genes, top_genes, title = title,
                                 fontsize_row = 12, fontsize = 16, legend_fontsize = 12)
  ggsave(file.path(OUT, paste0(stem, ".png")),
         plot_grid(h$plot + theme(plot.margin = margin(8, 2, 2, 2)),
                   ggdraw(h$legend), ncol = 2, rel_widths = c(4, 0.6)),
         width = 7.5, height = 6.6, dpi = 240, bg = "white")  # 1800 px wide
  message("Saved ", stem, ".png")
}
save_k5(load_precomputed("fit_std_elbowk_tcgacptac"),     "NMF, k = 5",    "bk_nmf_k5")
save_k5(load_precomputed("tar_fit_desurv_elbowk_tcgacptac"), "DeSurv, k = 5", "bk_desurv_k5")

# ── Shared setup for subtype + cutpoint (merged validation, lp stats) ───────
data_val <- load_precomputed("data_val_filtered_tcgacptac")
val_named <- data_val
names(val_named) <- vapply(val_named, function(x)
  if (!is.null(x$dataname) && nzchar(x$dataname)) x$dataname else "unknown", character(1))
val_surv <- merge_paca_au_datasets(val_named)
tar_fit   <- load_precomputed("tar_fit_desurv_tcgacptac")
std_fit   <- load_precomputed("fit_std_desurvk_tcgacptac")
lp_d      <- load_precomputed("desurv_lp_stats_tcgacptac")
lp_n      <- load_precomputed("std_desurvk_lp_stats_tcgacptac")
best      <- load_precomputed("tar_params_best_tcgacptac")

# ── B12: Risk-group subtype composition (DeSurv vs NMF), enlarged ───────────
# Local LP/risk-group (the DeSurv package isn't installed here; get_top_genes
# reproduces its contrast-based top-gene selection, as used in code/15).
compute_lp_local <- function(W, beta, X, ntop = NULL) {
  if (is.null(ntop)) return(drop((t(X) %*% W) %*% beta))
  genes <- unique(unlist(get_top_genes(W, as.integer(ntop))$top_genes, use.names = FALSE))
  idx <- match(genes, rownames(W)); idx <- idx[!is.na(idx)]
  if (!length(idx)) return(drop((t(X) %*% W) %*% beta))
  theta_sub <- (W %*% beta)[idx, , drop = FALSE]
  tn <- sqrt(sum(theta_sub^2))
  if (!is.finite(tn)) theta_sub[] <- 0 else if (tn > 0) theta_sub <- theta_sub / tn
  drop(t(X[idx, , drop = FALSE]) %*% theta_sub)
}
risk_group_local <- function(W, beta, X, ntop, lp_mean, lp_sd, z_cut) {
  z <- (compute_lp_local(W, beta, X, ntop) - lp_mean) / lp_sd
  factor(ifelse(z > z_cut, "High", "Low"), levels = c("Low", "High"))
}
collect_so <- function(W, beta, ntop, lp, dv) {
  rows <- list()
  for (nm in names(dv)) {
    vd <- dv[[nm]]; cg <- intersect(rownames(W), rownames(vd$ex))
    if (length(cg) < 2) next
    tt <- vd$sampInfo$time; ev <- vd$sampInfo$event
    idx <- which(is.finite(tt) & !is.na(ev) & tt > 0)
    if (length(idx) < 2) next
    g <- risk_group_local(W[cg, , drop = FALSE], beta, vd$ex[cg, idx, drop = FALSE],
                          ntop, lp$lp_mean, lp$lp_sd, lp$optimal_z_cutpoint)
    si <- vd$sampInfo[idx, ]
    if (!all(c("PurIST", "DeCAF") %in% names(si))) next
    d <- data.frame(group = g, PurIST = si$PurIST, DeCAF = si$DeCAF, stringsAsFactors = FALSE)
    d <- d[complete.cases(d), ]; if (nrow(d)) rows[[nm]] <- d
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}
bar <- function(df, col, colors, lab) {
  df[[col]] <- as.character(df[[col]]); df[[col]][df[[col]] == "permCAF"] <- "proCAF"
  tbl <- table(df$group, df[[col]]); fp <- tryCatch(fisher.test(tbl)$p.value, error = function(e) NA)
  pd <- as.data.frame(tbl); names(pd) <- c("group", "subtype", "count")
  tot <- tapply(pd$count, pd$group, sum); pd$prop <- pd$count / tot[pd$group]
  plab <- if (is.na(fp)) "Fisher p = NA" else if (fp < 0.001) sprintf("Fisher p = %.1e", fp) else sprintf("Fisher p = %.3f", fp)
  ggplot(pd, aes(group, prop, fill = subtype)) + geom_col(width = 0.7) +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5),
              size = 6, color = "white", fontface = "bold") +
    scale_fill_manual(values = colors, name = lab) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Risk group", y = "Proportion", subtitle = plab) +
    theme_bw(base_size = 19) + theme(legend.position = "bottom")
}
pc <- c(`Basal-like` = "orange", Classical = "blue"); dc <- c(proCAF = "#FF4DA6", restCAF = "#1B9E9E")
so_d <- collect_so(tar_fit$W, tar_fit$beta, best$ntop, lp_d, val_surv)
so_n <- collect_so(std_fit$W, std_fit$beta, NULL,      lp_n, val_surv)
row_lab <- function(txt) ggdraw() + draw_label(txt, fontface = "bold", size = 19, angle = 90)
desurv_row <- plot_grid(row_lab("DeSurv"),
                        bar(so_d, "PurIST", pc, "PurIST"), bar(so_d, "DeCAF", dc, "DeCAF"),
                        ncol = 3, rel_widths = c(0.08, 1, 1))
nmf_row <- plot_grid(row_lab("NMF k=3"),
                     bar(so_n, "PurIST", pc, "PurIST"), bar(so_n, "DeCAF", dc, "DeCAF"),
                     ncol = 3, rel_widths = c(0.08, 1, 1))
ggsave(file.path(OUT, "bk_subtype_overlap.png"),
       plot_grid(desurv_row, nmf_row, ncol = 1),
       width = 11, height = 9.6, dpi = 165, bg = "white")  # 1815 px wide, taller still
message("Saved bk_subtype_overlap.png")

# ── B11: Log-rank cutpoint-selection curve, enlarged ────────────────────────
cut_sum <- load_precomputed("desurv_cutpoint_summary_tcgacptac")
p_cut <- plot_cutpoint_curve_logrank(cut_sum, k = best$k, alpha = best$alpha,
                                     ntop = best$ntop, optimal_z = lp_d$optimal_z_cutpoint) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(size = 18), axis.title = element_text(size = 20),
        axis.text = element_text(size = 16, color = "black"))
ggsave(file.path(OUT, "bk_cutpoint_curve.png"), p_cut,
       width = 8, height = 5.6, dpi = 300, bg = "white")
message("Saved bk_cutpoint_curve.png")
message("=== Defense backup figure rebuilds complete ===")
