#!/usr/bin/env Rscript
# code/24_reference_program_multiplicity.R
# ---------------------------------------------------------------------------
# R5 reviewer refinement: quantify higher-rank NMF fragmentation using
# "reference-program multiplicity" from the SAME factor-to-reference
# correspondence analysis already used for main-text Fig. 2a,b and SI Fig.
# S6 (fig-nmf-k7-heatmap), at the prespecified r > 0.2 representation
# threshold, plus a companion COVERAGE metric so a lower multiplicity cannot
# be misread as a failure to represent reference programs.
#
# Reused routine (CALLED, not reimplemented):
#   R/figure_plot_helpers.R :: make_gene_overlap_heatmap(..., return_matrix = TRUE)
#   R/get_top_genes.R       :: get_top_genes()   (SI Appendix Section 8
#                              specificity-score top-50-per-factor selection)
# make_gene_overlap_heatmap(return_matrix = TRUE) returns the factor x
# reference-program Spearman correlation matrix (top 50 genes per factor,
# same reference-program universe, same r > 0.2 threshold) used to build
# Fig. 2a,b / Fig. S6 -- this script only summarizes that matrix; it does not
# recompute correlations or redefine the threshold, top-n, or reference
# universe.
#
# Metrics, over the FIXED reference-program universe (all programs in
# data/derv/cmbSubtypes_formatted.RData `top_genes`, after the identical
# renaming/dedup performed inside make_gene_overlap_heatmap()):
#   coverage     = # reference programs with >= 1 factor at r > 0.2
#   edges        = # (factor, reference-program) pairs with r > 0.2
#   multiplicity = edges / coverage (defined only when coverage > 0)
#
# Curves (rank-wise, at fixed ntop = 270 -- the same n_top the production
# figures use, from results/tar_params_best_tcgacptac.rds$ntop):
#   unsupervised: alpha = 0 at each k, from results/cv_grid/cv_grid_fit_list.rds
#   supervised:   at each k, the best alpha by mean CV C-index
#                 (results/cv_grid/cv_grid_best_alpha.csv,
#                  selection_method == "max_cindex"), the same convention
#                 already used for SI Fig. S4 (fig-cindex-by-k) in
#                 R/cv_grid_helpers.R::plot_cindex_by_k().
# These rank-wise points are SINGLE fits at fixed lambda/nu/n_top (not
# 100-restart consensus fits).
#
# Anchors (100-restart consensus fits, overlaid on the curves):
#   DeSurv k=3            results/tar_fit_desurv_tcgacptac.rds
#   Matched-rank Lee NMF  results/fit_std_desurvk_tcgacptac.rds
#   Unsupervised NMF k=7  results/tar_fit_desurv_alpha0_tcgacptac.rds
#                         (confirmed by dim(W)[2] == 7; this is the fit
#                         underlying the "Rank-optimized unsupervised control (k = 7, alpha = 0)"
#                         row and is already used for SI Fig. S6 in
#                         code/09a_figures.R / code/09b_si_figures.R.)
#
# Inputs are read-only cached objects; nothing upstream is recomputed. Run
# from the repository root as:
#   DESURV_RECOMPUTE=FALSE Rscript code/24_reference_program_multiplicity.R
#
# Outputs:
#   results/reference_program_multiplicity.rds
#   figures/fig_reference_program_multiplicity.pdf
#   figures/fig_reference_program_multiplicity.png (300 dpi)
# ---------------------------------------------------------------------------

message("=== Step 24: Reference-program multiplicity (R5) ===")
source("code/00_helpers.R")
suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

source("R/get_top_genes.R")
source("R/theme_nature.R")
source("R/figure_plot_helpers.R")   # make_gene_overlap_heatmap()

REF_THRESHOLD <- 0.2   # prespecified r > 0.2 representation threshold (unchanged)
TOP_N         <- 50    # top-50-genes-per-factor rule (unchanged; enforced inside
                        # make_gene_overlap_heatmap(), not re-specified here)

# ── Reference-program universe (same file/loader as code/09a_figures.R and
#    code/23_stroma_resolution.R) ─────────────────────────────────────────
top_genes_path <- "data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(top_genes_path)) {
  top_genes_path <- "../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData"
}
if (!file.exists(top_genes_path)) {
  stop("Reference gene list not found. Copy data/derv/cmbSubtypes_formatted.RData ",
       "from DeSurv-paper or set the path manually.")
}
load(top_genes_path)  # loads: top_genes, colors, subtypeList, etc.

# ── Cached model fits / grids (read-only; no recomputation) ──────────────
ntop_value <- load_precomputed("tar_params_best_tcgacptac")$ntop   # 270

tar_fit_desurv        <- load_precomputed("tar_fit_desurv_tcgacptac")        # D, k=3
fit_std_desurvk       <- load_precomputed("fit_std_desurvk_tcgacptac")       # matched-rank Lee NMF, k=3
tar_fit_desurv_alpha0 <- load_precomputed("tar_fit_desurv_alpha0_tcgacptac") # unsupervised NMF, BO k=7

stopifnot(
  "results/tar_fit_desurv_alpha0_tcgacptac.rds is not the k=7 unsupervised anchor" =
    ncol(tar_fit_desurv_alpha0$W) == 7
)

cv_fit_list <- readRDS(file.path(CV_GRID_DIR, "cv_grid_fit_list.rds"))
cv_best_alpha <- read.csv(file.path(CV_GRID_DIR, "cv_grid_best_alpha.csv"),
                           stringsAsFactors = FALSE)

# ── Core summary: call the existing correspondence routine, summarize its
#    output matrix into (coverage, edges, multiplicity) ──────────────────
# `fit_obj` needs only a $W field (make_gene_overlap_heatmap() only reads
# fit_desurv$W internally).
summarize_multiplicity <- function(W, ntop = ntop_value) {
  tops <- get_top_genes(W = W, ntop = ntop)$top_genes
  hm   <- make_gene_overlap_heatmap(
    fit_desurv    = list(W = W),
    tops          = tops,
    top_genes_ref = top_genes,
    return_matrix = TRUE
  )
  cor_mat <- hm$cor_mat   # reference-program x factor Spearman correlations

  hits         <- cor_mat > REF_THRESHOLD
  edges        <- sum(hits, na.rm = TRUE)
  covered_row  <- apply(hits, 1, any, na.rm = TRUE)
  coverage     <- sum(covered_row, na.rm = TRUE)
  multiplicity <- if (coverage > 0) edges / coverage else NA_real_

  list(coverage = coverage, edges = edges, multiplicity = multiplicity,
       n_factors = ncol(W), universe_size = nrow(cor_mat), cor_mat = cor_mat)
}

# Reference-program universe size is fixed by top_genes/the renaming+dedup
# inside make_gene_overlap_heatmap(); confirm it is identical across calls
# below rather than assuming it.
universe_sizes <- c()

# ── Rank-wise curves from results/cv_grid (single fits, fixed lambda/nu/ntop) ─
ntop_val_fn <- function(x) if (is.null(x$ntop) || length(x$ntop) == 0) NA_real_ else x$ntop
cv_ntop   <- vapply(cv_fit_list, ntop_val_fn, numeric(1))
cv_k      <- vapply(cv_fit_list, function(x) x$k, numeric(1))
cv_alpha  <- vapply(cv_fit_list, function(x) x$alpha, numeric(1))

k_grid <- sort(unique(cv_k[!is.na(cv_ntop) & cv_ntop == ntop_value]))

per_config_rows <- list()

for (k_i in k_grid) {
  # Unsupervised arm: alpha = 0 at this k, ntop = ntop_value
  idx_un <- which(!is.na(cv_ntop) & cv_ntop == ntop_value & cv_k == k_i & abs(cv_alpha - 0) < 1e-8)
  if (length(idx_un) != 1) {
    stop(sprintf("Expected exactly one unsupervised (alpha=0) fit at k=%d, ntop=%d; found %d",
                  k_i, ntop_value, length(idx_un)))
  }
  m_un <- summarize_multiplicity(cv_fit_list[[idx_un]]$fit$W)
  universe_sizes <- c(universe_sizes, m_un$universe_size)
  per_config_rows[[length(per_config_rows) + 1]] <- data.frame(
    k = k_i, alpha = 0, arm = "unsupervised",
    coverage = m_un$coverage, edges = m_un$edges, multiplicity = m_un$multiplicity
  )

  # Supervised arm: best alpha by mean CV C-index at this k (same convention
  # as SI Fig. S4 / plot_cindex_by_k(): selection_method == "max_cindex").
  best_row <- cv_best_alpha[
    cv_best_alpha$k == k_i &
    cv_best_alpha$selection_method == "max_cindex" &
    ((is.na(cv_best_alpha$ntop) & is.na(ntop_value)) |
     (!is.na(cv_best_alpha$ntop) & cv_best_alpha$ntop == ntop_value)),
  ]
  if (nrow(best_row) != 1) {
    stop(sprintf("Expected exactly one cv_grid_best_alpha row at k=%d, ntop=%d, max_cindex; found %d",
                  k_i, ntop_value, nrow(best_row)))
  }
  best_alpha_k <- best_row$best_alpha[1]

  idx_sup <- which(!is.na(cv_ntop) & cv_ntop == ntop_value & cv_k == k_i &
                     abs(cv_alpha - best_alpha_k) < 1e-6)
  if (length(idx_sup) != 1) {
    stop(sprintf("Expected exactly one supervised (alpha=%.2f) fit at k=%d, ntop=%d; found %d",
                  best_alpha_k, k_i, ntop_value, length(idx_sup)))
  }
  m_sup <- summarize_multiplicity(cv_fit_list[[idx_sup]]$fit$W)
  universe_sizes <- c(universe_sizes, m_sup$universe_size)
  per_config_rows[[length(per_config_rows) + 1]] <- data.frame(
    k = k_i, alpha = best_alpha_k, arm = "supervised",
    coverage = m_sup$coverage, edges = m_sup$edges, multiplicity = m_sup$multiplicity
  )
}

per_config <- do.call(rbind, per_config_rows)
rownames(per_config) <- NULL

# ── Anchors: 100-restart consensus fits ───────────────────────────────────
m_desurv_k3 <- summarize_multiplicity(tar_fit_desurv$W)
m_nmf_k3    <- summarize_multiplicity(fit_std_desurvk$W)
m_nmf_k7    <- summarize_multiplicity(tar_fit_desurv_alpha0$W)
universe_sizes <- c(universe_sizes, m_desurv_k3$universe_size,
                     m_nmf_k3$universe_size, m_nmf_k7$universe_size)

stopifnot("Reference-program universe size is not constant across configurations" =
            length(unique(universe_sizes)) == 1)

anchors <- data.frame(
  anchor_id    = c("desurv_k3", "nmf_k3_matched_rank", "nmf_k7_unsupervised"),
  label        = c("DeSurv (k = 3)", "Matched-rank Lee NMF (k = 3)",
                    "Rank-optimized unsupervised control (k = 7, alpha = 0)"),
  k            = c(ncol(tar_fit_desurv$W), ncol(fit_std_desurvk$W), ncol(tar_fit_desurv_alpha0$W)),
  n_factors    = c(m_desurv_k3$n_factors, m_nmf_k3$n_factors, m_nmf_k7$n_factors),
  coverage     = c(m_desurv_k3$coverage, m_nmf_k3$coverage, m_nmf_k7$coverage),
  edges        = c(m_desurv_k3$edges, m_nmf_k3$edges, m_nmf_k7$edges),
  multiplicity = c(m_desurv_k3$multiplicity, m_nmf_k3$multiplicity, m_nmf_k7$multiplicity),
  fit_source   = c("results/tar_fit_desurv_tcgacptac.rds",
                    "results/fit_std_desurvk_tcgacptac.rds",
                    "results/tar_fit_desurv_alpha0_tcgacptac.rds"),
  stringsAsFactors = FALSE
)

metadata <- list(
  threshold             = REF_THRESHOLD,
  top_n_genes_per_factor = TOP_N,
  reference_universe_size = unique(universe_sizes),
  ntop_for_get_top_genes  = ntop_value,
  alpha_selection_method  = "max_cindex",
  alpha_selection_source  = "results/cv_grid/cv_grid_best_alpha.csv",
  rank_wise_fits_are      = "single fits at fixed lambda/nu/n_top (results/cv_grid/cv_grid_fit_list.rds), NOT 100-restart consensus fits",
  anchor_fits_are         = "100-restart consensus fits",
  helper_reused           = "R/figure_plot_helpers.R::make_gene_overlap_heatmap(return_matrix = TRUE)",
  top_genes_helper_reused = "R/get_top_genes.R::get_top_genes()",
  top_genes_source        = top_genes_path
)

out <- list(per_config = per_config, anchors = anchors, metadata = metadata)
saveRDS(out, file.path(RESULTS_DIR, "reference_program_multiplicity.rds"))
message(sprintf("Saved %s", file.path(RESULTS_DIR, "reference_program_multiplicity.rds")))

# ── Figure: coverage / edges / multiplicity vs. k, both arms, anchors overlaid ─
plot_df <- per_config |>
  tidyr::pivot_longer(cols = c(coverage, edges, multiplicity),
                       names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    metric = factor(metric, levels = c("coverage", "edges", "multiplicity"),
                     labels = c("Coverage", "Edges", "Reference-program multiplicity")),
    arm    = factor(arm, levels = c("unsupervised", "supervised"),
                     labels = c("Unsupervised (alpha = 0)", "Supervised (best alpha by CV C-index)"))
  )

anchor_plot_df <- anchors |>
  dplyr::select(anchor_id, label, k, coverage, edges, multiplicity) |>
  tidyr::pivot_longer(cols = c(coverage, edges, multiplicity),
                       names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    metric = factor(metric, levels = c("coverage", "edges", "multiplicity"),
                     labels = c("Coverage", "Edges", "Reference-program multiplicity"))
  )

fig <- ggplot2::ggplot(plot_df, ggplot2::aes(x = k, y = value)) +
  ggplot2::geom_line(ggplot2::aes(linetype = arm), linewidth = 0.5) +
  ggplot2::geom_point(ggplot2::aes(shape = arm), size = 1.8) +
  ggplot2::geom_point(
    data = anchor_plot_df,
    ggplot2::aes(x = k, y = value),
    inherit.aes = FALSE, size = 3, stroke = 1.1, shape = 4
  ) +
  ggrepel::geom_text_repel(
    data = anchor_plot_df, ggplot2::aes(x = k, y = value, label = label),
    inherit.aes = FALSE, size = 2.4, max.overlaps = Inf,
    # The two k = 3 anchors sit at nearly the same y on the coverage and
    # edges panels, so their labels collided; repel along y only, push them
    # off the marker to the right, and always draw a leader segment.
    box.padding = 0.6, point.padding = 0.5, force = 6, direction = "y",
    nudge_x = 1.3, hjust = 0, min.segment.length = 0, segment.size = 0.3,
    seed = 1
  ) +
  ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  ggplot2::scale_x_continuous(breaks = seq(2, 12, by = 1)) +
  ggplot2::labs(
    x = "Factorization rank (k)", y = NULL, linetype = NULL, shape = NULL,
    title = "Reference-program coverage, edges and multiplicity by rank"
  ) +
  theme_nature(base_size = 9) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(FIGURE_DIR, "fig_reference_program_multiplicity.pdf"),
                 fig, width = 6, height = 8)
ggplot2::ggsave(file.path(FIGURE_DIR, "fig_reference_program_multiplicity.png"),
                 fig, width = 6, height = 8, dpi = 300)
message("Saved figures/fig_reference_program_multiplicity.pdf and .png")

message("=== Step 24 complete ===")
