#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Builds the rank-biserial enrichment source table for Supplementary Table
# (tab:rb-enrichment). For every reference program in the fixed 1,970-gene
# universe, reports the rank-biserial statistic r_RB = 2*AUC - 1 for each
# DeSurv factor (D1-D3) and each matched-rank standard NMF factor (N1-N3),
# together with a size-preserving gene-label permutation P value and a joint
# Benjamini-Hochberg adjusted P value computed across ALL factor x signature
# tests (both methods pooled). Saved to results/rb_enrichment_table.rds and
# a human-readable CSV; the SI table reads the cached object.
# ---------------------------------------------------------------------------
suppressMessages({ library(ggplot2) })
source("R/get_top_genes.R")
source("R/figure_plot_helpers.R")

B    <- 20000L   # permutation replicates
SEED <- 1L

fd <- readRDS("results/tar_fit_desurv_tcgacptac.rds")
fn <- readRDS("results/fit_std_desurvk_tcgacptac.rds")
pp <- readRDS("results/tar_params_best_tcgacptac.rds")

tp <- "data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(tp)) tp <- "../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(tp)) stop("Reference program file cmbSubtypes_formatted.RData not found")
load(tp)   # -> top_genes

td <- get_top_genes(W = fd$W, ntop = pp$ntop)$top_genes
tn <- get_top_genes(W = fn$W, ntop = pp$ntop)$top_genes

d_lab <- c("D1 Classical/restCAF", "D2 proCAF", "D3 Basal-like")
n_lab <- c("N1 Tumor", "N2 Exocrine", "N3 Microenviron.")

rd <- make_gene_overlap_heatmap(fd, td, top_genes, factor_labels = d_lab,
                                perm_B = B, perm_seed = SEED)
rn <- make_gene_overlap_heatmap(fn, tn, top_genes, factor_labels = n_lab,
                                perm_B = B, perm_seed = SEED)

# Format the raw "GROUP_Subtype" keys the same way the figure does.
fmt_label <- function(x) {
  idx <- regexpr("_", x)
  out <- ifelse(idx == -1L, x, {
    g <- substr(x, 1, idx - 1L)
    s <- substr(x, idx + 1L, nchar(x))
    s <- gsub("_", " ", s)
    s <- gsub("([a-z])([A-Z][a-z])", "\\1 \\2", s)
    paste0(g, ": ", s)
  })
  ov <- c("DECODER: Classical Tumor" = "DECODER: Classical tumor",
          "DECODER: Basal Tumor"     = "DECODER: Basal-like tumor",
          "PurIST: Basal Like"       = "PurIST: Basal-like",
          "Puleo: Pure Basal-like"   = "Puleo Basal-like",
          "Puleo: tumor Basal-like"  = "Puleo: Basal-like",
          "Puleo: tumor Classical"   = "Puleo: Immune Classical")
  h <- match(out, names(ov)); out[!is.na(h)] <- ov[h[!is.na(h)]]
  out
}

to_long <- function(rb, pv, method, labs) {
  sigs <- rownames(rb)
  do.call(rbind, lapply(seq_len(ncol(rb)), function(j) {
    data.frame(method = method, factor = labs[j],
               signature_raw = sigs, signature = fmt_label(sigs),
               r_RB = rb[, j], p_perm = pv[, j],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

tab <- rbind(to_long(rd$full_rb, rd$full_p, "DeSurv", d_lab),
             to_long(rn$full_rb, rn$full_p, "NMF",    n_lab))
tab <- tab[!is.na(tab$r_RB), ]
# Joint Benjamini-Hochberg across all factor x signature tests, both methods.
tab$q_BH <- p.adjust(tab$p_perm, method = "BH")

# Flag the rows shown in the main figure (the prespecified common panel).
# Shared source of truth with the figure builders.
source("R/fig2_display_panel.R")
tab$in_main_figure <- tab$signature_raw %in% fig2_display_sigs

attr(tab, "perm_B") <- B
attr(tab, "perm_seed") <- SEED
saveRDS(tab, "results/rb_enrichment_table.rds")
write.csv(tab, "results/rb_enrichment_table.csv", row.names = FALSE)
cat(sprintf("Wrote rb_enrichment_table: %d rows (%d signatures x %d factors x 2 methods), B=%d perms.\n",
            nrow(tab), length(unique(tab$signature_raw)), 3L, B))
cat(sprintf("  significant at q<0.05 (joint BH): %d / %d tests\n", sum(tab$q_BH < 0.05), nrow(tab)))
