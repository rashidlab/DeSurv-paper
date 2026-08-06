#!/usr/bin/env Rscript
# 20_program_stability.R
# ---------------------------------------------------------------------------
# Empirical PROGRAM stability of the DeSurv basis (Nature Cancer revision).
#
# Motivation. The manuscript invokes the non-identifiability of reconstruction-
# based factorization against unsupervised methods in three separate passages
# (Introduction, and Discussion twice), and the word "reproducible" appears in
# the title. The only stability evidence currently offered for DeSurv itself is
# Supplementary Fig. S1, which reads each fit's objective trace (`f$lossit`) and
# never inspects `f$W`. Monotone objective descent is numerical convergence, not
# basis stability: a non-identifiable objective can descend monotonically to many
# different optima. This script measures the missing quantity.
#
# Three parts, all computed from cached artifacts. Nothing is refitted.
#   (A) restart vs the reported consensus basis
#   (B) restart vs restart, over independent pairs. Fairer than (A), because the
#       consensus fit aggregates all 100 runs and so is not itself a draw from
#       the restart distribution.
#   (C) consensus reproducibility. Split the 100 cached restarts into two
#       DISJOINT blocks of 50, build the consensus gene partition from each
#       independently per the SI procedure, and compare by Adjusted Rand Index.
#
# IMPORTANT SCOPE LIMIT. Consensus yields an INITIALIZATION W0, after which the
# model is fitted. Part (C) therefore measures the stability of the starting
# point, not of the fitted endpoint. The decisive experiment -- fit the pipeline
# from two independent 100-restart consensus initializations and compare the
# final bases -- requires new fitting and is NOT covered here. See
# code/21_consensus_endpoint.R.
#
# In-repo only. No new data.
# ---------------------------------------------------------------------------

NTOP     <- 270L   # BO-selected signature size, matches tar_params_best$ntop
N_PAIRS  <- 30L    # restarts sampled for the pairwise comparison in part (B)
N_SPLITS <- 20L    # random 50/50 splits for part (C)

seed_obj <- readRDS("results/desurv_seed_fits_tcgacptac.rds")
ref      <- readRDS("results/tar_fit_desurv_tcgacptac.rds")
X        <- readRDS("results/tar_data_filtered_tcgacptac.rds")$ex

fits <- seed_obj$fits
Wref <- ref$W
gn   <- rownames(Wref)
k    <- ncol(Wref)
stopifnot(nrow(Wref) == nrow(X))

## --- helpers ---------------------------------------------------------------

# Top genes per factor by the factor-specificity score s_ij = W*_ij - max_{j'!=j} W*_ij',
# where W* is W scaled by its column maxima. Mirrors DeSurv::desurv_get_top_genes.
top_idx <- function(W, n = NTOP) {
  Ws <- sweep(W, 2, apply(W, 2, max), FUN = "/")
  lapply(seq_len(ncol(Ws)), function(j) {
    other <- if (ncol(Ws) > 1) apply(Ws[, -j, drop = FALSE], 1, max) else 0
    order(Ws[, j] - other, decreasing = TRUE)[seq_len(n)]
  })
}
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# Greedy one-to-one column matching on |Spearman| of loadings. At k = 3 greedy
# and optimal assignment coincide except in pathological cases, so this avoids a
# dependency on clue::solve_LSAP.
match_cols <- function(Wa, Wb) {
  C <- outer(seq_len(ncol(Wa)), seq_len(ncol(Wb)),
             Vectorize(function(i, j) abs(cor(Wa[, i], Wb[, j], method = "spearman"))))
  perm <- integer(ncol(Wa)); avail <- seq_len(ncol(Wb))
  for (i in order(-apply(C, 1, max))) {
    j <- avail[which.max(C[i, avail])]; perm[i] <- j; avail <- setdiff(avail, j)
  }
  perm
}

adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  ci <- sum(choose(rowSums(tab), 2)); cj <- sum(choose(colSums(tab), 2))
  idx <- sum(choose(tab, 2)); expected <- ci * cj / choose(n, 2)
  (idx - expected) / ((ci + cj) / 2 - expected)
}

summarise <- function(m) {
  t(apply(m, 2, function(v) {
    v <- v[!is.na(v)]
    c(median = median(v), q25 = unname(quantile(v, .25)),
      q75 = unname(quantile(v, .75)), min = min(v), max = max(v), n = length(v))
  }))
}

ok <- which(vapply(fits, function(f) !is.null(f$W) && ncol(f$W) == k, logical(1)))

## --- (A) restart vs the reported consensus basis ----------------------------

ref_top <- top_idx(Wref)
Zref    <- t(X) %*% Wref

A_load <- A_score <- A_jac <- matrix(NA_real_, length(fits), k)
for (b in ok) {
  Wb <- fits[[b]]$W; rownames(Wb) <- gn
  perm <- match_cols(Wref, Wb); tb <- top_idx(Wb); Zb <- t(X) %*% Wb
  for (j in seq_len(k)) {
    A_load[b, j]  <- abs(cor(Wref[, j], Wb[, perm[j]], method = "spearman"))
    A_score[b, j] <- abs(cor(Zref[, j], Zb[, perm[j]], method = "spearman"))
    A_jac[b, j]   <- jaccard(ref_top[[j]], tb[[perm[j]]])
  }
}

## --- (B) restart vs restart -------------------------------------------------

set.seed(1)
prs <- t(combn(sample(ok, min(N_PAIRS, length(ok))), 2))
B_load <- B_score <- B_jac <- matrix(NA_real_, nrow(prs), k)
for (r in seq_len(nrow(prs))) {
  Wa <- fits[[prs[r, 1]]]$W; Wb <- fits[[prs[r, 2]]]$W
  rownames(Wa) <- rownames(Wb) <- gn
  perm <- match_cols(Wa, Wb); ta <- top_idx(Wa); tb <- top_idx(Wb)
  Za <- t(X) %*% Wa; Zb <- t(X) %*% Wb
  for (j in seq_len(k)) {
    B_load[r, j]  <- abs(cor(Wa[, j], Wb[, perm[j]], method = "spearman"))
    B_score[r, j] <- abs(cor(Za[, j], Zb[, perm[j]], method = "spearman"))
    B_jac[r, j]   <- jaccard(ta[[j]], tb[[perm[j]]])
  }
}

## --- (C) consensus reproducibility across disjoint restart blocks -----------

# Consensus gene partition following the SI "Consensus initialization" section:
# gene-gene co-occurrence over top-gene sets across R runs, average-linkage
# clustering on 1 - co-occurrence/R, genes appearing in < ceiling(0.3R) runs dropped.
consensus_partition <- function(idx) {
  R <- length(idx)
  Co <- matrix(0, length(gn), length(gn)); freq <- integer(length(gn))
  for (b in idx) {
    W <- fits[[b]]$W; if (is.null(W) || ncol(W) != k) next
    for (g in top_idx(W)) { Co[g, g] <- Co[g, g] + 1; freq[g] <- freq[g] + 1 }
  }
  keep <- which(freq >= ceiling(0.3 * R))
  if (length(keep) < k * 5) return(NULL)
  cl <- cutree(hclust(as.dist(1 - Co[keep, keep] / R), method = "average"), k = k)
  setNames(cl, gn[keep])
}

set.seed(7)
C_ari <- rep(NA_real_, N_SPLITS)
for (s in seq_len(N_SPLITS)) {
  sh <- sample(ok); half <- floor(length(sh) / 2)
  pa <- consensus_partition(sh[seq_len(half)])
  pb <- consensus_partition(sh[(half + 1):length(sh)])
  if (is.null(pa) || is.null(pb)) next
  shared <- intersect(names(pa), names(pb))
  if (length(shared) > 20) C_ari[s] <- adj_rand(pa[shared], pb[shared])
}
C_ari <- C_ari[!is.na(C_ari)]

## --- predictive performance across the same restarts ------------------------

ci <- unlist(seed_obj$cindex)

## --- chance baseline for the top-gene Jaccard -------------------------------
# Two independent NTOP-of-p draws: E|A n B| ~ NTOP^2/p, E|A u B| ~ 2*NTOP - NTOP^2/p
p_genes    <- length(gn)
exp_int    <- NTOP^2 / p_genes
jac_chance <- exp_int / (2 * NTOP - exp_int)

program_stability_stats <- list(
  meta = list(n_restarts = length(ok), k = k, ntop = NTOP, n_genes = p_genes,
              n_pairs = nrow(prs), n_splits = length(C_ari),
              jaccard_chance = jac_chance),
  vs_consensus = list(loading = summarise(A_load), score = summarise(A_score),
                      top_jaccard = summarise(A_jac)),
  pairwise     = list(loading = summarise(B_load), score = summarise(B_score),
                      top_jaccard = summarise(B_jac)),
  consensus_ari = c(median = median(C_ari), q25 = unname(quantile(C_ari, .25)),
                    q75 = unname(quantile(C_ari, .75)),
                    min = min(C_ari), max = max(C_ari), n = length(C_ari)),
  cindex = c(median = median(ci), min = min(ci), max = max(ci), sd = sd(ci))
)
saveRDS(program_stability_stats, "results/program_stability_stats.rds")

## --- report -----------------------------------------------------------------

pr <- function(tag, m) {
  cat(sprintf("  %-14s", tag))
  for (j in seq_len(k)) cat(sprintf("  D%d %.3f", j, m[j, "median"]))
  cat("\n")
}
cat(sprintf("\nProgram stability across %d cached random restarts (k = %d, ntop = %d)\n",
            length(ok), k, NTOP))
cat("\n(A) restart vs reported consensus basis (medians)\n")
pr("loading r",  program_stability_stats$vs_consensus$loading)
pr("score r",    program_stability_stats$vs_consensus$score)
pr("top Jaccard",program_stability_stats$vs_consensus$top_jaccard)
cat(sprintf("\n(B) restart vs restart, %d independent pairs (medians)\n", nrow(prs)))
pr("loading r",  program_stability_stats$pairwise$loading)
pr("score r",    program_stability_stats$pairwise$score)
pr("top Jaccard",program_stability_stats$pairwise$top_jaccard)
cat(sprintf("\n     chance top-%d Jaccard for two independent draws from %d genes: %.3f\n",
            NTOP, p_genes, jac_chance))
cat(sprintf("\n(C) consensus partition ARI between disjoint 50-run blocks, %d splits\n",
            length(C_ari)))
cat(sprintf("     median %.3f   IQR %.3f-%.3f   range %.3f-%.3f\n",
            median(C_ari), quantile(C_ari, .25), quantile(C_ari, .75), min(C_ari), max(C_ari)))
cat(sprintf("\nTraining C-index across the same restarts: median %.4f, range %.4f-%.4f, SD %.4f\n",
            median(ci), min(ci), max(ci), sd(ci)))
cat("\nRead (B) against the C-index SD: near-identical predictive performance from\n")
cat("bases sharing little top-gene membership is the signature of non-identifiability.\n")
cat("\nSaved -> results/program_stability_stats.rds\n")
