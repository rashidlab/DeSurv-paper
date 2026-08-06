#!/usr/bin/env Rscript
# 21_consensus_endpoint.R
# ---------------------------------------------------------------------------
# Consensus ENDPOINT reproducibility (Nature Cancer revision).
#
# code/20 measures the dispersion of individual random restarts and the
# reproducibility of the consensus gene PARTITION. Neither characterises the
# object the manuscript actually reports. The reported basis is produced by
#     consensus over R restarts  ->  initialization W0  ->  desurv_fit(W0)
# so consensus supplies a starting point and the model is then fitted. The
# relevant question for the paper's "reproducible" claim is whether two
# independently derived consensus initializations converge to the SAME FITTED
# BASIS, which code/20 cannot answer.
#
# Design (split-half). Partition the 100 cached restarts in
# results/desurv_seed_fits_tcgacptac.rds into two DISJOINT blocks of 50, build a
# consensus initialization from each block with DeSurv::desurv_consensus_seed
# (min_frequency = ceiling(0.3 * 50), matching code/04), fit the full model from
# each initialization at the BO-selected hyperparameters, and compare the two
# fitted bases. Repeat over several random splits.
#
# This reuses the cached restarts, so each split costs two final fits, not 100.
# It is a split-half proxy for "two independent 100-restart consensus runs":
# strictly it halves R, which if anything makes consensus noisier than the
# reported R = 100 procedure, so a GOOD result here is conservative evidence and
# a BAD result is not automatically damning. Stated plainly rather than glossed.
#
# Requires the DeSurv package. It is installed under the R 4.5 library tree and
# loads cleanly under R 4.6; the .libPaths() call below is what code/16 already
# relies on.
# ---------------------------------------------------------------------------

N_SPLITS <- as.integer(Sys.getenv("ENDPOINT_SPLITS", "6"))
NTOP     <- 270L

.lib45 <- "/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.5"
if (dir.exists(.lib45)) .libPaths(c(.lib45, .libPaths()))
suppressMessages(library(DeSurv))

seed_obj <- readRDS("results/desurv_seed_fits_tcgacptac.rds")
ref      <- readRDS("results/tar_fit_desurv_tcgacptac.rds")
dat      <- readRDS("results/tar_data_filtered_tcgacptac.rds")
prm      <- readRDS("results/tar_params_best_tcgacptac.rds")

fits <- seed_obj$fits
gn   <- rownames(ref$W)
k    <- as.integer(prm$k)
X    <- dat$ex; y <- dat$sampInfo$time; d <- dat$sampInfo$event

RUN_TOL   <- 1e-7
RUN_MAXIT <- 3000L

## --- helpers (same conventions as code/20) ---------------------------------

top_idx <- function(W, n = NTOP) {
  Ws <- sweep(W, 2, apply(W, 2, max), FUN = "/")
  lapply(seq_len(ncol(Ws)), function(j) {
    other <- if (ncol(Ws) > 1) apply(Ws[, -j, drop = FALSE], 1, max) else 0
    order(Ws[, j] - other, decreasing = TRUE)[seq_len(n)]
  })
}
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
all_perms <- function(k) {
  if (k == 1) return(list(1L))
  out <- list()
  for (i in seq_len(k)) for (p in all_perms(k - 1)) {
    rest <- setdiff(seq_len(k), i); out[[length(out) + 1]] <- c(i, rest[p])
  }
  out
}
match_cols <- function(Wa, Wb) {
  C <- outer(seq_len(ncol(Wa)), seq_len(ncol(Wb)),
             Vectorize(function(i, j) abs(cor(Wa[, i], Wb[, j], method = "spearman"))))
  perms <- all_perms(ncol(Wa))
  as.integer(perms[[which.max(vapply(perms, function(p)
    sum(C[cbind(seq_len(ncol(Wa)), p)]), numeric(1)))]])
}

# consensus initialization -> final fit, mirroring code/04_fit_models.R
consensus_fit <- function(idx) {
  init <- DeSurv::desurv_consensus_seed(
    fits = fits[idx], X = X, ntop = NTOP, k = k,
    min_frequency = max(1L, as.integer(ceiling(0.3 * length(idx))))
  )
  desurv_fit(X = X, y = y, d = d, k = k, alpha = prm$alpha, lambda = prm$lambda,
             nu = prm$nu, lambdaW = 0, lambdaH = 0,
             W0 = init$W0, H0 = init$H0, beta0 = init$beta0,
             seed = NULL, tol = RUN_TOL / 100, maxit = RUN_MAXIT, verbose = FALSE)
}

ok <- which(vapply(fits, function(f) !is.null(f$W) && ncol(f$W) == k, logical(1)))

## --- run -------------------------------------------------------------------

set.seed(101)
load_m <- score_m <- jac_m <- matrix(NA_real_, N_SPLITS, k)
ci_a <- ci_b <- rep(NA_real_, N_SPLITS)
ref_load <- ref_jac <- matrix(NA_real_, 2 * N_SPLITS, k)

for (s in seq_len(N_SPLITS)) {
  sh <- sample(ok); half <- floor(length(sh) / 2)
  fa <- try(consensus_fit(sh[seq_len(half)]), silent = TRUE)
  fb <- try(consensus_fit(sh[(half + 1):length(sh)]), silent = TRUE)
  if (inherits(fa, "try-error") || inherits(fb, "try-error")) {
    message(sprintf("  split %d: fit failed, skipping", s)); next
  }
  Wa <- fa$W; Wb <- fb$W; rownames(Wa) <- rownames(Wb) <- gn
  ci_a[s] <- fa$cindex; ci_b[s] <- fb$cindex

  p <- match_cols(Wa, Wb); ta <- top_idx(Wa); tb <- top_idx(Wb)
  Za <- t(X) %*% Wa; Zb <- t(X) %*% Wb
  for (j in seq_len(k)) {
    load_m[s, j]  <- abs(cor(Wa[, j], Wb[, p[j]], method = "spearman"))
    score_m[s, j] <- abs(cor(Za[, j], Zb[, p[j]], method = "spearman"))
    jac_m[s, j]   <- jaccard(ta[[j]], tb[[p[j]]])
  }
  # each half-consensus fit against the REPORTED basis
  rt <- top_idx(ref$W)
  for (ii in seq_len(2)) {
    W <- if (ii == 1) Wa else Wb
    pr <- match_cols(ref$W, W); tt <- top_idx(W)
    for (j in seq_len(k)) {
      ref_load[2 * (s - 1) + ii, j] <- abs(cor(ref$W[, j], W[, pr[j]], method = "spearman"))
      ref_jac[2 * (s - 1) + ii, j]  <- jaccard(rt[[j]], tt[[pr[j]]])
    }
  }
  message(sprintf("  split %d done: C %.3f / %.3f, median top-%d Jaccard %.3f",
                  s, fa$cindex, fb$cindex, NTOP, median(jac_m[s, ])))
}

summarise <- function(m) t(apply(m, 2, function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(c(median = NA, q25 = NA, q75 = NA, min = NA, max = NA, n = 0))
  c(median = median(v), q25 = unname(quantile(v, .25)), q75 = unname(quantile(v, .75)),
    min = min(v), max = max(v), n = length(v))
}))

consensus_endpoint_stats <- list(
  meta = list(n_splits = N_SPLITS, block_size = floor(length(ok) / 2), k = k,
              ntop = NTOP, design = "split-half of 100 cached restarts",
              caveat = "halves R relative to the reported R=100 procedure"),
  half_vs_half = list(loading = summarise(load_m), score = summarise(score_m),
                      top_jaccard = summarise(jac_m)),
  vs_reported  = list(loading = summarise(ref_load), top_jaccard = summarise(ref_jac)),
  cindex = c(half_a = median(ci_a, na.rm = TRUE), half_b = median(ci_b, na.rm = TRUE),
             reported = ref$cindex)
)
saveRDS(consensus_endpoint_stats, "results/consensus_endpoint_stats.rds")

pr <- function(tag, m) {
  cat(sprintf("  %-14s", tag))
  for (j in seq_len(k)) cat(sprintf("  D%d %.3f", j, m[j, "median"]))
  cat("\n")
}
cat(sprintf("\nConsensus ENDPOINT reproducibility, %d split-half replicates\n", N_SPLITS))
cat("\nTwo independently seeded consensus fits vs each other (medians)\n")
pr("loading r",   consensus_endpoint_stats$half_vs_half$loading)
pr("score r",     consensus_endpoint_stats$half_vs_half$score)
pr("top Jaccard", consensus_endpoint_stats$half_vs_half$top_jaccard)
cat("\nEach half-consensus fit vs the REPORTED basis (medians)\n")
pr("loading r",   consensus_endpoint_stats$vs_reported$loading)
pr("top Jaccard", consensus_endpoint_stats$vs_reported$top_jaccard)
cat(sprintf("\nC-index: half A %.4f, half B %.4f, reported %.4f\n",
            median(ci_a, na.rm = TRUE), median(ci_b, na.rm = TRUE), ref$cindex))
cat("\nCompare against code/20: restart-vs-restart top-270 Jaccard was 0.13\n")
cat("(chance 0.074). Movement toward 1.0 here is the consensus step doing its job.\n")
cat("\nSaved -> results/consensus_endpoint_stats.rds\n")
