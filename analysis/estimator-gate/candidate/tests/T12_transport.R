# T12: transport of the stored prediction vector (spec section 4, amendment 2). v~ = P W diag(1/sigma) b over
# the training gene universe, centered before the ntop support G is taken; scores use X_G^T v~_G. Removing or
# reordering genes outside G in an external (simulated) matrix leaves the score unchanged; re-ranking within
# a reduced gene universe changes the input ranks and is reported descriptively.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("T12_transport.txt")
cat("T12 transport |", format(Sys.time()), "\n")
nchk <- 0L; npass <- 0L
ok <- function(cond, label) { nchk <<- nchk + 1L; npass <<- npass + isTRUE(cond); cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "pass" else "FAIL", label)) }
s <- sim_data(1201, 200, 160, beta = c(2, 1, 0)); tr <- 1:110; ex <- 111:160
rk <- palm_rank
D <- palm_fold_data(s$X_raw, s$y, s$d, tr)
fit <- palm_fit(D, list(alpha = 0.5, lambda = 0.1, nu = 0.5, rho = 1), seed = 1202L, nstarts = 4L, maxit = 300L)
m <- palm_model(fit, D, ntop = 25)
v_ind <- as.vector(fit$best$W %*% (fit$best$b / palm_scores(fit$best$W, D, fit$eps)$sigma))
ok(max(abs(m$v_tilde - (v_ind - mean(v_ind)))) < 1e-15 * max(abs(v_ind)) * 10 && identical(names(m$v_tilde), D$genes),
   "stored v~ = v - mean(v) over the full training gene universe (independent recomputation)")
ok(abs(sum(m$v_tilde)) < 1e-12 * sum(abs(m$v_tilde)), sprintf("v~ sums to 0 over the universe (|sum| %.1e)", abs(sum(m$v_tilde))))
cat(sprintf("  info: |G| = %d of %d genes; mean of v~ over G = %.3e (not recentered on G)\n", length(m$G), D$p, mean(m$v_tilde[m$G])))
Xe <- rk(s$X_raw[, ex]); s0 <- palm_score(Xe, m)
ok(identical(s0, as.vector(crossprod(Xe[m$G, ], m$v_tilde[m$G]))), "score = X_G^T v~_G")
set.seed(1203); nonG <- setdiff(rownames(Xe), m$G)
Xr1 <- Xe[sample(rownames(Xe)), ]
ok(identical(palm_score(Xr1, m), s0), "all genes reordered: score identical")
Xr2 <- Xe[c(m$G, sample(nonG, length(nonG) %/% 2)), ]; Xr2 <- Xr2[sample(nrow(Xr2)), ]
ok(identical(palm_score(Xr2, m), s0), "half of the non-G genes removed and rows reordered: score identical")
Xr3 <- Xe[rev(m$G), , drop = FALSE]
ok(identical(palm_score(Xr3, m), s0), "only G kept (reversed order): score identical")
ok(tryCatch({ palm_score(Xe[setdiff(rownames(Xe), m$G[1]), ], m); FALSE }, error = function(e) grepl("missing", conditionMessage(e))),
   "a scored gene missing: deterministic error, no silent rescoring")
cat("\ndescriptive: re-ranking within a reduced gene universe (rank construction, not centering)\n")
for (fr in c(0.75, 0.5, 0.25)) {
  keep <- c(m$G, sample(nonG, round(fr * length(nonG))))
  sR <- palm_score(rk(s$X_raw[keep, ex]), m)
  cat(sprintf("  universe %3.0f%% of non-G genes kept (%d genes): Pearson %.4f, Spearman %.4f with full-universe scores, mean shift %.3g (sd of full scores %.3g)\n",
              100 * fr, length(keep), cor(sR, s0), cor(sR, s0, method = "spearman"), mean(sR - s0), sd(s0)))
}
cat(sprintf("\nT12 TOTAL: %d checks, %d passed\n", nchk, npass))
verdict(nchk == npass)
sink()
