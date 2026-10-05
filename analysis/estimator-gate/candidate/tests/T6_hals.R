# T6: HALS H block. Every sweep must not increase the H part (1e-13 relative, I4); repeated sweeps must reach
# the exact regularized NNLS solution, computed independently per column by enumerating the 2^k supports.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("T6_hals.txt")
cat("T6 HALS |", format(Sys.time()), "\n")
s <- sim_data(16, 80, 50); D <- palm_data(s$X, s$y, s$d)
Hpart <- function(W, H, lH) sum((D$X - W %*% H)^2) / (2 * D$X2) + lH * sum(H^2) / (2 * D$n * ncol(W))
exact_H <- function(W, lH) {
  k <- ncol(W); cc <- lH * D$X2 / (D$n * k); sup <- as.matrix(expand.grid(rep(list(0:1), k)))[-1, , drop = FALSE]
  sapply(seq_len(D$n), function(i) {
    x <- D$X[, i]; best <- rep(0, k); bv <- sum(x^2)
    for (r in seq_len(nrow(sup))) {
      S <- which(sup[r, ] == 1); h <- numeric(k)
      h[S] <- solve(crossprod(W[, S, drop = FALSE]) + cc * diag(length(S)), crossprod(W[, S, drop = FALSE], x))
      if (all(h >= 0)) { v <- sum((x - W %*% h)^2) + cc * sum(h^2); if (v < bv) { bv <- v; best <- h } }
    }
    best
  })
}
run <- function(sweep_fun, W, lH, nsw = 5000L, seed = 1L) {
  set.seed(seed); H <- matrix(runif(3 * D$n, 0, 50), 3, D$n); hp <- list(lambdaH = lH, alpha = 0.3)
  f <- Hpart(W, H, lH); ninc <- 0L; maxinc <- 0
  for (i in seq_len(nsw)) {
    H <- sweep_fun(W, H, D, hp); fn <- Hpart(W, H, lH)
    maxinc <- max(maxinc, fn - f); if (fn - f > 1e-13 * max(1, abs(f))) ninc <- ninc + 1L; f <- fn
  }
  He <- exact_H(W, lH)
  kkt <- sqrt(sum(pmin(H, palm_grad_H(W, H, D, hp))^2)) / sqrt(sum(H^2))
  list(ninc = ninc, maxinc = maxinc, relerr = max(abs(H - He)) / max(abs(He)), kkt = kkt, nsw = nsw)
}
setups <- list()
for (lH in c(1e-6, 1e-2, 1)) for (wt in c("unit", "nonunit")) for (rep_ in 1:3) {
  set.seed(600 + rep_); W <- proj_S(matrix(runif(D$p * 3)^(1 + rep_), D$p, 3))
  if (wt == "nonunit") W <- sweep(W, 2, c(0.5, 2, 3), "*")
  setups[[length(setups) + 1]] <- list(lH = lH, wt = wt, rep = rep_, W = W)
}
res <- do.call(rbind, lapply(setups, function(st) {
  r <- run(hals_sweep, st$W, st$lH)
  data.frame(lambdaH = st$lH, W = st$wt, rep = st$rep, sweeps = r$nsw, increases = r$ninc, max_raw_increase = r$maxinc,
             relerr_vs_exact = r$relerr, kkt = r$kkt)
}))
print(res, digits = 3)
cat(sprintf("\nmonotonicity: %d sweeps checked, %d increases, largest raw change %.3e\n", sum(res$sweeps), sum(res$increases), max(res$max_raw_increase)))
cat(sprintf("exact NNLS: %d setups, %d with max rel error < 1e-8, worst %.3e; KKT residual worst %.3e (%d < 1e-10)\n",
            nrow(res), sum(res$relerr_vs_exact < 1e-8), max(res$relerr_vs_exact), max(res$kkt), sum(res$kkt < 1e-10)))

nT6 <- 2 * nrow(res); pT6 <- sum(res$increases == 0) + sum(res$relerr_vs_exact < 1e-8 & res$kkt < 1e-10)
cat(sprintf("T6 TOTAL: %d checks (monotone runs + exact-solution matches), %d passed\n", nT6, pT6))
verdict(nT6 == pT6)
# Breakage check (B8): wrong denominators must miss the exact solution.
cat("\nbreakage check (B8):\n")
bad_k <- function(W, H, D, hp) {             # lambdaH term without the 1/k
  WtX <- crossprod(W, D$X); WtW <- crossprod(W)
  for (j in 1:3) H[j, ] <- pmax(0, (WtX[j, ] - WtW[j, -j, drop = FALSE] %*% H[-j, , drop = FALSE]) / (WtW[j, j] + hp$lambdaH * D$X2 / D$n))
  H
}
bad_unit <- function(W, H, D, hp) {          # assumes ||w_j|| = 1
  WtX <- crossprod(W, D$X); WtW <- crossprod(W)
  for (j in 1:3) H[j, ] <- pmax(0, (WtX[j, ] - WtW[j, -j, drop = FALSE] %*% H[-j, , drop = FALSE]) / (1 + hp$lambdaH * D$X2 / (D$n * 3)))
  H
}
st <- setups[[which(sapply(setups, function(z) z$lH == 1 && z$wt == "unit" && z$rep == 1))]]
r1 <- run(bad_k, st$W, 1, 2000L)
st2 <- setups[[which(sapply(setups, function(z) z$lH == 1e-6 && z$wt == "nonunit" && z$rep == 1))]]
r2 <- run(bad_unit, st2$W, 1e-6, 2000L)
cat(sprintf("  missing 1/k (lambdaH = 1):         rel error vs exact %.2e -> %s\n", r1$relerr, if (r1$relerr >= 1e-8) "flagged" else "NOT flagged"))
cat(sprintf("  unit-norm assumed (non-unit W):     rel error vs exact %.2e -> %s\n", r2$relerr, if (r2$relerr >= 1e-8) "flagged" else "NOT flagged"))
sink()
