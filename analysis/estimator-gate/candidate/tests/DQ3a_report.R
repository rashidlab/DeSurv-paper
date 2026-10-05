# D-Q3a report, regenerated from tests/out/DQ3a_continuation.rds (supersedes the report block of DQ3a_continuation.R:
# adds the sweep-1,000 baseline from D-Q1 and tolerates a final trace row without a b-exit record).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
hp <- readRDS(file.path(CAND, "tests", "out", "superseded", "pre-amendment3", "T7_fit.rds"))$fit$hp
dq1 <- readRDS(file.path(OUT, "DQ1_block_replay.rds")); x <- readRDS(file.path(OUT, "DQ3a_continuation.rds"))
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
SS <- c(1L, 5L, 16L); u <- rep(1, nrow(D$X)) / sqrt(nrow(D$X))
row_of <- function(W, b) {
  sc <- palm_scores(W, D, hp$eps); C <- cor(sc$Zc)
  c(sigma = sc$sigma, beta = b / sc$sigma, b = b, maxcor = max(abs(C[upper.tri(C)])), cosu = as.vector(crossprod(W, u)))
}
line <- function(sweep, F, r, rW) cat(sprintf("  sweep %5d: F %.10f | sigma %s | beta %s | max|cor| %.5f | cos(w,u) %s | median rW (last 100) %.3e\n",
  sweep, F, paste(sprintf("%.4g", r[1:3]), collapse = "/"), paste(sprintf("%.4g", r[4:6]), collapse = "/"), r[["maxcor"]],
  paste(sprintf("%.4f", r[c("cosu1", "cosu2", "cosu3")]), collapse = "/"), rW))

sink(file.path(OUT, "DQ3a_continuation.txt"), split = TRUE)
cat("D-Q3a continuation (descriptive). Gate:", paste(SS, ifelse(x$gate, "PASS", "FAIL"), collapse = "; "),
    "| settings: maxit", x$settings$maxit, "hvp_cap", x$settings$hvp_cap, "tol", x$settings$tol, "\n")
cat("Sweep numbers are totals: D-Q1's 1,000 plus the continuation. Sweep 1000 = D-Q1 end point.\n")
for (r in x$res) {
  if (inherits(r, "try-error")) { cat("ERROR:", r, "\n"); next }
  f <- r$f; tr <- f$trace; rec <- r$rec; n <- nrow(tr); nr <- nrow(rec)
  e <- dq1$res[[match(r$start, SS)]]; ef <- e$f; eb <- e$blk[e$blk$block == "b", ]
  if (nr != n) cat(sprintf("  NOTE: %d b-exit records for %d sweeps (final sweep(s) ended before the b block)\n", nr, n))
  J <- which(abs(f$b) > 0); if (!length(J)) J <- 1:3
  cat(sprintf("\n=== start %d: %d continuation sweeps, status %s (%s), %.0f s, HVPs %d ===\n",
              r$start, n, f$status, f$reason, r$secs, f$hvp_total))
  line(1000L, ef$F, row_of(ef$W, ef$b), median(tail(eb$rW, 100)))
  at <- intersect(c(2000L, 5000L, 10000L, 20000L) - 1000L, seq_len(nr)); at <- unique(c(at, nr))
  for (k in at) line(k + 1000L, tr[k, "Fb"], rec[k, ], median(tr[max(1L, k - 99L):k, "rW"]))
  ls <- log(apply(rec[, J, drop = FALSE], 1, min)); sw <- seq_len(nr) + 1000
  sl <- function(lo, hi) { i <- which(sw >= lo & sw <= hi); if (length(i) < 10) NA else coef(lm(ls[i] ~ log(sw[i])))[2] }
  cat(sprintf("  log-log slope of min coefficient-carrying sigma: sweeps 1000-5000 %.3f; 5000-20000 %.3f\n",
              sl(1000, 5000), sl(5000, 20000)))
  Fb <- c(ef$F, tr[seq_len(nr), "Fb"]); blk <- split(seq_len(nr), ceiling(seq_len(nr) / 1000))
  cat("  F drop per 1,000 sweeps:", paste(sprintf("%.2e", vapply(blk, function(i) Fb[min(i)] - Fb[max(i) + 1], 0)), collapse = " "), "\n")
  dls <- -diff(c(log(min(ef$sigma[J])), ls)); long <- tr[seq_len(nr), "cg_iter"] >= 10
  cat(sprintf("  share of log-sigma drop from sweeps with >= 10 CG iterations: %.3f (%d such sweeps of %d)\n",
              sum(dls[long]) / sum(dls), sum(long), nr))
  sg <- sign(rec[, c("b1", "b2", "b3")])
  cat("  b sign pattern changes:", sum(rowSums(abs(diff(sg))) > 0), "; final b", paste(sprintf("%.4f", f$b), collapse = ", "), "\n")
  cat(sprintf("  final r: rW %.3e rH %.3e rb %.3e; face changes enter/leave %d/%d; negative curvature %d\n",
              f$r[["rW"]], f$r[["rH"]], f$r[["rb"]], sum(tr[, "n_enter"]), sum(tr[, "n_leave"]), sum(tr[, "negative_curvature"])))
}
sink()
