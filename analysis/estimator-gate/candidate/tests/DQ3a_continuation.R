# D-Q3a (docs/plans/2026-10-05-dq3-collapse-trajectory-prespec.md): continue starts 1, 5, 16 from their D-Q1 end
# points to 20,000 total sweeps with the unchanged amendment-3 solver; hvp_cap raised to 1e6 (diagnostic setting only).
# Descriptive. Per-sweep record at b-block exit: sigma, beta, b, max |cor| of centered scores, cos(w_j, u).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressMessages(library(parallel))

hp <- readRDS(file.path(CAND, "tests", "out", "superseded", "pre-amendment3", "T7_fit.rds"))$fit$hp
dq1 <- readRDS(file.path(OUT, "DQ1_block_replay.rds")); q1 <- readRDS(file.path(OUT, "Q1_warm.rds"))
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
SS <- c(1L, 5L, 16L); MAXIT <- as.integer(Sys.getenv("DQ3_MAXIT", "19000"))  # env override for the smoke test only

# gate: saved D-Q1 end points equal the committed Q1 end points
gate <- vapply(SS, function(i) {
  a <- dq1$res[[match(i, SS)]]$f; b <- q1$fits[[match(i, q1$tab$start)]]
  identical(a$W, b$W) && identical(a$H, b$H) && identical(a$b, b$b) && identical(a$F, b$F)
}, logical(1))
cat("GATE (D-Q1 end points == committed Q1 end points):", paste(SS, ifelse(gate, "PASS", "FAIL"), collapse = "; "), "\n")
if (!all(gate)) stop("D-Q3a gate failed; stopping as pre-specified")

solve_b_plain <- solve_b; b_accept_plain <- PALM_DEFS$b_accept
u <- rep(1, nrow(D$X)) / sqrt(nrow(D$X))

one <- function(i) {
  REC <- new.env(); REC$rows <- vector("list", MAXIT); REC$n <- 0L
  solve_b <<- function(b0, W, D, hp, L0 = 1, tol = 1e-10, maxit = 1e5) {
    res <- solve_b_plain(b0, W, D, hp, L0, tol, maxit); REC$pending <- list(W = W, b_old = b0, b_cand = res$b); res
  }
  PALM_DEFS$b_accept <<- function(F_new, F_old) {
    dec <- b_accept_plain(F_new, F_old); pd <- REC$pending
    bf <- if (isTRUE(dec)) pd$b_cand else pd$b_old
    sc <- palm_scores(pd$W, D, hp$eps); C <- cor(sc$Zc)
    REC$n <- REC$n + 1L
    REC$rows[[REC$n]] <- c(sigma = sc$sigma, beta = bf / sc$sigma, b = bf, maxcor = max(abs(C[upper.tri(C)])),
                           cosu = as.vector(crossprod(pd$W, u)))
    dec
  }
  e <- dq1$res[[match(i, SS)]]$f
  t0 <- proc.time()[["elapsed"]]
  f <- palm_fit_one(D, hp, e$W, e$H, e$b, maxit = MAXIT, trace = TRUE, hvp_cap = 1e6L)
  rec <- do.call(rbind, REC$rows[seq_len(REC$n)])
  list(start = i, f = f, rec = rec, secs = proc.time()[["elapsed"]] - t0)
}
res <- mclapply(SS, one, mc.cores = 3L, mc.preschedule = FALSE)
saveRDS(list(res = res, gate = gate, settings = list(maxit = MAXIT, hvp_cap = 1e6, tol = 1e-5)), file.path(OUT, if (MAXIT < 19000L) "DQ3a_smoke.rds" else "DQ3a_continuation.rds"))

sink(file.path(OUT, if (MAXIT < 19000L) "DQ3a_smoke.txt" else "DQ3a_continuation.txt"), split = TRUE)
cat("D-Q3a continuation (descriptive). Sweep numbers are totals: D-Q1's 1,000 plus the continuation.\n")
for (r in res) {
  if (inherits(r, "try-error")) { cat("ERROR:", r, "\n"); next }
  f <- r$f; tr <- f$trace; rec <- r$rec; n <- nrow(tr)
  if (nrow(rec) != n) cat(sprintf("  WARNING: %d b-exit records for %d sweeps\n", nrow(rec), n))
  J <- which(abs(f$b) > 0); if (!length(J)) J <- 1:3
  cat(sprintf("\n=== start %d: %d continuation sweeps, status %s (%s), %.0f s, HVPs %d ===\n",
              r$start, n, f$status, f$reason, r$secs, f$hvp_total))
  at <- intersect(c(1000L, 2000L, 5000L, 10000L, 20000L) - 1000L, seq_len(n)); at <- unique(c(pmax(at, 1L), n))
  for (k in at) {
    w <- max(1L, k - 99L):k
    cat(sprintf("  sweep %5d: F %.10f | sigma %s | beta %s | max|cor| %.5f | cos(w,u) %s | median rW (last 100) %.3e\n",
      k + 1000L, tr[k, "Fb"], paste(sprintf("%.4g", rec[k, 1:3]), collapse = "/"), paste(sprintf("%.4g", rec[k, 4:6]), collapse = "/"),
      rec[k, "maxcor"], paste(sprintf("%.4f", rec[k, c("cosu1", "cosu2", "cosu3")]), collapse = "/"), median(tr[w, "rW"])))
  }
  ls <- log(apply(rec[, J, drop = FALSE], 1, min)); sw <- seq_len(n) + 1000
  sl <- function(lo, hi) { i <- which(sw >= lo & sw <= hi); if (length(i) < 10) NA else coef(lm(ls[i] ~ log(sw[i])))[2] }
  cat(sprintf("  log-log slope of min coefficient-carrying sigma: sweeps 1000-5000 %.3f; 5000-20000 %.3f\n",
              sl(1000, 5000), sl(5000, 20000)))
  Fb <- c(dq1$res[[match(r$start, SS)]]$f$F, tr[, "Fb"]); blk <- split(seq_len(n), ceiling(seq_len(n) / 1000))
  cat("  F drop per 1,000 sweeps:", paste(sprintf("%.2e", vapply(blk, function(i) Fb[min(i)] - Fb[max(i) + 1], 0)), collapse = " "), "\n")
  dls <- -diff(c(log(min(dq1$res[[match(r$start, SS)]]$f$sigma[J])), ls))
  cat(sprintf("  share of log-sigma drop from sweeps with >= 10 CG iterations: %.3f (%d such sweeps of %d)\n",
              sum(dls[tr[, "cg_iter"] >= 10]) / sum(dls), sum(tr[, "cg_iter"] >= 10), n))
  sg <- sign(rec[, 7:9]); cat("  b sign pattern changes:", sum(rowSums(abs(diff(sg))) > 0), "; final b", paste(sprintf("%.4f", f$b), collapse = ", "), "\n")
  cat(sprintf("  final r: rW %.3e rH %.3e rb %.3e; face changes enter/leave %d/%d; negative curvature %d\n",
              f$r[["rW"]], f$r[["rH"]], f$r[["rb"]], sum(tr[, "n_enter"]), sum(tr[, "n_leave"]), sum(tr[, "negative_curvature"])))
}
sink()
