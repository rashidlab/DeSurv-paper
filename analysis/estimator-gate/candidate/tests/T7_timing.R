# T7: realistic size (p = 3000, k = 3) at n = 175 (inner-fold training size) and n = 218 (outer training size):
# one 30-start fit each at the midpoint of the frozen tuning box (logit alpha, log10 lambda, nu, log10 rho).
# Per-start time distribution, iterations, backtracking counts, FISTA iterations, statuses, r at stop, eps
# shared across starts, parallel efficiency; extrapolated CPU-hours for the gate.
# Amendment 3: W step = Newton-CG (w_solver "ncg", the palm_fit default); maxit 5000 -> the per-start budget (MAXIT
# is the amendment-3 sweep budget of 1,000 sweeps, plus 25,000 HVPs, not the old PG iteration count). The outcome is
# reported; Q1 (tests/Q1_warm.R), not T7, is the binding realistic-size criterion. Iterations below are sweeps.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
NC <- as.integer(Sys.getenv("T_NCORES", "15")); P <- 3000L; NS <- c(175L, 218L); NST <- 30L; MAXIT <- 1000L
if (QUICK) { P <- 300L; NS <- c(50L, 60L); NST <- 6L; MAXIT <- 50L; NC <- min(NC, 6L) }
hp <- list(alpha = plogis((qlogis(1e-4) + qlogis(0.9)) / 2), lambda = 1, nu = 0.475, rho = 10^0.5)
runs <- list()
for (n in NS) {
  s <- sim_data(700, P, n); D <- palm_data(s$X, s$y, s$d)
  t0 <- proc.time()[["elapsed"]]
  fit <- palm_fit(D, hp, k = 3L, seed = 123L, nstarts = NST, maxit = MAXIT, ncores = NC, keep_starts = (n == max(NS)))
  wall <- proc.time()[["elapsed"]] - t0
  runs[[as.character(n)]] <- list(n = n, ne = D$ne, ties = s$n_tied, fit = fit, wall = wall)
  if (n == max(NS)) saveRDS(list(fit = fit, W_true = s$W_true, wall = wall, ncores = NC), file.path(OUT, "T7_fit.rds"))
}
open_out("T7_timing.txt")
cat("T7 realistic size and timing |", format(Sys.time()), "| cores", NC, "| BLAS", basename(extSoftVersion()[["BLAS"]]), "\n")
cat(sprintf("hp (box midpoint): alpha=%.6g lambda=%g nu=%g rho=%.6g lambdaH=1e-6; p=%d, k=3, %d starts, w_solver ncg, per-start budget %d sweeps + 25000 HVPs, fit seed 123\n",
            hp$alpha, hp$lambda, hp$nu, hp$rho, P, NST, MAXIT))
cat("note: the machine is shared; other users' jobs were running, so seconds are upper-biased.\n")
pass <- TRUE; ms <- c()
for (z in runs) {
  st <- z$fit$starts; q <- quantile(st$secs, c(0, 0.25, 0.5, 0.75, 0.9, 1))
  shared <- length(unique(st$eps)) == 1 && st$eps[1] == z$fit$eps
  cat(sprintf("\nn = %d (events %d, tied times %d): fit status %s, returned start %s\n", z$n, z$ne, z$ties, z$fit$status, z$fit$best_start))
  cat(sprintf("  eps %.6e shared by all %d starts: %s\n", z$fit$eps, nrow(st), shared))
  cat(sprintf("  start status: converged %d, nonconverged %d, failure %d\n", sum(st$status == "converged"), sum(st$status == "nonconverged"), sum(st$status == "failure")))
  if (any(st$status != "converged")) { cat("  reasons: "); print(table(st$reason[st$status != "converged"])) }
  cat(sprintf("  HVPs per start: median %g, max %g (cap 25000); W exits summed: Newton %g, no-op %g, line-search failed %g\n",
              median(st$hvp_total), max(st$hvp_total), sum(st$w_newton), sum(st$w_noop), sum(st$w_ls_failed)))
  cat(sprintf("  seconds per start: mean %.1f; quantiles 0/25/50/75/90/100%%: %s\n", mean(st$secs), paste(sprintf("%.1f", q), collapse = " / ")))
  cat(sprintf("  iterations: median %g, max %g; seconds per iteration median %.4f\n", median(st$iter), max(st$iter), median(st$secs / st$iter)))
  cat(sprintf("  W halvings per start: median %g, max %g; FISTA iterations per start: median %g, max %g; rejected b steps total %d\n",
              median(st$w_backtracks), max(st$w_backtracks), median(st$b_iters), max(st$b_iters), sum(st$b_rejected)))
  cat(sprintf("  r at stop: median %.3e (r_W %.3e, r_H %.3e, r_b %.3e), r_original median %.3e\n", median(st$r), median(st$rW), median(st$rH), median(st$rb), median(st$r_original)))
  cat(sprintf("  per 30-start fit: summed CPU %.1f s, wall %.1f s on %d cores, parallel efficiency %.2f\n", sum(st$secs), z$wall, NC, sum(st$secs) / (z$wall * NC)))
  ms[as.character(z$n)] <- mean(st$secs)
  pass <- pass && shared && nrow(st) == NST && !any(st$status == "failure")
}
n_inner <- 20 * 180 * 5 * 30; n_final <- 20 * 30
cpu <- (n_inner * ms[[1]] + n_final * ms[[2]]) / 3600
cat(sprintf("\ngate extrapolation: %d inner starts (20 x 180 x 5 x 30) x %.1f s (n = %d) + %d final starts (20 x 30) x %.1f s (n = %d) = %.0f CPU-hours\n",
            n_inner, ms[[1]], NS[1], n_final, ms[[2]], NS[2], cpu))
cat("caveat: one hyperparameter point (box midpoint); BO evaluations span the box and time per start depends on alpha, lambda, rho and on how many starts run to maxit.\n")
for (z in runs) { cat(sprintf("\nper start, n = %d:\n", z$n)); print(z$fit$starts[, c("start", "status", "reason", "F", "iter", "r", "rW", "rH", "rb", "r_original", "hvp_total", "w_backtracks", "b_iters", "b_rejected", "eps", "secs")], digits = 4) }
verdict(pass)
sink()
