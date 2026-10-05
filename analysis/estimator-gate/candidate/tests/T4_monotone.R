# T4: full runs to the stopping rule (r < 1e-5, the amendment-3 per-start budget, or numerical failure) on 3
# simulated datasets x 9 corners x 5 starts; eps pooled over the 5 starts of each (dataset, corner) fit. Checks: no
# F increase after any block (counter, and an independent recount from the per-block trace), no numerical failure,
# all residuals recorded (r_W, r_H, r_b, r_original; r = max of the first three), accepted W steps t <= 1 on the
# retracted path with every accepted step satisfying the amendment-3 Armijo condition (recounted from the trace:
# slope <G_W, W+ - W> < 0, F_W <= F_H + 1e-4 slope, F_W < F_H), the b block never raises F. Reported: convergence
# fraction, iterations, residuals, rejected b steps, HVPs, W exits, time.
# Amendment 3 (Q-T table): "accepted W steps <= 2" -> "accepted t <= 1 and Armijo"; maxit 5000 -> per-start budget
# (1,000 sweeps, 25,000 HVPs; MAXIT is the amendment-3 sweep budget, not the old PG iteration count).
# A per-start wall-clock budget (test harness only) ends a start as "incomplete"; such starts are counted.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
NC <- as.integer(Sys.getenv("T_NCORES", "16")); BUDGET <- as.numeric(Sys.getenv("T4_BUDGET", "1800"))
corners <- rbind(expand.grid(alpha = c(1e-4, 0.9), rho = c(1e-2, 1e3), lnu = 1:2), data.frame(alpha = 0.3, rho = 1, lnu = 3))
LN <- list(c(1e-3, 0), c(1, 0.95), c(1, 0.5))
corners$lambda <- sapply(corners$lnu, function(i) LN[[i]][1]); corners$nu <- sapply(corners$lnu, function(i) LN[[i]][2])
NDS <- 3L; NST <- 5L; MAXIT <- 1000L
if (QUICK) { NDS <- 1L; NST <- 2L; MAXIT <- 150L; corners <- corners[c(4, 9), ]; BUDGET <- 120 }
dsets <- lapply(seq_len(NDS), function(i) { s <- sim_data(400 + i, 150, 60); list(s = s, D = palm_data(s$X, s$y, s$d)) })
jobs <- expand.grid(start = seq_len(NST), corner = seq_len(nrow(corners)), ds = seq_len(NDS))
epsv <- list()
for (di in seq_len(NDS)) for (ci in seq_len(nrow(corners)))
  epsv[[paste(di, ci)]] <- palm_eps(lapply(seq_len(NST), function(st) proj_S(palm_draw(dsets[[di]]$D, 3L, 7L, st))), dsets[[di]]$D)
one <- function(j) {
  di <- jobs$ds[j]; ci <- jobs$corner[j]; st <- jobs$start[j]; D <- dsets[[di]]$D
  hp <- list(alpha = corners$alpha[ci], rho = corners$rho[ci], lambda = corners$lambda[ci], nu = corners$nu[ci],
             lambdaH = 1e-6, eps = epsv[[paste(di, ci)]])
  W <- proj_S(palm_draw(D, 3L, 7L, st)); H <- palm_init_H(W, D, hp)
  f <- palm_fit_one(D, hp, W, H, rep(0, 3), maxit = MAXIT, tol = 1e-5, trace = TRUE, max_secs = BUDGET)
  tr <- f$trace; tol13 <- function(Fo) 1e-13 * pmax(1, abs(Fo))
  rc <- sum(tr[, "FH"] - tr[, "Fstart"] > tol13(tr[, "Fstart"]), na.rm = TRUE) +
        sum(tr[, "FW"] - tr[, "FH"] > tol13(tr[, "FH"]), na.rm = TRUE) + sum(tr[, "Fb"] - tr[, "FW"] > tol13(tr[, "FW"]), na.rm = TRUE)
  rok <- tr[!is.na(tr[, "r"]), , drop = FALSE]
  acc <- !is.na(tr[, "w_exit"]) & tr[, "w_exit"] == 1                                   # accepted Newton-CG steps
  arm <- !is.na(tr[, "slope"]) & tr[, "slope"] < 0 & tr[, "FW"] <= tr[, "FH"] + 1e-4 * tr[, "slope"] & tr[, "FW"] < tr[, "FH"]
  data.frame(ds = di, corner = ci, start = st, alpha = hp$alpha, rho = hp$rho, lambda = hp$lambda, nu = hp$nu, eps = hp$eps,
             status = f$status, reason = f$reason, F = f$F, iter = f$iter, r = f$r[["r"]], rW = f$r[["rW"]], rH = f$r[["rH"]],
             rb = f$r[["rb"]], r_original = f$r[["r_original"]], inc_counter = sum(f$ninc), inc_recount = rc,
             maxinc = max(f$maxinc), b_rejected = f$b_rejected, w_null = f$w_null_steps, w_backtracks = f$w_backtracks,
             b_iters = f$b_iters, max_bit = max(tr[, "bit"]), t_max = suppressWarnings(max(tr[, "t"], na.rm = TRUE)),
             b_raised = sum(tr[, "Fb"] > tr[, "FW"], na.rm = TRUE), hvp_total = f$hvp_total,
             w_newton = sum(acc), armijo_ok = sum(arm[acc]), w_noop = f$w_exits[["w_converged"]], w_ls_failed = f$w_exits[["line_search_failed"]],
             cg_curv = f$cg_exits[["curvature"]], cg_ceiling = f$cg_exits[["hvp_ceiling"]], nc1 = f$n_nc_first_iter,
             resid_recorded = all(is.finite(rok[, c("rW", "rH", "rb", "r_original")])) && all(rok[, "r"] == pmax(rok[, "rW"], rok[, "rH"], rok[, "rb"])),
             r_hit_1e3 = if (any(rok[, "r"] < 1e-3)) min(rok[rok[, "r"] < 1e-3, "iter"]) else NA,
             r_hit_1e4 = if (any(rok[, "r"] < 1e-4)) min(rok[rok[, "r"] < 1e-4, "iter"]) else NA,
             secs = f$secs, stringsAsFactors = FALSE)
}
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(jobs)), function(j) tryCatch(one(j), error = function(e) conditionMessage(e)), mc.cores = NC, mc.preschedule = FALSE)
bad <- !vapply(res, is.data.frame, TRUE)
x <- do.call(rbind, res[!bad])
saveRDS(x, file.path(OUT, "T4_monotone.rds"))
open_out("T4_monotone.txt")
cat("T4 monotone F |", format(Sys.time()), "| wall", format(Sys.time() - t0), "| cores", NC, "| per-start budget", BUDGET, "s\n")
cat(sprintf("datasets: p=150 n=60, events %s, tied times %s; starts per fit %d; per-start budget %d sweeps, 25000 HVPs (w_solver ncg)\n", paste(sapply(dsets, function(z) z$D$ne), collapse = ","),
            paste(sapply(dsets, function(z) z$s$n_tied), collapse = ","), NST, MAXIT))
if (any(bad)) { cat("ERRORS in", sum(bad), "runs:\n"); print(unique(unlist(res[bad]))) }
nit <- sum(x$iter)
cat(sprintf("runs: %d; outer iterations %d; block steps checked %d (H %d, W %d, b %d)\n", nrow(x), nit, 3 * nit, nit, nit, nit))
cat(sprintf("F increases: counter %d, independent recount from trace %d; largest raw block change %.3e\n", sum(x$inc_counter), sum(x$inc_recount), max(x$maxinc)))
cat(sprintf("status: converged %d, nonconverged %d, failure %d, incomplete (time budget) %d\n", sum(x$status == "converged"),
            sum(x$status == "nonconverged"), sum(x$status == "failure"), sum(x$status == "incomplete")))
if (any(x$status == "failure")) print(table(x$reason[x$status == "failure"]))
cat(sprintf("reached r < 1e-5: %d of %d (%.1f%%); r < 1e-4 at some iteration: %d; r < 1e-3: %d\n", sum(x$status == "converged"), nrow(x),
            100 * mean(x$status == "converged"), sum(!is.na(x$r_hit_1e4)), sum(!is.na(x$r_hit_1e3))))
cat(sprintf("residuals recorded every iteration with r = max(r_W, r_H, r_b): %d of %d runs\n", sum(x$resid_recorded), nrow(x)))
cat(sprintf("accepted W steps: %d; max accepted t %.3g (must be <= 1); satisfying the Armijo recount %d of %d; b blocks that raised F: %d; rejected b steps: %d\n",
            sum(x$w_newton), max(x$t_max), sum(x$armijo_ok), sum(x$w_newton), sum(x$b_raised), sum(x$b_rejected)))
cat(sprintf("W exits: Newton %d, no-op (r_W < 1e-5) %d, line-search failed %d; CG curvature exits %d (nc_first_iter %d), HVP-ceiling exits %d; HVPs per run median %g, max %g\n",
            sum(x$w_newton), sum(x$w_noop), sum(x$w_ls_failed), sum(x$cg_curv), sum(x$nc1), sum(x$cg_ceiling), median(x$hvp_total), max(x$hvp_total)))
if (any(x$status == "nonconverged")) print(table(nonconverged_reason = x$reason[x$status == "nonconverged"]))
cat(sprintf("eps shared within each fit: %d fits, %d with a single eps value\n", nrow(unique(x[, c("ds", "corner")])),
            sum(tapply(x$eps, paste(x$ds, x$corner), function(e) length(unique(e)) == 1))))
cat(sprintf("iterations: median %g, max %g; secs per run median %.1f, max %.1f; FISTA iterations per b block max %d\n",
            median(x$iter), max(x$iter), median(x$secs), max(x$secs), max(x$max_bit)))
cat("\nby corner (all datasets and starts): medians, and counts of statuses\n")
agg <- aggregate(cbind(iter, r, rW, rH, rb, r_original, b_iters, secs) ~ corner + alpha + rho + lambda + nu, x, median, na.rm = TRUE, na.action = na.pass)
agg$converged <- tapply(x$status == "converged", x$corner, sum)[as.character(agg$corner)]
agg$failure <- tapply(x$status == "failure", x$corner, sum)[as.character(agg$corner)]
agg$incomplete <- tapply(x$status == "incomplete", x$corner, sum)[as.character(agg$corner)]
print(agg, digits = 3)
cat("\nper run:\n"); print(x[, c("ds", "corner", "start", "status", "F", "iter", "r", "rW", "rH", "rb", "r_original", "inc_counter", "inc_recount", "b_rejected", "t_max", "w_newton", "armijo_ok", "hvp_total", "max_bit", "secs")], digits = 4)
pass <- !any(bad) && sum(x$inc_counter) == 0 && sum(x$inc_recount) == 0 && !any(x$status == "failure") && all(x$resid_recorded) &&
  max(x$t_max) <= 1 && all(x$armijo_ok == x$w_newton) && sum(x$b_raised) == 0
cat(sprintf("\nT4 TOTAL: %d runs; checks per run: increases (counter, recount), status not failure, residuals recorded, accepted t <= 1 and Armijo, b never raises F\n", nrow(x)))
verdict(pass)
sink()
