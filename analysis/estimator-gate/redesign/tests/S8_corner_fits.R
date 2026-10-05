# S8: full fits at the corners of the spec section 3 box, p = 80, n = 50 with tied times (one sim_data draw).
# alpha {1e-4, 0.3, 0.9} x lambda {1e-3, 1, 1e3} x rho {1e-2, 1e3} x nu {0, 0.95} = 36 configurations, 3 starts
# each (fit seed 801), maxit 5,000 sweeps, tol 1e-5, lambdaH 1e-6. Pass: zero F increases and zero numerical
# failures over all 108 starts. Convergence counts are REPORTED only. Configurations run in parallel (16 cores);
# each start is deterministic.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S8_corner_fits.txt")
cat("S8 corner fits |", format(Sys.time()), "\n")
s <- sim_data(21, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d s0=%.6f lambda_1=%.6e\n", D$p, D$n, D$ne, s$n_tied, D$s0, D$lam1))
G8 <- expand.grid(alpha = c(1e-4, 0.3, 0.9), lambda = c(1e-3, 1, 1e3), rho = c(1e-2, 1e3), nu = c(0, 0.95))
runs <- parallel::mclapply(seq_len(nrow(G8)), function(i) {
  hp <- c(as.list(G8[i, ]), lambdaH = 1e-6)
  palm_fit(D, hp, k = 3L, seed = 801L, nstarts = 3L, maxit = 5000L, keep_starts = TRUE)
}, mc.cores = 16L)
stopifnot(all(vapply(runs, is.list, TRUE)))
saveRDS(list(grid = G8, runs = runs), file.path(OUT, "S8_corner_fits.rds"))
rows <- lapply(seq_len(nrow(G8)), function(i) {
  st <- runs[[i]]$starts
  data.frame(cfg = i, G8[i, ], fit = runs[[i]]$status, conv = sum(st$status == "converged"),
             noncv = sum(st$status == "nonconverged"), fail = sum(st$status == "failure"),
             finc = sum(st$ninc, na.rm = TRUE), numfail = sum(st$status == "failure" & st$ninc == 0),
             th_rej = sum(st$theta_rejected), iter_max = max(st$iter), r_max = signif(max(st$r), 3),
             reasons = paste(unique(st$reason[nzchar(st$reason)]), collapse = "; "), row.names = NULL)
})
tb <- do.call(rbind, rows)
op <- options(width = 250); print(tb, row.names = FALSE); options(op)
cat(sprintf("\nstarts: %d; converged %d, nonconverged %d, failure %d\n", 3 * nrow(G8), sum(tb$conv), sum(tb$noncv), sum(tb$fail)))
cat(sprintf("fits: converged %d, nonconverged %d, failed %d (reported only)\n", sum(tb$fit == "converged"),
            sum(tb$fit == "nonconverged"), sum(tb$fit == "failed")))
cat(sprintf("theta steps rejected (F would increase): %d\n", sum(tb$th_rej)))
incs <- do.call(rbind, lapply(seq_along(runs), function(i) do.call(rbind, lapply(runs[[i]]$fits, function(f)
  if (!is.null(f$f_increases)) cbind(cfg = i, start = f$start, f$f_increases)))))
if (!is.null(incs)) { cat("\nF increases recorded:\n"); print(incs) }
cat("\n")
tally("F increases per start (108 starts)", unlist(lapply(runs, function(r) r$starts$ninc)), 0)
tally("numerical failures per start (1 = failure without F increase)",
      unlist(lapply(runs, function(r) as.numeric(r$starts$status == "failure" & r$starts$ninc == 0))), 0)
total("S8")
sink()
