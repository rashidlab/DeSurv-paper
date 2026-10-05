# S9: end to end at p = 300. One sim_data draw (seed 29, n = 320) split by column: train = columns 1..120,
# test = 121..320 (within-sample ranks are per column, so the split ranks equal palm_rank of each part). Centre
# hp alpha 0.05, lambda 1, nu 0.5, rho 1, lambdaH 1e-6, ntop 50; 5 starts, fit seed 901, maxit 5,000, recording on.
# Checks: the fit runs and returns a start (status not "failed"); the model scores the test set with finite scores;
# a rerun with the same seed is bit-identical (fit and model, wall-clock secs fields excluded); recording present
# (column statistics every 10th sweep, checkpoints 10, 50, 100, 500, 3750 and the endpoint). Statuses and the
# held-out C are reported, descriptive only.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S9_end_to_end.txt")
cat("S9 end to end |", format(Sys.time()), "\n")
s <- sim_data(29, 300, 320); train <- 1:120; test <- 121:320
D <- palm_data(s$X[, train], s$y[train], s$d[train])
cat(sprintf("train: p=%d n=%d events=%d tied times=%d s0=%.6f lambda_1(XCX^T)=%.6e | test n=%d events=%d\n",
            D$p, D$n, D$ne, sum(duplicated(D$y)), D$s0, D$lam1, length(test), sum(s$d[test])))
hp <- list(alpha = 0.05, lambda = 1, nu = 0.5, rho = 1, lambdaH = 1e-6)
run <- function() {
  f <- palm_fit(D, hp, k = 3L, seed = 901L, nstarts = 5L, maxit = 5000L, ncores = 5L, keep_starts = TRUE, record = TRUE)
  m <- palm_model(f, D, ntop = 50L)
  list(fit = f, model = m)
}
a <- run(); b <- run()
strip <- function(r) {
  r$fit$starts$secs <- NULL; if (!is.null(r$fit$best)) r$fit$best$secs <- NULL
  r$fit$fits <- lapply(r$fit$fits, function(x) { x$secs <- NULL; x }); r
}
saveRDS(a, file.path(OUT, "S9_end_to_end.rds"))
op <- options(width = 250); print(a$fit$starts[, setdiff(names(a$fit$starts), "secs")], row.names = FALSE); options(op)
cat(sprintf("\nfit status %s, returned start %s, final F %.10e\n", a$fit$status, a$fit$best_start, a$fit$best$F))
cat(sprintf("theta = (%s); beta = theta/s0 = (%s)\n", paste(sprintf("%.6g", a$model$theta), collapse = ", "),
            paste(sprintf("%.6g", a$model$beta), collapse = ", ")))
sc <- palm_score(s$X[, test], a$model)
cc <- harrell(sc, s$y[test], s$d[test]); rl <- cor(sc, s$linpred[test])
cat(sprintf("support |G| = %d genes; held-out C = %.4f; cor(score, true linear predictor) = %.4f (descriptive only)\n",
            length(a$model$G), cc, rl))
rec <- a$fit$best$record
cat(sprintf("recording: %d column-statistic rows (every 10th sweep of %d sweeps); checkpoints stored: %s\n",
            if (is.null(rec$colstats)) 0L else nrow(rec$colstats), a$fit$best$iter, paste(names(rec$checkpoints), collapse = ", ")))
cat("\n")
tally("fit returned a start (1 = status failed)", as.numeric(a$fit$status == "failed"), 0)
tally("test scores finite (1 = any nonfinite)", as.numeric(any(!is.finite(sc))), 0)
tally("rerun bit-identical (1 = differs)", as.numeric(!identical(strip(a), strip(b))), 0)
tally("recording complete (1 = missing checkpoint or column rows)",
      as.numeric(!all(c("10", "50", "100", "500", "3750", "end") %in% names(rec$checkpoints)) ||
                 (a$fit$best$iter >= 10 && (is.null(rec$colstats) || nrow(rec$colstats) != a$fit$best$iter %/% 10))), 0)
cat(sprintf("  F increases over the 5 starts (reported only): %d\n", sum(a$fit$starts$ninc)))
total("S9")
sink()
