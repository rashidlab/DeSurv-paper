# May (published) BO reference values for the k = 3 ladder: per-fold held-out and training C of
# May evals 136 and 129 (mean over the 30 starts, as desurv_cv's summary_fold), from the cached BO
# diagnostics. Checks that ladder_eval.R's hard-coded configurations equal the cached history to
# full precision and that the per-fold means reproduce the history score.
# Read-only on results/. Usage (repo root): Rscript analysis/k3-ladder/may_reference.R
b <- readRDS("results/desurv_bo_results_tcgacptac.rds")
run <- b$runs[[1L]]; h <- run$history; dg <- run$diagnostics
stopifnot(identical(run$seed, 123L), run$fixed$n_starts == 30, run$fixed$tol == 1e-5, run$fixed$maxit == 4000,
          run$fixed$ngene == 3000, run$fixed$nfolds == 5, identical(run$fixed$engine, "warmstart"))
expect <- list(`136` = c(k = 3, alpha = 0.442207523710734685, lambda = 0.240956814604169955,
                         nu = 0.201480919147608795, ntop = 299),
               `129` = c(k = 7, alpha = 0.333935556177224480, lambda = 0.349448737941329357,
                         nu = 0.055651216804573779, ntop = 270))
rows <- list()
for (e in names(expect)) {
  he <- h[h$eval_id == as.integer(e), ]
  got <- c(k = he$k_grid, alpha = he$alpha_grid, lambda = he$lambda_grid, nu = he$nu_grid, ntop = he$ntop)
  if (!identical(unname(got), unname(expect[[e]]))) stop("eval ", e, " params differ from ladder_eval.R")
  de <- dg[dg$eval_id == as.integer(e), ]
  stopifnot(nrow(de) == 150L, all(table(de$fold) == 30L))
  f <- stats::aggregate(cbind(val_C = val_cindex, train_C = train_cindex) ~ fold, data = de, FUN = mean)
  stopifnot(abs(mean(f$val_C) - he$mean_cindex) < 1e-12)
  rows[[e]] <- data.frame(may_eval = as.integer(e), config = if (e == "136") "c136_k3" else "pub_k7",
                          rep = 1L, fold = f$fold, val_C = f$val_C, train_C = f$train_C,
                          history_mean_cindex = he$mean_cindex)
}
ref <- do.call(rbind, rows)
utils::write.csv(ref, "analysis/k3-ladder/may_reference.csv", row.names = FALSE)
print(ref, digits = 6)
