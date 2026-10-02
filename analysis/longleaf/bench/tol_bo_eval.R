# Step 4 diagnostic (2026-10-02): does convergence tolerance cause the low-k overfitting under the new evaluator?
# One config (BO start point 46: k 3, alpha 0.2925, lambda 0.0363, nu 0.5643, ntop 149) scored 2 x 2:
# estimator {reported (100 seed fits -> consensus -> final fit), mean_starts (30 starts, May)} x tol {1e-5, 1e-7}.
# Same CV partitions and seed bank as the rev 3 BO. Reports val and train C per cell.
suppressMessages(library(DeSurv)); stopifnot(packageVersion("DeSurv") >= "1.2.0")
dir <- Sys.getenv("BENCH_DIR"); stopifnot(nzchar(dir)); ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "32"))
td <- readRDS(file.path(dir, "tar_data_tcgacptac.rds"))
base <- list(X = td$ex, y = td$sampInfo$time, d = td$sampInfo$event, dataset = td$sampInfo$dataset,
  samp_keeps = td$samp_keeps, preprocess = TRUE, ngene = 3000L, method_trans_train = "rank", nfolds = 5L,
  maxit = 4000L, engine = "warmstart", cv_seeds = c(20261101L, 20261102L, 20261103L),
  lambdaW_grid = 0, lambdaH_grid = 0, parallel_grid = TRUE, ncores_grid = ncores, cv_only = TRUE, verbose = FALSE,
  k_grid = 3L, alpha_grid = 0.29254, lambda_grid = 0.03631, nu_grid = 0.5643, ntop = 149L)
cells <- expand.grid(estimator = c("reported", "mean_starts"), tol = c(1e-5, 1e-7), stringsAsFactors = FALSE)
out <- list()
for (i in seq_len(nrow(cells))) {
  e <- cells$estimator[i]; tl <- cells$tol[i]
  a <- c(base, list(tol = tl, estimator = e))
  if (e == "reported") a <- c(a, list(n_seed_fits = 100L, seed_bank_base = 4000000L,
                                      reported_control = list(tol_init = 1e-5, imaxit = 4000L)))
  else a <- c(a, list(n_starts = 30L))
  t0 <- Sys.time(); res <- try(do.call(desurv_cv, a))
  d <- if (inherits(res, "try-error")) NULL else res$diagnostics
  v <- if (is.null(d)) NA else mean(d$val_cindex, na.rm = TRUE); tr <- if (is.null(d)) NA else mean(d$train_cindex, na.rm = TRUE)
  message(sprintf("TOL estimator=%s tol=%g val=%.4f train=%.4f minutes=%.1f %s", e, tl, v, tr,
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), if (inherits(res, "try-error")) as.character(res) else ""))
  out[[i]] <- list(cell = cells[i, ], res = res); saveRDS(out, file.path(dir, "tol_results.rds"))
}
