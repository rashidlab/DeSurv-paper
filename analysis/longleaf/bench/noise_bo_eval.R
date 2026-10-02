# Step 4 diagnostic (2026-10-02): seed-bank noise of the reported-estimator CV score; adapted from of the "reported" estimator at the full
# settings code/03 uses (3 CV partitions x 5 folds x NINIT_FULL = 100 seed fits per fold).
# Standalone: reads tar_data_tcgacptac.rds from BENCH_DIR, writes noise_results.rds there.
suppressMessages(library(DeSurv))
stopifnot(packageVersion("DeSurv") >= "1.2.0")
dir <- Sys.getenv("BENCH_DIR"); stopifnot(nzchar(dir)); ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "32"))
td <- readRDS(file.path(dir, "tar_data_tcgacptac.rds"))
# code/03 BO_COMMON / BO_FIXED / R/reported_settings.R (full mode)
common <- list(X = td$ex, y = td$sampInfo$time, d = td$sampInfo$event,
  dataset = td$sampInfo$dataset, samp_keeps = td$samp_keeps,
  preprocess = TRUE, ngene = 3000L, method_trans_train = "rank", nfolds = 5L,
  tol = 1e-7, maxit = 4000L, engine = "warmstart", estimator = "reported",
  cv_seeds = c(20261101L, 20261102L, 20261103L), n_seed_fits = 100L,
  reported_control = list(tol_init = 1e-5, imaxit = 4000L),
  lambdaW_grid = 0, lambdaH_grid = 0, parallel_grid = TRUE, ncores_grid = ncores,
  cv_only = TRUE, verbose = TRUE)
# current reported point (k=3, alpha 0.334, lambda 0.349, nu 0.0557, ntop 270) and the same
# point at k = 7 and at the BO upper bound k = 12 (cost grows with k)
cfgs <- data.frame(label = c("A_new_incumbent_k6","B_published_k3","C_oldBObest_k3"), k = c(6L, 3L, 3L), alpha = c(0.003361216, 0.334, 0.442), lambda = c(0.03079838, 0.349, 0.2410), nu = c(0.9694315, 0.0557, 0.2015), ntop = c(174L, 270L, 299L))
banks <- c(4000000L, 5000000L, 6000000L)   # 4,000,000 = the BO bank (reproduction check); others independent
cfgs <- merge(cfgs, data.frame(bank = banks)); cfgs <- cfgs[order(cfgs$label, cfgs$bank), ]
out <- vector("list", nrow(cfgs))
for (i in seq_len(nrow(cfgs))) {
  cf <- cfgs[i, ]
  t0 <- Sys.time()
  res <- try(do.call(desurv_cv, c(common, list(k_grid = cf$k, alpha_grid = cf$alpha,
    lambda_grid = cf$lambda, nu_grid = cf$nu, ntop = cf$ntop, seed_bank_base = cf$bank))))
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  out[[i]] <- list(cfg = cf, minutes = el, ok = !inherits(res, "try-error"),
                   res = if (inherits(res, "try-error")) as.character(res) else res)
  message(sprintf("NOISE %s bank=%d C=%.4f minutes=%.1f", cf$label, cf$bank, if (out[[i]]$ok) mean(res$cv_results$cindex, na.rm = TRUE) else NA, el)) #cf$k, el, out[[i]]$ok, ncores))
  saveRDS(list(out = out, ncores = ncores, node = Sys.info()[["nodename"]],
               desurv = as.character(packageVersion("DeSurv"))), file.path(dir, "noise_results.rds"))
}
