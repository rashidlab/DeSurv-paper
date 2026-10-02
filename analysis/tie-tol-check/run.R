# Diagnostic (2026-10-02): tie fix (DeSurv 1.0.1 vs 1.2.0) x convergence tolerance (1e-5 vs 1e-7) on single-start
# supervised fits (May-evaluator style: ninit = 1). Fixed 80/20 split of tar_data_filtered; train C from the fit,
# held-out C from top-ntop scoring (each factor's top ntop genes, W zeroed elsewhere; Z = X'W, lp = Z beta).
# Run once per library: R_LIBS=<lib> Rscript analysis/tie-tol-check/run.R <label>
suppressMessages(library(DeSurv)); lab <- commandArgs(TRUE)[1]
td <- readRDS("results/tar_data_filtered_tcgacptac.rds"); X <- td$ex; y <- td$sampInfo$time; d <- td$sampInfo$event
set.seed(20261002); n <- ncol(X); va <- sort(sample(n, round(0.2 * n))); tr <- setdiff(seq_len(n), va)
cfg <- list(weak_eval46 = c(k = 3, alpha = 0.29254, lambda = 0.03631, nu = 0.5643, ntop = 149),
            published   = c(k = 3, alpha = 0.334,   lambda = 0.349,   nu = 0.0557, ntop = 270))
cidx <- function(lp, t, e) survival::concordance(survival::Surv(t, e) ~ lp, reverse = TRUE)$concordance
rows <- list()
for (cn in names(cfg)) for (tl in c(1e-5, 1e-7)) for (s in 1:4) {
  p <- cfg[[cn]]
  f <- desurv_fit(X[, tr], y[tr], d[tr], k = p[["k"]], alpha = p[["alpha"]], lambda = p[["lambda"]], nu = p[["nu"]],
                  lambdaW = 0, lambdaH = 0, seed = s, tol = tl, tol_init = tl, maxit = 4000, imaxit = 4000,
                  ninit = 1, verbose = FALSE)
  W <- f$W; for (j in seq_len(ncol(W))) { o <- order(W[, j], decreasing = TRUE); W[-o[seq_len(p[["ntop"]])], j] <- 0 }
  lp_tr <- drop(t(X[, tr]) %*% W %*% f$beta); lp_va <- drop(t(X[, va]) %*% W %*% f$beta)
  rows[[length(rows) + 1]] <- data.frame(pkg = lab, version = as.character(packageVersion("DeSurv")), config = cn,
    tol = tl, seed = s, iters = length(f$lossit), fit_train_C = as.numeric(f$cindex)[1],
    top_train_C = cidx(lp_tr, y[tr], d[tr]), top_val_C = cidx(lp_va, y[va], d[va]))
}
res <- do.call(rbind, rows); write.csv(res, sprintf("analysis/tie-tol-check/res_%s.csv", lab), row.names = FALSE)
print(aggregate(cbind(iters, fit_train_C, top_train_C, top_val_C) ~ config + tol, res, function(x) round(mean(x), 3)))
