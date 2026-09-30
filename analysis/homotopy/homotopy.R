# Homotopy experiment, prespecified in docs/plans/2026-09-30-homotopy-prespec.md (commit 20abbe4).
# Existing desurv_alpha_warmstart, 10 paths (seeds 5,000,001..5,000,010), alpha grid seq(0, alpha*, 8).
suppressPackageStartupMessages(library(DeSurv))
stopifnot(packageVersion("DeSurv") == "1.1.0", exists("cox_stats_cpp", asNamespace("DeSurv")))
DATA <- Sys.getenv("H_DATA_DIR"); OUT <- Sys.getenv("H_OUT_DIR"); dir.create(OUT, FALSE, TRUE)
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
p  <- readRDS(file.path(DATA, "tar_params_best_tcgacptac.rds"))
dat <- readRDS(file.path(DATA, "tar_data_filtered_tcgacptac.rds"))
X <- dat$ex; y <- dat$sampInfo$time; d <- dat$sampInfo$event
grid <- seq(0, p$alpha, length.out = 8)
t0 <- Sys.time()
ws <- desurv_alpha_warmstart(X = X, y = y, d = d, k = p$k, alpha_grid = grid,
                             lambda = p$lambda, nu = p$nu, lambdaW = 0, lambdaH = 0,
                             n_starts = 10, seed = 5000001L, tol = 1e-7, maxit = 4000,
                             parallel = TRUE, ncores = min(NC, 10L), verbose = FALSE)
saveRDS(ws, file.path(OUT, "homotopy_warmstart.rds"))
nrm <- function(W) sweep(W, 2, sqrt(colSums(W^2)), "/")
Wat <- function(i) lapply(ws$fits, function(path) nrm(path[[i]]$W))
pairwise <- function(Ws) {
  k <- ncol(Ws[[1]]); cm <- utils::combn(length(Ws), 2); out <- matrix(NA, ncol(cm), k)
  for (c in seq_len(ncol(cm))) {
    A <- Ws[[cm[1, c]]]; B <- Ws[[cm[2, c]]]
    Wb <- B[, DeSurv:::.dc_align(A, B), drop = FALSE]
    out[c, ] <- diag(cor(A, Wb))
  }
  out
}
stable <- function(M) all(apply(M, 2, median) >= 0.98) && all(apply(M, 2, min) >= 0.95)
r0 <- pairwise(Wat(1)); rS <- pairwise(Wat(length(grid)))
s0 <- stable(r0); sS <- stable(rS)
decision <- if (!s0) "STOP (alpha = 0 not stable)" else if (sS) "PROCEED (alpha = 0 and alpha* stable)" else "STOP (alpha* not stable)"
ref <- readRDS(file.path(DATA, "tar_fit_desurv_tcgacptac.rds"))
Wref <- nrm(ref$W)
toref <- t(vapply(Wat(length(grid)), function(W) diag(cor(Wref, W[, DeSurv:::.dc_align(Wref, W), drop = FALSE])), numeric(p$k)))
fmt <- function(M) paste(sprintf("f%d median %.4f min %.4f", seq_len(ncol(M)), apply(M, 2, median), apply(M, 2, min)), collapse = " | ")
L <- c(sprintf("Homotopy experiment (prespec 20abbe4). DeSurv %s. %s", packageVersion("DeSurv"), format(Sys.time())),
       sprintf("alpha grid: %s", paste(round(grid, 4), collapse = ", ")),
       sprintf("DECISION: %s", decision),
       sprintf("alpha = 0 (45 pairs): %s -> %s", fmt(r0), if (s0) "stable" else "NOT stable"),
       sprintf("alpha* (45 pairs):   %s -> %s", fmt(rS), if (sS) "stable" else "NOT stable"),
       "[descriptive] alpha* endpoints vs reported fit W (per path, per factor):",
       apply(round(toref, 3), 1, paste, collapse = " "),
       sprintf("[descriptive] training C at alpha*: %s", paste(round(vapply(ws$fits, function(pth) pth[[length(grid)]]$cindex, 0), 4), collapse = " ")),
       sprintf("Elapsed %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
writeLines(L, file.path(OUT, "homotopy_summary.txt")); cat(L, sep = "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
