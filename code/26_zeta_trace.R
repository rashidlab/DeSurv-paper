#!/usr/bin/env Rscript
# code/26_zeta_trace.R
# ---------------------------------------------------------------------------
# Empirical behaviour of the gradient-scaling factor zeta in the W update, and
# the distance between the consensus initialization W0 and the final fitted W.
#
# READ-ONLY with respect to every existing cached object. Loads cached results
# via paper/load_precomputed.R and writes exactly one NEW file,
# results/zeta_trace.rds. Never calls cache_or_compute() and never sources
# code/00_helpers.R, so it cannot overwrite the objects that anchor the paper.
#
#   Run as:  DESURV_RECOMPUTE=FALSE Rscript code/26_zeta_trace.R
#
# Background (2026-09-06 audit). In src/functions.cpp of DeSurv v1.0.1 (commit
# afb00d5) the W update forms
#     g_N = (1 - alpha) * grad_W L_NMF,   g_C = alpha * grad_W L_Cox,
#     zeta = ||g_N||_F / ||g_C||_F  (plain Frobenius norms; capped at 1e6),
# and adds alpha * zeta * g_C to the multiplicative numerator. Away from the cap
# and the zero-gradient case the added term therefore has norm alpha * ||g_N||_F.
# The package does not log zeta, so this script reconstructs the production
# optimization path from the cached seed fits and the consensus rule of
# code/04, ASSERTS that the reconstruction agrees with the cached final fit to
# declared tolerances, and evaluates zeta after every iteration by re-running
# the deterministic fit with maxit = 1, 2, ... (same init, same settings).
#
# Outputs (results/zeta_trace.rds)
#   trace         data.frame: iter, zeta, gn, gc, capped, cindex (after each
#                 beta update along the production path)
#   final         zeta at the cached final iterate; norms of the scaled Cox term
#                 and of the NMF numerator term
#   init          the state at initialization: beta0 = 0, Cox gradient zero,
#                 ratio capped, added term identically zero
#   endpoints     zeta at each of the 100 single-run (random-start) endpoints
#   w0            per-program correlation and top-ntop Jaccard between W0 and
#                 the final W (columns correspond by position; the identity
#                 alignment is checked against the full correlation matrix)
#   reconstruction  agreement of the reconstructed path with the cached fit
#   meta
# ---------------------------------------------------------------------------
suppressMessages(library(DeSurv))
source("paper/load_precomputed.R")
stopifnot(identical(as.character(packageVersion("DeSurv")), "1.0.1"))

fit  <- read_result("tar_fit_desurv_tcgacptac")
sf   <- read_result("desurv_seed_fits_tcgacptac")
data <- read_result("tar_data_filtered_tcgacptac")
pb   <- read_result("tar_params_best_tcgacptac")

ntop     <- as.integer(round(pb$ntop))
min_freq <- max(1L, as.integer(ceiling(0.3 * length(sf$fits))))   # as in code/04
hy <- fit$hyper
alpha <- hy$alpha

## ---- 1. zeta as implemented -------------------------------------------------
# Breslow gradient of the partial log-likelihood, g_z, as in cox_suffix_stats_full.
zeta_at <- function(W, H, beta, X, y, d) {
  beta <- as.numeric(beta); y <- as.numeric(y); d <- as.numeric(d)
  lp <- as.numeric(t(X) %*% W %*% beta); e <- exp(lp - max(lp))
  S0  <- vapply(seq_along(y), function(i) sum(e[y >= y[i]]), numeric(1))
  cum <- vapply(seq_along(y), function(i) sum(d[y <= y[i]] / S0[y <= y[i]]), numeric(1))
  g <- d - e * cum
  gN <- ((1 - alpha) / sum(X^2)) * (W %*% H - X) %*% t(H)          # (1-alpha) grad_W L_NMF
  gC <- (2 * alpha / sum(d)) * (X %*% g) %*% t(beta)                # alpha grad_W L_Cox
  gn <- sqrt(sum(gN^2)) + 1e-12; gc <- sqrt(sum(gC^2)) + 1e-12
  zeta <- gn / gc
  c(gn = gn, gc = gc, zeta = zeta, capped = as.numeric(zeta >= 1e6),
    cox_term_norm = alpha * min(zeta, 1e6) * gc,
    nmf_num_norm = sqrt(sum((((1 - alpha) / sum(X^2)) * X %*% t(H))^2)))
}

## ---- 2. reconstruct the production path and assert agreement ----------------
init <- DeSurv::desurv_consensus_seed(fits = sf$fits, X = data$ex, ntop = ntop,
                                      k = pb$k, min_frequency = min_freq)
refit <- function(maxit) desurv_fit(
  X = data$ex, y = data$sampInfo$time, d = data$sampInfo$event,
  k = pb$k, alpha = pb$alpha, lambda = pb$lambda, nu = pb$nu,
  lambdaW = hy$lambdaW, lambdaH = hy$lambdaH,
  W0 = init$W0, H0 = init$H0, beta0 = init$beta0,
  seed = NULL, tol = hy$tol, maxit = maxit, verbose = FALSE)
full <- suppressWarnings(refit(hy$maxit))
TOL <- list(max_abs_dW = 1e-3, max_abs_dbeta = 1e-5, abs_dcindex = 1e-3)
recon <- list(
  n_iter_reconstructed = length(full$lossit), n_iter_cached = length(fit$lossit),
  max_abs_dW = max(abs(full$W - fit$W)),
  max_abs_dbeta = max(abs(as.numeric(full$beta) - as.numeric(fit$beta))),
  abs_dcindex = abs(full$cindex - fit$cindex), tolerances = TOL)
recon$pass <- with(recon, n_iter_reconstructed == n_iter_cached &&
                     max_abs_dW <= TOL$max_abs_dW && max_abs_dbeta <= TOL$max_abs_dbeta &&
                     abs_dcindex <= TOL$abs_dcindex)
if (!recon$pass) stop("Reconstructed production path does not match the cached fit: ",
                      paste(names(recon)[1:5], unlist(recon[1:5]), sep = "=", collapse = ", "))

## ---- 3. zeta along the production path ---------------------------------------
n_it <- recon$n_iter_cached
trace <- do.call(rbind, lapply(seq_len(n_it), function(m) {
  f <- suppressWarnings(refit(m))
  z <- zeta_at(f$W, f$H, f$beta, f$data$X, f$data$y, f$data$d)
  data.frame(iter = m, zeta = z[["zeta"]], gn = z[["gn"]], gc = z[["gc"]],
             capped = z[["capped"]] == 1, cindex = f$cindex)
}))
z_init  <- zeta_at(init$W0, init$H0, init$beta0, data$ex, data$sampInfo$time, data$sampInfo$event)
z_final <- zeta_at(fit$W, fit$H, fit$beta, fit$data$X, fit$data$y, fit$data$d)

## ---- 4. single-run endpoints ---------------------------------------------------
ep <- t(sapply(sf$fits, function(s) zeta_at(s$W, s$H, s$beta, s$data$X, s$data$y, s$data$d)))

## ---- 5. W0 -> W ----------------------------------------------------------------
cm <- cor(fit$W, init$W0)                                # rows: final W; cols: W0
aligned_identity <- all(apply(cm, 1, which.max) == seq_len(nrow(cm)))
tg_W  <- desurv_get_top_genes(fit$W,  ntop)$top_genes
tg_W0 <- desurv_get_top_genes(init$W0, ntop)$top_genes
jacc <- vapply(seq_len(pb$k), function(j) {
  a <- tg_W[[j]]; b <- tg_W0[[j]]; length(intersect(a, b)) / length(union(a, b)) }, numeric(1))

res <- list(
  trace = trace,
  final = list(zeta = z_final[["zeta"]], cox_term_norm = z_final[["cox_term_norm"]],
               nmf_num_norm = z_final[["nmf_num_norm"]],
               cox_to_nmf_ratio = z_final[["cox_term_norm"]] / z_final[["nmf_num_norm"]]),
  init = list(beta0_all_zero = all(init$beta0 == 0), cox_grad_norm = z_init[["gc"]] - 1e-12,
              capped = z_init[["capped"]] == 1, cox_term_norm = z_init[["cox_term_norm"]]),
  endpoints = list(zeta = unname(ep[, "zeta"]), n = nrow(ep), n_capped = sum(ep[, "capped"]),
                   min = min(ep[, "zeta"]), median = median(ep[, "zeta"]), max = max(ep[, "zeta"])),
  w0 = list(r = unname(diag(cm)), jaccard = jacc, ntop = ntop, aligned_identity = aligned_identity,
            cor_matrix = cm),
  reconstruction = recon,
  meta = list(package = "DeSurv 1.0.1 (afb00d5)", alpha = alpha, zeta_cap = 1e6,
              definition = "zeta = ||(1-alpha) grad_W L_NMF||_F / ||alpha grad_W L_Cox||_F; the W-update numerator adds alpha * min(zeta, 1e6) * alpha grad_W L_Cox",
              trace_convention = "zeta evaluated at (W, H, beta) after each completed iteration, i.e. after that iteration's beta update; iteration 0 is the consensus initialization with beta0 = 0",
              consensus = list(ntop = ntop, min_frequency = min_freq)))
saveRDS(res, "results/zeta_trace.rds")

cat(sprintf("zeta: iter1 %.3f -> final %.4f; capped after init: %d of %d iterations; init capped: %s (Cox term norm %.3g)\n",
            trace$zeta[1], z_final[["zeta"]], sum(trace$capped), n_it, res$init$capped, res$init$cox_term_norm))
cat(sprintf("endpoints: n=%d min %.4f median %.4f max %.4f capped %d\n", nrow(ep), min(ep[,"zeta"]), median(ep[,"zeta"]), max(ep[,"zeta"]), sum(ep[,"capped"])))
cat(sprintf("W0->W: r = %s; Jaccard(top-%d) = %s; identity alignment: %s\n", paste(round(diag(cm),3), collapse=" "), ntop, paste(round(jacc,3), collapse=" "), aligned_identity))
cat(sprintf("reconstruction: iters %d/%d, max|dW| %.2e, max|dbeta| %.2e, |dC| %.2e -> PASS\n", recon$n_iter_reconstructed, recon$n_iter_cached, recon$max_abs_dW, recon$max_abs_dbeta, recon$abs_dcindex))
