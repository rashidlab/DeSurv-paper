# D-Q3b (docs/plans/2026-10-05-dq3-collapse-trajectory-prespec.md): one-dimensional collapse path from each D-Q1 end
# point toward the uniform-vector limit, holding the score direction fixed. w_j(s) = (c_j u + s d_j)/||.|| for
# j in J = {b_j != 0}. Variant A: H fixed, b re-solved. Variant B: H profiled by HALS, then b re-solved. Descriptive.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))

hp <- readRDS(file.path(CAND, "tests", "out", "superseded", "pre-amendment3", "T7_fit.rds"))$fit$hp
dq1 <- readRDS(file.path(OUT, "DQ1_block_replay.rds"))
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
SS <- c(1L, 5L, 16L); S_GRID <- c(1, 0.9, 0.75, 0.5, 0.25, 0.1, 0.03, 0.01, 0.003, 0.001, 3e-4, 1e-4)
u <- rep(1, nrow(D$X)) / sqrt(nrow(D$X))
xu <- as.vector(crossprod(D$X, u))

profile_H <- function(W, H) {
  for (it in seq_len(5000L)) {
    Hn <- hals_sweep(W, H, D, hp)
    ch <- max(abs(Hn - H) / pmax(abs(H), 1e-300)); H <- Hn
    if (ch < 1e-12) return(list(H = H, sweeps = it, capped = FALSE))
  }
  list(H = H, sweeps = 5000L, capped = TRUE)
}

one <- function(i) {
  e <- dq1$res[[match(i, SS)]]$f; W0 <- e$W; J <- which(e$b != 0)
  cj <- as.vector(crossprod(W0, u)); Dv <- W0 - outer(u, cj)
  eta_of <- function(W, b) { sc <- palm_scores(W, D, hp$eps); list(eta = as.vector(sc$Zc %*% (b / sc$sigma)), Zs = sweep(sc$Zc, 2, sc$sigma, "/"), sigma = sc$sigma) }
  rows <- list(); ref <- list()
  for (v in c("A", "B")) {
    b <- e$b; H <- e$H
    for (sv in S_GRID) {
      W <- W0
      for (j in J) { w <- cj[j] * u + sv * Dv[, j]; W[, j] <- w / sqrt(sum(w^2)) }
      feas <- all(W >= 0); hs <- NA_integer_; hcap <- NA
      if (v == "B") { ph <- profile_H(W, H); H <- ph$H; hs <- ph$sweeps; hcap <- ph$capped }
      sb <- solve_b(b, W, D, hp); b <- sb$b
      pp <- palm_parts(W, H, b, D, hp); pf <- palm_parts(W, H, e$b, D, hp)   # pf: survival parts at fixed end-point b
      et <- eta_of(W, b)
      if (sv == 1) ref[[v]] <- et
      rr <- palm_resid(W, H, b, D, hp)
      rows[[length(rows) + 1L]] <- data.frame(start = i, variant = v, s = sv, feasible = feas, H_sweeps = hs, H_capped = hcap,
        F = pp$F, nmf = pp$nmf, hreg = pp$hreg, surv = pp$surv, l1 = pp$l1, l2 = pp$l2, ridge = pp$ridge,
        surv_fixb = pf$surv, ridge_fixb = pf$ridge,
        sigma1 = et$sigma[1], sigma2 = et$sigma[2], sigma3 = et$sigma[3],
        beta1 = pp$beta[1], beta2 = pp$beta[2], beta3 = pp$beta[3], b1 = b[1], b2 = b[2], b3 = b[3],
        cor_eta = cor(et$eta, ref[[v]]$eta), min_cor_Zs = min(diag(cor(et$Zs, ref[[v]]$Zs))[J]), rW = rr[["rW"]])
    }
  }
  do.call(rbind, rows)
}
res <- do.call(rbind, lapply(SS, one))
saveRDS(list(res = res, sd_Xu = sd(xu), S_GRID = S_GRID), file.path(OUT, "DQ3b_collapse_path.rds"))

sink(file.path(OUT, "DQ3b_collapse_path.txt"), split = TRUE)
cat(sprintf("D-Q3b collapse path (descriptive). Premise check: sd over samples of X'u = %.3e (0 means the uniform direction gives a constant score).\n", sd(xu)))
op <- options(width = 250)
for (i in SS) for (v in c("A", "B")) {
  d <- res[res$start == i & res$variant == v, ]
  cat(sprintf("\n=== start %d, variant %s ===\n", i, v))
  print(data.frame(s = d$s, F = sprintf("%.10f", d$F), dF = sprintf("%+.3e", d$F - d$F[1]),
    nmf = sprintf("%.8f", d$nmf), surv = sprintf("%.8f", d$surv), l1 = sprintf("%.6f", d$l1), l2 = sprintf("%.6f", d$l2),
    ridge = sprintf("%.6f", d$ridge), d_surv_fixb = sprintf("%+.2e", d$surv_fixb - d$surv_fixb[1]),
    d_ridge_fixb = sprintf("%+.2e", d$ridge_fixb - d$ridge_fixb[1]),
    sigma = sprintf("%.3g/%.3g/%.3g", d$sigma1, d$sigma2, d$sigma3), beta = sprintf("%.3g/%.3g/%.3g", d$beta1, d$beta2, d$beta3),
    cor_eta = sprintf("%.8f", d$cor_eta), min_cor_Zs = sprintf("%.8f", d$min_cor_Zs), rW = sprintf("%.2e", d$rW),
    H_sw = d$H_sweeps, feas = d$feasible), row.names = FALSE)
}
options(op); sink()
