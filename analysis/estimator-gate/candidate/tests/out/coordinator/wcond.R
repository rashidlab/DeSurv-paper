# Coordinator diagnostic (PI ruling 2026-10-05, step 1): is the T7 stall a W-block conditioning problem?
# Simulated data only. Operator: Riemannian Hessian of Psi in W (H, b fixed) on the active face of S_+^k,
# i.e. tangent to each unit sphere and supported on {W > 0}: Hess[V]_j = P_j((d2Psi V)_j) - (w_j' g_j) V_j.
# Euclidean Hessian-vector products by central differences of palm_grad_W (h = 1e-5 / ||V||).
# Lanczos with full reorthogonalization: run A from a random face vector (extreme eigenvalues), run B from the
# support part of the stalled residual R_W (Gauss-quadrature weights = spectral measure of R_W).
#
# Decision rule, fixed before the first run. SUPPORTED only if all five hold at the p = 3000 end points
# (majority of the four starts for S1-S4; all four for S5):
#   S1  W-only projected-gradient steps at fixed (H, b) (same w_step, 3000 steps) do not bring r_W below 1e-5:
#       the W block alone reproduces the stall (otherwise the stall is block coupling, not W conditioning).
#   S2  kappa = lambda_max / lambda_min (Ritz; a lower bound on the true ratio) >= 1e3 at p = 3000 and >= 10x
#       the kappa of the converged p = 150 control.
#   S3  >= 50% of the spectral weight of R_W (support part) lies on Ritz values <= 1e-2 lambda_max.
#   S4  the boundary part of R_W (W = 0 entries with negative gradient) carries < 50% of ||R_W||^2
#       (otherwise: active-set identification, not conditioning).
#   S5  lambda_min >= -1e-3 lambda_max (otherwise: negative curvature / saddle, not conditioning).
CAND <- "/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/estimator-gate/candidate"
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
for (f in c("kkt_oracle.R", "palm_lib.R", "sim.R")) source(file.path(CAND, f))
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
M_LANCZOS <- 300L; N_WONLY <- 3000L

face_proj <- function(V, W) {
  V[W <= 0] <- 0
  V - sweep(W, 2, colSums(W * V), "*")
}
make_op <- function(W, gradfun, mu) function(V) {
  h <- 1e-5 / sqrt(sum(V^2))
  E <- (gradfun(W + h * V) - gradfun(W - h * V)) / (2 * h)
  face_proj(E, W) - sweep(V, 2, mu, "*")
}
lanczos <- function(op, v0, m, shape) {
  N <- length(v0); Q <- matrix(0, N, m); a <- numeric(m); bt <- numeric(m)
  q <- v0 / sqrt(sum(v0^2)); qprev <- 0; bprev <- 0
  for (j in seq_len(m)) {
    Q[, j] <- q
    w <- as.vector(op(matrix(q, shape[1], shape[2]))) - bprev * qprev
    a[j] <- sum(w * q); w <- w - a[j] * q
    for (rep in 1:2) w <- w - Q[, 1:j, drop = FALSE] %*% crossprod(Q[, 1:j, drop = FALSE], w)
    bt[j] <- sqrt(sum(w^2)); if (bt[j] < 1e-12) { m <- j; break }
    qprev <- q; bprev <- bt[j]; q <- as.vector(w) / bt[j]
  }
  Tm <- diag(a[1:m], m); if (m > 1) { Tm[cbind(1:(m - 1), 2:m)] <- bt[1:(m - 1)]; Tm[cbind(2:m, 1:(m - 1))] <- bt[1:(m - 1)] }
  e <- eigen(Tm, symmetric = TRUE)
  list(theta = e$values, wt = e$vectors[1, ]^2, res = abs(bt[m] * e$vectors[m, ]), S = e$vectors, Q = Q[, 1:m, drop = FALSE])
}

analyse <- function(lab, D, hp, W, H, b, t_last) {
  set.seed(1)
  gtot <- function(Wx) palm_grad_W(Wx, H, b, D, hp)
  grec <- function(Wx) (1 - hp$alpha) * grad_W_terms(Wx, H, b, D, hp)$recon
  gsur <- function(Wx) { g <- grad_W_terms(Wx, H, b, D, hp); hp$alpha * (g$cox + g$ridge) }
  G <- gtot(W); mu <- colSums(W * G)
  R <- ifelse(W > 0, G - sweep(W, 2, mu, "*"), pmin(G, 0))
  Rsup <- R; Rsup[W <= 0] <- 0; Rbd <- R - Rsup
  op <- make_op(W, gtot, mu)
  # FD accuracy check: same vector at h scaled 10x
  V <- face_proj(matrix(rnorm(length(W)), nrow(W)), W); V <- V / sqrt(sum(V^2))
  hv1 <- op(V); h2 <- function(V) { h <- 1e-4 / sqrt(sum(V^2)); face_proj((gtot(W + h * V) - gtot(W - h * V)) / (2 * h), W) - sweep(V, 2, mu, "*") }
  fd_rel <- sqrt(sum((hv1 - h2(V))^2)) / sqrt(sum(hv1^2))
  A <- lanczos(op, as.vector(V), M_LANCZOS, dim(W))
  B <- lanczos(op, as.vector(Rsup), M_LANCZOS, dim(W))
  lmax <- max(A$theta, B$theta); lmin <- min(A$theta, B$theta)
  # per-term curvature at the 3 smallest and 3 largest Ritz vectors of run A
  mu_r <- colSums(W * grec(W)); mu_s <- colSums(W * gsur(W))
  op_r <- make_op(W, grec, mu_r); op_s <- make_op(W, gsur, mu_s)
  ix <- c(order(A$theta)[1:3], order(-A$theta)[1:3])
  terms <- t(sapply(ix, function(i) {
    x <- matrix(A$Q %*% A$S[, i], nrow(W)); x <- x / sqrt(sum(x^2))
    c(theta = A$theta[i], ritz_res = A$res[i], recon = sum(x * op_r(x)), surv = sum(x * op_s(x)))
  }))
  # W-only continuation at fixed (H, b)
  Wc <- W; tp <- t_last; tr <- NULL; t0 <- proc.time()[["elapsed"]]
  for (s in seq_len(N_WONLY)) {
    ws <- w_step(Wc, H, b, D, hp, tp); if (ws$fail) break
    Wc <- ws$W; tp <- ws$t
    if (s %% 250 == 0 || s == 1) {
      Gc <- gtot(Wc); rw <- PALM_DEFS$resid_W(Wc, Gc, D, hp)
      tr <- rbind(tr, c(step = s, rW = rw, Psi = palm_Psi(Wc, H, b, D, hp), t = tp))
      if (rw < 1e-5) break
    }
  }
  list(lab = lab, p = D$p, n = D$n, support = colSums(W > 0), rW = sqrt(sum(R^2)),
       frac_boundary = sum(Rbd^2) / sum(R^2), fd_rel = fd_rel, t_last = t_last,
       lmax = lmax, lmin = lmin, kappa = lmax / lmin, A = A[c("theta", "wt", "res")],
       wt_slow = sum(B$wt[B$theta <= 1e-2 * lmax]), wt_slow3 = sum(B$wt[B$theta <= 1e-3 * lmax]),
       rq_R = sum(B$wt * B$theta), B = B[c("theta", "wt", "res")], terms = terms, wonly = tr,
       wonly_secs = proc.time()[["elapsed"]] - t0)
}

x <- readRDS(file.path(CAND, "tests/out/T7_fit.rds")); fit <- x$fit; hp <- fit$hp
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
ss <- c(1L, 5L, 9L, 16L)
# control: p = 150 (dataset 401, n = 60) at the same hyperparameters, where the frozen solver converges
s150 <- sim_data(401, 150L, 60L); D150 <- palm_data(s150$X, s150$y, s150$d)
hp150 <- hp[c("alpha", "lambda", "nu", "rho", "lambdaH")]
f150 <- palm_fit(D150, hp150, k = 3L, seed = 123L, nstarts = 5L, maxit = 5000L, ncores = 5L, keep_starts = TRUE)
c150 <- which(vapply(f150$fits, function(f) isTRUE(f$conv), TRUE))
cat("control fit: converged starts", paste(c150, collapse = ","), "iters", paste(f150$starts$iter[c150], collapse = ","), "\n")
hp150$eps <- f150$eps
jobs <- c(lapply(ss, function(i) list(lab = paste0("p3000_n218_start", i), D = D, hp = hp, g = fit$fits[[i]])),
          lapply(head(c150, 2), function(i) list(lab = paste0("p150_n60_start", i, "_converged"), D = D150, hp = hp150, g = f150$fits[[i]])))
res <- parallel::mclapply(jobs, function(j) analyse(j$lab, j$D, j$hp, j$g$W, j$g$H, j$g$b, j$g$t_last), mc.cores = length(jobs))
saveRDS(res, file.path(HERE, "wcond.rds"))

sink(file.path(HERE, "wcond.txt"), split = TRUE)
cat("W-block conditioning diagnostic |", format(Sys.time()), "| Lanczos steps", M_LANCZOS, "| W-only steps", N_WONLY, "\n")
cat(sprintf("hp: alpha=%.6g lambda=%g nu=%g rho=%.6g eps=%.4g\n", hp$alpha, hp$lambda, hp$nu, hp$rho, hp$eps))
for (r in res) {
  if (!is.list(r) || is.null(r$lab)) { cat("\nERROR:", as.character(r), "\n"); next }
  cat(sprintf("\n== %s (p %d, n %d) support per column %s\n", r$lab, r$p, r$n, paste(r$support, collapse = "/")))
  cat(sprintf("  ||R_W|| %.3e; boundary share of ||R_W||^2 %.3f; FD rel. diff (h 1e-5 vs 1e-4) %.2e\n", r$rW, r$frac_boundary, r$fd_rel))
  cat(sprintf("  lambda_max %.4e  lambda_min %.4e  kappa %.3e  (t_last %.4g, t*lambda_max %.3f)\n", r$lmax, r$lmin, r$kappa, r$t_last, r$t_last * r$lmax))
  cat(sprintf("  run A smallest 5 Ritz: %s\n  run A largest 3 Ritz: %s\n", paste(sprintf("%.3e", sort(r$A$theta)[1:5]), collapse = " "),
              paste(sprintf("%.3e", sort(r$A$theta, TRUE)[1:3]), collapse = " ")))
  cat(sprintf("  R_W spectral weight on theta <= 1e-2 lmax: %.3f; <= 1e-3 lmax: %.3f; Rayleigh quotient of R_W %.3e (= %.2e lmax)\n",
              r$wt_slow, r$wt_slow3, r$rq_R, r$rq_R / r$lmax))
  cat("  per-term curvature at Ritz vectors (3 smallest, 3 largest):\n"); print(signif(r$terms, 4))
  cat(sprintf("  W-only continuation at fixed (H, b), %.0f s:\n", r$wonly_secs)); print(signif(r$wonly, 4))
}
p3 <- Filter(function(r) is.list(r) && r$p == 3000, res); p1 <- Filter(function(r) is.list(r) && r$p == 150, res)
k150 <- max(vapply(p1, `[[`, 0, "kappa"))
S1 <- vapply(p3, function(r) min(r$wonly[, "rW"]) >= 1e-5, TRUE)
S2 <- vapply(p3, function(r) r$kappa >= 1e3 && r$kappa >= 10 * k150, TRUE)
S3 <- vapply(p3, function(r) r$wt_slow >= 0.5, TRUE)
S4 <- vapply(p3, function(r) r$frac_boundary < 0.5, TRUE)
S5 <- vapply(p3, function(r) r$lmin >= -1e-3 * r$lmax, TRUE)
cat("\nDecision rule (fixed before the run): S1 W-only stall, S2 kappa, S3 R_W on slow modes, S4 not boundary, S5 no saddle\n")
print(rbind(S1 = S1, S2 = S2, S3 = S3, S4 = S4, S5 = S5))
maj <- function(v) sum(v) > length(v) / 2
sup <- maj(S1) && maj(S2) && maj(S3) && maj(S4) && all(S5)
cat(sprintf("control kappa (max over converged p=150 starts) %.3e\nCONDITIONING HYPOTHESIS: %s\n", k150, if (sup) "SUPPORTED" else "NOT SUPPORTED"))
sink()
