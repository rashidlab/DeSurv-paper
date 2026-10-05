# T9: constructed edge cases (simulated or hand-built data only).
# W/H/b boundary cases; Breslow tie patterns, row order, one event, large contrasts and offsets against
# coxph and an independent log-sum-exp reference; invalid input fails deterministically with a clear error;
# sample-permutation invariance, Cox offset invariance, seed replay; amendment-2 behaviour: sphere KKT
# residual, curvature-scaled H residual, W step cap, rejected b step, start statuses and the fit status rule.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressMessages(library(survival))
open_out("T9_edge_cases.txt")
cat("T9 edge cases |", format(Sys.time()), "\n")
nchk <- 0L; npass <- 0L
ok <- function(cond, label) {
  nchk <<- nchk + 1L; npass <<- npass + isTRUE(cond)
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "pass" else "FAIL", label))
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)
hp0 <- list(alpha = 0.3, lambda = 1, nu = 0.5, rho = 1, lambdaH = 1e-6)
s <- sim_data(19, 60, 40); D <- palm_data(s$X, s$y, s$d)
W0 <- proj_S(palm_draw(D, 3L, 1L, 1L)); hp <- hp0; hp$eps <- palm_eps(list(W0), D)
fd <- function(f, h = 1e-4) (f(h) - f(-h)) / (2 * h)

cat("\n1. W, H, b boundary cases\n")
Wu <- W0; Wu[, 2] <- 1 / sqrt(D$p)                                  # uniform column: X^T w constant (ranks)
sc <- palm_scores(Wu, D, hp$eps)
ok(abs(sc$sigma[2] - hp$eps) / hp$eps < 1e-3, sprintf("uniform column: sigma_2 = %.3e vs eps %.3e (sigma at the floor)", sc$sigma[2], hp$eps))
Hu <- palm_init_H(Wu, D, hp); bu <- c(0.2, 1e-9, -0.1)
ok(all(is.finite(c(palm_F(Wu, Hu, bu, D, hp), palm_grad_W(Wu, Hu, bu, D, hp), palm_resid(Wu, Hu, bu, D, hp)))),
   "uniform column: F, grad_W and residuals finite")
med_u <- median(palm_scores(matrix(1 / sqrt(D$p), D$p, 3), D, 0)$sigma)
ok(isTRUE(all.equal(palm_eps(list(matrix(1 / sqrt(D$p), D$p, 3)), D), 1e-6 * med_u)) && med_u < 1e-9,
   sprintf("all-uniform starts on rank data: pooled median sd is rounding-level (%.2e), not 0, so palm_eps returns 1e-6 x it", med_u))
Wn <- Wu; Wn[, 2] <- proj_S(matrix(1 + 1e-7 * runif(D$p)))[, 1]       # nearly uniform
fn <- palm_fit_one(D, hp, Wn, palm_init_H(Wn, D, hp), rep(0, 3), maxit = 30L)
ok(fn$status %in% c("nonconverged", "converged") && sum(fn$ninc) == 0, sprintf("nearly uniform column: 30 iterations, status %s, %d increases", fn$status, sum(fn$ninc)))
Wd <- W0; Wd[, 3] <- Wd[, 1]                                        # duplicate columns
fdup <- palm_fit_one(D, hp, Wd, palm_init_H(Wd, D, hp), rep(0, 3), maxit = 30L)
ok(fdup$status != "failure" && sum(fdup$ninc) == 0, sprintf("duplicate columns: 30 iterations, status %s, no increase", fdup$status))
H0 <- palm_init_H(W0, D, hp); G0 <- grad_W_terms(W0, H0, rep(0, 3), D, hp)
ok(max(abs(G0$cox)) == 0 && max(abs(G0$ridge)) == 0, "b = 0: Cox and rho-ridge W gradients exactly 0")
dv <- matrix(rnorm(D$p * 3), D$p, 3)
ok(rel(sum(palm_grad_W(W0, H0, rep(0, 3), D, hp) * dv), fd(function(e) palm_Psi(W0 + e * dv, H0, rep(0, 3), D, hp))) < 1e-6,
   "b = 0: directional derivative matches finite difference")
Ws <- W0; Ws[sample(D$p * 3, round(0.8 * D$p * 3))] <- 0; Ws <- proj_S(Ws)   # sparse boundary W
Hs <- palm_init_H(Ws, D, hp); Hs[2, ] <- 0; Hs[1, 1:10] <- 0                  # zero row and zeros in H
r_s <- palm_resid(Ws, Hs, c(0.1, 0, -0.2), D, hp)
ok(all(is.finite(r_s)), sprintf("sparse W (%.0f%% zeros) and H with a zero row: residuals finite (r = %.3e)", 100 * mean(Ws == 0), r_s[["r"]]))
fsp <- palm_fit_one(D, hp, Ws, Hs, c(0.1, 0, -0.2), maxit = 30L)
ok(fsp$status != "failure" && sum(fsp$ninc) == 0, sprintf("sparse start: 30 iterations, status %s, no increase", fsp$status))
GH0 <- palm_grad_H(W0, matrix(0, 3, D$n), D, hp)
ok(is.finite(PALM_DEFS$resid_H(matrix(0, 3, D$n), GH0, W0, D, hp)), "r_H defined at H = 0 (denominator max(1, ||H||))")

cat("\n2. Breslow: tie patterns, row order, one event, contrasts, offsets\n")
cox_ll <- function(eta, y, d) coxph(Surv(y, d) ~ offset(eta), ties = "breslow")$loglik[1]
lse_ref <- function(eta, y, d) {                                      # independent reference, log-sum-exp per risk set
  ll <- 0; g <- d
  for (i in which(d == 1)) { rs <- y >= y[i]; m <- max(eta[rs]); ls <- m + log(sum(exp(eta[rs] - m))); ll <- ll + eta[i] - ls
    g[rs] <- g[rs] - exp(eta[rs] - ls) }
  list(ll = ll, g = g)
}
mk <- function(y, d) list(y = y, d = d, D = palm_data(matrix(runif(2 * length(y)) + 1, 2), y, d))
set.seed(91)
tie_sets <- list(
  event_event = mk(c(1, 2, 2, 2, 3, 4, 5, 5), c(1, 1, 1, 1, 0, 1, 1, 1)),
  event_censor = mk(c(1, 2, 2, 3, 3, 4, 5, 6), c(1, 1, 0, 0, 1, 1, 0, 1)),
  censor_only_group = mk(c(1, 2, 3, 3, 3, 4, 5, 6), c(1, 0, 0, 0, 0, 1, 0, 1)),
  unsorted_rows = mk(c(5, 1, 3, 3, 2, 6, 2, 4), c(1, 0, 1, 1, 1, 0, 1, 1)),
  one_event = mk(c(3, 1, 2, 5, 4, 6), c(0, 0, 1, 0, 0, 0)),
  all_tied = mk(rep(2, 6), c(1, 0, 1, 1, 0, 1)))
for (nm in names(tie_sets)) {
  z <- tie_sets[[nm]]; eta <- rnorm(length(z$y))
  bz <- breslow(eta, z$D); rf <- lse_ref(eta, z$y, z$d)
  ok(abs(bz$ll - cox_ll(eta, z$y, z$d)) < 1e-12 * max(1, abs(bz$ll)) && rel(bz$g, rf$g) < 1e-12,
     sprintf("%-18s ll vs coxph %.1e, gradient vs reference %.1e", nm, abs(bz$ll - cox_ll(eta, z$y, z$d)), rel(bz$g, rf$g)))
}
z <- tie_sets$event_censor
for (sc_ in c(50, 500, 2000)) {
  eta <- c(sc_, -sc_, 0.5 * sc_, -0.3 * sc_, sc_ / 7, -sc_ / 3, 0, sc_ / 2)
  bz <- breslow(eta, z$D); rf <- lse_ref(eta, z$y, z$d)
  ok(is.finite(bz$ll) && all(is.finite(bz$g)) && abs(bz$ll - rf$ll) < 1e-10 * max(1, abs(rf$ll)) && max(abs(bz$g - rf$g)) < 1e-12,
     sprintf("contrast range %5g: ll %.6g vs reference %.6g, max grad diff %.1e", 2 * sc_, bz$ll, rf$ll, max(abs(bz$g - rf$g))))
}
eta <- rnorm(8)
for (off in c(1e3, -1e3, 1e5)) {
  b1 <- breslow(eta, z$D); b2 <- breslow(eta + off, z$D)
  ok(abs(b1$ll - b2$ll) < 1e-9 * max(1, abs(b1$ll)) && max(abs(b1$g - b2$g)) < 1e-9,
     sprintf("common offset %+g: ll change %.1e, gradient change %.1e", off, abs(b1$ll - b2$ll), max(abs(b1$g - b2$g))))
}
s1 <- sim_data(23, 40, 12); s1$d[] <- 0; s1$d[3] <- 1
D1 <- palm_data(s1$X, s1$y, s1$d); f1 <- palm_fit(D1, hp0, nstarts = 2L, maxit = 20L)
ok(D1$ne == 1 && f1$status %in% c("converged", "nonconverged"), sprintf("one-event dataset: 2 starts x 20 iterations, fit status %s", f1$status))

cat("\n3. invalid input fails deterministically with a clear error\n")
err <- function(expr) tryCatch({ force(expr); "no error" }, error = function(e) conditionMessage(e))
Xb <- s$X
cases <- list(
  n_lt_2 = list(quote(palm_data(Xb[, 1, drop = FALSE], s$y[1], s$d[1])), "n >= 2"),
  no_events = list(quote(palm_data(Xb, s$y, 0 * s$d)), "no events"),
  X_NA = list(quote({ X2 <- Xb; X2[3, 4] <- NA; palm_data(X2, s$y, s$d) }), "non-finite"),
  X_NaN = list(quote({ X2 <- Xb; X2[3, 4] <- NaN; palm_data(X2, s$y, s$d) }), "non-finite"),
  X_Inf = list(quote({ X2 <- Xb; X2[3, 4] <- Inf; palm_data(X2, s$y, s$d) }), "non-finite"),
  X_zero = list(quote(palm_data(0 * Xb, s$y, s$d)), "\\|\\|X\\|\\| = 0"),
  y_Inf = list(quote({ y2 <- s$y; y2[2] <- Inf; palm_data(Xb, y2, s$d) }), "y has non-finite"),
  y_NA = list(quote({ y2 <- s$y; y2[2] <- NA; palm_data(Xb, y2, s$d) }), "y has non-finite"),
  d_NA = list(quote({ d2 <- s$d; d2[2] <- NA; palm_data(Xb, s$y, d2) }), "d has non-finite"),
  d_NaN = list(quote({ d2 <- as.numeric(s$d); d2[2] <- NaN; palm_data(Xb, s$y, d2) }), "d has non-finite"),
  sigma_zero = list(quote(palm_fit(palm_data(matrix(rep(1:5, 6), 5, 6), 1:6, c(1, 0, 1, 0, 1, 1)), hp0, nstarts = 2L, maxit = 2L)), "pooled median"))
for (nm in names(cases)) {
  m1 <- err(eval(cases[[nm]][[1]])); m2 <- err(eval(cases[[nm]][[1]]))
  ok(grepl(cases[[nm]][[2]], m1) && identical(m1, m2), sprintf("%-11s -> \"%s\"", nm, m1))
}

cat("\n4. invariances and replay\n")
pm <- sample(D$n)
Dp <- palm_data(s$X[, pm], s$y[pm], s$d[pm])
fa <- palm_fit_one(D, hp, W0, palm_init_H(W0, D, hp), rep(0, 3), maxit = 40L)
fb <- palm_fit_one(Dp, hp, W0, palm_init_H(W0, Dp, hp), rep(0, 3), maxit = 40L)
ok(rel(fb$W, fa$W) < 1e-8 && rel(fb$b, fa$b) < 1e-8 && abs(fb$F - fa$F) < 1e-10 * abs(fa$F) && rel(fb$H, fa$H[, pm]) < 1e-8,
   sprintf("sample permutation, 40 iterations: W rel diff %.1e, b %.1e, F %.1e, H (permuted) %.1e",
           rel(fb$W, fa$W), rel(fb$b, fa$b), abs(fb$F - fa$F) / abs(fa$F), rel(fb$H, fa$H[, pm])))
Xo <- s$X + 0                                                       # Cox offset via F: shift all scores by a constant
Wt <- W0; bt <- c(0.3, -0.2, 0.1)
F1 <- palm_parts(Wt, H0, bt, D, hp)$loglik
eta <- as.vector(palm_scores(Wt, D, hp$eps)$Zc %*% (bt / palm_scores(Wt, D, hp$eps)$sigma))
ok(abs(breslow(eta + 37, D)$ll - F1) < 1e-10 * abs(F1), "loglik invariant to a common offset of the centered scores")
r1 <- palm_fit(D, hp0, seed = 77L, nstarts = 3L, maxit = 25L); r2 <- palm_fit(D, hp0, seed = 77L, nstarts = 3L, maxit = 25L)
r3 <- palm_fit(D, hp0, seed = 77L, nstarts = 3L, maxit = 25L, ncores = 3L)
strip <- function(f) { f$best$secs <- NULL; f$starts$secs <- NULL; f }
ok(identical(strip(r1)$best[c("W", "H", "b", "F")], strip(r2)$best[c("W", "H", "b", "F")]) &&
   identical(strip(r1)$starts, strip(r2)$starts), "seed replay: two serial runs identical")
ok(identical(strip(r1)$best[c("W", "H", "b", "F")], strip(r3)$best[c("W", "H", "b", "F")]) &&
   identical(strip(r1)$starts, strip(r3)$starts), "seed replay: serial vs 3-core run identical")

cat("\n5. amendment 2: residuals, step cap, b acceptance, statuses\n")
p <- 6; We <- matrix(0, p, 2); We[1, 1] <- 1; We[3, 2] <- 1
Grad <- matrix(0, p, 2); Grad[1, 1] <- 5; Grad[3, 2] <- 3; Grad[5, 2] <- 0.7     # radial (outward) + positive on a zero
rWe <- PALM_DEFS$resid_W(We, Grad, NULL, NULL); rOe <- sqrt(sum((We - proj_S(We - Grad))^2))
ok(rWe == 0 && rOe > 1, sprintf("stationary point w = e_i with radial gradient (|G| > 1): r_W = %g, original W residual = %.3f", rWe, rOe))
Grad2 <- Grad; Grad2[2, 1] <- -0.4
ok(abs(PALM_DEFS$resid_W(We, Grad2, NULL, NULL) - 0.4) < 1e-15, "negative gradient on a zero coordinate: r_W = |min(G, 0)| = 0.4")
Grad3 <- Grad; Grad3[1, 1] <- 5; We3 <- We; We3[1:2, 1] <- c(0.6, 0.8); Grad3[1:2, 1] <- c(1, 0)
mu <- sum(We3[, 1] * Grad3[, 1]); Rref <- Grad3[1:2, 1] - mu * We3[1:2, 1]
ok(abs(PALM_DEFS$resid_W(We3, Grad3, NULL, NULL) - sqrt(sum(Rref^2))) < 1e-14,
   "tangential residual on the support matches G - (w^T G) w (positive G on a zero coordinate adds 0)")
Hh <- matrix(c(2, 0, 0.5, 3), 2); GHh <- matrix(c(-1, 0.3, 4, 0), 2)
aH <- (1 - hp$alpha) * (1 / D$X2 + hp$lambdaH / (D$n * 2))
ok(abs(PALM_DEFS$resid_H(Hh, GHh, matrix(0, 3, 2), D, hp) - sqrt(sum(pmin(Hh, GHh / aH)^2)) / max(1, sqrt(sum(Hh^2)))) < 1e-15,
   "r_H = ||min(H, G_H / a_H)|| / max(1, ||H||)")
ok(PALM_DEFS$step_init(1) == 2 && PALM_DEFS$step_init(1.7) == 2 && PALM_DEFS$step_init(0.25) == 0.5, "step trial min(2 t_prev, 2)")
ft <- palm_fit_one(D, hp, W0, H0, rep(0, 3), maxit = 60L, trace = TRUE)
ok(max(ft$trace[, "t"]) <= 2, sprintf("accepted W steps never exceed 2 (max %.3g over 60 iterations)", max(ft$trace[, "t"])))
solve_b_good <- solve_b
solve_b <- function(b0, W, D, hp, L0 = 1, tol = 1e-10, maxit = 1e5) { z <- solve_b_good(b0, W, D, hp, L0, tol, maxit); z$b <- z$b + 1; z }
fr <- palm_fit_one(D, hp, W0, H0, rep(0, 3), maxit = 5L)
ok(fr$b_rejected == 5 && identical(fr$b, rep(0, 3)) && fr$status == "nonconverged" && sum(fr$ninc) == 0,
   sprintf("b+ that raises F is rejected: %d of 5 rejected, b kept, status %s, increases %d", fr$b_rejected, fr$status, sum(fr$ninc)))
solve_b <- function(b0, W, D, hp, L0 = 1, tol = 1e-10, maxit = 1e5) solve_b_good(b0, W, D, hp, L0, tol, maxit = 2)
fc <- palm_fit_one(D, hp, W0, H0, rep(0, 3), maxit = 5L)
ok(fc$status == "failure" && grepl("FISTA", fc$reason) && fc$iter == 1, sprintf("FISTA cap -> numerical failure at iteration %d (%s)", fc$iter, fc$reason))
solve_b <- solve_b_good
# The floor needs a step that never becomes exactly null: as t shrinks, W - t G rounds to W on positive
# coordinates and the null step satisfies the inequality with equality; only zero coordinates with negative
# gradient keep W+ != W. Psi is mutated to rise by 1 for any W other than the start, on a sparse start.
Gs0 <- palm_grad_W(Ws, Hs, rep(0, 3), D, hp)
psi_good <- palm_Psi
palm_Psi <- function(W, H, b, D, hp) psi_good(W, H, b, D, hp) + if (identical(W, Ws)) 0 else 1
fs <- palm_fit_one(D, hp, Ws, Hs, rep(0, 3), maxit = 5L)
palm_Psi <- psi_good
ok(min(Gs0[Ws == 0]) < 0 && fs$status == "failure" && grepl("floor", fs$reason) && identical(fs$W, Ws) && fs$iter == 1,
   sprintf("W step floor -> numerical failure at iteration %d, W unchanged (%s)", fs$iter, fs$reason))
grad_good <- palm_grad_W                                         # observation (not counted): ascent direction
palm_grad_W <- function(W, H, b, D, hp) -grad_good(W, H, b, D, hp)
fasc <- palm_fit_one(D, hp, W0, H0, rep(0, 3), maxit = 5L, trace = TRUE)
palm_grad_W <- grad_good
cat(sprintf("  [obs ] ascent-direction W gradient on a dense start: status %s after %d iterations (%s); accepted t %s; W moved by %.1e; F change over the W blocks %.1e\n",
            fasc$status, fasc$iter, if (nzchar(fasc$reason)) fasc$reason else "no failure", paste(signif(fasc$trace[, "t"], 2), collapse = ","),
            max(abs(fasc$W - W0)), max(fasc$trace[, "FW"] - fasc$trace[, "FH"], na.rm = TRUE)))
hals_good <- hals_sweep
hals_sweep <- function(W, H, D, hp) hals_good(W, H, D, hp) * 1.5           # an H block that raises F
fi <- palm_fit_one(D, hp, W0, H0, rep(0, 3), maxit = 5L)
fF <- palm_fit(D, hp0, nstarts = 2L, maxit = 5L)
hals_sweep <- hals_good
ok(fi$status == "failure" && grepl("increase", fi$reason) && fi$iter == 1, sprintf("F increase -> numerical failure, start ends at iteration %d (%s)", fi$iter, fi$reason))
fn5 <- palm_fit_one(D, hp, W0, H0, rep(0, 3), maxit = 3L)
ok(fn5$status == "nonconverged" && fn5$iter == 3, "maxit reached with finite values -> nonconverged")
fcv <- palm_fit_one(D, hp, W0, H0, rep(0, 3), maxit = 5000L, tol = 1e3)
ok(fcv$status == "converged" && fcv$iter == 1, "r < tol -> converged, start ends at once")
ok(identical(pick_start(c(1, 0.5, 2), c("converged", "nonconverged", "converged")), list(i = 1L, status = "converged")),
   "fit rule: lowest-F converged start returned although a nonconverged start has lower F")
ok(identical(pick_start(c(3, 0.5, 2), c("nonconverged", "failure", "nonconverged")), list(i = 3L, status = "nonconverged")),
   "fit rule: no converged start -> lowest-F nonconverged, failure never returned")
ok(identical(pick_start(c(NA, 0.5), c("failure", "failure")), list(i = NA_integer_, status = "failed")), "fit rule: all failures -> fit failed")
ok(fF$status == "failed" && is.null(fF$best) && tryCatch({ palm_model(fF, D, 10); FALSE }, error = function(e) TRUE),
   "fit with every start failed: status failed, no model can be built (never scored)")

cat(sprintf("\nT9 TOTAL: %d checks, %d passed\n", nchk, npass))
verdict(nchk == npass)
sink()
