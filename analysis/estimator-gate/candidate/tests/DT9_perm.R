# D-T9a / D-T9b (post-failure diagnostics, PI approval 2026-10-05). Descriptive only; changes nothing about
# the recorded amendment-3 FAIL (0c44cca) or the T9 criterion. palm_lib.R is read unmodified; no existing
# test is changed. This script reproduces the T9 permutation pair (tests/T9_edge_cases.R section 4) by
# running the identical code T9 runs before `pm <- sample(D$n)` (sections 1-3, verbatim, output suppressed),
# so the RNG stream reaching `pm <- sample(D$n)` is the same path, then compares the pair under the frozen
# D-T9a/b protocol.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressMessages(library(survival))

# ---- verbatim replica of tests/T9_edge_cases.R, sections 1-3 (output suppressed; this is what builds the
# RNG state T9 carries into `pm <- sample(D$n)`; nchk/npass/ok/rel/hp0/s/D/W0/hp/fd/H0/Ws/Hs are produced here) ----
invisible(capture.output({
  nchk <- 0L; npass <- 0L
  ok <- function(cond, label) { nchk <<- nchk + 1L; npass <<- npass + isTRUE(cond)
    cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "pass" else "FAIL", label)) }
  rel <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)
  hp0 <- list(alpha = 0.3, lambda = 1, nu = 0.5, rho = 1, lambdaH = 1e-6)
  s <- sim_data(19, 60, 40); D <- palm_data(s$X, s$y, s$d)
  W0 <- proj_S(palm_draw(D, 3L, 1L, 1L)); hp <- hp0; hp$eps <- palm_eps(list(W0), D)
  fd <- function(f, h = 1e-4) (f(h) - f(-h)) / (2 * h)
  Wu <- W0; Wu[, 2] <- 1 / sqrt(D$p)
  sc <- palm_scores(Wu, D, hp$eps)
  ok(abs(sc$sigma[2] - hp$eps) / hp$eps < 1e-3, "uniform column sigma")
  Hu <- palm_init_H(Wu, D, hp); bu <- c(0.2, 1e-9, -0.1)
  ok(all(is.finite(c(palm_F(Wu, Hu, bu, D, hp), palm_grad_W(Wu, Hu, bu, D, hp), palm_resid(Wu, Hu, bu, D, hp)))), "uniform column finite")
  med_u <- median(palm_scores(matrix(1 / sqrt(D$p), D$p, 3), D, 0)$sigma)
  ok(isTRUE(all.equal(palm_eps(list(matrix(1 / sqrt(D$p), D$p, 3)), D), 1e-6 * med_u)) && med_u < 1e-9, "all-uniform palm_eps")
  Wn <- Wu; Wn[, 2] <- proj_S(matrix(1 + 1e-7 * runif(D$p)))[, 1]
  fn <- palm_fit_one(D, hp, Wn, palm_init_H(Wn, D, hp), rep(0, 3), maxit = 30L)
  ok(fn$status %in% c("nonconverged", "converged") && sum(fn$ninc) == 0, "nearly uniform column")
  Wd <- W0; Wd[, 3] <- Wd[, 1]
  fdup <- palm_fit_one(D, hp, Wd, palm_init_H(Wd, D, hp), rep(0, 3), maxit = 30L)
  ok(fdup$status != "failure" && sum(fdup$ninc) == 0, "duplicate columns")
  H0 <- palm_init_H(W0, D, hp); G0 <- grad_W_terms(W0, H0, rep(0, 3), D, hp)
  ok(max(abs(G0$cox)) == 0 && max(abs(G0$ridge)) == 0, "b = 0 zero gradients")
  dv <- matrix(rnorm(D$p * 3), D$p, 3)
  ok(rel(sum(palm_grad_W(W0, H0, rep(0, 3), D, hp) * dv), fd(function(e) palm_Psi(W0 + e * dv, H0, rep(0, 3), D, hp))) < 1e-6, "b = 0 directional derivative")
  Ws <- W0; Ws[sample(D$p * 3, round(0.8 * D$p * 3))] <- 0; Ws <- proj_S(Ws)
  Hs <- palm_init_H(Ws, D, hp); Hs[2, ] <- 0; Hs[1, 1:10] <- 0
  r_s <- palm_resid(Ws, Hs, c(0.1, 0, -0.2), D, hp)
  ok(all(is.finite(r_s)), "sparse W residuals finite")
  fsp <- palm_fit_one(D, hp, Ws, Hs, c(0.1, 0, -0.2), maxit = 30L)
  ok(fsp$status != "failure" && sum(fsp$ninc) == 0, "sparse start")
  GH0 <- palm_grad_H(W0, matrix(0, 3, D$n), D, hp)
  ok(is.finite(PALM_DEFS$resid_H(matrix(0, 3, D$n), GH0, W0, D, hp)), "r_H at H = 0")

  cox_ll <- function(eta, y, d) coxph(Surv(y, d) ~ offset(eta), ties = "breslow")$loglik[1]
  lse_ref <- function(eta, y, d) { ll <- 0; g <- d
    for (i in which(d == 1)) { rs <- y >= y[i]; m <- max(eta[rs]); ls <- m + log(sum(exp(eta[rs] - m))); ll <- ll + eta[i] - ls
      g[rs] <- g[rs] - exp(eta[rs] - ls) }
    list(ll = ll, g = g) }
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
    ok(abs(bz$ll - cox_ll(eta, z$y, z$d)) < 1e-12 * max(1, abs(bz$ll)) && rel(bz$g, rf$g) < 1e-12, nm)
  }
  z <- tie_sets$event_censor
  for (sc_ in c(50, 500, 2000)) {
    eta <- c(sc_, -sc_, 0.5 * sc_, -0.3 * sc_, sc_ / 7, -sc_ / 3, 0, sc_ / 2)
    bz <- breslow(eta, z$D); rf <- lse_ref(eta, z$y, z$d)
    ok(is.finite(bz$ll) && all(is.finite(bz$g)) && abs(bz$ll - rf$ll) < 1e-10 * max(1, abs(rf$ll)) && max(abs(bz$g - rf$g)) < 1e-12, "contrast")
  }
  eta <- rnorm(8)
  for (off in c(1e3, -1e3, 1e5)) {
    b1 <- breslow(eta, z$D); b2 <- breslow(eta + off, z$D)
    ok(abs(b1$ll - b2$ll) < 1e-9 * max(1, abs(b1$ll)) && max(abs(b1$g - b2$g)) < 1e-9, "offset")
  }
  s1 <- sim_data(23, 40, 12); s1$d[] <- 0; s1$d[3] <- 1
  D1 <- palm_data(s1$X, s1$y, s1$d); f1 <- palm_fit(D1, hp0, nstarts = 2L, maxit = 20L)
  ok(D1$ne == 1 && f1$status %in% c("converged", "nonconverged"), "one event")

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
    ok(grepl(cases[[nm]][[2]], m1) && identical(m1, m2), nm)
  }
}))

# ---- the T9 permutation pair, exactly as tests/T9_edge_cases.R section 4 builds it ----
pm <- sample(D$n)
Dp <- palm_data(s$X[, pm], s$y[pm], s$d[pm])
fa <- palm_fit_one(D, hp, W0, palm_init_H(W0, D, hp), rep(0, 3), maxit = 40L)
fb <- palm_fit_one(Dp, hp, W0, palm_init_H(W0, Dp, hp), rep(0, 3), maxit = 40L)

rel_ <- function(a, b) max(abs(a - b)) / max(abs(b), 1e-300)
Wrel <- rel_(fb$W, fa$W); brel <- rel_(fb$b, fa$b); Frel <- abs(fb$F - fa$F) / abs(fa$F); Hrel <- rel_(fb$H, fa$H[, pm])

# ---- comparison helpers ----
eta_of <- function(f, Dx) { sc <- palm_scores(f$W, Dx, hp$eps); as.vector(sc$Zc %*% (f$b / sc$sigma)) }
match_perm3 <- function(C) { perms <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  perms[[which.max(sapply(perms, function(q) sum(C[cbind(1:3, q)])))]] }
colcos <- function(Wref, Wother) { C <- crossprod(Wref, Wother); q <- match_perm3(C); list(perm = q, cos = C[cbind(1:3, q)]) }
match_cos_truth <- function(Wt, W) { C <- crossprod(Wt, W); q <- match_perm3(C); C[cbind(1:3, q)] }
TOPN <- 5L   # = markers_per_factor used by sim_data for p = 60 (max(5, round(150*60/3000)) = 5); chosen before looking at results
overlap_topn <- function(Wref, Wother, perm, n = TOPN) {
  sapply(1:3, function(j) { ta <- order(Wref[, j], decreasing = TRUE)[1:n]; tb <- order(Wother[, perm[j]], decreasing = TRUE)[1:n]
    length(intersect(ta, tb)) / n })
}
compare_pair <- function(fa, fb, Dp, pm, label) {
  ra <- palm_resid(fa$W, fa$H, fa$b, D, hp); rb <- palm_resid(fb$W, fb$H, fb$b, Dp, hp)
  cm <- colcos(fa$W, fb$W)
  eta_a <- eta_of(fa, D); eta_b <- eta_of(fb, Dp)
  pear <- cor(eta_a[pm], eta_b); sdr <- sd(eta_b) / sd(eta_a[pm])
  ov <- overlap_topn(fa$W, fb$W, cm$perm)
  rec_a <- match_cos_truth(s$W_true, fa$W); rec_b <- match_cos_truth(s$W_true, fb$W)
  cat(sprintf("\n[%s]\n", label))
  cat(sprintf("  F: a %.10g, b %.10g, rel diff %.3e\n", fa$F, fb$F, abs(fa$F - fb$F) / abs(fa$F)))
  cat(sprintf("  KKT residuals a: r %.3e rW %.3e rH %.3e rb %.3e\n", ra[["r"]], ra[["rW"]], ra[["rH"]], ra[["rb"]]))
  cat(sprintf("  KKT residuals b: r %.3e rW %.3e rH %.3e rb %.3e\n", rb[["r"]], rb[["rW"]], rb[["rH"]], rb[["rb"]]))
  cat(sprintf("  Hungarian-matched (k=3, exhaustive over 3! perms) W-column cosines: %s (perm %s)\n",
              paste(sprintf("%.4f", cm$cos), collapse = ", "), paste(cm$perm, collapse = ",")))
  cat(sprintf("  risk-score (eta) Pearson correlation (b's order, a matched via pm): %.6f; SD ratio (b/a): %.4f\n", pear, sdr))
  cat(sprintf("  top-%d-loading overlap per matched column (markers_per_factor for p=%d): %s\n", TOPN, D$p, paste(sprintf("%.2f", ov), collapse = ", ")))
  cat(sprintf("  recovery of simulated truth (T8_recovery.R match_cos convention) a: %s\n", paste(sprintf("%.4f", rec_a), collapse = ", ")))
  cat(sprintf("  recovery of simulated truth (T8_recovery.R match_cos convention) b: %s\n", paste(sprintf("%.4f", rec_b), collapse = ", ")))
  invisible(list(ra = ra, rb = rb, cm = cm, pear = pear, sdr = sdr, ov = ov, rec_a = rec_a, rec_b = rec_b))
}

open_out("DT9_perm.txt")
options(width = 220)
cat("D-T9a/b permutation-pair diagnostic |", format(Sys.time()), "\n")

cat("\n=== D-T9a: seeded replay of the T9 permutation pair (40 iterations, as T9 runs it) ===\n")
pr <- sprintf("%.1e", c(Wrel, brel, Frel, Hrel)); committed <- c("7.3e-03", "5.4e-05", "5.3e-07", "5.2e-03")
cat(sprintf("GATE (vs committed tests/out/T9_edge_cases.txt line 45): replay W %s, b %s, F %s, H %s | committed W %s, b %s, F %s, H %s\n",
            pr[1], pr[2], pr[3], pr[4], committed[1], committed[2], committed[3], committed[4]))
gate_ok <- identical(pr, committed)
cat(sprintf("GATE: %s\n", if (gate_ok) "PASS (matches committed T9 line 45 to printed precision)" else "FAIL"))

out_obj <- list(gate_ok = gate_ok, Wrel = Wrel, brel = brel, Frel = Frel, Hrel = Hrel, fa = fa, fb = fb, pm = pm)
if (!gate_ok) {
  cat("\nD-T9a/b STOPPED: the replay did not reproduce the committed T9 permutation-pair differences. No further comparison run.\n")
  saveRDS(out_obj, file.path(OUT, "DT9_perm.rds"))
  verdict(FALSE)
} else {
  r40 <- compare_pair(fa, fb, Dp, pm, "D-T9a: 40-iteration (truncated) pair")

  cat("\n\n=== D-T9b: same pair, run to convergence (SEPARATE, NONBINDING); settings frozen in the brief: ===\n")
  cat("w_solver = ncg, tol = 1e-5, per-start budget 1,000 sweeps / 25,000 HVPs (palm_fit_one defaults)\n")
  fca <- palm_fit_one(D, hp, W0, palm_init_H(W0, D, hp), rep(0, 3), w_solver = "ncg", tol = 1e-5)
  fcb <- palm_fit_one(Dp, hp, W0, palm_init_H(W0, Dp, hp), rep(0, 3), w_solver = "ncg", tol = 1e-5)
  cat(sprintf("  member a: status %s (%s), sweeps %d, HVPs %d, final F %.10g, r %.3e\n", fca$status, fca$reason, fca$iter, fca$hvp_total, fca$F, fca$r[["r"]]))
  cat(sprintf("  member b: status %s (%s), sweeps %d, HVPs %d, final F %.10g, r %.3e\n", fcb$status, fcb$reason, fcb$iter, fcb$hvp_total, fcb$F, fcb$r[["r"]]))
  rconv <- compare_pair(fca, fcb, Dp, pm, "D-T9b: converged-path pair")

  cat("\n--- interpretation (fixed in advance; no threshold change) ---\n")
  both_conv <- fca$status == "converged" && fcb$status == "converged"
  Wrel40 <- Wrel; Wrelc <- rel_(fcb$W, fca$W)
  cat(sprintf("both converged: %s. truncated W rel diff %.3e; converged-path W rel diff %.3e\n", both_conv, Wrel40, Wrelc))
  if (!both_conv) {
    cat("At least one member did not converge within the frozen budget. That is the outcome; no rerun with other settings.\n")
  } else {
    i_fit <- Wrelc < 1e-4 && rconv$pear > 0.999
    ii_fit <- Wrelc > 1e-2
    iii_fit <- Wrelc > 1e-3 && rconv$pear > 0.999 && abs(fca$F - fcb$F) / abs(fca$F) < 1e-6
    cat(sprintf("(i) truncated differ, converged agree: %s\n", i_fit))
    cat(sprintf("(ii) converged differ materially: %s\n", ii_fit))
    cat(sprintf("(iii) parameters differ but objective/predictions/recovery agree (flat/nonidentifiable): %s\n", iii_fit))
  }
  out_obj$r40 <- r40; out_obj$rconv <- rconv; out_obj$fca <- fca; out_obj$fcb <- fcb
  saveRDS(out_obj, file.path(OUT, "DT9_perm.rds"))
  verdict(gate_ok)
}
sink()
