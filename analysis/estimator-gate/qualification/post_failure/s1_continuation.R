# Post-failure diagnostic (PI ruling 19). Cannot alter the Stage 1 FAIL or any criterion. Question: do S1 starts that
# end nonconverged at 5,000 sweeps converge slowly or stall? Inputs fixed before any run: S1 replicate 1, centre
# configuration at n = 175 (g = 5), the first three starts (by index) whose Stage 1 status is "nonconverged".
# The grid outputs keep no endpoint H, so each start is rerun deterministically from its frozen seed to 25,000
# sweeps with trace on. Gate: F at sweep 5,000 identical to the Stage 1 stored final F of that start.
# Usage (repo root, QUAL_LIB set): QUAL_MODE=qual NCORES=3 Rscript analysis/estimator-gate/qualification/post_failure/s1_continuation.R
QUAL <- normalizePath("analysis/estimator-gate/qualification"); source(file.path(QUAL, "qual_common.R"))
g1 <- readRDS(file.path(QUAL, "out", "grid", "S1_r01_g05.rds"))
ss <- head(g1$starts$start[g1$starts$status == "nonconverged"], 3L); stopifnot(length(ss) == 3L)
z <- qual_data("S1", 1L); D <- palm_data(z$X[, 1:175], z$y[1:175], z$d[1:175])
hp <- list(alpha = 0.05, lambda = 1, nu = 0.5, rho = 1, lambdaH = LAMBDA_H)
seed <- seed_grid(1L, 1L, 5L)
Ws <- lapply(seq_len(30L), function(s) proj_S(palm_draw(D, 3L, seed, s)))
one <- function(s) {
  W <- Ws[[s]]; H <- palm_init_H(W, D, hp)
  f <- palm_fit_one(D, hp, W, H, rep(0, 3), maxit = 25000L, tol = QUAL_TOL, trace = TRUE)
  tr <- f$trace; if (is.list(tr) && !is.matrix(tr)) tr <- do.call(rbind, tr)
  list(start = s, status = f$status, iter = f$iter, r_end = f$r, trace = tr)
}
res <- parallel::mclapply(ss, one, mc.cores = 3L, mc.preschedule = FALSE)
for (x in res) {
  tr <- x$trace; F5000 <- tr[tr[, "iter"] == 5000, "Ftheta"]
  Fst <- g1$starts$F[g1$starts$start == x$start]
  at <- function(k) if (any(tr[, "iter"] == k)) tr[tr[, "iter"] == k, "r"] else NA
  cat(sprintf("start %d: gate F(5000) == Stage 1 F: %s | status %s at %s sweeps | r at 1k %.2e, 5k %.2e, 10k %.2e, 15k %.2e, 20k %.2e, 25k %.2e\n",
              x$start, isTRUE(all.equal(F5000, Fst, tolerance = 0)), x$status, x$iter, at(1000), at(5000), at(10000),
              at(15000), at(20000), at(25000)))
}
saveRDS(list(starts = ss, res = res), file.path(QUAL, "post_failure", "s1_continuation.rds"))
