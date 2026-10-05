# S4: HALS H block (ruling 12 denominator 1 + lambdaH). (a) After one sweep from a random H, each row h_j
# satisfies its row KKT at the state right after its update (rows 1..j from the new H, rows j+1..k from the old H,
# which is the H at which row j was minimized): |min(h_j, G_hj / a_H)|_inf / max(1, |h_j|_inf) <= 1e-12, with
# G_hj = row j of palm_grad_H and a_H = (1 - alpha)(1 + lambdaH)/X2 (the row-Hessian, so both arguments are on the
# H scale). (b) F non-increasing over 50 sweeps at fixed (W, theta): F_new <= F_old + 1e-13 max(1, |F_old|).
# 9 configurations: alpha {1e-4, 0.3, 0.9} x 3 (W, H, theta): reps 1-2 random; rep 3 has a very large row h_3 in H0
# (1e4 + H0), so the row-1 and row-2 updates are clipped at zero (exercises the h_j = 0 branch of the KKT).
# p = 80, n = 50 with tied times.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S4_hals.txt")
cat("S4 HALS |", format(Sys.time()), "\n")
s <- sim_data(14, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d\n", D$p, D$n, D$ne, s$n_tied))
kk <- numeric(0); mono <- numeric(0); rawinc <- 0L
for (alpha in c(1e-4, 0.3, 0.9)) for (rep in 1:3) {
  set.seed(400 + round(10 * alpha) + rep)
  hp <- list(alpha = alpha, lambda = 1, nu = 0.5, rho = 1, lambdaH = 1e-6)
  W <- rand_W(D$p); H0 <- matrix(runif(3 * D$n, 0, 50), 3, D$n); th <- rnorm(3)
  if (rep == 3) H0[3, ] <- H0[3, ] + 1e4
  H1 <- hals_sweep(W, H0, D, hp); aH <- (1 - alpha) * (1 + hp$lambdaH) / D$X2
  for (j in 1:3) {
    Hs <- H0; Hs[1:j, ] <- H1[1:j, ]
    g <- palm_grad_H(W, Hs, D, hp)[j, ]
    kk <- c(kk, max(abs(pmin(Hs[j, ], g / aH))) / max(1, max(abs(Hs[j, ]))))
  }
  H <- H0; f <- palm_F(W, H, th, D, hp)
  for (i in 1:50) {
    H <- hals_sweep(W, H, D, hp); fn <- palm_F(W, H, th, D, hp)
    mono <- c(mono, max(0, fn - f) / max(1, abs(f))); if (fn > f) rawinc <- rawinc + 1L; f <- fn
  }
  cat(sprintf("  alpha %.0e rep %d: zero entries after sweep 1 = %d of %d\n", alpha, rep, sum(H1 == 0), length(H1)))
}
cat("\n")
tally("row KKT after one sweep (9 configs x 3 rows)", kk, 1e-12)
tally("F increase over 50 sweeps (relative, 9 x 50)", mono, 1e-13)
cat(sprintf("  strict raw increases (fn > f, any size): %d\n", rawinc))
total("S4")
sink()
