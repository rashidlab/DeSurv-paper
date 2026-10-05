# S2: theta_lib F vs the independent evaluator ref_F_theta (survival::coxph at init = theta / s0, iter.max = 0)
# at 20 random points: random W in S_+^k, H > 0, theta (two points with a zero component), hp drawn from the
# spec section 3 box. p = 80, n = 50 with tied times. Pass: relative difference <= 1e-10. Also checks s0 agrees.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S2_ref_objective.txt")
cat("S2 reference objective |", format(Sys.time()), "\n")
s <- sim_data(12, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d\n", D$p, D$n, D$ne, s$n_tied))
TOL <- 1e-10
rF <- numeric(20); rs0 <- numeric(20)
for (i in 1:20) {
  set.seed(200 + i)
  hp <- rand_hp(); W <- rand_W(D$p); H <- matrix(runif(3 * D$n, 0, 50), 3, D$n); th <- rnorm(3, 0, 3)
  if (i %in% c(5, 15)) th[2] <- 0
  a <- palm_F(W, H, th, D, hp)
  b <- ref_F_theta(W, H, th, D$X, D$y, D$d, hp$alpha, hp$lambda, hp$nu, hp$rho, hp$lambdaH)
  rF[i] <- abs(a - b$F) / max(abs(a), abs(b$F)); rs0[i] <- abs(D$s0 - b$s0) / D$s0
  cat(sprintf("  point %2d: alpha %.3e lambda %.3e nu %.3f rho %.3e  F lib %.15e ref %.15e rel %.2e\n",
              i, hp$alpha, hp$lambda, hp$nu, hp$rho, a, b$F, rF[i]))
}
cat("\n")
tally("F: theta_lib vs ref_F_theta", rF, TOL)
tally("s0: theta_lib vs ref_F_theta", rs0, TOL)
total("S2")
sink()
