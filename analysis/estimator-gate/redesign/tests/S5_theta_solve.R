# S5: theta block. (a) solve_theta from theta = 0 at 10 random (W, hp): converged flag and independent-oracle
# KKT <= 1e-10. (b) theta_lib theta-subproblem gradient vs the oracle's own gradient at the returned theta and at
# 3 random theta per configuration: max|g_lib - g_oracle| / max(1, max|g_oracle|) <= 1e-10.
# p = 80, n = 50 with tied times; hp drawn from the spec section 3 box.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S5_theta_solve.txt")
cat("S5 theta solve |", format(Sys.time()), "\n")
s <- sim_data(15, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d\n", D$p, D$n, D$ne, s$n_tied))
kk <- numeric(10); conv <- logical(10); gd <- numeric(0)
for (i in 1:10) {
  set.seed(500 + i); hp <- rand_hp(); W <- rand_W(D$p)
  st <- solve_theta(rep(0, 3), W, D, hp)
  kk[i] <- kkt_oracle_theta(st$theta, W, D$oracle, hp$lambda, hp$nu, hp$rho); conv[i] <- st$conv
  for (th in c(list(st$theta), replicate(3, rnorm(3, 0, 3), simplify = FALSE))) {
    gl <- tsub_grad(W, th, D, hp); go <- kkt_oracle_theta_grad(th, W, D$oracle, hp$lambda, hp$nu, hp$rho)
    gd <- c(gd, max(abs(gl - go)) / max(1, max(abs(go))))
  }
  cat(sprintf("  cfg %2d: alpha %.2e lambda %.2e nu %.3f rho %.2e  conv %s it %d  theta = (%s)  oracle KKT %.2e\n",
              i, hp$alpha, hp$lambda, hp$nu, hp$rho, st$conv, st$it, paste(sprintf("%.4g", st$theta), collapse = ", "), kk[i]))
}
cat("\n")
tally("solve_theta converged flag (1 = not converged)", as.numeric(!conv), 0)
tally("oracle KKT after solve_theta", kk, 1e-10)
tally("theta_lib vs oracle gradient (10 x 4 points)", gd, 1e-10)
total("S5")
sink()
