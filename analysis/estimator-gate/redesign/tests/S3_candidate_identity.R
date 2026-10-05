# S3: identity with the frozen candidate. With eps = 0 and theta = s0 b / sigma(W), the survival term and the rho
# term of theta_lib equal those of candidate/ref_objective.R ref_F (eps = 0). ref_F returns loglik and F only:
# survival term = -(2/ne) loglik; rho term = ref_F(alpha = 1, lambda = 0, rho) $F - ref_F(alpha = 1, lambda = 0,
# rho = 0) $F (both share the identical survival term). 10 random points, rho log-uniform on [1e-2, 1e3].
# p = 80, n = 50 with tied times. Pass: relative difference <= 1e-12.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S3_candidate_identity.txt")
cat("S3 identity with the frozen candidate |", format(Sys.time()), "\n")
C <- load_candidate()
s <- sim_data(13, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d\n", D$p, D$n, D$ne, s$n_tied))
TOL <- 1e-12
rs <- numeric(10); rr <- numeric(10)
for (i in 1:10) {
  set.seed(300 + i)
  W <- rand_W(D$p); H <- matrix(runif(3 * D$n, 0, 50), 3, D$n); b <- rnorm(3, 0, 1); rho <- 10^runif(1, -2, 3)
  sig <- palm_scores(W, D)$sigma
  th <- D$s0 * b / sig
  pr <- palm_parts(W, H, th, D, list(alpha = 0.5, lambda = 1, nu = 0.5, rho = rho, lambdaH = 1e-6))
  f1 <- C$ref_F(W, H, b, D$X, D$y, D$d, alpha = 1, lambda = 0, nu = 0.5, rho = rho, lambdaH = 1e-6, eps = 0)
  f0 <- C$ref_F(W, H, b, D$X, D$y, D$d, alpha = 1, lambda = 0, nu = 0.5, rho = 0, lambdaH = 1e-6, eps = 0)
  surv_c <- -(2 / D$ne) * f1$loglik; ridge_c <- f1$F - f0$F
  rs[i] <- abs(pr$surv - surv_c) / max(abs(pr$surv), abs(surv_c))
  rr[i] <- abs(pr$ridge - ridge_c) / max(abs(pr$ridge), abs(ridge_c))
  cat(sprintf("  point %2d: rho %.3e  surv %.15e vs %.15e (rel %.2e)  ridge %.15e vs %.15e (rel %.2e)\n",
              i, rho, pr$surv, surv_c, rs[i], pr$ridge, ridge_c, rr[i]))
}
cat("\n")
tally("survival term: theta_lib vs candidate ref_F", rs, TOL)
tally("rho term: theta_lib vs candidate ref_F", rr, TOL)
total("S3")
sink()
