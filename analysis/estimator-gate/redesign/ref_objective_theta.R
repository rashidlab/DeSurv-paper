# Independent evaluator of the redesigned objective F (spec 2026-10-05 revision 4, section 2), written
# WITHOUT theta_lib.R. Pattern of candidate/ref_objective.R. Cox log partial likelihood from survival::coxph at
# fixed coefficients beta = theta / s0 on the unscaled centered scores Zc (Breslow ties, no iterations).
ref_F_theta <- function(W, H, theta, X, y, d, alpha, lambda, nu, rho, lambdaH) {
  n <- ncol(X); p <- nrow(X); k <- ncol(W)
  s0 <- sqrt(mean(apply(X, 1, var)))
  scores <- t(X) %*% W
  Zc <- scores - matrix(colMeans(scores), n, k, byrow = TRUE)
  fit <- survival::coxph(survival::Surv(y, d) ~ Zc, init = theta / s0, ties = "breslow",
                         control = survival::coxph.control(iter.max = 0))
  ll <- fit$loglik[1]
  Pw <- W - matrix(colMeans(W), p, k, byrow = TRUE)
  v <- Pw %*% theta
  X2 <- sum(X^2)
  recon <- sum((X - W %*% H)^2) / (2 * X2) + lambdaH * sum(H^2) / (2 * X2)
  surv <- -(2 / sum(d)) * ll; ridge <- rho / 2 * sum(v^2)
  sup <- surv + lambda * (nu * sum(abs(theta)) + (1 - nu) / 2 * sum(theta^2)) + ridge
  list(F = (1 - alpha) * recon + alpha * sup, loglik = ll, surv = surv, ridge = ridge, s0 = s0, beta = theta / s0)
}
