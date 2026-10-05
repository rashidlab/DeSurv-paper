# Independent evaluator of the frozen objective F (spec section 4), written without palm_lib.R.
# Cox log partial likelihood from survival::coxph at fixed coefficients (Breslow ties, no iterations).
ref_F <- function(W, H, b, X, y, d, alpha, lambda, nu, rho, lambdaH, eps) {
  n <- ncol(X); p <- nrow(X); k <- ncol(W)
  scores <- t(X) %*% W
  cen <- scores - matrix(colMeans(scores), n, k, byrow = TRUE)
  sig <- sqrt(apply(cen, 2, function(z) sum(z^2)) / (n - 1) + eps^2)
  Zs <- cen / matrix(sig, n, k, byrow = TRUE)
  fit <- survival::coxph(survival::Surv(y, d) ~ Zs, init = b, ties = "breslow",
                         control = survival::coxph.control(iter.max = 0))
  ll <- fit$loglik[1]
  Pw <- W - matrix(colMeans(W), p, k, byrow = TRUE)
  v <- Pw %*% diag(1 / sig, k) %*% b
  recon <- sum((X - W %*% H)^2) / (2 * sum(X^2)) + lambdaH * sum(H^2) / (2 * n * k)
  sup <- -(2 / sum(d)) * ll + lambda * (nu * sum(abs(b)) + (1 - nu) / 2 * sum(b^2)) +
    rho * mean(apply(X, 1, var)) / 2 * sum(v^2)
  list(F = (1 - alpha) * recon + alpha * sup, loglik = ll, beta = b / sig)
}
