# Independent KKT evaluator for the alpha-free theta subproblem (spec 2026-10-05 revision 4, sections 2 and 4):
#   f(theta) = -(2/n_e) l(Zc theta / s0) + lambda (nu |theta|_1 + (1-nu)/2 |theta|^2) + (rho/2) |P W theta|^2
# Written WITHOUT theta_lib.R: own s0, own scores, own Breslow gradient (dense risk-set matrix, O(n^2)).
# Pattern of candidate/kkt_oracle.R. Returns ||KKT||_inf with KKT_j = max(|g_j| - lambda nu, 0) if theta_j = 0,
# |g_j + lambda nu sign(theta_j)| otherwise, g the gradient of the smooth part (Cox term plus both quadratics).
kkt_oracle_theta_setup <- function(X, y, d) {
  list(X = X, d = as.numeric(d), ne = sum(d), n = ncol(X), atrisk = outer(y, y, function(yi, yj) yj >= yi),
       s0 = sqrt(mean(apply(X, 1, function(r) sum((r - mean(r))^2) / (length(r) - 1)))))
}
kkt_oracle_theta_grad <- function(theta, W, O, lambda, nu, rho) {
  n <- O$n; k <- ncol(W)
  sc <- t(O$X) %*% W
  cen <- sc - matrix(colMeans(sc), n, k, byrow = TRUE)           # Zc = C X^T W, no sigma
  eta <- as.vector(cen %*% theta) / O$s0
  E <- matrix(eta, n, n, byrow = TRUE); E[!O$atrisk] <- -Inf      # E[i, j] = eta_j for j in the risk set of i
  mi <- apply(E, 1, max); lS <- mi + log(rowSums(exp(E - mi)))    # per-risk-set log-sum-exp
  P <- exp(E - lS)                                                # P[i, j] = exp(eta_j) / S_i, 0 outside risk set
  dl <- O$d - as.vector(crossprod(P, O$d))                        # d l / d eta
  Wc <- W - matrix(colMeans(W), nrow(W), k, byrow = TRUE)         # P W
  -(2 / O$ne) * as.vector(crossprod(cen, dl)) / O$s0 + lambda * (1 - nu) * theta +
    rho * as.vector(crossprod(Wc, Wc %*% theta))
}
kkt_oracle_theta <- function(theta, W, O, lambda, nu, rho) {
  g <- kkt_oracle_theta_grad(theta, W, O, lambda, nu, rho); l1 <- lambda * nu
  max(ifelse(theta == 0, pmax(abs(g) - l1, 0), abs(g + l1 * sign(theta))))
}
