# Independent KKT evaluator for the alpha-free b subproblem (spec section 4, amendment 2):
#   -(2/n_e) l(Zs b) + lambda (nu |b|_1 + (1-nu)/2 |b|^2) + (rho sigma_bar^2 / 2) |P W diag(1/sigma) b|^2
# Written without palm_lib.R: own scores, own Breslow gradient (dense risk-set matrix, O(n^2)).
# Returns ||KKT||_inf with KKT_j = max(|g_j| - lambda nu, 0) if b_j = 0, |g_j + lambda nu sign(b_j)| otherwise,
# g the gradient of the smooth part (Cox term plus both quadratic terms).
kkt_oracle_setup <- function(X, y, d) {
  list(X = X, d = as.numeric(d), ne = sum(d), n = ncol(X), atrisk = outer(y, y, function(yi, yj) yj >= yi),
       sb2 = mean(apply(X, 1, function(r) sum((r - mean(r))^2) / (length(r) - 1))))
}
kkt_oracle_grad <- function(b, W, O, eps, lambda, nu, rho) {
  n <- O$n; k <- ncol(W)
  sc <- t(O$X) %*% W
  cen <- sc - matrix(colMeans(sc), n, k, byrow = TRUE)
  sig <- sqrt(colSums(cen * cen) / (n - 1) + eps^2)
  Zs <- cen / matrix(sig, n, k, byrow = TRUE)
  eta <- as.vector(Zs %*% b); w <- exp(eta - max(eta))
  Ssum <- as.vector(O$atrisk %*% w)                         # risk-set sums at each sample's time
  dl <- O$d - w * as.vector(crossprod(O$atrisk, O$d / Ssum))  # d l / d eta
  Wc <- W - matrix(colMeans(W), nrow(W), k, byrow = TRUE)
  v <- Wc %*% (b / sig)
  -(2 / O$ne) * as.vector(crossprod(Zs, dl)) + lambda * (1 - nu) * b + rho * O$sb2 * as.vector(crossprod(Wc, v)) / sig
}
kkt_oracle <- function(b, W, O, eps, lambda, nu, rho) {
  g <- kkt_oracle_grad(b, W, O, eps, lambda, nu, rho); l1 <- lambda * nu
  max(ifelse(b == 0, pmax(abs(g) - l1, 0), abs(g + l1 * sign(b))))
}
