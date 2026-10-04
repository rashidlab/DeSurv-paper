# Independent base-R implementation of the stated DeSurv objective (Methods Eq. 1 /
# SI L) and stationarity checks. No DeSurv package code is used here.

breslow <- function(eta, y, d) {
  # Breslow log partial likelihood and its gradient wrt eta (ties handled exactly)
  R <- outer(y, y, function(ti, tj) tj >= ti)     # R[i, j] = 1 if j at risk at t_i
  m <- max(eta); e <- exp(eta - m)
  S <- as.vector(R %*% e)                          # S_i (scaled)
  ll <- sum(d * (eta - m - log(S)))
  g <- d - e * as.vector(crossprod(R, d / S))      # dl/deta_m
  list(ll = ll, g = g)
}

obj_parts <- function(W, H, beta, X, y, d, hp) {
  n <- ncol(X); p <- nrow(X); k <- ncol(W); ne <- sum(d)
  Xnorm <- sum(X^2)
  Z <- crossprod(X, W); s <- apply(Z, 2, sd)
  b <- as.vector(beta) * s
  P <- (1 - hp$nu) / 2 * sum(b^2) + hp$nu * sum(abs(b))
  bl <- breslow(as.vector(Z %*% beta), y, d)
  nmf <- sum((X - W %*% H)^2) / (2 * Xnorm)
  surv <- 2 / ne * bl$ll - hp$lambda * P
  L <- (1 - hp$alpha) * (nmf + hp$lambdaH * sum(H^2) / (2 * n * k)) - hp$alpha * surv +
       hp$lambdaW * sum(W^2) / (2 * p * k)
  list(L = L, nmf = nmf, surv = surv, P = P, ll = bl$ll, b = b, s = s)
}
objL <- function(W, H, beta, X, y, d, hp) obj_parts(W, H, beta, X, y, d, hp)$L

gradW <- function(W, H, beta, X, y, d, hp) {
  # full gradient of L wrt W, including the dependence of sd(Z) in the beta penalty
  n <- ncol(X); p <- nrow(X); k <- ncol(W); ne <- sum(d); Xnorm <- sum(X^2)
  Z <- crossprod(X, W); s <- apply(Z, 2, sd); beta <- as.vector(beta)
  g <- breslow(as.vector(Z %*% beta), y, d)$g
  G <- (1 - hp$alpha) / Xnorm * (W %*% H - X) %*% t(H)
  G <- G - hp$alpha * (2 / ne) * (X %*% g) %*% t(beta)
  Zc <- sweep(Z, 2, colMeans(Z))
  for (j in seq_len(k)) {
    dPds <- (1 - hp$nu) * beta[j]^2 * s[j] + hp$nu * abs(beta[j])
    G[, j] <- G[, j] + hp$alpha * hp$lambda * dPds * as.vector(X %*% Zc[, j]) / ((n - 1) * s[j])
  }
  G + hp$lambdaW * W / (p * k)
}

# beta subproblem in standardized coefficients b, W fixed:
#   f(b) = -(2/ne) l(Zs b) + lam [(1-nu)/2 |b|^2 + nu |b|_1],  Zs = Z / sd
beta_sub <- function(b, Zs, y, d, lam, nu, ne) {
  bl <- breslow(as.vector(Zs %*% b), y, d)
  list(f = -(2 / ne) * bl$ll + lam * ((1 - nu) / 2 * sum(b^2) + nu * sum(abs(b))),
       gs = -(2 / ne) * as.vector(crossprod(Zs, bl$g)) + lam * (1 - nu) * b)  # smooth-part gradient
}
kkt <- function(b, Zs, y, d, lam, nu, ne) {
  gs <- beta_sub(b, Zs, y, d, lam, nu, ne)$gs
  ifelse(b != 0, gs + lam * nu * sign(b), pmax(abs(gs) - lam * nu, 0))  # 0 at optimum
}
solve_b <- function(b0, Zs, y, d, lam, nu, ne, tol = 1e-13, maxit = 1e5) {
  # FISTA with backtracking; exact optimum of the convex beta subproblem
  soft <- function(v, t) sign(v) * pmax(abs(v) - t, 0)
  fs <- function(b) beta_sub(b, Zs, y, d, lam, nu, ne)$f - lam * nu * sum(abs(b))
  b <- b0; v <- b0; tk <- 1; L <- 1
  for (it in seq_len(maxit)) {
    gv <- beta_sub(v, Zs, y, d, lam, nu, ne)$gs; fv <- fs(v)
    repeat {
      bn <- soft(v - gv / L, lam * nu / L)
      if (fs(bn) <= fv + sum(gv * (bn - v)) + L / 2 * sum((bn - v)^2) + 1e-15) break
      L <- L * 2
    }
    tn <- (1 + sqrt(1 + 4 * tk^2)) / 2
    v <- bn + (tk - 1) / tn * (bn - b)
    if (max(abs(bn - b)) < tol) { b <- bn; break }
    b <- bn; tk <- tn
  }
  b
}

cidx <- function(risk, y, d) survival::concordance(survival::Surv(y, d) ~ risk, reverse = TRUE)$concordance
