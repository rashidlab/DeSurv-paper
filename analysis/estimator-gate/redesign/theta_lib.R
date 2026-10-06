# Redesigned DeSurv estimator (spec 2026-10-05 DRAFT revision 4, sections 1, 2, 4, 5, 6): theta
# parameterization, eta = Zc(W) theta / s0, no 1/sigma(W) and no eps anywhere. Block order H (one HALS sweep),
# W (projected gradient with backtracking), theta (FISTA to an independent KKT certificate). Base R.
# X is p x n (within-sample ranks), y times, d event indicators. Built from candidate/palm_lib.R (the failure
# record, never edited); see REPORT.md for the function-by-function breakage inventory.
# Requires kkt_oracle_theta.R (the independent theta-block KKT evaluator) sourced first.
if (!exists("kkt_oracle_theta")) stop("source kkt_oracle_theta.R before theta_lib.R", call. = FALSE)

# ---- data object --------------------------------------------------------------------------------
# Invalid input (spec section 2) fails deterministically: n < 2, n_e = 0, nonfinite X, y or d, ||X|| = 0, s0 = 0.
# s0 = sqrt(sb2), sb2 = mean over genes of the row variance (denominator n - 1), from X only, once per fit.
# lam1 = lambda_1(X C X^T), computed as the largest eigenvalue of the n x n matrix C X^T X C (same nonzero spectrum).
palm_data <- function(X, y, d) {
  if (!is.matrix(X) && !is.data.frame(X)) stop("palm_data: X must be a p x n matrix", call. = FALSE)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("palm_data: X must be numeric", call. = FALSE)
  storage.mode(X) <- "double"
  n <- ncol(X); p <- nrow(X)
  if (n < 2L) stop("palm_data: need n >= 2 samples (columns of X)", call. = FALSE)
  if (p < 1L) stop("palm_data: need p >= 1 genes (rows of X)", call. = FALSE)
  if (any(!is.finite(X))) stop("palm_data: X has non-finite entries (NA, NaN or Inf)", call. = FALSE)
  if (sum(X^2) == 0) stop("palm_data: ||X|| = 0", call. = FALSE)
  if (length(y) != n || length(d) != n) stop("palm_data: y and d must have length ncol(X)", call. = FALSE)
  if (!is.numeric(y) || any(!is.finite(y))) stop("palm_data: y has non-finite or non-numeric entries", call. = FALSE)
  if (!is.numeric(d) && !is.logical(d) || any(!is.finite(as.numeric(d)))) stop("palm_data: d has non-finite entries", call. = FALSE)
  d <- as.numeric(d)
  if (!all(d %in% c(0, 1))) stop("palm_data: d must be 0/1", call. = FALSE)
  if (sum(d) < 1) stop("palm_data: no events (n_e = 0)", call. = FALSE)
  sb2 <- mean(apply(X, 1, var)); s0 <- sqrt(sb2)
  if (!is.finite(s0) || s0 == 0) stop("palm_data: s0 = 0 (every gene is constant across samples)", call. = FALSE)
  if (is.null(rownames(X))) rownames(X) <- paste0("g", seq_len(p))
  if (anyDuplicated(rownames(X))) stop("palm_data: duplicated gene names", call. = FALSE)
  o <- order(y); ys <- y[o]
  gs <- seq_len(n); for (i in 2:n) if (ys[i] == ys[i - 1]) gs[i] <- gs[i - 1]
  ge <- seq_len(n); for (i in (n - 1):1) if (ys[i] == ys[i + 1]) ge[i] <- ge[i + 1]
  Xc <- X - rowMeans(X)                                                   # X C
  lam1 <- max(eigen(crossprod(Xc), symmetric = TRUE, only.values = TRUE)$values)
  list(X = X, y = as.numeric(y), d = d, p = p, n = n, ne = sum(d), X2 = sum(X^2), genes = rownames(X),
       sb2 = sb2, s0 = s0, lam1 = lam1, maxX = max(X), o = o, gs = gs, ge = ge,
       oracle = kkt_oracle_theta_setup(X, y, d))
}

# Breslow log partial likelihood and its gradient in eta (exact ties). Unchanged from the candidate.
breslow <- function(eta, D) {
  o <- D$o; es <- eta[o]; ds <- D$d[o]; m <- max(es); e <- exp(es - m)
  S <- rev(cumsum(rev(e)))[D$gs]
  if (min(S[ds > 0]) > 1e-200) {
    ll <- sum((ds * (es - m - log(S)))[ds > 0])
    h <- cumsum(ifelse(ds > 0, ds / S, 0))[D$ge]
    g <- numeric(length(eta)); g[o] <- ds - e * h
    if (is.finite(ll) && all(is.finite(g))) return(list(ll = ll, g = g))   # else fall through to log-sum-exp
  }
  n <- length(es); lsuf <- numeric(n); acc <- -Inf
  for (i in n:1) { mx <- max(acc, es[i]); acc <- mx + log(exp(acc - mx) + exp(es[i] - mx)); lsuf[i] <- acc }
  logS <- lsuf[D$gs]
  ll <- sum(ds * (es - logS))
  lh <- numeric(n); acc <- -Inf
  for (i in 1:n) {
    if (ds[i] > 0) { v <- -logS[i]; mx <- max(acc, v); acc <- mx + log(exp(acc - mx) + exp(v - mx)) }
    lh[i] <- acc
  }
  g <- numeric(length(eta)); g[o] <- ds - exp(es + lh[D$ge])
  list(ll = ll, g = g)
}

# ---- scores and objective -----------------------------------------------------------------------
# Zc = C X^T W (no sigma). sigma_j = sd(Zc_j) is returned for reporting and diagnostics only.
palm_scores <- function(W, D) {
  Z <- crossprod(D$X, W)
  Zc <- sweep(Z, 2, colMeans(Z))
  list(Zc = Zc, sigma = sqrt(colSums(Zc^2) / (D$n - 1)))
}

# All parts of F. Psi = smooth part (F without lambda*nu*||theta||_1 and without the indicators).
palm_parts <- function(W, H, theta, D, hp) {
  a <- hp$alpha; sc <- palm_scores(W, D)
  eta <- as.vector(sc$Zc %*% theta) / D$s0
  ll <- breslow(eta, D)$ll
  u <- as.vector(sweep(W, 2, colMeans(W)) %*% theta)              # P W theta
  nmf <- sum((D$X - W %*% H)^2) / (2 * D$X2)
  hreg <- hp$lambdaH * sum(H^2) / (2 * D$X2)                        # ruling 12
  surv <- -(2 / D$ne) * ll
  l1 <- hp$lambda * hp$nu * sum(abs(theta)); l2 <- hp$lambda * (1 - hp$nu) / 2 * sum(theta^2)
  ridge <- hp$rho / 2 * sum(u^2)
  Psi <- (1 - a) * (nmf + hreg) + a * (surv + l2 + ridge)
  list(F = Psi + a * l1, Psi = Psi, nmf = nmf, hreg = hreg, loglik = ll, surv = surv, l1 = l1, l2 = l2,
       ridge = ridge, beta = theta / D$s0, sigma = sc$sigma)
}
palm_F <- function(W, H, theta, D, hp) palm_parts(W, H, theta, D, hp)$F
palm_Psi <- function(W, H, theta, D, hp) palm_parts(W, H, theta, D, hp)$Psi

# Reported raw-score Cox coefficient: beta = theta / s0
palm_beta <- function(theta, D) theta / D$s0

# ---- gradients of Psi ---------------------------------------------------------------------------
palm_grad_H <- function(W, H, D, hp)
  (1 - hp$alpha) * (-crossprod(W, D$X - W %*% H) + hp$lambdaH * H) / D$X2

# Per-term W gradients (unweighted): reconstruction ||X-WH||^2/(2||X||^2), Cox -(2/ne) l(Zc theta/s0),
# rho-ridge (rho/2)||P W theta||^2. eta is linear in W: no chain rule through sigma.
grad_W_terms <- function(W, H, theta, D, hp) {
  sc <- palm_scores(W, D)
  g <- breslow(as.vector(sc$Zc %*% theta) / D$s0, D)$g
  Cg <- g - mean(g)
  u <- as.vector(sweep(W, 2, colMeans(W)) %*% theta)
  list(recon = -(D$X - W %*% H) %*% t(H) / D$X2,
       cox = outer(-(2 / D$ne) * as.vector(D$X %*% Cg) / D$s0, theta),
       ridge = outer(hp$rho * u, theta))
}
palm_grad_W <- function(W, H, theta, D, hp) {
  gt <- grad_W_terms(W, H, theta, D, hp)
  (1 - hp$alpha) * gt$recon + hp$alpha * (gt$cox + gt$ridge)
}

# Smooth gradient of the alpha-free theta subproblem; grad_theta Psi = alpha * this.
# f(theta) = -(2/ne) l(Zs theta) + lambda (1-nu)/2 ||theta||^2 + (rho/2) theta^T M theta, Zs = Zc / s0,
# M = (P W)^T (P W).
tsub_M <- function(W) { Wc <- sweep(W, 2, colMeans(W)); crossprod(Wc) }
tsub_terms <- function(theta, Zs, M, D, hp) {
  br <- breslow(as.vector(Zs %*% theta), D); Mt <- as.vector(M %*% theta)
  list(f_cox = -(2 / D$ne) * br$ll, g_cox = -(2 / D$ne) * as.vector(crossprod(Zs, br$g)),
       f_tridge = hp$lambda * (1 - hp$nu) / 2 * sum(theta^2), g_tridge = hp$lambda * (1 - hp$nu) * theta,
       f_rho = hp$rho / 2 * sum(theta * Mt), g_rho = hp$rho * Mt)
}
tsub_smooth <- function(theta, Zs, M, D, hp) {
  t <- tsub_terms(theta, Zs, M, D, hp)
  list(f = t$f_cox + t$f_tridge + t$f_rho, g = t$g_cox + t$g_tridge + t$g_rho)
}
tsub_setup <- function(W, D) list(Zs = palm_scores(W, D)$Zc / D$s0, M = tsub_M(W))
tsub_grad <- function(W, theta, D, hp) { ts <- tsub_setup(W, D); tsub_smooth(theta, ts$Zs, ts$M, D, hp)$g }
palm_grad_theta <- function(W, theta, D, hp) hp$alpha * tsub_grad(W, theta, D, hp)

# minimum-norm subgradient of the theta subproblem (smooth gradient g plus lambda*nu*d|theta|_1); unchanged form
kkt_theta <- function(theta, g, lam1) ifelse(theta != 0, g + lam1 * sign(theta), sign(g) * pmax(abs(g) - lam1, 0))

# ---- W-block Hessian-vector product (Euclidean, fixed H and theta) ---------------------------------
# Breslow Hessian-vector product, unchanged from the candidate: breslow_hvp(eta, D, v) = d g(eta)[v], the
# derivative of the log-likelihood gradient g (log-likelihood Hessian, negative semidefinite). The NEGATIVE
# log-likelihood HVP of the brief is Hl[v] = -breslow_hvp(eta, D, v) (sign checked by finite differences in S6).
lse2 <- function(a, b) { mx <- max(a, b); if (mx == -Inf) -Inf else mx + log(exp(a - mx) + exp(b - mx)) }
breslow_hvp <- function(eta, D, v, force_lse = FALSE) {
  o <- D$o; es <- eta[o]; ds <- D$d[o]; vs <- v[o]; m <- max(es); e <- exp(es - m)
  S <- rev(cumsum(rev(e)))[D$gs]
  if (!force_lse && min(S[ds > 0]) > 1e-200) {
    ll <- sum((ds * (es - m - log(S)))[ds > 0])
    h <- cumsum(ifelse(ds > 0, ds / S, 0))[D$ge]
    if (is.finite(ll) && all(is.finite(ds - e * h))) {
      T <- rev(cumsum(rev(e * vs)))[D$gs]
      kk <- cumsum(ifelse(ds > 0, ds * (T / S) / S, 0))[D$ge]   # (T/S)/S: S^2 may underflow where S does not
      out <- numeric(length(eta)); out[o] <- -e * vs * h + e * kk
      return(out)
    }
  }
  n <- length(es); lsuf <- numeric(n); lTp <- numeric(n); lTn <- numeric(n); acc <- -Inf; ap <- -Inf; an <- -Inf
  for (i in n:1) {
    acc <- lse2(acc, es[i]); lsuf[i] <- acc
    if (vs[i] > 0) ap <- lse2(ap, es[i] + log(vs[i])); if (vs[i] < 0) an <- lse2(an, es[i] + log(-vs[i]))
    lTp[i] <- ap; lTn[i] <- an
  }
  logS <- lsuf[D$gs]; lTp <- lTp[D$gs]; lTn <- lTn[D$gs]
  lh <- numeric(n); lkp <- numeric(n); lkn <- numeric(n); acc <- -Inf; ap <- -Inf; an <- -Inf
  for (i in 1:n) {
    if (ds[i] > 0) {
      acc <- lse2(acc, -logS[i]); ap <- lse2(ap, lTp[i] - 2 * logS[i]); an <- lse2(an, lTn[i] - 2 * logS[i])
    }
    lh[i] <- acc; lkp[i] <- ap; lkn[i] <- an
  }
  out <- numeric(length(eta))
  out[o] <- -vs * exp(es + lh[D$ge]) + exp(es + lkp[D$ge]) - exp(es + lkn[D$ge])
  out
}
# Per-term W-block Hessian applied to V (p x k), unweighted:
# recon V H H^T / X2; Cox (2/ne) X C Hl[C X^T V theta / s0] theta^T / s0; ridge rho P V theta theta^T.
hvp_W_terms <- function(W, H, theta, D, hp, V) {
  Zc <- palm_scores(W, D)$Zc
  eta <- as.vector(Zc %*% theta) / D$s0
  Y <- crossprod(D$X, V); Yc <- sweep(Y, 2, colMeans(Y))
  deta <- as.vector(Yc %*% theta) / D$s0
  hl <- -breslow_hvp(eta, D, deta); Chl <- hl - mean(hl)
  Vc <- sweep(V, 2, colMeans(V))
  list(recon = V %*% tcrossprod(H) / D$X2,
       cox = outer((2 / D$ne) * as.vector(D$X %*% Chl) / D$s0, theta),
       ridge = outer(hp$rho * as.vector(Vc %*% theta), theta))
}
palm_hvp_W <- function(W, H, theta, D, hp, V) {
  ht <- hvp_W_terms(W, H, theta, D, hp, V)
  (1 - hp$alpha) * ht$recon + hp$alpha * (ht$cox + ht$ridge)
}
# Sigma-free curvature bound (spec section 2):
# L_W <= (1 - alpha) ||H||_2^2 / X2 + alpha ||theta||^2 (lambda_1(X C X^T) / s0^2 + rho)
palm_LW <- function(H, theta, D, hp)
  (1 - hp$alpha) * norm(H, "2")^2 / D$X2 + hp$alpha * sum(theta^2) * (D$lam1 / D$s0^2 + hp$rho)

# Free face and tangent projection: unchanged from the candidate. Not used by the projected-gradient W step;
# kept as utilities (breakage inventory).
w_face <- function(W, GW) (W > 0) | (W == 0 & GW < 0)
face_tangent <- function(V, W, free) { V <- V * free; V - sweep(W, 2, colSums(W * V), "*") }

# ---- projection onto S_+^k (unchanged) -----------------------------------------------------------
proj_S <- function(V) {
  for (j in seq_len(ncol(V))) {
    v <- V[, j]; vp <- pmax(v, 0); nv <- sqrt(sum(vp^2))
    if (nv > 0) V[, j] <- vp / nv else { e <- numeric(length(v)); e[which.max(v)] <- 1; V[, j] <- e }
  }
  V
}

# ---- swappable definitions ----------------------------------------------------------------------
# resid_W: W stationarity residual (KKT on the nonnegative unit sphere, unchanged); resid_H: H residual with
# a_H = (1 - alpha)(1 + lambdaH)/X2 (ruling 12); step_init: W-step trial size (unchanged); theta_accept: whether
# theta+ replaces theta (the candidate's b_accept, renamed); score: risk score from a stored model (unchanged).
PALM_DEFS <- list(
  resid_W = function(W, GW, D, hp) {
    mu <- colSums(W * GW)                                   # w_j^T G_j
    R <- ifelse(W > 0, GW - sweep(W, 2, mu, "*"), pmin(GW, 0))
    sqrt(sum(R^2))
  },
  resid_H = function(H, GH, W, D, hp) {
    aH <- (1 - hp$alpha) * (1 + hp$lambdaH) / D$X2
    sqrt(sum(pmin(H, GH / aH)^2)) / max(1, sqrt(sum(H^2)))
  },
  step_init = function(t_prev) min(2 * t_prev, 2),
  theta_accept = function(F_new, F_old) F_new <= F_old,
  score = function(Xnew, model) {
    G <- model$G
    if (!all(G %in% rownames(Xnew))) stop("palm_score: scored genes missing from the new data", call. = FALSE)
    as.vector(crossprod(Xnew[G, , drop = FALSE], model$v_tilde[G]))
  }
)
# The original (pre-amendment) residual, recorded as r_original and binding nothing (unchanged form).
resid_original <- function(W, H, GW, GH, rt)
  max(sqrt(sum((W - proj_S(W - GW))^2)), sqrt(sum(pmin(H, GH)^2)) / sqrt(sum(H^2)), rt)

# ---- H block: one HALS sweep (ruling 12) ---------------------------------------------------------
# h_j = max(0, w_j^T R_j / (1 + lambdaH)), R_j = X - sum_{l != j} w_l h_l, rows j = 1..k in order.
hals_sweep <- function(W, H, D, hp) {
  k <- ncol(W); WtX <- crossprod(W, D$X); WtW <- crossprod(W)
  den <- 1 + hp$lambdaH
  for (j in seq_len(k)) {
    num <- WtX[j, ] - as.vector(WtW[j, -j, drop = FALSE] %*% H[-j, , drop = FALSE])
    H[j, ] <- pmax(0, num / den)
  }
  H
}

# ---- W block: projected gradient step with backtracking (unchanged mechanics) -------------------
# t = min(2 t_prev, 2); halve until Psi(W+) <= Psi(W) + <G, W+-W> + ||W+-W||^2/(4t). Below t_min the
# step is not taken and the start is a numerical failure.
w_step <- function(W, H, theta, D, hp, t_prev, t_min = 1e-20, defs = PALM_DEFS) {
  G <- palm_grad_W(W, H, theta, D, hp); P0 <- palm_Psi(W, H, theta, D, hp)
  t <- defs$step_init(t_prev); nbt <- 0L
  repeat {
    Wn <- proj_S(W - t * G); Dl <- Wn - W
    if (isTRUE(palm_Psi(Wn, H, theta, D, hp) <= P0 + sum(G * Dl) + sum(Dl^2) / (4 * t))) break
    t <- t / 2; nbt <- nbt + 1L
    if (t < t_min) return(list(W = W, t = t_prev, fail = TRUE, nbt = nbt))
  }
  list(W = Wn, t = t, fail = FALSE, nbt = nbt)
}

# ---- theta block: FISTA with adaptive restart, stop certified by the independent KKT evaluator ---
# Internal check (kkt_theta) proposes a stop; the stop happens only when kkt_oracle_theta() <= tol. Cap = failure.
solve_theta <- function(theta0, W, D, hp, L0 = 1, tol = 1e-10, maxit = 1e5) {
  ts <- tsub_setup(W, D); Zs <- ts$Zs; M <- ts$M
  lam1 <- hp$lambda * hp$nu
  orc <- function(th) kkt_oracle_theta(th, W, D$oracle, hp$lambda, hp$nu, hp$rho)
  soft <- function(v, s) sign(v) * pmax(abs(v) - s, 0)
  x <- theta0; fx <- tsub_smooth(x, Zs, M, D, hp)
  norc <- 0L
  if (max(abs(kkt_theta(x, fx$g, lam1))) <= tol) {
    norc <- 1L; ko <- orc(x)
    if (ko <= tol) return(list(theta = x, it = 0L, kkt = ko, L = L0, conv = TRUE, n_oracle = norc))
  }
  y <- x; fy <- fx; tk <- 1; L <- L0; ko <- NA_real_
  for (it in seq_len(maxit)) {
    repeat {
      xn <- soft(y - fy$g / L, lam1 / L); fxn <- tsub_smooth(xn, Zs, M, D, hp)
      if (isTRUE(sum((fxn$g - fy$g) * (xn - y)) <= L / 2 * sum((xn - y)^2))) break   # non-finite: shrink step
      L <- 2 * L
      if (!is.finite(L)) return(list(theta = x, it = it, kkt = NA_real_, L = L, conv = FALSE, n_oracle = norc))
    }
    if (sum((y - xn) * (xn - x)) > 0) tk <- 1                  # gradient-mapping restart
    tn <- (1 + sqrt(1 + 4 * tk^2)) / 2
    y <- xn + (tk - 1) / tn * (xn - x); tk <- tn
    x <- xn; fx <- fxn
    if (max(abs(kkt_theta(x, fx$g, lam1))) <= tol) {
      norc <- norc + 1L; ko <- orc(x)
      if (isTRUE(ko <= tol)) return(list(theta = x, it = it, kkt = ko, L = L, conv = TRUE, n_oracle = norc))
    }
    fy <- if (identical(y, x)) fx else tsub_smooth(y, Zs, M, D, hp)
  }
  list(theta = x, it = maxit, kkt = if (is.na(ko)) orc(x) else ko, L = L, conv = FALSE, n_oracle = norc)
}

# ---- residuals (spec section 6) -----------------------------------------------------------------
palm_resid <- function(W, H, theta, D, hp, defs = PALM_DEFS) {
  GW <- palm_grad_W(W, H, theta, D, hp); GH <- palm_grad_H(W, H, D, hp)
  rW <- defs$resid_W(W, GW, D, hp)
  rH <- defs$resid_H(H, GH, W, D, hp)
  rt <- kkt_oracle_theta(theta, W, D$oracle, hp$lambda, hp$nu, hp$rho)
  c(r = max(rW, rH, rt), rW = rW, rH = rH, rtheta = rt, r_original = resid_original(W, H, GW, GH, rt))
}

# ---- initialization (spec section 5) ------------------------------------------------------------
# DeSurv 1.0.1 warmstart draw for start s of a fit with seed S: set.seed(S + s - 1); W0 then H0 from
# runif(0, max X); H0 is drawn so the stream matches and then discarded. Unchanged.
palm_draw <- function(D, k, seed, s) {
  set.seed(as.integer(seed + s - 1L))
  W0 <- matrix(runif(D$p * k, 0, D$maxX), D$p, k)
  H0 <- matrix(runif(D$n * k, 0, D$maxX), k, D$n)
  W0
}
palm_init_H <- function(W, D, hp, nsweep = 20L) {
  H <- matrix(0, ncol(W), D$n)
  for (i in seq_len(nsweep)) H <- hals_sweep(W, H, D, hp)
  H
}

# ---- diagnostics recording (qualification prespec section 6; binds nothing) ----------------------
# Per column j: sigma_j = sd(Zc_j), r_j = ||P w_j||, g_j = sigma_j / r_j, theta_j.
palm_colstats <- function(W, theta, D) {
  sig <- palm_scores(W, D)$sigma; r <- sqrt(colSums(sweep(W, 2, colMeans(W))^2))
  c(setNames(sig, paste0("sigma", seq_along(sig))), setNames(r, paste0("r", seq_along(r))),
    setNames(sig / r, paste0("g", seq_along(r))), setNames(theta, paste0("theta", seq_along(theta))))
}
PALM_CHECKPOINTS <- c(10L, 50L, 100L, 500L, 3750L, 4375L)   # 4375: Q3 window midpoint (ruling 14)
PALM_Q4_SWEEPS <- c(1L, 10L, 100L)   # (W, H, theta) for the Q4 curvature check (prespec section 7); recording only

# ---- one start ----------------------------------------------------------------------------------
# Status: converged (r < tol), nonconverged (maxit sweeps reached; finite, no failure) or failure (nonfinite
# value, W step floor, FISTA cap, or any F increase F_new > F_old + 1e-13 max(1, |F_old|) after a block).
# The start ends at the first of these. trace = TRUE records F after each block. record = TRUE records the
# prespec section 6 diagnostics: column statistics every 10th sweep, (W, theta) at PALM_CHECKPOINTS and the
# endpoint (a start that stops earlier stores its endpoint for later checkpoints). F increases are always recorded.
# max_secs is a test-harness wall-clock budget (default Inf, not part of the spec).
palm_fit_one <- function(D, hp, W, H, theta, maxit = 5000L, tol = 1e-5, trace = FALSE, record = FALSE,
                         defs = PALM_DEFS, max_secs = Inf) {
  inc_tol <- function(Fo) 1e-13 * max(1, abs(Fo))
  ninc <- c(H = 0L, W = 0L, theta = 0L); maxinc <- c(H = 0, W = 0, theta = 0)
  t_prev <- 1; Lt <- 1; nbt <- 0L; tits <- 0L; nrej <- 0L; norc <- 0L; nnull <- 0L
  status <- NA_character_; reason <- ""
  fail <- function(why) { status <<- "failure"; reason <<- why }
  Fcur <- palm_F(W, H, theta, D, hp); F0 <- Fcur
  r <- c(r = NA_real_, rW = NA_real_, rH = NA_real_, rtheta = NA_real_, r_original = NA_real_)
  if (!is.finite(Fcur)) fail("nonfinite F at start")
  tr <- if (trace) vector("list", maxit) else NULL
  incs <- list()
  rec_cols <- if (record) list() else NULL; rec_ck <- if (record) list() else NULL; rec_q4 <- if (record) list() else NULL
  inc <- function(Fn, Fo, blk) {
    di <- Fn - Fo
    if (is.finite(di) && di > maxinc[blk]) maxinc[blk] <<- di
    if (is.finite(di) && di > inc_tol(Fo)) {
      ninc[blk] <<- ninc[blk] + 1L
      incs[[length(incs) + 1L]] <<- data.frame(sweep = it, block = blk, F_old = Fo, F_new = Fn, increase = di)
      fail(paste("F increase in", blk, "block"))
    }
  }
  it <- 0L; t0 <- proc.time()[["elapsed"]]
  while (is.na(status) && it < maxit) {
    if (proc.time()[["elapsed"]] - t0 > max_secs) { status <- "incomplete"; reason <- "test time budget"; break }
    it <- it + 1L; Fit <- Fcur; FW <- NA_real_; Ft <- NA_real_; ws <- list(nbt = 0L, fail = FALSE); stit <- 0L; rej <- FALSE
    H <- hals_sweep(W, H, D, hp); FH <- palm_F(W, H, theta, D, hp)
    if (!is.finite(FH)) { fail("nonfinite F after H block"); break }
    inc(FH, Fcur, "H"); Fcur <- FH
    if (is.na(status)) {
      ws <- w_step(W, H, theta, D, hp, t_prev, defs = defs); nbt <- nbt + ws$nbt
      if (ws$fail) fail("W step floor") else {
        if (identical(ws$W, W)) nnull <- nnull + 1L                      # accepted null step (recorded only)
        W <- ws$W; t_prev <- ws$t; FW <- palm_F(W, H, theta, D, hp)
        if (!is.finite(FW)) fail("nonfinite F after W block") else { inc(FW, Fcur, "W"); Fcur <- FW }
      }
    }
    if (is.na(status)) {
      st <- solve_theta(theta, W, D, hp, L0 = max(Lt / 4, 1e-8)); tits <- tits + st$it; stit <- st$it; norc <- norc + st$n_oracle
      if (!st$conv) fail("FISTA cap or nonfinite in theta block") else {
        Lt <- st$L; Fn <- palm_F(W, H, st$theta, D, hp)
        if (isTRUE(defs$theta_accept(Fn, Fcur))) theta <- st$theta else { nrej <- nrej + 1L; rej <- TRUE }
        Ft <- palm_F(W, H, theta, D, hp)
        if (!is.finite(Ft)) fail("nonfinite F after theta block") else { inc(Ft, Fcur, "theta"); Fcur <- Ft }
      }
    }
    if (is.na(status)) {   # finiteness checked on binding residuals only; r_original binds nothing
      r <- palm_resid(W, H, theta, D, hp, defs)
      if (any(!is.finite(r[c("r", "rW", "rH", "rtheta")]))) fail("nonfinite residual") else if (r[["r"]] < tol) status <- "converged"
    }
    if (record) {
      if (it %% 10L == 0L) rec_cols[[length(rec_cols) + 1L]] <- c(sweep = it, palm_colstats(W, theta, D))
      if (it %in% PALM_CHECKPOINTS) rec_ck[[as.character(it)]] <- list(W = W, theta = theta, sweep = it)
      if (it %in% PALM_Q4_SWEEPS) rec_q4[[as.character(it)]] <- list(W = W, H = H, theta = theta, sweep = it)
    }
    if (trace) tr[[it]] <- c(iter = it, Fstart = Fit, FH = FH, FW = FW, Ftheta = Ft, r, t = t_prev, nbt = ws$nbt,
                             theta_it = stit, theta_rejected = rej, failed = identical(status, "failure"))
  }
  if (is.na(status)) status <- "nonconverged"
  if (record) {
    for (ck in PALM_CHECKPOINTS) if (is.null(rec_ck[[as.character(ck)]])) rec_ck[[as.character(ck)]] <- list(W = W, theta = theta, sweep = it)
    rec_ck[["end"]] <- list(W = W, theta = theta, sweep = it)
    for (q in PALM_Q4_SWEEPS) if (is.null(rec_q4[[as.character(q)]])) rec_q4[[as.character(q)]] <- list(W = W, H = H, theta = theta, sweep = it)
    rec_q4[["end"]] <- list(W = W, H = H, theta = theta, sweep = it)
  }
  secs <- proc.time()[["elapsed"]] - t0
  pr <- palm_parts(W, H, theta, D, hp)
  list(W = W, H = H, theta = theta, beta = pr$beta, sigma = pr$sigma, F = Fcur, F0 = F0, iter = it, status = status,
       reason = reason, conv = identical(status, "converged"), r = r, ninc = ninc, maxinc = maxinc,
       f_increases = if (length(incs)) do.call(rbind, incs) else NULL,
       theta_rejected = nrej, theta_iters = tits, n_oracle = norc, w_backtracks = nbt, w_null_steps = nnull,
       t_last = t_prev, secs = secs,
       record = if (record) list(colstats = if (length(rec_cols)) do.call(rbind, rec_cols) else NULL, checkpoints = rec_ck, q4 = rec_q4) else NULL,
       trace = if (trace) do.call(rbind, tr[seq_len(it)]) else NULL)
}

# ---- fit: 30 starts, status rule ---------------------------------------------------------------
# Returns the lowest-F converged start; else the lowest-F nonconverged start (fit status "nonconverged");
# else fit status "failed" with best = NULL. hp: alpha, lambda, nu, rho (lambdaH defaults to 1e-6).
pick_start <- function(Fs, status) {
  for (st in c("converged", "nonconverged")) {
    ok <- which(status == st)
    if (length(ok)) return(list(i = ok[which.min(Fs[ok])], status = st))
  }
  list(i = NA_integer_, status = "failed")
}
palm_fit <- function(D, hp, k = 3L, seed = 123L, nstarts = 30L, maxit = 5000L, tol = 1e-5,
                     ncores = 1L, keep_starts = FALSE, defs = PALM_DEFS, trace = FALSE, record = FALSE, max_secs = Inf) {
  if (is.null(hp$lambdaH)) hp$lambdaH <- 1e-6
  Ws <- lapply(seq_len(nstarts), function(s) proj_S(palm_draw(D, k, seed, s)))
  one <- function(s) {
    W <- Ws[[s]]; H <- palm_init_H(W, D, hp)
    f <- tryCatch(palm_fit_one(D, hp, W, H, rep(0, k), maxit = maxit, tol = tol, defs = defs, trace = trace,
                               record = record, max_secs = max_secs),
                  error = function(e) list(status = "failure", reason = paste("error:", conditionMessage(e)), F = NA_real_,
                                           iter = NA_integer_, r = c(r = NA, rW = NA, rH = NA, rtheta = NA, r_original = NA),
                                           ninc = c(H = NA, W = NA, theta = NA), theta_rejected = NA, w_backtracks = NA,
                                           theta_iters = NA, secs = NA))
    f$start <- s; f
  }
  fits <- if (ncores > 1L) parallel::mclapply(seq_len(nstarts), one, mc.cores = ncores) else lapply(seq_len(nstarts), one)
  if (!all(vapply(fits, is.list, TRUE))) stop("palm_fit: a worker returned no result", call. = FALSE)
  g <- function(nm, typ) vapply(fits, function(f) { v <- f[[nm]]; if (is.null(v)) typ[NA] else v }, typ)
  rr <- function(nm) vapply(fits, function(f) unname(f$r[nm]), 0)
  summ <- data.frame(start = seq_len(nstarts), status = g("status", ""), reason = g("reason", ""), F = g("F", 0),
                     iter = g("iter", 0L), r = rr("r"), rW = rr("rW"), rH = rr("rH"), rtheta = rr("rtheta"),
                     r_original = rr("r_original"), ninc = vapply(fits, function(f) sum(f$ninc), 0),
                     theta_rejected = g("theta_rejected", 0), w_backtracks = g("w_backtracks", 0),
                     w_null_steps = g("w_null_steps", 0), theta_iters = g("theta_iters", 0), secs = g("secs", 0),
                     stringsAsFactors = FALSE)
  pk <- pick_start(summ$F, summ$status)
  list(status = pk$status, best = if (is.na(pk$i)) NULL else fits[[pk$i]], best_start = pk$i, starts = summ,
       s0 = D$s0, lam1 = D$lam1, hp = hp, seed = seed, fits = if (keep_starts) fits else NULL)
}

# ---- stored model and scoring (spec section 2, prediction vector) -------------------------------
# beta = theta / s0, v = W beta, v~ = v - mean(v) over the training gene universe; G = unique union of
# desurv_get_top_genes(W, ntop)$top_genes. Centering happens before the support restriction; v~ is
# never recentered on another dataset's genes. score = X_G^T v~_G.
palm_model <- function(fit, D, ntop, top_genes_fun = DeSurv::desurv_get_top_genes) {
  if (is.null(fit$best)) stop("palm_model: fit has no returned start (status failed)", call. = FALSE)
  W <- fit$best$W; rownames(W) <- D$genes
  beta <- fit$best$theta / D$s0
  v <- as.vector(W %*% beta); names(v) <- D$genes
  vt <- v - mean(v)
  G <- unique(unlist(top_genes_fun(W, ntop)$top_genes))
  list(genes = D$genes, v_tilde = vt, G = G, ntop = ntop, W = W, theta = fit$best$theta, beta = beta,
       s0 = D$s0, sigma = fit$best$sigma, hp = fit$hp, fit_status = fit$status)
}
palm_score <- function(Xnew, model, defs = PALM_DEFS) defs$score(Xnew, model)

# ---- fold-local preprocessing (unchanged) --------------------------------------------------------
palm_rank <- function(Xraw) { R <- apply(Xraw, 2, rank, ties.method = "average"); dimnames(R) <- dimnames(Xraw); R }
palm_fold_data <- function(Xraw, y, d, train) palm_data(palm_rank(Xraw[, train, drop = FALSE]), y[train], d[train])
