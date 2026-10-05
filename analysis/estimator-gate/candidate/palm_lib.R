# Candidate DeSurv estimator (spec 2026-10-04 section 4, amendments 1 and 2): block prox-linear (PALM)
# reference implementation. Base R. X is p x n (within-sample ranks), y times, d event indicators.
# Definitions that a spec amendment may change are isolated in PALM_DEFS so they can be swapped without
# touching the optimizer. Requires kkt_oracle.R (the independent b-block KKT evaluator) sourced first.
if (!exists("kkt_oracle")) stop("source kkt_oracle.R before palm_lib.R", call. = FALSE)

# ---- data object --------------------------------------------------------------------------------
# Invalid input (spec section 4) fails deterministically: n < 2, n_e = 0, nonfinite X, y or d, ||X|| = 0.
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
  if (is.null(rownames(X))) rownames(X) <- paste0("g", seq_len(p))
  if (anyDuplicated(rownames(X))) stop("palm_data: duplicated gene names", call. = FALSE)
  o <- order(y); ys <- y[o]
  gs <- seq_len(n); for (i in 2:n) if (ys[i] == ys[i - 1]) gs[i] <- gs[i - 1]
  ge <- seq_len(n); for (i in (n - 1):1) if (ys[i] == ys[i + 1]) ge[i] <- ge[i + 1]
  list(X = X, y = as.numeric(y), d = d, p = p, n = n, ne = sum(d), X2 = sum(X^2), genes = rownames(X),
       sb2 = mean(apply(X, 1, var)), maxX = max(X), o = o, gs = gs, ge = ge, oracle = kkt_oracle_setup(X, y, d))
}

# Breslow log partial likelihood and its gradient in eta (exact ties). Fast path O(n log n) with one shift;
# when a risk-set sum underflows under the single shift (predictor contrasts beyond ~460), a log-sum-exp
# recursion over the sorted samples is used instead.
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
# Zc = C X^T W, sigma_j(W) = sqrt(||Zc_j||^2/(n-1) + eps^2)
palm_scores <- function(W, D, eps) {
  Z <- crossprod(D$X, W)
  Zc <- sweep(Z, 2, colMeans(Z))
  list(Zc = Zc, sigma = sqrt(colSums(Zc^2) / (D$n - 1) + eps^2))
}

# All parts of F. Psi = smooth part (F without lambda*nu*||b||_1 and without the indicators).
palm_parts <- function(W, H, b, D, hp) {
  a <- hp$alpha; sc <- palm_scores(W, D, hp$eps); beta <- b / sc$sigma
  eta <- as.vector(sc$Zc %*% beta)
  ll <- breslow(eta, D)$ll
  u <- as.vector(sweep(W, 2, colMeans(W)) %*% beta)               # P W diag(1/sigma) b
  nmf <- sum((D$X - W %*% H)^2) / (2 * D$X2)
  hreg <- hp$lambdaH * sum(H^2) / (2 * D$n * ncol(W))
  surv <- -(2 / D$ne) * ll
  l1 <- hp$lambda * hp$nu * sum(abs(b)); l2 <- hp$lambda * (1 - hp$nu) / 2 * sum(b^2)
  ridge <- hp$rho * D$sb2 / 2 * sum(u^2)
  Psi <- (1 - a) * (nmf + hreg) + a * (surv + l2 + ridge)
  list(F = Psi + a * l1, Psi = Psi, nmf = nmf, hreg = hreg, loglik = ll, surv = surv, l1 = l1, l2 = l2,
       ridge = ridge, beta = beta, sigma = sc$sigma)
}
palm_F <- function(W, H, b, D, hp) palm_parts(W, H, b, D, hp)$F
palm_Psi <- function(W, H, b, D, hp) palm_parts(W, H, b, D, hp)$Psi

# Reported coefficients: beta_j = b_j / sigma_j(W) (floored scale, never the raw sd)
palm_beta <- function(W, b, D, eps) b / palm_scores(W, D, eps)$sigma

# ---- gradients of Psi ---------------------------------------------------------------------------
palm_grad_H <- function(W, H, D, hp)
  (1 - hp$alpha) * (-crossprod(W, D$X - W %*% H) / D$X2 + hp$lambdaH * H / (D$n * ncol(W)))

# Per-term W gradients (unweighted): reconstruction ||X-WH||^2/(2||X||^2), Cox -(2/ne) l, rho-ridge
# (rho sb2/2)||P W diag(1/sigma) b||^2. With beta = b/sigma(W), a term phi(W, beta(W)) whose partials are
# dphi/dW|beta = q beta^T and dphi/dbeta = W^T q has gradient q beta^T - X Zc diag(c beta/((n-1) sigma^2)),
# c = W^T q (sigma chain rule).
grad_W_terms <- function(W, H, b, D, hp) {
  sc <- palm_scores(W, D, hp$eps); beta <- b / sc$sigma
  g <- breslow(as.vector(sc$Zc %*% beta), D)$g
  u <- as.vector(sweep(W, 2, colMeans(W)) %*% beta)
  XZ <- D$X %*% sc$Zc
  chain <- function(q) outer(q, beta) - sweep(XZ, 2, as.vector(crossprod(W, q)) * beta / ((D$n - 1) * sc$sigma^2), "*")
  list(recon = -(D$X - W %*% H) %*% t(H) / D$X2,
       cox = chain(-(2 / D$ne) * as.vector(D$X %*% g)),
       ridge = chain(hp$rho * D$sb2 * u))
}
palm_grad_W <- function(W, H, b, D, hp) {
  gt <- grad_W_terms(W, H, b, D, hp)
  (1 - hp$alpha) * gt$recon + hp$alpha * (gt$cox + gt$ridge)
}

# smooth b-gradient of the alpha-free b subproblem (I3); grad_b Psi = alpha * this
bsub_M <- function(W, sigma) { Wc <- sweep(W, 2, colMeans(W)); crossprod(Wc) / outer(sigma, sigma) }
bsub_terms <- function(b, Zs, M, D, hp) {
  br <- breslow(as.vector(Zs %*% b), D); Mb <- as.vector(M %*% b)
  list(f_cox = -(2 / D$ne) * br$ll, g_cox = -(2 / D$ne) * as.vector(crossprod(Zs, br$g)),
       f_bridge = hp$lambda * (1 - hp$nu) / 2 * sum(b^2), g_bridge = hp$lambda * (1 - hp$nu) * b,
       f_rho = hp$rho * D$sb2 / 2 * sum(b * Mb), g_rho = hp$rho * D$sb2 * Mb)
}
bsub_smooth <- function(b, Zs, M, D, hp) {
  t <- bsub_terms(b, Zs, M, D, hp)
  list(f = t$f_cox + t$f_bridge + t$f_rho, g = t$g_cox + t$g_bridge + t$g_rho)
}
bsub_setup <- function(W, D, hp) {
  sc <- palm_scores(W, D, hp$eps)
  list(Zs = sweep(sc$Zc, 2, sc$sigma, "/"), M = bsub_M(W, sc$sigma))
}
palm_grad_b <- function(W, b, D, hp) {
  bs <- bsub_setup(W, D, hp)
  hp$alpha * bsub_smooth(b, bs$Zs, bs$M, D, hp)$g
}

# minimum-norm subgradient of the b subproblem (smooth gradient g plus lambda*nu*d|b|_1)
kkt_b <- function(b, g, lam1) ifelse(b != 0, g + lam1 * sign(b), sign(g) * pmax(abs(g) - lam1, 0))

# ---- W-block Hessian-vector products (amendment 3, section 3; analytic, no finite differences) -----
# Breslow Hessian-vector product Hll[v] = d g(eta)[v] (exact ties, same sort and tie groups as breslow()).
# Hll[v]_l = -e_l v_l h_l + e_l sum_{events i: y_i <= y_l} d_i T_i / S_i^2, T_i = sum_{m: y_m >= y_i} e_m v_m.
# Fast path under exactly breslow()'s condition (one shift; e h and e T/S^2 are shift invariant); otherwise
# log-sum-exp, with T split into its positive and negative parts. force_lse = TRUE is for testing only.
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

# Directional derivatives along V (p x k) of the per-term W gradients of grad_W_terms (unweighted).
# With s_j = c_j beta_j/((n-1) sigma_j^2), G = q beta^T - X Zc diag(s), the derivative is
# dq beta^T + q dbeta^T - X Yc diag(s) - X Zc diag(ds), Yc = C X^T V, dsigma_j = <Zc_j, Yc_j>/((n-1) sigma_j),
# dbeta = -beta dsigma/sigma, dc = V^T q + W^T dq, ds = (dc beta + c dbeta)/((n-1) sigma^2) - 2 s dsigma/sigma.
hvp_W_terms <- function(W, H, b, D, hp, V) {
  sc <- palm_scores(W, D, hp$eps); Zc <- sc$Zc; sig <- sc$sigma; beta <- b / sig; n1 <- D$n - 1
  Y <- crossprod(D$X, V); Yc <- sweep(Y, 2, colMeans(Y))
  dsig <- colSums(Zc * Yc) / (n1 * sig); dbeta <- -beta * dsig / sig
  XZ <- D$X %*% Zc; XY <- D$X %*% Yc
  dchain <- function(q, dq) {
    cc <- as.vector(crossprod(W, q)); dc <- as.vector(crossprod(V, q)) + as.vector(crossprod(W, dq))
    s <- cc * beta / (n1 * sig^2); dsv <- (dc * beta + cc * dbeta) / (n1 * sig^2) - 2 * s * dsig / sig
    outer(dq, beta) + outer(q, dbeta) - sweep(XY, 2, s, "*") - sweep(XZ, 2, dsv, "*")
  }
  eta <- as.vector(Zc %*% beta); deta <- as.vector(Yc %*% beta + Zc %*% dbeta)
  qc <- -(2 / D$ne) * as.vector(D$X %*% breslow(eta, D)$g)
  dqc <- -(2 / D$ne) * as.vector(D$X %*% breslow_hvp(eta, D, deta))
  Wc <- sweep(W, 2, colMeans(W)); Vc <- sweep(V, 2, colMeans(V))
  qr <- hp$rho * D$sb2 * as.vector(Wc %*% beta); dqr <- hp$rho * D$sb2 * as.vector(Vc %*% beta + Wc %*% dbeta)
  list(recon = V %*% tcrossprod(H) / D$X2, cox = dchain(qc, dqc), ridge = dchain(qr, dqr))
}
palm_hvp_W <- function(W, H, b, D, hp, V) {
  ht <- hvp_W_terms(W, H, b, D, hp, V)
  (1 - hp$alpha) * ht$recon + hp$alpha * (ht$cox + ht$ridge)
}

# Free face: W > 0, or W = 0 with G_W < 0 (the two cases of resid_W). face_tangent: zero off the face, then
# remove the w_j component of each column (tangent to the unit sphere; ||Pface(G_W)||_F = resid_W).
w_face <- function(W, GW) (W > 0) | (W == 0 & GW < 0)
face_tangent <- function(V, W, free) { V <- V * free; V - sweep(W, 2, colSums(W * V), "*") }
# Riemannian Hessian on the tangent face: Hess[V] = Pface(d2Psi V) - V diag(mu), mu_j = <w_j, G_W[, j]>.
riem_hess_W <- function(W, H, b, D, hp, V, free, GW)
  face_tangent(palm_hvp_W(W, H, b, D, hp, V), W, free) - sweep(V, 2, colSums(W * GW), "*")

# ---- projection onto S_+^k ----------------------------------------------------------------------
proj_S <- function(V) {
  for (j in seq_len(ncol(V))) {
    v <- V[, j]; vp <- pmax(v, 0); nv <- sqrt(sum(vp^2))
    if (nv > 0) V[, j] <- vp / nv else { e <- numeric(length(v)); e[which.max(v)] <- 1; V[, j] <- e }
  }
  V
}

# ---- swappable definitions (spec amendment 2) ---------------------------------------------------
# resid_W: W stationarity residual (KKT on the nonnegative unit sphere); resid_H: H residual;
# step_init: W-step trial size; b_accept: whether b+ replaces b; score: risk score from a stored model.
PALM_DEFS <- list(
  resid_W = function(W, GW, D, hp) {
    mu <- colSums(W * GW)                                   # w_j^T G_j
    R <- ifelse(W > 0, GW - sweep(W, 2, mu, "*"), pmin(GW, 0))
    sqrt(sum(R^2))
  },
  resid_H = function(H, GH, W, D, hp) {
    aH <- (1 - hp$alpha) * (1 / D$X2 + hp$lambdaH / (D$n * ncol(W)))
    sqrt(sum(pmin(H, GH / aH)^2)) / max(1, sqrt(sum(H^2)))
  },
  step_init = function(t_prev) min(2 * t_prev, 2),
  b_accept = function(F_new, F_old) F_new <= F_old,
  score = function(Xnew, model) {
    G <- model$G
    if (!all(G %in% rownames(Xnew))) stop("palm_score: scored genes missing from the new data", call. = FALSE)
    as.vector(crossprod(Xnew[G, , drop = FALSE], model$v_tilde[G]))
  }
)
# The original (pre-amendment) residual, recorded as r_original and binding nothing.
resid_original <- function(W, H, GW, GH, rb)
  max(sqrt(sum((W - proj_S(W - GW))^2)), sqrt(sum(pmin(H, GH)^2)) / sqrt(sum(H^2)), rb)

# ---- H block: one HALS sweep (I6) ---------------------------------------------------------------
hals_sweep <- function(W, H, D, hp) {
  k <- ncol(W); WtX <- crossprod(W, D$X); WtW <- crossprod(W)
  den0 <- hp$lambdaH * D$X2 / (D$n * k)
  for (j in seq_len(k)) {
    num <- WtX[j, ] - as.vector(WtW[j, -j, drop = FALSE] %*% H[-j, , drop = FALSE])
    H[j, ] <- pmax(0, num / (WtW[j, j] + den0))
  }
  H
}

# ---- W block: projected gradient step with backtracking -----------------------------------------
# t = min(2 t_prev, 2); halve until Psi(W+) <= Psi(W) + <G, W+-W> + ||W+-W||^2/(4t). Below t_min the
# step is not taken and the start is a numerical failure.
w_step <- function(W, H, b, D, hp, t_prev, t_min = 1e-20, defs = PALM_DEFS) {
  G <- palm_grad_W(W, H, b, D, hp); P0 <- palm_Psi(W, H, b, D, hp)
  t <- defs$step_init(t_prev); nbt <- 0L
  repeat {
    Wn <- proj_S(W - t * G); Dl <- Wn - W
    if (isTRUE(palm_Psi(Wn, H, b, D, hp) <= P0 + sum(G * Dl) + sum(Dl^2) / (4 * t))) break
    t <- t / 2; nbt <- nbt + 1L
    if (t < t_min) return(list(W = W, t = t_prev, fail = TRUE, nbt = nbt))
  }
  list(W = Wn, t = t, fail = FALSE, nbt = nbt)
}

# ---- W block (amendment 3, section 3): active-face Newton-CG with a retracted Armijo search -------
# Retraction W(t) = N([W + t d]_+): negatives set to zero, columns renormalized. NULL when a column is
# annihilated (N undefined). For d in the tangent face this cannot happen in exact arithmetic (<w_j, w_j + t d_j> = 1).
w_retract <- function(W, d, t) {
  V <- pmax(W + t * d, 0); nv <- sqrt(colSums(V^2))
  if (any(nv == 0)) return(NULL)
  sweep(V, 2, nv, "/")
}
# Monotone Armijo search on the retracted path at fixed (H, b). Trials t = 1, 1/2, ..., 2^-60 (the trial after the
# 60th halving is the last). A column-annihilating trial is rejected and halved. Accept the first nonzero trial with
# slope = <G_W, W(t) - W> < 0, F(W(t)) <= F(W) + 1e-4 slope and F(W(t)) < F(W) (projected-search Armijo on the actual
# displacement). Fails with kind "no_change" when W(t) equals W in every entry, "halvings" after 60 halvings without
# acceptance, "nonfinite" on a nonfinite trial point or F (a numerical failure, not a line-search failure).
ncg_line_search <- function(W, d, GW, H, b, D, hp, F0, max_halv = 60L) {
  t <- 1; nh <- 0L; nann <- 0L
  repeat {
    Wt <- w_retract(W, d, t)
    if (is.null(Wt)) nann <- nann + 1L else {
      if (any(!is.finite(Wt))) return(list(W = W, ok = FALSE, kind = "nonfinite", t = t, nbt = nh, n_annih = nann))
      if (all(Wt == W)) return(list(W = W, ok = FALSE, kind = "no_change", t = t, nbt = nh, n_annih = nann))
      slope <- sum(GW * (Wt - W))
      if (slope < 0) {
        Ft <- palm_F(Wt, H, b, D, hp)
        if (!is.finite(Ft)) return(list(W = W, ok = FALSE, kind = "nonfinite", t = t, nbt = nh, n_annih = nann))
        if (Ft <= F0 + 1e-4 * slope && Ft < F0) return(list(W = Wt, ok = TRUE, kind = "accepted", t = t, nbt = nh, n_annih = nann))
      }
    }
    if (nh >= max_halv) return(list(W = W, ok = FALSE, kind = "halvings", t = t, nbt = nh, n_annih = nann))
    t <- t / 2; nh <- nh + 1L
  }
}
# One Newton-CG W step at fixed (H, b). r_W < tol at entry: no-op, exit "w_converged", 0 HVPs. Otherwise plain
# truncated CG on Hess d = -g over the tangent face (g = projected gradient, ||g|| = r_W), d0 = 0, no
# preconditioner, stopped by the forcing test on the recursively updated residual (PI ruling B3.1),
# ||r|| <= min(0.5, sqrt(||g||)) ||g||; by max_hvp HVPs in this solve; by the per-start HVP total (hvp_left; PI
# ruling B3.2: the accumulated direction is still searched and the sweep completes); or by nonpositive or flat
# curvature q_j = p'Hess p / p'p <= 1e-12 max(1, |q_ref|), q_ref the first finite positive q (j = 1: d = -g,
# nc_first_iter; j > 1: accumulated d, negative_curvature). Nonfinite q: numerical failure. exit: "w_converged",
# "newton", "line_search_failed" or "numerical_failure". GW and d are returned for diagnostics only.
w_step_ncg <- function(W, H, b, D, hp, hvp_left, tol = 1e-5, defs = PALM_DEFS, max_hvp = 100L) {
  GW <- palm_grad_W(W, H, b, D, hp); rW <- defs$resid_W(W, GW, D, hp)
  out <- list(W = W, exit = "w_converged", GW = GW, rW = rW, face = NA_integer_, cg_iter = 0L, hvp = 0L, cg_exit = "none",
              negative_curvature = FALSE, nc_first_iter = FALSE, t = NA_real_, nbt = 0L, n_annih = 0L, n_enter = 0L,
              n_leave = 0L, ls_kind = "", d = NULL, res_rec = NA_real_, forcing_thr = NA_real_)
  if (!is.finite(rW)) { out$exit <- "numerical_failure"; return(out) }
  if (rW < tol) return(out)                                    # convergence guard: W kept, no CG, no search
  free <- w_face(W, GW); g <- face_tangent(GW, W, free); ng <- sqrt(sum(g^2))
  thr <- min(0.5, sqrt(ng)) * ng
  d <- 0 * g; r <- g; p <- -g; rr <- sum(r^2); qref <- NA_real_; j <- 0L; cg_exit <- ""
  repeat {
    if (j >= max_hvp) { cg_exit <- "hvp_ceiling"; break }
    if (j >= hvp_left) { cg_exit <- "hvp_total_cap"; break }
    Hp <- riem_hess_W(W, H, b, D, hp, p, free, GW); j <- j + 1L
    pHp <- sum(p * Hp); q <- pHp / sum(p * p)
    if (!is.finite(q)) { out[c("exit", "face", "cg_iter", "hvp")] <- list("numerical_failure", sum(free), j, j); return(out) }
    if (is.na(qref) && q > 0) qref <- q
    if (q <= 1e-12 * max(1, abs(if (is.na(qref)) q else qref))) {
      if (j == 1L) { d <- -g; out$nc_first_iter <- TRUE } else out$negative_curvature <- TRUE
      cg_exit <- "curvature"; break
    }
    a <- rr / pHp; d <- d + a * p; r <- r + a * Hp; rrn <- sum(r^2)
    if (sqrt(rrn) <= thr) { rr <- rrn; cg_exit <- "forcing"; break }
    p <- -r + (rrn / rr) * p; rr <- rrn
  }
  out[c("face", "cg_iter", "hvp", "cg_exit", "d", "res_rec", "forcing_thr")] <- list(sum(free), j, j, cg_exit, d, sqrt(rr), thr)
  if (any(!is.finite(d))) { out$exit <- "numerical_failure"; return(out) }
  ls <- ncg_line_search(W, d, GW, H, b, D, hp, palm_F(W, H, b, D, hp))
  out[c("t", "nbt", "n_annih", "ls_kind")] <- list(if (ls$ok) ls$t else NA_real_, ls$nbt, ls$n_annih, ls$kind)
  if (ls$kind == "nonfinite") { out$exit <- "numerical_failure"; return(out) }
  if (!ls$ok) { out$exit <- "line_search_failed"; return(out) }
  out$n_enter <- sum(W == 0 & ls$W > 0); out$n_leave <- sum(W > 0 & ls$W == 0)   # face changes
  out$W <- ls$W; out$exit <- "newton"
  out
}

# ---- b block: FISTA with adaptive restart, stop certified by the independent KKT evaluator -------
# Internal check (kkt_b) proposes a stop; the stop happens only when kkt_oracle() <= tol. Cap = failure.
solve_b <- function(b0, W, D, hp, L0 = 1, tol = 1e-10, maxit = 1e5) {
  bs <- bsub_setup(W, D, hp); Zs <- bs$Zs; M <- bs$M
  lam1 <- hp$lambda * hp$nu
  orc <- function(b) kkt_oracle(b, W, D$oracle, hp$eps, hp$lambda, hp$nu, hp$rho)
  soft <- function(v, s) sign(v) * pmax(abs(v) - s, 0)
  x <- b0; fx <- bsub_smooth(x, Zs, M, D, hp)
  norc <- 0L
  if (max(abs(kkt_b(x, fx$g, lam1))) <= tol) {
    norc <- 1L; ko <- orc(x)
    if (ko <= tol) return(list(b = x, it = 0L, kkt = ko, L = L0, conv = TRUE, n_oracle = norc))
  }
  y <- x; fy <- fx; tk <- 1; L <- L0; ko <- NA_real_
  for (it in seq_len(maxit)) {
    repeat {
      xn <- soft(y - fy$g / L, lam1 / L); fxn <- bsub_smooth(xn, Zs, M, D, hp)
      # f convex: f(xn) - f(y) - <g(y), xn - y> <= <g(xn) - g(y), xn - y>, so this implies the descent
      # lemma inequality and stays free of cancellation in f near the solution (I10)
      if (isTRUE(sum((fxn$g - fy$g) * (xn - y)) <= L / 2 * sum((xn - y)^2))) break   # non-finite: shrink step
      L <- 2 * L
      if (!is.finite(L)) return(list(b = x, it = it, kkt = NA_real_, L = L, conv = FALSE, n_oracle = norc))
    }
    # gradient-mapping restart (O'Donoghue and Candes 2015): no comparison of objective values, which
    # near the solution differ by less than rounding and would reset momentum at every step (I16)
    if (sum((y - xn) * (xn - x)) > 0) tk <- 1
    tn <- (1 + sqrt(1 + 4 * tk^2)) / 2
    y <- xn + (tk - 1) / tn * (xn - x); tk <- tn
    x <- xn; fx <- fxn
    if (max(abs(kkt_b(x, fx$g, lam1))) <= tol) {
      norc <- norc + 1L; ko <- orc(x)
      if (isTRUE(ko <= tol)) return(list(b = x, it = it, kkt = ko, L = L, conv = TRUE, n_oracle = norc))
    }
    fy <- if (identical(y, x)) fx else bsub_smooth(y, Zs, M, D, hp)
  }
  list(b = x, it = maxit, kkt = if (is.na(ko)) orc(x) else ko, L = L, conv = FALSE, n_oracle = norc)
}

# ---- residuals (spec section 4, amendment 2) ----------------------------------------------------
palm_resid <- function(W, H, b, D, hp, defs = PALM_DEFS) {
  GW <- palm_grad_W(W, H, b, D, hp); GH <- palm_grad_H(W, H, D, hp)
  rW <- defs$resid_W(W, GW, D, hp)
  rH <- defs$resid_H(H, GH, W, D, hp)
  rb <- kkt_oracle(b, W, D$oracle, hp$eps, hp$lambda, hp$nu, hp$rho)
  c(r = max(rW, rH, rb), rW = rW, rH = rH, rb = rb, r_original = resid_original(W, H, GW, GH, rb))
}

# ---- initialization ----------------------------------------------------------------------------
# DeSurv 1.0.1 warmstart draw for start s of a fit with seed S: set.seed(S + s - 1); W0 then H0 from
# runif(0, max X); H0 is drawn so the stream matches and then discarded.
palm_draw <- function(D, k, seed, s) {
  set.seed(as.integer(seed + s - 1L))
  W0 <- matrix(runif(D$p * k, 0, D$maxX), D$p, k)
  H0 <- matrix(runif(D$n * k, 0, D$maxX), k, D$n)
  W0
}
palm_init_H <- function(W, D, hp, nsweep = 20L) {
  H <- matrix(0, ncol(W), D$n)                               # I11
  for (i in seq_len(nsweep)) H <- hals_sweep(W, H, D, hp)
  H
}
# eps: 1e-6 x median of sigma_j (eps = 0) over j and over all starts' projected initial W; median <= 0 fails
palm_eps <- function(Ws, D) {
  m <- median(unlist(lapply(Ws, function(W) palm_scores(W, D, 0)$sigma)))
  if (!is.finite(m) || m <= 0) stop("palm_eps: pooled median initial score sd is <= 0 or nonfinite", call. = FALSE)
  1e-6 * m
}
# ---- one start ----------------------------------------------------------------------------------
# Status: converged (r < tol), nonconverged (work cap reached, or a line-search failure under "ncg"; finite, no
# failure) or failure (nonfinite value, W step floor ("pg"), FISTA cap, or any F increase
# F_new > F_old + 1e-13 max(1, |F_old|) after a block). The start ends at the first of these. trace = TRUE records
# F after each block, and under "ncg" the amendment-3 section 5 per-sweep diagnostics, for an independent recount.
# w_solver = "ncg" (amendment 3): Newton-CG W step (w_step_ncg); maxit is the amendment-3 per-start sweep budget
# (default 1,000 sweeps), NOT the old PG iteration count (PI ruling B3.3), and hvp_cap the per-start W-block HVP
# total (25,000). Reaching either cap -> nonconverged, reason "work_cap_reached"; the HVP cap is applied after the
# sweep in which it is reached completes and the certificate is checked (PI ruling B3.2). Line-search failure ->
# nonconverged, reason "line_search_failed", W kept. w_solver = "pg": the frozen projected-gradient step (w_step,
# PALM_DEFS$step_init), maxit default 5,000 iterations, unchanged.
# max_secs is a test-harness wall-clock budget (default Inf, not part of the spec); exceeding it ends the start
# with status "incomplete", which tests report separately.
palm_fit_one <- function(D, hp, W, H, b, maxit = NULL, tol = 1e-5, trace = FALSE, defs = PALM_DEFS, max_secs = Inf,
                         w_solver = c("ncg", "pg"), hvp_cap = 25000L) {
  w_solver <- match.arg(w_solver); ncg <- w_solver == "ncg"
  if (is.null(maxit)) maxit <- if (ncg) 1000L else 5000L        # ncg: sweep budget; pg: frozen maxit
  inc_tol <- function(Fo) 1e-13 * max(1, abs(Fo))
  ninc <- c(H = 0L, W = 0L, b = 0L); maxinc <- c(H = 0, W = 0, b = 0)
  t_prev <- 1; Lb <- 1; nbt <- 0L; bits <- 0L; nrej <- 0L; norc <- 0L; nnull <- 0L
  hvp_tot <- 0L; wex <- c(w_converged = 0L, newton = 0L, line_search_failed = 0L, numerical_failure = 0L)
  cgx <- c(forcing = 0L, hvp_ceiling = 0L, hvp_total_cap = 0L, curvature = 0L); nneg <- 0L; nnc1 <- 0L
  status <- NA_character_; reason <- ""
  fail <- function(why) { status <<- "failure"; reason <<- why }
  Fcur <- palm_F(W, H, b, D, hp); F0 <- Fcur
  r <- c(r = NA_real_, rW = NA_real_, rH = NA_real_, rb = NA_real_, r_original = NA_real_)
  if (!is.finite(Fcur)) fail("nonfinite F at start")
  tr <- if (trace) vector("list", maxit) else NULL
  inc <- function(Fn, Fo, blk) {
    di <- Fn - Fo
    if (is.finite(di) && di > maxinc[blk]) maxinc[blk] <<- di
    if (is.finite(di) && di > inc_tol(Fo)) { ninc[blk] <<- ninc[blk] + 1L; fail(paste("F increase in", blk, "block")) }
  }
  it <- 0L; t0 <- proc.time()[["elapsed"]]
  while (is.na(status) && it < maxit) {
    if (proc.time()[["elapsed"]] - t0 > max_secs) { status <- "incomplete"; reason <- "test time budget"; break }
    it <- it + 1L; Fit <- Fcur; FW <- NA_real_; Fb <- NA_real_; ws <- list(nbt = 0L, fail = FALSE); sbit <- 0L; rej <- FALSE
    wd <- if (ncg) c(rW_start = NA, w_exit = NA, face = NA, cg_iter = NA, hvp = NA, cg_exit = NA, negative_curvature = NA,
                     nc_first_iter = NA, n_annih = NA, n_enter = NA, n_leave = NA, slope = NA, res_rec = NA, forcing_thr = NA,
                     hvp_total = hvp_tot) else NULL                 # trace columns when the W block does not run
    H <- hals_sweep(W, H, D, hp); FH <- palm_F(W, H, b, D, hp)
    if (!is.finite(FH)) { fail("nonfinite F after H block"); break }
    inc(FH, Fcur, "H"); Fcur <- FH
    if (is.na(status) && ncg) {
      ws <- w_step_ncg(W, H, b, D, hp, hvp_left = hvp_cap - hvp_tot, tol = tol, defs = defs)
      hvp_tot <- hvp_tot + ws$hvp; nbt <- nbt + ws$nbt; wex[ws$exit] <- wex[ws$exit] + 1L
      if (ws$cg_exit %in% names(cgx)) cgx[ws$cg_exit] <- cgx[ws$cg_exit] + 1L
      nneg <- nneg + ws$negative_curvature; nnc1 <- nnc1 + ws$nc_first_iter
      slope <- NA_real_
      if (ws$exit == "numerical_failure") fail("nonfinite in W block (Newton-CG)") else
      if (ws$exit == "line_search_failed") { status <- "nonconverged"; reason <- "line_search_failed" } else {
        if (ws$exit == "newton") { slope <- sum(ws$GW * (ws$W - W)); t_prev <- ws$t }   # recorded for the Armijo recount
        W <- ws$W; FW <- palm_F(W, H, b, D, hp)
        if (!is.finite(FW)) fail("nonfinite F after W block") else { inc(FW, Fcur, "W"); Fcur <- FW }
      }
      wd <- c(rW_start = ws$rW, w_exit = match(ws$exit, names(wex)) - 1L, face = ws$face, cg_iter = ws$cg_iter, hvp = ws$hvp,
              cg_exit = match(ws$cg_exit, names(cgx), nomatch = 0L), negative_curvature = ws$negative_curvature,
              nc_first_iter = ws$nc_first_iter, n_annih = ws$n_annih, n_enter = ws$n_enter, n_leave = ws$n_leave,
              slope = slope, res_rec = ws$res_rec, forcing_thr = ws$forcing_thr, hvp_total = hvp_tot)
    }
    if (is.na(status) && !ncg) {
      ws <- w_step(W, H, b, D, hp, t_prev, defs = defs); nbt <- nbt + ws$nbt
      if (ws$fail) fail("W step floor") else {
        if (identical(ws$W, W)) nnull <- nnull + 1L                      # accepted null step (recorded only)
        W <- ws$W; t_prev <- ws$t; FW <- palm_F(W, H, b, D, hp)
        if (!is.finite(FW)) fail("nonfinite F after W block") else { inc(FW, Fcur, "W"); Fcur <- FW }
      }
    }
    if (is.na(status)) {
      sb <- solve_b(b, W, D, hp, L0 = max(Lb / 4, 1e-8)); bits <- bits + sb$it; sbit <- sb$it; norc <- norc + sb$n_oracle
      if (!sb$conv) fail("FISTA cap or nonfinite in b block") else {
        Lb <- sb$L; Fn <- palm_F(W, H, sb$b, D, hp)
        if (isTRUE(defs$b_accept(Fn, Fcur))) b <- sb$b else { nrej <- nrej + 1L; rej <- TRUE }
        Fb <- palm_F(W, H, b, D, hp)
        if (!is.finite(Fb)) fail("nonfinite F after b block") else { inc(Fb, Fcur, "b"); Fcur <- Fb }
      }
    }
    if (is.na(status)) {   # finiteness checked on binding residuals only; r_original binds nothing
      r <- palm_resid(W, H, b, D, hp, defs)
      if (any(!is.finite(r[c("r", "rW", "rH", "rb")]))) fail("nonfinite residual") else if (r[["r"]] < tol) status <- "converged"
    }
    if (is.na(status) && ncg && hvp_tot >= hvp_cap) { status <- "nonconverged"; reason <- "work_cap_reached" }
    if (trace) tr[[it]] <- c(iter = it, Fstart = Fit, FH = FH, FW = FW, Fb = Fb, r, t = if (ncg) { if (is.null(ws$t)) NA else ws$t } else t_prev, nbt = ws$nbt,
                             bit = sbit, b_rejected = rej, failed = identical(status, "failure"), wd)
  }
  if (is.na(status)) { status <- "nonconverged"; if (ncg) reason <- "work_cap_reached" }
  secs <- proc.time()[["elapsed"]] - t0
  pr <- palm_parts(W, H, b, D, hp)
  list(W = W, H = H, b = b, beta = pr$beta, sigma = pr$sigma, F = Fcur, F0 = F0, iter = it, status = status,
       reason = reason, conv = identical(status, "converged"), r = r, ninc = ninc, maxinc = maxinc,
       b_rejected = nrej, b_iters = bits, n_oracle = norc, w_backtracks = nbt, w_null_steps = nnull, t_last = t_prev, secs = secs,
       eps = hp$eps, w_solver = w_solver, hvp_total = hvp_tot, w_exits = wex, cg_exits = cgx, n_negative_curvature = nneg,
       n_nc_first_iter = nnc1, trace = if (trace) do.call(rbind, tr[seq_len(it)]) else NULL)
}

# ---- fit: 30 starts, status rule ---------------------------------------------------------------
# Returns the lowest-F converged start; else the lowest-F nonconverged start (fit status "nonconverged");
# else fit status "failed" with best = NULL. hp: alpha, lambda, nu, rho (lambdaH defaults to 1e-6).
# w_solver (default "ncg") and hvp_cap pass through to palm_fit_one; maxit = NULL takes palm_fit_one's default for
# the solver (ncg: the amendment-3 sweep budget of 1,000 sweeps, not the old PG iteration count; pg: 5,000).
pick_start <- function(Fs, status) {
  for (st in c("converged", "nonconverged")) {
    ok <- which(status == st)
    if (length(ok)) return(list(i = ok[which.min(Fs[ok])], status = st))
  }
  list(i = NA_integer_, status = "failed")
}
palm_fit <- function(D, hp, k = 3L, seed = 123L, nstarts = 30L, maxit = NULL, tol = 1e-5,
                     ncores = 1L, keep_starts = FALSE, defs = PALM_DEFS, trace = FALSE, max_secs = Inf,
                     w_solver = c("ncg", "pg"), hvp_cap = 25000L) {
  w_solver <- match.arg(w_solver)
  if (is.null(hp$lambdaH)) hp$lambdaH <- 1e-6
  Ws <- lapply(seq_len(nstarts), function(s) proj_S(palm_draw(D, k, seed, s)))
  hp$eps <- palm_eps(Ws, D)
  one <- function(s) {
    W <- Ws[[s]]; H <- palm_init_H(W, D, hp)
    f <- tryCatch(palm_fit_one(D, hp, W, H, rep(0, k), maxit = maxit, tol = tol, defs = defs, trace = trace, max_secs = max_secs,
                                       w_solver = w_solver, hvp_cap = hvp_cap),
                  error = function(e) list(status = "failure", reason = paste("error:", conditionMessage(e)), F = NA_real_,
                                           iter = NA_integer_, r = c(r = NA, rW = NA, rH = NA, rb = NA, r_original = NA),
                                           ninc = c(H = NA, W = NA, b = NA), b_rejected = NA, w_backtracks = NA, b_iters = NA,
                                           secs = NA, eps = hp$eps))
    f$start <- s; f
  }
  fits <- if (ncores > 1L) parallel::mclapply(seq_len(nstarts), one, mc.cores = ncores) else lapply(seq_len(nstarts), one)
  if (!all(vapply(fits, is.list, TRUE))) stop("palm_fit: a worker returned no result", call. = FALSE)
  g <- function(nm, typ) vapply(fits, function(f) { v <- f[[nm]]; if (is.null(v)) typ[NA] else v }, typ)
  rr <- function(nm) vapply(fits, function(f) unname(f$r[nm]), 0)
  wx <- function(nm) vapply(fits, function(f) if (is.null(f$w_exits)) NA_real_ else unname(f$w_exits[[nm]]) + 0, 0)
  summ <- data.frame(start = seq_len(nstarts), status = g("status", ""), reason = g("reason", ""), F = g("F", 0),
                     iter = g("iter", 0L), r = rr("r"), rW = rr("rW"), rH = rr("rH"), rb = rr("rb"),
                     r_original = rr("r_original"), ninc = vapply(fits, function(f) sum(f$ninc), 0),
                     b_rejected = g("b_rejected", 0), w_backtracks = g("w_backtracks", 0), w_null_steps = g("w_null_steps", 0), b_iters = g("b_iters", 0),
                     hvp_total = g("hvp_total", 0L), w_newton = wx("newton"), w_noop = wx("w_converged"),
                     w_ls_failed = wx("line_search_failed"), eps = g("eps", 0), secs = g("secs", 0), stringsAsFactors = FALSE)
  pk <- pick_start(summ$F, summ$status)
  list(status = pk$status, best = if (is.na(pk$i)) NULL else fits[[pk$i]], best_start = pk$i, starts = summ,
       eps = hp$eps, hp = hp, seed = seed, w_solver = w_solver, fits = if (keep_starts) fits else NULL)
}

# ---- stored model and scoring (spec section 4, prediction vector) -------------------------------
# v = W diag(1/sigma) b and v~ = v - mean(v) over the training gene universe; G = unique union of
# desurv_get_top_genes(W, ntop)$top_genes. Centering happens before the support restriction; v~ is
# never recentered on another dataset's genes.
palm_model <- function(fit, D, ntop, top_genes_fun = DeSurv::desurv_get_top_genes) {
  if (is.null(fit$best)) stop("palm_model: fit has no returned start (status failed)", call. = FALSE)
  W <- fit$best$W; rownames(W) <- D$genes
  v <- as.vector(W %*% (fit$best$b / fit$best$sigma)); names(v) <- D$genes
  vt <- v - mean(v)
  G <- unique(unlist(top_genes_fun(W, ntop)$top_genes))
  list(genes = D$genes, v_tilde = vt, G = G, ntop = ntop, W = W, b = fit$best$b, beta = fit$best$beta,
       sigma = fit$best$sigma, eps = fit$eps, hp = fit$hp, fit_status = fit$status)
}
palm_score <- function(Xnew, model, defs = PALM_DEFS) defs$score(Xnew, model)

# ---- fold-local preprocessing ------------------------------------------------------------------
# Within-sample ranks (average ties) of the training columns only, then palm_data. Gene filtering is not
# done here (the gate's preprocessing step applies the published ngene filter on training data).
palm_rank <- function(Xraw) { R <- apply(Xraw, 2, rank, ties.method = "average"); dimnames(R) <- dimnames(Xraw); R }
palm_fold_data <- function(Xraw, y, d, train) palm_data(palm_rank(Xraw[, train, drop = FALSE]), y[train], d[train])
