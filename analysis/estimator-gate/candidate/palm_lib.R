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
    ll <- sum(ds * (es - m - log(S)))
    h <- cumsum(ds / S)[D$ge]
    g <- numeric(length(eta)); g[o] <- ds - e * h
    return(list(ll = ll, g = g))
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
# Status: converged (r < tol), nonconverged (maxit, finite, no failure) or failure (nonfinite value, W step
# floor, FISTA cap, or any F increase F_new > F_old + 1e-13 max(1, |F_old|) after a block). The start ends at
# the first of these. trace = TRUE records F after each block for an independent recount.
# max_secs is a test-harness wall-clock budget (default Inf, not part of the spec); exceeding it ends the start
# with status "incomplete", which tests report separately.
palm_fit_one <- function(D, hp, W, H, b, maxit = 5000L, tol = 1e-5, trace = FALSE, defs = PALM_DEFS, max_secs = Inf) {
  inc_tol <- function(Fo) 1e-13 * max(1, abs(Fo))
  ninc <- c(H = 0L, W = 0L, b = 0L); maxinc <- c(H = 0, W = 0, b = 0)
  t_prev <- 1; Lb <- 1; nbt <- 0L; bits <- 0L; nrej <- 0L; norc <- 0L; nnull <- 0L
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
    H <- hals_sweep(W, H, D, hp); FH <- palm_F(W, H, b, D, hp)
    if (!is.finite(FH)) { fail("nonfinite F after H block"); break }
    inc(FH, Fcur, "H"); Fcur <- FH
    if (is.na(status)) {
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
    if (is.na(status)) {
      r <- palm_resid(W, H, b, D, hp, defs)
      if (any(!is.finite(r))) fail("nonfinite residual") else if (r[["r"]] < tol) status <- "converged"
    }
    if (trace) tr[[it]] <- c(iter = it, Fstart = Fit, FH = FH, FW = FW, Fb = Fb, r, t = t_prev, nbt = ws$nbt,
                             bit = sbit, b_rejected = rej, failed = identical(status, "failure"))
  }
  if (is.na(status)) status <- "nonconverged"
  secs <- proc.time()[["elapsed"]] - t0
  pr <- palm_parts(W, H, b, D, hp)
  list(W = W, H = H, b = b, beta = pr$beta, sigma = pr$sigma, F = Fcur, F0 = F0, iter = it, status = status,
       reason = reason, conv = identical(status, "converged"), r = r, ninc = ninc, maxinc = maxinc,
       b_rejected = nrej, b_iters = bits, n_oracle = norc, w_backtracks = nbt, w_null_steps = nnull, t_last = t_prev, secs = secs,
       eps = hp$eps, trace = if (trace) do.call(rbind, tr[seq_len(it)]) else NULL)
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
                     ncores = 1L, keep_starts = FALSE, defs = PALM_DEFS, trace = FALSE, max_secs = Inf) {
  if (is.null(hp$lambdaH)) hp$lambdaH <- 1e-6
  Ws <- lapply(seq_len(nstarts), function(s) proj_S(palm_draw(D, k, seed, s)))
  hp$eps <- palm_eps(Ws, D)
  one <- function(s) {
    W <- Ws[[s]]; H <- palm_init_H(W, D, hp)
    f <- tryCatch(palm_fit_one(D, hp, W, H, rep(0, k), maxit = maxit, tol = tol, defs = defs, trace = trace, max_secs = max_secs),
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
  summ <- data.frame(start = seq_len(nstarts), status = g("status", ""), reason = g("reason", ""), F = g("F", 0),
                     iter = g("iter", 0L), r = rr("r"), rW = rr("rW"), rH = rr("rH"), rb = rr("rb"),
                     r_original = rr("r_original"), ninc = vapply(fits, function(f) sum(f$ninc), 0),
                     b_rejected = g("b_rejected", 0), w_backtracks = g("w_backtracks", 0), w_null_steps = g("w_null_steps", 0), b_iters = g("b_iters", 0),
                     eps = g("eps", 0), secs = g("secs", 0), stringsAsFactors = FALSE)
  pk <- pick_start(summ$F, summ$status)
  list(status = pk$status, best = if (is.na(pk$i)) NULL else fits[[pk$i]], best_start = pk$i, starts = summ,
       eps = hp$eps, hp = hp, seed = seed, fits = if (keep_starts) fits else NULL)
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
