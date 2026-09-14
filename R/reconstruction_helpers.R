# R/reconstruction_helpers.R — reconstruction-based factor metrics (GitHub issue #7)
#
# Reframes per-factor contribution as RECONSTRUCTION of the expression matrix
# rather than "variance explained". NMF seeks a low-rank reconstruction of X;
# many factorizations reconstruct comparably well and the unsupervised
# objective gives no incentive to favor a prognostic direction.
#
# Provides:
#   recon_shapley(W, H, X)        — per-factor reconstruction Shapley share
#                                   (sums to 100% by the efficiency axiom;
#                                    replaces the non-partitioning centered-
#                                    variance metric from issue #6).
#   build_recon_surv_df(W,H,X,...) — Fig 2C data: x = Shapley share,
#                                    y = Delta partial log-likelihood.
#   projection_r2(src_W, basis_W) — OLS R^2 of each src factor loading vector
#                                   onto the other method's basis.
#   overall_recon_r2(W, H, X)     — whole-matrix reconstruction R^2.

# ── §1 Reconstruction Shapley value ─────────────────────────────────────────
# v(S) = ||X||^2 - ||X - sum_{k in S} w_k h_k^T||^2 ; phi_k = avg marginal gain
# across all orderings. Efficiency: sum_k phi_k = v(N), so normalized shares
# sum to exactly 1.
recon_shapley <- function(W, H, X) {
  W <- as.matrix(W); H <- as.matrix(H); X <- as.matrix(X)
  Xnorm2 <- sum(X^2)
  K  <- ncol(W)
  Rk <- lapply(seq_len(K), function(j) W[, j] %o% H[j, ])
  v  <- function(S) if (!length(S)) 0 else Xnorm2 - sum((X - Reduce(`+`, Rk[S]))^2)
  subs <- do.call(c, lapply(0:K, function(m) utils::combn(K, m, simplify = FALSE)))
  phi <- numeric(K)
  for (k in seq_len(K)) for (S in Filter(function(z) !(k %in% z), subs)) {
    w <- factorial(length(S)) * factorial(K - length(S) - 1) / factorial(K)
    phi[k] <- phi[k] + w * (v(c(S, k)) - v(S))
  }
  list(phi = phi, share = phi / sum(phi), vN = v(seq_len(K)))
}

# ── Fig 2C data: Shapley x-axis + Type III partial-likelihood y-axis ────────
build_recon_surv_df <- function(W, H, X, time, event, method) {
  W <- as.matrix(W); H <- as.matrix(H); X <- as.matrix(X)
  k <- ncol(W)
  share <- recon_shapley(W, H, X)$share

  XtW     <- t(X) %*% W
  ll_full <- survival::coxph(survival::Surv(time, event) ~ XtW)$loglik[2]
  delta_ll <- vapply(seq_len(k), function(j) {
    XtW_mj <- XtW[, -j, drop = FALSE]
    reduced <- if (ncol(XtW_mj) == 0L) survival::coxph(survival::Surv(time, event) ~ 1)
               else survival::coxph(survival::Surv(time, event) ~ XtW_mj)
    ll_full - reduced$loglik[2]
  }, numeric(1))

  # NOTE: the `variance_explained` column holds the reconstruction Shapley
  # share, not a variance fraction. The name is retained for API
  # compatibility with code/09a_figures.R, code/09b_si_figures.R, and the
  # paper Rmd setup chunk, which read $variance_explained as the Fig 2C
  # x-axis. (Issue #7 changed the metric, not the column name.)
  data.frame(method = method, factor = seq_len(k),
             variance_explained = share, delta_loglik = delta_ll,
             stringsAsFactors = FALSE)
}

# ── §5 Cross-method projection R^2 ──────────────────────────────────────────
# How reconstructable is each src factor's gene-loading vector from the other
# method's basis (OLS R^2). Low R^2 => the factor is novel to src.
projection_r2 <- function(src_W, basis_W) {
  src_W <- as.matrix(src_W); basis_W <- as.matrix(basis_W)
  vapply(seq_len(ncol(src_W)),
         function(j) summary(stats::lm(src_W[, j] ~ basis_W))$r.squared,
         numeric(1))
}

# ── §6 Overall whole-matrix reconstruction R^2 ──────────────────────────────
# NOTE ON THE DENOMINATOR. This uses the UNCENTERED total sum of squares,
# sum(X^2). X is the within-sample rank matrix, whose entries run 1..p with
# grand mean p/2, so the grand mean alone accounts for about 75% of sum(X^2)
# and any nonnegative low-rank fit scores high by construction. The uncentered
# value is not comparable to an ordinary regression R^2 and must not be read as
# "fraction of expression variation explained". Report it alongside the
# centered version below rather than on its own.
overall_recon_r2 <- function(W, H, X) {
  W <- as.matrix(W); H <- as.matrix(H); X <- as.matrix(X)
  1 - sum((X - W %*% H)^2) / sum(X^2)
}

# Conventional mean-centered R^2: residual sum of squares against the total sum
# of squares about the grand mean. This is what a reader assumes on seeing
# "R^2", and it is the definition to quote when stating how much reconstruction
# fidelity is given up.
overall_recon_r2_centered <- function(W, H, X) {
  W <- as.matrix(W); H <- as.matrix(H); X <- as.matrix(X)
  1 - sum((X - W %*% H)^2) / sum((X - mean(X))^2)
}
