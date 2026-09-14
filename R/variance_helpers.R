# R/variance_helpers.R — Per-factor cross-sample variance vs survival
# contribution helpers, used by code/09a_figures.R (main text Fig 2C) and
# code/09b_si_figures.R (SI k=5 variance/survival panel).
#
# Cross-sample variance fraction: ||(W_j H_j^T)^(c)||_F^2 / ||X^(c)||_F^2
# where (c) denotes gene-wise centering. The uncentered analogue
# ||W_j H_j^T||^2 / ||X||^2 is dominated by the baseline of rank-normalized
# X and under-represents factors that encode genuine cross-sample variation
# (GitHub issue #6: the exocrine factor N2 appeared at 6.9% under the
# uncentered metric, contradicting the paper's own thesis; the corrected
# centered metric puts N2 at 39.3% in PDAC NMF at k=3). For non-orthogonal
# factorizations such as NMF the per-factor centered shares do not partition
# to 100% (in PDAC NMF k=3 they sum to approximately 1.97); see Results
# paragraph footnote and SI Section 10 ("Variance and Survival Contribution
# per Factor").

build_var_surv_df <- function(W, H, X, time, event, method) {
  W <- as.matrix(W); H <- as.matrix(H); X <- as.matrix(X)
  k <- ncol(W)

  # Cross-sample variance fraction (gene-centered)
  Xc <- X - rowMeans(X, na.rm = TRUE)
  var_den <- sum(Xc^2, na.rm = TRUE)
  var_exp <- vapply(seq_len(k), function(j) {
    comp  <- W[, j] %o% H[j, ]
    compc <- comp - rowMeans(comp, na.rm = TRUE)
    sum(compc^2, na.rm = TRUE) / var_den
  }, numeric(1))

  # Per-factor survival contribution: Δ partial log-likelihood when factor j
  # is dropped from the full k-factor Cox model. Uses XtW = t(X) %*% W per
  # the original convention (R/figure_targets.R, since deleted).
  XtW     <- t(X) %*% W
  ll_full <- survival::coxph(survival::Surv(time, event) ~ XtW)$loglik[2]
  delta_ll <- vapply(seq_len(k), function(j) {
    XtW_mj <- XtW[, -j, drop = FALSE]
    reduced <- if (ncol(XtW_mj) == 0L) {
      survival::coxph(survival::Surv(time, event) ~ 1)
    } else {
      survival::coxph(survival::Surv(time, event) ~ XtW_mj)
    }
    ll_full - reduced$loglik[2]
  }, numeric(1))

  data.frame(method = method, factor = seq_len(k),
             variance_explained = var_exp, delta_loglik = delta_ll,
             stringsAsFactors = FALSE)
}
