# R/consensus_settings.R — one source for the consensus-DeSurv estimator settings
#
# Robust-rerun plan rev 2.1 (docs/plans/2026-09-30-robust-rerun-plan.md), section 2:
# "The BO tunes the estimator that is reported, not a cheaper approximation of it."
# The BO's CV folds (code/03, via DeSurv::desurv_cv(estimator = "consensus")) and
# the final fit (fit_consensus_desurv() below) therefore use the SAME entry point,
# DeSurv::desurv_consensus_fit(), with the SAME constituent-run settings, taken
# from this file only. The values are the code/04 seed-fit settings
# (code/04_fit_models.R:30-38, 78-79: tol = RUN_TOL / 100, tol_init = RUN_TOL,
# maxit = imaxit = RUN_MAXIT).
#
# The section 3 rule parameters (outlier and competing-mode rules) are frozen at
# desurv_consensus_fit()'s defaults and are deliberately not set anywhere here.
#
# Pure definitions: sourcing this file reads and writes nothing.

DESURV_ESTIMATOR_ENTRY <- "desurv_consensus_fit"

# The namespace binding that DeSurv::desurv_cv() calls inside every CV fold.
# Looked up in the namespace (not with `DeSurv::`) because, until Track 1 lands,
# the package carries an unexported TRACK1-STUB of this function; once Track 1's
# exported version is merged this is the same object as DeSurv::desurv_consensus_fit.
desurv_estimator_entry <- function() {
  get(DESURV_ESTIMATOR_ENTRY, envir = asNamespace("DeSurv"), inherits = FALSE)
}

desurv_consensus_settings <- function(quick = FALSE) {
  run_tol   <- if (isTRUE(quick)) 1e-4 else 1e-5
  run_maxit <- if (isTRUE(quick)) 500L else 4000L
  list(
    tol      = run_tol / 100,
    tol_init = run_tol,
    maxit    = run_maxit,
    imaxit   = run_maxit
  )
}

# Final-fit entry point: the consensus estimator at a configuration from
# tar_params_best, with the shared settings. (code/04 is not yet rewired to call
# this; that is part of the full rerun, plan section 5.)
fit_consensus_desurv <- function(X, y, d, params, seeds, settings,
                                 ncores = 1L, verbose = FALSE) {
  get_num <- function(nm, default) {
    v <- params[[nm]]
    if (is.null(v) || length(v) != 1L || is.na(v)) default else as.numeric(v)
  }
  desurv_estimator_entry()(
    X = X, y = y, d = d,
    k = as.integer(params$k), alpha = get_num("alpha", NA_real_),
    lambda = get_num("lambda", NA_real_), nu = get_num("nu", NA_real_),
    lambdaW = get_num("lambdaW", 0), lambdaH = get_num("lambdaH", 0),
    seeds = seeds,
    tol = settings$tol, tol_init = settings$tol_init,
    maxit = settings$maxit, imaxit = settings$imaxit,
    ncores = ncores, verbose = verbose
  )
}
