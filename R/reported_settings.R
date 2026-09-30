# R/reported_settings.R — one source for the reported DeSurv estimator settings
#
# Rerun plan rev 3 (docs/plans/2026-09-30-rerun-plan-rev3.md), section 2: the
# BO objective in each CV fold is the REPORTED estimator, the code/04 recipe
# (seed fits -> desurv_consensus_seed -> one final desurv_fit), packaged as
# DeSurv::desurv_reported_fit(). The BO's CV folds (code/03, via
# DeSurv::desurv_cv(estimator = "reported")) and the final fit
# (fit_reported_desurv() below) therefore use the SAME entry point with the SAME
# seed-fit settings, taken from this file only. The values are the code/04
# settings (code/04_fit_models.R:30-38, 78-79, 112: NINIT_FULL seed fits,
# tol = RUN_TOL / 100, tol_init = RUN_TOL, maxit = imaxit = RUN_MAXIT).
#
# Pure definitions: sourcing this file reads and writes nothing.

DESURV_ESTIMATOR_ENTRY <- "desurv_reported_fit"

# The function DeSurv::desurv_cv(estimator = "reported") calls in every fold.
desurv_estimator_entry <- function() {
  DeSurv::desurv_reported_fit
}

desurv_reported_settings <- function(quick = FALSE) {
  run_tol   <- if (isTRUE(quick)) 1e-4 else 1e-5
  run_maxit <- if (isTRUE(quick)) 500L else 4000L
  list(
    n_seed_fits = if (isTRUE(quick)) 3L else 100L,   # code/04 NINIT_FULL
    tol         = run_tol / 100,
    tol_init    = run_tol,
    maxit       = run_maxit,
    imaxit      = run_maxit
  )
}

# Final-fit entry point: the reported estimator at a configuration from
# tar_params_best, with the shared settings. (code/04 is not yet rewired to call
# this; that is part of the full rerun, plan section 5.)
fit_reported_desurv <- function(X, y, d, params, seeds, settings,
                                ncores = 1L) {
  get_num <- function(nm, default) {
    v <- params[[nm]]
    if (is.null(v) || length(v) != 1L || is.na(v)) default else as.numeric(v)
  }
  desurv_estimator_entry()(
    X = X, y = y, d = d,
    k = as.integer(params$k), alpha = get_num("alpha", NA_real_),
    lambda = get_num("lambda", NA_real_), nu = get_num("nu", NA_real_),
    lambdaW = get_num("lambdaW", 0), lambdaH = get_num("lambdaH", 0),
    ntop = as.integer(round(get_num("ntop", 100))),   # code/04 ntop_value rule
    seeds = seeds,
    tol = settings$tol, tol_init = settings$tol_init,
    maxit = settings$maxit, imaxit = settings$imaxit,
    ncores = ncores
  )
}
