# R/bo_setup.R: shared BO configuration for code/03 (all stages), code/03b (rank-stratified
# BO, rev 3 amendment 1) and code/03c (its cross-rank selection).
#
# Moved verbatim from code/03_bayesian_optimization.R (lines 28-139 at 8a587d0) so that every
# BO in the rerun uses one set of values; the constructed BO_BOUNDS, BO_FIXED and BO_COMMON
# were checked identical() before and after the move (docs/plans/2026-10-02-amendment1-impl-report.md).
#
# Requires, already loaded by the caller: CONFIG (code/00_helpers.R) and
# desurv_reported_settings() (R/reported_settings.R). run_bo_finalists() reads a global
# `tar_data` when called. Defines objects only; reads and writes no files.

# ── Configuration ─────────────────────────────────────────────────────────
# Rerun plan rev 3 (docs/plans/2026-09-30-rerun-plan-rev3.md), section 2:
# every BO evaluation scores the reported estimator (the code/04 recipe: seed
# fits -> consensus seed -> final fit, DeSurv::desurv_reported_fit) in each
# training fold, with code/04's seed-fit count and tolerances, on 3 repeated CV
# partitions with the same partitions and per-fold seed banks for every
# configuration; the incumbent is the GP posterior-mean argmax; k is chosen
# among re-evaluated finalists by the 1-SE rule on the empirical paired SE.
BO_ESTIMATOR    <- "reported"           # "mean_starts" = pre-rev-2 evaluator
BO_INCUMBENT    <- "posterior_mean"     # "observed" = pre-rev-2 argmax
REPORTED        <- desurv_reported_settings(CONFIG$quick)    # = final-fit settings
BO_N_SEED_FITS  <- REPORTED$n_seed_fits # code/04 NINIT_FULL: 100 full, 3 quick
CV_SEED_BANK_BASE <- 4000000L           # plan section 9: per-fold bank 4,000,001 onward
if (CONFIG$quick) {
  NINIT          <- 2L                  # used only if BO_ESTIMATOR == "mean_starts"
  BO_N_INIT      <- 4L
  BO_N_ITER      <- 4L
  NCORES_GRID    <- 1L
  PARALLEL       <- FALSE
  STD_NMF_NRUN   <- 2L
  STD_NMF_K_GRID <- 2:6
  BO_K_UPPER     <- 6L
  BO_NGENE       <- 500L
  BO_CV_SEEDS    <- 20261101L           # one partition in quick mode
  FINALISTS_PER_K <- 2L
} else {
  NINIT          <- 30L                 # used only if BO_ESTIMATOR == "mean_starts"
  BO_N_INIT      <- 50L
  BO_N_ITER      <- 100L
  NCORES_GRID    <- CONFIG$ncores
  PARALLEL       <- CONFIG$ncores > 1
  STD_NMF_NRUN   <- 30L
  STD_NMF_K_GRID <- 2:12
  BO_K_UPPER     <- 12L
  BO_NGENE       <- 3000L
  BO_CV_SEEDS    <- c(20261101L, 20261102L, 20261103L)   # plan section 9
  FINALISTS_PER_K <- 3L
}
# Seed-fit and final-fit tolerances inside CV = the final fit's (one source).
BO_TOL   <- REPORTED$tol
BO_MAXIT <- REPORTED$maxit

BO_BOUNDS <- list(
  k_grid     = list(lower = 2L, upper = BO_K_UPPER, type = "integer"),
  alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
  nu_grid    = list(lower = 0, upper = 1, type = "continuous")
)

BO_FIXED <- list(
  ngene        = BO_NGENE,
  lambdaW_grid = 0,
  lambdaH_grid = 0
)
if (identical(BO_ESTIMATOR, "mean_starts")) BO_FIXED$n_starts <- NINIT

# ── ntop handling ──────────────────────────────────────────────────────────
BO_BOUNDS$ntop <- list(lower = CONFIG$ntop_lower, upper = CONFIG$ntop_upper,
                       type = "integer")
message("  ntop BO bounds: [", CONFIG$ntop_lower, ", ", CONFIG$ntop_upper, "]")

BO_COMMON <- list(
  preprocess        = TRUE,
  method_trans_train = "rank",
  engine            = "warmstart",
  nfolds            = 5L,
  tol               = BO_TOL,
  maxit             = BO_MAXIT,
  max_refinements   = 0L,
  tol_gain          = 0.002,
  plateau           = 1L,
  top_k             = 10L,
  shrink_base       = 0.3,
  importance_gain   = 0.1,
  coarse_control    = list(n_init = BO_N_INIT, n_iter = BO_N_ITER, candidate_pool = 4000L,
                           exploration_weight = 0.01, seed = 123),
  refine_control    = NULL,
  verbose           = TRUE,
  parallel_grid     = PARALLEL,
  ncores_grid       = NCORES_GRID,
  estimator         = BO_ESTIMATOR,
  incumbent         = BO_INCUMBENT,
  cv_seeds          = BO_CV_SEEDS,
  n_seed_fits       = BO_N_SEED_FITS,
  seed_bank_base    = CV_SEED_BANK_BASE,
  reported_control  = list(tol_init = REPORTED$tol_init, imaxit = REPORTED$imaxit)
)

# Finalists: top FINALISTS_PER_K configurations per k by GP posterior mean,
# re-evaluated on the identical repeated CV (same partitions and seed banks).
run_bo_finalists <- function(bo_results) {
  DeSurv::desurv_bo_finalists(
    bo_results,
    X = tar_data$ex,
    y = tar_data$sampInfo$time,
    d = tar_data$sampInfo$event,
    dataset = tar_data$sampInfo$dataset,
    samp_keeps = tar_data$samp_keeps,
    n_per_k = FINALISTS_PER_K,
    cv_overrides = list(parallel_grid = PARALLEL, ncores_grid = NCORES_GRID)
  )
}

# ngene is fixed in BO_FIXED (not tuned), so the filtered data use that value.
bo_ngene <- function(params) {
  if (!is.null(params$ngene) && !is.na(params$ngene)) {
    as.integer(round(params$ngene))
  } else {
    as.integer(BO_FIXED$ngene)
  }
}

