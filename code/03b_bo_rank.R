#!/usr/bin/env Rscript
# code/03b_bo_rank.R — rank-stratified BO, one rank per run (rev 3 amendment 1)
#
# docs/plans/2026-10-02-rev3-amendment-1-rank-stratified-bo.md (APPROVED). Only the search
# allocation changes: k is fixed to DESURV_BO_K (2..12) and alpha, lambda, nu, ntop are
# searched with the main BO's estimator, partitions, seed banks, bounds and acquisition
# (R/bo_setup.R, shared with code/03). The BO is initialized from every valid main-BO
# evaluation at this k; k = 3 also receives the published and May-best configurations,
# evaluated earlier with the identical recipe, partitions and seed bank 4,000,000
# (docs/plans/2026-10-02-step4-bo-diagnosis.md, job 3500640).
#
# Env:     DESURV_BO_K (required), plus the usual DESURV_RESULTS_DIR / DESURV_QUICK / DESURV_NCORES.
# Inputs:  tar_data_tcgacptac, desurv_bo_results_tcgacptac (the main BO; copied into the
#          job's results directory by slurm/run_rev3_bo_rank.sh)
# Outputs: desurv_bo_results_k{K}_tcgacptac   desurv_cv_bayesopt() result; attr "rank_bo"
#          desurv_bo_finalists_k{K}_tcgacptac desurv_bo_finalists() result; configs gain
#                                             train_cindex_bo (see section 4)

message("=== Step 3b: rank-stratified BO (rev 3 amendment 1) ===")
source("code/00_helpers.R")
library(DeSurv)
source("R/bo_helpers.R")
source("R/reported_settings.R")
source("R/bo_setup.R")

BO_K <- as.integer(Sys.getenv("DESURV_BO_K", NA))
if (is.na(BO_K) || BO_K < 2L || BO_K > BO_K_UPPER) {
  stop("Set DESURV_BO_K to an integer in 2..", BO_K_UPPER, " (got '",
       Sys.getenv("DESURV_BO_K"), "').")
}
message("  k = ", BO_K)

# ── Budget (amendment 1) ──────────────────────────────────────────────────
# Every rank: RANK_N_INIT new Latin-hypercube points, then EI iterations until at least
# RANK_MIN_TOTAL valid evaluations, counting inherited ones: n_iter = max(0, 40 - existing - 10).
if (CONFIG$quick) {
  RANK_N_INIT    <- 2L
  RANK_MIN_TOTAL <- NA_integer_    # quick: n_iter below, budget rule not applied
  # Quick only: at least 1 EI iteration and at least 5 evaluations in total. With 4 points
  # in 4 dimensions the package's final GP fit failed (local test, k = 2: 1 inherited + 2 + 1),
  # and desurv_bo_finalists() then stops with "No GP posterior mean".
  RANK_QUICK_MIN_TOTAL <- 5L
} else {
  RANK_N_INIT    <- 10L
  RANK_MIN_TOTAL <- 40L
}

# Seed observations for k = 3 (amendment 1; values as in analysis/longleaf/bench/noise_bo_eval.R,
# CV C from job 3500640 at bank 4,000,000).
RANK_SEEDS <- data.frame(
  label       = c("published", "may_best"),
  k           = 3L,
  alpha_grid  = c(0.334, 0.442),
  lambda_grid = c(0.349, 0.2410),
  nu_grid     = c(0.0557, 0.2015),
  ntop        = c(270, 299),
  mean_cindex = c(0.6258, 0.5803),
  stringsAsFactors = FALSE
)

# ── Inputs ────────────────────────────────────────────────────────────────
tar_data  <- load_precomputed("tar_data_tcgacptac")
main_bo   <- load_precomputed("desurv_bo_results_tcgacptac")
main_run  <- if (inherits(main_bo, "desurv_cv_bo")) main_bo else main_bo$runs[[length(main_bo$runs)]]

# The main BO must have been run on these data with this script's CV settings.
fp <- DeSurv:::.desurv_bo_data_fingerprint(tar_data$ex, tar_data$sampInfo$time,
                                           tar_data$sampInfo$event, tar_data$sampInfo$dataset,
                                           tar_data$samp_keeps)
if (!isTRUE(all.equal(fp, main_run$data_fingerprint, tolerance = 1e-12))) {
  stop("tar_data differs from the data the main BO was run on.")
}

# Bounds and fixed values: k dropped from the search and pinned (code/03 elbow-k idiom).
bounds_rank <- BO_BOUNDS
bounds_rank$k_grid <- NULL
bo_fixed_rank <- c(BO_FIXED, list(k_grid = BO_K))
param_names <- names(bounds_rank)

# ── 1. History at this k ─────────────────────────────────────────────────
mh <- as.data.frame(main_run$history, stringsAsFactors = FALSE)
keep <- mh$status == "ok" & is.finite(mh$mean_cindex) &
  as.integer(round(as.numeric(mh$k_grid))) == BO_K
hist <- data.frame(
  main_eval_id = mh$eval_id[keep],
  mh[keep, param_names, drop = FALSE],
  mean_cindex  = mh$mean_cindex[keep],
  status       = "ok",
  stage        = mh$stage[keep],
  stringsAsFactors = FALSE
)
n_inherited <- nrow(hist)
n_seeds <- 0L
if (BO_K == 3L) {
  seeds <- data.frame(main_eval_id = NA_integer_, RANK_SEEDS[, param_names],
                      mean_cindex = RANK_SEEDS$mean_cindex, status = "ok", stage = "seed",
                      stringsAsFactors = FALSE)
  hist <- rbind(hist, seeds)
  n_seeds <- nrow(seeds)
}
rownames(hist) <- NULL
message(sprintf("  history: %d inherited main-BO evaluations at k = %d, %d seed rows",
                n_inherited, BO_K, n_seeds))

# Check the package imports every row (it drops rows outside the bounds and duplicate
# parameter keys and renumbers eval_id 1..n in row order). The rank BO's eval_id j
# (j <= n_existing) is then row j of `hist`; the train-C lookup below relies on that.
bound_info <- DeSurv:::.desurv_bo_parse_bounds(bounds_rank)
prep <- DeSurv:::.desurv_bo_prepare_history(if (nrow(hist)) hist else NULL, bound_info)
if (prep$n_loaded != nrow(hist)) {
  stop(sprintf("Package would import %d of %d history rows at k = %d (discarded %d); not proceeding.",
               prep$n_loaded, nrow(hist), BO_K, prep$n_discarded))
}
n_existing <- nrow(hist)

if (CONFIG$quick) {
  n_iter_rank <- max(1L, RANK_QUICK_MIN_TOTAL - n_existing - RANK_N_INIT)
} else {
  n_iter_rank <- max(0L, RANK_MIN_TOTAL - n_existing - RANK_N_INIT)
}
message(sprintf("  budget: %d existing + %d LHS + %d EI = %d", n_existing, RANK_N_INIT,
                n_iter_rank, n_existing + RANK_N_INIT + n_iter_rank))

# ── 2. Fixed-k BO ─────────────────────────────────────────────────────────
# The arguments desurv_cv_bayesopt_refine() passes to its coarse run (DeSurv
# R/desurv_bo_refine.R:282-298): data, bo_fixed and BO_COMMON minus the refine-only
# entries, modified by coarse_control, plus bo_bounds. Here also bo_history and this
# rank's n_init / n_iter; candidate_pool, exploration_weight and seed stay as in BO_COMMON.
REFINE_ONLY <- c("max_refinements", "tol_gain", "plateau", "top_k", "shrink_base",
                 "importance_gain", "coarse_control", "refine_control")
rank_control <- utils::modifyList(BO_COMMON$coarse_control,
                                  list(n_init = RANK_N_INIT, n_iter = n_iter_rank))

desurv_bo_results_k <- cache_or_compute(sprintf("desurv_bo_results_k%d_tcgacptac", BO_K), {
  args <- c(
    list(
      X = tar_data$ex,
      y = tar_data$sampInfo$time,
      d = tar_data$sampInfo$event,
      dataset = tar_data$sampInfo$dataset,
      samp_keeps = tar_data$samp_keeps,
      bo_fixed = bo_fixed_rank
    ),
    BO_COMMON[setdiff(names(BO_COMMON), REFINE_ONLY)]
  )
  args <- utils::modifyList(args, rank_control)
  args$bo_bounds  <- bounds_rank
  args$bo_history <- if (nrow(hist)) hist else NULL
  res <- do.call(DeSurv::desurv_cv_bayesopt, args)
  attr(res, "rank_bo") <- list(
    k = BO_K, n_inherited = n_inherited, n_seeds = n_seeds, n_existing = n_existing,
    n_init = RANK_N_INIT, n_iter = n_iter_rank,
    imported = hist,                         # row j = rank-BO eval_id j
    main_eval_ids = hist$main_eval_id,
    seeds = if (n_seeds) RANK_SEEDS else NULL
  )
  res
})

# The scoring arguments must equal the main BO's (k and execution-only settings aside):
# same estimator, partitions, seed banks, tolerances.
local({
  drop <- c("k_grid", "parallel_grid", "ncores_grid")
  a <- desurv_bo_results_k$cv_args; b <- main_run$cv_args
  a <- a[setdiff(names(a), drop)]; b <- b[setdiff(names(b), drop)]
  if (!identical(a[sort(names(a))], b[sort(names(b))])) {
    stop("Rank BO cv_args differ from the main BO's (k_grid, parallel_grid, ncores_grid excluded).")
  }
})

h <- desurv_bo_results_k$history
message(sprintf("  rank BO history: %d rows (%d ok); stages: %s", nrow(h),
                sum(h$status == "ok" & is.finite(h$mean_cindex)),
                paste(names(table(h$stage)), table(h$stage), sep = "=", collapse = ", ")))

# ── 3. Finalists: top FINALISTS_PER_K by this rank's GP posterior mean ──────
# desurv_bo_finalists() ranks every ok row of the history (inherited, seeds and new)
# by history$gp_mean, the final GP fitted on all of them; k is taken from cv_args$k_grid.
desurv_bo_finalists_k <- cache_or_compute(sprintf("desurv_bo_finalists_k%d_tcgacptac", BO_K), {
  fin <- run_bo_finalists(desurv_bo_results_k)

  # ── 4. Training C for the finalists (diagnostic only) ──
  # desurv_bo_finalists() keeps only held-out fold C, so training C is taken from the
  # BO evaluation of the same configuration (same partitions and seed bank): rank-BO
  # diagnostics for new evaluations, main-BO diagnostics (by main eval_id) for inherited
  # ones, NA for the two k = 3 seeds (evaluated outside both BOs).
  train_by_eval <- function(diag) {
    if (is.null(diag) || !"train_cindex" %in% names(diag)) return(numeric())
    tapply(as.numeric(diag$train_cindex), diag$eval_id, mean, na.rm = TRUE)
  }
  tr_rank <- train_by_eval(desurv_bo_results_k$diagnostics)
  tr_main <- train_by_eval(main_run$diagnostics)
  fin$configs$train_cindex_bo <- vapply(fin$configs$eval_id, function(id) {
    if (id <= n_existing) {
      mid <- hist$main_eval_id[id]
      if (is.na(mid)) NA_real_ else unname(tr_main[as.character(mid)])
    } else {
      unname(tr_rank[as.character(id)])
    }
  }, numeric(1))
  inherited <- fin$configs$eval_id <= n_existing
  fin$configs$source <- "rank_bo"
  fin$configs$source[inherited] <- ifelse(hist$stage[fin$configs$eval_id[inherited]] == "seed",
                                          "seed", "main_bo")
  fin$configs$main_eval_id <- NA_integer_
  fin$configs$main_eval_id[inherited] <- hist$main_eval_id[fin$configs$eval_id[inherited]]
  fin
})

print(desurv_bo_finalists_k$configs)
message("=== Step 3b complete (k = ", BO_K, ") ===")
