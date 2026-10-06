# One tuning run (prespec section 5 "Tuning run"; spec section 7 candidate procedure with the reduced budget).
# Stage (prespec section 9.1): S3 = stage 2, S1 = stage 3, S2 = stage 4; in qual mode the job refuses to run unless
# the previous stage's evaluation reports a cumulative PASS (and, for S3, the stage 2 realized-ratio check passed).
# Usage: QUAL_MODE=qual|plumbing NCORES=<n> Rscript qual_tune_job.R <scenario> <replicate>
# 218 training columns of the replicate's sim_qual draw. Inner split: 5 folds stratified by event status
# (DeSurv 1.0.1 .desurv_make_folds_stratified, the comparator's fold function, dataset = NULL), seed = seed_tune.
# BO: bo_engine_copy (verbatim DeSurv 1.0.1 engine; candidate pool 4,000, exploration weight 0.01, BO seed =
# seed_tune) over the spec section 3 box: logit(alpha) in [logit(1e-4), logit(0.9)] (searched as the linear-scale
# parameter logit_alpha, because the engine has only linear and log10 scales), lambda [1e-3, 1e3] log10,
# nu [0, 0.95] linear, rho [1e-2, 1e3] log10, ntop [50, 300] integer. Budget n_init 20, n_iter 40 (plumbing 2, 1).
# Each evaluation: 5 inner-fold fits, 30 starts each (plumbing 3), seeds seed_inner(f) shared by every evaluation,
# s0 per inner-training set; inner C per fold = comparator scorer on the held-out fold with the fold's v~ on its
# ntop support; feasible only if all five fits are converged; infeasible = NA score + flag (excluded from the GP
# and from selection by the engine's status/NA filter). Selection: feasible evaluation with the highest mean
# inner-CV C, full precision. Final fit: n = 218, 30 starts, seed_final, recording on, scored on the test set.
# No feasible configuration = no tuned configuration = recorded Q1 failure. Writes one .rds under <OUT_DIR>/tune/.
QUAL <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(QUAL, "qual_common.R")); source(file.path(QUAL, "manifest.R")); source(file.path(QUAL, "bo_engine_copy.R"))
a <- commandArgs(TRUE)
if (length(a) != 2L) stop("usage: qual_tune_job.R <scenario> <replicate>", call. = FALSE)
scen <- a[1]; r <- as.integer(a[2])
NC <- as.integer(Sys.getenv("NCORES", "1"))
FOLD_PAR <- as.integer(Sys.getenv("QUAL_FOLD_PAR", "1")); stopifnot(FOLD_PAR %in% c(1L, 5L))
man <- qual_manifest(); bud <- man$budget
job <- man$tune[man$tune$scenario == scen & man$tune$r == r, ]
if (nrow(job) != 1L) stop(sprintf("tuning job %s r%d is not in the %s manifest", scen, r, QUAL_MODE), call. = FALSE)
if (file.exists(job$out)) stop(sprintf("output exists, not overwritten: %s", job$out), call. = FALSE)
require_stage_pass(job$stage); require_precheck(scen, job$stage)
t0 <- proc.time()[["elapsed"]]
z <- qual_data(scen, r); sc <- job$sc
D218 <- palm_data(z$X, z$y, z$d)
folds <- DeSurv:::.desurv_make_folds_stratified(z$d, dataset = NULL, nfolds = 5L, seed = seed_tune(sc, r))
stopifnot(attr(folds, "nfolds") == 5L)
FD <- lapply(1:5, function(f) palm_data(z$X[, folds != f], z$y[folds != f], z$d[folds != f]))
cat(sprintf("[%s] tuning job %s r%d: read 1 dataset (seed %d, p %d, train n %d, test %d); fold sizes %s; events per fold %s; cores %d\n",
            QUAL_MODE, scen, r, z$seed, D218$p, D218$n, ncol(z$X_test), paste(tabulate(folds), collapse = "/"),
            paste(tapply(z$d, folds, sum), collapse = "/"), NC))

BOUNDS <- list(logit_alpha = list(lower = qlogis(1e-4), upper = qlogis(0.9), type = "continuous", scale = "linear"),
               lambda = list(lower = 1e-3, upper = 1e3, scale = "log10"),
               nu = list(lower = 0, upper = 0.95, type = "continuous"),
               rho = list(lower = 1e-2, upper = 1e3, scale = "log10"),
               ntop = list(lower = 50, upper = 300, type = "integer"))
EVALS <- new.env(); EVALS$log <- list()
hp_of <- function(v) list(alpha = plogis(v$logit_alpha), lambda = v$lambda, nu = v$nu, rho = v$rho)
objective <- function(args) {
  v <- args[names(BOUNDS)]; hp <- hp_of(v); ntop <- as.integer(v$ntop)
  fits <- vector("list", 5L); cf <- rep(NA_real_, 5L)
  # QUAL_FOLD_PAR > 1 runs the 5 inner-fold fits in parallel (NC %/% QUAL_FOLD_PAR cores each for the starts); each fit
  # seeds its own starts, so the result does not depend on it (default 1 = folds in sequence).
  fits <- if (FOLD_PAR > 1L) parallel::mclapply(1:5, function(f) qual_fit(FD[[f]], hp, seed_inner(sc, r, f),
            ncores = max(1L, NC %/% FOLD_PAR), keep = "light"), mc.cores = FOLD_PAR) else
          lapply(1:5, function(f) qual_fit(FD[[f]], hp, seed_inner(sc, r, f), ncores = NC, keep = "light"))
  stopifnot(length(fits) == 5L, all(vapply(fits, is.list, TRUE)))
  for (f in 1:5) {
    rl <- fits[[f]]$ret_light
    if (!is.null(rl)) {
      va <- folds == f
      cf[f] <- qual_cindex(palm_score(z$X[, va, drop = FALSE], iterate_model(rl$W, rl$theta, FD[[f]], hp, ntop)), z$y[va], z$d[va])
    }
  }
  st <- vapply(fits, `[[`, "", "status")
  feasible <- all(st == "converged")
  score <- if (feasible) { x <- cf[is.finite(cf)]; if (length(x)) mean(x) else NA_real_ } else NA_real_
  i <- length(EVALS$log) + 1L
  EVALS$log[[i]] <- list(eval = i, values = v, hp = hp, ntop = ntop, feasible = feasible, infeasible_flag = !feasible,
                         score = score, fold_status = st, fold_C = cf, fits = fits)
  cat(sprintf("  eval %d: alpha %.4g lambda %.4g nu %.3f rho %.4g ntop %d | fold status %s | %s\n", i, hp$alpha, hp$lambda,
              hp$nu, hp$rho, ntop, paste(substr(st, 1, 4), collapse = ","),
              if (feasible) sprintf("mean inner C %.4f", score) else "INFEASIBLE (NA)"))
  list(summary = data.frame(mean_cindex = score), diagnostics = NULL)
}
bo <- tryCatch(bo_engine_copy(cv_objective = objective, X = NULL, y = NULL, d = NULL, bo_bounds = BOUNDS,
                              n_init = bud$n_init, n_iter = bud$n_iter, candidate_pool = 4000L, exploration_weight = 0.01,
                              seed = seed_tune(sc, r), verbose = FALSE),
               error = function(e) e)
log <- EVALS$log
ev <- data.frame(eval = seq_along(log), feasible = vapply(log, `[[`, TRUE, "feasible"), score = vapply(log, `[[`, 0, "score"))
ok <- ev$feasible & !is.na(ev$score)
res <- list(job = job, kind = "tune", mode = QUAL_MODE, hashes = HASHES, data = data_summary(z), folds = folds,
            budget = bud, bounds = BOUNDS, n_evals = length(log), evals = log,
            bo_error = if (inherits(bo, "error")) conditionMessage(bo) else NULL,
            bo_history = if (inherits(bo, "error")) NULL else bo$history)
if (!any(ok)) {
  res$tuned <- FALSE; res$q1_failure <- "no feasible configuration: the tuning run has no tuned configuration"
  cat(sprintf("NO FEASIBLE CONFIGURATION in %d evaluations: recorded as a Q1 failure (prespec section 5)\n", length(log)))
} else {
  sel <- which(ok)[which.max(ev$score[ok])]                                  # first maximum, as the engine
  bp <- bo$best$params
  stopifnot(isTRUE(all(unlist(log[[sel]]$values)[names(bp)] == bp)))        # engine's best = this selection
  res$tuned <- TRUE; res$selected <- sel; res$selected_values <- log[[sel]]$values
  res$selected_hp <- log[[sel]]$hp; res$selected_ntop <- log[[sel]]$ntop; res$selected_score <- log[[sel]]$score
  res$inner_fits_selected <- log[[sel]]$fits
  cat(sprintf("selected eval %d: alpha %.17g lambda %.17g nu %.17g rho %.17g ntop %d, mean inner C %.6f (%d of %d feasible)\n",
              sel, res$selected_hp$alpha, res$selected_hp$lambda, res$selected_hp$nu, res$selected_hp$rho,
              res$selected_ntop, res$selected_score, sum(ok), length(log)))
  res$final <- qual_fit(D218, res$selected_hp, seed_final(sc, r), ncores = NC, keep = "full", z = z, ntop = res$selected_ntop)
  cat(sprintf("final fit (n %d): status %s, returned start %s; %s\n", D218$n, res$final$status, res$final$best_start,
              fit_counts(res$final)))
}
res$wall <- proc.time()[["elapsed"]] - t0
cat(sprintf("evaluations %d (feasible %d), inner-fold fits %d, wall %.0f s\n", length(log), sum(ok), 5L * length(log), res$wall))
save_job(res, job$out)
