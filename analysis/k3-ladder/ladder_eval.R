# k = 3 reconciliation ladder (rev 3 amendment 1c, section 4): one CV evaluation of one configuration
# under one ladder step, through DeSurv::desurv_cv exactly as the BO calls it (cv_only = TRUE), so the
# score is comparable to a BO history value. No cache_or_compute, no results/ or figures/ writes.
#
# Usage:  LADDER_DATA=<tar_data_tcgacptac.rds> LADDER_OUT=<dir> Rscript analysis/k3-ladder/ladder_eval.R <config> <step>
#   (or LADDER_CONFIG / LADDER_STEP in the environment; or a Slurm array index into LADDER_MANIFEST)
# config: c136_k3 (May eval 136), pub_k3 (May eval 129 hyperparameters at k 3), pub_k7 (May eval 129)
# step (one ingredient at a time; S3p is the swapped order, S4' = S4 by construction):
#   S0  DeSurv 1.0.1, mean_starts (30 random starts, mean held-out C per fold), tol 1e-5, May partition (seed 123)
#   S1  as S0 on the rev 3 partitions cv_seeds 20261101-3 (start seeds derive from cv_seeds; the bank is unused)
#   S2  as S1 under DeSurv 1.2.0 (Breslow tie fix + deterministic top-gene tie-break)
#   S3  as S2 at tol 1e-7
#   S4  reported estimator (100 seed fits from bank 4,000,000 -> consensus -> final fit), tol 1e-7, tol_init 1e-5
#   S3p reported estimator at tol 1e-5 (S2 -> S3p -> S4 is the swapped path)
# Env: LADDER_NCORES (default SLURM_CPUS_PER_TASK, else 8); LADDER_SMOKE=1 -> 2 starts / 2 seed fits / 1 partition.
# Output: <LADDER_OUT>/<config>__<step>[__smoke].rds (refuses to overwrite).

suppressMessages(library(DeSurv))
t_start <- Sys.time()

args <- commandArgs(TRUE)
manifest <- Sys.getenv("LADDER_MANIFEST")
task <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (length(args) >= 2L) {
  cfg_name <- args[1L]; step <- args[2L]
} else if (nzchar(manifest) && nzchar(task)) {
  m <- utils::read.csv(manifest, stringsAsFactors = FALSE)
  row <- m[m$task == as.integer(task), ]
  stopifnot("array task id not in manifest" = nrow(row) == 1L)
  cfg_name <- row$config; step <- row$step
} else {
  cfg_name <- Sys.getenv("LADDER_CONFIG"); step <- Sys.getenv("LADDER_STEP")
}

# Exact May values (results/desurv_bo_results_tcgacptac.rds, runs[[1]]$history, eval_id 136 and 129;
# checked against the cached object by may_reference.R).
CONFIGS <- list(
  c136_k3 = list(k = 3L, alpha = 0.442207523710734685, lambda = 0.240956814604169955,
                 nu = 0.201480919147608795, ntop = 299L, may_eval = 136L),
  pub_k3  = list(k = 3L, alpha = 0.333935556177224480, lambda = 0.349448737941329357,
                 nu = 0.055651216804573779, ntop = 270L, may_eval = NA_integer_),
  pub_k7  = list(k = 7L, alpha = 0.333935556177224480, lambda = 0.349448737941329357,
                 nu = 0.055651216804573779, ntop = 270L, may_eval = 129L)
)
STEPS <- list(
  S0  = list(pkg = "1.0.1", estimator = "mean_starts", tol = 1e-5, partitions = "may"),
  S1  = list(pkg = "1.0.1", estimator = "mean_starts", tol = 1e-5, partitions = "rev3"),
  S2  = list(pkg = "1.2.0", estimator = "mean_starts", tol = 1e-5, partitions = "rev3"),
  S3  = list(pkg = "1.2.0", estimator = "mean_starts", tol = 1e-7, partitions = "rev3"),
  S4  = list(pkg = "1.2.0", estimator = "reported",    tol = 1e-7, partitions = "rev3"),
  S3p = list(pkg = "1.2.0", estimator = "reported",    tol = 1e-5, partitions = "rev3")
)
if (!cfg_name %in% names(CONFIGS)) stop("config must be one of ", paste(names(CONFIGS), collapse = ", "))
if (!step %in% names(STEPS)) stop("step must be one of ", paste(names(STEPS), collapse = ", "))
cf <- CONFIGS[[cfg_name]]; st <- STEPS[[step]]

# Library check: the step decides the DeSurv version; never run a step under the wrong one.
ver <- as.character(packageVersion("DeSurv")); pkg_path <- normalizePath(find.package("DeSurv"))
if (st$pkg == "1.0.1") {
  if (ver != "1.0.1") stop(step, " needs DeSurv 1.0.1; loaded ", ver, " from ", pkg_path)
} else if (packageVersion("DeSurv") < "1.2.0") {
  stop(step, " needs DeSurv >= 1.2.0; loaded ", ver, " from ", pkg_path)
}

smoke <- identical(Sys.getenv("LADDER_SMOKE"), "1")
ncores <- as.integer(Sys.getenv("LADDER_NCORES", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
data_path <- Sys.getenv("LADDER_DATA"); out_dir <- Sys.getenv("LADDER_OUT")
if (!nzchar(data_path) || !file.exists(data_path)) stop("set LADDER_DATA to tar_data_tcgacptac.rds")
if (!nzchar(out_dir)) stop("set LADDER_OUT")
abs_path <- function(p) {  # absolute, symlinks resolved for the existing part; works for not-yet-existing dirs
  if (!startsWith(p, "/")) p <- file.path(getwd(), p)
  tail_parts <- character(0)
  while (!file.exists(p)) { tail_parts <- c(basename(p), tail_parts); p <- dirname(p) }
  do.call(file.path, as.list(c(normalizePath(p), tail_parts)))
}
out_dir <- abs_path(out_dir)
repo <- abs_path(Sys.getenv("LADDER_REPO", getwd()))
for (bad in file.path(repo, c("results", "figures"))) {
  if (startsWith(paste0(out_dir, "/"), paste0(bad, "/"))) stop("LADDER_OUT inside ", bad, "; refusing")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("%s__%s%s.rds", cfg_name, step, if (smoke) "__smoke" else ""))
if (file.exists(out_file)) stop(out_file, " exists; refusing to overwrite")

N_STARTS    <- if (smoke) 2L else 30L     # May BO_FIXED$n_starts (public:code/03:58-63)
N_SEED_FITS <- if (smoke) 2L else 100L    # rev 3 R/reported_settings.R n_seed_fits
CV_SEEDS    <- c(20261101L, 20261102L, 20261103L)
if (smoke) CV_SEEDS <- CV_SEEDS[1L]
MAY_SEED    <- 123L                       # desurv_cv default `seed`; the May BO passed none

td <- readRDS(data_path)
common <- list(X = td$ex, y = td$sampInfo$time, d = td$sampInfo$event, dataset = td$sampInfo$dataset,
  samp_keeps = td$samp_keeps, preprocess = TRUE, ngene = 3000L, method_trans_train = "rank",
  nfolds = 5L, maxit = 4000L, engine = "warmstart", lambdaW_grid = 0, lambdaH_grid = 0,
  parallel_grid = ncores > 1L, ncores_grid = ncores, cv_only = TRUE, verbose = FALSE,
  k_grid = cf$k, alpha_grid = cf$alpha, lambda_grid = cf$lambda, nu_grid = cf$nu, ntop = cf$ntop,
  tol = st$tol)

# One desurv_cv call -> per (rep, fold) unit rows. For mean_starts the unit value is the mean over
# starts (as the BO summary), train C likewise; for reported there is one fit per unit.
units_of <- function(res, rep_id) {
  dg <- res$diagnostics
  if (!"rep" %in% names(dg)) dg$rep <- rep_id
  dg$val_cindex[!is.finite(dg$val_cindex)] <- NA
  u <- stats::aggregate(cbind(val_C = val_cindex, train_C = train_cindex) ~ rep + fold, data = dg,
                        FUN = function(x) mean(x, na.rm = TRUE), na.action = stats::na.pass)
  n_fits <- stats::aggregate(val_cindex ~ rep + fold, data = dg, FUN = length, na.action = stats::na.pass)
  u$n_fits <- n_fits$val_cindex[match(paste(u$rep, u$fold), paste(n_fits$rep, n_fits$fold))]
  u
}

calls <- list()
if (st$estimator == "mean_starts") {
  seeds_use <- if (st$partitions == "may") MAY_SEED else CV_SEEDS
  if (st$pkg == "1.0.1") {
    # 1.0.1 has no cv_seeds: one call per partition with seed = cv_seed, which is what 1.2.0's
    # run_mean_starts() does per partition (cv_warm.R and cv_helpers.R are unchanged 1.0.1 -> 1.2.0).
    for (r in seq_along(seeds_use)) {
      res <- do.call(desurv_cv, c(common, list(n_starts = N_STARTS, seed = seeds_use[r])))
      calls[[r]] <- list(res = res, units = units_of(res, r), score = res$summary$mean_cindex)
    }
  } else {
    a <- c(common, list(n_starts = N_STARTS, estimator = "mean_starts"))
    if (st$partitions == "may") a$seed <- MAY_SEED else a$cv_seeds <- seeds_use
    res <- do.call(desurv_cv, a)
    calls[[1L]] <- list(res = res, units = units_of(res, 1L), score = res$summary$mean_cindex)
  }
} else {
  res <- do.call(desurv_cv, c(common, list(estimator = "reported", cv_seeds = CV_SEEDS,
    n_seed_fits = N_SEED_FITS, seed_bank_base = 4000000L,
    reported_control = list(tol_init = 1e-5, imaxit = 4000L))))
  calls[[1L]] <- list(res = res, units = units_of(res, 1L), score = res$summary$mean_cindex)
}

units <- do.call(rbind, lapply(calls, `[[`, "units"))
units <- units[order(units$rep, units$fold), ]
mean_C <- mean(units$val_C, na.rm = TRUE)
score_pkg <- if (length(calls) == 1L) calls[[1L]]$score else NA_real_
# A single desurv_cv call's own score must equal the unit mean (all units finite).
if (length(calls) == 1L && all(is.finite(units$val_C))) stopifnot(abs(score_pkg - mean_C) < 1e-12)
minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

out <- list(
  config = cfg_name, step = step, params = cf, mean_C = mean_C, desurv_cv_score = score_pkg,
  se_unit = stats::sd(units$val_C, na.rm = TRUE) / sqrt(sum(is.finite(units$val_C))),
  mean_train_C = mean(units$train_C, na.rm = TRUE), units = units,
  folds = lapply(calls, function(x) x$res$folds),
  settings = c(st, list(n_starts = if (st$estimator == "mean_starts") N_STARTS else NA,
    n_seed_fits = if (st$estimator == "reported") N_SEED_FITS else NA,
    tol_init = if (st$estimator == "reported") 1e-5 else NA, imaxit = if (st$estimator == "reported") 4000L else NA,
    seed_bank_base = if (st$estimator == "reported") 4000000L else NA,
    seeds = if (st$partitions == "may") MAY_SEED else CV_SEEDS, maxit = 4000L, ngene = 3000L, nfolds = 5L,
    engine = "warmstart", smoke = smoke)),
  desurv_version = ver, desurv_path = pkg_path, r_version = R.version.string,
  data_path = normalizePath(data_path), data_md5 = unname(tools::md5sum(data_path)),
  ncores = ncores, node = Sys.info()[["nodename"]], job = Sys.getenv("SLURM_JOB_ID", "local"),
  minutes = minutes, finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
saveRDS(out, out_file)
message(sprintf("LADDER %s %s DeSurv %s C=%.4f train=%.4f units=%d minutes=%.1f -> %s", cfg_name, step, ver,
                mean_C, out$mean_train_C, nrow(units), minutes, out_file))
