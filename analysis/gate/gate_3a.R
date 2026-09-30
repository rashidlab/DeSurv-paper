# Track 3a convergence probe, fixed-R rule.
# Governing documents: plan rev 2.1 (frozen 960b76d), section 3 clarifications (0dc8aff), and the
# section 3 amendment docs/plans/2026-09-30-section3-amendment-fixedR.md (d120e56), which supersedes the
# R ladder, the smallest-passing-R rule and the 2R rule.
#
# At the CURRENT configuration (tar_params_best_tcgacptac.rds, training data
# tar_data_filtered_tcgacptac.rds; no validation-cohort data is read):
#   Stage R = 50: four independent 50-run consensus-DeSurv ensembles from disjoint blocks
#     E1 = 1,000,001-050, E2 = 051-100, E3 = 101-150, E4 = 151-200.
#     For each of the 6 ensemble pairs: Hungarian alignment on correlation, per-factor Pearson r.
#     PASS       = every factor: median over the 6 pairs >= 0.98 and min >= 0.95; mode_pass TRUE in all 4;
#                  no ensemble fit error.
#     BORDERLINE = every factor median >= 0.95, mode_pass TRUE in all 4, no fit error, not PASS
#                  -> the same job runs stage R = 100 (1,000,201-300, 301-400, 401-500, 501-600),
#                  same PASS rule; a non-PASS at R = 100 is a FAIL.
#     Anything else = FAIL.
#   Descriptive only (never pass/fail): half-vs-half desurv_consensus_gate within each ensemble, beta,
#   training C, top-gene overlap, occupancy, tightness, outlier counts, and the current-estimator
#   reference (code/04 recipe on seeds 1,100,001-100 / 101-200 / 201-300).
#
# Environment:
#   GATE_DATA_DIR  directory holding the two rds inputs (read-only)
#   GATE_OUT_DIR   output directory (created; must not be inside any results/ directory)
#   SLURM_CPUS_PER_TASK  cores (default 1)
#   GATE_DRY=TRUE  correctness dry run: 4 ensembles of 6 runs (R = 12 stage: 4 x 12), 2 half-splits,
#                  one 4-seed reference set, at most 4 cores
#   GATE_FORCE_BORDERLINE=TRUE  dry run only: treat the first-stage verdict as BORDERLINE to exercise
#                  the second stage (refused unless GATE_DRY=TRUE)
#
# Never run code/NN_*.R or make from this job. This script only reads GATE_DATA_DIR.

t_start <- Sys.time()
suppressMessages(library(DeSurv))
stopifnot(packageVersion("DeSurv") == "1.1.0", exists("cox_stats_cpp", asNamespace("DeSurv")))

env_or_stop <- function(v) {
  x <- Sys.getenv(v, "")
  if (!nzchar(x)) stop(sprintf("Environment variable %s must be set.", v), call. = FALSE)
  x
}
DATA_DIR <- normalizePath(env_or_stop("GATE_DATA_DIR"), mustWork = TRUE)
OUT_DIR  <- env_or_stop("GATE_OUT_DIR")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_DIR  <- normalizePath(OUT_DIR, mustWork = TRUE)
if (grepl("(^|/)results(/|$)", OUT_DIR)) {
  stop("GATE_OUT_DIR must not be inside a results/ directory: ", OUT_DIR, call. = FALSE)
}
if (identical(OUT_DIR, DATA_DIR) || startsWith(OUT_DIR, paste0(DATA_DIR, "/"))) {
  stop("GATE_OUT_DIR must not be inside GATE_DATA_DIR.", call. = FALSE)
}
RUN_DIR <- file.path(OUT_DIR, "runs")
dir.create(RUN_DIR, showWarnings = FALSE)

DRY <- identical(toupper(Sys.getenv("GATE_DRY", "FALSE")), "TRUE")
FORCE_BORDERLINE <- identical(toupper(Sys.getenv("GATE_FORCE_BORDERLINE", "FALSE")), "TRUE")
if (FORCE_BORDERLINE && !DRY) stop("GATE_FORCE_BORDERLINE is allowed only with GATE_DRY=TRUE.", call. = FALSE)
NCORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(NCORES) || NCORES < 1L) NCORES <- 1L
if (DRY) NCORES <- min(NCORES, 4L)

# ---- Frozen settings (plan section 9 + fixed-R amendment) -----------------------------------------
N_ENS      <- 4L
MEDIAN_R   <- 0.98   # section 3 thresholds, unchanged by the amendment
MIN_R      <- 0.95
BORDER_R   <- 0.95   # amendment: borderline needs every factor median >= 0.95
SPLIT_SEED <- 20260930L   # descriptive half-vs-half diagnostic only
if (!DRY) {
  STAGES   <- list(list(R = 50L,  seeds = 1000001L:1000200L),
                   list(R = 100L, seeds = 1000201L:1000600L))
  N_SPLITS <- 10L
  REF_SETS <- list(set1 = 1100001L:1100100L, set2 = 1100101L:1100200L, set3 = 1100201L:1100300L)
  CHUNK    <- 50L
} else {
  STAGES   <- list(list(R = 6L,  seeds = 1000001L:1000024L),
                   list(R = 12L, seeds = 1000025L:1000072L))
  N_SPLITS <- 2L
  REF_SETS <- list(set1 = 1100001L:1100004L)
  CHUNK    <- 6L
}
for (st in STAGES) stopifnot(length(st$seeds) == N_ENS * st$R, st$R %% CHUNK == 0L)

# code/04 full-mode seed-fit settings (code/04_fit_models.R:35-37, 74-80, 109-114)
RUN_TOL   <- 1e-5
RUN_MAXIT <- 4000L
FIT_TOL      <- RUN_TOL / 100   # tol for seed fits, consensus refit and code/04 final fit
FIT_TOL_INIT <- RUN_TOL
FIT_MAXIT    <- RUN_MAXIT
FIT_IMAXIT   <- RUN_MAXIT

logf <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf(...), "\n", sep = "")
  flush.console()
}

# ---- Session info ---------------------------------------------------------------------------------
si <- capture.output(sessionInfo())
writeLines(c(sprintf("DeSurv %s from %s", packageVersion("DeSurv"), find.package("DeSurv")),
             sprintf("GATE_DRY=%s FORCE_BORDERLINE=%s NCORES=%d DATA_DIR=%s OUT_DIR=%s",
                     DRY, FORCE_BORDERLINE, NCORES, DATA_DIR, OUT_DIR),
             si), file.path(OUT_DIR, "sessionInfo.txt"))
cat(si, sep = "\n")
logf("DeSurv %s at %s; dry=%s; force_borderline=%s; cores=%d",
     packageVersion("DeSurv"), find.package("DeSurv"), DRY, FORCE_BORDERLINE, NCORES)

# ---- Inputs (training data only) ------------------------------------------------------------------
dat <- readRDS(file.path(DATA_DIR, "tar_data_filtered_tcgacptac.rds"))
par <- readRDS(file.path(DATA_DIR, "tar_params_best_tcgacptac.rds"))
extract_param <- function(params, name, default) {       # code/04_fit_models.R:55-58
  val <- params[[name]]
  if (is.null(val) || is.na(val)) default else val
}
X <- dat$ex; y <- dat$sampInfo$time; d <- dat$sampInfo$event
HP <- list(k = as.integer(par$k), alpha = par$alpha, lambda = par$lambda, nu = par$nu,
           lambdaW = extract_param(par, "lambdaW", 0), lambdaH = extract_param(par, "lambdaH", 0))
NTOP <- if (!is.null(par$ntop) && !is.na(par$ntop)) as.integer(round(par$ntop)) else 100L  # code/04:62-66
logf("Data: %d genes x %d samples (%s). Params: k=%d alpha=%.6g lambda=%.6g nu=%.6g lambdaW=%g lambdaH=%g ntop=%d",
     nrow(X), ncol(X), paste(names(table(dat$sampInfo$dataset)), table(dat$sampInfo$dataset), collapse = ", "),
     HP$k, HP$alpha, HP$lambda, HP$nu, HP$lambdaW, HP$lambdaH, NTOP)
input_md5 <- tools::md5sum(file.path(DATA_DIR, c("tar_data_filtered_tcgacptac.rds", "tar_params_best_tcgacptac.rds")))
rm(dat); invisible(gc())

# ---- One constituent run: identical call to code/04 run_seed_fits and to desurv_consensus_fit ------
fit_one <- function(s) {
  f <- tryCatch(
    desurv_fit(X = X, y = y, d = d, k = HP$k, alpha = HP$alpha, lambda = HP$lambda, nu = HP$nu,
               lambdaW = HP$lambdaW, lambdaH = HP$lambdaH,
               seed = s, tol = FIT_TOL, tol_init = FIT_TOL_INIT,
               maxit = FIT_MAXIT, imaxit = FIT_IMAXIT, ninit = 1L,
               parallel_init = FALSE, verbose = FALSE),
    error = function(e) structure(conditionMessage(e), class = "gate_run_error"))
  if (inherits(f, "desurv_fit")) {
    f$data <- NULL      # a copy of the training data (~4 MB per run); not used by either estimator
    f$seed <- s
  }
  f
}
run_ok <- function(f) {
  inherits(f, "desurv_fit") && !isTRUE(f$nan_flag) && is.matrix(f$W) && all(is.finite(f$W))
}

# Fit a seed bank in chunks; each chunk is saved atomically and skipped on restart.
# Worker crashes (mclapply returns try-error/NULL) are retried once sequentially; a run that still
# fails is recorded and returned as a failure (caller decides: the gate bank stops, code/04 drops).
fit_bank <- function(seeds, tag) {
  chunks <- split(seeds, ceiling(seq_along(seeds) / CHUNK))
  out <- vector("list", length(chunks))
  for (ci in seq_along(chunks)) {
    s <- chunks[[ci]]
    fn <- file.path(RUN_DIR, sprintf("%s_%d-%d.rds", tag, s[1], s[length(s)]))
    if (file.exists(fn)) {
      out[[ci]] <- readRDS(fn)
      ok <- vapply(out[[ci]], run_ok, TRUE)
      stopifnot(length(out[[ci]]) == length(s),
                identical(vapply(out[[ci]][ok], function(f) as.integer(f$seed), integer(1)), s[ok]))
      logf("%s chunk %d/%d (seeds %d-%d): loaded checkpoint", tag, ci, length(chunks), s[1], s[length(s)])
      next
    }
    t0 <- Sys.time()
    res <- if (NCORES > 1L) {
      parallel::mclapply(s, fit_one, mc.cores = min(NCORES, length(s)), mc.preschedule = FALSE)
    } else lapply(s, fit_one)
    crashed <- which(!vapply(res, function(f) inherits(f, "desurv_fit") || inherits(f, "gate_run_error"), TRUE))
    for (i in crashed) {
      logf("  seed %d: worker crash, refitting sequentially", s[i])
      res[[i]] <- fit_one(s[i])
    }
    tmp <- paste0(fn, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, fn)
    out[[ci]] <- res
    logf("%s chunk %d/%d (seeds %d-%d): %d ok, %d failed, %.1f min", tag, ci, length(chunks), s[1], s[length(s)],
         sum(vapply(res, run_ok, TRUE)), sum(!vapply(res, run_ok, TRUE)),
         as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  res <- unlist(out, recursive = FALSE)
  stopifnot(length(res) == length(seeds))
  res
}

align_perm <- function(A, B) DeSurv:::.dc_align(A, B)   # B[, perm] matches A (Hungarian on Pearson r)
top_genes <- function(W) {
  tg <- desurv_get_top_genes(W, NTOP)$top_genes
  lapply(seq_len(ncol(W)), function(j) tg[[j]])
}

# ---- One gate stage: fit the bank, 4 ensembles, pairwise agreement, verdict ------------------------
cons_args <- c(HP, list(tol = FIT_TOL, maxit = FIT_MAXIT))
run_stage <- function(st) {
  R <- st$R; seeds <- st$seeds
  logf("Stage R = %d: fitting %d runs (seeds %d-%d)", R, length(seeds), seeds[1], seeds[length(seeds)])
  bank <- fit_bank(seeds, "bank")
  bad <- !vapply(bank, run_ok, TRUE)
  if (any(bad)) {
    # A failed constituent run is never dropped; it is an ensemble fit error (gate FAIL), recorded.
    msgs <- vapply(bank[bad], function(f) if (is.character(f)) f[1] else "invalid fit", "")
    writeLines(sprintf("seed %d: %s", seeds[bad], msgs), file.path(OUT_DIR, sprintf("bank_failures_R%d.txt", R)))
    logf("  %d constituent run(s) failed; see bank_failures_R%d.txt", sum(bad), R)
  }
  blocks <- split(seq_along(seeds), rep(seq_len(N_ENS), each = R))
  names(blocks) <- paste0("E", seq_len(N_ENS))

  # Consensus fit and descriptive half-vs-half diagnostic per ensemble (parallel over ensembles)
  one_ens <- function(e) {
    idx <- blocks[[e]]
    if (any(bad[idx])) {
      msg <- sprintf("constituent run(s) failed: seeds %s", paste(seeds[idx][bad[idx]], collapse = ", "))
      return(list(fit = structure(msg, class = "gate_cons_error"), hvh = NULL))
    }
    fit <- tryCatch(
      do.call(desurv_consensus_fit, c(list(X = X, y = y, d = d, run_fits = bank[idx], seeds = seeds[idx]),
                                      cons_args)),
      error = function(e) structure(conditionMessage(e), class = "gate_cons_error"))
    hvh <- tryCatch(
      do.call(desurv_consensus_gate, c(list(run_fits = bank[idx], X = X, y = y, d = d, R_values = R,
                                            n_splits = N_SPLITS, split_seed = SPLIT_SEED), cons_args)),
      error = function(e) structure(conditionMessage(e), class = "gate_cons_error"))
    list(fit = fit, hvh = hvh)
  }
  t0 <- Sys.time()
  ens <- if (NCORES > 1L) {
    parallel::mclapply(names(blocks), one_ens, mc.cores = min(NCORES, N_ENS), mc.preschedule = FALSE)
  } else lapply(names(blocks), one_ens)
  names(ens) <- names(blocks)
  for (e in names(ens)) {   # worker crash: redo sequentially
    if (!is.list(ens[[e]]) || is.null(ens[[e]]$fit)) ens[[e]] <- one_ens(e)
  }
  logf("  consensus fits + half-vs-half diagnostics: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

  fits <- lapply(ens, `[[`, "fit")
  fit_err <- vapply(fits, function(f) !inherits(f, "desurv_consensus_fit"), TRUE)
  mode_ok <- vapply(fits, function(f) inherits(f, "desurv_consensus_fit") && isTRUE(f$consensus$mode_pass), TRUE)
  k <- HP$k

  # Factor identity: each ensemble mapped to E1 by Hungarian; per pair, the pair's own Hungarian
  # alignment gives the partner factor (the amendment's per-pair alignment).
  pairs <- NULL; to_ref <- NULL; tg_over <- NULL
  if (!any(fit_err)) {
    W <- lapply(fits, `[[`, "W")
    to_ref <- lapply(W, function(Wb) align_perm(W[[1]], Wb))   # ensemble column to_ref[[e]][j] ~ E1 factor j
    tg <- lapply(W, top_genes)
    cmb <- utils::combn(names(W), 2L)
    rows <- list()
    for (i in seq_len(ncol(cmb))) {
      a <- cmb[1, i]; b <- cmb[2, i]
      p_ab <- align_perm(W[[a]], W[[b]])                      # b column p_ab[ia] ~ a column ia
      for (j in seq_len(k)) {
        ia <- to_ref[[a]][j]; ib <- p_ab[ia]
        rows[[length(rows) + 1L]] <- data.frame(
          pair = paste(a, b, sep = "-"), factor = j, col_a = ia, col_b = ib,
          consistent_with_E1_map = (ib == to_ref[[b]][j]),
          pearson = cor(W[[a]][, ia], W[[b]][, ib]),
          spearman = cor(W[[a]][, ia], W[[b]][, ib], method = "spearman"),
          beta_a = fits[[a]]$beta[ia], beta_b = fits[[b]]$beta[ib],
          topgene_overlap = length(intersect(tg[[a]][[ia]], tg[[b]][[ib]])) / NTOP,
          stringsAsFactors = FALSE)
      }
    }
    pairs <- do.call(rbind, rows)
  }

  per_factor <- if (!is.null(pairs)) {
    data.frame(factor = seq_len(k),
               median_r = vapply(seq_len(k), function(j) median(pairs$pearson[pairs$factor == j]), 0),
               min_r = vapply(seq_len(k), function(j) min(pairs$pearson[pairs$factor == j]), 0),
               median_topgene_overlap = vapply(seq_len(k), function(j)
                 median(pairs$topgene_overlap[pairs$factor == j]), 0))
  } else NULL

  reasons <- character(0)
  if (any(fit_err)) reasons <- c(reasons, sprintf("ensemble fit error in %s", paste(names(fits)[fit_err], collapse = ", ")))
  if (any(!mode_ok & !fit_err)) reasons <- c(reasons, sprintf("competing mode (mode_pass FALSE) in %s",
                                                               paste(names(fits)[!mode_ok & !fit_err], collapse = ", ")))
  verdict <- if (length(reasons)) {
    "FAIL"
  } else if (all(per_factor$median_r >= MEDIAN_R) && all(per_factor$min_r >= MIN_R)) {
    "PASS"
  } else if (all(per_factor$median_r >= BORDER_R)) {
    reasons <- sprintf("agreement below PASS: factor(s) %s",
                       paste(per_factor$factor[per_factor$median_r < MEDIAN_R | per_factor$min_r < MIN_R], collapse = ", "))
    "BORDERLINE"
  } else {
    reasons <- sprintf("factor median r < %.2f for factor(s) %s", BORDER_R,
                       paste(per_factor$factor[per_factor$median_r < BORDER_R], collapse = ", "))
    "FAIL"
  }
  logf("  Stage R = %d verdict: %s%s", R, verdict, if (length(reasons)) paste0(" (", paste(reasons, collapse = "; "), ")") else "")
  list(R = R, seeds = seeds, blocks = lapply(blocks, function(i) seeds[i]), verdict = verdict, reasons = reasons,
       per_factor = per_factor, pairs = pairs, to_E1 = to_ref, fit_error = fit_err, mode_pass = mode_ok,
       ensembles = fits, hvh = lapply(ens, `[[`, "hvh"),
       bank_diag = data.frame(seed = seeds, ok = !bad,
                              converged = vapply(bank, function(f) run_ok(f) && isTRUE(f$convergence), TRUE),
                              iters = vapply(bank, function(f) if (run_ok(f)) length(f$lossit) else NA_integer_, 1L),
                              objective = vapply(bank, function(f) if (run_ok(f)) as.numeric(f$loss$loss)[1] else NA_real_, 0),
                              cindex = vapply(bank, function(f) if (run_ok(f)) as.numeric(f$cindex)[1] else NA_real_, 0)))
}

# ---- Gate ----------------------------------------------------------------------------------------
stage_res <- list()
stage_res[["R_first"]] <- run_stage(STAGES[[1]])
first_verdict <- stage_res[["R_first"]]$verdict
if (FORCE_BORDERLINE) {
  logf("DRY RUN: forcing first-stage verdict %s -> BORDERLINE to exercise the second stage", first_verdict)
  first_verdict <- "BORDERLINE"
}
if (first_verdict == "BORDERLINE") {
  stage_res[["R_second"]] <- run_stage(STAGES[[2]])
}
names(stage_res) <- paste0("R", vapply(stage_res, `[[`, 1L, "R"))
s1 <- stage_res[[1]]
VERDICT <- if (first_verdict == "PASS") {
  sprintf("PASS at R = %d", s1$R)
} else if (first_verdict == "BORDERLINE") {
  s2 <- stage_res[[2]]
  if (s2$verdict == "PASS") sprintf("BORDERLINE at R = %d -> PASS at R = %d", s1$R, s2$R)
  else sprintf("FAIL: BORDERLINE at R = %d (%s), then %s at R = %d (%s)", s1$R, paste(s1$reasons, collapse = "; "),
               s2$verdict, s2$R, paste(s2$reasons, collapse = "; "))
} else {
  sprintf("FAIL at R = %d: %s", s1$R, paste(s1$reasons, collapse = "; "))
}
if (FORCE_BORDERLINE) VERDICT <- paste("[DRY RUN, borderline FORCED; not a gate result]", VERDICT)
logf("VERDICT: %s", VERDICT)
saveRDS(stage_res, file.path(OUT_DIR, "gate_stages.rds"))

# ---- Current-estimator reference (code/04 recipe; descriptive) ---------------------------------------
# code/04 consensus_and_fit (code/04_fit_models.R:104-117); failed seed fits dropped as in run_seed_fits.
code04_final <- function(fits) {
  init_vals <- DeSurv::desurv_consensus_seed(
    fits = fits, X = X, ntop = NTOP, k = HP$k,
    min_frequency = max(1L, as.integer(ceiling(0.3 * length(fits)))))
  f <- desurv_fit(X = X, y = y, d = d, k = HP$k, alpha = HP$alpha, lambda = HP$lambda, nu = HP$nu,
                  lambdaW = HP$lambdaW, lambdaH = HP$lambdaH,
                  W0 = init_vals$W0, H0 = init_vals$H0, beta0 = init_vals$beta0,
                  seed = NULL, tol = FIT_TOL, maxit = FIT_MAXIT, verbose = FALSE)
  f$data <- NULL
  f
}
logf("Reference: %d seed set(s) through the code/04 recipe", length(REF_SETS))
ref <- lapply(names(REF_SETS), function(nm) {
  s <- REF_SETS[[nm]]
  fits <- fit_bank(s, paste0("ref_", nm))
  ok <- vapply(fits, run_ok, TRUE)
  if (!any(ok)) stop("No successful seed fits for reference ", nm, call. = FALSE)
  logf("  %s: %d/%d successful seed fits", nm, sum(ok), length(s))
  fin <- code04_final(fits[ok])
  logf("  %s: final fit cindex = %.4f", nm, fin$cindex)
  list(seeds = s, n_ok = sum(ok), failed_seeds = s[!ok],
       seed_cindex = vapply(fits[ok], function(f) as.numeric(f$cindex), 0),
       final = fin)
})
names(ref) <- names(REF_SETS)
ref_pairs <- if (length(ref) >= 2L) {
  cmb <- utils::combn(names(ref), 2L)
  do.call(rbind, lapply(seq_len(ncol(cmb)), function(i) {
    A <- ref[[cmb[1, i]]]$final$W; B <- ref[[cmb[2, i]]]$final$W
    p <- align_perm(A, B)
    tA <- top_genes(A); tB <- top_genes(B)
    data.frame(pair = paste(cmb[1, i], cmb[2, i], sep = "-"), factor = seq_len(ncol(A)), matched_col = p,
               pearson = vapply(seq_len(ncol(A)), function(j) cor(A[, j], B[, p[j]]), 0),
               spearman = vapply(seq_len(ncol(A)), function(j) cor(A[, j], B[, p[j]], method = "spearman"), 0),
               topgene_overlap = vapply(seq_len(ncol(A)), function(j) length(intersect(tA[[j]], tB[[p[j]]])) / NTOP, 0))
  }))
} else NULL

# ---- Save everything ---------------------------------------------------------------------------------
results <- list(
  verdict = VERDICT,
  stages = stage_res,
  reference = ref,
  reference_pairs = ref_pairs,
  settings = list(dry = DRY, force_borderline = FORCE_BORDERLINE, ncores = NCORES, hyper = HP, ntop = NTOP,
                  tol = FIT_TOL, tol_init = FIT_TOL_INIT, maxit = FIT_MAXIT, imaxit = FIT_IMAXIT,
                  stages = STAGES, n_ensembles = N_ENS, median_r = MEDIAN_R, min_r = MIN_R, border_r = BORDER_R,
                  hvh_n_splits = N_SPLITS, hvh_split_seed = SPLIT_SEED, ref_sets = REF_SETS,
                  rule = "docs/plans/2026-09-30-section3-amendment-fixedR.md (d120e56)",
                  input_md5 = input_md5,
                  desurv_version = as.character(packageVersion("DeSurv")),
                  desurv_path = find.package("DeSurv"),
                  started = t_start, finished = Sys.time()),
  sessionInfo = sessionInfo()
)
saveRDS(results, file.path(OUT_DIR, "gate3a_results.rds"))

# ---- Text summary ------------------------------------------------------------------------------------
fmt <- function(x, d = 3) paste(formatC(x, format = "f", digits = d), collapse = " ")
L <- c(sprintf("Track 3a convergence probe (fixed-R rule)%s. DeSurv %s. %s", if (DRY) " (DRY RUN)" else "",
               packageVersion("DeSurv"), format(Sys.time())),
       sprintf("Config: k=%d alpha=%.6g lambda=%.6g nu=%.6g lambdaW=%g lambdaH=%g ntop=%d; tol=%g tol_init=%g maxit=%d",
               HP$k, HP$alpha, HP$lambda, HP$nu, HP$lambdaW, HP$lambdaH, NTOP, FIT_TOL, FIT_TOL_INIT, FIT_MAXIT),
       "",
       sprintf("VERDICT: %s", VERDICT),
       sprintf("Rule: per factor, median over the 6 ensemble pairs >= %.2f and min >= %.2f; mode_pass in all 4; no fit error. Borderline: all medians >= %.2f.",
               MEDIAN_R, MIN_R, BORDER_R))
for (st in stage_res) {
  L <- c(L, "", sprintf("=== Stage R = %d: %s%s ===", st$R, st$verdict,
                        if (length(st$reasons)) paste0(" (", paste(st$reasons, collapse = "; "), ")") else ""))
  L <- c(L, sprintf("Runs: %d fitted (seeds %d-%d), %d failed, %d converged, median iters %s",
                    nrow(st$bank_diag), st$seeds[1], st$seeds[length(st$seeds)], sum(!st$bank_diag$ok),
                    sum(st$bank_diag$converged), format(median(st$bank_diag$iters, na.rm = TRUE))))
  if (!is.null(st$per_factor)) {
    L <- c(L, "Per factor (factor ids = E1 columns): median r | min r | median top-gene overlap",
           sprintf("  factor %d: %.4f | %.4f | %.3f", st$per_factor$factor, st$per_factor$median_r,
                   st$per_factor$min_r, st$per_factor$median_topgene_overlap),
           "Pair table: pair factor col_a col_b pearson spearman topgene_overlap beta_a beta_b consistent_with_E1_map")
    pt <- st$pairs
    L <- c(L, sprintf("  %s f%d %d %d %.4f %.4f %.3f %.4f %.4f %s", pt$pair, pt$factor, pt$col_a, pt$col_b,
                      pt$pearson, pt$spearman, pt$topgene_overlap, pt$beta_a, pt$beta_b, pt$consistent_with_E1_map))
  }
  L <- c(L, "Per ensemble (columns in the ensemble's own canonical order):")
  for (e in names(st$ensembles)) {
    cf <- st$ensembles[[e]]
    if (!inherits(cf, "desurv_consensus_fit")) { L <- c(L, sprintf("  %s: FIT ERROR: %s", e, cf)); next }
    cs <- cf$consensus; occ <- cs$occupancy_summary
    L <- c(L, sprintf("  %s (seeds %d-%d): mode_pass=%s; map to E1 factors=%s; beta=%s; training C=%.4f",
                      e, st$blocks[[e]][1], st$blocks[[e]][length(st$blocks[[e]])], cs$mode_pass,
                      paste(match(seq_len(HP$k), st$to_E1[[e]]), collapse = " "), fmt(cf$beta, 4), cf$cindex),
              sprintf("      tightness=%s; retained=%s; n_members=%s; n_outliers=%s",
                      fmt(cs$tightness), fmt(cs$tightness_retained),
                      paste(cs$n_members, collapse = " "), paste(cs$n_outliers, collapse = " ")),
              sprintf("      occupancy (runs with 0/1/>1 factors in family): %s",
                      paste(sprintf("f%d %d/%d/%d", seq_len(ncol(occ)), occ["n0", ], occ["n1", ], occ["n_gt1", ]),
                            collapse = "; ")),
              sprintf("      competing-mode check: medoid r=%s; sub-cluster sizes=%s",
                      fmt(cs$competing_mode$medoid_r), paste(sprintf("%d/%d", cs$competing_mode$size_1,
                                                                     cs$competing_mode$size_2), collapse = " ")))
    h <- st$hvh[[e]]
    if (inherits(h, "desurv_consensus_gate")) {
      L <- c(L, sprintf("      [descriptive] half-vs-half (%d splits): per-factor median r %s; p05 r %s; half-level mode flags %s; errors %d",
                        N_SPLITS, fmt(h$per_factor$median_r, 4), fmt(h$per_factor$p05_r, 4),
                        format(h$summary$half_mode_fails), h$summary$n_errors))
    } else if (!is.null(h)) {
      L <- c(L, sprintf("      [descriptive] half-vs-half: ERROR %s", h))
    }
  }
}
L <- c(L, "", "CURRENT-ESTIMATOR REFERENCE (descriptive; code/04 recipe: seed fits -> desurv_consensus_seed -> one final fit)")
for (nm in names(ref)) {
  L <- c(L, sprintf("  %s (seeds %d-%d): %d/%d seed fits ok; final training C=%.4f; beta=%s",
                    nm, ref[[nm]]$seeds[1], ref[[nm]]$seeds[length(ref[[nm]]$seeds)], ref[[nm]]$n_ok,
                    length(ref[[nm]]$seeds), ref[[nm]]$final$cindex, fmt(ref[[nm]]$final$beta, 4)))
}
if (!is.null(ref_pairs)) {
  for (p in unique(ref_pairs$pair)) {
    z <- ref_pairs[ref_pairs$pair == p, ]
    L <- c(L, sprintf("  %s: pearson %s | spearman %s | top-gene overlap %s", p, fmt(z$pearson, 4),
                      fmt(z$spearman, 4), fmt(z$topgene_overlap, 3)))
  }
} else L <- c(L, "  (one reference set only: no pairwise correlations)")
L <- c(L, "", sprintf("Elapsed: %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
writeLines(L, file.path(OUT_DIR, "gate3a_summary.txt"))
cat(L, sep = "\n")
logf("Done. Outputs in %s", OUT_DIR)
