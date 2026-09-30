# Track 3a convergence probe (plan rev 2.1, FROZEN at 960b76d; section 3 + 2026-09-30 clarifications).
#
# At the CURRENT configuration (tar_params_best_tcgacptac.rds, training data
# tar_data_filtered_tcgacptac.rds; no validation-cohort data is read):
#   A. Gate run bank: seeds 1,000,001..1,000,400, one desurv_fit per seed (code/04 seed-fit recipe),
#      checkpointed in chunks so a killed job resumes.
#   B. desurv_consensus_gate on the bank prefixes (split_seed 20260930, 10 splits, package-default
#      frozen thresholds), plus desurv_consensus_fit at every R (objects saved).
#   C. Current-estimator reference: code/04 recipe (seed fits -> desurv_consensus_seed -> one final
#      desurv_fit) on seed sets 1,100,001-100 / 101-200 / 201-300; pairwise W correlations of the three
#      endpoints (factors aligned by Hungarian on Pearson r; Pearson and Spearman reported).
#
# Environment:
#   GATE_DATA_DIR  directory holding the two rds inputs (read-only)
#   GATE_OUT_DIR   output directory (created; must not be inside any results/ directory)
#   SLURM_CPUS_PER_TASK  cores (default 1)
#   GATE_DRY=TRUE  correctness dry run: seeds 1,000,001-012, R 4 and 6, 2 splits, one 4-seed
#                  reference set, at most 4 cores
#   GATE_R_VALUES  optional override of the gate R values (comma-separated even integers); see README
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
NCORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(NCORES) || NCORES < 1L) NCORES <- 1L
if (DRY) NCORES <- min(NCORES, 4L)

# ---- Frozen settings (plan section 9) -------------------------------------------------------------
SPLIT_SEED <- 20260930L
if (!DRY) {
  BANK_SEEDS  <- 1000001L:1000400L
  R_DESCRIBE  <- c(25L, 50L, 100L, 200L, 400L)   # desurv_consensus_fit at each (400 = 2R check at 200)
  R_GATE      <- c(50L, 100L, 200L, 400L)        # desurv_consensus_gate requires even R (see README)
  N_SPLITS    <- 10L
  REF_SETS    <- list(set1 = 1100001L:1100100L, set2 = 1100101L:1100200L, set3 = 1100201L:1100300L)
  CHUNK       <- 50L
} else {
  BANK_SEEDS  <- 1000001L:1000012L
  R_DESCRIBE  <- c(4L, 6L)
  R_GATE      <- c(4L, 6L)
  N_SPLITS    <- 2L
  REF_SETS    <- list(set1 = 1100001L:1100004L)
  CHUNK       <- 4L
}
if (nzchar(Sys.getenv("GATE_R_VALUES", ""))) {
  R_GATE <- as.integer(strsplit(Sys.getenv("GATE_R_VALUES"), ",")[[1]])
  R_DESCRIBE <- sort(unique(c(R_DESCRIBE, R_GATE)))
}
stopifnot(!anyNA(R_GATE), all(R_GATE %% 2L == 0L), max(c(R_GATE, R_DESCRIBE)) <= length(BANK_SEEDS))

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
             sprintf("GATE_DRY=%s NCORES=%d DATA_DIR=%s OUT_DIR=%s", DRY, NCORES, DATA_DIR, OUT_DIR),
             si), file.path(OUT_DIR, "sessionInfo.txt"))
cat(si, sep = "\n")
logf("DeSurv %s at %s; dry=%s; cores=%d", packageVersion("DeSurv"), find.package("DeSurv"), DRY, NCORES)

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
      stopifnot(identical(vapply(out[[ci]], function(f) if (inherits(f, "desurv_fit")) f$seed else NA_integer_,
                                 integer(1))[vapply(out[[ci]], run_ok, TRUE)], s[vapply(out[[ci]], run_ok, TRUE)]))
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

# ---- A. Gate run bank -------------------------------------------------------------------------------
logf("A. Gate run bank: %d seeds (%d-%d)", length(BANK_SEEDS), BANK_SEEDS[1], BANK_SEEDS[length(BANK_SEEDS)])
bank <- fit_bank(BANK_SEEDS, "bank")
bad <- !vapply(bank, run_ok, TRUE)
if (any(bad)) {
  msgs <- vapply(bank[bad], function(f) if (is.character(f)) f[1] else "invalid fit", "")
  writeLines(sprintf("seed %d: %s", BANK_SEEDS[bad], msgs), file.path(OUT_DIR, "bank_failures.txt"))
  stop(sprintf("%d gate-bank run(s) failed (seeds %s). No run is dropped; see bank_failures.txt.",
               sum(bad), paste(head(BANK_SEEDS[bad], 10), collapse = ", ")), call. = FALSE)
}
bank_conv <- vapply(bank, function(f) isTRUE(f$convergence), TRUE)
bank_iters <- vapply(bank, function(f) length(f$lossit), integer(1))
logf("Bank complete: %d runs, %d converged (non-converged runs are kept, per Track 1 design)",
     length(bank), sum(bank_conv))

# ---- B. Consensus fits at each R and the gate --------------------------------------------------------
cons_args <- c(HP, list(tol = FIT_TOL, maxit = FIT_MAXIT))
cfit_R <- function(R) {
  tryCatch(
    do.call(desurv_consensus_fit, c(list(X = X, y = y, d = d, run_fits = bank[seq_len(R)],
                                         seeds = BANK_SEEDS[seq_len(R)]), cons_args)),
    error = function(e) structure(conditionMessage(e), class = "gate_cons_error"))
}
logf("B1. desurv_consensus_fit at R = %s", paste(R_DESCRIBE, collapse = ", "))
t0 <- Sys.time()
cons_fits <- if (NCORES > 1L) {
  parallel::mclapply(R_DESCRIBE, cfit_R, mc.cores = min(NCORES, length(R_DESCRIBE)), mc.preschedule = FALSE)
} else lapply(R_DESCRIBE, cfit_R)
names(cons_fits) <- paste0("R", R_DESCRIBE)
for (i in which(!vapply(cons_fits, function(z) inherits(z, "desurv_consensus_fit") ||
                                        inherits(z, "gate_cons_error"), TRUE))) {
  cons_fits[[i]] <- cfit_R(R_DESCRIBE[i])   # worker crash: refit sequentially
}
saveRDS(cons_fits, file.path(OUT_DIR, "consensus_fits.rds"))
logf("B1 done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

logf("B2. desurv_consensus_gate: R = %s, n_splits = %d, split_seed = %d, thresholds = package defaults",
     paste(R_GATE, collapse = ", "), N_SPLITS, SPLIT_SEED)
t0 <- Sys.time()
gate <- do.call(desurv_consensus_gate,
                c(list(run_fits = bank, X = X, y = y, d = d, R_values = R_GATE,
                       n_splits = N_SPLITS, split_seed = SPLIT_SEED), cons_args))
saveRDS(gate, file.path(OUT_DIR, "gate.rds"))
logf("B2 done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

# ---- C. Current-estimator reference (code/04 recipe) -----------------------------------------------
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
logf("C. Reference: %d seed set(s) through the code/04 recipe", length(REF_SETS))
ref <- lapply(names(REF_SETS), function(nm) {
  s <- REF_SETS[[nm]]
  fits <- fit_bank(s, paste0("ref_", nm))
  ok <- vapply(fits, run_ok, TRUE)
  if (!any(ok)) stop("No successful seed fits for reference ", nm, call. = FALSE)
  logf("  %s: %d/%d successful seed fits", nm, sum(ok), length(s))
  fin <- code04_final(fits[ok])
  logf("  %s: final fit cindex = %.4f, beta = %s", nm, fin$cindex, paste(sprintf("%.3f", fin$beta), collapse = ", "))
  list(seeds = s, n_ok = sum(ok), failed_seeds = s[!ok],
       seed_cindex = vapply(fits[ok], function(f) as.numeric(f$cindex), 0),
       final = fin)
})
names(ref) <- names(REF_SETS)

align_cor <- function(A, B) {
  # B's columns aligned to A by Hungarian on Pearson r (the estimator's own rule)
  perm <- DeSurv:::.dc_align(A, B)
  B <- B[, perm, drop = FALSE]
  data.frame(factor = seq_len(ncol(A)), matched_col = perm,
             pearson = vapply(seq_len(ncol(A)), function(j) cor(A[, j], B[, j]), 0),
             spearman = vapply(seq_len(ncol(A)), function(j) cor(A[, j], B[, j], method = "spearman"), 0))
}
ref_pairs <- if (length(ref) >= 2L) {
  cmb <- utils::combn(names(ref), 2L)
  do.call(rbind, lapply(seq_len(ncol(cmb)), function(i) {
    a <- cmb[1, i]; b <- cmb[2, i]
    cbind(pair = paste(a, b, sep = "_vs_"), align_cor(ref[[a]]$final$W, ref[[b]]$final$W))
  }))
} else NULL

# ---- Descriptive cross-R agreement of the consensus endpoints -----------------------------------------
ok_cons <- names(cons_fits)[vapply(cons_fits, inherits, TRUE, "desurv_consensus_fit")]
cons_vs_top <- if (length(ok_cons) >= 2L) {
  top <- ok_cons[length(ok_cons)]
  do.call(rbind, lapply(setdiff(ok_cons, top), function(nm)
    cbind(pair = paste(nm, top, sep = "_vs_"), align_cor(cons_fits[[top]]$W, cons_fits[[nm]]$W))))
} else NULL

# ---- Save everything ---------------------------------------------------------------------------------
results <- list(
  gate = gate,
  consensus_fits = cons_fits,
  consensus_vs_largest_R = cons_vs_top,
  reference = ref,
  reference_pairs = ref_pairs,
  bank = list(seeds = BANK_SEEDS, converged = bank_conv, iters = bank_iters,
              objective = vapply(bank, function(f) as.numeric(f$loss$loss)[1], 0),
              cindex = vapply(bank, function(f) as.numeric(f$cindex)[1], 0),
              run_files = list.files(RUN_DIR, pattern = "^bank_")),
  settings = list(dry = DRY, ncores = NCORES, hyper = HP, ntop = NTOP,
                  tol = FIT_TOL, tol_init = FIT_TOL_INIT, maxit = FIT_MAXIT, imaxit = FIT_IMAXIT,
                  R_describe = R_DESCRIBE, R_gate = R_GATE, n_splits = N_SPLITS, split_seed = SPLIT_SEED,
                  ref_sets = REF_SETS, input_md5 = input_md5,
                  desurv_version = as.character(packageVersion("DeSurv")),
                  desurv_path = find.package("DeSurv"),
                  started = t_start, finished = Sys.time()),
  sessionInfo = sessionInfo()
)
saveRDS(results, file.path(OUT_DIR, "gate3a_results.rds"))

# ---- Text summary ------------------------------------------------------------------------------------
fmt <- function(x, d = 3) paste(formatC(x, format = "f", digits = d), collapse = " ")
L <- c(sprintf("Track 3a convergence probe%s. DeSurv %s. %s", if (DRY) " (DRY RUN)" else "",
               packageVersion("DeSurv"), format(Sys.time())),
       sprintf("Config: k=%d alpha=%.6g lambda=%.6g nu=%.6g lambdaW=%g lambdaH=%g ntop=%d; tol=%g tol_init=%g maxit=%d",
               HP$k, HP$alpha, HP$lambda, HP$nu, HP$lambdaW, HP$lambdaH, NTOP, FIT_TOL, FIT_TOL_INIT, FIT_MAXIT),
       sprintf("Bank: %d runs (seeds %d-%d), %d converged, median iters %d",
               length(bank), BANK_SEEDS[1], BANK_SEEDS[length(BANK_SEEDS)], sum(bank_conv), as.integer(median(bank_iters))),
       "",
       sprintf("GATE (half-vs-half: median >= %.2f and p05 >= %.2f per factor; no competing mode). converged_R = %s",
               gate$settings$median_r, gate$settings$p05_r, format(gate$converged_R)))
for (i in seq_len(nrow(gate$summary))) {
  s <- gate$summary[i, ]
  pf <- gate$per_factor[gate$per_factor$R == s$R, , drop = FALSE]
  L <- c(L, sprintf("  R=%d: %s (hvh_pass=%s, mode_pass=%s, errors=%d, half-level mode flags=%s, median |dC_train| halves=%s)",
                    s$R, if (isTRUE(s$pass)) "PASS" else "FAIL", s$hvh_pass, s$mode_pass, s$n_errors,
                    format(s$half_mode_fails), fmt(s$median_abs_dC_train, 4)))
  if (nrow(pf)) {
    L <- c(L, sprintf("        per-factor median r: %s", fmt(pf$median_r, 4)),
              sprintf("        per-factor p05 r:    %s", fmt(pf$p05_r, 4)))
    sp <- gate$splits[gate$splits$R == s$R, , drop = FALSE]
    L <- c(L, sprintf("        beta half1-half2, median |diff| per factor: %s",
                      fmt(vapply(pf$factor, function(j) median(abs(sp$beta_half1[sp$factor == j] -
                                                                   sp$beta_half2[sp$factor == j])), 0), 4)))
  }
}
if (length(gate$errors)) L <- c(L, "  Gate errors:", paste0("    ", gate$errors))
L <- c(L, "", "CONSENSUS FIT per R (descriptive; R=25 is not gated, see README)")
for (nm in names(cons_fits)) {
  cf <- cons_fits[[nm]]
  if (!inherits(cf, "desurv_consensus_fit")) { L <- c(L, sprintf("  %s: ERROR %s", nm, cf)); next }
  cs <- cf$consensus
  occ <- cs$occupancy_summary
  L <- c(L, sprintf("  %s: mode_pass=%s; beta=%s; training C=%.4f", nm, cs$mode_pass, fmt(cf$beta), cf$cindex),
            sprintf("      tightness=%s; tightness_retained=%s; n_members=%s; n_outliers=%s",
                    fmt(cs$tightness), fmt(cs$tightness_retained),
                    paste(cs$n_members, collapse = " "), paste(cs$n_outliers, collapse = " ")),
            sprintf("      occupancy (runs with 0 / 1 / >1 factors in family): %s",
                    paste(sprintf("f%d %d/%d/%d", seq_len(ncol(occ)), occ["n0", ], occ["n1", ], occ["n_gt1", ]),
                          collapse = "; ")),
            sprintf("      competing-mode medoid r=%s, sub-cluster sizes=%s",
                    fmt(cs$competing_mode$medoid_r), paste(sprintf("%d/%d", cs$competing_mode$size_1,
                                                                   cs$competing_mode$size_2), collapse = " ")))
}
if (!is.null(cons_vs_top)) {
  L <- c(L, "  Consensus W vs largest R (aligned): ")
  for (p in unique(cons_vs_top$pair)) {
    z <- cons_vs_top[cons_vs_top$pair == p, ]
    L <- c(L, sprintf("    %s: pearson %s | spearman %s", p, fmt(z$pearson, 4), fmt(z$spearman, 4)))
  }
}
L <- c(L, "", "CURRENT-ESTIMATOR REFERENCE (code/04 recipe: seed fits -> desurv_consensus_seed -> one final fit)")
for (nm in names(ref)) {
  L <- c(L, sprintf("  %s (seeds %d-%d): %d/%d seed fits ok; final training C=%.4f; beta=%s",
                    nm, ref[[nm]]$seeds[1], ref[[nm]]$seeds[length(ref[[nm]]$seeds)], ref[[nm]]$n_ok,
                    length(ref[[nm]]$seeds), ref[[nm]]$final$cindex, fmt(ref[[nm]]$final$beta)))
}
if (!is.null(ref_pairs)) {
  for (p in unique(ref_pairs$pair)) {
    z <- ref_pairs[ref_pairs$pair == p, ]
    L <- c(L, sprintf("  %s: pearson %s | spearman %s", p, fmt(z$pearson, 4), fmt(z$spearman, 4)))
  }
} else L <- c(L, "  (one reference set only: no pairwise correlations)")
L <- c(L, "", sprintf("Elapsed: %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
writeLines(L, file.path(OUT_DIR, "gate3a_summary.txt"))
cat(L, sep = "\n")
logf("Done. Outputs in %s", OUT_DIR)
