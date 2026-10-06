# Shared setup for the simulation-qualification harness (brief docs/plans/2026-10-05-qualification-harness-brief.md).
# Authority: FROZEN spec docs/plans/2026-10-05-estimator-redesign-spec.md (revision 4) section 9 and FROZEN
# docs/plans/2026-10-05-qualification-prespec.md. Estimator: analysis/estimator-gate/redesign/theta_lib.R (sourced,
# never modified here). Data: sim_qual.R (frozen). Simulated data only; no real data is read anywhere.
#
# QUAL_MODE (environment, required, no default):
#   "qual"     : the qualification itself. p = 3,000, outputs under out/. Seeds of prespec section 4.
#   "plumbing" : pipe test only. p = 300 generator (see plumbing_install_generator), reduced manifest, seeds
#                49,300,001+, outputs under out_plumbing/. Its PASS/FAIL values mean nothing.
# Order: (1) frozen-file hashes are checked BEFORE anything is sourced or generated; (2) the frozen files and the
# estimator are sourced; (3) seeds, data, fit, scoring and recording helpers are defined.
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
if (!exists("QUAL")) stop("qual_common.R: set QUAL (path of analysis/estimator-gate/qualification) first", call. = FALSE)
QUAL <- normalizePath(QUAL)
ROOT <- normalizePath(file.path(QUAL, "..", "..", ".."))
RED <- file.path(ROOT, "analysis", "estimator-gate", "redesign")
QUAL_MODE <- Sys.getenv("QUAL_MODE", "")
if (!QUAL_MODE %in% c("qual", "plumbing")) stop("set QUAL_MODE to 'qual' or 'plumbing' (no default)", call. = FALSE)
OUT_DIR <- file.path(QUAL, if (QUAL_MODE == "qual") "out" else "out_plumbing")
QUAL_P <- if (QUAL_MODE == "qual") 3000L else 300L

# ---- 1. frozen-file hashes (prespec section 3.3) and implementation hashes ----------------------
# Git blob hashes, first 12 characters. FROZEN_HASHES are the prespec section 3.3 list. IMPL_HASHES are theta_lib.R
# (6f0c862 + the Q3 recording point of ruling 14 + the recording-only Q4 snapshot addition of this harness, sweeps
# 1, 10, 100) and kkt_oracle_theta.R as they stand when this harness was built; every job records the values it
# found and stops on a mismatch. If theta_lib.R is committed unchanged its blob hash stays the same.
FROZEN_HASHES <- c(
  "analysis/estimator-gate/qualification/sim_qual.R" = "2379bb00de3f",
  "analysis/estimator-gate/qualification/s3_calibrate.R" = "0f2017a0ac58",
  "analysis/estimator-gate/qualification/s3_calibration.rds" = "d8e7d1f9f534",
  "R/simulation_functions/simulate_data.R" = "301920a23600",
  "R/simulation_functions/simulate_W.R" = "f21b1af5d47a",
  "R/simulation_functions/simulate_H.R" = "77a66ef5fb8e",
  "R/simulation_functions/simulate_X.R" = "d0e55aede746",
  "R/simulation_functions/simulate_survival.R" = "ac52c17b1e01",
  "R/simulation_functions/scenario_defaults.R" = "dca30349342d")
IMPL_HASHES <- c(
  "analysis/estimator-gate/redesign/theta_lib.R" = "b3dbc0cf4edc",
  "analysis/estimator-gate/redesign/kkt_oracle_theta.R" = "782e86656faf")
RECORDED_ONLY <- c("docs/plans/2026-10-05-estimator-redesign-spec.md", "docs/plans/2026-10-05-qualification-prespec.md",
                   "analysis/estimator-gate/qualification/bo_engine_copy.R", "analysis/estimator-gate/qualification/qual_common.R")
blob12 <- function(rel) substr(system2("git", c("-C", ROOT, "hash-object", rel), stdout = TRUE), 1, 12)
qual_check_hashes <- function() {
  want <- c(FROZEN_HASHES, IMPL_HASHES)
  got <- vapply(names(want), blob12, "")
  bad <- names(want)[got != want]
  if (length(bad)) stop(sprintf("frozen-file hash mismatch (stopping before any data): %s",
                                paste(sprintf("%s found %s expected %s", bad, got[bad], want[bad]), collapse = "; ")), call. = FALSE)
  head <- tryCatch(system2("git", c("-C", ROOT, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_)
  list(checked = got, recorded = vapply(RECORDED_ONLY, blob12, ""), git_head = head,
       desurv = c(version = as.character(utils::packageVersion("DeSurv")),
                  sha = utils::packageDescription("DeSurv")$RemoteSha %||% NA_character_))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
HASHES <- qual_check_hashes()

# ---- 2. frozen files and estimator ---------------------------------------------------------------
source(file.path(QUAL, "sim_qual.R"))
source(file.path(RED, "kkt_oracle_theta.R")); source(file.path(RED, "theta_lib.R"))
suppressMessages({ library(survival) })
CAL <- readRDS(file.path(QUAL, "s3_calibration.rds"))
stopifnot(CAL$S3$cA == 0.1911, CAL$S3b$cA == 0.474)                    # prespec section 3.3 constants
SCEN_IDX <- c(S1 = 1L, S2 = 2L, S3 = 3L, S3b = 4L)
SCEN_CVEC <- list(S1 = c(1, 1, 1), S2 = c(1, 1, 1), S3 = c(CAL$S3$cA, 1, 1), S3b = c(CAL$S3b$cA, 1, 1))
QUAL_NSTARTS <- if (QUAL_MODE == "qual") 30L else 3L
QUAL_MAXIT <- 5000L; QUAL_TOL <- 1e-5; LAMBDA_H <- 1e-6

# Plumbing only: sim_qual() calls gen_scaled() with its default p = 3,000 and fixed marker / background sizes
# (150 per factor, 500 background), which do not fit in 300 genes. The plumbing session rebinds gen_scaled with the
# candidate/sim.R scaling (markers max(5, round(150 p / 3000)), background round(500 p / 3000)). In memory only;
# sim_qual.R is not edited, and qual mode never calls this.
plumbing_install_generator <- function(p = 300L) {
  stopifnot(QUAL_MODE == "plumbing")
  assign("gen_scaled", function(seed, beta, cvec, N, p0 = p) {
    calls <- 0L
    assign("simulate_H", function(...) { calls <<- calls + 1L; sweep(.simulate_H_orig(...), 2, cvec, "*") }, envir = globalenv())
    on.exit(assign("simulate_H", .simulate_H_orig, envir = globalenv()))
    s <- suppressWarnings(simulate_desurv_scenario(scenario = "R0", G = p0, N = N, K = 3L, beta = beta,
           markers_per_factor = max(5L, round(150 * p0 / 3000)), B_size = round(500 * p0 / 3000), seed = seed))
    stopifnot(calls == 1L)
    s
  }, envir = globalenv())
}
if (QUAL_MODE == "plumbing") plumbing_install_generator(QUAL_P)

# ---- 3. seeds (prespec section 4: one function per row of the seed table) ------------------------
# sc = 1..4 (S1, S2, S3, S3b), r = replicate, g = fit index 1..12, f = inner fold 1..5. Start s of a fit with
# seed S uses set.seed(S + s - 1) inside palm_draw. Plumbing seeds: 49,300,001 + 100,000 i + 10,000 sc + 1,000 r
# + offset (i = row 0..4), inside 49,300,001-49,999,999 and outside every qualification namespace.
.pl <- function(i, sc, r, off = 0) 49300001 + 100000 * i + 10000 * sc + 1000 * r + off
seed_dataset <- function(sc, r) if (QUAL_MODE == "qual") 40000000 + 100000 * sc + r else .pl(0, sc, r)
seed_grid <- function(sc, r, g) if (QUAL_MODE == "qual") 41000000 + 100000 * sc + 1000 * r + 50 * (g - 1) else .pl(1, sc, r, 50 * (g - 1))
seed_tune <- function(sc, r) if (QUAL_MODE == "qual") 42000000 + 100000 * sc + r else .pl(2, sc, r)   # folds and BO
seed_inner <- function(sc, r, f) if (QUAL_MODE == "qual") 43000000 + 100000 * sc + 1000 * r + 100 * (f - 1) else .pl(3, sc, r, 100 * (f - 1))
seed_final <- function(sc, r) if (QUAL_MODE == "qual") 44000000 + 100000 * sc + 1000 * r else .pl(4, sc, r)

# Disjointness over the full job list: every seed interval [S, S + nstarts - 1] (single seeds for datasets and the
# tuning partition/BO seed) is pairwise disjoint, and none touches a used or excluded range of prespec section 4.
SEED_EXCLUDED <- rbind(c(300101, 300834), c(10001100, 10004529), c(3000, 3429), c(20000000, 20000029),
                       c(11, 23), c(123, 123), c(600, 600), c(700, 700), c(1000, 5000), c(31337, 31337),
                       c(424242, 424243), c(700000, 799999), c(20261005, 20261005),
                       c(49000001, 49000020), c(49000101, 49000200), c(49100001, 49100330), c(49200000, 49299999))
assert_seeds_disjoint <- function(man) {
  iv <- list()
  add <- function(lo, len, what) iv[[length(iv) + 1L]] <<- data.frame(lo = lo, hi = lo + len - 1, what = what)
  for (i in seq_len(nrow(man$grid))) {
    j <- man$grid[i, ]
    add(seed_grid(j$sc, j$r, j$g), QUAL_NSTARTS, sprintf("grid %s r%d g%d", j$scenario, j$r, j$g))
  }
  for (i in seq_len(nrow(man$tune))) {
    j <- man$tune[i, ]
    add(seed_tune(j$sc, j$r), 1, sprintf("tune-partition %s r%d", j$scenario, j$r))
    for (f in 1:5) add(seed_inner(j$sc, j$r, f), QUAL_NSTARTS, sprintf("inner %s r%d f%d", j$scenario, j$r, f))
    add(seed_final(j$sc, j$r), QUAL_NSTARTS, sprintf("final %s r%d", j$scenario, j$r))
  }
  ds <- unique(rbind(man$grid[, c("scenario", "sc", "r")], man$tune[, c("scenario", "sc", "r")]))
  for (i in seq_len(nrow(ds))) add(seed_dataset(ds$sc[i], ds$r[i]), 1, sprintf("dataset %s r%d", ds$scenario[i], ds$r[i]))
  iv <- do.call(rbind, iv); iv <- iv[order(iv$lo), ]
  ov <- which(iv$lo[-1] <= iv$hi[-nrow(iv)])
  if (length(ov)) stop(sprintf("seed overlap: %s vs %s", iv$what[ov[1]], iv$what[ov[1] + 1]), call. = FALSE)
  for (e in seq_len(nrow(SEED_EXCLUDED))) {
    hit <- which(iv$lo <= SEED_EXCLUDED[e, 2] & iv$hi >= SEED_EXCLUDED[e, 1])
    if (length(hit)) stop(sprintf("seed in an excluded range: %s", iv$what[hit[1]]), call. = FALSE)
  }
  if (QUAL_MODE == "plumbing" && any(iv$lo < 49300001 | iv$hi > 49999999)) stop("plumbing seed outside 49,300,001-49,999,999", call. = FALSE)
  if (QUAL_MODE == "qual" && any(iv$lo < 40000000 | iv$hi > 44999999)) stop("qualification seed outside 40,000,000-44,999,999", call. = FALSE)
  cat(sprintf("seed check: %d intervals (%d seeds), pairwise disjoint, no excluded range touched\n",
              nrow(iv), sum(iv$hi - iv$lo + 1)))
  invisible(iv)
}

# ---- 4. data ----------------------------------------------------------------------------------------
# One sim_qual draw per (scenario, replicate). S3/S3b run the wrapper self-check (check = TRUE). S3: the replicate's
# realized ranked masked ratio must lie in [0.18, 0.37] (prespec section 3.4), checked in qual mode only (the
# calibration does not apply to the p = 300 plumbing generator).
qual_data <- function(scenario, r) {
  z <- sim_qual(seed_dataset(SCEN_IDX[[scenario]], r), scenario, cvec = SCEN_CVEC[[scenario]])
  if (QUAL_MODE == "qual") stopifnot(nrow(z$X) == 3000L, ncol(z$X) == 218L, ncol(z$X_test) == 1000L)
  if (QUAL_MODE == "qual" && scenario == "S3" && !(z$ratio_rank_train >= 0.18 && z$ratio_rank_train <= 0.37))
    stop(sprintf("S3 r%d realized SD ratio %.4f outside [0.18, 0.37]: harness defect, stop", r, z$ratio_rank_train), call. = FALSE)
  z
}
# Qual mode, fail-fast order (prespec section 9.1): a stage-k job refuses to run unless (a) for S3, the realized-ratio
# check of its stage (qual_precheck.R <k>) passed, and (b) for k > 1, the stage k-1 evaluation reports a cumulative
# PASS. Plumbing mode records but does not enforce (its values mean nothing).
precheck_file <- function(stage) file.path(OUT_DIR, "precheck", sprintf("ratio_check_stage%d.rds", stage))
eval_file <- function(stage, ext = "rds") file.path(OUT_DIR, sprintf("QUALIFICATION_RESULT_stage%d.%s", stage, ext))
require_precheck <- function(scenario, stage) {
  if (QUAL_MODE != "qual" || scenario != "S3") return(invisible(TRUE))
  f <- precheck_file(stage)
  if (!file.exists(f) || !isTRUE(readRDS(f)$pass)) stop(sprintf("stage %d S3 realized-ratio precheck missing or failed: run qual_precheck.R %d first", stage, stage), call. = FALSE)
  invisible(TRUE)
}
require_controls <- function() {
  if (QUAL_MODE != "qual") return(invisible(TRUE))
  f <- file.path(OUT_DIR, "controls", "controls.rds")
  if (!file.exists(f) || !isTRUE(readRDS(f)$pass_all)) stop("verifier controls missing or not all passed (prespec section 7): no qualification fit may run", call. = FALSE)
  invisible(TRUE)
}
require_stage_pass <- function(stage) {
  if (QUAL_MODE != "qual" || stage <= 1L) return(invisible(TRUE))
  f <- eval_file(stage - 1L)
  if (!file.exists(f) || !isTRUE(readRDS(f)$cumulative_pass)) stop(sprintf("stage %d may not run: %s missing or not a cumulative PASS", stage, f), call. = FALSE)
  invisible(TRUE)
}
data_summary <- function(z) list(seed = z$seed, cvec = z$cvec, ratio_rank_train = z$ratio_rank_train,
  ratio_latent_train = z$ratio_latent_train, ratio_full_rank_train = z$ratio_full_rank_train,
  cor_AB_test = z$cor_AB_test, event_rate = z$event_rate, n_tied = z$n_tied, p = nrow(z$X))

# ---- 5. scoring (gate spec section 5: the comparator's CV scorer, cvwrapr::getCindex on Surv(y, d)) ----------
qual_cindex <- function(lp, y, d) {
  if (sum(d) == 0 || any(!is.finite(lp))) return(NA_real_)
  cvwrapr::getCindex(lp, survival::Surv(y, d))
}
# Model of one iterate (W, theta), through palm_model unchanged: own v~ and own ntop support G.
iterate_model <- function(W, theta, D, hp, ntop)
  palm_model(list(best = list(W = W, theta = theta, sigma = NULL), hp = hp, status = "iterate"), D, ntop)
iterate_testC <- function(W, theta, D, hp, ntop, z) {
  m <- iterate_model(W, theta, D, hp, ntop)
  qual_cindex(palm_score(z$X_test, m), z$y_test, z$d_test)
}

# ---- 6. Q4: power iteration on the W-block HVP restricted to the free coordinates (prespec section 7) ----------
# Start vector: ones on {W > 0}, normalized. Estimate ||H v|| for unit v (from below for PSD H). Stop at relative
# change < 1e-8 or 500 iterations. Pass: estimate <= L_W (1 + 1e-8), L_W = palm_LW at the same (H, theta).
# lambda_min (shifted power iteration on L_W I - H, same rule) is recorded for the reported condition number.
q4_power <- function(W, H, theta, D, hp, tol = 1e-8, maxit = 500L) {
  free <- W > 0
  LW <- palm_LW(H, theta, D, hp)
  hv <- function(v) palm_hvp_W(W, H, theta, D, hp, v) * free
  pw <- function(op) {
    v <- free * 1; nv <- sqrt(sum(v^2)); if (nv == 0) return(list(est = 0, it = 0L, conv = TRUE))
    v <- v / nv; est <- NA_real_
    for (i in seq_len(maxit)) {
      w <- op(v); e <- sqrt(sum(w^2))
      if (!is.finite(e)) return(list(est = NA_real_, it = i, conv = FALSE))
      if (e == 0) return(list(est = 0, it = i, conv = TRUE))
      if (i > 1L && abs(e - est) <= tol * e) return(list(est = e, it = i, conv = TRUE))
      est <- e; v <- w / e
    }
    list(est = est, it = maxit, conv = FALSE)
  }
  mx <- pw(hv)
  mn <- pw(function(v) LW * v * free - hv(v))
  lmin <- if (is.finite(mn$est)) LW - mn$est else NA_real_
  list(lmax = mx$est, iters = mx$it, conv = mx$conv, LW = LW, pass = isTRUE(mx$est <= LW * (1 + 1e-8)),
       lmin = lmin, lmin_iters = mn$it, lmin_conv = mn$conv, n_free = sum(free),
       min_sigma = min(palm_scores(W, D)$sigma))
}

# ---- 7. Q3: window statistics (needs X) and the detector (pure comparisons) -----------------------------
# Window sweeps 3,751-5,000, halves split at 4,375 (spec section 9, rulings 11 and 14). Inputs: W and theta at the
# window points a = 3,750, m = 4,375, e = 5,000 and theta at every recorded sweep of the window. Per column:
# sigma_j = sd(Zc_j), r_j = ||P w_j||, g_j = sigma_j / r_j (log sigma = log r + log g), standardized score
# Zc_j / ||Zc_j||, net angle 2 asin(||u - v|| / 2) (as D-Q4).
q3_point <- function(W, D) {
  sc <- palm_scores(W, D); r <- sqrt(colSums(sweep(W, 2, colMeans(W))^2))
  nz <- sqrt(colSums(sc$Zc^2))
  list(lsig = log(sc$sigma), lr = log(r), lg = log(sc$sigma / r), Zs = sweep(sc$Zc, 2, nz, "/"))
}
q3_window_stats <- function(Wa, Wm, We, tha, thm, the, D) {
  a <- q3_point(Wa, D); m <- q3_point(Wm, D); e <- q3_point(We, D)
  ang <- function(u, v) 2 * asin(pmin(1, sqrt(colSums((u - v)^2)) / 2))
  seg <- function(x, y, tx, ty) list(dlsig = y$lsig - x$lsig, dlr = y$lr - x$lr, dlg = y$lg - x$lg,
                                     dlth = log(abs(ty)) - log(abs(tx)), angle = ang(x$Zs, y$Zs))
  list(full = seg(a, e, tha, the), h1 = seg(a, m, tha, thm), h2 = seg(m, e, thm, the))
}
# theta_rec: matrix (rows = every recorded sweep of the window incl. the three window points, cols = columns j).
q3_detect <- function(st, theta_rec, thr = 1e-3) {
  k <- ncol(theta_rec); carry <- logical(k); trans <- logical(k); rec <- logical(k); cond <- vector("list", k)
  for (j in seq_len(k)) {
    th <- theta_rec[, j]
    carry[j] <- all(abs(th) > thr) && length(unique(sign(th))) == 1L
    trans[j] <- !carry[j] && any(abs(th) > thr)          # crosses the threshold, hits zero or changes sign
    if (!carry[j]) next
    f <- st$full; h1 <- st$h1; h2 <- st$h2
    c1 <- f$dlsig[j] < 0; c2 <- abs(f$dlg[j]) > abs(f$dlr[j]); c3 <- f$angle[j] < abs(f$dlsig[j]); c4 <- f$dlth[j] > 0
    half <- function(h) h$dlsig[j] < 0 && abs(h$dlg[j]) > abs(h$dlr[j]) && h$dlth[j] > 0
    c5 <- half(h1) && half(h2)
    cond[[j]] <- c(i = c1, ii = c2, iii = c3, iv = c4, v = c5)
    rec[j] <- isTRUE(c1 && c2 && c3 && c4 && c5)
  }
  list(recurrence = any(rec), carrying = which(carry), transitions = which(trans), recurrent_cols = which(rec), cond = cond)
}
# Scope (spec section 9 Q3): converged = negative by definition; numerical failure = reported separately;
# finite nonconverged at the 5,000-sweep work limit = evaluated; anything else = not evaluable.
q3_scope <- function(status, iter, maxit = QUAL_MAXIT) {
  if (identical(status, "converged")) "negative" else if (identical(status, "failure")) "numerical_failure" else
    if (identical(status, "nonconverged") && isTRUE(iter == maxit)) "evaluate" else "not_evaluable"
}
q3_theta_rec <- function(colstats, ck) {
  th <- colstats[colstats[, "sweep"] > 3750 & colstats[, "sweep"] <= 5000, paste0("theta", 1:3), drop = FALSE]
  rbind(ck[["3750"]]$theta, th, ck[["4375"]]$theta, ck[["end"]]$theta)
}

# ---- 8. V3: loading and program recovery (reported only) ---------------------------------------------
v3_recovery <- function(W, theta, z, ntop) {
  Wn <- sweep(W, 2, pmax(sqrt(colSums(W^2)), 1e-300), "/")
  cosm <- crossprod(Wn, z$W_true)                                         # fitted j x true m
  asg <- as.integer(clue::solve_LSAP(cosm - min(cosm), maximum = TRUE))  # fitted j -> true asg[j]
  rownames(W) <- rownames(z$W_true)
  tg <- DeSurv::desurv_get_top_genes(W, ntop)$top_genes
  ov <- vapply(seq_len(ncol(W)), function(j) length(intersect(tg[[j]], z$marker_sets[[asg[j]]])) / length(tg[[j]]), 0)
  jmax <- which.max(abs(theta)); tgt <- which.max(abs(z$beta_true))     # A (S3, S3b) or the predictive factor (S1)
  list(cos_matched = cosm[cbind(seq_len(ncol(W)), asg)], match = asg, top_overlap = ov,
       maxtheta_col = jmax, maxtheta_matched_to_target = if (all(z$beta_true == 0)) NA else asg[jmax] == tgt)
}

# ---- 9. one fit: palm_fit + per-start records (prespec section 6) -------------------------------------
# keep = "full" (grid and tuned final fits): per-start status, reason, sweeps, F, residuals, F increases,
# column statistics every 10th sweep, (W, theta) at every checkpoint and the end, test C of converged starts,
# Q3 window statistics of evaluable starts; returned start: all checkpoint iterates with their test C, the Q4 power iteration at sweeps
# 1, 10, 100 and the end, and its test-set predictor for Q7/V1-V3.
# keep = "light" (inner-fold fits of every BO evaluation; brief item 3 keeps only the selected configuration's five
# for Q1): status table, F increases, Q3 window statistics, returned start's end iterate.
qual_fit <- function(D, hp, seed, ncores, keep = c("full", "light"), z = NULL, ntop = NULL) {
  keep <- match.arg(keep)
  hp$lambdaH <- LAMBDA_H
  fit <- palm_fit(D, hp, k = 3L, seed = seed, nstarts = QUAL_NSTARTS, maxit = QUAL_MAXIT, tol = QUAL_TOL,
                  ncores = ncores, keep_starts = TRUE, record = TRUE)
  fs <- fit$fits; st <- fit$starts
  st$maxinc <- vapply(fs, function(f) if (is.null(f$maxinc)) NA_real_ else max(f$maxinc), 0)
  st$q3_scope <- vapply(seq_along(fs), function(s) q3_scope(fs[[s]]$status, fs[[s]]$iter), "")
  q3 <- lapply(seq_along(fs), function(s) {
    f <- fs[[s]]; if (st$q3_scope[s] != "evaluate") return(NULL)
    ck <- f$record$checkpoints
    list(stats = q3_window_stats(ck[["3750"]]$W, ck[["4375"]]$W, ck[["end"]]$W, ck[["3750"]]$theta,
                                 ck[["4375"]]$theta, ck[["end"]]$theta, D),
         theta_rec = q3_theta_rec(f$record$colstats, ck), sweeps = c(ck[["3750"]]$sweep, ck[["4375"]]$sweep, ck[["end"]]$sweep))
  })
  out <- list(status = fit$status, best_start = fit$best_start, starts = st, s0 = fit$s0, lam1 = fit$lam1, hp = fit$hp,
              seed = seed, nstarts = QUAL_NSTARTS, maxit = QUAL_MAXIT, n = D$n, p = D$p, ne = D$ne,
              f_increases = do.call(rbind, lapply(fs, function(f) if (!is.null(f$f_increases)) cbind(start = f$start, f$f_increases))),
              q3 = q3)
  if (keep == "light") {   # returned start's end iterate, for the inner-fold validation C
    if (!is.na(fit$best_start)) out$ret_light <- list(start = fit$best_start, W = fit$best$W, theta = fit$best$theta,
                                                      status = fit$best$status, iter = fit$best$iter)
    return(out)
  }
  conv <- which(st$status == "converged")
  out$starts$testC <- NA_real_
  for (s in conv) out$starts$testC[s] <- iterate_testC(fs[[s]]$W, fs[[s]]$theta, D, hp, ntop, z)
  out$end <- lapply(fs, function(f) if (is.null(f$W)) NULL else list(W = f$W, theta = f$theta, sigma = f$sigma, F = f$F))
  out$colstats <- lapply(fs, function(f) f$record$colstats)
  out$checkpoints <- lapply(fs, function(f) f$record$checkpoints)   # (W, theta) at 10, 50, 100, 500, 3750, 4375, end
  b <- fit$best_start
  if (!is.na(b)) {
    fb <- fs[[b]]; ck <- fb$record$checkpoints
    out$ret <- list(
      start = b, status = fb$status, iter = fb$iter, W = fb$W, H = fb$H, theta = fb$theta, sigma = fb$sigma,
      checkpoints = ck,
      ck_testC = vapply(names(ck), function(t) iterate_testC(ck[[t]]$W, ck[[t]]$theta, D, hp, ntop, z), 0),
      ck_sweep = vapply(ck, function(x) as.numeric(x$sweep), 0),
      q4 = lapply(fb$record$q4, function(x) c(q4_power(x$W, x$H, x$theta, D, hp), list(sweep = x$sweep))),
      testC = iterate_testC(fb$W, fb$theta, D, hp, ntop, z))
    m <- iterate_model(fb$W, fb$theta, D, hp, ntop)
    out$ret$lp_test <- palm_score(z$X_test, m); out$ret$G_size <- length(m$G)
    out$ret$v3 <- v3_recovery(fb$W, fb$theta, z, ntop)
  }
  out$truth <- list(linpred_test = z$linpred_test, scores_test = z$scores_test[, 1:2, drop = FALSE],
                    oracleC = qual_cindex(z$linpred_test, z$y_test, z$d_test))
  out
}
fit_counts <- function(o) sprintf("starts %d: converged %d, nonconverged %d, failure %d, F-increase starts %d",
  nrow(o$starts), sum(o$starts$status == "converged"), sum(o$starts$status == "nonconverged"),
  sum(o$starts$status == "failure"), sum(o$starts$ninc > 0, na.rm = TRUE))
save_job <- function(obj, path) {
  if (QUAL_MODE == "plumbing" && grepl(paste0(.Platform$file.sep, "out", .Platform$file.sep), path, fixed = TRUE))
    stop("plumbing mode must not write into out/", call. = FALSE)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp"); saveRDS(obj, tmp); file.rename(tmp, path)
  cat(sprintf("wrote 1 file: %s\n", path))
}
