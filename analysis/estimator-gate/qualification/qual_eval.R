# Qualification evaluator, one fail-fast stage at a time (spec section 9; prespec sections 3.4, 7 and 9.1, ruling 15).
# Usage: QUAL_MODE=qual|plumbing Rscript qual_eval.R <stage 1-4>
#   stage 1: verifier controls (precondition), realized S3 ratio on grid replicates 1-5 (each in [0.18, 0.37]), and the
#            180 grid and reference fits: Q1, Q3, Q4 on them; Q5 at the centre (n = 175) and reference (n = 218) cells.
#   stage 2: S3 tuning (10): realized-ratio median check; Q1 (final fits, inner-fold fits at the selected configuration,
#            runs without a feasible configuration), Q0-class F increases in other BO evaluations, Q3, Q4, Q5 tuned
#            cell, Q6 S3 part, Q7 (R = 10, t_{0.95, 9}).
#   stage 3: S1 tuning (3): Q1, Q0-class, Q3, Q4, Q5 tuned cell, Q6 S1 part.
#   stage 4: S2 tuning (3): Q1, Q0-class, Q3, Q4, Q5 tuned cell (Q6 S2 reported).
# Reads only <OUT_DIR> (job outputs, the controls file, and earlier stages' result files). Refuses to run unless every
# output of stages 1..k in the manifest exists and, for k > 1, the stage k-1 result file exists: it never evaluates a
# subset. Writes <OUT_DIR>/QUALIFICATION_RESULT_stage<k>.md and .rds with one PASS/FAIL line per binding criterion of
# the stage, the stage line, and the cumulative overall line (PASS only if every stage so far passed; qualification
# passes only if stage 4 reports a cumulative PASS).
# Quantities that need the training X (Q3 window statistics, Q4 power iterations, test-set predictors) are computed
# in the jobs from the stored iterates (qual_common.R sections 6-7) and stored there; this script applies the rules.
QUAL <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(QUAL, "qual_common.R")); source(file.path(QUAL, "manifest.R")); source(file.path(QUAL, "qual_q7.R"))
STAGE <- as.integer(commandArgs(TRUE)[1])
if (!isTRUE(STAGE %in% 1:4)) stop("usage: qual_eval.R <stage 1-4>", call. = FALSE)
STAGE_SCEN <- c(`2` = "S3", `3` = "S1", `4` = "S2")
man <- qual_manifest()
need <- stage_files(man, seq_len(STAGE))
miss <- need[!file.exists(need)]
if (length(miss)) stop(sprintf("REFUSING TO EVALUATE stage %d: %d of %d expected outputs of stages 1-%d missing (first: %s)",
                               STAGE, length(miss), length(need), STAGE, miss[1]), call. = FALSE)
prev <- NULL
if (STAGE > 1L) {
  if (!file.exists(eval_file(STAGE - 1L))) stop(sprintf("REFUSING TO EVALUATE stage %d: %s missing", STAGE, eval_file(STAGE - 1L)), call. = FALSE)
  prev <- readRDS(eval_file(STAGE - 1L))
}
PLUMB <- QUAL_MODE == "plumbing"
if (STAGE == 1L) { GR <- lapply(man$grid$out, readRDS); TU <- list() } else {
  tj <- man$tune[man$tune$stage == STAGE, ]; GR <- list(); TU <- lapply(tj$out, readRDS)
}
cat(sprintf("[%s] stage %d: read %d grid outputs and %d tuning outputs (manifest for this stage: %d + %d)\n", QUAL_MODE, STAGE,
            length(GR), length(TU), sum(man$grid$stage == STAGE), sum(man$tune$stage == STAGE)))

# ---- consistency of the outputs -----------------------------------------------------------------
allo <- c(GR, TU)
stopifnot(all(vapply(allo, function(o) identical(o$mode, QUAL_MODE), TRUE)))
hk <- unique(lapply(allo, function(o) o$hashes$checked))
if (length(hk) != 1L || !identical(hk[[1]], HASHES$checked)) stop("outputs carry different frozen/implementation hashes", call. = FALSE)
if (!is.null(prev) && !identical(prev$hashes$checked, HASHES$checked)) stop("earlier stage was evaluated under different hashes", call. = FALSE)
rec_hash <- unique(lapply(allo, function(o) o$hashes$recorded))
nst_ok <- all(vapply(GR, function(o) o$nstarts == QUAL_NSTARTS, TRUE)) &&
  all(vapply(TU, function(o) !isTRUE(o$tuned) || o$final$nstarts == QUAL_NSTARTS, TRUE))
if (!nst_ok) stop("an output has the wrong number of starts", call. = FALSE)
if (STAGE == 1L) for (i in seq_along(GR)) stopifnot(identical(GR[[i]]$job$out, man$grid$out[i])) else
  for (i in seq_along(TU)) stopifnot(identical(TU[[i]]$job$out, tj$out[i]), identical(TU[[i]]$job$scenario, STAGE_SCEN[[as.character(STAGE)]]))

# ---- fit units of this stage --------------------------------------------------------------------
lab_g <- function(o) sprintf("grid %s r%d g%d", o$job$scenario, o$job$r, o$job$g)
lab_t <- function(o) sprintf("tuned %s r%d", o$job$scenario, o$job$r)
UNITS <- c(lapply(GR, function(o) list(label = lab_g(o), kind = "grid", scen = o$job$scenario, r = o$job$r, g = o$job$g, fit = o)),
           unlist(lapply(TU, function(o) if (isTRUE(o$tuned)) list(list(label = lab_t(o), kind = "final", scen = o$job$scenario,
                                                                           r = o$job$r, g = NA, fit = o$final))), recursive = FALSE))
INNER <- unlist(lapply(TU, function(o) if (isTRUE(o$tuned)) lapply(1:5, function(f) list(label = sprintf("%s inner f%d", lab_t(o), f),
                                                                                                fit = o$inner_fits_selected[[f]]))), recursive = FALSE)
LINES <- list()
line <- function(id, pass, detail, binding = TRUE) {
  v <- if (is.na(pass)) "NOT EVALUABLE" else if (pass) "PASS" else "FAIL"
  if (!binding) v <- paste(v, "(reported)")
  LINES[[length(LINES) + 1L]] <<- data.frame(id = id, result = v, binding = binding, detail = detail, stringsAsFactors = FALSE)
}

# ---- stage preconditions: verifier controls and the realized S3 ratio ------------------------------
if (STAGE == 1L) {
  cf <- file.path(OUT_DIR, "controls", "controls.rds")
  ctl <- if (file.exists(cf)) readRDS(cf) else NULL
  pz <- function(x) isTRUE(!is.null(x) && "pass" %in% names(x) && x[["pass"]] == 1)
  line("Verifier controls (precondition, prespec section 7)", if (is.null(ctl)) FALSE else isTRUE(ctl$pass_all),
       if (is.null(ctl)) "controls file missing" else sprintf("Q7 %s, Q6 %s, Q3 negative %s, Q3 positive %s",
         if (pz(ctl$Q7)) "PASS" else "FAIL", if (pz(ctl$Q6)) "PASS" else "FAIL/NOT RUN",
         if (pz(ctl$Q3neg)) "PASS" else "FAIL/NOT RUN", if (pz(ctl$Q3pos)) "PASS" else "BLOCKED/FAIL"))
  r1 <- unique(do.call(rbind, lapply(Filter(function(o) o$job$scenario == "S3", GR), function(o) data.frame(r = o$job$r, ratio = o$data$ratio_rank_train))))
  line("H1 (realized S3 SD ratio, grid replicates, each in [0.18, 0.37])", nrow(r1) > 0 && all(r1$ratio >= 0.18 & r1$ratio <= 0.37),
       sprintf("S3 replicates %d, range %.3f-%.3f, median %.3f (median reported at stage 1)%s", nrow(r1), min(r1$ratio), max(r1$ratio),
               median(r1$ratio), if (PLUMB) " [p = 300 generator: calibration does not apply]" else ""))
}
if (STAGE == 2L) {
  s3r <- vapply(TU, function(o) o$data$ratio_rank_train, 0)
  line("H2 (realized S3 SD ratio, tuning replicates: each in [0.18, 0.37], median in [0.22, 0.29])",
       all(s3r >= 0.18 & s3r <= 0.37) && median(s3r) >= 0.22 && median(s3r) <= 0.29,
       sprintf("S3 replicates %d, range %.3f-%.3f, median %.3f%s", length(s3r), min(s3r), max(s3r), median(s3r),
               if (PLUMB) " [p = 300 generator: calibration does not apply]" else ""))
}

# ---- Q1 (+ Q0-class F increases elsewhere) -------------------------------------------------------
q1_one <- function(f) {
  st <- f$starts; need <- ceiling(0.9 * nrow(st)); nc <- sum(st$status == "converged")
  data.frame(conv = nc, need = need, finc = sum(st$ninc > 0, na.rm = TRUE), fail = sum(st$status == "failure"),
             pass = nc >= need && all(st$ninc == 0, na.rm = TRUE) && !any(st$status == "failure"))
}
q1 <- rbind(if (length(UNITS)) do.call(rbind, lapply(UNITS, function(u) cbind(unit = u$label, scope = u$kind, q1_one(u$fit)))),
            if (length(INNER)) do.call(rbind, lapply(INNER, function(u) cbind(unit = u$label, scope = "inner-selected", q1_one(u$fit)))))
if (is.null(q1)) q1 <- data.frame(unit = character(), scope = character(), conv = integer(), need = numeric(), finc = integer(),
                                  fail = integer(), pass = logical())   # e.g. every tuning run of the stage without a feasible configuration
untuned <- vapply(Filter(function(o) !isTRUE(o$tuned), TU), lab_t, "")
line("Q1 (ordinary convergence)", (is.null(q1) || all(q1$pass)) && !length(untuned),
     sprintf("fits %d (grid %d, tuned final %d, inner-fold at selected %d): failing %d; tuning runs without a feasible configuration %d%s",
             if (is.null(q1)) 0L else nrow(q1), sum(q1$scope == "grid"), sum(q1$scope == "final"), sum(q1$scope == "inner-selected"),
             sum(!q1$pass), length(untuned), if (length(untuned)) paste0(" (", paste(untuned, collapse = ", "), ")") else ""))
other_inc <- NULL
if (STAGE > 1L) {
  other_inc <- do.call(rbind, lapply(TU, function(o) do.call(rbind, lapply(seq_along(o$evals), function(e) {
    if (isTRUE(o$tuned) && e == o$selected) return(NULL)
    data.frame(run = lab_t(o), eval = e, finc_starts = sum(vapply(o$evals[[e]]$fits, function(f) sum(f$starts$ninc > 0, na.rm = TRUE), 0)))
  }))))
  n_other <- if (is.null(other_inc)) 0 else sum(other_inc$finc_starts)
  line("Q0-class (F increase in any other qualification fit)", n_other == 0,
       sprintf("non-selected BO evaluations %d (inner-fold fits %d): starts with an F increase %d",
               if (is.null(other_inc)) 0L else nrow(other_inc), 5L * (if (is.null(other_inc)) 0L else nrow(other_inc)), n_other))
}

# ---- Q3 ---------------------------------------------------------------------------------------
q3_rows <- function(units) do.call(rbind, lapply(units, function(u) {
  f <- u$fit; st <- f$starts
  do.call(rbind, lapply(seq_len(nrow(st)), function(s) {
    sc <- st$q3_scope[s]; rec <- NA; car <- NA_character_; tra <- NA_character_
    if (sc == "evaluate") {
      q <- f$q3[[s]]; d <- q3_detect(q$stats, q$theta_rec)
      rec <- d$recurrence; car <- paste(d$carrying, collapse = ","); tra <- paste(d$transitions, collapse = ",")
    }
    data.frame(unit = u$label, start = s, scope = sc, recurrence = rec, carrying = car, transitions = tra, stringsAsFactors = FALSE)
  }))
}))
q3 <- q3_rows(UNITS)
q3_inner <- if (length(INNER)) q3_rows(INNER) else NULL
nev <- sum(q3$scope == "evaluate"); nrec <- sum(q3$recurrence %in% TRUE)
ntrans <- sum(nzchar(q3$transitions[q3$scope == "evaluate"]))
line("Q3 (no score-preserving collapse with compensation)", if (!length(UNITS)) NA else nrec == 0,
     sprintf("starts %d: evaluated %d, negative (converged) %d, not evaluable %d, numerical failure %d (reported separately, never passing); recurrence %d; evaluated starts with an active-set transition %d",
             if (is.null(q3)) 0L else nrow(q3), nev, sum(q3$scope == "negative"), sum(q3$scope == "not_evaluable"),
             sum(q3$scope == "numerical_failure"), nrec, ntrans))
q3_desc <- do.call(rbind, lapply(UNITS, function(u) do.call(rbind, lapply(which(u$fit$starts$q3_scope == "evaluate"), function(s) {
  q <- u$fit$q3[[s]]; d <- q3_detect(q$stats, q$theta_rec); if (!length(d$carrying)) return(NULL)
  f <- q$stats$full
  data.frame(unit = u$label, start = s, col = d$carrying, dlog_sigma = f$dlsig[d$carrying], dlog_g = f$dlg[d$carrying],
             dlog_r = f$dlr[d$carrying], dlog_abs_theta = f$dlth[d$carrying], net_angle_rad = f$angle[d$carrying])
}))))

# ---- Q4 ---------------------------------------------------------------------------------------
q4 <- do.call(rbind, lapply(UNITS, function(u) {
  rt <- u$fit$ret; if (is.null(rt)) return(NULL)
  do.call(rbind, lapply(names(rt$q4), function(t) { x <- rt$q4[[t]]
    data.frame(unit = u$label, point = t, sweep = x$sweep, lmax = x$lmax, LW = x$LW, ratio = x$lmax / x$LW, iters = x$iters,
               conv = x$conv, lmin = x$lmin, cond = x$lmax / x$lmin, min_sigma = x$min_sigma, pass = x$pass) }))
}))
n_q4_units <- sum(vapply(UNITS, function(u) !is.null(u$fit$ret), TRUE))
line("Q4 (W-block curvature within the analytic bound)", if (!length(UNITS)) NA else !is.null(q4) && all(q4$pass) && n_q4_units == length(UNITS),
     sprintf("returned starts %d of %d units, estimates %d (sweeps 1, 10, 100, end): above L_W(1 + 1e-8) %d; power iterations not converged in 500 %d",
             n_q4_units, length(UNITS), if (is.null(q4)) 0L else nrow(q4), if (is.null(q4)) 0L else sum(!q4$pass),
             if (is.null(q4)) 0L else sum(!q4$conv)))

# ---- Q5 ---------------------------------------------------------------------------------------
q5_fit <- function(f) {
  st <- f$starts; b <- f$best_start; if (is.na(b) || is.null(f$ret)) return(c(m = NA, n_other = 0))
  oth <- which(st$status == "converged" & seq_len(nrow(st)) != b)
  c(m = if (length(oth)) median(st$testC[oth] - f$ret$testC) else NA_real_, n_other = length(oth))
}
cells <- if (STAGE == 1L) expand.grid(cell = c("centre n175", "reference n218"), scenario = c("S1", "S2", "S3"), stringsAsFactors = FALSE) else
  data.frame(cell = "tuned", scenario = STAGE_SCEN[[as.character(STAGE)]], stringsAsFactors = FALSE)
q5 <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  s <- cells$scenario[i]; cl <- cells$cell[i]
  us <- Filter(function(u) u$scen == s && ((cl == "centre n175" && isTRUE(u$g == 5)) || (cl == "reference n218" && isTRUE(u$g == 12)) ||
                                           (cl == "tuned" && u$kind == "final")), UNITS)
  expected <- if (cl == "tuned") sum(man$tune$scenario == s) else sum(man$grid$scenario == s & man$grid$g == (if (cl == "centre n175") 5 else 12))
  ms <- do.call(rbind, lapply(us, function(u) q5_fit(u$fit)))
  m <- if (is.null(ms)) numeric() else ms[, "m"]
  mm <- if (length(m) && any(!is.na(m))) mean(m, na.rm = TRUE) else NA_real_
  data.frame(scenario = s, cell = cl, fits = length(us), expected = expected, m_na = sum(is.na(m)), mean_m = mm,
             pass = if (is.na(mm)) NA else mm <= 0.01)
}))
q5_pass <- if (any(is.na(q5$pass))) NA else all(q5$pass) && all(q5$fits == q5$expected)
line("Q5 (multi-start stability, one-sided, cell mean of per-fit median <= 0.01)", q5_pass,
     paste(sprintf("%s %s: %s (%s, fits %d/%d, fits without a converged non-returned start %d)", q5$scenario, q5$cell,
                   ifelse(is.na(q5$pass), "NOT EVALUABLE", ifelse(q5$pass, "PASS", "FAIL")),
                   ifelse(is.na(q5$mean_m), "mean m NA", sprintf("mean m %.4f", q5$mean_m)), q5$fits, q5$expected, q5$m_na), collapse = "; "))
q5_extra <- do.call(rbind, lapply(UNITS, function(u) {
  f <- u$fit; st <- f$starts; b <- f$best_start; if (is.na(b) || is.null(f$ret)) return(NULL)
  ok <- which(st$status == "converged" & seq_len(nrow(st)) != b)
  Wb <- f$ret$W; carry <- which(abs(f$ret$theta) > 1e-3)
  lr <- vapply(ok, function(s) {
    W <- f$end[[s]]$W; if (!length(carry)) return(NA_real_)
    cm <- cor(Wb, W); a <- as.integer(clue::solve_LSAP(cm - min(cm), maximum = TRUE)); median(cm[cbind(carry, a[carry])])
  }, 0)
  Fc <- st$F[st$status == "converged"]; Fmin <- if (length(Fc)) min(Fc) else NA_real_
  data.frame(unit = u$label, conv_nonreturned = length(ok), median_loading_r = if (length(lr)) median(lr, na.rm = TRUE) else NA,
             share_within_1e4_relF = if (length(Fc)) mean(abs(Fc - Fmin) <= 1e-4 * abs(Fmin)) else NA_real_)
}))

# ---- Q6 (per scenario at its stage: S3 at stage 2, S1 at stage 3, S2 reported at stage 4) ------------------
TT <- c("10", "50", "100", "500")
q6 <- NULL
if (STAGE > 1L) {
  s <- STAGE_SCEN[[as.character(STAGE)]]
  q6 <- do.call(rbind, lapply(Filter(function(u) u$kind == "final", UNITS), function(u) {
    rt <- u$fit$ret
    data.frame(unit = u$label, returned_converged = identical(rt$status, "converged"), t(rt$testC - rt$ck_testC[TT]))
  }))
  R6 <- sum(man$tune$scenario == s)
  ev6 <- !is.null(q6) && nrow(q6) == R6
  mt <- if (ev6) colMeans(as.matrix(q6[, paste0("X", TT)])) else rep(NA_real_, 4)
  line(sprintf("Q6 %s part (convergence vs checkpoints, mean >= -0.01)", s), if (ev6) all(mt >= -0.01) else NA,
       sprintf("%s (R = %d of %d, returned start converged in %d): mean C_conv - C_t at t = 10/50/100/500: %s", s,
               if (is.null(q6)) 0L else nrow(q6), R6, if (is.null(q6)) 0L else sum(q6$returned_converged),
               paste(sprintf("%.4f", mt), collapse = "/")), binding = s != "S2")
}

# ---- Q7 (stage 2) -------------------------------------------------------------------------------
q7 <- NULL
if (STAGE == 2L) {
  q7 <- do.call(rbind, lapply(Filter(function(u) u$kind == "final", UNITS), function(u) {
    f <- u$fit; o <- TU[[which(vapply(TU, lab_t, "") == u$label)]]
    d <- q7_dab(f$ret$lp_test, f$truth$scores_test[, 1], f$truth$scores_test[, 2])
    data.frame(unit = u$label, D_AB = d[["joint"]], D_AB_marginal = d[["marginal"]],
               cor_AB_generator = o$data$cor_AB_test[["generator"]], cor_AB_masked_rank = o$data$cor_AB_test[["masked_rank"]])
  }))
  R7 <- sum(man$tune$scenario == "S3")
  if (!is.null(q7) && nrow(q7) == R7 && R7 >= 2 && all(is.finite(q7$D_AB))) {
    lb <- mean(q7$D_AB) - qt(0.95, R7 - 1) * sd(q7$D_AB) / sqrt(R7)
    line("Q7 (effect on the predictive direction, S3)", lb > 0,
         sprintf("R = %d: mean D_AB %.4f, SD %.4f, lower bound mean - t(0.95, %d) SD / sqrt(R) = %.4f (no top-up, ruling 15)",
                 R7, mean(q7$D_AB), sd(q7$D_AB), R7 - 1, lb))
  } else line("Q7 (effect on the predictive direction, S3)", NA,
              sprintf("S3 tuned final fits %d of %d; replicates with a constant predictor (D_AB undefined) %d",
                      if (is.null(q7)) 0L else nrow(q7), R7, if (is.null(q7)) 0L else sum(!is.finite(q7$D_AB))))
}

# ---- Q2 and V1-V3 (reported) --------------------------------------------------------------------
q2 <- do.call(rbind, lapply(UNITS, function(u) {
  f <- u$fit; st <- f$starts
  do.call(rbind, lapply(which(st$status == "converged"), function(s) {
    e <- f$end[[s]]; cs <- f$colstats[[s]]; carry <- which(abs(e$theta) > 1e-3); if (!length(carry)) return(NULL)
    s50 <- if (!is.null(cs) && any(cs[, "sweep"] == 50)) cs[cs[, "sweep"] == 50, paste0("sigma", carry)] else e$sigma[carry]
    data.frame(unit = u$label, scen = u$scen, start = s, col = carry, dlog_sigma_end_vs_50 = log(e$sigma[carry] / s50))
  }))
}))
v <- do.call(rbind, lapply(UNITS, function(u) {
  f <- u$fit; rt <- f$ret; if (is.null(rt)) return(NULL)
  size <- if (u$kind == "final") "tuned n218" else if (isTRUE(u$g == 12)) "reference n218" else sprintf("grid g%d n175", u$g)
  cc <- function(m) if (sd(rt$lp_test) > 0) cor(rt$lp_test, f$truth$linpred_test, method = m) else NA_real_
  data.frame(unit = u$label, scen = u$scen, size = size, testC = rt$testC, oracleC = f$truth$oracleC, oracle_gap = f$truth$oracleC - rt$testC,
             pearson_lp = cc("pearson"), spearman_lp = cc("spearman"), mean_cos_matched = mean(rt$v3$cos_matched),
             maxtheta_on_target = rt$v3$maxtheta_matched_to_target, mean_top_overlap = mean(rt$v3$top_overlap))
}))

# ---- write -------------------------------------------------------------------------------------
L <- do.call(rbind, LINES)
stage_pass <- all(L$result[L$binding] == "PASS")
cum_pass <- stage_pass && (STAGE == 1L || isTRUE(prev$cumulative_pass))
md_table <- function(d, digits = 4) {
  if (is.null(d) || !nrow(d)) return("(none)\n")
  d <- as.data.frame(d); for (j in seq_along(d)) if (is.numeric(d[[j]])) d[[j]] <- formatC(d[[j]], digits = digits, format = "g")
  paste0("| ", paste(names(d), collapse = " | "), " |\n|", paste(rep("---", ncol(d)), collapse = "|"), "|\n",
         paste(apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), collapse = "\n"), "\n")
}
qs <- function(x) if (is.null(x) || !length(x)) "none" else paste(sprintf("%.3g", quantile(x, c(0, 0.5, 1), na.rm = TRUE)), collapse = " / ")
stage_name <- c("1 (controls, grid and reference fits)", "2 (S3 tuning)", "3 (S1 tuning)", "4 (S2 tuning)")[STAGE]
hdr <- if (PLUMB) "# PLUMBING TEST, NOT QUALIFICATION\n\nThis run checks that the pipes connect and that the evaluator counts correctly (p = 300 generator, reduced manifest, 3 starts per fit, BO n_init 2 n_iter 1, seeds 49,300,001+). **Its PASS/FAIL values mean nothing and must not be read as qualification outcomes.**\n\n" else
  "# Qualification result (spec 2026-10-05 revision 4, section 9; prespec 2026-10-05 with ruling 15)\n\n"
md <- c(hdr, sprintf("## Stage %s\n\n", stage_name),
  sprintf("Generated %s. Mode %s. Outputs read for this stage: grid %d, tuning %d (all expected by the manifest; stages 1-%d complete). Frozen and implementation hashes (identical in every output): %s. Recorded-only hash sets: %d. git HEAD at the first job: %s. DeSurv %s (%s).\n\n",
          format(Sys.time()), QUAL_MODE, length(GR), length(TU), STAGE, paste(sprintf("%s %s", basename(names(HASHES$checked)), HASHES$checked), collapse = ", "),
          length(rec_hash), allo[[1]]$hashes$git_head, HASHES$desurv[["version"]], HASHES$desurv[["sha"]]),
  "### Binding criteria of this stage\n\n",
  paste(sprintf("- **%s: %s.** %s", L$id, L$result, L$detail), collapse = "\n"),
  sprintf("\n\n**STAGE %d: %s.**\n\n**CUMULATIVE OVERALL (stages 1-%d): %s**%s%s\n\n", STAGE, if (stage_pass) "PASS" else "FAIL", STAGE,
          if (cum_pass) "PASS" else "FAIL", if (STAGE == 4L) " (qualification result)" else " (qualification passes only if stage 4 reports a cumulative PASS)",
          if (PLUMB) " [plumbing: meaningless]" else ""),
  if (STAGE == 1L) "Q0 (mathematical checks, smoke tests S1-S7) is rerun at the start of stage 1 and is not evaluated by this script.\n\n" else "",
  "### Q1 fits failing\n\n", md_table(q1[!q1$pass, , drop = FALSE]),
  "\n### Q3 descriptive (evaluated starts, carrying columns; D-Q4 reference: log sigma loss 6 to 9, net angle at most 0.0023 rad)\n\n", md_table(q3_desc),
  sprintf("\nInner-fold fits at the selected configurations (reported, not in the binding scope): evaluated %d, recurrence %d.\n",
          if (is.null(q3_inner)) 0L else sum(q3_inner$scope == "evaluate"), if (is.null(q3_inner)) 0L else sum(q3_inner$recurrence %in% TRUE)),
  "\n### Q4 realized curvature (reported)\n\n",
  sprintf("lambda_max / L_W min / median / max: %s. Condition number (lambda_max / lambda_min on the free face) min / median / max: %s. Slope of log lambda_max on log min_j sigma_j: %s.\n\n",
          qs(q4$ratio), qs(q4$cond[is.finite(q4$cond)]),
          if (!is.null(q4) && nrow(q4) > 2 && all(q4$lmax > 0)) sprintf("%.3f", coef(lm(log(lmax) ~ log(min_sigma), data = q4))[2]) else "n/a"),
  md_table(q4),
  "\n### Q5 cells and reported extras\n\n", md_table(q5), "\n", md_table(q5_extra),
  if (!is.null(q6)) c("\n### Q6 per fit (C_conv - C_t)\n\n", md_table(q6)) else "",
  if (!is.null(q7)) c("\n### Q7 per replicate (S3)\n\n", md_table(q7)) else "",
  "\n### Q2 (reported): log(sigma_end / sigma_50), carrying columns of converged starts\n\n",
  md_table(if (is.null(q2)) NULL else do.call(rbind, lapply(split(q2, q2$scen), function(d) data.frame(scenario = d$scen[1], columns = nrow(d),
           min = min(d$dlog_sigma_end_vs_50), median = median(d$dlog_sigma_end_vs_50), max = max(d$dlog_sigma_end_vs_50))))),
  "\n### V1-V3 (reported)\n\n", md_table(v))
writeLines(md, eval_file(STAGE, "md"))
saveRDS(list(stage = STAGE, lines = L, stage_pass = stage_pass, cumulative_pass = cum_pass, mode = QUAL_MODE, q1 = q1, q0_other = other_inc,
             q3 = q3, q3_inner = q3_inner, q3_desc = q3_desc, q4 = q4, q5 = q5, q5_extra = q5_extra, q6 = q6, q7 = q7,
             q2 = q2, v = v, hashes = HASHES, recorded_hashes = rec_hash), eval_file(STAGE))
cat(sprintf("wrote 2 files: %s and .md\n", eval_file(STAGE)))
for (i in seq_len(nrow(L))) cat(sprintf("%s: %s\n", L$id[i], L$result[i]))
cat(sprintf("STAGE %d: %s\nCUMULATIVE OVERALL (stages 1-%d): %s%s\n", STAGE, if (stage_pass) "PASS" else "FAIL", STAGE,
            if (cum_pass) "PASS" else "FAIL", if (PLUMB) " (PLUMBING: meaningless)" else ""))
