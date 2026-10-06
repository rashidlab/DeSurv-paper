# Verifier controls (prespec section 7; run after Q0 and before any qualification fit; a failure stops the run and
# is reported, never answered by adjusting a criterion). Usage: QUAL_MODE=qual|plumbing Rscript qual_controls.R
#  Q7 controls   S3 replicate 1's test set: eta = true A gives D_AB > 0.9, eta = true B gives D_AB < -0.9,
#                eta = A + B gives |D_AB| < 0.1 (exact inputs, no fit).
#  Q6 control    a start converged before sweep 10 has every checkpoint difference exactly 0. Built by restarting
#                palm_fit_one (recording on) from the converged endpoint of one control start (control seeds
#                49,900,001 / 49,900,101, outside every qualification and plumbing namespace), scored through the same
#                iterate_testC used by the jobs. Needs no qualification output, so it runs before any qualification fit.
#  Q3 negative   q3_scope routes the converged control start and its restart to "negative" (the scope filter).
#  Q3 positive   the detector, unchanged, on the failed candidate's stalled T7 starts. Needs W (and b) at the window
#                points of the exact-gated D-Q4 replay; checks whether they are stored and, if not, reports the rerun
#                that would produce them (no p = 3,000 candidate job is run here).
QUAL <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(QUAL, "qual_common.R")); source(file.path(QUAL, "manifest.R")); source(file.path(QUAL, "qual_q7.R"))
RES <- list()
cat(sprintf("[%s] verifier controls\n", QUAL_MODE))

# ---- Q7 controls ----
z <- qual_data("S3", 1L)
A <- z$scores_test[, 1]; B <- z$scores_test[, 2]
dA <- q7_dab(A, A, B)[["joint"]]; dB <- q7_dab(B, A, B)[["joint"]]; dAB <- q7_dab(A + B, A, B)[["joint"]]
RES$Q7 <- c(D_A = dA, D_B = dB, D_AplusB = dAB, pass = dA > 0.9 && dB < -0.9 && abs(dAB) < 0.1)
cat(sprintf("Q7 controls (S3 r1 test set, %d columns, generator cor(A, B) %.3f): D_AB(A) = %.6f (> 0.9), D_AB(B) = %.6f (< -0.9), D_AB(A + B) = %.6f (|.| < 0.1) => %s\n",
            length(A), cor(A, B), dA, dB, dAB, if (RES$Q7[["pass"]] == 1) "PASS" else "FAIL"))

# ---- Q6 control and Q3 negative control ----
# Control namespace (outside every qualification and plumbing seed): dataset sim_qual(49,900,001, "S1"), three
# starts of the centre configuration at n = 175 with seed 49,900,101 (starts use 49,900,101-103). The first converged
# start (if any) is the control start; palm_fit_one is restarted
# from its endpoint (W, H, theta) with recording on; a converged-before-sweep-10 restart must give C_conv - C_t
# exactly 0 at t = 10, 50, 100, 500 (endpoint stands in), and q3_scope must route it to "negative".
CTRL_DATA_SEED <- 49900001; CTRL_FIT_SEED <- 49900101
iv <- assert_seeds_disjoint(qual_manifest())
stopifnot(!any(iv$lo <= CTRL_FIT_SEED + 2 & iv$hi >= CTRL_DATA_SEED))   # control seeds outside the manifest's seeds
zz <- sim_qual(CTRL_DATA_SEED, "S1"); cols <- 1:175
D <- palm_data(zz$X[, cols], zz$y[cols], zz$d[cols])
hp <- list(alpha = 0.05, lambda = 1, nu = 0.5, rho = 1, lambdaH = LAMBDA_H); ntop <- 150L
NC <- as.integer(Sys.getenv("NCORES", "1"))
f0 <- palm_fit(D, hp, k = 3L, seed = CTRL_FIT_SEED, nstarts = 3L, maxit = QUAL_MAXIT, tol = QUAL_TOL, keep_starts = TRUE, ncores = NC)
cat(sprintf("control fit (p %d, n %d, 3 starts): statuses %s, sweeps %s\n", D$p, D$n, paste(f0$starts$status, collapse = "/"),
            paste(f0$starts$iter, collapse = "/")))
ci <- which(f0$starts$status == "converged")[1]
if (is.na(ci)) {
  RES$Q6 <- NA; RES$Q3neg <- NA
  cat("Q6 / Q3-negative controls: INCONCLUSIVE, no control start converged (no converged iterate to restart from)\n")
} else {
  g0 <- f0$fits[[ci]]
  f <- palm_fit_one(D, hp, g0$W, g0$H, g0$theta, maxit = QUAL_MAXIT, tol = QUAL_TOL, record = TRUE)
  ck <- f$record$checkpoints
  cC <- vapply(c("10", "50", "100", "500"), function(t) iterate_testC(ck[[t]]$W, ck[[t]]$theta, D, hp, ntop, zz), 0)
  cE <- iterate_testC(f$W, f$theta, D, hp, ntop, zz)
  dif <- cE - cC
  early <- identical(f$status, "converged") && f$iter < 10L
  RES$Q6 <- c(iter = f$iter, converged = f$status == "converged", max_abs_diff = max(abs(dif)), pass = early && all(dif == 0))
  cat(sprintf("Q6 control: restart converged at sweep %d (%s); C_conv - C_t at t = 10/50/100/500: %s => %s\n",
              f$iter, f$status, paste(format(dif), collapse = "/"),
              if (!early) "INCONCLUSIVE (not converged before sweep 10)" else if (all(dif == 0)) "PASS" else "FAIL"))
  sc_restart <- q3_scope(f$status, f$iter); sc_orig <- q3_scope(g0$status, g0$iter)
  RES$Q3neg <- c(restart_negative = sc_restart == "negative", original_negative = sc_orig == "negative",
                 pass = sc_restart == "negative" && sc_orig == "negative")
  cat(sprintf("Q3 negative control: converged control start -> scope '%s', its restart -> '%s' (never routed to the detector) => %s\n",
              sc_orig, sc_restart, if (RES$Q3neg[["pass"]] == 1) "PASS" else "FAIL"))
}

# ---- Q3 positive control (ruling 18) ----
# Evaluable inputs, fixed from the archived trajectories independently of the detector: the D-Q3a continuation
# starts 1 and 5 (complete final-25% window with the documented scale collapse; exact replays, gates PASS). Both must
# be flagged. Not evaluable, reasons retained: start 9 (converged after 9 sweeps, no stalled collapse trajectory);
# start 16 (numerical termination after the collapse, flat final-quarter window). Specificity record: the PG-crawl
# replay of starts 1, 5, 9, 16 (late drift r-driven), not flagged. Results from q3pos/q3pos_dq3a_eval.R and
# q3pos/q3pos_eval.R (committed .rds; no p = 3,000 candidate job runs here).
q3d <- readRDS(file.path(QUAL, "q3pos", "dq3a_result.rds")); q3pg <- readRDS(file.path(QUAL, "q3pos", "q3pos_result.rds"))
stopifnot(identical(sort(names(q3d)), c("1", "5")))
flag <- vapply(q3d, function(x) isTRUE(x$detect$recurrence), TRUE)
pg_flag <- vapply(q3pg$res, function(x) isTRUE(x$detect$recurrence), TRUE)
RES$Q3pos <- c(start1 = flag[["1"]], start5 = flag[["5"]], pass = all(flag))
RES$Q3pos_record <- list(not_evaluable = c(`9` = "converged after 9 sweeps; no stalled collapse trajectory",
                                           `16` = "numerical termination after the collapse; flat final-quarter window"),
                         pg_specificity_flags = pg_flag)
cat(sprintf("Q3 positive control (ruling 18): D-Q3a start 1 %s, start 5 %s => %s; not evaluable: 9, 16 (reasons recorded); PG-crawl specificity record: flagged %d of %d\n",
            flag[["1"]], flag[["5"]], if (all(flag)) "PASS" else "FAIL", sum(pg_flag), length(pg_flag)))
pass_of <- function(x) isTRUE(!is.null(x) && "pass" %in% names(x) && x[["pass"]] == 1)
ok <- c(Q7 = pass_of(RES$Q7), Q6 = pass_of(RES$Q6), Q3neg = pass_of(RES$Q3neg), Q3pos = pass_of(RES$Q3pos))
cat(sprintf("CONTROLS: Q7 %s, Q6 %s, Q3 negative %s, Q3 positive %s%s\n", if (ok[["Q7"]]) "PASS" else "FAIL",
            if (ok[["Q6"]]) "PASS" else "FAIL/NOT RUN", if (ok[["Q3neg"]]) "PASS" else "FAIL/NOT RUN", if (ok[["Q3pos"]]) "PASS" else "FAIL",
            if (QUAL_MODE == "plumbing") " (plumbing data)" else ""))
RES$pass_all <- all(ok); RES$mode <- QUAL_MODE; RES$hashes <- HASHES
save_job(RES, file.path(OUT_DIR, "controls", "controls.rds"))
if (QUAL_MODE == "qual" && !all(ok)) { cat("controls not all passed: the qualification must not start (prespec section 7)\n"); quit(status = 1) }
