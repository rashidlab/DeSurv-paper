# Q3 positive control, diagnosis (PI ruling 17: diagnose only the failed control). Exact replay of the D-Q3a
# continuation (candidate/tests/DQ3a_continuation.R, unchanged solver and settings) for starts 1 and 5, with a
# recording-only hook keeping W and b at the b-block exits of the final-25% window points (a, m, n) and every 10th
# sweep in between. Gate: final F identical to the committed DQ3a_continuation.rds. Then the unchanged detector.
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
CAND <- normalizePath("analysis/estimator-gate/candidate"); OUTC <- file.path(CAND, "tests", "out")
source(file.path(CAND, "kkt_oracle.R")); source(file.path(CAND, "palm_lib.R")); source(file.path(CAND, "sim.R"))
hp <- readRDS(file.path(OUTC, "superseded", "pre-amendment3", "T7_fit.rds"))$fit$hp
dq1 <- readRDS(file.path(OUTC, "DQ1_block_replay.rds")); ref <- readRDS(file.path(OUTC, "DQ3a_continuation.rds"))
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
SS <- c(1L, 5L, 16L); solve_b_plain <- solve_b
one <- function(i) {
  n_ref <- nrow(ref$res[[match(i, SS)]]$rec); a <- floor(0.75 * n_ref); m <- floor((a + n_ref) / 2)
  keep <- sort(unique(c(a, m, n_ref, seq(a, n_ref, by = 10L))))
  REC <- new.env(); REC$n <- 0L; REC$W <- list(); REC$b <- list()
  solve_b <<- function(b0, W, D, hp, L0 = 1, tol = 1e-10, maxit = 1e5) {
    res <- solve_b_plain(b0, W, D, hp, L0, tol, maxit); REC$n <- REC$n + 1L
    if (REC$n %in% keep) { REC$W[[as.character(REC$n)]] <- W }
    REC$pend <- list(b0 = b0, b1 = res$b); res
  }
  acc_plain <- PALM_DEFS$b_accept
  PALM_DEFS$b_accept <<- function(F_new, F_old) { dec <- acc_plain(F_new, F_old)
    if (REC$n %in% keep) REC$b[[as.character(REC$n)]] <- if (isTRUE(dec)) REC$pend$b1 else REC$pend$b0; dec }
  e <- dq1$res[[match(i, SS)]]$f
  f <- palm_fit_one(D, hp, e$W, e$H, e$b, maxit = 19000L, trace = TRUE, hvp_cap = 1e6L)
  list(start = i, a = a, m = m, n = n_ref, nb = REC$n, W = REC$W, b = REC$b, eps = hp$eps,
       gate = identical(f$F, ref$res[[match(i, SS)]]$f$F))
}
res <- parallel::mclapply(c(1L, 5L), one, mc.cores = 2L, mc.preschedule = FALSE)
saveRDS(res, file.path(normalizePath("analysis/estimator-gate/qualification/q3pos"), "dq3a_replay.rds"))
for (r in res) cat(sprintf("start %d: b-blocks %d (ref %d), gate %s, snapshots %d\n", r$start, r$nb, r$n, r$gate, length(r$W)))
