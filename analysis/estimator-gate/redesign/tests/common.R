# Shared setup for the redesign smoke tests tests/S*.R (pattern of candidate/tests/common.R).
# Run from any directory: Rscript tests/S1_gradients.R. Simulated data only (candidate/sim.R), p <= 300.
.args <- commandArgs(FALSE)
RED <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])), ".."))
CAND <- normalizePath(file.path(RED, "..", "candidate"))   # read-only: sim.R, and the frozen reference code in S3/S10
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))   # R 4.6 user library (DeSurv 1.0.1)
source(file.path(RED, "kkt_oracle_theta.R")); source(file.path(RED, "theta_lib.R"))
source(file.path(CAND, "sim.R")); source(file.path(RED, "ref_objective_theta.R"))
OUT <- file.path(RED, "tests", "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
open_out <- function(name) { f <- file.path(OUT, name); sink(f, split = TRUE); invisible(f) }
verdict <- function(pass) cat(sprintf("VERDICT: %s\n", if (isTRUE(pass)) "PASS" else "FAIL"))
# check tally: every test prints checks run, passed, max error
TALLY <- new.env(); TALLY$n <- 0L; TALLY$pass <- 0L; TALLY$maxerr <- 0
tally <- function(label, err, tol) {
  ok <- is.finite(err) & err <= tol
  TALLY$n <- TALLY$n + length(err); TALLY$pass <- TALLY$pass + sum(ok)
  TALLY$maxerr <- max(TALLY$maxerr, err, na.rm = TRUE)
  cat(sprintf("  %-55s checks %4d  passed %4d  max err %.3e  (tol %.0e)\n", label, length(err), sum(ok), max(err), tol))
  invisible(all(ok))
}
total <- function(id) {
  cat(sprintf("\n%s TOTAL: checks run %d, passed %d, max error %.3e\n", id, TALLY$n, TALLY$pass, TALLY$maxerr))
  verdict(TALLY$n > 0 && TALLY$n == TALLY$pass)
}
# The frozen candidate code, loaded into its own environment so nothing clashes with theta_lib.R.
load_candidate <- function() {
  e <- new.env()
  sys.source(file.path(CAND, "kkt_oracle.R"), envir = e); sys.source(file.path(CAND, "palm_lib.R"), envir = e)
  sys.source(file.path(CAND, "ref_objective.R"), envir = e)
  e
}
rand_W <- function(p, k = 3L) proj_S(matrix(runif(p * k), p, k))           # interior point of S_+^k (all W > 0)
rand_hp <- function() list(alpha = plogis(runif(1, qlogis(1e-4), qlogis(0.9))), lambda = 10^runif(1, -3, 3),
                           nu = runif(1, 0, 0.95), rho = 10^runif(1, -2, 3), lambdaH = 1e-6)
harrell <- function(score, y, d) {                         # simulated data only (copied from candidate T11)
  num <- 0; den <- 0
  for (i in which(d == 1)) { j <- y > y[i]; den <- den + sum(j); num <- num + sum(score[i] > score[j]) + 0.5 * sum(score[i] == score[j]) }
  num / den
}
