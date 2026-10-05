# Shared setup for tests/T*.R. Run from any directory: Rscript tests/T1_gradients.R
# MUTANT=<id> sources tests/mutants.R after the library (T10) and redirects output to tests/out/mutants/<id>/.
# T_QUICK=1 shrinks the long tests (T4, T7) for mutation runs.
.args <- commandArgs(FALSE)
CAND <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])), ".."))
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))   # R 4.6 user library (DeSurv 1.0.1)
source(file.path(CAND, "kkt_oracle.R")); source(file.path(CAND, "palm_lib.R")); source(file.path(CAND, "sim.R")); source(file.path(CAND, "ref_objective.R"))
MUTANT <- Sys.getenv("MUTANT", ""); QUICK <- Sys.getenv("T_QUICK", "") == "1"
OUT <- file.path(CAND, "tests", "out")
if (nzchar(MUTANT)) { source(file.path(CAND, "tests", "mutants.R")); apply_mutant(MUTANT); OUT <- file.path(OUT, "mutants", MUTANT) }
if (QUICK && !nzchar(MUTANT)) OUT <- file.path(OUT, "quick")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
open_out <- function(name) { f <- file.path(OUT, name); sink(f, split = TRUE); if (nzchar(MUTANT)) cat("MUTANT", MUTANT, "\n"); f }
verdict <- function(pass) cat(sprintf("VERDICT: %s\n", if (isTRUE(pass)) "PASS" else "FAIL"))
GRID <- expand.grid(alpha = c(1e-4, 0.3, 0.9), rho = c(1e-2, 1, 1e3), lambda = c(1e-3, 1, 1e3), nu = c(0, 0.95))
