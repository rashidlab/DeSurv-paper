# Test: bo_engine_copy() proposes exactly the same points as the installed DeSurv 1.0.1 desurv_cv_bayesopt() for a
# fixed toy objective and seed. The package function is run unchanged, with its desurv_cv() replaced by the toy
# (function environment rebound to a child of the DeSurv namespace); the copy gets the same toy as cv_objective.
# The toy returns NA (an infeasible evaluation) on part of the box, so the NA-exclusion path is exercised.
# Run: Rscript analysis/estimator-gate/qualification/bo_engine_test.R   (no data, no fit)
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
QUAL <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(QUAL, "bo_engine_copy.R"))
sha <- utils::packageDescription("DeSurv")$RemoteSha
cat(sprintf("DeSurv %s, RemoteSha %s\n", utils::packageVersion("DeSurv"), sha))
bounds <- list(logit_alpha = list(lower = qlogis(1e-4), upper = qlogis(0.9), type = "continuous", scale = "linear"),
               lambda = list(lower = 1e-3, upper = 1e3, scale = "log10"),
               nu = list(lower = 0, upper = 0.95, type = "continuous"),
               rho = list(lower = 1e-2, upper = 1e3, scale = "log10"),
               ntop = list(lower = 50, upper = 300, type = "integer"))
toy_value <- function(a) {
  if (a$nu > 0.85) return(NA_real_)                                   # infeasible region
  0.7 - (plogis(a$logit_alpha) - 0.1)^2 - (log10(a$lambda) - 0.5)^2 / 20 - (a$nu - 0.4)^2 -
    (log10(a$rho))^2 / 30 - ((a$ntop - 140) / 300)^2
}
toy_result <- function(a) list(summary = data.frame(mean_cindex = toy_value(a)), diagnostics = NULL)
run_args <- list(X = NULL, y = NULL, d = NULL, bo_bounds = bounds, n_init = 8L, n_iter = 10L, candidate_pool = 4000L,
                 exploration_weight = 0.01, seed = 4242L, verbose = FALSE)
pkg <- DeSurv::desurv_cv_bayesopt
e <- new.env(parent = asNamespace("DeSurv")); e$desurv_cv <- function(...) toy_result(list(...))
environment(pkg) <- e
a <- do.call(pkg, run_args)
b <- do.call(bo_engine_copy, c(list(cv_objective = toy_result), run_args))
pc <- c("eval_id", "stage", "iteration", names(bounds), "mean_cindex", "status")
same_hist <- identical(a$history[, pc], b$history[, pc])
same_best <- identical(a$best, b$best)
cat(sprintf("evaluations: package %d, copy %d; NA (infeasible) evaluations %d; stages %s\n", nrow(a$history), nrow(b$history),
            sum(is.na(a$history$mean_cindex)), paste(names(table(a$history$stage)), table(a$history$stage), collapse = ", ")))
cat(sprintf("proposed points and scores identical: %s; best identical: %s\n", same_hist, same_best))
# the copy's body differs from the package body only in the two documented edits
bp <- deparse(body(DeSurv::desurv_cv_bayesopt)); bc <- deparse(body(bo_engine_copy))
nd <- sum(bp != bc)
cat(sprintf("body lines: package %d, copy %d, differing lines %d (expected 1: the objective call)\n", length(bp), length(bc), nd))
cat(sprintf("formals: copy = c('cv_objective', package formals): %s\n",
            identical(names(formals(bo_engine_copy)), c("cv_objective", names(formals(DeSurv::desurv_cv_bayesopt))))))
ok <- same_hist && same_best && length(bp) == length(bc) && nd == 1L && identical(sha, "afb00d5fdf2621159200b3966e22d9ea12e16179")
cat(sprintf("VERDICT: %s\n", if (ok) "PASS" else "FAIL"))
