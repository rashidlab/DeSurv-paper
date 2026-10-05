# Q2 (amendment 3 section 6): estimator equivalence and non-regression at p = 150. On the p = 150 control settings,
# both solvers (w_solver "pg", the frozen projected-gradient step with maxit 5,000; w_solver "ncg", the amendment-3
# Newton-CG step with its per-start budget) run the frozen 30-start fit from the same initializations (palm_fit draws
# every start from the fit seed, so both solvers start from identical W0, H0 and share eps). Pass iff all five:
#   1. both fits converge under the same r < 1e-5 certificate (fit status "converged" for both);
#   2. the same start is selected;
#   3. |F_new - F_old| / max(1, |F_old|) <= 1e-6;
#   4. after optimal column matching (exhaustive over the k! permutations, which is the Hungarian optimum), every
#      W-column cosine >= 0.999;
#   5. fitted standardized risk scores eta~ (training-sample eta = Zc diag(1/sigma) b of the selected start, divided by
#      its sd): Pearson r >= 0.9999 and RMSE(eta~_new - eta~_old) / SD(eta~_old) <= 0.01.
# H and b differences are reported diagnostically only.
#
# SETTINGS NOT PINNED. The frozen amendment names "the p = 150 control settings" without fixing dataset, fit seed and
# hyperparameters, and the PI ruled (Stage B) that Q2 is not run, not even exploratory, until the settings are fixed
# prospectively. Candidate definitions are listed in AMENDMENT3_REPORT.md, Stage B. Set Q2_SETTINGS below (dataset
# seed, p, n, fit seed, hp) only after that ruling; until then the script stops.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
Q2_SETTINGS <- NULL   # e.g. list(sim_seed = , p = 150L, n = , fit_seed = , nstarts = 30L, hp = list(alpha = , lambda = , nu = , rho = , lambdaH = 1e-6))
if (is.null(Q2_SETTINGS)) stop("Q2: the p = 150 control settings are not pinned by the frozen documents; awaiting PI ruling (see AMENDMENT3_REPORT.md, Stage B)", call. = FALSE)
NC <- as.integer(Sys.getenv("T_NCORES", "15"))
S <- Q2_SETTINGS
s <- sim_data(S$sim_seed, S$p, S$n); D <- palm_data(s$X, s$y, s$d)
t0 <- proc.time()[["elapsed"]]
f_old <- palm_fit(D, S$hp, k = 3L, seed = S$fit_seed, nstarts = S$nstarts, ncores = NC, keep_starts = TRUE, w_solver = "pg")
t_old <- proc.time()[["elapsed"]] - t0; t0 <- proc.time()[["elapsed"]]
f_new <- palm_fit(D, S$hp, k = 3L, seed = S$fit_seed, nstarts = S$nstarts, ncores = NC, keep_starts = TRUE, w_solver = "ncg")
t_new <- proc.time()[["elapsed"]] - t0
perms <- function(k) if (k == 1) list(1L) else do.call(c, lapply(seq_len(k), function(i) lapply(perms(k - 1), function(q) c(i, setdiff(seq_len(k), i)[q]))))
match_cols <- function(Wo, Wn) {
  C <- crossprod(Wo, Wn); P <- perms(ncol(Wo)); best <- P[[which.max(sapply(P, function(q) sum(C[cbind(seq_along(q), q)])))]]
  list(perm = best, cos = C[cbind(seq_along(best), best)])
}
eta_std <- function(f) { z <- f$best; e <- as.vector(palm_scores(z$W, D, f$eps)$Zc %*% (z$b / z$sigma)); e / sd(e) }
c1 <- f_old$status == "converged" && f_new$status == "converged" && isTRUE(f_old$best$conv) && isTRUE(f_new$best$conv)
c2 <- isTRUE(f_old$best_start == f_new$best_start)
relF <- abs(f_new$best$F - f_old$best$F) / max(1, abs(f_old$best$F)); c3 <- relF <= 1e-6
mc <- match_cols(f_old$best$W, f_new$best$W); c4 <- all(mc$cos >= 0.999)
eo <- eta_std(f_old); en <- eta_std(f_new); pr <- cor(en, eo); rm <- sqrt(mean((en - eo)^2)) / sd(eo); c5 <- pr >= 0.9999 && rm <= 0.01
saveRDS(list(settings = S, f_old = f_old, f_new = f_new), file.path(OUT, "Q2_p150.rds"))
open_out("Q2_p150.txt")
options(width = 200)
cat("Q2 p = 150 estimator equivalence |", format(Sys.time()), "\n"); str(S)
cat(sprintf("wall: pg %.1f s, ncg %.1f s\n", t_old, t_new))
for (nm in c("pg", "ncg")) { f <- if (nm == "pg") f_old else f_new
  cat(sprintf("%s: fit status %s, selected start %s, F %.12g, start statuses converged %d / nonconverged %d / failure %d\n", nm, f$status,
              f$best_start, f$best$F, sum(f$starts$status == "converged"), sum(f$starts$status == "nonconverged"), sum(f$starts$status == "failure"))) }
cat(sprintf("1 both converge: %s\n2 same start selected: %s (pg %s, ncg %s)\n3 relative F difference %.3e (<= 1e-6): %s\n", c1, c2,
            f_old$best_start, f_new$best_start, relF, c3))
cat(sprintf("4 matched W-column cosines %s (all >= 0.999): %s\n5 eta~ Pearson %.8f (>= 0.9999), RMSE/SD %.3e (<= 0.01): %s\n",
            paste(sprintf("%.6f", mc$cos), collapse = ", "), c4, pr, rm, c5))
cat(sprintf("diagnostic: max |H diff| (matched rows) %.3e, b (matched) %s vs %s\n", max(abs(f_new$best$H[mc$perm, ] - f_old$best$H)),
            paste(signif(f_old$best$b, 4), collapse = ","), paste(signif(f_new$best$b[mc$perm], 4), collapse = ",")))
cat("\nper start (pg, ncg):\n"); print(cbind(f_old$starts[, c("start", "status", "F", "iter", "r")], f_new$starts[, c("status", "F", "iter", "hvp_total", "r")]), digits = 6)
verdict(c1 && c2 && c3 && c4 && c5)
sink()
