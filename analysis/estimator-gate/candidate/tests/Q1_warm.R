# Q1 (amendment 3 section 6): binding realistic-size qualification. Warm start from the four stalled p = 3000,
# n = 218 end points of the frozen projected-gradient T7 fit (starts 1, 5, 9, 16; fit$fits[[i]] carrying W, H, b,
# eps and hp), full sweeps H, W, b with the Newton-CG W step at the per-start budget (1,000 sweeps, 100 HVPs per
# solve, 25,000 HVPs). Pass iff all four (prespecified 4 of 4): (1) r < 1e-5; (2) no F increase under the frozen
# scale-aware definition (counter, and a recount from the per-block trace, including the first sweep against the
# stalled end point; F_final <= F_stalled is checked as part of it); (3) within the per-start budget (no cap
# reached). A line-search or numerical failure fails Q1; a CG solve truncated at 100 HVPs does not. Loading and
# prediction agreement with the stalled end point are descriptive only.
# The stalled end points are read from tests/out/superseded/pre-amendment3/T7_fit.rds: the copy of T7_fit.rds made
# before the amendment-3 T7 rerun overwrote tests/out/T7_fit.rds.
# Residual-drift diagnostic at every CG exit (tests/drift_hook.R; its HVPs do not count toward the budget).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
source(file.path(CAND, "tests", "drift_hook.R"))
SRC <- file.path(CAND, "tests", "out", "superseded", "pre-amendment3", "T7_fit.rds")
x <- readRDS(SRC); fit <- x$fit; hp <- fit$hp
SS <- c(1L, 5L, 9L, 16L)
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
pre <- vapply(SS, function(i) { g <- fit$fits[[i]]; g$iter == 5000L && g$status == "nonconverged" && identical(g$eps, hp$eps) &&
                                  isTRUE(all.equal(palm_F(g$W, g$H, g$b, D, hp), g$F, tolerance = 1e-14)) }, TRUE)
one <- function(i) {
  g <- fit$fits[[i]]; drift_reset(paste0("start", i))
  f <- palm_fit_one(D, hp, g$W, g$H, g$b, trace = TRUE)          # w_solver "ncg", maxit = sweep budget 1000, hvp_cap 25000
  tr <- f$trace; tol13 <- function(Fo) 1e-13 * pmax(1, abs(Fo))
  rc <- sum(tr[, "FH"] - tr[, "Fstart"] > tol13(tr[, "Fstart"]), na.rm = TRUE) +
        sum(tr[, "FW"] - tr[, "FH"] > tol13(tr[, "FH"]), na.rm = TRUE) + sum(tr[, "Fb"] - tr[, "FW"] > tol13(tr[, "FW"]), na.rm = TRUE)
  eta <- function(W, b) as.vector(palm_scores(W, D, hp$eps)$Zc %*% (b / palm_scores(W, D, hp$eps)$sigma))
  list(f = f, drift = drift_rows(), row = data.frame(start = i, status = f$status, reason = f$reason, sweeps = f$iter,
    hvp_total = f$hvp_total, r = f$r[["r"]], rW = f$r[["rW"]], rH = f$r[["rH"]], rb = f$r[["rb"]], F_stalled = g$F, F_final = f$F,
    F_improve = g$F - f$F, inc_counter = sum(f$ninc), inc_recount = rc, first_sweep_vs_stalled = tr[1, "Fstart"] - g$F,
    cg_forcing = f$cg_exits[["forcing"]], cg_ceiling = f$cg_exits[["hvp_ceiling"]], cg_total_cap = f$cg_exits[["hvp_total_cap"]],
    cg_curvature = f$cg_exits[["curvature"]], negative_curvature = f$n_negative_curvature, nc_first_iter = f$n_nc_first_iter,
    w_noop = f$w_exits[["w_converged"]], w_newton = f$w_exits[["newton"]], w_ls_failed = f$w_exits[["line_search_failed"]],
    max_t = suppressWarnings(max(tr[, "t"], na.rm = TRUE)), halvings = f$w_backtracks, enter = sum(tr[, "n_enter"], na.rm = TRUE),
    leave = sum(tr[, "n_leave"], na.rm = TRUE), cos_min = min(colSums(f$W * g$W)), eta_pearson = cor(eta(f$W, f$b), eta(g$W, g$b)),
    secs = f$secs, stringsAsFactors = FALSE))
}
t0 <- Sys.time()
res <- parallel::mclapply(SS, function(i) tryCatch(one(i), error = function(e) conditionMessage(e)), mc.cores = length(SS), mc.preschedule = FALSE)
bad <- !vapply(res, is.list, TRUE)
tab <- do.call(rbind, lapply(res[!bad], `[[`, "row"))
dr <- do.call(rbind, lapply(res[!bad], `[[`, "drift"))
saveRDS(list(tab = tab, drift = dr, fits = lapply(res[!bad], `[[`, "f"), src = SRC, hp = hp), file.path(OUT, "Q1_warm.rds"))
open_out("Q1_warm.txt")
options(width = 250)
cat("Q1 warm-started p = 3000 qualification |", format(Sys.time()), "| wall", format(Sys.time() - t0), "\n")
cat(sprintf("source %s; hp alpha=%.6g lambda=%g nu=%g rho=%.6g lambdaH=%g eps=%.6e; p=%d n=%d events %d\n", SRC, hp$alpha, hp$lambda,
            hp$nu, hp$rho, hp$lambdaH, hp$eps, D$p, D$n, D$ne))
cat(sprintf("precondition (stalled end points: 5000 PG iterations, nonconverged, fit eps, stored F reproduced): %s\n", paste(pre, collapse = ",")))
if (any(bad)) { cat("ERRORS:\n"); print(unlist(res[bad])) }
cat("\nper start:\n"); print(tab[, c("start", "status", "reason", "sweeps", "hvp_total", "r", "rW", "rH", "rb", "F_stalled", "F_final", "F_improve",
                                    "inc_counter", "inc_recount", "first_sweep_vs_stalled")], digits = 6)
cat("\nW-step and CG detail per start:\n"); print(tab[, c("start", "w_newton", "w_noop", "w_ls_failed", "cg_forcing", "cg_ceiling", "cg_total_cap",
                                                         "cg_curvature", "negative_curvature", "nc_first_iter", "max_t", "halvings", "enter", "leave", "secs")])
cat("\ndescriptive agreement with the stalled end point (not criteria): min column cosine, Pearson of fitted eta\n")
print(tab[, c("start", "cos_min", "eta_pearson")], digits = 6)
c1 <- tab$status == "converged" & tab$r < 1e-5
c2 <- tab$inc_counter == 0 & tab$inc_recount == 0 & tab$first_sweep_vs_stalled <= 1e-13 * pmax(1, abs(tab$F_stalled)) & tab$F_final <= tab$F_stalled
c3 <- tab$sweeps <= 1000 & tab$hvp_total <= 25000 & tab$reason != "work_cap_reached"
cat(sprintf("\ncriteria per start (1 certificate, 2 monotone incl. F_final <= F_stalled, 3 budget): %s\n",
            paste(sprintf("start %d: %s/%s/%s", tab$start, c1, c2, c3), collapse = "; ")))
npass <- sum(c1 & c2 & c3)
ntc <- if (is.null(dr)) 0L else drift_summary(dr, "Q1 (every CG solve)")$n_tier_c
cat(sprintf("\nQ1 TOTAL: %d of 4 starts pass all three criteria (required 4 of 4); residual-drift tier (c) events %d (binding, PI ruling)\n", npass, ntc))
verdict(!any(bad) && all(pre) && nrow(tab) == 4 && npass == 4 && ntc == 0)
sink()
