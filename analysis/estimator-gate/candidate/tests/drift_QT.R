# Residual-drift diagnostic on a seeded sample of Q-T solves (PI ruling, Stage B; diagnostic only, tier (c) binding).
# Sample: set.seed(20261005), 12 of the 135 T4 runs (3 datasets x 9 corners x 5 starts), rebuilt exactly as in
# tests/T4_monotone.R (same data, corners, eps pooling, draws) and run to the same stopping rule with every CG exit
# instrumented by tests/drift_hook.R. Writes tests/out/drift_QT.txt and .rds.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
source(file.path(CAND, "tests", "drift_hook.R"))
NC <- as.integer(Sys.getenv("T_NCORES", "12")); SEED <- 20261005L; NSAMP <- 12L
corners <- rbind(expand.grid(alpha = c(1e-4, 0.9), rho = c(1e-2, 1e3), lnu = 1:2), data.frame(alpha = 0.3, rho = 1, lnu = 3))
LN <- list(c(1e-3, 0), c(1, 0.95), c(1, 0.5))
corners$lambda <- sapply(corners$lnu, function(i) LN[[i]][1]); corners$nu <- sapply(corners$lnu, function(i) LN[[i]][2])
NDS <- 3L; NST <- 5L
dsets <- lapply(seq_len(NDS), function(i) { s <- sim_data(400 + i, 150, 60); list(s = s, D = palm_data(s$X, s$y, s$d)) })
jobs <- expand.grid(start = seq_len(NST), corner = seq_len(nrow(corners)), ds = seq_len(NDS))
epsv <- list()
for (di in seq_len(NDS)) for (ci in seq_len(nrow(corners)))
  epsv[[paste(di, ci)]] <- palm_eps(lapply(seq_len(NST), function(st) proj_S(palm_draw(dsets[[di]]$D, 3L, 7L, st))), dsets[[di]]$D)
set.seed(SEED); pick <- sort(sample(nrow(jobs), NSAMP))
one <- function(j) {
  di <- jobs$ds[j]; ci <- jobs$corner[j]; st <- jobs$start[j]; D <- dsets[[di]]$D
  hp <- list(alpha = corners$alpha[ci], rho = corners$rho[ci], lambda = corners$lambda[ci], nu = corners$nu[ci],
             lambdaH = 1e-6, eps = epsv[[paste(di, ci)]])
  drift_reset(sprintf("T4 job %d (ds %d, corner %d, start %d)", j, di, ci, st))
  W <- proj_S(palm_draw(D, 3L, 7L, st)); H <- palm_init_H(W, D, hp)
  f <- palm_fit_one(D, hp, W, H, rep(0, 3), maxit = 1000L, tol = 1e-5)
  list(drift = drift_rows(), row = data.frame(job = j, ds = di, corner = ci, start = st, status = f$status, reason = f$reason,
                                              sweeps = f$iter, hvp_total = f$hvp_total, r = f$r[["r"]], stringsAsFactors = FALSE))
}
res <- parallel::mclapply(pick, one, mc.cores = NC, mc.preschedule = FALSE)
bad <- !vapply(res, is.list, TRUE)
tab <- do.call(rbind, lapply(res[!bad], `[[`, "row")); dr <- do.call(rbind, lapply(res[!bad], `[[`, "drift"))
saveRDS(list(seed = SEED, pick = pick, tab = tab, drift = dr), file.path(OUT, "drift_QT.rds"))
open_out("drift_QT.txt")
options(width = 200)
cat("Residual-drift diagnostic, seeded Q-T sample |", format(Sys.time()), "\n")
cat(sprintf("sample: set.seed(%d), %d of %d T4 runs: jobs %s\n", SEED, NSAMP, nrow(jobs), paste(pick, collapse = ",")))
if (any(bad)) { cat("ERRORS:\n"); print(unlist(res[bad])) }
print(tab, digits = 4)
ds <- drift_summary(dr, "Q-T sample (T4 runs)")
sink()
