# Coordinator diagnostic: continue T7 (n = 218) starts past maxit 5000; is r_W still falling or at a floor?
CAND <- "/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/estimator-gate/candidate"
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
for (f in c("kkt_oracle.R", "palm_lib.R", "sim.R")) source(file.path(CAND, f))
x <- readRDS(file.path(CAND, "tests/out/T7_fit.rds")); fit <- x$fit
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
hp <- fit$hp; stopifnot(isTRUE(all.equal(hp$eps, fit$eps)))
ss <- c(1, 5, 9, 16)
res <- parallel::mclapply(ss, function(i) {
  g <- fit$fits[[i]]
  f <- palm_fit_one(D, hp, g$W, g$H, g$b, maxit = 15000L, tol = 1e-5, trace = TRUE)
  tr <- as.data.frame(f$trace); ix <- unique(c(1, seq(500, nrow(tr), 500), nrow(tr)))
  cbind(start = i, tr[ix, intersect(c("iter", "Fstart", "r", "rW", "rH", "t", "nbt"), names(tr))])
}, mc.cores = 4)
saveRDS(res, file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "cont.rds"))
for (r in res) if (is.data.frame(r)) print(signif(r, 6)) else print(r)
