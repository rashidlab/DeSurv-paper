# Coordinator diagnostic: does r_W keep falling past maxit at the box midpoint (p = 150, n = 60)?
CAND <- "/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/estimator-gate/candidate"
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
for (f in c("kkt_oracle.R", "palm_lib.R", "sim.R")) source(file.path(CAND, f))
s <- sim_data(401, 150, 60); D <- palm_data(s$X, s$y, s$d)
eps <- palm_eps(lapply(1:5, function(st) proj_S(palm_draw(D, 3L, 7L, st))), D)
hp <- list(alpha = plogis((qlogis(1e-4) + qlogis(0.9)) / 2), rho = 10^0.5, lambda = 1, nu = 0.475, lambdaH = 1e-6, eps = eps)
W <- proj_S(palm_draw(D, 3L, 7L, 1)); H <- palm_init_H(W, D, hp)
t0 <- Sys.time()
f <- palm_fit_one(D, hp, W, H, rep(0, 3), maxit = as.integer(commandArgs(TRUE)[1]), tol = 1e-5, trace = TRUE)
tr <- as.data.frame(f$trace); saveRDS(f, file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "long_fit.rds"))
cat("status", f$status, "iter", f$iter, "secs", round(as.numeric(Sys.time() - t0, units = "secs")), "\n")
ix <- unique(round(exp(seq(0, log(nrow(tr)), length.out = 25))))
print(signif(tr[ix, intersect(c("iter", "Fstart", "r", "rW", "rH", "rb", "t", "nbt"), names(tr))], 5))
