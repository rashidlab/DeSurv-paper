# Q3 positive control, phase 1 (prespec section 7, implementation correction recorded in q3pos/README.md).
# Replays the FAILED candidate's original projected-gradient trajectory for one T7 start at p = 3,000, n = 218
# (sim_data(700, 3000, 218), T7 hyperparameters, fit seed 123, eps pooled over the 30 T7 starts, w_solver "pg",
# maxit 5,000), with an instrumentation-only hook on w_step that stores (W after the W step, b entering that sweep)
# at every 10th sweep of the Q3 window and at sweeps 3,750, 4,375 and 5,000. The hook changes no computation.
# Usage (repo root): Rscript analysis/estimator-gate/qualification/q3pos/q3pos_replay.R <start>
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
s <- as.integer(commandArgs(TRUE)[1]); stopifnot(s %in% c(1L, 5L, 9L, 16L))
CAND <- normalizePath("analysis/estimator-gate/candidate")
source(file.path(CAND, "kkt_oracle.R")); source(file.path(CAND, "palm_lib.R")); source(file.path(CAND, "sim.R"))
z <- sim_data(700, 3000L, 218L); D <- palm_data(z$X, z$y, z$d)
hp <- list(alpha = plogis((qlogis(1e-4) + qlogis(0.9)) / 2), lambda = 1, nu = 0.475, rho = 10^0.5, lambdaH = 1e-6)
Ws <- lapply(1:30, function(i) proj_S(palm_draw(D, 3L, 123L, i)))
hp$eps <- palm_eps(Ws, D)
KEEP <- c(3750L, seq(3760L, 5000L, by = 10L), 4375L)
SNAP <- new.env(); SNAP$sweep <- 0L; SNAP$W <- list(); SNAP$b <- list()
w_step_plain <- w_step
w_step <- function(W, H, b, D, hp, t_prev, t_min = 1e-20, defs = PALM_DEFS) {
  z <- w_step_plain(W, H, b, D, hp, t_prev, t_min = t_min, defs = defs)
  SNAP$sweep <- SNAP$sweep + 1L
  if (SNAP$sweep %in% KEEP) { k <- as.character(SNAP$sweep); SNAP$W[[k]] <- z$W; SNAP$b[[k]] <- b }
  z
}
W <- Ws[[s]]; H <- palm_init_H(W, D, hp)
t0 <- proc.time()[["elapsed"]]
f <- palm_fit_one(D, hp, W, H, rep(0, 3), maxit = 5000L, tol = 1e-5, w_solver = "pg")
out <- list(start = s, status = f$status, iter = f$iter, w_sweeps_seen = SNAP$sweep, W = SNAP$W, b = SNAP$b,
            eps = hp$eps, hp = hp, secs = proc.time()[["elapsed"]] - t0)
saveRDS(out, file.path(CAND, "..", "qualification", "q3pos", sprintf("replay_start%02d.rds", s)))
cat(sprintf("start %d: status %s, sweeps %s, W steps seen %d, snapshots %d, %.0f s\n", s, f$status, f$iter, SNAP$sweep,
            length(SNAP$W), out$secs))
