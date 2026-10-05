# Timing-only run (qualification prespec section 8, frozen under ruling 14). One dataset, sim_qual(49,100,001, "S1"),
# outside every qualification namespace. Fits: centre configuration at n = 175 (g = 5) and the reference fit at
# n = 218 (g = 12), 30 starts each, plus one inner-fold BO evaluation (5 inner-fold fits of 30 starts, 5-fold split
# of the 218 training columns stratified by event). Recording ON, as in qualification.
# Written per start: wall seconds and sweeps per second only. Status, F, residuals, iteration counts and C are
# dropped before anything is written; the output has no field for them.
# Run from the repo root: NCORES=15 Rscript analysis/estimator-gate/qualification/timing/timing_run.R
.libPaths(c("/home/naimrashid/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
QUAL <- normalizePath("analysis/estimator-gate/qualification"); RED <- normalizePath("analysis/estimator-gate/redesign")
source(file.path(QUAL, "sim_qual.R")); source(file.path(RED, "kkt_oracle_theta.R")); source(file.path(RED, "theta_lib.R"))
NC <- as.integer(Sys.getenv("NCORES", "15"))
hp <- list(alpha = 0.05, lambda = 1, nu = 0.5, rho = 1, lambdaH = 1e-6)
z <- sim_qual(49100001L, "S1")
fold_ids <- function(d, K, seed) { set.seed(seed); f <- integer(length(d))
  for (lev in c(0, 1)) { i <- which(d == lev); f[i] <- sample(rep_len(seq_len(K), length(i))) }; f }
fid <- fold_ids(z$d, 5L, 49100002L)
jobs <- c(list(list(name = "centre_n175", cols = 1:175, seed = 49100101L),
               list(name = "reference_n218", cols = 1:218, seed = 49100201L)),
          lapply(1:5, function(f) list(name = sprintf("inner_fold_%d", f), cols = which(fid != f), seed = 49100301L + 100L * (f - 1L))))
rows <- list(); t_all <- proc.time()[["elapsed"]]
for (jb in jobs) {
  D <- palm_data(z$X[, jb$cols], z$y[jb$cols], z$d[jb$cols])
  t0 <- proc.time()[["elapsed"]]
  fit <- palm_fit(D, hp, k = 3L, seed = jb$seed, nstarts = 30L, maxit = 5000L, ncores = NC, record = TRUE)
  wall <- proc.time()[["elapsed"]] - t0
  st <- fit$starts
  rows[[jb$name]] <- data.frame(job = jb$name, n = length(jb$cols), start = st$start, secs = st$secs,
                                sweeps_per_sec = st$iter / st$secs, job_wall = wall, cores = NC)
  rm(fit, st); gc(verbose = FALSE)
  cat(sprintf("%s done: wall %.0f s\n", jb$name, wall))
}
out <- do.call(rbind, rows); rownames(out) <- NULL
saveRDS(out, file.path(QUAL, "timing", "timing_wall.rds"))
sink(file.path(QUAL, "timing", "timing_wall.txt"))
cat(sprintf("Timing-only run, %s, cores %d, total wall %.0f s\n", format(Sys.time()), NC, proc.time()[["elapsed"]] - t_all))
print(do.call(rbind, lapply(split(out, out$job), function(d) data.frame(job = d$job[1], n = d$n[1],
  mean_secs = mean(d$secs), median_secs = median(d$secs), max_secs = max(d$secs),
  median_sweeps_per_sec = median(d$sweeps_per_sec), job_wall = d$job_wall[1]))), row.names = FALSE)
sink()
