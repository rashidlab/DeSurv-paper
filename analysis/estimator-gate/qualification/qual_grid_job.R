# One grid or reference fit (prespec section 5 table, g = 1-12; replicates 1-5 of S1, S2, S3).
# Usage: QUAL_MODE=qual|plumbing NCORES=<n> Rscript qual_grid_job.R <scenario> <replicate> <g>
# Training columns 1..n of the replicate's sim_qual draw (n = 175 for g = 1-11, 218 for g = 12); test = the draw's
# 1,000 test columns. 30 starts (plumbing: 3), maxit 5,000, recording on. Writes one .rds under <OUT_DIR>/grid/.
QUAL <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(QUAL, "qual_common.R")); source(file.path(QUAL, "manifest.R"))
a <- commandArgs(TRUE)
if (length(a) != 3L) stop("usage: qual_grid_job.R <scenario> <replicate> <g>", call. = FALSE)
scen <- a[1]; r <- as.integer(a[2]); g <- as.integer(a[3])
NC <- as.integer(Sys.getenv("NCORES", "1"))
man <- qual_manifest()
job <- man$grid[man$grid$scenario == scen & man$grid$r == r & man$grid$g == g, ]
if (nrow(job) != 1L) stop(sprintf("job %s r%d g%d is not in the %s manifest", scen, r, g, QUAL_MODE), call. = FALSE)
if (file.exists(job$out)) stop(sprintf("output exists, not overwritten: %s", job$out), call. = FALSE)
require_controls(); require_precheck(scen, 1L)   # stage 1 (prespec section 9.1)
t0 <- proc.time()[["elapsed"]]
z <- qual_data(scen, r)
cols <- seq_len(job$n)
D <- palm_data(z$X[, cols], z$y[cols], z$d[cols])
cat(sprintf("[%s] grid job %s r%d g%d (%s): read 1 dataset (seed %d, p %d, train n %d of %d, test %d); s0 %.6f; cores %d\n",
            QUAL_MODE, scen, r, g, job$role, z$seed, D$p, D$n, ncol(z$X), ncol(z$X_test), D$s0, NC))
hp <- list(alpha = job$alpha, lambda = job$lambda, nu = job$nu, rho = job$rho)
res <- qual_fit(D, hp, seed_grid(job$sc, r, g), ncores = NC, keep = "full", z = z, ntop = job$ntop)
res$job <- job; res$kind <- "grid"; res$mode <- QUAL_MODE; res$hashes <- HASHES; res$data <- data_summary(z)
res$wall <- proc.time()[["elapsed"]] - t0
cat(sprintf("fit status %s, returned start %s; %s; wall %.0f s\n", res$status, res$best_start, fit_counts(res), res$wall))
save_job(res, job$out)
