# Job manifest (prespec sections 2, 5 and 9.1, PI ruling 15). qual: grid 3 scenarios x 5 replicates x 12 fits = 180;
# tuning S1 3 + S2 3 + S3 10 = 16 (S3b dropped). Fail-fast stages (prespec section 9.1, each evaluated before the
# next is submitted): stage 1 = the 180 grid and reference fits (with the verifier controls and the realized-ratio
# check on S3 replicates 1-5); stage 2 = S3 tuning r1-10; stage 3 = S1 tuning r1-3; stage 4 = S2 tuning r1-3.
# plumbing (compute-limited pipe test, 3 starts per fit): stage 1 = grid S1, S2, S3 replicate 1 at g = 5 (centre)
# and g = 12 (reference), 6 fits; stage 2 = S3 tuning r1-2 (two replicates so Q7's SD exists); stage 3 = S1 r1;
# stage 4 = S2 r1. Requires qual_common.R (QUAL_MODE, OUT_DIR).
# Rscript manifest.R prints the counts per stage and runs the seed-disjointness assertion for the full list.
GRID_CONFIGS <- data.frame(
  g      = 1:12,
  alpha  = c(1e-3, 1e-3, 1e-3, 0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.05, 0.05, 0.05),
  lambda = c(1e-2, 1, 1e2, 1e-2, 1, 1e2, 1e-2, 1, 1e2, 1, 1, 1),
  rho    = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1e-2, 1e3, 1),
  nu     = 0.5, ntop = 150L,
  n      = c(rep(175L, 11), 218L),
  role   = c(rep("grid", 4), "centre", rep("grid", 6), "reference"),
  stringsAsFactors = FALSE)
tune_budget <- function() if (QUAL_MODE == "qual") list(n_init = 20L, n_iter = 40L) else list(n_init = 2L, n_iter = 1L)
qual_manifest <- function() {
  grid_reps <- if (QUAL_MODE == "qual") list(S1 = 1:5, S2 = 1:5, S3 = 1:5) else list(S1 = 1L, S2 = 1L, S3 = 1L)
  tune_reps <- if (QUAL_MODE == "qual") list(S3 = 1:10, S1 = 1:3, S2 = 1:3) else list(S3 = 1:2, S1 = 1L, S2 = 1L)
  grid <- do.call(rbind, lapply(names(grid_reps), function(s) do.call(rbind, lapply(grid_reps[[s]], function(r)
    cbind(data.frame(scenario = s, sc = SCEN_IDX[[s]], r = r, stringsAsFactors = FALSE), GRID_CONFIGS)))))
  if (QUAL_MODE == "plumbing") grid <- grid[grid$g %in% c(5L, 12L), ]
  grid$out <- file.path(OUT_DIR, "grid", sprintf("%s_r%02d_g%02d.rds", grid$scenario, grid$r, grid$g))
  tune <- do.call(rbind, lapply(names(tune_reps), function(s)
    data.frame(scenario = s, sc = SCEN_IDX[[s]], r = tune_reps[[s]], stringsAsFactors = FALSE)))
  tune$out <- file.path(OUT_DIR, "tune", sprintf("%s_r%02d.rds", tune$scenario, tune$r))
  grid$stage <- 1L; tune$stage <- unname(c(S3 = 2L, S1 = 3L, S2 = 4L)[tune$scenario])
  rownames(grid) <- NULL; rownames(tune) <- NULL
  if (QUAL_MODE == "qual") stopifnot(nrow(grid) == 180L, nrow(tune) == 16L, identical(as.vector(table(tune$stage)), c(10L, 3L, 3L)))
  list(grid = grid, tune = tune, mode = QUAL_MODE, nstarts = QUAL_NSTARTS, budget = tune_budget(), stages = 1:4)
}
# Output files of one stage (the evaluator's completeness check uses stages 1..k).
stage_files <- function(man, k) c(man$grid$out[man$grid$stage %in% k], man$tune$out[man$tune$stage %in% k])
if (sys.nframe() == 0L) {   # Rscript manifest.R
  QUAL <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  source(file.path(QUAL, "qual_common.R"))
  man <- qual_manifest()
  cat(sprintf("mode %s: grid jobs %d, tuning jobs %d, starts per fit %d, BO n_init %d n_iter %d\n", man$mode,
              nrow(man$grid), nrow(man$tune), man$nstarts, man$budget$n_init, man$budget$n_iter))
  print(table(man$grid$scenario)); print(table(stage = man$tune$stage, scenario = man$tune$scenario))
  assert_seeds_disjoint(man)
  cat(sprintf("outputs present: grid %d of %d, tuning %d of %d\n", sum(file.exists(man$grid$out)), nrow(man$grid),
              sum(file.exists(man$tune$out)), nrow(man$tune)))
}
