# Simulated data for the candidate tests: the repo generator (R/simulation_functions/, scenario R0
# defaults, marker and background sizes scaled to p), within-sample ranks (average ties) as in the
# pipeline, survival times rounded up to integers to create tied times. Simulated data only.
# Requires CAND (path of analysis/estimator-gate/candidate) to be set by the caller.
local({
  sd <- file.path(CAND, "..", "..", "..", "R", "simulation_functions")
  for (f in c("simulate_W.R", "simulate_H.R", "simulate_X.R", "simulate_survival.R",
              "scenario_defaults.R", "simulate_data.R")) source(file.path(sd, f), local = globalenv())
})

sim_data <- function(seed, p, n, k = 3L, beta = c(2, 1, 0)[seq_len(k)], tie_unit = 1) {
  s <- suppressWarnings(simulate_desurv_scenario(
    scenario = "R0", G = p, N = n, K = k, beta = beta,
    markers_per_factor = max(5L, round(150 * p / 3000)), B_size = round(500 * p / 3000), seed = seed))
  X <- apply(s$X, 2, rank, ties.method = "average"); rownames(X) <- rownames(s$X)
  y <- ceiling(s$time / tie_unit) * tie_unit
  Wt <- sweep(s$W, 2, sqrt(colSums(s$W^2)), "/")
  list(X = X, X_raw = s$X, y = y, d = s$status, W_true = Wt, beta_true = beta, linpred = s$linpred,
       n_tied = sum(duplicated(y)))
}
