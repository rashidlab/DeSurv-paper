# Qualification data generator for the estimator-redesign spec (docs/plans/2026-10-05-qualification-prespec.md,
# section 3). One code path for S1, S2, S3 and S3b: the repo generator (R/simulation_functions/, scenario R0) with the
# true H columns multiplied by cvec between simulate_H and X = W H'. The multiplication draws no random numbers, so
# W, marker sets, H draws, survival and noise draws are identical to the unscaled call with the same seed.
# Training and test samples come from ONE generator call (N = n_train + n_test, same W, gene baselines and score
# standardization) and are split by index. Simulated data only.
# Requires QUAL (path of analysis/estimator-gate/qualification) to be set by the caller.
local({
  sd <- file.path(QUAL, "..", "..", "..", "R", "simulation_functions")
  for (f in c("simulate_W.R", "simulate_H.R", "simulate_X.R", "simulate_survival.R",
              "scenario_defaults.R", "simulate_data.R")) source(file.path(sd, f), local = globalenv())
})
.simulate_H_orig <- simulate_H

QUAL_SCEN <- list(S1  = list(beta = c(2, 1, 0)), S2 = list(beta = c(0, 0, 0)),
                  S3  = list(beta = c(2, 0, 0), target = 0.25), S3b = list(beta = c(2, 0, 0), target = 0.50))

# Raw generator call with H columns scaled by cvec. Asserts the override ran exactly once and was restored.
gen_scaled <- function(seed, beta, cvec, N, p = 3000L) {
  calls <- 0L
  assign("simulate_H", function(...) { calls <<- calls + 1L; sweep(.simulate_H_orig(...), 2, cvec, "*") },
         envir = globalenv())
  on.exit(assign("simulate_H", .simulate_H_orig, envir = globalenv()))
  s <- suppressWarnings(simulate_desurv_scenario(scenario = "R0", G = p, N = N, K = 3L, beta = beta,
         markers_per_factor = 150L, B_size = 500L, seed = seed))
  stopifnot(calls == 1L)
  s
}

# Unit-norm marker-masked true loadings: the generator's own score definition (simulate_survival.R, Wtilde with
# survival_gene_n = NULL keeps exactly the marker genes), normalized for a scale-free SD ratio.
masked_loadings <- function(s) {
  Wm <- matrix(0, nrow(s$W), ncol(s$W), dimnames = dimnames(s$W))
  for (j in seq_len(ncol(s$W))) { g <- match(s$marker_sets[[j]], rownames(s$W)); Wm[g, j] <- s$W[g, j] }
  sweep(Wm, 2, sqrt(colSums(Wm^2)), "/")
}

# Realized SD ratio (factor 1 over factor 2) of masked scores on the RANKED training X: the binding S3 quantity.
sd_ratio_rank <- function(Xr_train, Wm) { sc <- crossprod(Xr_train, Wm); sd(sc[, 1]) / sd(sc[, 2]) }

sim_qual <- function(seed, scenario, cvec = c(1, 1, 1), n_train = 218L, n_test = 1000L, tie_unit = 1, check = TRUE) {
  sc <- QUAL_SCEN[[scenario]]; N <- n_train + n_test
  s <- gen_scaled(seed, sc$beta, cvec, N)
  if (check && any(cvec != 1)) {   # wrapper self-check: same stream, H scaled exactly, nothing else touched
    s1 <- gen_scaled(seed, sc$beta, c(1, 1, 1), N)
    stopifnot(identical(s$W, s1$W), identical(s$marker_sets, s1$marker_sets),
              isTRUE(all.equal(s$H, sweep(s1$H, 2, cvec, "*"), tolerance = 1e-15)))
  }
  Xr <- apply(s$X, 2, rank, ties.method = "average"); rownames(Xr) <- rownames(s$X)
  tr <- seq_len(n_train); te <- n_train + seq_len(n_test)
  y <- ceiling(s$time / tie_unit) * tie_unit
  Wm <- masked_loadings(s); Wt <- sweep(s$W, 2, sqrt(colSums(s$W^2)), "/")
  lat <- crossprod(s$W %*% t(s$H), Wm)   # generator-scale masked scores (pre-standardization), N x 3
  list(X = Xr[, tr], y = y[tr], d = s$status[tr],
       X_test = Xr[, te], y_test = y[te], d_test = s$status[te],
       linpred_test = s$linpred[te], scores_test = s$scores[te, , drop = FALSE],   # generator's standardized scores
       W_true = Wt, W_mask = Wm, marker_sets = s$marker_sets, beta_true = sc$beta, cvec = cvec, seed = seed,
       ratio_rank_train = sd_ratio_rank(Xr[, tr], Wm),
       ratio_latent_train = sd(lat[tr, 1]) / sd(lat[tr, 2]),
       ratio_full_rank_train = { z <- crossprod(Xr[, tr], Wt); sd(z[, 1]) / sd(z[, 2]) },
       cor_AB_test = c(generator = cor(s$scores[te, 1], s$scores[te, 2]),
                       masked_rank = cor(crossprod(Xr[, te], Wm))[1, 2]),
       n_tied = sum(duplicated(y[tr])), event_rate = mean(s$status[tr]))
}
