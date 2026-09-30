# Tracks 1 and 2: shared contract and briefs (2026-09-30)

Governing plan: `docs/plans/2026-09-30-robust-rerun-plan.md`, rev 2.1, FROZEN at commit 960b76d on
branch `recompute/consensus-desurv`. Sections 2-3 are not open for redesign. Where this brief and the
plan disagree, the plan wins; report the conflict, do not improvise.

## Repos and worktrees (each agent works ONLY in its own)

| Track | Repo worktree | Branch | Base |
|---|---|---|---|
| 1 (estimator) | `~/Downloads/DeSurv-track1` | `feat/consensus-estimator` | 72d8495 (1.0.2 tie fix) |
| 2 (package BO part) | `~/Downloads/DeSurv-track2` | `feat/bo-consensus` | 72d8495 |
| 2 (pipeline part) | `~/Downloads/DeSurv-paper-recompute` | `recompute/consensus-desurv` | 960b76d |

Never touch `~/Downloads/DeSurv-clean` (fix branch), `~/Downloads/DeSurv-paper-clean` (the canonical
manuscript tree; its render depends on the site-library DeSurv 1.0.1), or the site/user R libraries.
Install only into a private library: `R CMD INSTALL --preclean --library=<track lib>`; deps live in
`~/R/x86_64-pc-linux-gnu-library/4.6` and `~/R/library` (put them after the private lib in R_LIBS).
After any install, verify `exists("cox_stats_cpp", asNamespace("DeSurv"))`.
Never run `code/NN_*.R`, `make` targets or renders (plain runs RECOMPUTE and OVERWRITE results).
No commits, pushes or tags: the coordinator commits. Machine is shared: at most 8 cores.

## Interface contract (Track 1 implements; Track 2 codes against it)

```r
desurv_consensus_fit(
  X, y, d, k, alpha, lambda, nu, lambdaW = 0, lambdaH = 0,
  seeds,                      # integer vector, one independent random-start run per seed (prefix rule)
  tol, tol_init, maxit, imaxit,   # passed to each constituent desurv_fit (code/04 seed-fit settings)
  ncores = 1L,
  outlier_L_frac = 0.3, outlier_nmad = 3,        # section 3 outlier rule (frozen defaults)
  mode_min_frac = 0.20, mode_r = 0.95,           # section 3 competing-mode rule (frozen defaults)
  run_fits = NULL,            # optional precomputed constituent fits (skip fitting; used by the gate)
  verbose = FALSE
)
```

Returns an object of class `c("desurv_consensus_fit", "desurv_fit")` that works unchanged with
`predict()`, `coef()` and the top-gene scoring path (`.desurv_lp_with_top_genes`, `desurv_get_top_genes`),
with W = consensus W (L2-normalized columns, as the scoring path expects), H and beta from the fixed-W
refit, plus `$consensus`:
- `seeds`, `run_objective`, `run_cindex` (recorded, never used for selection)
- `family` (R x k integer: family of each run's factor), `occupancy` (R x k counts of 0/1/>1)
- `tightness` (per family median pairwise r), `n_outliers` (per family)
- `competing_mode` (per family logical + sub-cluster sizes and medoid r), `mode_pass` (all families)
- `settings` (all rule parameters used)

Algorithm = plan section 2, Track 1, steps 1-8, in that order. Behaviour when a competing mode is found:
the family's members are NOT outlier-filtered; the median is taken over all its members; `mode_pass` is
FALSE. The object is still returned (callers decide; the gate treats it as a fail). This follows the
frozen rule "a competing mode is never filtered away" and is not a design change.

```r
desurv_consensus_gate(
  run_fits, X, y, d, R_values, n_splits = 10L, split_seed,
  median_r = 0.98, p05_r = 0.95, ...   # section 3 half-versus-half thresholds (frozen defaults)
)
```

Returns, per R in `R_values`: per-factor median and 5th-percentile half-versus-half r over the splits,
`mode_pass`, pass/fail, plus descriptive (non-gating) beta and training C-index agreement.
Half-splits: `set.seed(split_seed)`, then `n_splits` successive `sample(R)`; halves = first and second
R/2. Factor alignment across halves by the Hungarian method on correlation.

## Track 1 brief (package: the estimator)

1. Fixed-W refit: prefer a flag on the existing C++ optimizer (`optimize_loss_cpp`) that skips the W
   proposal and W normalization, so H and beta use exactly the same update code as a normal fit. If the
   code makes that unsafe (e.g. zeta or normalization couples W and beta), report and propose the
   alternative before building it.
2. `desurv_consensus_fit` and `desurv_consensus_gate` per the contract. Clustering of the pooled R x k
   L2-normalized factors into k families by correlation distance (average-linkage hierarchical or PAM;
   state which and why in the report; it must be deterministic). Competing-mode check with PAM(2)
   (`cluster::pam`, deterministic build/swap). Median, renormalize, fixed-W refit, then top genes with the
   deterministic tie-break (score, then gene name); also fix `desurv_get_top_genes` / consensus_seed tie
   handling if they break ties by row index (FP15).
3. `init()` best-by-training-C selection (pkg R/init.R:180): do not delete; leave as is and report its
   callers. The consensus path must not call it for selection.
4. Tests (testthat, tiny synthetic data, fast): the plan's section 4 tests 1-7. Each test must fail if
   the property it guards is broken (show this for at least tests 2, 3, 6 by a temporary mutation, then
   revert).
5. Version 1.1.0 in DESCRIPTION, NEWS entry. Exports and roxygen docs for the two new functions.
6. Report: design choices, file:line of changes, test output before/after, full-suite result,
   `R CMD check --no-manual` NOTE/WARNING delta, anything in the contract you could not honour.

## Track 2 brief (tuning correctness)

Package part (`DeSurv-track2`):
1. CV scoring with the consensus estimator: an option on `desurv_cv` (and through `desurv_cv_bayesopt`)
   `estimator = c("mean_starts", "consensus")`, default unchanged; "consensus" fits
   `desurv_consensus_fit` in each training fold with a per-fold seed bank and scores the held-out fold
   through the same top-gene scoring path with `ntop`. Code against the contract with a local stub; mark
   every stub use so the coordinator can swap in Track 1's implementation.
2. Repeated CV: a list of fold-assignment seeds (the plan: 20261101-3), identical for every
   configuration; per-fold run banks identical across configurations (paired).
3. BO incumbent by GP posterior mean (not max observed).
4. Finalists: a function that takes the BO history, picks the top 2-3 configurations per k by posterior
   mean, re-evaluates them on the identical repeated CV, and returns per-configuration fold-level
   C-indices, so the 1-SE rule uses the empirical paired SE. k selection: the smallest k whose best
   finalist is within 1 empirical SE of the overall best; report the runner-up k.
Paper part (`DeSurv-paper-recompute`):
5. code/03:117-123 and 160-167: all hyperparameters from the selected k's own row
   (`select_bo_params_1se`, R/bo_helpers.R:152-208, or its replacement from item 4). Wire the new options,
   with quick-mode branches for every new knob (repo quick-mode contract). Do not run the script.
6. code/06:49-51 and code/19 (ntop 270) read `tar_params_best`, no typed values.
7. Tests (create `tests/` in the paper repo, testthat, runnable by `Rscript -e 'testthat::test_dir("tests")'`
   without touching results/): plan section 4 tests 8 (BO estimator identity) and 9 (selected-k
   integrity on a toy BO table with distinct values per k).
8. Report: file:line of every change, how item 1 uses the stub, test output, and every place where the
   old estimator path is still reachable.
