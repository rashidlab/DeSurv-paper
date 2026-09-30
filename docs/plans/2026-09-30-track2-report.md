# Track 2 report: tuning correctness (2026-09-30)

Brief: `docs/plans/2026-09-30-tracks-1-2-briefs.md` ("Track 2 brief"). Plan: rev 2.1, sections 2-4 (frozen).
Nothing committed. Nothing in `results/` read or written. No `code/NN_*.R`, `make` or render was run.
Private library: `/home/naimrashid/research/desurv-libs/track2` (DeSurv 1.0.2 + Track 2 changes;
`exists("cox_stats_cpp", asNamespace("DeSurv"))` is TRUE after the final install).

## 1. What changed

### Package (`~/Downloads/DeSurv-track2`, branch `feat/bo-consensus`, uncommitted)

| File | Change |
|---|---|
| `R/cv_consensus.R` (new, 346 lines) | `desurv_cv_seed_bank()` (L30, exported); `.desurv_check_consensus_control()` (L56); `.desurv_cv_summarise_repeated()` (L78); `.desurv_cv_grid_consensus()` (L133), the consensus CV engine; the call to `desurv_consensus_fit` is at L279 |
| `R/desurv_cv.R` | roxygen for 6 new args (L172-197); signature `estimator`, `cv_seeds`, `consensus_runs`, `consensus_seed_base`, `consensus_seed_stride`, `consensus_control` (L257-263); argument checks (L269-281); dispatch: consensus engine, unchanged single-partition path, repeated `mean_starts` path (L309-403) |
| `R/desurv_cv_bayesopt.R` | args `estimator`, `incumbent = c("observed", "posterior_mean")` (L178-187); `estimator` passed to every `desurv_cv` call (L229); EI incumbent value = max GP latent mean when `posterior_mean` (L451-459); final incumbent by GP posterior mean, `history$gp_mean` (L499-527); new outputs `best$gp_mean`, `best$eval_id`, `km_fit_final`, `incumbent`, `design_units`, `cv_args`, `data_fingerprint` (L558-571); docs L48-60, L89-98 |
| `R/desurv_bo_refine.R` | L360-366: with `posterior_mean`, `overall_best` = last run's incumbent (its GP covers the imported history) instead of max observed across runs |
| `R/desurv_bo_finalists.R` (new, 302 lines) | `.desurv_bo_gp_latent_mean()` (L14), `.desurv_bo_gp_posterior_at_design()` (L28), `desurv_bo_finalists()` (L115, exported), `desurv_select_k_paired_1se()` (L230, exported), `print.desurv_bo_finalists` |
| `R/zzz_TRACK1_STUB_desurv_consensus_fit.R` (new) | TRACK1-STUB, see section 2 |
| `NAMESPACE` | exports `desurv_bo_finalists`, `desurv_cv_seed_bank`, `desurv_select_k_paired_1se`; `S3method(print, desurv_bo_finalists)` (hand-edited; identical to what roxygen generates) |
| `man/` | regenerated `desurv_cv.Rd`, `desurv_cv_bayesopt.Rd`; new `desurv_bo_finalists.Rd`, `desurv_cv_seed_bank.Rd`, `desurv_select_k_paired_1se.Rd` (roxygen 8.0.0 run in a scratch copy; only these 5 files copied, to avoid RoxygenNote churn in unrelated Rd files) |
| `tests/testthat/test-cv-consensus.R`, `test-bo-finalists.R` (new) | see section 3 |

Item-by-item:
1. **CV with the consensus estimator.** `desurv_cv(estimator = "consensus", cv_only = TRUE, consensus_runs = R, consensus_control = list(tol_init, imaxit))` fits `desurv_consensus_fit()` once per training fold (per partition, hyperparameter row and alpha) and scores the held-out fold with `.desurv_lp_with_top_genes(fit, X_val, ntop)`, the same top-gene scoring path as before. `tol`/`maxit` come from `desurv_cv`'s own args. `consensus_control` accepts only `tol_init`, `imaxit`, `ncores`; the section 3 rule parameters are refused ("frozen"). Default `estimator = "mean_starts"` is the untouched historical path (regression-tested). `desurv_cv_bayesopt(estimator = ...)` passes it through; the other consensus args pass through `...`/`bo_fixed` (also through `desurv_cv_bayesopt_refine`).
2. **Repeated CV, paired.** `cv_seeds` = one stratified partition per seed (identical for every configuration). Per-fold seed bank `desurv_cv_seed_bank(rep, fold, nfolds, R)` = `4,000,000 + ((rep-1)*nfolds + (fold-1))*10,000 + (1..R)`: first seed 4,000,001, prefix rule, independent of every hyperparameter, disjoint blocks (R <= 10,000 enforced). Summary: `mean_cindex` over all partition x fold units; `se_cindex` = sd(units)/sqrt(n); `se_cindex_rep` = SE of partition means; `n_units`. Works for both estimators.
3. **Incumbent by GP posterior mean.** Found that `DiceKriging::predict.km` interpolates the noisy observations at design points even with an estimated nugget (verified: max |y - mean| = 1.7e-15), so "posterior mean at evaluated points" via `predict.km` would just be the observed max again. `.desurv_bo_gp_latent_mean()` computes the latent (noise-free) mean with the nugget excluded from the cross-covariance; it equals `predict.km` off the design (tested to 1e-8). The final GP is refitted on all successful evaluations under the BO seed (RNG restored). Used both for the EI incumbent value and for `best`.
4. **Finalists and k.** `desurv_bo_finalists(bo, X, y, d, dataset, samp_keeps, n_per_k)`: top `n_per_k` per k by `gp_mean` (ties: observed mean, then eval_id), re-evaluated with `desurv_cv` under exactly the BO's stored `cv_args` (only `parallel_grid`, `ncores_grid`, `verbose` may be overridden; data checked against a fingerprint). Returns per-finalist fold-level C-indices (`rep`, `fold`) and `desurv_select_k_paired_1se()`: per k, its best finalist vs the overall best on paired units, `se_diff = sd(C_best - C_k)/sqrt(n_pairs)`; k within 1 SE if `mean_diff <= se_diff`; selected k = smallest such k; runner-up k = the other k with the highest best-finalist mean.

### Pipeline (`~/Downloads/DeSurv-paper-recompute`, branch `recompute/consensus-desurv`, uncommitted)

| File | Change |
|---|---|
| `R/consensus_settings.R` (new) | one source of the constituent-run settings = code/04 seed-fit settings (`tol = RUN_TOL/100`, `tol_init = RUN_TOL`, `maxit = imaxit = RUN_MAXIT`; quick and full); `desurv_estimator_entry()` (namespace binding the CV folds call); `fit_consensus_desurv()`, the final-fit entry point for code/04 (not yet wired, see section 5) |
| `R/bo_helpers.R` | `select_bo_params_finalists()` (L216): every hyperparameter from the selected finalist's own row, stops if that row's k differs from the selected k; `k_selection_from_finalists()` (L247): `tar_k_selection`-shaped list; `compute_bo_eval_se()` groups by `rep` when present (L351-354) |
| `code/03_bayesian_optimization.R` | L14-18 header; L27 sources `R/consensus_settings.R`; L29-40 knobs `BO_ESTIMATOR = "consensus"` (L36), `BO_INCUMBENT = "posterior_mean"`, `CONSENSUS`, `CV_SEED_BANK_BASE = 4000000L`; quick/full branches for every new knob (L41-72: `BO_CV_SEEDS` quick 20261101 / full 20261101-3, `BO_CONSENSUS_R` quick 3 / full from env `DESURV_BO_CONSENSUS_R`, which stops if unset, `FINALISTS_PER_K` 2 / 3); `BO_TOL`/`BO_MAXIT` = `CONSENSUS$tol`/`$maxit` (L74-75); `n_starts` only for `mean_starts` (L89); `BO_COMMON` gains estimator/incumbent/cv_seeds/consensus args (L115-121); `run_bo_finalists()` (L125) and `bo_ngene()` (L139); main BO: finalists, `tar_k_selection`, `tar_params_best` from the selected finalist row + `stopifnot(k)` (L167-186); alpha = 0 control: same finalists selection, alpha forced 0 (L210-224); `ngene` fallbacks read `BO_FIXED$ngene` via `bo_ngene()` (L185, L224, L270) and L322 |
| `code/06_cv_grid.R` | L5-6 header; L50-56: `lambda`, `nu`, `ntop` from `load_precomputed("tar_params_best_tcgacptac")` (was typed 0.349 / 0.056 / 270) |
| `code/19_desurv_vs_supervised_tuned.R` | L15-16 header; L45-53: `NTOP` from `results/tar_params_best_tcgacptac.rds`, and a hard stop if k != 3 (the D1/D2/D3 layout is structural); L77-78 `get_top_genes(W, NTOP)`; L151 method label built with sprintf (same string when ntop 270, k 3) |
| `tests/helper-paths.R`, `tests/test-bo-estimator-identity.R`, `tests/test-selected-k-integrity.R` (new) | plan tests 8 and 9 |

## 2. How item 1 uses the stub

`R/zzz_TRACK1_STUB_desurv_consensus_fit.R` defines `desurv_consensus_fit` with the contract's exact signature, only if no other file in the namespace already defined it (`if (!exists(..., envir = topenv(), inherits = FALSE))`; the file sorts last, there is no Collate field). Not exported. It fits one `desurv_fit(ninit = 1, seed = s)` per seed (or uses `run_fits`), aligns factors greedily to run 1, takes the gene-wise median, L2-normalizes, keeps H and beta from run 1 (no fixed-W refit), and returns class `c("desurv_consensus_fit", "desurv_fit")` with the contract's `$consensus` field names plus `stub = TRUE`. Every line is marked TRACK1-STUB.

Stub uses (all call the plain name, so no call site changes when the stub is deleted):
- `R/cv_consensus.R:279`, the only package call site.
- Paper repo `R/consensus_settings.R` `desurv_estimator_entry()` looks the binding up with `get(..., asNamespace("DeSurv"))` because the stub is unexported; once Track 1's exported function is merged it is the same object as `DeSurv::desurv_consensus_fit` (the paper test asserts this identity).
- Tests that run the estimator for real (package `test-cv-consensus.R` tests 5-6, paper test 8b) run the stub and check wiring only. Diagnostics carry `estimator_stub` so a stubbed CV result can never be mistaken for a real one.

Swap-in: delete the stub file, merge Track 1, reinstall, rerun both test suites. Consensus diagnostics read `fit$cindex`, `fit$convergence`, `fit$nan_flag`, `fit$consensus$mode_pass`, `fit$consensus$n_outliers` with NA fallbacks (the first three are not in the contract).

## 3. Tests

Package (`NOT_CRAN=true R_LIBS=<track2 lib>:~/R/x86_64-pc-linux-gnu-library/4.6:~/R/library Rscript -e 'testthat::test_dir("tests/testthat", package = "DeSurv", load_package = "installed")'`):
- Before (existing 14 files): 107 test blocks, 0 failed, 0 skipped.
- After: 16 files, 118 test blocks, 384 expectations, **0 failed, 0 errors, 0 skipped**. New: `test-cv-consensus.R` (6 blocks, 36 expectations: seed bank origin/prefix/disjointness; default path identical to `.desurv_cv_grid_warmstart` output; repeated-CV pooling; one `desurv_consensus_fit` per partition x fold with identical seed banks across two different configurations and the right settings; CV scores exactly the returned consensus object through the top-ntop rule; determinism and guards) and `test-bo-finalists.R` (5 blocks, 35 expectations: latent mean vs predict.km; posterior-mean incumbent; observed default unchanged; finalists pick top-n per k by gp_mean and re-evaluate with identical CV args, reproducing BO values exactly; data-fingerprint and override guards; paired 1-SE selection and runner-up).
- Mutation check: making the seed bank depend on (k, alpha) fails the pairing test (2 failures); reverted and reinstalled.
- `R CMD check --no-manual --ignore-vignettes`: base (HEAD 72d8495) 1 WARNING, 3 NOTEs; after: the same 1 WARNING, 3 NOTEs, with identical text for the R-code NOTE and the codoc WARNING. **No delta.** (Check runs no tests: the package has no `tests/testthat.R`.)

Paper repo (`R_LIBS=<track2 lib first> Rscript -e 'testthat::test_dir("tests")'` from the repo root; reads and writes nothing in `results/`, verified by `find results -newer`):
- 2 files, 5 test blocks, 57 expectations, **0 failed, 0 skipped**.
- Test 8 (`test-bo-estimator-identity.R`): (a) static parse of code/03: `BO_ESTIMATOR == "consensus"`, `BO_INCUMBENT == "posterior_mean"`, `BO_TOL`/`BO_MAXIT`/`consensus_control` all from `desurv_consensus_settings()`, all four BOs go through `BO_COMMON`, full-mode fold seeds 20261101-3; (b) dynamic: `desurv_cv(estimator = "consensus")` fold fits and `fit_consensus_desurv()` hit the same mocked namespace binding with equal (k, alpha, lambda, nu, lambdaW, lambdaH, tol, tol_init, maxit, imaxit), same class, and CV scores exactly the returned objects. Mutation `BO_ESTIMATOR <- "mean_starts"`: fails; reverted.
- Test 9 (`test-selected-k-integrity.R`): toy BO table with distinct values per k; `select_bo_params_1se` returns the whole k = 3 row and the pre-rev-2 recipe is shown to fail; `select_bo_params_finalists` returns the whole selected finalist row (runner-up 7, `lcb_threshold` NA) and refuses a row whose k disagrees; `compute_bo_eval_se` treats partition x fold as units. Mutations (take the overall best row): `select_bo_params_1se` 4 failures, `select_bo_params_finalists` errors on the k guard; reverted.
- With only the site DeSurv 1.0.1 on the path, the dynamic parts (8b, 9b) skip.

## 4. Where the old estimator path is still reachable

- Package defaults: `desurv_cv(estimator = "mean_starts")`, `desurv_cv_bayesopt(estimator = "mean_starts", incumbent = "observed")`, and `desurv_cv_bayesopt_refine` inheriting them. Kept for backward compatibility (brief: "default unchanged"); code/03 sets the new values explicitly.
- `desurv_cv(cv_only = FALSE)` refit (`mean_starts` only) still calls `desurv_fit()` with default `ninit = 20`, i.e. `init()` best-by-training-C (pkg `R/init.R:180`, Track 1's item 3).
- code/03 knobs `BO_ESTIMATOR <- "mean_starts"` / `BO_INCUMBENT <- "observed"` restore the old evaluator.
- code/04 final fits: still `desurv_consensus_seed()` + one `desurv_fit` (code/04:104-118). `fit_consensus_desurv()` is ready; rewiring code/04 is outside this brief (plan section 5).
- code/29_final_config_cv.R: re-runs the old BO evaluator (n_starts 30, tol 1e-5, seed 123) against the cached history.
- code/06 + `R/cv_grid_helpers.R`: its own estimator (30 short starts, best by training C via `init()`; `SEED = 123` single partition; `NGENE = 3000L` typed). Only lambda/nu/ntop were un-typed per the brief; rerun/re-scope/drop is owner decision 3.
- code/08:107 (`n_starts = 30`, seed 125), simulations (`R/sim_helpers.R:1300-1345`, `R/simulation_functions/select_bo_k_lcb.R`): old BO/CV path.
- `select_bo_k_by_cv_se()` and `select_bo_params_1se()` remain in `R/bo_helpers.R`, now unused by code/03.

## 5. Mismatches, ambiguities and decisions for the coordinator

1. **Seed-bank layout.** The plan fixes the origin (4,000,001) and "identical for every configuration" but not the layout across partitions and folds. Chosen: block per (partition, fold), width 10,000, prefix rule within a block. Max seed used at 3 x 5 blocks is 4,150,000, disjoint from all other banks.
2. **"Empirical paired SE".** Implemented as the SE of per-unit paired differences over all 15 partition x fold units. Units from different partitions share training patients, so this likely understates the SE (the rule becomes stricter, favouring larger k). The per-partition SE (`se_cindex_rep`, 3 values) is also returned but unused. A ruling may be wanted; I did not change the definition.
3. **Runner-up k** is defined as the other k whose best finalist has the highest CV mean.
4. **Finalist re-evaluation on identical CV** reproduces the BO's own numbers exactly (the evaluation is deterministic). It contributes the fold-level values and a determinism check (`delta_vs_bo`, 0 in tests), not new information. If the plan intended more partitions or a larger R for finalists, that is not written; not implemented.
5. **BO constituent tolerances** are now the final fit's (tol 1e-7, tol_init 1e-5, maxit = imaxit = 4000) per the contract, replacing tol 1e-5 / maxit 4000 with 30 averaged starts. Cost per BO evaluation rises with R_BO x 15 fold fits; the plan's runtime benchmark before the BO applies.
6. **R_BO** has no default in full mode: code/03 stops unless `DESURV_BO_CONSENSUS_R` is set (plan: the smallest R that passed Track 3a).
7. **Controls.** alpha = 0 (k tuned) uses the same finalists + paired 1-SE selection (brief item 5 allows `select_bo_params_1se` or the item-4 replacement; I used the replacement for symmetry). elbow-k and a0k3 (k fixed) use the posterior-mean incumbent without finalists. All four BOs use the consensus estimator via `BO_COMMON`.
8. **Manuscript fields.** `tar_k_selection` keeps `k_selected`, `k_best`, `best_mean`, `best_se` (best finalist's unit SE) and adds `runner_up_k`, `per_k`; `lcb_threshold` is NA because the paired rule has one threshold per k. Methods and SI read `.ks$lcb_threshold` (SI 4 sites, Methods 1): text work for plan section 6.
9. **code/19** hard-stops if k != 3 (its D1/D2/D3 decomposition and "D2 (proCAF factor)" label are structural), matching the plan's "if k != 3 ... rulings before anything downstream".
10. **Quick-mode ngene.** Previously quick BO used ngene 500 but the filtered data fell back to a typed 3000; `bo_ngene()` now uses `BO_FIXED$ngene` (quick 500). Full mode unchanged (3000).
11. **Version and NEWS** not touched (Track 1 owns 1.1.0 and the NEWS entry; avoids a merge conflict). Suggested NEWS lines: "desurv_cv(): `estimator = "consensus"` scores desurv_consensus_fit() in each training fold with a per-fold seed bank (desurv_cv_seed_bank()); `cv_seeds` for repeated, paired CV. desurv_cv_bayesopt(): `incumbent = "posterior_mean"` (GP latent mean). New desurv_bo_finalists() and desurv_select_k_paired_1se()."
12. Consensus CV parallelism: outer `mclapply` over (partition, fold, config) jobs with `ncores_grid`; spare cores go to the constituent runs (`floor(ncores_grid / outer)`), unless `consensus_control$ncores` is given.
13. The R CMD check was run with `--library` pointing at the track2 lib, so it reinstalled into that lib; I reinstalled from the worktree afterwards and re-verified `cox_stats_cpp` and the new exports.

Nothing in the frozen sections 2-3 was found impossible or contradictory in the code.
