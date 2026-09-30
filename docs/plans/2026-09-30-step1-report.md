# Rev 3 Step 1 report: "reported" fold estimator, tie-break, tests (2026-09-30)

Brief: `docs/plans/2026-09-30-step1-brief.md`. Authority: plan rev 3 sections 1, 2, 8 (step 1).
Status: DONE. All inventory items pass. No commits made in either repo; nothing staged.
Private library: `/tmp/claude-1000/-home-naimrashid-Downloads-DeSurv-paper-clean/da9c6c81-47b4-4511-8cfc-1c71451b090a/scratchpad/lib-rev3`
(final install: DeSurv 1.2.0, `cox_stats_cpp` present, `find.package` resolves to that library).

## Files changed

### Package `~/Downloads/DeSurv-clean` (branch `feat/reported-estimator`, uncommitted on 200a891)

| File | Change |
|---|---|
| `R/reported_fit.R` (new) | `desurv_reported_fit(X, y, d, k, alpha, lambda, nu, lambdaW = 0, lambdaH = 0, ntop, seeds, tol, tol_init, maxit, imaxit, ncores = 1L)`, exported. code/04 `run_seed_fits` + `consensus_and_fit` recipe; returns the final `desurv_fit` plus `$reported` (seeds, seeds_ok, n_seed_fits, n_successful, seed_cindex, min_frequency, ntop). |
| `R/cv_consensus.R` -> `R/cv_reported.R` | Renamed (plain `mv`). Engine `.desurv_cv_grid_reported` calls `desurv_reported_fit(ntop = ntop, seeds = bank, ...)`; args `consensus_runs` -> `n_seed_fits`, `consensus_seed_base/stride` -> `seed_bank_base/stride`, `consensus_control` -> `reported_control` (tol_init, imaxit, ncores; validation kept). Diagnostics: `n_runs` -> `n_seed_fits`, added `n_successful`; dropped `mode_pass`, `n_outliers`, `estimator_stub`. Stops if `ntop` is NULL. Seed bank function and layout unchanged. |
| `R/zzz_TRACK1_STUB_desurv_consensus_fit.R` | Deleted. |
| `R/desurv_cv.R` | `estimator = c("mean_starts", "reported")`, renamed args, `cv_obj$reported_settings` (entry_point `"desurv_reported_fit"`), docs. `cv_only = TRUE` requirement kept. |
| `R/desurv_cv_bayesopt.R` | `estimator = c("mean_starts", "reported")`, docs. Default incumbent `observed` unchanged. |
| `R/get_top_genes.R` | Line 91: `order(-diff_vector, gene_names, seq_along(diff_vector), method = "radix")` (score desc, name asc in C locale, then row index); doc sentence added. |
| `DESCRIPTION` | Version 1.0.2 -> 1.2.0 (nothing else). |
| `NEWS.md` | 1.2.0 entry: reported estimator, repeated paired CV + seed bank, posterior_mean incumbent, finalists, paired 1-SE k, tie-break. |
| `NAMESPACE`, `man/desurv_cv.Rd`, `man/desurv_cv_bayesopt.Rd`, `man/desurv_cv_seed_bank.Rd`, `man/desurv_get_top_genes.Rd`, `man/desurv_reported_fit.Rd` (new) | Regenerated with roxygen2 (man/ is roxygen-generated). |
| `tests/testthat/test-cv-consensus.R` -> `test-cv-reported.R` | Rewritten for "reported" (6 blocks; seed-bank and default-path blocks kept verbatim). |
| `tests/testthat/test-reported-fit.R` (new) | Numeric identity vs verbatim code/04 copy (serial and NCORES = 2), failure dropping, exact desurv_fit argument sets. |
| `tests/testthat/test-top-genes-tiebreak.R` (new) | Tie at ntop boundary under 4 row permutations (+ scoring lp invariance); untied 3000 x 3 W equals `order(diff, decreasing = TRUE)` for ntop 1/50/300/3000. |

### Paper `~/Downloads/DeSurv-paper-recompute` (branch `recompute/consensus-desurv`, uncommitted on 883b39e)

| File | Change |
|---|---|
| `R/consensus_settings.R` -> `R/reported_settings.R` | Renamed (trivial: 3 references). `desurv_reported_settings(quick)` = code/04 values incl. `n_seed_fits` (100 full, 3 quick); `desurv_estimator_entry()` returns `DeSurv::desurv_reported_fit` (no `get(asNamespace)`); `fit_reported_desurv()` passes `ntop` by code/04's `ntop_value` rule (round, fallback 100). |
| `code/03_bayesian_optimization.R` | `BO_ESTIMATOR <- "reported"`; `REPORTED <- desurv_reported_settings(CONFIG$quick)`; `BO_N_SEED_FITS <- REPORTED$n_seed_fits`; `BO_COMMON` passes `n_seed_fits`, `seed_bank_base`, `reported_control`; `DESURV_BO_CONSENSUS_R` stop and `BO_CONSENSUS_R` removed. Incumbent, fold seeds, finalists, 1-SE, quick branches untouched. code/04, 06, 19 untouched. |
| `tests/helper-paths.R` | Comment path updated. |
| `tests/test-bo-estimator-identity.R` | Test 8 re-pointed: static asserts `BO_ESTIMATOR == "reported"` and that settings equal code/04's `NINIT_FULL`, `RUN_TOL/100`, `RUN_TOL`, `RUN_MAXIT` in both modes (parsed from code/04's `if (CONFIG$quick)` block); dynamic asserts CV folds and `fit_reported_desurv()` call the same exported `DeSurv::desurv_reported_fit` object with identical args, and CV scores exactly the returned objects. Test 9 unchanged, passes. |

## Breakage inventory

1. **Recipe drift: none found.** `test-reported-fit.R`: `desurv_reported_fit()` vs verbatim inline copy of code/04:70-118 (seeds 1:5, RUN_TOL 1e-4, RUN_MAXIT 60, ntop 6) is `all.equal` on W, H, beta, cindex, the whole fit object minus `$reported`, kept seeds, seed C-indices and min_frequency; same again with NCORES = 2 (mclapply). Argument-name test: seed fits get exactly code/04's 17 args (incl. `ninit = 1`, `parallel_init = FALSE`); the final fit gets W0/H0/beta0, `seed = NULL`, tol, maxit, verbose and no tol_init/imaxit/ninit. Failure test: seeds that error are dropped (seeds_ok 1,3,5; min_frequency = ceiling(0.3 x 3) = 1); all failing -> error.
2. **Pairing: intact.** Mutation check in a scratch copy and scratch library: seed-bank base made to depend on k (`+ (hp$k - 2) * 100000`); `test-cv-reported.R:78` (`expect_identical(seeds1, seeds2)`) FAILED as required, the only failure. Mutant copy and library deleted; real source never mutated; $LIB reinstalled from the real tree.
3. **Tie-break on untied data: unchanged.** Old (200a891 build) vs new `desurv_get_top_genes` on 20 random 3000 x 3 W x ntop {10, 100, 500, 3000}: identical outputs 80/80. Tied W: old gives `gA gD` / `gA gC` depending on row order; new gives `gA gB` for both. Single ordering site confirmed: grep of R/ and src/ finds no other top-gene ordering; `desurv_consensus_seed` reaches it via `R/consensus_seed.R:227`.
4. **Tie fix intact.** `cox_stats_cpp` exists in the installed 1.2.0; test-cox-ties 5 blocks pass; `git diff 72d8495 -- src/` = 0 lines.
5. **Consensus residue.** Package grep (`R man NAMESPACE tests`): 0 lines. Paper grep (`code R tests`): 3 lines, all in `code/21_consensus_endpoint.R` (86, 108, 109): a local helper named `consensus_fit` in the existing program-stability analysis, unrelated to the removed estimator. Kept deliberately (brief says do not touch other scripts). No `BO_CONSENSUS_R`, `consensus_runs`, `consensus_control` or `consensus_settings` anywhere outside docs/.
6. **Default path unchanged.** "default estimator path is unchanged (single partition)" passes (block kept verbatim). Default `estimator = "mean_starts"`, `incumbent = "observed"`.
7. **Counts** (files / test blocks / expectations / failures / errors / skips; warnings are the pre-existing non-convergence warnings in test-top-genes and others):

   | Suite | Before (200a891 as installed) | After |
   |---|---|---|
   | Package | 16 / 118 / 384 / 0 / 0 / 0 (8 warnings) | 18 / 124 / 458 / 0 / 0 / 0 (8 warnings) |
   | Paper | 2 / 5 / 57 / 0 / 0 / 0 | 2 / 5 / 68 / 0 / 0 / 0 |

   `find results -newer <marker>` in the paper worktree: empty for both paper runs.
8. **R CMD check** (`--no-manual --ignore-vignettes`, throwaway `--library`, built tarballs): 72d8495 `1 WARNING, 3 NOTEs`; new `1 WARNING, 3 NOTEs`. Item text (codoc mismatch; hidden `.zenodo.json`; License stub; `setNames`/`best_so_far` globals) diffed after stripping timings: identical. Reinstalled into $LIB afterwards.

## Mismatches with the brief

None where the brief conflicts with the code. Decisions where the brief was silent:

- **ntop in CV.** `desurv_consensus_seed` needs a positive `ntop`, but `desurv_cv` normalises `ntop <= 0` or NULL to NULL. For `estimator = "reported"` the engine now stops if ntop is NULL (code/03 always tunes ntop, so this does not bind). `fit_reported_desurv()` follows code/04's fallback (100 if params lack ntop).
- **Worker crash.** If an mclapply worker dies (not an R error inside `desurv_fit`, which `try` catches), `desurv_reported_fit` stops with the worker's message. code/04 would also stop, with an obscure `[[` error. Failure dropping is otherwise the same as code/04.
- **roxygen version.** Installed roxygen2 is 8.0.0; RoxygenNote is 7.3.3. Roxygen 8 also rewrote `DESCRIPTION` (RoxygenNote -> Config/roxygen2/version), `man/DeSurv-package.Rd` and `man/preprocess_data.Rd` (pre-existing drift, likely the codoc WARNING). Those three changes were reverted so the only DESCRIPTION change is the version. Only the Rd files for edited functions changed.
- **Return object.** `desurv_reported_fit` returns the final `desurv_fit` (class unchanged) with an added `$reported` list; the identity test compares every other element to code/04's fit.
- **File renames** (`cv_consensus.R`, `test-cv-consensus.R`, paper `consensus_settings.R`) are plain renames on disk; the coordinator should `git add -A` so git sees delete + new.
- Build artefacts `src/*.o`, `src/DeSurv.so` reappear after `R CMD INSTALL` in the source tree; they are git-ignored (as at 72d8495).
