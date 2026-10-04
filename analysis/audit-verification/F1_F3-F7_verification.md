# Verification of audit findings F1, F3, F4, F5, F6, F7 (2026-10-04)

Audit: `docs/HANDOFF-independent-scientific-audit-2026-10-04.md` (DeSurv-paper-clean, branch revision/desurv-step6, HEAD fb731e7). F2 not examined (coordinator).

Read-only. `git status --short` in DeSurv-paper-clean was empty before and after; HEAD fb731e7 unchanged. (DeSurv-paper-recompute shows three pre-existing untracked files under `analysis/tie-tol-check/` and `analysis/whatif-d1-own-rank/`; not created by this session.) All scripts and compiled objects live in `scratchpad/sub/` (`cox101.cpp`, `f1.R`, `f1b.R`, `f5.R`, `f6.R`, `lp_sorted.rds`, `rcppcache/`).

Versions: site library DeSurv 1.0.1, RemoteSha afb00d5fdf26, RemoteRef submission/pnas-2026. DeSurv 1.2.0 from the private lib-rev3 library (source 1f3d3f0). The adjacent `../DeSurv` checkout (HEAD 370c88a, Version 1.0.1) has `src/functions.cpp` byte-identical to the extracted 1.0.1 source; afb00d5 is not an object in that repo.

## Summary table

| Finding | Verdict | Revised severity | Disclosed in manuscript | Remediation class |
|---|---|---|---|---|
| F1 tie derivative | CONFIRMED (defect); downstream CV-C impact measured negligible; W/beta impact unmeasured | MEDIUM (was HIGH) | SI line 111 discloses Breslow ties for the likelihood; not the derivative defect | code correction (done in 1.0.2/1.2.0) + deterministic verification of fit impact |
| F3 hybrid final config below threshold; CV estimator differs from consensus fit | CONFIRMED (all numbers exact) | LOW-MEDIUM (fully disclosed) | Yes, Methods line 86 and SI 811-879, Table final-config-cv | manuscript-only (already done); optional rerun |
| F4 grid preprocessing outside folds | CONFIRMED (gene counts exact) | LOW (sensitivity analysis only) | No | manuscript-only disclosure, or rerun grid fold-local |
| F5 superpc correlation surrogate | CONFIRMED (all C values exact) | LOW-MEDIUM | No; Methods line 119 says "the fitted model was frozen and applied", which is inaccurate | code correction (native predict) + manuscript number update |
| F6 zero-fill of absent assay genes | CONFIRMED as mechanism; impact NARROWED to negligible | LOW | No; Methods line 163 says "No missing-data imputation was performed" | manuscript-only |
| F7 provenance gaps | CONFIRMED for items checked | MEDIUM (reproducibility, not correctness) | Partly (SI table note discloses unrecorded BO version) | deterministic verification + manifest; manuscript fix for session-info promise |

## F1. Tied-time Cox derivative

Verdict: CONFIRMED for the defect. Impact narrowed.

Evidence. 1.0.1 `src/functions.cpp:cox_suffix_stats_full` lines 14-72. `h[i]` and `B[i]` (lines 47-56) accumulate `d[j]/S0` over j <= i in storage order, so a subject at the start of a tie group excludes later tied events whose risk sets contain it. `g_z = d - eeta % h` and `w = eeta % h - eeta^2 % B` (lines 62-63) inherit the error. The log-likelihood (lines 66-70) uses `logS0_full[gstart[i]]` and is correct. 1.2.0 (`~/Downloads/DeSurv-clean/src/functions.cpp` lines 36-67) adds a backward `gend` pass and reads `h[i] = hcum[gend[i]]`, `B[i] = Bcum[gend[i]]`, which is the Breslow derivative. Callers: lines 124, 195 (W step), 270, 328 (coordinate descent in `cdfit_cox_dh_one_lambda`), so the W gradient and the beta coordinate-descent weights both use it.

Experiment (`sub/f1.R`): compiled the 1.0.1 function unchanged via Rcpp::sourceCpp with a thin wrapper; 1.2.0 via `DeSurv:::cox_stats_cpp`; independent base-R Breslow log-likelihood with central finite differences.
- Audit toy (t = 1,1,2,3,4,5; d = 1,1,1,0,1,0; eta = 0): 1.0.1 gradient (0.8333, 0.6667, 0.4167, -0.5833, -0.0833, -1.0833), sum 0.1667. 1.2.0 and finite differences give first entry 0.6667, sum -4.4e-16. Matches the audit exactly.
- Random tied data (n = 40, 12 unique times, seed 42): max |gradient - FD| 0.508 (1.0.1) vs 7.2e-9 (1.2.0); max |diag weight - FD| 0.443 vs 1.8e-6; gradient sum (common-shift invariance, should be 0) 3.16 vs 1.2e-14; 26/40 subjects wrong in 1.0.1. Log-likelihood equal to independent value to 4e-14 in both.
- Untied data: 1.0.1 and 1.2.0 identical (difference 0), both match FD to 6e-9. So the defect is tie-only.
- Swapping two tied subjects changes per-subject gradients by 0.004 in 1.0.1, 0 in 1.2.0.
- Audit reproducer on installed 1.0.1 `update_beta_cpp`: 0, -1.74085 after swap, coxph Breslow -1.5275. Reproduced exactly.

Manuscript training data (`results/tar_fit_desurv_tcgacptac.rds$data`, identical times in `tar_data_filtered_tcgacptac.rds`): n 273, 139 events, 242 unique times, 29 distinct tied time values, 60 subjects at tied times, 30 events at tied times, 22 tie groups containing at least one event, 7 groups with two or more events. Matches the audit. At the archived fitted linear predictor (sorted by time), 1.0.1 vs 1.2.0 gradient differs for 18 subjects, max absolute difference 0.0295, relative L2 difference 0.56 percent; diagonal weight relative difference 0.52 percent; 1.0.1 gradient sum 0.190 (should be 0).

Downstream: recompute `docs/plans/2026-10-03-k3-ladder-report.md` second "## 8. Results" section (line 130) exists; S1 to S2 (1.0.1 to 1.2.0, all else fixed) moves CV C 0.6150 to 0.6147, 0.6123 to 0.6122, 0.6317 to 0.6320, i.e. -0.0003, -0.0001, +0.0003 as stated in the overlap note; owner reading 2 there: "Tie fix ~0.000".

Proven vs conjectured: the derivative defect is proven. A material effect on CV C is refuted by the ladder. Effect on the archived W, beta, programs and support is not measured; given a 0.5 percent gradient perturbation it is likely small but not shown.

Disclosure: the manuscript describes the Breslow likelihood with ties; it does not mention the derivative defect. Revised severity MEDIUM (implementation defect, measured CV impact negligible). Smallest decisive check: refit the published configuration (consensus init, same seeds) under 1.2.0 with everything else fixed and report W loading correlation, top-270 support Jaccard, and beta, versus the archived fit. Remediation: code correction (already in 1.0.2/1.2.0) plus that deterministic comparison.

## F3. Final hyperparameters outside recomputed one-SE eligibility; CV estimator differs from final estimator

Verdict: CONFIRMED; already disclosed.

Evidence. `code/03_bayesian_optimization.R:112-123`: `tar_k_selection` from `select_bo_k_by_cv_se()`, then `tar_params_best` takes `desurv_bo_results$overall_best$params` and overwrites only `k`. Matches amendment 1c section 1 (lines 10-23): published configuration is the eval-129 (k 7) optimum run at k 3; eval 136 admitted k 3.

Artifacts (full precision): `tar_k_selection_tcgacptac.rds` best_mean 0.6545342172, best_se 0.01098411963, lcb_threshold 0.6435500975, k_selected 3. `tar_params_best_tcgacptac.rds` k 3, alpha 0.3339355562, lambda 0.3494487379, nu 0.0556512168, ntop 270. `final_config_cv_tcgacptac.rds` (desurv_version 1.0.1): k7 best 0.654822868 (SE 0.0106712), best k3 0.6466307793 (SE 0.008651577, alpha 0.4422075, lambda 0.2409568, nu 0.2014809, ntop 299), final 0.6419248411 (SE 0.006918191). Same-build threshold 0.6441516682; final below by 0.002226827. All match the audit. Ladder S0 also shows 0.6419 is below the original May threshold 0.64355.

Disclosure: Methods line 86 states the hybrid, that it was not evaluated in the search, the re-evaluated value, the threshold and the gap, and "The consensus initialization described next was used only for the final model and is not reproduced within each fold." SI lines 811-879 and Table final-config-cv repeat this. Revised severity LOW-MEDIUM: a selection-procedure weakness, not misreporting. Note the BO fold SE is fold-level, descriptive; the manuscript does not present it as a post-selection CI. Smallest decisive check: none needed for honesty; for substance, the rev 3 reported-estimator CV already evaluates the full recipe in folds (ladder S4: published k3 0.6362). Remediation: manuscript-only (done); a principled fix would select the best row conditional on k, which is a rerun decision for the owner.

## F4. Sensitivity grid preprocesses outside folds

Verdict: CONFIRMED.

Evidence. `code/06_cv_grid.R:80-101` runs `DeSurv::preprocess_data(..., ngene = 3000)` once on the pooled TCGA+CPTAC training data, then the caller at line 146 passes the result to `run_cv_grid_point()` (`R/cv_grid_helpers.R:79`), which subsets columns per fold (lines 121-131) and fits `desurv_fit(ninit = 30, tol = 1e-5, maxit = 3000)` on fold training columns. `desurv_fit` picks a single training-selected initialization, unlike BO's averaging over starts. BO uses fold-local preprocessing (`code/03`, `preprocess = TRUE`).

Experiment (read-only load_data + preprocess_data, folds `.desurv_make_folds_stratified(event, dataset, 5, 123)`): full-data filter gives 1970 genes, identical to cached `tar_data_filtered` (genes and matrix all.equal TRUE). Fold-local filters: 1945, 1957, 1944, 1944, 1967 genes; full-only 73, 80, 86, 95, 78; fold-only 48, 67, 60, 69, 75. Exactly the audit's numbers.

Grid (`results/cv_grid/cv_grid_summary.csv`, ntop 270): best k 9 alpha 0.25, C 0.67059006, SE 0.008928048, threshold 0.66166202; best k3 alpha 0.40 C 0.66298174 qualifies; k3 alpha 0.35 C 0.65646761 (SE 0.00996) does not. Matches the audit and the SI's own stopifnot guards (si_appendix.Rmd:752-754).

Proven vs conjectured: held-out expression influences unsupervised gene filtering (about 4 percent of genes per fold); outcome leakage is not involved (the filter uses mean and variance only, and within-sample ranks do not cross samples). Effect on grid CV C unmeasured, plausibly small. Primary BO is unaffected.

Disclosure: SI "Cross-validated C-index across factorization rank" (line 2014) and SI 830-855 describe the grid as a sensitivity analysis and already say the final configuration's neighbor fails the grid's own threshold; neither discloses the outside-fold filter or the single-init estimator. Revised severity LOW. Smallest decisive check: rerun the k 3 and k 9 grid rows fold-locally (10 points) and compare. Remediation: manuscript-only disclosure sentence, or targeted rerun.

## F5. Tuned supervised PCA is a correlation-weight surrogate

Verdict: CONFIRMED.

Evidence. `code/19_desurv_vs_supervised_tuned.R:83-87` `superpc_dir()`: native `superpc.predict` on training, components combined by training Cox, then gene weights `w[keepg] <- cor(gene, scr)`; external scoring at line 58-66 (`transfer`) projects external rank expression through `w`. Native predictor (training-fitted component structure and centering) is never applied to validation data. Orientation (`al`, line 95) flips each comparator to correlate positively with D1.

Experiment (`sub/f5.R`, same seed 1, 20 thresholds, 5 folds): selected 3 components, threshold 2.476, 5 genes; native vs surrogate training score correlation 0.9569627. External C (same PACA dedup, protective orientation for both):

| Cohort | Native | Surrogate |
|---|---:|---:|
| Dijk | 0.6208054 | 0.6239029 |
| Moffitt | 0.5682978 | 0.5699981 |
| PACA array (17) | 0.5877193 | 0.5877193 |
| PACA RNA-seq | 0.6282195 | 0.6315789 |
| Puleo | 0.5932278 | 0.6055994 |
| Event-weighted | 0.5962774 | 0.6033275 |

Identical to the audit. Cached `desurv_vs_supervised_tuned.rds$prediction` pooled supervised PCA 0.603 is the surrogate.

Consequence: the surrogate flatters supervised PCA by 0.007 pooled; native 0.596 is further below DeSurv full LP (0.613), so the Results conclusion (comparable discrimination, different biology) is not reversed. The gene-set attribution (selected genes, enrichment) is unaffected; the |r| with D2 (Results line 615) is computed on the surrogate score and would need recomputing.

Disclosure: not disclosed; Methods line 119 states "the fitted model was frozen and applied to rank-transformed validation cohorts without refitting", which the surrogate does not do. Revised severity LOW-MEDIUM. Remediation: code correction (use native `superpc.predict` with training-derived Cox combination) and update the pooled value and the D2 |r|. Smallest check is already done above.

## F6. Absent assay genes zero-filled before ranking

Verdict: CONFIRMED as a mechanism and as a text mismatch; impact narrowed to negligible.

Evidence. Actual caller `code/05_external_validation.R:46-56` `load_val_datasets()` calls `preprocess_validation_data(..., zero_fill_missing = TRUE)` (line 114); `R/preprocess_helpers.R:69-81` rbinds zero rows for missing genes before `DeSurv::preprocess_data` ranks. `code/06_cv_grid.R:116` and code/19 (via `data_val_filtered`) inherit the same objects.

Counts (raw cohort gene lists vs 1970 namespace and 810 union support): Dijk 25/15, Moffitt 73/27, PACA array 11/6, PACA RNA-seq 28/16, Puleo 78/28. Inserted ranks: Dijk 18.5-95.5 (33 distinct values), Moffitt 37, PACA array 6, PACA RNA-seq 15-30, Puleo 39.5. All match the audit.

Narrowing experiment (`sub/f6.R`): recomputed the primary union-support score from cached `data_val_filtered` (max |difference| vs cached `val_latent` risk_score 6e-14), then the alternative policy the treated cohorts use (rank only measured namespace genes, drop missing support genes). When inserted ranks are constant across samples (Moffitt, PACA array, Puleo), the two scores differ by a constant and C is identical. C zero-fill vs drop-missing: Dijk 0.6073826 vs 0.6076407 (+0.00026), Moffitt, PACA array, PACA RNA-seq, Puleo identical to 7 decimals; score correlation at least 0.999999 in every cohort. The audit's "Impact NEEDS VERIFICATION" is resolved: negligible for primary C and, given correlation 0.999999, for the HRs.

Disclosure: not disclosed. Methods line 163 "No missing-data imputation was performed" refers to complete cases, and SI lines 936-946 say validation sets were "restricted to this same gene set" and ranked "by exactly the same per-sample procedure", which omits the zero fill. SI line 2392 discloses the different drop rule for treated cohorts. Revised severity LOW. Remediation: manuscript-only (one sentence stating that absent genes enter at the bottom rank, with the sensitivity numbers above).

## F7. Provenance and self-contained rebuild

Verdict: CONFIRMED for each item checked; not a correctness finding.

- `code/01_install.R:13-22`: pins `submission/pnas-2026` but installs `../DeSurv` unpinned when it exists (with a warning message only). `../DeSurv` is HEAD 370c88a, reports Version 1.0.1, and does not contain afb00d5; its functions.cpp equals the 1.0.1 extract, so the C++ is the same, but R-level identity with afb00d5 is not established here.
- `paper/paper.Rmd:150` promises "the exact session-info file"; no session-info file exists in the tree (find and git ls-files return nothing).
- `paper/paper.Rmd` versioning chunk (around line 175) is `eval=FALSE` and reads `targets::tar_read("version_info")`, a store this repo does not use.
- SI table note (si_appendix.Rmd:870-872) discloses that the cached BO history came from an unrecorded DeSurv version; best k3 0.6466 now vs historical 0.6462.
- `code/09a_figures.R:147-155` (also code/23 and code/24) fall back to `../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData`; the local `data/derv/` file is absent; no generating script found.
- `R/load_data_internal.R:79` consumes curated `*.survival_data.rds`; no raw-to-curated endpoint derivation in the repo.

Revised severity MEDIUM for reproducibility. Remediation: deterministic verification (hash manifest of inputs and results, sessionInfo file, the missing formatted subtype object or its generator), manuscript fix of the session-info promise, and pin enforcement in 01_install.R.

## Positive findings spot-check (from cached .rds, independent code)

- Per-cohort C from `val_latent_desurv_tcgacptac.rds`: Dijk 0.6073826, Moffitt 0.5616852, PACA array 0.6490455, PACA RNA-seq 0.6338186, Puleo 0.6359692. Match.
- After `dedup_combined`: 570 patients, 388 events; PACA-AU combined 0.6303550; event-weighted pooled 0.6134885. Match.
- Stored raw beta (-0.001002651, 0, 0.00009359101); union support 810. Match.
- Full-data preprocessing equals cached object exactly (F4 experiment). Match.

## Mismatches with the audit

None numerical: every number checked reproduced to the printed precision. Narrowings: F1 severity HIGH to MEDIUM (CV impact measured at about 0.0003; W and beta impact still open); F6 impact resolved as negligible (largest C change +0.00026, Dijk only), because constant inserted ranks only shift scores by a constant. F5 surrogate inflates supervised PCA, so correcting it strengthens rather than weakens DeSurv's relative position.
