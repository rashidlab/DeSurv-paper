# Robust full rerun of the DeSurv analysis (plan, rev 2.1, 2026-09-30)

Owner ask (2026-09-30, verbatim): "rerun the entire paper analysis ... I want the rerun version to have
more random restarts, more consensus stability fits, and other proposed changes that would make the final
result less immune to the choice and number of restarts and etc. That is, we want to make sure that we
rerun things so that the results are less likely to change with more runs." The paper is not submitted:
the rerun uses the fixed package and the manuscript carries no correction note.

Rev 2 replaces rev 1's A1-A5/B1-B6/C1-C5/D1-D4 structure after an external review of rev 1, accepted by
the owner: fix one estimator instability, fix two tuning errors, show independent reproducibility. Rev 1
is not preserved (docs/ is gitignored); its evidence section is kept below in short form.

STATUS: rev 2.1, revised per a second external review (three edits: W-only gate, mode check before
outlier removal, ntop defined as a scoring-rule parameter; plus the R_BO principle). The review
recommends no further redesign. Awaiting the owner's approval and the section 3 freeze. Section 3 (prespecified rules) must be confirmed or edited
by the owner and frozen BEFORE the first convergence probe runs.

## 1. Evidence (details in the files named)

- `2026-09-30-tie-fix-sizing-report.md`: the Cox tie fix alone barely moves the reported fit (same
  initialization: W r >= 0.998, pooled HR 1.42 -> 1.43). Rerunning the seed fits under the fix moves the
  consensus start slightly (r ~0.987) and the final fit into a different, higher-objective basin
  (D1 r 0.90).
- `2026-09-30-robustness-code-audit.md` (FP1-FP15, file:line):
  - the consensus keeps genes by a hard count cut (`min_frequency`; 308 genes within +-5 of it) and a
    `cutree`, so dropping ONE of 100 restarts changes the final D1/D3 (r 0.95-0.97)
  - the BO scores the mean of 30 single random starts, not the consensus estimator that is reported
  - `code/03:117-123` takes alpha/lambda/nu/ntop from the best k = 7 row and overwrites only k, so the
    reported configuration was never evaluated (`select_bo_params_1se` exists, unused)
  - one fold split throughout; fold noise (SE 0.011) is ~8x restart noise
  - code/21's "r 0.97-0.99" is not independent of the reported fit
- `2026-09-30-robustness-probe.md`:
  - 100, 300 and 1,000-start ensembles with disjoint seeds never agree (5, 5, 3 distinct endpoints;
    worst D1 r 0.92, 0.81, 0.95)
  - D3 is significant in validation in 8 of 13 endpoints, so Results L443/L463 are basin-dependent
  - lower objective tracks worse validation (the global minimum found: training C 0.90, validation 0.57)
  - neither multi-start around the consensus nor a filtered consensus fixes it

Diagnosis: a multimodal, non-identifiable factorization plus a discontinuous consensus rule. More
restarts search the same landscape harder; they do not give a unique answer. Precedent for the remedy:
consensus NMF on gene programs (Kotliar et al., cNMF): normalize the factors from many runs, cluster them,
drop outlying components, take gene-wise medians, refit usages with the consensus programs fixed.

## 2. Architecture: three tracks, one scientific gate

```
                 ┌─ Track 1: consensus estimator + smoke tests ─┐
start ───────────┤                                              ├─> Track 3a: convergence probe (current config)
                 └─ Track 2: tuning correctness fixes + tests ──┘            │
                                                              ┌─ GATE passes ─┴─> Track 2 run: repeated-CV BO
                                                              │                    │
                                                              │                    └─> Track 3b: same gate at FINAL config
                                                              │                              │
                                                              │                              └─ passes ─> full rerun, sims,
                                                              │                                          replicate stability, text
                                                              └─ GATE fails ──> STOP; only then consider section 7
```

Tracks 1 and 2 start in parallel now. Smoke tests are engineering gates only: passing them unlocks the
convergence probe, never the expensive BO or the full rerun.

### Track 1: the consensus estimator ("consensus-DeSurv"), a DeSurv package change

At fixed (k, alpha, lambda, nu). ntop is not an estimator parameter: `desurv_fit` never uses it. It
defines the scoring rule (the risk score uses each factor's top-ntop genes, `.desurv_lp_with_top_genes`,
pkg R/cv_warm.R:293), so it is applied to the consensus W and is still tuned by the BO because it changes
the cross-validated C-index.
1. R independent DeSurv fits from independent random starts (seeds from a prospectively defined seed
   bank, so larger R contains smaller R as a prefix). Each fit runs to convergence. No fit is selected by
   objective or training C-index, anywhere (this also retires `init()`'s best-by-training-C rule).
2. L2-normalize every W column.
3. Pool the R x k factor vectors and cluster them into k families (correlation distance). No
   one-factor-per-run constraint.
4. Check each family for a competing mode (section 3) BEFORE any outlier removal. A competing mode
   fails the gate; it is never filtered away.
5. Only if no competing mode: remove isolated outlying factors by the prespecified geometric rule
   (section 3). Outlier status never depends on survival, C-index or validation performance.
6. Consensus W: gene-wise median of each family's remaining members, renormalized.
7. Refit H and beta with the consensus W held FIXED. No free refinement of W afterwards.
8. Only then derive the top-ntop genes of the consensus W for the scoring rule, with a deterministic
   tie-break (score, then gene name).

Report per run: family occupancy (0, 1 or more than 1 factors per family); per family: within-family
tightness (median pairwise r) and the number of factors removed as outliers.

Consequence for the paper: the reported W is no longer a minimizer of the DeSurv objective; it is the
consensus of an ensemble of DeSurv solutions, with H and beta fitted conditional on it. Methods and SI
say so explicitly (the optimization/convergence passages, including those edited in TX-0016, are
rewritten).

Package version: build on `fix/cox-breslow-ties` (72d8495, v1.0.2 local, unpushed) and release as 1.1.0.

### Track 2: tuning correctness (paper repo + package BO code)

- The BO objective inside each CV fold is consensus-DeSurv at R_BO runs. Principle: R_BO is the smallest R
  that passed the gate in Track 3a, unless a runtime benchmark of a few BO evaluations shows that is
  impractical (then the owner decides). The BO tunes the estimator that is reported, not a cheaper
  approximation of it.
- 3 repeated CV partitions, the SAME partitions and seed banks for every configuration (paired
  comparisons).
- BO incumbent chosen by the GP posterior mean, not the best single noisy evaluation.
- Finalists: the top 2-3 configurations per k are re-evaluated on the identical repeated CV; k by the
  1-SE rule on the EMPIRICAL paired CV SE (not the GP's posterior SE). Report the runner-up k.
- All hyperparameters come from the selected k's own row (`select_bo_params_1se`; code/03:117-123 and
  160-167 for the alpha = 0 control).
- No typed hyperparameters downstream: code/06:49-51 and code/19 (ntop 270) read `tar_params_best`.
- Stability is reported at the chosen k and the runner-up; it does not choose k. The chosen k must still
  pass the gate.

### Track 3: convergence (the gate)

- 3a, now, at the current configuration: R = 25, 50, 100, 200 as needed (prefixes of one seed bank).
- 3b, mandatory, after tuning: the identical check at the final configuration.
- Current-estimator reference: rerun the existing pipeline on a few disjoint seed sets to reproduce the
  instability under the same scripts, so the comparison is like for like.

## 3. Prespecified rules (owner confirms or edits, then FROZEN before 3a runs)

Proposed values; the point is that they are written down before any output is seen. They are
engineering tolerances, not scientific constants: freeze them and never tune them after seeing results.

- **Outlier rule (step 4):** for each L2-normalized factor, the mean correlation distance (1 - Pearson r)
  to its L = ceil(0.3 R) nearest factors. Outlier if that statistic exceeds median + 3 x MAD of all
  factors' statistics.
- **Half-versus-half agreement:** 10 random splits of the R runs into two disjoint halves; consensus W on
  each half; align factors; per-factor Pearson r. Pass if, for every factor, the median over splits is
  >= 0.98 and the 5th percentile is >= 0.95.
- **Competing mode:** within each family, split members into two sub-clusters (k-medoids, 2). A competing
  mode exists if the smaller sub-cluster holds >= 20% of the family's members AND the two sub-cluster
  medoids correlate at r < 0.95. Any competing mode FAILS the gate; it is never filtered away as an outlier.
- **The gate is about W only:** pass = the half-versus-half rule and no competing mode. Beta, training
  C-index, HRs and top-gene overlap are reported descriptively as downstream stability outcomes; they do
  not decide pass or fail. Validation-cohort quantities never enter the gate, the outlier rule or any
  estimator choice.
- **Order of operations:** families, then the competing-mode check, then (only if it passes) outlier
  removal, then the median.
- **Convergence in R:** the smallest R at which the rules above pass, and they still pass at 2R. Report
  that R; the final fit uses at least it.

## 4. Tests

Estimator smoke tests (package, tiny synthetic data, run in CI and `make quick`):
1. Determinism: same seeds and inputs give identical W, H, beta.
2. Run-order invariance: shuffling the constituent fits gives the same consensus W.
3. Factor-label invariance: permuting factor columns within constituent fits gives the same consensus
   after alignment.
4. One aberrant fit: adding a deliberately aberrant fit to near-identical fits leaves the consensus
   essentially unchanged (outlier rule plus median).
5. More runs: adding runs from the same neighborhood keeps the consensus at r ~1.
6. Fixed W: after the H/beta refit, W is bit-identical to the consensus W.
7. Top genes after aggregation: ntop selection happens on the consensus W, never inside constituent fits.

Pipeline tests (paper repo):
8. BO estimator identity: the object scored inside CV comes from the same entry point as the final fit.
9. Selected-k integrity: on a toy BO table with distinct values per k, selecting k = 3 returns the whole
   k = 3 row.

Preflight / environment check (every job): the installed DeSurv is the pinned version and contains the
fix (`exists("cox_stats_cpp", asNamespace("DeSurv"))` plus tie test T1). Tracked stale objects already
made a checkout install run 1.0.1 code once. Existing `test-cox-ties.R` carries over.

## 5. Full rerun (after the gate passes at the final configuration)

- Hygiene: seed `NMF::nmf` (code/03:183-187) and `cv.glmnet` (R/fit_cox_model.R:10-15); clean run into an
  empty results directory on a new branch, no mixed cache hits; DESURV_RECOMPUTE set explicitly in every
  command; canonical `results/` untouched until the owner accepts (before side = commit fb731e7).
- Tie-fix plan items still open (`2026-09-30-cox-tie-fix-plan.md` B0a): code/26 version hard-stop,
  code/01 DESURV_REF, paper.Rmd L150 typed "v1.0.1", Slurm coverage for 11, 14-18, 20-26, 29, CV-grid
  walltime, cmbSubtypes file on Longleaf, a producing script for `paper/gene_lists_top270_k3.csv`.
- Order: the inventory's topological list (`2026-09-30-tie-fix-recompute-inventory.md`), staged on
  Longleaf, with an OWNER STOP after the BO (new k and hyperparameters, runner-up k, empirical SEs). If
  k != 3, program labels and "three programs" wording need rulings before anything downstream.
- Simulations (step 7): rerun with consensus-DeSurv. Their claims compare DeSurv with other methods on
  program recovery, so they test the reported estimator. Before rerunning, list each simulation claim
  and confirm it concerns the estimator; any claim about a lower-level optimization property can keep
  its current run. About 6-12 h.
- code/06 CV grid: decide at the owner stop whether it is rerun with the new estimator, re-scoped or
  dropped (it is 24-48 h with the old estimator and more with the new one).
- Independent replicate stability: 10 complete consensus-DeSurv fits at the final configuration with
  disjoint seed banks. Report per-factor r (Pearson and Spearman), top-gene overlap, beta and HR ranges.
  This replaces code/21's non-independent check.
- Claims: every Results claim is checked against the replicate distribution. Claims that vary across
  replicates are stated as frequencies or dropped (currently L443, L463 and D3's significance).

## 6. Numbers and text

The inline-hook diff of every manuscript number (tie-fix plan Part C), typed values (abstract HRs via
check_abstract_hr.R, "seven-factor", k), then paper-system transactions for Results, Methods and SI,
including the estimator definition (section 2, Track 1, "Consequence for the paper"). Release gate and
owner read after that.

## 7. Deferred: only if the gate fails

alpha-homotopy (forward and reverse, to detect path dependence), NNDSVD or other deterministic
initialization, patient bootstrap stability, constrained factor matching, stronger NMF constraints
(minimum-volume, orthogonal, sparse), Bayesian NMF. None of these is built or run unless the
consensus estimator fails the gate at the final configuration.

## 8. Owner decisions

1. Confirm the section 3 values as written (or edit once); then they are frozen, with the date recorded here.
2. R_BO target and the Longleaf compute budget (cores x hours).
3. code/06 CV grid: rerun, re-scope or drop (can wait until the owner stop after the BO).
4. Push `fix/cox-breslow-ties` now as 1.0.2, or release only as 1.1.0 with the estimator.

Rough elapsed time: Tracks 1-2 build and tests 2-3 days; probe 3a under a day; BO 1-3 days of cluster
time; gate 3b hours; full rerun plus simulations 3-5 days; numbers and text 2-3 days. About 2-3 weeks.

## 9. FROZEN (2026-09-30)

Owner (2026-09-30, verbatim): "Approve Rev 2.1. Freeze Section 3 exactly as written. Record the
commit/hash and gate seeds. Start Tracks 1 and 2 in parallel. No further methodological redesign unless
the convergence gate fails."

Section 3 is frozen exactly as written above. Nothing in sections 2-3 changes unless the gate fails.
Freeze record: this file is committed with `git add -f` on the local branch `recompute/consensus-desurv`
of the paper repo (worktree `~/Downloads/DeSurv-paper-recompute`). The commit hash and this file's
sha256 are recorded in `docs/HANDOFF-2026-09-30.md` (a file cannot contain its own commit hash).

Gate seeds (prospective; disjoint from every seed used so far: cached seed fits 1-100, probe <= 3,000,
BO/CV small integers):
- Track 3a run bank: seeds 1,000,001 to 1,000,400 in order (R = 25, 50, 100, 200 are prefixes; 400 covers
  the 2R check at R = 200). Each run: `desurv_fit(..., seed = s)` from its own random start.
- Track 3a half-split draws: `set.seed(20260930)`, then 10 successive `sample(R)` permutations; halves
  are the first and second R/2 of each permutation.
- Track 3a current-estimator reference: three disjoint 100-seed sets, 1,100,001-1,100,100,
  1,100,101-1,100,200, 1,100,201-1,100,300, through the existing code/04 recipe.
- Track 3b run bank: 2,000,001 onward (same prefix rule); half-split draws `set.seed(20261001)`.
- Track 2 repeated-CV partitions: fold-assignment seeds 20261101, 20261102, 20261103, identical for every
  configuration; per-fold run bank 4,000,001 onward, identical for every configuration (paired).
- Final independent replicates (section 5): replicate j = 1..10 uses 3,000,000 + (j - 1) x 10,000 + (1..R).
- Competing-mode check: k-medoids with the deterministic PAM build/swap (no random seed involved).
If R must exceed 200, banks extend upward in the same order; no seed is ever reused across tracks.
