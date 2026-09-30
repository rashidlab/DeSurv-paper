# DeSurv rerun: plan rev 3 (2026-09-30)

Supersedes rev 2.1 (frozen at 960b76d) for everything after its gate. Rev 2.1's consensus-of-final-fits
estimator FAILED its prospective gate (Slurm 3109627: competing modes in all 12 families at R = 50), and the
one follow-up, alpha = 0 homotopy, FAILED its prespecified test (Slurm 3112383: alpha = 0 not stable; min pair
r 0.85-0.90; endpoints r ~0.25-0.28 with the reported fit's first factor, training C ~0.57). Per the hard stop,
method development for a unique W is CLOSED for this paper.

Owner scope (2026-09-30): "correct the implementation -> correctly retune the existing estimator -> rerun the
full analysis -> independently characterize solution sensitivity -> revise claims accordingly."

STATUS: PROPOSAL. Owner approves this file and the section 4 threshold before execution.

## 1. What is kept, what is removed

Kept:
- Cox tie fix, DeSurv 1.0.2 (`fix/cox-breslow-ties`, 4c4809b + 72d8495 untracking stale objects). Mandatory.
- The existing estimator, unchanged: seed fits (code/04 `run_seed_fits`) -> `desurv_consensus_seed` -> one final
  `desurv_fit` (code/04 `consensus_and_fit`).
- Track 2 correctness pieces only (section 2).
- Generic hygiene that is not estimator-specific: deterministic top-gene tie-break (score, then gene name; FP15)
  (owner to confirm; it changes top-gene sets only on exact score ties).

Removed (not shipped, not run again): `desurv_consensus_fit` / `desurv_consensus_gate`, the `fix_W` flag, the
R = 50/100 gate, alpha homotopy, NNDSVD, alternative NMF formulations. Branches `feat/consensus-estimator`,
`feat/consensus-integrated` stay as a record, never merged or released. The two failed experiments are
documented results (summaries in `docs/plans/gate3a_summary_3109627.txt` and the homotopy summary); whether the
paper mentions them (e.g. one SI paragraph on estimator sensitivity) is an owner call at the text stage.

## 2. Retuning the existing estimator (Track 2, re-pointed)

- The BO objective in each CV fold is the REPORTED estimator: the code/04 recipe (seed fits -> consensus seed ->
  final fit) on the training fold, scored on the held-out fold through the existing top-ntop scoring rule. The
  `estimator = "consensus"` option is replaced by this recipe (call it `estimator = "reported"`).
- Seed fits per fold: the same count as the final fit (NINIT_FULL = 100) unless a runtime benchmark of a few BO
  evaluations shows that is impractical; then the owner decides (e.g. 30), recorded before the BO runs.
- 3 repeated CV partitions (fold seeds 20261101-3), identical for every configuration; per-fold seed banks
  identical across configurations (paired).
- BO incumbent by GP posterior mean; top 2-3 configurations per k re-evaluated on the identical repeated CV;
  k by the 1-SE rule on the empirical paired CV SE; runner-up k reported.
- All hyperparameters from the selected k's own row; no typed hyperparameters downstream (code/06, code/19).
- Implementation: a new package branch from 72d8495 carrying the tie fix plus the Track 2 commits minus the
  consensus-estimator path (0408f59 is based on 72d8495 and does not contain Track 1; its stub file and the
  "consensus" option are replaced by the "reported" recipe). Paper-repo Track 2 changes (59c3d4e) kept.
  Tests 8 and 9 kept; test 8 now asserts the fold estimator is the code/04 recipe.

## 3. Primary fit (after retuning; never informed by the sensitivity runs)

Hyperparameters frozen from section 2. Primary fit = code/04 recipe with the pipeline's own seed convention
(seed fits 1..NINIT_FULL), as before. This, and only this, is the reported model.

## 4. Sensitivity analysis (after the primary fit is frozen)

- 10 complete, independent reruns of the frozen estimator (seed fits -> consensus seed -> final fit) with
  disjoint seed banks: rerun j uses seeds 3,000,000 + (j - 1) x 10,000 + (1..NINIT_FULL).
- Per rerun, factors aligned to the primary fit (Hungarian on correlation), reporting: per-factor Pearson and
  Spearman r with the primary fit, top-ntop gene overlap, beta, validation HRs and C-indices, and whether each
  claim in the claims inventory holds.
- Claims inventory: before the reruns, list every substantive Results/Discussion claim that depends on the fit
  (e.g. L443 "D1 the smallest reconstruction share", L463 "D1 the only program with a significant combined
  association", the D3 per-cohort pattern, the pooled and adjusted HRs being significant), each with an
  operational test. Committed before any rerun output exists.
- Rule (PROPOSED, owner confirms): a claim that holds in at least 9 of 10 reruns stays as an ordinary finding;
  a claim holding in 5-8 of 10 is described as sensitivity-dependent with its frequency; below 5 of 10 it is
  removed. Validation data enter only this descriptive characterization, never selection.
- The sensitivity runs never change the hyperparameters, the primary fit, or which solution is reported.

## 5. Full rerun (Longleaf, staged)

- Pre-run: code/26 version hard-stop; code/01 DESURV_REF; paper.Rmd L150 typed version -> inline; Slurm coverage
  for 11, 14-18, 20-26, 29; CV-grid walltime; cmbSubtypes file; a producing script for
  `paper/gene_lists_top270_k3.csv`; seed `NMF::nmf` (code/03:183-187) and `cv.glmnet` (R/fit_cox_model.R);
  fix `preprocessCore` on Longleaf (7 quantile-preprocessing test errors) if step 2/5 paths use it; the
  `setup_env.sh` test invocation.
- Order: inventory's topological list; OWNER STOP after the BO (k, hyperparameters, runner-up k, empirical
  SEs). If k != 3: labels and "three programs" wording need rulings.
- Clean run into an empty results directory on the recompute branch; DESURV_RECOMPUTE explicit; canonical
  `results/` untouched until the owner accepts (before side = fb731e7).
- Simulations: the estimator is unchanged and simulated times have no ties, so they are rerun only if their
  hyperparameters come from the BO (check code/07 at the owner stop).
- code/06 CV grid: rerun, re-scope or drop (owner, at the stop).
- code/21 stability section is replaced by section 4.

## 6. Numbers and text

Inline-hook diff of every manuscript number; typed values (abstract HRs, "seven-factor", k); claims revised per
the section 4 rule; Methods/SI describe the tie fix, the tuning procedure and the sensitivity analysis; text
through paper-system transactions; release gate and owner read afterwards.

## 7. Owner decisions

1. Approve rev 3.
2. The section 4 claim rule (9/10 ordinary; 5-8/10 sensitivity-dependent with frequency; <5/10 removed).
3. Keep the deterministic top-gene tie-break (section 1)?
4. Longleaf compute budget; the BO seed-fits-per-fold count if the benchmark shows 100 is impractical.

## 8. APPROVED (2026-09-30)

Owner (verbatim): "Approve Rev 3. Approve the 9/10 / 5-8/10 / <5/10 claims rule with frequencies reported in
the sensitivity analysis. Keep the deterministic top-gene tie-break. Benchmark NINIT_FULL=100 and use 100 unless
measured runtime makes it genuinely impractical; no smaller count is approved prospectively yet."

Clarification of the section 4 rule (owner): a claim meeting 9/10 can remain an ordinary finding in the main
text, but its full sensitivity frequency is reported in the sensitivity results/SI. Claims meeting 5-8/10 must
explicitly be described as solution-sensitive; claims below 5/10 are not asserted. The operational definition of
each claim (exact test, CI or comparison criterion) is frozen in the claims inventory before any rerun. No
intervals around the 9/10 cutoff; these are transparent sensitivity categories, not inferential thresholds.

BO seed fits per fold: 100 (NINIT_FULL). If the benchmark shows 100 is genuinely impractical, the measured
runtime goes to the owner before any change. The primary fit (seeds 1..NINIT_FULL) is frozen before the ten
sensitivity reruns and is never replaced by a rerun. Methodological design is CLOSED.

Execution order:
1. Package: new branch from 72d8495 (tie fix) + the Track 2 correctness commits, re-pointed to the "reported"
   estimator (code/04 recipe per fold); deterministic top-gene tie-break; tests (incl. 8, 9); install on Longleaf.
2. Paper repo: section 5 pre-run fixes; claims inventory committed (before step 6).
3. Benchmark a few BO evaluations at NINIT_FULL = 100 on Longleaf; report runtime.
4. BO with 3 repeated CV partitions -> OWNER STOP (k, hyperparameters, runner-up, empirical SEs).
5. Primary fit + full downstream rerun (staged); owner accepts before results/ is replaced.
6. Ten sensitivity reruns; claims scored against the frozen inventory.
7. Numbers diff and text transactions; release gate.
