# Track 3a gate script: report (2026-09-30, revised for the fixed-R amendment)

This revision implements `docs/plans/2026-09-30-section3-amendment-fixedR.md` (d120e56) on top of the committed script (e4195a2).
- Nothing is committed.
- Nothing has run on Longleaf.
- Only the three files in `analysis/gate/` and this report changed.

The first version's open point (R = 25 cannot be gated because R must be even) is moot. The amendment supersedes it.

## Files

- `analysis/gate/gate_3a.R`: the probe. Rewritten around the fixed-R rule. The R ladder, the 2R rule and `GATE_R_VALUES` are removed.
- `analysis/gate/gate_3a.sbatch`: unchanged resources (general, 32 cpus, 64G, 12 h, `r/4.5.0`, briefed paths). The header comment names the amendment, and the script unsets `GATE_FORCE_BORDERLINE`.
- `analysis/gate/README.md`: rewritten. It covers staging (the same two rds files), submitting, and reading the output.

## Implementation of the amendment

- **Stage R = 50.** Fits only seeds 1,000,001-1,000,200, in checkpointed chunks of 50. Each chunk is one ensemble block.
  - One `desurv_consensus_fit(run_fits = block)` per ensemble E1-E4, run in parallel.
  - For each of the 6 pairs: Hungarian alignment (`DeSurv:::.dc_align`, the estimator's own rule) and per-factor Pearson r.
  - Factor identity across pairs comes from mapping each ensemble to E1. Each pair's partner factor comes from that pair's own Hungarian alignment, as the amendment says. The column `consistent_with_E1_map` flags any disagreement between the two alignments.
- **Verdict**, in this order:
  1. Any ensemble fit error, or any `mode_pass` FALSE: FAIL, with the reason.
  2. Otherwise, per-factor median >= 0.98 and min >= 0.95: PASS.
  3. Otherwise, all medians >= 0.95: BORDERLINE.
  4. Otherwise: FAIL.
- **Stage R = 100.** On BORDERLINE the same job fits seeds 1,000,201-1,000,600 (4 blocks of 100) and applies the same PASS rule. Anything but PASS is a FAIL.
- **Failed constituent runs.** A failed run is not dropped and does not stop the job. Its ensemble is recorded as a fit error, which makes the stage FAIL, and the seeds go to `bank_failures_R<R>.txt`.
- **Descriptive outputs**, never used for the verdict:
  - `desurv_consensus_gate` half-vs-half within each ensemble (split seed 20260930, 10 splits);
  - beta, training C and top-gene overlap (ntop 270) for each pair;
  - occupancy, tightness, n_members, n_outliers and competing-mode medoid r for each ensemble;
  - the code/04 reference on 1,100,001-300 (pairwise Pearson, Spearman and top-gene overlap).
- **Kept from the first version:** preflight, sessionInfo, input md5 sums, the guard against writing under `results/`, the code/04 fit settings (tol 1e-7, tol_init 1e-5, maxit 4000) and training data only.

## Local dry run (correctness only; Track 1 library, paper-clean results read-only, 4 cores)

1. **`GATE_DRY=TRUE`** (4 ensembles of 6 runs, 2 half-splits, one 4-seed reference set).
   - Exit 0, 0.6 min.
   - Verdict: FAIL at R = 6 (competing mode in all 4). At this size that is expected and not a result.
   - The pair table, per-ensemble lines and reference all printed.
2. **`GATE_FORCE_BORDERLINE=TRUE` without `GATE_DRY`.** Refused at startup, as intended.
3. **`GATE_DRY=TRUE GATE_FORCE_BORDERLINE=TRUE`**, run in the same output directory.
   - All 5 existing checkpoints were loaded.
   - The stage-1 pair table and ensemble W were identical to run 1.
   - The second stage fitted only the 48 new seeds (1,000,025-072).
   - The verdict line is labelled `[DRY RUN, borderline FORCED; not a gate result]`.
4. `git status results/` in paper-clean is clean.

## Not exercised by the dry run

At dry-run sizes every ensemble has a competing mode. So the numeric PASS and BORDERLINE branches (the threshold comparisons) were reached only through the forced flag, not through real numbers. The logic is short and inline (`run_stage`, verdict block). A one-line code read is the cheapest check before submission.
