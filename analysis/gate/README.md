# Track 3a convergence probe (fixed-R rule)

These documents govern the probe:
- plan rev 2.1 (frozen 960b76d);
- the section 3 clarifications (0dc8aff);
- the section 3 amendment `docs/plans/2026-09-30-section3-amendment-fixedR.md` (d120e56). It replaces the R ladder and the 2R rule.

The script requires DeSurv 1.1.0 (`feat/consensus-estimator`, b42d7a9). It stops unless `packageVersion("DeSurv") == "1.1.0"` and `DeSurv:::cox_stats_cpp` exists.

## What it does

1. **Stage R = 50.** Fits seeds 1,000,001-1,000,200 and builds four independent consensus-DeSurv ensembles from them:
   - E1 = 001-050
   - E2 = 051-100
   - E3 = 101-150
   - E4 = 151-200

   For each of the 6 ensemble pairs, it aligns the factors (Hungarian on Pearson r) and computes per-factor Pearson r. The verdict:
   - **PASS:** every factor has median >= 0.98 and min >= 0.95 over the 6 pairs, `mode_pass` is TRUE in all 4 ensembles, and no ensemble fit errors.
   - **BORDERLINE:** every factor median is >= 0.95, `mode_pass` holds in all 4, there are no errors, but PASS is not met.
   - **FAIL:** anything else.
2. **Stage R = 100.** This runs only if R = 50 is BORDERLINE, and it runs in the same job. It fits seeds 1,000,201-1,000,600 (blocks of 100) and applies the same PASS rule. Anything other than PASS here is a FAIL.
3. **Descriptive outputs (never pass/fail):**
   - a half-vs-half `desurv_consensus_gate` within each ensemble (split seed 20260930, 10 splits);
   - beta, training C and top-gene overlap (ntop = 270);
   - occupancy, tightness and outlier counts;
   - the current-estimator reference: the code/04 recipe on seeds 1,100,001-100, 101-200 and 201-300, with pairwise aligned W correlations.

A failed constituent run is never dropped. It makes its ensemble a fit error, so the stage FAILs, and the job writes `bank_failures_R<R>.txt`.

## Stage before submitting

1. `GATE_DATA_DIR` (`/work/users/n/u/nur2/desurv-rerun/data`) must hold exactly these two files. Copy them from the canonical `results/`. The job reads nothing else, and no validation cohort.
   - `tar_data_filtered_tcgacptac.rds` (TCGA + CPTAC training data, 1970 genes x 273 samples)
   - `tar_params_best_tcgacptac.rds` (k = 3, alpha, lambda, nu, ntop = 270)

   Their md5 sums are recorded in `settings$input_md5`. Compare them with the local copies.
2. Install DeSurv 1.1.0 (plus `cluster`) in `R_LIBS=/work/users/n/u/nur2/desurv-rerun/lib/r4.5`.
3. Put `gate_3a.R` and `gate_3a.sbatch` in one directory on Longleaf.

## Submit

```bash
cd <dir with gate_3a.R>
sbatch gate_3a.sbatch
```

- Constituent runs are checkpointed in chunks of 50 to `$GATE_OUT_DIR/runs/`.
- If the job is killed, resubmit the same command. Finished chunks are loaded, not refitted.
- Delete `runs/` only if you want to start over.
- Expected runtime is well under 12 h. On one local core a fit takes about 7 s, so even the worst case (600 gate fits + 300 reference fits) is short on 32 cores.

## Outputs (`$GATE_OUT_DIR`)

- `gate3a_summary.txt`: read this first. It contains:
  - the `VERDICT:` line (`PASS at R = 50`, `BORDERLINE at R = 50 -> PASS at R = 100`, or `FAIL ...` with the reason);
  - for each stage, the per-factor median, min and top-gene overlap, plus the full 6-pair table;
  - for each ensemble: mode_pass, beta, training C, tightness, occupancy, outliers, competing-mode medoid r, and the descriptive half-vs-half line;
  - the reference comparison.
- `gate3a_results.rds`: everything. `gate_stages.rds` is also saved on its own.
- `sessionInfo.txt`.

**Factor ids.** Each ensemble is mapped to E1 by Hungarian alignment. Each pair's partner factor comes from that pair's own Hungarian alignment. The column `consistent_with_E1_map` flags any pair whose own alignment disagrees with the E1-induced map.

## Local dry run

```bash
GATE_DRY=TRUE GATE_DATA_DIR=<results dir, read-only> GATE_OUT_DIR=<scratch> SLURM_CPUS_PER_TASK=4 Rscript gate_3a.R
```

The dry run shrinks the job:
- stage 1 = 4 ensembles of 6 runs;
- stage 2 = 4 ensembles of 12 runs;
- 2 half-splits;
- one 4-seed reference set.

Add `GATE_FORCE_BORDERLINE=TRUE` to force the second stage. The script refuses that flag unless `GATE_DRY=TRUE`, and it labels the verdict as forced.
