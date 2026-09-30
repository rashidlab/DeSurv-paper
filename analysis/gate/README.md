# Track 3a convergence probe

Plan rev 2.1 (frozen at 960b76d), section 2 Track 3, section 3 rules, section 9 seeds, and the
2026-09-30 section 3 clarifications. Requires DeSurv 1.1.0 (`feat/consensus-estimator`, b42d7a9); the
script stops unless `packageVersion("DeSurv") == "1.1.0"` and `DeSurv:::cox_stats_cpp` exists.

## Stage before submitting

1. `GATE_DATA_DIR` (`/work/users/n/u/nur2/desurv-rerun/data`) must hold exactly these two files, copied
   from the canonical `results/` (read-only; nothing else is read, no validation cohort):
   - `tar_data_filtered_tcgacptac.rds` (TCGA + CPTAC training data, 1970 genes x 273 samples)
   - `tar_params_best_tcgacptac.rds` (k = 3, alpha, lambda, nu, ntop = 270)
   Their md5 sums are recorded in the output (`settings$input_md5`); compare with the local copies.
2. DeSurv 1.1.0 installed in `R_LIBS=/work/users/n/u/nur2/desurv-rerun/lib/r4.5` (plus `cluster`).
3. `gate_3a.R` and `gate_3a.sbatch` in one directory on Longleaf.

## Submit

```bash
cd <dir with gate_3a.R>
sbatch gate_3a.sbatch
```

Constituent runs are checkpointed in chunks of 50 to `$GATE_OUT_DIR/runs/`. If the job is killed,
resubmit the same command: finished chunks are loaded, not refitted. Delete `runs/` only to start over.
The local dry run (12 + 4 fits, all stages) finished in a few minutes on 4 cores, so the 12 h limit is
generous.

## Outputs (`$GATE_OUT_DIR`)

- `gate3a_summary.txt`: read this first. Per gated R: PASS/FAIL, per-factor median and p05 half-vs-half
  r, mode_pass, errors, median |delta training C| and per-factor median |delta beta| between halves
  (descriptive). Then per R the consensus fit's mode_pass, beta, training C, tightness, n_members,
  n_outliers, occupancy (runs with 0/1/>1 factors per family), competing-mode medoid r; then the
  code/04 reference: each set's final fit and the pairwise aligned W correlations (Pearson, Spearman).
- `gate3a_results.rds`: everything (gate object, consensus fits, reference fits, bank diagnostics,
  settings, sessionInfo). `gate.rds` and `consensus_fits.rds` are also saved on their own.
- `sessionInfo.txt`, `bank_failures.txt` (only if a gate-bank run failed; the job then stops).

The gate result is `converged_R` in the summary header: the smallest R that passes and whose 2R also
passes. `NA` means no converged R in the tested range (plan: STOP, section 7).

## Open point: R = 25

`desurv_consensus_gate` requires even R (halves are the first and second R/2 of each permutation, as
in section 9), so R = 25 from section 2 cannot be gated as written. The script gates
R = 50, 100, 200, 400 and fits the R = 25 consensus descriptively only. This needs an owner ruling
before submission. If the ruling substitutes an even R (e.g. 24 or 26), set
`export GATE_R_VALUES=24,50,100,200,400` in the sbatch file; no code change is needed.
