# Track 3a gate script: report (2026-09-30)

Written, not run on Longleaf. Nothing committed. Only `analysis/gate/` was created in this worktree.

## Files

- `analysis/gate/gate_3a.R`: the probe.
- `analysis/gate/gate_3a.sbatch`: Longleaf job (general, 1 node, 32 cpus, 64G, 12 h, `r/4.5.0`, paths as briefed; BLAS threads pinned to 1).
- `analysis/gate/README.md`: what to stage, how to submit, what to read.

## What the script does

1. **Preflight.** Checks `packageVersion("DeSurv") == "1.1.0"` and that `cox_stats_cpp` exists. Then it reads the env vars, refuses any `GATE_OUT_DIR` inside a `results/` directory or inside `GATE_DATA_DIR`, and writes `sessionInfo.txt`.
2. **Inputs.** Reads only `tar_data_filtered_tcgacptac.rds` and `tar_params_best_tcgacptac.rds` and records their md5 sums.
   - Hyperparameters: lambdaW and lambdaH default to 0 (code/04:55-61). ntop = 270.
   - Fit settings: the code/04 full-mode ones: tol 1e-7, tol_init 1e-5, maxit = imaxit = 4000, ninit 1.
3. **Gate bank.** Seeds 1,000,001-1,000,400, fitted in order with the exact call `desurv_consensus_fit` uses.
   - Runs are saved in atomic chunks of 50 to `runs/`, and a resubmitted job loads finished chunks instead of refitting them.
   - A worker crash is retried once, sequentially.
   - A genuinely failed seed stops the job and writes `bank_failures.txt`. No run is dropped.
   - The `data` element (a copy of X, about 4 MB) is removed from each stored fit. Neither estimator reads it.
4. **Consensus fits.** `desurv_consensus_fit(run_fits = prefix)` at R = 25, 50, 100, 200, 400. Each object is saved.
5. **Gate.** `desurv_consensus_gate` on R = 50, 100, 200, 400 (see Blocker 1), with split_seed 20260930, 10 splits and the package-default thresholds.
6. **Reference.** The code/04 recipe on the three 100-seed sets:
   - `run_seed_fits` semantics: failed fits are dropped and the count is recorded.
   - `desurv_consensus_seed` with ntop and `min_frequency = ceil(0.3 n)`.
   - One final `desurv_fit`.
   - Pairwise W correlations of the three endpoints, with factors aligned by `DeSurv:::.dc_align` (Hungarian on Pearson r), reported as Pearson and Spearman.
7. **Outputs.** `gate3a_results.rds` (everything), `gate.rds`, `consensus_fits.rds` and `gate3a_summary.txt`. The summary gives, per R:
   - pass/fail;
   - per-factor median and p05 r;
   - mode_pass;
   - occupancy;
   - tightness;
   - n_outliers;
   - beta and training C agreement between halves.

   It also includes the reference correlations.

## Local dry run (correctness only)

Run with `GATE_DRY=TRUE` on 4 cores, using the Track 1 library and the paper-clean `results/` read-only.
- It covered 12 bank seeds, R 4 and 6, 2 splits and one 4-seed reference set.
- Exit 0 in 0.7 min. Each fit took about 7 s per core.
- A second invocation loaded every checkpoint and reproduced the gate, the consensus W and the reference W identically.
- `git status results/` in paper-clean is clean.

At the full size, the whole job should take well under 1 h on 32 cores.

The dry-run numbers are meaningless: at R ≤ 6 every factor fails and one half emptied a family. They show that the error path is recorded as a FAIL, per clarification 4.

## Blocker (needs an owner ruling before submission)

1. **R = 25 cannot be gated.** Section 2 lists R = 25, but `desurv_consensus_gate` requires even R, and section 9 defines halves as "the first and second R/2".
   - What I did: the script gates 50/100/200/400 and fits R = 25 only as a descriptive consensus.
   - Options for the owner: accept this, substitute 24 or 26 (`GATE_R_VALUES=24,50,100,200,400`, no code change), or allow unequal halves (an API change).
   - It matters only if R = 25 would have been the smallest passing R.

## Minor observations (no action taken)

- When a split errors, the gate's per-factor median is NA (it has no `na.rm`), while p05 uses `na.rm = TRUE`. That R already fails, so the gate outcome is unaffected. Only the reported median is uninformative.
- In the dry run the constituent fits converged in a median of 37 recorded iterations at tol 1e-7, as the code/04 recipe does.
- The reference recipe runs under 1.1.0, whose top-gene tie-break (FP15) differs from 1.0.1 only on exact ties at ntop.
