# Amendment 1b implementation report (2026-10-03)

Spec: docs/plans/2026-10-03-rev3-amendment-1b-fresh-rescoring.md (APPROVED). Worktree
~/Downloads/DeSurv-paper-recompute at 564f6f1 (plus the uncommitted amendment edit). Nothing committed, no
sbatch or ssh, package untouched (DeSurv 1.2.0 from lib-rev3, `packageVersion == "1.2.0"` checked).
Every Rscript ran with DESURV_QUICK=TRUE and DESURV_RESULTS_DIR under the session scratchpad
(`.../4d6f63f8-.../scratchpad/amend1b/`); `git status --short results figures` is empty after all tests.

## Files

| File | Status | What |
|---|---|---|
| R/bo_rescore.R | new | `rescore_finalists(fin, tar_data, cv_seeds, seed_bank_base)`: the per-finalist call of `desurv_bo_finalists()` (DeSurv R/desurv_bo_finalists.R:176-193) with only those two cv_args replaced; returns configs, long fold C, per-fold seed diagnostics, the cv_args used |
| code/03d_bo_rank_rescore.R | new | Rank mode (DESURV_BO_K) with hard checks (a)-(c); select mode (DESURV_BO_RESCORE_SELECT=1) writing the ordinary selection names, tagged, plus the A/B/C side-by-side table |
| code/03c_bo_rank_select.R | modified | tar_k_selection and tar_params_best get `attr(, "selection_rule") <- "bo_partition_paired_1se"` (+7/-2 lines) |
| code/04_fit_models.R | modified | Guard right after `tar_params_best` is loaded (line 48); see below |
| slurm/submit_rev3_downstream.sh | modified | `export DESURV_REQUIRE_SELECTION_RULE="${...:-fresh_partition_paired_1se}"` and echo |
| slurm/run_rev3_rescore.sh | new | Runner: `<K>` -> `$STEP5A_ROOT/rescore-k<K>/results`, `select` -> `$STEP5A_ROOT/rescore-select/results`; refuses a non-empty results dir; copies inputs; MD5SUMS |
| slurm/submit_rev3_rescore.sh | new | RANKS, DEPS ("K=jobid" -> afterok), SELECT_RANKS; validates everything before the first sbatch; appends PROVENANCE.txt |
| docs/plans/2026-10-03-amendment1b-impl-report.md | new | This report |

R/bo_rescore.R is an addition to the brief: the rescoring call lives in one sourced function so that
the original-seed reproduction test runs the identical code path that 03d uses.

## Design points

- Seeds: full `FRESH_CV_SEEDS = c(20261111L, 20261112L, 20261113L)`, quick `20261111L` (one partition, as
  bo_setup.R's quick branch); `FRESH_SEED_BANK_BASE = 5000000L` in both modes. Quick finalists per rank stay 2
  (inherited from 03b's FINALISTS_PER_K).
- No cv_overrides. parallel_grid/ncores_grid come from the stored cv_args, i.e. the 03b job's (32 cores on
  Longleaf). desurv_bo_finalists() passed the same values as overrides in that job, so the call is the same.
  The rescore job is also 32 cpus.
- 03d rank mode also checks the tar_data fingerprint against `desurv_bo_results_k<K>$data_fingerprint`. That is
  why the rank job copies the BO results, and why it refuses fresh seeds or a bank that coincide with the BO's.
- Select mode pools as 03c does: finalist_id is renumbered 1..n across ranks (rank order). desurv_select_k_paired_1se()
  keys finalists and paired units by finalist_id, so the ids must be unique (confirmed in the package source).
  A and B use the same ids. B's pooled object sets `mean_cindex` to the fresh value and keeps
  `mean_cindex_bo_partition`.
- Select mode writes the same object set as 03c in its own directory: `desurv_bo_finalists_tcgacptac` (pooled,
  FRESH fold C), `tar_k_selection_tcgacptac`, `tar_params_best_tcgacptac` (both tagged
  `fresh_partition_paired_1se`), `tar_data_filtered_tcgacptac`, and `bo_rank_selection_side_by_side_tcgacptac`
  (.rds and .csv). It also writes `bo_rank_selection_side_by_side_finalists_tcgacptac.csv`, with all finalists'
  A C, B C, B minus A, parameters and train C. The rds carries attributes `finalists`, `selection`
  (k chosen by B, A and C) and `note`. tar_data_filtered does not depend on the rule (ngene is fixed).
- Side-by-side per k: B_* (best fresh finalist: rank_finalist_id, eval_id, cv, mean_diff, se_diff, n_pairs,
  within_1se, selected, alpha/lambda/nu/ntop, train_cindex_bo), the same A_* columns,
  C_best_gp_mean / C_top, B_minus_A_best (best under B minus best under A; the two may be different finalists),
  B_minus_A_same_finalist (B's best finalist, fresh minus BO-partition), and selected_by.
- Hard check (a) is stricter than the brief. Every diagnostic `cv_seed` must be in the fresh set and the fresh
  set must be fully covered. Every `seed_first` must be >= 5,000,000 and also equal
  `desurv_cv_seed_bank(rep, fold, nfolds, n_seed_fits, base = 5e6)[1]`. Non-ok fold status is counted in
  the log but does not stop the job, matching desurv_bo_finalists().

## Quick-mode test

Inputs: the amendment-1 quick main BO (`.../885b7638-.../scratchpad/amend1/qmain/quick`: tar_data and the
quick main desurv_bo_results). R/bo_setup.R, code/03 and code/03b have not changed since e6080ff, so I reused
it rather than rerunning the main BO. Then, as in the amendment-1 report: 03b K=2,3,4 (all exit 0) -> 03c
ranks 2,3,4 (exit 0) -> 03d rank K=2,3,4 (exit 0; 8.5 s, 12.3 s, 19.9 s) -> 03d select ranks 2,3,4 (exit 0).
Script: `scratchpad/amend1b/run_quick.sh`, logs alongside.

### Hard checks, rank mode (actual log lines; identical for K = 2, 3, 4)

```
check (a) PASS: 10 diagnostic rows; cv_seed in {20261111}; seed_first 5000001..5040001 (all = bank start); status not ok: 0
check (b) PASS: cv_args identical to the stored ones apart from cv_seeds and seed_bank_base
check (c) PASS: 5/5 fold rows per finalist (expected 5)
```
(Quick: 2 finalists x 1 partition x 5 folds = 10 diagnostic rows. Full mode requires 15 per finalist.)

Negative tests: a doctored cached rescore object was loaded through 03d with DESURV_RECOMPUTE=FALSE. Each one stops:
- one fold's seed_first set to 4,020,001: `HARD CHECK (a) FAILED`
- cv_args$tol x 10: (a) PASS, `HARD CHECK (b) FAILED`
- one fold row dropped: (a), (b) PASS, `HARD CHECK (c) FAILED: fold rows per finalist 4/5, expected 5.`
- DESURV_BO_K and DESURV_BO_RESCORE_SELECT both set: stops with "Set either ... not both."

### Fresh vs BO-partition C (quick; per finalist)

| k | finalist (eval) | BO-partition C | fresh C | fresh minus BO | train C (BO) | train C (fresh) |
|---|---|---|---|---|---|---|
| 2 | 1 (2) | 0.6055 | 0.5764 | -0.0291 | 0.723 | 0.767 |
| 2 | 2 (4) | 0.5808 | 0.5666 | -0.0142 | 0.753 | 0.727 |
| 3 | 1 (3, published seed) | 0.5905 | 0.6226 | +0.0321 | NA | 0.715 |
| 3 | 2 (5) | 0.6088 | 0.6025 | -0.0062 | 0.675 | 0.736 |
| 4 | 1 (3) | 0.6147 | 0.5889 | -0.0258 | 0.694 | 0.819 |
| 4 | 2 (5) | 0.5776 | 0.5730 | -0.0046 | 0.878 | 0.928 |

The k = 3 seed's BO-partition value is its quick-recipe rescore (0.5905), not the full-recipe 0.6258 (same caveat as
the amendment-1 report). Quick numbers say nothing about the full run.

### Select mode (quick, ranks 2,3,4)

| k | B cv | B diff / SE | B within | A cv | A diff / SE | A within | C best GP mean | B-A best | B-A same finalist |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.5764 | 0.0462 / 0.0170 | no | 0.6055 | 0.0092 / 0.0057 | no | 0.6055 | -0.0291 | -0.0291 |
| 3 | 0.6226 | 0 / 0 | yes, **selected** | 0.6088 | 0.0060 / 0.0119 | yes, **selected** | 0.6258 (top) | +0.0138 | +0.0321 |
| 4 | 0.5889 | 0.0337 / 0.0262 | no | 0.6147 | 0 / 0 | yes (best) | 0.6147 | -0.0258 | -0.0258 |

B: k_selected 3, k_best 3, runner-up 4. A: k_selected 3, k_best 4, runner-up 4. C: top 3. B's tar_params_best:
k 3L, alpha 0.334, lambda 0.349, nu 0.0557, ntop 270L, with `identical(tar_params_best$k, k_selected)` TRUE.
B's tar_k_selection and tar_params_best have the same names as 03c's.
**A as computed in 03d (from the finalists' own fold_cindex) equals 03c's pooled selection per_k**
(`all.equal` TRUE, k_selected 3 in both). Attributes: 03c objects `bo_partition_paired_1se`, 03d objects
`fresh_partition_paired_1se`.

### Original-seed reproduction (scratchpad/amend1b/test_repro.R)

The same `rescore_finalists()` was given the finalists object's ORIGINAL cv_seeds (20261101) and bank (4,000,000):

```
k=2 finalist 1 (eval 2): stored mean_cindex 0.6055005323 | original-seed rerun 0.6055005323 | identical TRUE | fold C identical TRUE | args identical to stored TRUE
   original-seed diag cv_seed {20261101}, seed_first 4000001..4040001
   fresh mean_cindex 0.5763682406 (differs: TRUE; diff -0.0291)
k=4 finalist 1 (eval 3): stored mean_cindex 0.6147383132 | original-seed rerun 0.6147383132 | identical TRUE | fold C identical TRUE | args identical to stored TRUE
   original-seed diag cv_seed {20261101}, seed_first 4000001..4040001
   fresh mean_cindex 0.5889055522 (differs: TRUE; diff -0.0258)
```
The reproduction is exact under `identical()`, for the mean and for all 5 fold C values. So the call path is equivalent to
desurv_bo_finalists() apart from the seeds, and the fresh values differ.

### Guard in code/04 (actual code/04 run, scratch results dirs holding only tar_data_filtered,
tar_params_best and, where stated, tar_k_selection; DESURV_RECOMPUTE=FALSE)

| Inputs | DESURV_REQUIRE_SELECTION_RULE | Result |
|---|---|---|
| B (03d) both objects | fresh_partition_paired_1se | both checks print the rule, then the next load fails (`fit_std_tcgacptac` not found), so the guard PASSED |
| A (03c, tagged) | fresh_partition_paired_1se | `Error: tar_params_best_tcgacptac has selection_rule 'bo_partition_paired_1se', but ... requires 'fresh_partition_paired_1se'.` |
| A from the amendment-1 test (untagged, like job 3529511's) | same | `... has selection_rule '<none>' ...` stop |
| B tar_params_best + A tar_k_selection | same | tar_params_best passes; `tar_k_selection_tcgacptac has selection_rule 'bo_partition_paired_1se'` stop |
| B tar_params_best only | same | passes; logs "tar_k_selection_tcgacptac not in ...; not checked" |
| A (03c) | unset | no check; proceeds to the fit_std load (no behavior change) |

code/04 loads only tar_params_best. The guard checks it, then checks tar_k_selection only if
`precomputed_path("tar_k_selection_tcgacptac")` exists, as the brief allows.

### Slurm scripts (bash -n OK; dry runs with a fake `sbatch` on PATH, scratch STEP5A_ROOT)

- submit_rev3_rescore.sh, RANKS=2:4, DEPS="3=3529501 4=3529502": k2 `--time=04:00:00` (no dep, finalists present),
  k3 `06:00:00 --dependency=afterok:3529501`, k4 `06:00:00 --dependency=afterok:3529502`, select
  `--cpus-per-task=4 --mem=16G --time=01:00:00 --dependency=afterok:<k2>:<k3>:<k4>` with DESURV_BO_RANKS=2:4 in its
  environment. PROVENANCE.txt was appended.
- With RANKS=2,3 and no DEPS entry or finalists for k3, it aborts before any sbatch (0 calls). With SELECT_RANKS
  covering a rank that was neither rescored here nor finished, it aborts.
- submit_rev3_downstream.sh: prints `DESURV_REQUIRE_SELECTION_RULE=fresh_partition_paired_1se`, and the fake sbatch
  sees the variable in each job's environment. A preset value (`xyz`) is respected.
- run_rev3_rescore.sh was not executed (it needs `module` and Longleaf paths). Its md5 line was checked under
  `set -euo pipefail` with no csv present (rank mode) and does not fail.

## Spec mismatches and open points

1. **submit_rev3_downstream.sh does not run code/04.** It submits code/14-29 only. No rev 3 submitter for
   code/04 exists yet in slurm/ (run_fair_comparator.sh runs 03/04/05 for its own purpose). The export works and
   reaches every job it submits, but none of those jobs contains the guard. The Step 5 submission that runs
   code/04 must set the variable itself (or source the same default), or the safeguard never fires. I made no change
   beyond the brief.
2. Locally, `submit_rev3_downstream.sh` fails at its first `cd` when CDPATH is set (this machine's profile
   sets it: `cd` echoes the path into `$(...)`). This predates the change and is harmless on Longleaf if CDPATH is unset there.
   The dry run used `env -u CDPATH`. Not changed.
3. The guard checks tar_k_selection only if present (as briefed). If Step 5's directory lacks it, only
   tar_params_best is checked.
4. B_minus_A has two forms: best-vs-best, which may compare different finalists, and the same configuration for
   B's best. The amendment's "B minus A per rank" is ambiguous between them, so both are reported.
5. Walltimes are the brief's values. The rescore job's cores equal the stored ncores_grid only if the rank BO also
   ran on 32 cpus (it did, per run_rev3_bo_rank.sh).
6. No other mismatches: cv_args carries cv_seeds and seed_bank_base, desurv_cv diagnostics expose cv_seed and
   seed_first, and the finalists object has configs, fold_cindex, param_names and cv_args.
