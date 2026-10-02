# Amendment 1 implementation report (2026-10-02)

Brief: docs/plans/2026-10-02-amendment1-impl-brief.md. Worktree ~/Downloads/DeSurv-paper-recompute,
branch recompute/consensus-desurv, base 8a587d0. Nothing committed, no sbatch or ssh, package untouched
(DeSurv 1.2.0 from lib-rev3). No writes to results/ or figures/ (`git status --short results figures`: 0 lines).
Every code/NN run used DESURV_RESULTS_DIR under the session scratchpad (`.../scratchpad/amend1/`).

## Files

| File | Status | What |
|---|---|---|
| R/bo_setup.R | new | Shared BO setup moved verbatim from code/03 lines 28-139 (knobs, BO_BOUNDS, BO_FIXED, BO_COMMON, run_bo_finalists(), bo_ngene()) plus a header |
| code/03_bayesian_optimization.R | modified | Lines 28-139 replaced by `source("R/bo_setup.R")` (3 added, 110 removed; nothing else changed) |
| code/03b_bo_rank.R | new | One fixed-k BO per run (DESURV_BO_K), then its finalists |
| code/03c_bo_rank_select.R | new | Pools the per-rank finalists, paired 1-SE, writes the code/03 main-stage selection objects and the summary |
| slurm/run_rev3_bo_rank.sh | new | Runner: `<K>` (2..12) or `select` |
| slurm/submit_rev3_bo_rank.sh | new | 11 rank jobs + select job (afterok on all 11), PROVENANCE.txt |
| docs/plans/2026-10-02-amendment1-impl-report.md | new | This report |

## Verification

**Shared setup unchanged.** HEAD's code/03 lines 28-139 and R/bo_setup.R were each evaluated in four
configurations (quick TRUE/FALSE x ncores 1/32). All 23 objects matched under `identical()` in every
configuration (92/92): the 21 knobs including BO_BOUNDS, BO_FIXED, BO_COMMON, plus the deparsed bodies of
run_bo_finalists and bo_ngene. The whole captured list was also `identical()`. A `diff` of HEAD
code/03 lines 28-139 against the R/bo_setup.R body is empty.

**Parse/syntax.** `parse()` passes on R/bo_setup.R, code/03, code/03b and code/03c. `bash -n` passes on both
new shell scripts.

**Quick-mode test** (rev3 lib, DESURV_QUICK=TRUE, DESURV_NCORES=1, scratch dirs):

| Run | Wall time | Files in results dir (incl. copied inputs) |
|---|---|---|
| code/03 main stage (quick main BO, 8 evals, k 2-6) | 123 s | 6 |
| 03b K=2 | 24 s | 4 (2 written) |
| 03b K=3 | 26 s | 4 (2 written) |
| 03b K=4 | 43 s | 4 (2 written) |
| 03c on ranks 2,3,4 (DESURV_BO_RANKS=2,3,4) | 2 s | 12 (5 written) |

- History import, counted. K=2: 1 inherited + 2 LHS + 2 EI = 5 rows. K=3: **2 inherited + 2 seeds** + 2 LHS + 1 EI = 7 rows. K=4: 2 + 2 + 1 = 5 rows.
  For every K the first n_existing rows of the rank-BO history equal the supplied frame (parameters and
  mean_cindex, all.equal TRUE). For K=3, eval_id 3-4 have stage "seed" with the exact seed values and CV C
  (0.6258, 0.5803). The package-side import count (`n_loaded`) equalled nrow(history) for every rank.
- Scoring arguments. The rank BO's cv_args are `identical()` to the main BO's, ignoring k_grid,
  parallel_grid and ncores_grid. 03b checks this at run time and stops on a mismatch. bo_fixed k_grid equals K.
- Finalists. In every rank, the finalists equal the top n_per_k by history$gp_mean over the full history
  (inherited + seeds + new). gp_mean was finite for every row: 5/5, 7/7, 5/5. At K=3 the **published seed was
  picked as a finalist** (expected in quick mode). Its rescore delta was -0.035, because the seed CV C comes
  from the full recipe and the rescore used the quick recipe. In the full run this delta is a reproduction
  check of 0.6258. For new evaluations delta_vs_bo = 0 exactly, so rescoring is deterministic.
- 03c pooled 6 finalists with unique ids 1-6 and 30 fold rows (6 x 5). param_names are
  k_grid, alpha_grid, lambda_grid, nu_grid, ntop. Selection: k_selected 3, k_best 4, runner-up 4.
  tar_params_best has k = 3, alpha, lambda, nu, ntop (int), and `identical(tar_params_best$k, k_selected)` holds.
  tar_data_filtered was written, and bo_rank_summary was written as list(per_k, finalists, note).

## What desurv_bo_finalists does on a rank BO (confirmed by reading DeSurv R/desurv_bo_finalists.R:115-227 and the test)

It takes the run's ok rows, ranks them by `history$gp_mean`, and keeps the top n_per_k. gp_mean is the latent
posterior mean of the BO's final GP, fitted on all ok rows (inherited, seeds and new). The ties-break order is
higher observed mean, then earlier eval_id. With k not tuned, k is read from `cv_args$k_grid`. Each finalist is
re-scored with the BO's cv_args (partitions, seed bank, estimator) plus parallel_grid/ncores_grid overrides
from code/03. If the final GP fit fails, the function refits it; if that also fails, it stops.

## Mismatches and decisions to review

1. **Training C is not in the finalist object.** desurv_bo_finalists() keeps only held-out fold C. 03b
   therefore adds `configs$train_cindex_bo`, the mean train_cindex over partition x fold from the BO
   evaluation of the same configuration. New evaluations use the rank-BO diagnostics. Inherited ones use the
   main-BO diagnostics, mapped through the main eval_id. The mapping relies on the package renumbering
   imported rows 1..n in order, and 03b stops if any row would be dropped. The two k = 3 seeds have NA:
   they were evaluated outside both BOs, and the diagnosis note has their train C from job 3500640 (0.80 and
   0.91), which is not loaded. This matches the finalist's own train C only because rescoring is
   deterministic (delta_vs_bo = 0 above). The inherited-row lookup returned the right values in a direct
   check, but no inherited main-BO row became a finalist in the quick test, so that path was not exercised
   end to end.
2. **Quick-mode budget deviates from the brief (quick only).** With n_init 2 and n_iter 1, K=2 had 4
   points in 4 dimensions. The package's final GP fit failed ("Final GP fit failed; incumbent falls back"),
   and desurv_bo_finalists() then stopped with "No GP posterior mean". Quick mode now uses
   n_iter = max(1, 5 - n_existing - 2), which gives K=2 two EI iterations. Full mode is unchanged:
   n_init 10, n_iter = max(0, 40 - n_existing - 10).
3. **Same LHS design in every rank.** With seed 123 (as code/03) and a 4-dimensional search, every rank draws
   identical 10 LHS points in (alpha, lambda, nu, ntop). This is visible in the quick test, where k = 2, 3
   and 4 share the same new configurations. It follows from the brief, but the owner may want a per-rank
   seed such as 123 + K. Not changed.
4. **n_existing.** The brief's "history rows incl. seeds" is used. 03b stops rather than proceed if the
   package would import fewer rows (out of bounds or duplicate keys), because the budget rule and train-C
   mapping depend on it. Before submitting, the coordinator may want to check per-rank counts on the real
   main BO (4/6+2/6/7/18/29/17/22/19/13/9 expected). The 03b log prints them on its first lines.
5. **Inputs path for tar_data.** The brief says "from base". The runner uses
   `$STEP4_ROOT/base/results/tar_data_tcgacptac.rds` with STEP4_ROOT defaulting to
   `/work/users/n/u/nur2/desurv-rerun/out/step4-dd54827-r2`, and `main/results/desurv_bo_results_tcgacptac.rds`
   for the main BO. 03b also checks the tar_data fingerprint against the main BO's data_fingerprint.
6. **03c inputs.** The summary needs evaluation counts, so 03c loads the 11 `desurv_bo_results_k{K}` objects
   as well as the 11 finalists. The select job copies both from each `k{K}/results`.
7. **Pooled object shape.** select_bo_params_finalists() takes hyperparameters from `param_names`. The rank
   finalists lack k_grid, so 03c adds `configs$k_grid = k` and prepends "k_grid" to param_names, matching the
   main BO's finalists. eval_id stays rank-local (unique only within k). It is used only as a within-k
   tie-break, and `rank_finalist_id` keeps the original id. cv_args in the pooled object drop k_grid. 03c checks
   that all ranks' cv_args are identical apart from k_grid, parallel_grid and ncores_grid.
8. **Walltime.** Each rank gets (planned new evaluations + 3 finalist rescorings) x min/eval, using 5 at k=2-3
   rising linearly to 19 at k=12. That is doubled, rounded up, with a 12 h floor: k2-9 12 h, k10 13 h, k11 18 h,
   k12 22 h. The select job gets 4 cores, 16G and 1 h. Note that the diagnosis measured 12.8 min for one
   k = 3 evaluation at an overfitting setting, so per-evaluation time varies by more than k alone. At that
   rate k = 3 would need about 7.5 h, which fits within 12 h.
