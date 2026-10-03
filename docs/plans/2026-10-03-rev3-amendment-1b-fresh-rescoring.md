# Rev 3 amendment 1b: cross-rank selection on fresh-partition rescoring (APPROVED, owner, 2026-10-03 00:05)

Amends the selection step of amendment 1 (docs/plans/2026-10-02-rev3-amendment-1-rank-stratified-bo.md). The
per-rank BOs (jobs 3516485, 3529501-10) and their finalist choice are unchanged.

## Problem
Finalist rescoring was meant to correct the winner's curse of picking each rank's best of ~40 noisy CV
evaluations. It does not: `desurv_bo_finalists` re-evaluates on the BO's own partitions and seed banks, so the
rescored C equals the BO value exactly (k = 2: delta_vs_bo 0.0000 for all three finalists). The cross-rank paired
1-SE rule therefore compares selection-biased maxima. Sizes seen: k = 2 best raw 0.6285 vs GP mean 0.6228;
k = 10 single-evaluation peak 0.6430 against neighbouring evaluations 0.602-0.627.

## Rule (primary)
For every rank k = 2..12, take the same 3 finalists per rank (top GP posterior mean, as now) and re-evaluate each
with the identical estimator and settings (`cv_args` stored in the finalists object) on FRESH CV: 3 new repeated
5-fold partitions, partition seeds 20261111-20261113, per-fold seed banks from base 5,000,000 (BO uses 20261101-3
and 4,000,000+; sensitivity reruns 3,000,000+; 20261201-2 are used elsewhere). Apply the unchanged paired 1-SE
rule (`DeSurv::desurv_select_k_paired_1se`) to the fresh fold C-indices across all 33 finalists. The selected k,
its best finalist's parameters, and the runner-up come from this.

## Reported side by side (owner stop)
| | Statistic per rank | SE | Role |
|---|---|---|---|
| B | fresh-partition C of the best finalist | paired, 15 fresh units | PRIMARY selection |
| A | BO-partition C of the best finalist (job 3529511, amendment 1 as approved) | paired, 15 BO units | comparison |
| C | best GP posterior mean | none | descriptive only (no paired SE; GP fits differ by rank) |
Plus per rank: finalists' alpha, lambda, nu, ntop, training C, and B minus A: the change induced by fresh-partition
evaluation. For one rank it mixes optimism (selection bias) with fresh-CV sampling variation; its systematic
direction and size across the 11 ranks is the informative diagnostic.
If A and B select different k, B stands and the disagreement is reported.

## Order of events (for the record)
- 2026-10-02 ~21:30: winner's curse identified on k = 2 (raw 0.6285 vs GP mean 0.6228; rescoring delta 0).
- 2026-10-02 23:55: k = 10 single peak 0.6430 observed in the running log (no parameters, no other rank finished).
- 2026-10-03 00:05: owner approves this amendment. No rank except k = 2 has finished; no cross-rank selection
  under any rule has been computed.

## Implementation
Paper-side only; DeSurv stays 1.2.0 (its `cv_overrides` deliberately refuses partition changes, so the new step
calls `DeSurv::desurv_cv` with the stored `cv_args`, replacing only `cv_seeds` and `seed_bank_base`).
- `code/03d_bo_rank_rescore.R`: one rank per job (DESURV_BO_K) -> `desurv_bo_rescore_k<K>_tcgacptac.rds`
  (configs + fresh fold_cindex), and a select mode that applies the 1-SE rule across ranks and writes the
  ORDINARY downstream names `tar_k_selection_tcgacptac.rds` and `tar_params_best_tcgacptac.rds` into its own
  directory `$STEP5A_ROOT/rescore-select/results`, plus the A/B/C side-by-side table. Both objects carry
  `attr(, "selection_rule") = "fresh_partition_paired_1se"`.

## Which selection is authoritative (amendment to the review fix, owner 2026-10-03)
- The DIRECTORY, not the filename, defines authority. A and B both write unsuffixed names.
  - `$STEP5A_ROOT/rescore-select/results` (B): PRIMARY. Step 5 takes `tar_k_selection` and `tar_params_best`
    from here and nowhere else.
  - `$STEP5A_ROOT/select/results` (A, job 3529511): comparison only; never an input to any later step.
- Safeguard: Step 5 fails immediately unless its selection inputs identify themselves as the fresh rule.
  `code/04_fit_models.R` checks `attr(tar_k_selection, "selection_rule")` and `attr(tar_params_best,
  "selection_rule")` against `DESURV_REQUIRE_SELECTION_RULE` when that variable is set, and
  `slurm/rev3_env.sh` (sourced by every rev 3 job, including whatever runs code/04 for Step 5) and `slurm/submit_rev3_downstream.sh` set it to `fresh_partition_paired_1se` by default. A's objects carry no
  such attribute (or `bo_partition_paired_1se` if 03c is ever rerun), so feeding A to Step 5 stops.
- Longleaf: one rescoring job per rank (afterok on that rank's BO job; k2 can start now), then one select job.
  Separate output directory; job 3529511 (A) runs untouched.
- Cost: 3 finalists x (15 fold fits) per rank, about 1-5 h per rank at the observed pace.
- Quick-mode branch required (1 partition, 2 finalists), per the repo's quick-mode contract.
