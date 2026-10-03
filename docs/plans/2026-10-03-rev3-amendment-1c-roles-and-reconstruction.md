# Rev 3 amendment 1c: analysis roles, published-procedure reconstruction, k = 3 ladder (APPROVED, owner, 2026-10-03)

Supersedes the role statements of amendment 1b (B "primary", B's directory authoritative for Step 5, the
global DESURV_REQUIRE_SELECTION_RULE default). The rescoring design of 1b (fresh partitions 20261111-3, seed bank
5,000,000, paired 1-SE, hard checks) is unchanged.

Fixed BEFORE any B (fresh-rescoring) result exists and before ranks k5, k11, k12 finish (jobs 3529503, 3529509,
3529510 running at the time of this commit; no cross-rank selection under any rule computed on rev 3).

## 1. What the published selection did (code trace, public:code/03:111-122, public:R/bo_helpers.R)
Cached May BO (results/desurv_bo_results_tcgacptac.rds, 150 evaluations over k = 2..12):
1. Rank rule (`select_bo_k_by_cv_se`): best raw CV C per k; overall best k = 7 (eval 129, C 0.6545); threshold =
   0.6545 minus that evaluation's unpaired fold SE (0.011) = 0.644; smallest k whose best evaluation clears it:
   k = 3 (eval 136, C 0.6462; alpha 0.442, lambda 0.241, nu 0.201, ntop 299).
2. Parameter extraction: `tar_params_best` = the OVERALL best evaluation's hyperparameters (eval 129, k = 7:
   alpha 0.334, lambda 0.349, nu 0.056, ntop 270) with k replaced by 3.
So the published configuration is k = 7's optimum run at k = 3. It was never cross-validated at k = 3 in May; the
0.646 that admitted k = 3 belongs to eval 136. Whether step 2 was an intended convention (k as structural
complexity, other hyperparameters as transferable nuisance parameters) or an implementation artifact (global-best
row instead of best row conditional on the selected k) is undocumented; describe it as the latter unless
documentation is found: "the 1-SE procedure selected k = 3, but the parameter-extraction step used the overall
BO optimum rather than the optimum conditional on k = 3."

## 2. Roles (no analysis may overwrite another's)
| Analysis | Question | Role |
|---|---|---|
| Historical May rule | Did the published procedure select k = 3? (yes) | validation / reproduction |
| Step 5 fit | The prespecified published configuration (k 3, alpha 0.3339356, lambda 0.3494487, nu 0.05565122, ntop 270) re-fit under the rev 3 estimator | MAIN model-reproduction analysis |
| B (amendment 1b, fresh-partition rescoring + paired 1-SE) | What rank does fresh-partition evaluation support? | prospectively specified robustness / sensitivity analysis of rank selection |
| A (job 3529511, BO-partition paired 1-SE) | same, without fresh partitions | comparison |
The Step 5 object is called "the prespecified published model configuration re-fit under the rev 3 estimator",
not "the selected rev 3 model". B selecting a different k (e.g. 6) does not change Step 5; it is a finding about
the sensitivity of rank choice to the modernized estimator and validation procedure.

Mechanics:
- B writes `tar_k_selection_fresh_tcgacptac` / `tar_params_best_fresh_tcgacptac` (tag selection_rule
  "fresh_partition_paired_1se") in `$STEP5A_ROOT/rescore-select/results`. Not a Step 5 input.
- Step 5's `tar_params_best_tcgacptac` is written by `code/03e_published_config.R` from the constants above
  (full precision of May eval 129; value-identical to the published results/tar_params_best, checked) with
  `attr(, "model_role") = "published_configuration_k3"` and `attr(, "source") =
  "published_May_global_BO_optimum_with_k_replaced_by_3"`. code/04 reads no tar_k_selection, so no k object
  is written; 03e refuses to run if tar_params_best or tar_k_selection already exists in the directory.
- `slurm/run_rev3_step5.sh` is the only rev 3 path to code/04: it refuses a directory holding any
  selection or rescoring artifact (tar_params_best*, tar_k_selection*, desurv_bo_finalists*,
  desurv_bo_rescore*, bo_rank_selection_*), exports DESURV_REQUIRE_MODEL_ROLE, runs 03e then 04 in the same job.
  `slurm/run_rev3_step.sh` and `slurm/run_main_pipeline.sh` refuse 03e and 04; `slurm/submit_all.sh` (whose
  chain starts with 04) refuses outright. Staging code/04's other inputs belongs to the Step 5 plan.
- Guard: code/04 stops unless `attr(tar_params_best, "model_role")` equals DESURV_REQUIRE_MODEL_ROLE when that
  variable is set. ONLY the Step 5 submission sets it (to "published_configuration_k3"); it is not in
  rev3_env.sh or any other submitter, so other flows (run_main_pipeline.sh, submit_all.sh) are unaffected. A and
  B objects carry no model_role and cannot reach Step 5.

## 3. Rank-selection comparison (computed from existing outputs after k12 finishes; no new BO, no reruns)
| Method | Rank rule | Non-k hyperparameters | Purpose |
|---|---|---|---|
| May historical | published 1-SE (raw per-k best, unpaired SE of the overall best eval) | global-best row, k replaced | what produced the published model (re-run `select_bo_k_by_cv_se` on the cached May BO; must return k 3, threshold 0.644) |
| Rev 3 published-procedure reconstruction | same rule on the pooled rev 3 rank-BO histories (k2-k12) | global-best rev 3 row, k replaced | direct procedural replication |
| Rev 3 conditional-k | same rule | best rev 3 row within the selected k | effect of the parameter-extraction step |
| B (fresh) | fresh partitions + paired 1-SE | fresh-rescored finalists | modern robustness analysis |
A and best GP mean per rank are reported alongside. Reading:
- reconstruction selects k 3 -> the conflict with B comes from the newer rule;
- reconstruction selects another k (e.g. 2 or 6) -> the change comes from the rev 3 performance surface;
- reconstruction and conditional-k differ materially -> the parameter transplant was consequential.
Caveat to state: the May BO was one search over k (uneven evaluations per k, 5-26); rev 3 pools 11 fixed-k
searches of 40 evaluations each (including evaluations inherited from the step 4 main BO).
The reconstruction's unpaired SE must be computed exactly as `compute_bo_eval_se` did (public:R/bo_helpers.R).

## 4. k = 3 reconciliation ladder (no BO; 3 configurations)
Configurations: (i) eval 136 at k 3 (the point that admitted k 3; May 0.6462 -> rev 3 ~0.585); (ii) the published
hyperparameters at k 3 (May value never computed; rev 3 0.6258); (iii) the published hyperparameters at k 7
(eval 129, May 0.6545; control with a known May value).
Steps, one ingredient at a time: May recipe (mean_starts 30 starts, tol 1e-5, DeSurv 1.0.1, May partitions if
recoverable from the cached cv_args) -> rev 3 partitions/seed banks -> DeSurv 1.2.0 (tie fix) -> tol 1e-7 ->
reported estimator (100 seed fits -> consensus -> final). Also run the last two in the other order (the step 4
diagnosis found a tolerance x estimator interaction). Questions: (A) why eval 136 fell from 0.646 to ~0.585;
(B) how the published hyperparameters behave at k 3 under the May recipe. Do not expand beyond these three
configurations unless the ladder fails to explain the discrepancy.
