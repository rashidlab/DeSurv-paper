# Brief: implement rev 3 amendment 1 (rank-stratified BO), 2026-10-02

Read first: docs/plans/2026-10-02-rev3-amendment-1-rank-stratified-bo.md (APPROVED) and
docs/plans/2026-10-02-step4-bo-diagnosis.md. Worktree ~/Downloads/DeSurv-paper-recompute (branch
recompute/consensus-desurv). Package: DeSurv 1.2.0, local lib
/tmp/claude-1000/-home-naimrashid-Downloads-DeSurv-paper-clean/da9c6c81-47b4-4511-8cfc-1c71451b090a/scratchpad/lib-rev3
(source ~/Downloads/DeSurv-clean branch feat/reported-estimator). Do NOT change the package.

## Rules
- Never run a code/NN script without DESURV_RESULTS_DIR pointing at a scratch directory; never write to results/
  or figures/. DESURV_RECOMPUTE defaults to TRUE.
- No git commits, no sbatch, no ssh. Coordinator commits and submits.
- Report mismatches between this brief and the code; do not improvise design. Base R; match code/03 style.
- Report to docs/plans/2026-10-02-amendment1-impl-report.md; return at most 150 words.

## What to build
1. `code/03b_bo_rank.R` (one rank per run; env DESURV_BO_K = 2..12):
   - Reuse code/03's setup verbatim where possible (source the same helpers; BO_BOUNDS, BO_FIXED, BO_COMMON,
     quick/full knobs). Best: factor the shared setup out of code/03 into e.g. R/bo_setup.R sourced by both,
     WITHOUT changing any value code/03 uses (verify by comparing the constructed BO_COMMON/BO_BOUNDS/BO_FIXED
     objects before and after with identical(); report).
   - Inputs: the main BO result `desurv_bo_results_tcgacptac` via `load_precomputed` (the job copies it into its
     results dir), and `tar_data_tcgacptac`.
   - History: rows of the main BO's last run history with status ok and k_grid == K, columns alpha_grid,
     lambda_grid, nu_grid, ntop, mean_cindex, status (+ stage). For K = 3 append two seed rows (stage "seed"):
     published alpha 0.334, lambda 0.349, nu 0.0557, ntop 270, mean_cindex 0.6258; May-best alpha 0.442,
     lambda 0.2410, nu 0.2015, ntop 299, mean_cindex 0.5803. Check how .desurv_bo_prepare_history
     (DeSurv R/desurv_bayesopt.R:617) needs the frame (unit coords, eval ids, required columns).
   - Call DeSurv::desurv_cv_bayesopt directly with the args desurv_cv_bayesopt_refine would pass
     (R/desurv_bo_refine.R:282-298: common args + coarse_control + bo_bounds), i.e. BO_COMMON minus the
     refine-only entries, bo_bounds = BO_BOUNDS without k_grid, bo_fixed = c(BO_FIXED, list(k_grid = K)),
     bo_history = the filtered history, n_init = 10, n_iter = max(0, 40 - n_existing - 10) (n_existing counts
     history rows incl. seeds), candidate_pool 4000, exploration_weight 0.01, seed 123, incumbent and estimator
     as code/03. Quick mode: n_init 2, n_iter 1 (keep quick cheap) and the quick knobs code/03 already uses.
   - Then DeSurv::desurv_bo_finalists(<result>, ..., n_per_k = 3, cv_overrides as code/03). Confirm it picks the 3
     best by this rank's GP posterior mean over the FULL history (inherited + new); report what it does.
   - Save via cache_or_compute: `desurv_bo_results_k{K}_tcgacptac`, `desurv_bo_finalists_k{K}_tcgacptac`.
     Also save training C for the 3 finalists (from finalist fold diagnostics if present; else report).
2. `code/03c_bo_rank_select.R`: load the 11 finalist objects; renumber finalist ids uniquely; rbind configs and
   fold_cindex; apply DeSurv::desurv_select_k_paired_1se(fold_cindex, configs); build an object accepted by the
   existing helpers k_selection_from_finalists() / select_bo_params_finalists() (R/bo_helpers.R) and write
   `desurv_bo_finalists_tcgacptac`, `tar_k_selection_tcgacptac`, `tar_params_best_tcgacptac`,
   `tar_data_filtered_tcgacptac` exactly as code/03's main stage does after its finalists. Also write a table
   `bo_rank_summary_tcgacptac` (per k: n evals, best CV C, finalists' CV C and train C).
3. Slurm: `slurm/run_rev3_bo_rank.sh <K>` and `slurm/submit_rev3_bo_rank.sh` (11 jobs + a select job with
   afterok), modeled on slurm/run_rev3_bo_stage.sh / submit_rev3_bo.sh (immutable per-job dirs under
   $STEP5A_ROOT/k{K}/results, refuse non-empty, copy inputs in: tar_data from base, main BO result from
   /work/users/n/u/nur2/desurv-rerun/out/step4-dd54827-r2/main/results; MD5SUMS; PROVENANCE.txt). Walltime per
   rank from the amendment's evaluation counts and 5-19 min per eval (k 3-12); 32 cores, 64G.
4. Quick-mode test locally (rev3 lib, DESURV_QUICK=TRUE, scratch dirs): run code/03's main stage quick into a
   scratch dir to get a quick main BO, then 03b for K = 2, 3, 4 and 03c on those three (03c must tolerate the
   available ranks via an env list, default 2:12). Report run times, objects written, and that K = 3 imported its
   history plus 2 seeds. Parse-check all edited R; bash -n all shell.
