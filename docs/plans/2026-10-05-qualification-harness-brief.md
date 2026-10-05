# Brief: qualification harness for the frozen redesign (2026-10-05)

Authority: FROZEN spec `docs/plans/2026-10-05-estimator-redesign-spec.md` (revision 4, ruling 14, tag
`spec-redesign-rev4-frozen`) section 9 and FROZEN `docs/plans/2026-10-05-qualification-prespec.md` (all sections).
The estimator is `analysis/estimator-gate/redesign/theta_lib.R` (do not change its mathematics). Data:
`analysis/estimator-gate/qualification/sim_qual.R` (frozen; do not edit). Where this brief and the frozen files
disagree, the frozen files win: stop on that item and report the mismatch. Do not improvise. Do not commit.

## Deliverables (all under `analysis/estimator-gate/qualification/`)

1. `qual_common.R`: sources the frozen files; checks the git blob hashes listed in prespec section 3.3 ("Frozen
   mechanics") and the hashes of theta_lib.R and kkt_oracle_theta.R recorded at job start; stops on mismatch; seed
   functions of prespec section 4 as code (one function per row of the seed table) with an assertion that the
   namespaces are disjoint for the full job list.
2. `qual_grid_job.R` (args: scenario, replicate, g): one grid or reference fit (prespec section 5 table, rows g = 1-12,
   replicates 1-5 of S1, S2, S3). 30 starts, recording on, keep per-start endpoints and checkpoints. Saves one `.rds`
   per job under `out/grid/`. Contents: per-start status, reason, sweeps, F, residuals, F-increase count, column
   statistics, checkpoints; per-start test C (converged starts, own v~ and G, gate-spec section 5 scorer); test C of
   the returned start's checkpoint iterates; s0, lambda_1; the frozen-file hashes; the realized SD-ratio fields from
   sim_qual for S3/S3b.
3. `qual_tune_job.R` (args: scenario, replicate): one tuning run (prespec section 5 "Tuning run"): 5-fold inner split
   of the 218 training columns stratified by event (seed: prespec section 4), reduced BO n_init = 20, n_iter = 40 over
   the spec section 3 box (alpha logit, lambda log10, nu linear, rho log10, ntop integer), feasibility = all five
   inner-fold fits converged, infeasible = NA + flag excluded from surrogate and selection, selection = highest mean
   inner-CV C among feasible at full precision, then the final n = 218 fit (30 starts, recording on) with the
   per-start and checkpoint test C as in item 2, plus the five inner-fold fits at the selected configuration kept for
   Q1. Use the SAME BO engine as the gate spec (DeSurv 1.0.1 `desurv_cv_bayesopt` internals: candidate pool 4,000,
   exploration weight 0.01, its GP and acquisition code). First find whether DeSurv exposes the BO core so that a
   custom objective (our inner CV) can be plugged in. If it does not, report exactly what would have to be copied and
   from which file:line, and implement the copy verbatim in `bo_engine_copy.R` with a test that it proposes the same
   points as the package for a fixed toy objective and seed. If no configuration is feasible: record the Q1 failure
   (prespec section 5). Saves under `out/tune/`.
4. **Q4 records.** Q4 needs the W-block Hessian at sweeps 1, 10, 100 and the end of the RETURNED start, which is only
   known at the end. Store (W, H, theta) at sweeps 1, 10, 100 for every start (add these to the recording without
   changing any mathematics), and compute Q4 post hoc in the evaluator by power iteration on `palm_hvp_W` restricted
   to free coordinates (prespec section 7: start vector ones on free coordinates, normalized; stop at relative change
   < 1e-8 or 500 iterations; pass estimate <= palm_LW (1 + 1e-8)). Report if adding storage requires touching anything
   beyond the recording code.
5. `qual_eval.R`: reads only `out/`; computes every binding criterion exactly as prespec section 7 and spec section 9
   (Q1, Q3 with conditions (i)-(v), scope, not-evaluable, numerical-failure and active-set-transition rules; Q4; Q5
   one-sided, cell = (scenario, configuration, size), mean over replicates; Q6; Q7 with t_{0.95,19}) plus the realized
   SD-ratio harness check (prespec section 3.4), Q2 and V1-V3 as reported tables. Writes `QUALIFICATION_RESULT.md`
   and `.rds` with one PASS/FAIL line per binding criterion and an overall line. It must refuse to run if any expected
   output file is missing (counts from the job manifest), never silently on a subset.
6. `qual_controls.R`: the verifier controls of prespec section 7: Q3 positive control on the failed candidate's four
   stalled T7 starts (1, 5, 9, 16; `sim_data(700, 3000, 218)`), using W at the window ends of the exact-gated D-Q4
   replay (`candidate/tests/DQ4_step_decomposition.R` and its inputs; find which recorded trajectory gives the window
   and say so), theta-equivalent s0 b_j / sigma_j(W); must flag all four. Q3 negative control on a converged start.
   Q7 controls on S3 replicate 1's test set (true A > 0.9, true B < -0.9, A + B |D_AB| < 0.1). Q6 control. This uses
   the candidate code read-only and existing candidate outputs; if a needed trajectory is not stored, report what
   rerun would produce it rather than running a p = 3,000 candidate job.
7. `manifest.R`: the full job list (grid: 3 x 5 x 12 = 180; tuning: 5 + 5 + 20 + 5 = 35) with output paths.
8. Slurm (Longleaf) scripts in `slurm/qual/`: array jobs for grid and tuning, following `slurm/run_rev3_bo_stage.sh`
   and `slurm/rev3_env.sh` conventions (modules, library path, logs in `slurm/logs/`). Tuning runs parallelize the 30
   starts of each inner-fold fit across cores; size requests from the timing file
   `analysis/estimator-gate/qualification/timing/timing_wall.txt` if it exists when you get there, otherwise leave
   placeholders marked TODO.
9. **Plumbing test (not qualification):** run grid, tuning (with n_init = 2, n_iter = 1, 3 starts) and the evaluator
   end to end at p = 300 on seeds 49,300,001+ (outside every qualification namespace), outputs under `out_plumbing/`.
   It shows the pipes connect and the evaluator counts correctly; its PASS/FAIL values mean nothing and must be
   labelled so. Do NOT run anything at p = 3,000 and do NOT write into `out/`.

Rules: base R (data.table allowed), style of the candidate and redesign code. Every script prints counts of what it
read and wrote. Final reply to the coordinator at most 150 words: what was built, plumbing result, every mismatch.
