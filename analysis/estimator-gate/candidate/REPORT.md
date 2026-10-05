# Candidate optimizer: implementation preconditions report (simulated data only)

Spec: `docs/plans/2026-10-04-estimator-gate-spec.md` section 4 at e1502f4 (amendments 1 and 2 and their rulings).
Plan: `docs/plans/2026-10-04-candidate-optimizer-plan.md` at bb0d967. Written by the build agent (Opus) and saved
by the coordinator, because the harness blocked the agent from writing this file. The coordinator's verification
is in the last section.

Environment: R 4.6.0, reference BLAS, DeSurv 1.0.1, glmnet 5.0, survival 3.8.6. The machine was shared (load 20 to
34 on 20 cores), so times are upper-biased. Nothing reads `results/`, `data/` or any cohort file.

## Verdicts and counts

| Test | Checks | Passed | Worst value | Verdict |
|---|---|---|---|---|
| T1 gradients | 1,554 (part A 810 = 54 corners x 5 points x 3 blocks; part B 744 per-term, 420 coordinate directions) | 1,554 (1 via Richardson); 0 fail; 0 inconclusive | A: rel 3.6e-8; B: one direction 2.1e-6 plain, 8.2e-7 Richardson | PASS |
| T2 objective vs coxph | 810 (F, loglik, beta at 270 points) | 810 | 4.4e-16 | PASS |
| T3 b block | 360 (120 independent KKT, 120 enumeration solver, 90 x 1,000 perturbations, 30 glmnet) | 360 | KKT 9.99e-11; coef 2.6e-9 vs enumeration; 0 of 90,000 perturbations lower; glmnet coef within 6.2e-7, never lower in objective by > 1e-12 | PASS |
| T4 monotone F | 135 runs, 1,204,041 block steps | 0 increases (counter and independent recount); 0 failures; residuals recorded; W step <= 2; b never raised F | largest raw increase 6.9e-18 | PASS |
| T5 projection and draws | 9 | 9 | W0 bitwise identical to DeSurv 1.0.1 for 30 of 30 starts | PASS |
| T6 HALS | 36 (90,000 sweeps; 18 exact-solution matches) | 36 | 0 increases; 2.2e-15 from exact | PASS |
| T7 timing | 2 fits x 30 starts | eps shared in both; 0 failures | 0 of 60 starts converged | PASS (recorded) |
| T8 recovery | descriptive | n/a | cosines with true W 0.38, 0.34, 0.36 | n/a |
| T9 edge cases | 54 + 1 observation | 54 | see below | PASS |
| T10 mutation | 10 controls, 28 mutant runs | 10 of 10 controls; 24 of 28 mutants caught | 4 misses, all T4 | FAIL as defined |
| T11 fold-local | 7 | 7 | fit identical after changing test outcomes and expression; simulated C = 0.858 | PASS |
| T12 transport | 7 | 7 | score identical after removing or reordering genes outside G | PASS |

**T4 (p = 150, n = 60; 3 datasets x 9 corners x 5 starts; 1,200 s budget per start).** 67 converged (49.6%), 67
nonconverged, 0 failed, 1 stopped by the time budget. Corners: (alpha 1e-4, rho 1e-2, lambda 1, nu 0.95),
(1e-4, 1e3, 1, 0.95) and (0.9, 1e3, 1, 0.95) converged 15 of 15; (1e-4, 1e3, 1e-3, 0) 14 of 15; the box centre
(alpha 0.3 linear, rho 1, lambda 1, nu 0.5; median r 6e-4), (0.9, 1e3, 1e-3, 0) (4.6e-4) and (0.9, 1e-2, 1e-3, 0)
(0.18) 0 of 15. Seven rejected b steps (one re-run: rounding, F up 3.2e-16 relative). Largest b block 29,201 FISTA
iterations (cap 1e5).

**T7 (p = 3000, k = 3, logit-scale box midpoint alpha = 0.0291, lambda 1, nu 0.475, rho 3.16; 15 cores).**

| n | Converged | Seconds per start mean (min / median / max) | s per iteration | median r at stop | parallel efficiency |
|---|---|---|---|---|---|
| 175 | 0 of 30 | 568 (507 / 565 / 646) | 0.113 | 3.0e-4 (r_W) | 0.98 |
| 218 | 0 of 30 | 635 (572 / 642 / 697) | 0.128 | 3.7e-4 (r_W) | 0.97 |

Every start ran 5,000 iterations (about 5,000 W backtracks, 120,000 FISTA iterations). One 30-start fit: 4.7
CPU-hours (n = 175), 5.3 (n = 218). Gate extrapolation: 540,000 x 568 s + 600 x 635 s = 85,268 CPU-hours in base R.

**T10 misses**, all T4 runs whose trigger never occurs on correct code: B7 (floor never reached; T9 catches it), B9
(no increase to miss; with an F-raising H block, T4's recount catches it), B18 and B19-b (no increase occurs; T9
catches both). All other rows fail their named test as intended (`tests/out/T10_mutation.txt`).

**T9 covers** uniform, near-uniform and duplicate W columns, b = 0, sparse W and zero H rows; Breslow for every tie
pattern, a one-event dataset, large offsets and contrasts (all within 1e-11 abs, mostly 1e-14); 11 invalid-input
cases failing deterministically; sample-permutation invariance and seed replay (serial and 3 cores); r_W = 0 at a
stationary point where r_original = 2.0; the W cap; rejected b steps leaving status unchanged; each status
transition (FISTA cap, step floor, F increase, maxit, converged); the fit status rule.

**Grep** for `zeta`, `results/`, `concordance`, `cidx`: one hit, `README.md:5` (the sentence saying nothing reads
`results/`). The only C-index code is `harrell()` in `tests/T11_fold_local.R`, on simulated test samples.

## Conflicts and problems

1. **Convergence (most important).** r_W stalls near 3e-4 at p = 3000; the box centre does not converge at p = 150
   either. H and b residuals are far below tolerance. Under the frozen spec this makes BO configurations
   infeasible and fails Gate A's 90% precondition. No change made to optimizer or criterion.
2. **W-step floor effectively unreachable (spec gap).** With a deliberately uphill W gradient, once t reaches
   rounding scale (1.8e-15 to 3.5e-18) W+ differs from W only by rounding, the sufficient-decrease test is decided
   by noise, steps are accepted, and the 1e-20 floor is never reached. A broken W block shows as nonconverged, not
   as a numerical failure. Implemented as written; exactly-null steps counted in `w_null_steps` (I14).
3. **Agent's own FISTA restart defect, fixed.** The first build restarted momentum on any objective rise; near the
   solution that fired on rounding and 15 starts hit the 1e5 cap. Newton reached KKT 6e-14 on the same
   subproblem (no precision floor). A gradient-based restart (I16) certifies it in 1,837 iterations. Superseded
   runs in `tests/out/superseded/`.
4. **Uniform columns on rank data.** A uniform W column gives score sd about 4e-14 from rounding, not 0, so the
   "pooled median sigma <= 0" check does not fire. Random starts never produce one; observation only.
5. **Re-ranking in a reduced gene universe (T12, descriptive).** Dropping 25% to 75% of genes outside G before
   ranking shifts external scores by 2.7 to 7.7 training-score sd, Pearson >= 0.993. Within-cohort C unaffected;
   anything on an absolute scale is.
6. **Poor recovery in T8 (descriptive).** Fitted W columns correlate 0.78 to 0.81 with each other.
7. **Seed streams.** BO evaluations in DeSurv 1.0.1 use the warmstart stream, matched bitwise; the comparator's
   final fit uses `init()`. The seed namespaces keep these separate; not a conflict.
8. **Extra input validation.** d must be 0/1, gene names unique. A nonnegativity check on X was removed (not in the
   spec).
9. No formula in spec section 4 disagreed with the agent's derivation.

## Interpretations added beyond I1-I9 (details in README.md)

I10 FISTA step test on gradient differences; I11 H = 0 before the 20 HALS sweeps; I12 b block warm-starts L0 =
previous L / 4, solver KKT proposes the stop and the independent evaluator certifies it; I13 finite-difference
method in T1; I14 accepted null W steps counted, not failures; I15 per-start time budget in the T4 harness only;
I16 FISTA momentum restart on the gradient-mapping test.

## Coordinator verification (2026-10-05)

- T1 and T2 rerun independently: identical totals (1,554 of 1,554; 810 of 810), VERDICT PASS.
- T4 output: counter 0 and independent recount 0, largest raw block change 6.9e-18. Confirmed.
- Correction: T4's "box midpoint" is alpha = 0.3 on the linear scale, not the logit-scale midpoint T7 uses. At
  p = 150 (dataset 401, start 1) the logit midpoint (alpha 0.0291) converges in 417 iterations.
- T4 medians by corner: nonconvergence concentrates at alpha 0.3 and 0.9 with weak lambda (r 0.07 to 0.18 at
  alpha 0.9, rho 1e-2: not near convergence); small alpha converges in 120 to 1,370 iterations.
- T7 continuation (`tests/out/coordinator/cont.R`, output `cont.txt`): four n = 218 starts (1, 5, 9, 16) continued
  from their iteration-5,000 state for 15,000 more iterations (20,000 total). r_W does not trend down: it
  oscillates between 2.4e-4 and 4.8e-4 in every start (at 15,000: 3.3e-4, 3.4e-4, 2.4e-4, 2.4e-4). F is still
  decreasing, by about 1e-5 to 6e-5 over the 15,000 iterations, and accepted W steps are small (t 0.008 to 0.125,
  so local curvature about 8 to 128 against the cap of 2). r_H 1e-6 to 1e-5 and r_b about 1e-10 are not binding.
  Reading: a slow, ill-conditioned crawl of the projected-gradient W step (consistent with very different
  curvature between the rank-scale survival term and the 1/||X||^2-scaled reconstruction term), not a rounding
  floor (F still moves at about 1e-9 per iteration). The frozen algorithm does not reach r < 1e-5 at this point
  in 20,000 iterations, four times maxit, and shows no trend toward it.
- Consequence under the frozen spec: if real folds behave like this simulation at moderate alpha, most inner-BO
  configurations are infeasible and the outer final fits nonconverged, so Gate A would be ineligible
  (decision table: specification failure). This is a simulated-data finding; no real data or outcome was used.
- Commit-review fixes (codex review of b0d6b62), applied by the coordinator: (1) the fast Breslow path returned
  NaN when censored risk sets underflowed (y = 1:3, d = (1,0,0), eta = (0,-1000,-1000)); it now falls through to
  the log-sum-exp path unless ll and g are finite. (2) The KKT oracle used one global exponential shift; it now
  uses per-risk-set log-sum-exp. (3) A start was marked a numerical failure when the non-binding r_original was
  nonfinite (H = 0); finiteness is now checked on r, r_W, r_H, r_b only. After the fixes T1, T2, T3, T5, T6, T9,
  T11 and T12 rerun with identical totals and VERDICT PASS. T4 and T7 were run on the pre-fix code; the fixes touch
  only branches those runs never entered (0 numerical failures, no H = 0 state), so they stand.
