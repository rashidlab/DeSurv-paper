# Handoff 2026-10-05: candidate optimizer built; frozen algorithm does not converge at realistic size

Branch recompute/consensus-desurv, all unpushed: 26b9ee7, bb0d967, e1502f4 (spec amendment 2 and rulings; window
CLOSED), b0d6b62 (candidate optimizer + tests), ac8d805 (commit-review fixes).

## State

- Spec: `docs/plans/2026-10-04-estimator-gate-spec.md` (amendment 2, window closed by the PI 2026-10-04).
- Plan: `docs/plans/2026-10-04-candidate-optimizer-plan.md`; codex review and prompt beside it.
- Code: `analysis/estimator-gate/candidate/` (base R reference; `REPORT.md` has all counts and the coordinator
  verification). Simulated data only; nothing reads real data.
- Preconditions: T1-T6, T9, T11, T12 PASS (verified by rerun); T4 0 F increases in 1.2M block steps (recount
  confirmed); T10 24 of 28 mutants caught by the named test (4 misses are T4 rows whose trigger cannot occur on
  correct code; T9 catches them).

## The blocking finding

At p = 3000 (n = 175 and 218, logit-midpoint hyperparameters) 0 of 60 starts reach r < 1e-5 in 5,000 iterations.
Continuing four starts to 20,000 iterations: r_W oscillates at 2.4e-4 to 4.8e-4 with no downward trend while F
still creeps down (accepted t 0.008 to 0.125). Reading: ill-conditioned projected-gradient W step, not a rounding
floor. At p = 150 the same point converges (417 iterations); large alpha with weak lambda does not converge even
at p = 150 (r 0.07 to 0.18). Under the frozen spec this predicts Gate A ineligibility (decision table: specification
failure). Base-R gate cost estimate 85,268 CPU-hours.

## Decision needed from the PI (not taken)

The amendment window is closed. Options: (a) run as frozen and record the expected ineligibility (large compute
for a predictable outcome); (b) declare the frozen candidate a specification failure now, from simulated-data
evidence, per the decision table's last row; (c) reopen the window for an optimizer-only change that keeps F and
the r < 1e-5 certificate fixed (e.g. a preconditioned or second-order W block, outer extrapolation, larger maxit).
No real data or outcome has been seen, so (c) is not outcome-driven, but it reverses the PI's closure.

## Next

PI decision on (a)/(b)/(c). Optional diagnostic before deciding: Hessian-vector Lanczos at a T7 end point to
measure the W-block condition number. Do not push without the PI. Do not start the C++ port before the decision.

## PI ruling 2026-10-05 (supersedes "Decision needed" and "Next" above)

Option (c), gated on the diagnostic, as a one-time numerical-feasibility amendment. Not (a): 85k CPU-hours to
confirm what T7 shows. Not (b): nonconvergence of one solver is not evidence the objective lacks a usable
stationary solution. A better optimizer can land in different basins of a nonconvex problem, so this is not "the
same fit, faster": it needs a documented amendment and full simulated requalification.

1. **Diagnostic first** (simulated data): W-block condition number at a T7 end point (Hessian-vector Lanczos);
   check whether the slow direction aligns with the stalled KKT residual r_W. Proceed only if it supports the
   conditioning hypothesis.
2. **Amendment 3, written before implementing.** Unchanged and frozen: objective, initialization, 30-start seed
   schedule, KKT residuals, r < 1e-5 certificate, start-selection rule, BO design and tuning space, regularizer,
   Gate A/B rules, biology endpoints. Changed: only the W-block solver. State that the amendment exists because
   prespecified simulated qualification (T7) showed computational nonconvergence at p = 3000 before any real-data
   gate run. Close the window again afterwards.
3. **Solver choice:** a principled preconditioned or second-order / proximal block method with monotone line
   search and the same feasibility constraints (S_+^k projection). Not "increase maxit". Extrapolation only behind
   a preconditioned method, not as the fix.
4. **Requalification:** rerun T1-T12 (gradients, zero material increases, same KKT oracle, recovery, mutation
   tests) plus: (i) realistic-size convergence at p = 3000 with a criterion fixed in amendment 3 BEFORE running,
   e.g. >= 90% of realistic-size final fits reach r < 1e-5 within a per-fit budget set from what the gate can
   genuinely afford (not from what the new solver achieves); (ii) solver-equivalence test on the p = 150 settings
   where the current solver converges: compare final F, KKT residuals, loadings and predictions, and selected
   start, to see whether the new solver reaches the same stationary points or different basins.
5. **Chronology wording** (for the spec and any write-up): "The statistical specification was frozen before
   implementation. Prespecified simulated-data qualification subsequently showed that the reference W-block
   optimizer failed to reach the frozen KKT tolerance at realistic dimensionality despite monotone objective
   decrease. Before any real-data gate analysis, the specification was amended once to replace only the W-block
   numerical solver; the objective, convergence certificate, tuning procedure, and gate rules were unchanged."

Open for amendment 3: the affordable per-fit compute budget (PI to set from real capacity).

## Step 1 result 2026-10-05: conditioning hypothesis SUPPORTED

`analysis/estimator-gate/candidate/tests/out/coordinator/wcond.R` (decision rule S1-S5 written in its header before
the run), output `wcond.txt`. Simulated data only. Operator: Riemannian Hessian of Psi in W on the active face of
S_+^k (H, b fixed), Hessian-vector products by central differences (h 1e-5 vs 1e-4 agree to 2e-6), Lanczos 300 steps
with full reorthogonalization. Four T7 end points (n = 218, starts 1, 5, 9, 16) and a converged p = 150 control
(dataset 401, n = 60, same hyperparameters; all 5 control starts converge in 270 to 556 iterations).

| | lambda_max | lambda_min | kappa | R_W weight on theta <= 1e-2 lambda_max |
|---|---|---|---|---|
| p = 3000, four starts | 47 to 75 | 9.5e-5 to 2.0e-4 | 2.4e5 to 5.0e5 | 0.92 to 0.99 |
| p = 150 control, two starts | 4.0 to 5.7 | 0.022 to 0.030 | 183 to 190 | 0.96 to 0.99 |

S1-S5 all TRUE at all four end points:
- S1: 3,000 W-only steps at fixed (H, b) leave r_W at 2.5e-4 to 5.4e-4 (no trend). The W block alone reproduces
  the stall, so the cause is not W-H-b coupling.
- S4: the boundary (W = 0, negative gradient) share of ||R_W||^2 is 0.000, so the cause is not active-set
  identification. Support is essentially full (2,982 to 3,000 of 3,000 per column).
- S5: no negative curvature, so the point is not a saddle.

Mechanism (per-term Rayleigh quotients at Ritz vectors): the stiff directions are entirely survival curvature
(3 outliers, 13 to 75). The bulk is reconstruction curvature at about 2e-4, roughly isotropic, which is
(1 - alpha) H H' / ||X||^2 on the face. The backtracking step settles at t·lambda_max of 1.5 to 2.9, at the edge of
stability on the stiff modes. That explains why t cycles between 0.008 and 0.125 and r_W oscillates. The bulk then
contracts by 1 - t·lambda_min, about 1 - 6e-6 per iteration, which is 1e5 to 1e6 iterations per e-fold. The
residual sits almost entirely in that slow bulk.

Caveat: the control differs in both p and n (150 x 60 against 3000 x 218), so the 1,300x to 2,600x kappa ratio is
not a pure p effect. Lanczos Ritz extremes bound the true kappa from below.

Solver implication for amendment 3 (a proposal, not decided): the spectrum is a near-isotropic bulk plus a few
stiff outliers. That is the case where Newton-CG / truncated Newton on the active face converges in roughly
(outliers + a few) CG steps per outer iteration. A diagonal or k x k (H H') preconditioner alone would not separate
the outliers from the bulk.

Next: write amendment 3 BEFORE implementing (scope per the PI ruling). Open: the per-fit compute budget for the
realistic-size convergence criterion (PI to set from real capacity). Do not push.
