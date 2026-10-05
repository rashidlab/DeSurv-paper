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
