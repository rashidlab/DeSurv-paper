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
- S5: no substantial negative curvature on the W face (smallest Ritz value positive at all four points). This does
  not exclude negative curvature in the H or b blocks or in coupled directions, which the diagnostic does not probe.

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

## PI ruling 2026-10-05, second (amendment 3 scope and budget; supersedes "Open" above)

Proceed with amendment 3. Basis: the W block alone reproduces the stall, and 92-99% of the residual sits in the
lowest-curvature directions. That defines a specific numerical problem; it is not licence to redesign the estimator.

- **Budget (algorithmic, not wall clock):** at most 1,000 Newton-CG outer iterations per fit; at most 100
  Hessian-vector products per W update; at most 25,000 W-block HVPs per fit. Convergence must occur under the frozen
  KKT/residual criterion before a cap; hitting a cap is failure, not convergence. Record wall time, outer
  iterations, HVP count, objective path and final KKT residual for every benchmark fit (wall clock is not a pass/fail
  criterion: it depends on hardware, BLAS and parallelization).
- **Realistic-size qualification:** on the frozen p = 3000, n = 218 benchmark end points, the amended optimizer
  converges all four previously stalled fits (starts 1, 5, 9, 16) within the budget, without worsening the objective
  or changing the fitted solution beyond a prespecified numerical-equivalence tolerance.
- **Non-regression control:** p = 150. The new solver reaches the same solution and does not materially degrade a
  case the current algorithm already solves.
- **Amendment scope:** (1) problem: pathological W-block conditioning at realistic p, established prospectively by
  S1-S5; (2) permitted change: the W update only, replaced by active-face Newton-CG; (3) frozen: objective,
  penalties, H update, b update, normalization, initialization, stopping criterion, estimator definition;
  (4) acceptance: numerical equivalence on small/converged controls plus convergence of all four stalled end points
  within the budget; (5) failure rule: if Newton-CG does not pass, do not tune it against T7 outcomes; diagnose and
  write another amendment before changing the algorithm.
- **No preconditioner in amendment 3.** Test plain active-face Newton-CG first. A preconditioner needs evidence
  that CG itself is now the bottleneck. Do not prespecify "a handful of CG steps" (a hypothesis from the spectrum,
  not established).

### Items the amendment-3 draft must settle (raised by the coordinator, for the PI)

1. "Per fit": per start or per 30-start fit? The budget figures read as per start (palm_fit_one); state it.
2. Cap = failure changes status semantics. The frozen spec has "nonconverged" (maxit reached, finite, no failure)
   as a separate status, and the frozen start-selection rule falls back to the lowest-F nonconverged start. Does
   cap = failure apply only to the qualification benchmark, or does it replace "nonconverged" in the estimator? The
   latter touches the frozen start-selection rule.
3. Equivalence at p = 3000 is "relative to what"? The stalled end points are not stationary, so the solution is
   expected to move off them. Proposal: at p = 3000 require F_new <= F_endpoint and report loading correlation and
   prediction agreement against the end point descriptively; apply the numerical-equivalence tolerance (F, KKT,
   loadings, predictions, selected start) at p = 150 only, where both solvers reach r < 1e-5.
4. Benchmark start: warm start from the stalled end points (tests the local problem, which is what S1-S5
   diagnosed) or from the frozen initialization (tests the full path)? The ruling says "end points"; the
   earlier ruling's item 4(i) said "realistic-size final fits". Proposal: end points binding, full 30-start fits
   from initialization recorded as secondary.
5. Specification details the draft fixes: what counts as one outer iteration (one H-W-b sweep with the W block a
   Newton-CG step); the HVP source (analytic, or central differences of the gradient as in wcond.R); handling of
   negative curvature in CG (Steihaug truncation, so the step stays a descent direction); the monotone line search
   along the projected Newton path with the same S_+^k projection; the CG relative-residual stopping rule.

## PI ruling 2026-10-05, third (settles items 1-5 above; frozen for the amendment-3 draft)

Proposals for items 2-4 adopted; item 3 refined; item 5 constrained on the HVP source.

1. **Budget is per start.** A 30-start fit has up to 30 independent per-start budgets (1,000 outer iterations,
   25,000 W-block HVPs, 100 HVPs per W solve, each per start). 25,000 HVPs is never a fit-wide pool. The 25,000
   total is the binding global work cap per start; 100 is only the ceiling for one W solve. Do not relax 25,000
   on the ground that 1,000 x 100 is larger.
2. **Cap = failure applies to amendment-3 qualification only, not to the estimator's frozen status logic.** Do not
   change the existing "nonconverged" status or the start-selection fallback. In production, a start that hits
   the Newton-CG work cap is nonconverged, with a diagnostic reason recorded (e.g. `work_cap_reached`); the frozen
   rule may still select the lowest-F nonconverged start if no start converges. In the amendment-3 benchmark, a
   start that hits the cap means the new solver fails qualification. This keeps the amendment W-solver-only.
3. **p = 3000: no numerical equivalence to the stalled end point** (that end point is what the solver must escape).
   Binding requirements, per warm-started end point: (a) satisfy the already-frozen convergence/KKT criterion;
   (b) F_final <= F_stalled_endpoint; (c) stay within the per-start work budget. No prespecified minimum
   improvement in F (improvement is expected, not required). Loading and prediction agreement with the stalled end
   point are reported descriptively only, never as acceptance criteria. Numerical equivalence (F, KKT, loadings,
   predictions, selected start) applies at p = 150 only, where both algorithms reach stationary solutions and so
   share a meaningful target.
4. **Warm-started four end points (starts 1, 5, 9, 16) are the binding realistic-size qualification.** They test
   the pathology S1-S5 established without mixing in initialization or multi-start selection. The full frozen
   30-start p = 3000 fits are a secondary end-to-end check, not the amendment-3 pass/fail gate. Before the amended
   solver is used for the definitive recomputation, those full fits must still show the unchanged multi-start
   pipeline behaves sensibly under the new W solver. Two distinct questions:
   - amendment-3 qualification: did Newton-CG fix the diagnosed W-block stall?
   - end-to-end qualification: does the unchanged estimator pipeline work with the corrected W solver?
5. **Mechanics are the draft's to specify, with one constraint: analytic Hessian-vector products, not finite
   differences**, unless a compelling implementation obstacle arises (if one does, report it to the PI; do not
   substitute differences silently). Rationale: the point is to fix conditioning cleanly; finite-difference noise
   in the curvature operator makes both convergence and the HVP budget harder to interpret. The draft still fixes
   the outer-iteration definition, Steihaug negative-curvature truncation, the monotone projected line search, and
   the CG relative-residual stopping rule.

**Frozen summary:** per-start budgets; cap failure only for qualification; p = 150 equivalence, p = 3000
convergence plus F_final <= F_endpoint (not endpoint equivalence); four warm end points binding, full 30-start
fits secondary/end-to-end.

**Next:** write amendment 3 text from the second and third rulings; implement analytic W-block HVP.
