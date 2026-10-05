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

## PI ruling 2026-10-05, fourth: amendment 3 APPROVED and frozen; window closed

Amendment 3 approved with items 1, 2, 4, 5a, 6 as proposed (with refinements) and items 3, 5b, 7, 8 changed
(retracted path instead of step cap; invariant p = 150 criteria; scale-aware curvature test; W no-op once
r_W < 1e-5). The p = 3000 criterion is certificate, proportion, monotonicity and budget; lower F is evidence, not
the criterion. Full record: `docs/plans/2026-10-05-amendment3-newton-cg-w-solver.md` section 9. Window closed;
bug fixes only.

**Next:** implementation brief for the analytic W-block HVP and the Newton-CG W step (Q0 first), against the
frozen text.

**Closed 2026-10-05 (PI).** T4/T9/T10 remapping confirmed as a consistency correction; Armijo acceptance requires
<G_W, Delta W(t)> < 0; B19 replacement mutant is B19_retraction. Amendment 3 closed. Next session implements
against the frozen document, analytic HVP and Q0 first; no further solver-design changes.

## Stage A DONE 2026-10-05: analytic HVP built; amended Q0 PASS; Stage B unblocked

- Brief: `docs/plans/2026-10-05-amendment3-implementation-brief.md` (PI rulings on A5/B3 and the Q0 result
  appended; all six coordinator choices approved; added Stage B diagnostic: direct ||Hess d + g|| vs recursive
  CG residual at CG exit, every Q1 solve and a seeded Q-T sample, not a gate).
- 045f593: HVP code in `palm_lib.R` (72 lines added, 0 deleted). Original Q0 FAILED on the order clause only
  (Euclidean 1120/2700, Riemannian 2624/2700; accuracy held on all 5,400). Retained.
- 68dc7e9: test-only Q0 amendment, amendment 3 section 10 (half-decade grid 1e-1..1e-8, order by eps ratio,
  two consecutive p_k >= 1.5), committed before running. 2c9e29c: amended Q0 PASS 2700/2700 both, recomputed
  independently by the coordinator. Report: `analysis/estimator-gate/candidate/AMENDMENT3_REPORT.md`.

**Next (new session):** dispatch Stage B of the brief (Newton-CG W step, breakage inventory B2, Q-T, Q1, Q2, Q3 report)
to a subagent; coordinator verifies Q-T verdicts, Q1 4/4, Q2 criteria, and the residual-drift diagnostic before
committing. Do not touch the HVP. All commits on recompute/consensus-desurv unpushed; do not push without the PI.

## Stage B DONE 2026-10-05: amendment 3 FAILS; diagnostics localize the stall to a W-b zig-zag

- b0330ab: PI interpretation ruling recorded before any Stage B result (drift tier (c) binding, other drift raw only).
- 0c44cca: Stage B failure record. Q1 1/4 (starts 1, 5, 16 hit the 1,000-sweep cap at r_W 1.5e-3 to 3e-3, using
  3.1k-3.7k of 25k HVPs); T9 55/56 (permutation W diff 7.3e-3 vs 1e-8); T7, T8, T10 not rerun; Q2 not run
  (settings undefined, candidates in the report); drift tier (c) = 0.
- 08ef292 + 86f0e88: PI rulings on the outcome (T7 mandatory; T10 untriggerable mutants to be replaced with
  branch-forcing tests; Q1 budget strict <; Q2 moot; T9 not relaxed; no amendment 4 before diagnostics).
- fad6dd3: diagnostics approved, settings frozen before running. 09801d6 + 450e3f0: results.
  - D-Q1 (exact-replay gate PASS 3/3): r_W after H -> W -> b = 3.1e-3 -> 1.3e-4 -> 3.1e-3 (start 1), 2.3e-3 ->
    8.7e-5 -> 2.3e-3 (5), 1.5e-3 -> 4.9e-5 -> 1.5e-3 (16). The W step pushes r_b from ~1e-10 to 1e-4..5e-4 and the
    b step undoes W's stationarity. H is inert. |dW| ~1e-8, |db| ~2e-6 per sweep; sigma and beta medians do not move
    at printed precision. Outlier sweeps: r_W after W reaches 0.32 (start 1) and 0.16 (start 16) in a few of the
    last 100 sweeps (not yet examined).
  - D-T9a (gate PASS): the truncated pair is the same fit on every invariant (cosines 1.0000, eta r 1.000000,
    top-5 overlap 1.00, recovery within 1e-4). D-T9b: both converge (188/190 sweeps), W diff 7.2e-5, F diff 4.5e-11:
    pre-stated interpretation (i). The raw-W T9 criterion is the wrong invariant for an NCG path; replacing it needs a
    PI ruling.

**Next (new session):** PI reads the D-Q1/D-T9 results, then amendment 4 design targeting W-b coordination (e.g. a
joint W-b step or safeguarded acceleration of the full sweep map), plus the prospective test-contract fixes (T10
mutants, T9 invariant criterion, T7 in scope). Do not touch the HVP or Newton-CG. All commits on
recompute/consensus-desurv are unpushed; do not push without the PI.

## PI ruling 2026-10-05, Amendment 4 boundary (Option A)

Option A (variable-projection W step) is allowed only through Amendment 4. Amendment 3's `H_WW[V]` is frozen and
stays byte-for-byte unchanged. Amendment 4 may build a separately derived reduced-objective operator around it,
H_red[V] = H_WW[V] - H_Wb H_bb^{-1} H_bW[V], and explicitly supersedes "do not touch the HVP" for this wrapper only.
Three protections are required before it becomes the solver:
1. Qualify H_red independently, like Q0: compare it to a finite difference of the reduced gradient
   grad_W F{W, b*(W)}, not the fixed-b gradient. The existing H_WW Q0 stays valid and is not rerun.
2. b has an L1 term, so the Schur formula holds only on a locally fixed b active/sign set. Record that set at the
   stalled end points and along the proposed perturbations. If it flips, stop and rethink the derivative.
3. The line search evaluates the profiled objective phi(W) = F{W, b*(W)}, re-solving b with the frozen solver at
   every trial W(t).
Inspect the r_W = 0.32 / 0.16 outlier sweeps and b active-set behaviour before freezing any D-Q2 decision threshold.

## Outlier and active-set inspection 2026-10-05 (descriptive; nothing frozen, nothing run beyond reading saved fits)

Full record: docs/plans/2026-10-05-dq1-outliers-and-score-collapse.md. In short: the b sign/active set never changes
(protection 2 is met at the end points). The outlier sweeps are not a separate mechanism; they are the long-CG steps
that carry the end points further along a path where the score scale sigma of the coefficient-carrying factors
shrinks with no plateau, and every end point has collinear score columns. This may be the statistical-design
evidence the PI said would be needed to reopen the estimator. D-Q2 threshold NOT frozen; PI decision required.

## D-Q3 DONE 2026-10-05 (prespec 50a79be; results 7b15de6). Amendment 4 ON HOLD for the PI

- 70bb55e + 9964e70: D-Q1 follow-up record. Every end point has collinear score columns; the coefficient-carrying
  sigma shrinks with no plateau.
- d01d2b0: D-Q3b. The one-dimensional shrink-to-uniform path is NOT a descent direction. The supervised terms are
  exactly invariant along it, and reconstruction pays.
- 7b15de6: D-Q3a.
  - All three continuations ended on `line_search_failed`, at 1,581-4,604 sweeps.
  - sigma kept falling (slopes -0.31 to -0.50), scores stayed collinear, F kept decreasing.
  - Post-hoc: curvature along the residual grew from about 1 to between 1e6 and 4e6, almost all from the rho ridge
    term (about sigma^-4). From these end points, r_W < 1e-5 cannot be reached in double precision.
- Full record: `docs/plans/2026-10-05-dq1-outliers-and-score-collapse.md`, section 6.

**Open for the PI:**
1. Category for D-Q3a: "collapse continues", or mixed, given the unanticipated early numerical termination.
2. Whether to pre-specify a step-decomposition diagnostic. It would split the long-CG W steps that carry the sigma
   drop into a shrink part and a rotation part relative to the score direction, using the saved runs only. The aim
   is to name the descent direction before the estimator is redesigned; the likely targets are the sigma
   standardization and the rho term's 1/sigma.

**Standing rulings:**
- No real-data diagnostics until the estimator mechanics are settled.
- No Amendment 4 draft. Option A, if ever pursued, is a separately qualified Schur wrapper around the frozen H_WW.
- Nothing pushed.

## PI ruling 2026-10-05, fifth and sixth (settle the two open items above)

1. D-Q3a = "collapse continues", with the numerical-endpoint caveat (wording in the D-Q1 record, section 7). Option A
   stops. The burden of proof is on any fourth solver amendment: it would have to explain why a better optimizer
   helps on a trajectory whose curvature is becoming singular.
2. Step-decomposition diagnostic (D-Q4) approved: pre-specify, then an instrumentation-only deterministic replay of
   the D-Q3a runs (starts 1, 5, 16), gated on exact reproduction of the saved trajectory, decomposition computed
   online at every accepted W step, summaries only saved. Primary decomposition in centered score space; PW space
   secondary. Every accepted step, stratified afterwards by recorded CG count. Descriptive, non-gating.
3. No Amendment 4 draft. 4. No real data.

**D-Q4 draft:** `docs/plans/2026-10-05-dq4-step-decomposition-prespec.md`, NOT frozen. Open for the PI:
- Q-A. In score space the orthogonal part can only raise ||z|| (Pythagoras), so all sigma loss is assigned to shrink
  by construction and "mostly rotation" is unobservable. Proposal: attribute sigma loss with the exact loading
  identity sigma~ = ||P w|| g(v), keep score space as the description of direction change.
- Q-B. Cutoffs S >= 0.75 / S <= 0.25 (coordinator proposal).
- Coverage note (section 6 of the draft): the D-Q3a window holds only about 5-22% of the log-sigma drop; most of it
  is in Q1 sweeps 1-1,000. Optional addition, same gate pattern.

## PI ruling 2026-10-05, seventh: D-Q4 FROZEN (settles Q-A, Q-B, coverage above)

- Q-A: the primary attribution is the exact loading identity, Delta log sigma = Delta log ||Pw|| + Delta log g(v).
  Standardized-score angle is a separate direction-change diagnostic. No score-angle cutoff for the null-space case:
  report the three pieces separately, and describe g-driven loss with a near-zero score angle as orientation into a
  low-variance or null-space direction with the score direction preserved, not as rotation.
- Q-B: no cutoffs. Continuous quantities only; "predominantly norm-driven" or "predominantly g(v)-driven" only when
  one part is clearly larger in absolute magnitude.
- Coverage: Q1 sweeps 1-1,000 added, gated exactly against `Q1_warm.rds` (plus D-Q1 per-sweep sigma). Q1 and D-Q3a
  analyzed separately and combined. Solver settings unchanged.

**NEXT (fresh session):** write `tests/DQ4_step_decomposition.R` per the frozen prespec, smoke-test at 30 sweeps per
segment (discard outputs), run, commit results as observed. Still no real data, no Amendment 4, nothing pushed.

**D-Q4 DONE (2026-10-05).** `tests/DQ4_step_decomposition.R`; all six start-segment gates exact; results in
`tests/out/DQ4_step_decomposition.{rds,txt}` and section 8 of `docs/plans/2026-10-05-dq1-outliers-and-score-collapse.md`.
Reading: sigma loss is predominantly g(v)-driven in every start and in both segments, with the standardized-score
direction preserved (net angle <= 0.0023 rad). This is orientation into low-variance or null-space directions of X,
a rescaling of the same score, not loading shrinkage. Nearly all of it comes from the 1-5% of steps with >= 10 CG
iterations. Not pushed. NEXT: PI reads section 8; estimator redesign scoping is a PI call (no Amendment 4 yet).

## PI ruling 2026-10-05, eighth (after D-Q4): reopen the estimator

1. No Amendment 4 solver fix. Stop Newton-CG, Schur-complement and outer-acceleration work; a better optimizer moves
   faster along the pathological direction. Amendment 3 stands as a legitimate failure record.
2. Primary target: the sigma(W)-standardized parameterization. With z = C X' w, sigma proportional to ||z|| and
   eta = z b / sigma, moving w toward small-g(v) directions shrinks z and sigma together, so z / sigma (the survival
   predictor) barely changes; the rho term carries the same 1/sigma and produces the 1e6-level curvature. A
   parameterization problem, not a numerical one.
3. Preferred redesign: eta = Z_c beta with beta optimized and penalized directly (no beta = b / sigma(W) inside the
   objective). Keep the column normalization of W. The rho penalty is written in terms of P W beta. No 1/sigma
   singularity. sigma may still be used after fitting to report standardized Cox effects. No variance floor on sigma
   or g(v) as primary; a denominator-floor formulation is a fallback candidate only, not co-primary.
4. Treat this as a new estimator specification, not an optimizer amendment. Freeze before implementation: k = 3;
   W/H normalization; objective in the beta parameterization; penalties and tuning ranges; form of the rho term;
   initialization; convergence/KKT definition; nested tuning; Gate A (-0.02) unchanged; Gate B hidden until Gate A
   is committed. The current implemented estimator remains the comparator. No real data.
5. Qualification on simulated data before Gate A must test the pathology directly: ordinary convergence without
   sigma -> 0; no preserved-score/null-space collapse; bounded curvature (no sigma^-4 growth); stable multi-start
   behavior; predictive recovery under known simulation truth. Only then spend compute on the nested Gate A run.

**NEXT (fresh session):** draft the estimator-redesign specification (docs/plans/), status DRAFT, open items listed
for the PI. No implementation until it is frozen.

Open items the spec draft must list for the PI (coordinator notes, not decisions):
- F2 link: exact optimization of the published objective overfit (CV 0.638 -> 0.571; stall acted as implicit
  regularization). Qualification should therefore include held-out discrimination on the simulations at
  convergence, not only clean convergence, since a cleanly converging estimator that overfits still fails Gate A.
- Penalizing raw beta is scale-dependent: a factor with larger score variance pays less penalty per unit of effect.
  Tuning ranges for lambda (and nu) must be stated relative to the score scale (rank X, p = 3,000), and whether to
  rescale is a spec choice.
- The rho term moves from P W diag(1/sigma) b to P W beta; its scaling constant (currently sb2 / 2) needs restating.
- Comparator: confirm "current implemented estimator" means the production DeSurv fitting procedure (k = 3, tuned
  by fresh inner BO in every outer fold per docs/plans/2026-10-04-estimator-gate-spec.md), not the PALM candidate
  and not the published hyperparameters.

## PI ruling 2026-10-05, ninth (on redesign spec DRAFT revision 1, 8772b86)

Beta redesign approved in principle, not frozen. Optimize theta = s0 beta with s0 = sigma_bar (fixed, from training X
only, once per fit); binding F2 safeguard = converged vs frozen early checkpoints (mean paired dC >= -0.01 at each),
oracle gap descriptive; rho and range kept; comparator DeSurv 1.0.1 confirmed (1.2.0 sensitivity only); PALM with
maxit 5,000 as a work limit; unit sphere kept; qualification thresholds tied to the pathology, S3 = predictive
low-variance vs nonpredictive high-variance; no convergence-theorem claim, PALM assumptions checked explicitly
before implementation; score and loading recovery both assessed. Full record and the revision-2 open items 9a-9d:
`docs/plans/2026-10-05-estimator-redesign-spec.md` section 10.

**NEXT:** PI rules 9a-9d, then freeze. No implementation, real data or compute before freeze.

## PI ruling 2026-10-05, tenth (on spec DRAFT revision 2, 593a4d1)

9a: S3 bias test binding at predictor level (Q7, lower 95% MC bound of mean D_AB > 0; factor match descriptive).
9b: fail on mechanism recurrence, not historical magnitudes: raw sigma and factor-100 curvature descriptive; analytic
curvature bound and KKT binding; Q5 0.01 predictive binding, 0.80 loading descriptive. 9c: 0.25 primary, 0.50
optional reported. 9d: keep HALS/FISTA, decline a classical PALM theorem claim. Spec revision 3 applies these.

**NEXT:** independent Codex review of revision 3 with `docs/plans/2026-10-05-estimator-redesign-codex-prompt.md`. Attempted 2026-10-05: Codex usage limit, available again 2026-10-09 17:11. Re-run from the repo root:
`codex exec --sandbox read-only -c model_reasoning_effort=medium -o docs/plans/2026-10-05-estimator-redesign-codex-review.md "$(cat docs/plans/2026-10-05-estimator-redesign-codex-prompt.md)"`.
Then PI freeze, then implementation.

## PI ruling 2026-10-05, eleventh (on spec revision 3)

Wait for the independent Codex review (Oct 9); no freeze or implementation before it. Q3 reworded to bind only on the
mechanism (non-converged start with g-driven scale loss and a score preserved relative to the scale change, all
comparisons); historical D-Q4 magnitudes reported descriptively. No further substantive edits before the review.
