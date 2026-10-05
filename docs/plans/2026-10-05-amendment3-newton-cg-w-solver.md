# Amendment 3: active-face Newton-CG W solver (APPROVED, PI, 2026-10-05; FROZEN)

Amends `2026-10-04-estimator-gate-spec.md` section 4, "Optimizer", step 2, and the maxit clause of "Stopping".
Approved and frozen before any Newton-CG code exists and before any Newton-CG output has been seen. Sources: PI
rulings 2026-10-05 (first, second, third) in `docs/HANDOFF-2026-10-05-candidate-optimizer.md`, the PI's
mechanics specification of the same day, and the PI's ruling on the eight [PI] items of the draft (commit
c2f5482), recorded in section 9.

**Guiding principle (PI).** Make the numerical solver mathematically coherent without changing the objective,
the statistical estimator, the tuning space or the convergence certificate.

**After approval.** The amendment window is closed again. Bug fixes that make code conform to this text are
permitted. No change to the numerical mechanics is permitted on the basis of whether the candidate looks
favorable; any such change requires a new amendment under the failure rule (section 7).

## 1. Hierarchy

1. S1-S5 (`tests/out/coordinator/wcond.R`, commit c4b6db6) established the numerical pathology: the W block
   alone reproduces the stall at p = 3000; the Riemannian W-block Hessian on the active face has kappa 2.4e5 to
   5.0e5 (about 190 at the p = 150 control); 92-99% of the residual's spectral weight sits in the lowest-curvature
   directions.
2. Amendment 3 changes only the W optimizer.
3. Four warm-started p = 3000 end points are the binding qualification.
4. p = 150 is the estimator-equivalence and non-regression control.
5. Thirty-start realistic fits are the end-to-end check before the definitive recomputation.

## 2. Frozen (unchanged by this amendment)

Objective Psi and F, penalties, the H update (one HALS sweep), the b update (FISTA to KKT_b <= 1e-10 with
F-nonincrease acceptance), sweep order H, W, b, the feasible set S_+^k = {W >= 0, ||w_j|| = 1}, initialization and
seed schedule, the stopping certificate r = max(r_W, r_H, r_b) < 1e-5 under the amendment-2 definitions, the
F-increase definition, start statuses (converged, nonconverged, numerical failure) and the start-selection rule
(lowest-F converged start, else lowest-F nonconverged start, flagged), and every element of BO, scoring, Gate A
and Gate B.

## 3. The amended W step

Replaces spec section 4, "Optimizer", step 2. Everything is evaluated at fixed (H, b). G_W = grad_W Psi
(Euclidean).

**Sweep counter.** One sweep is one pass H, W, b. Each sweep counts once toward the per-start budget, whether the
W step is a Newton-CG update or a no-op. This is the only per-start sweep counter (section 4).

**Convergence guard.** At the start of the W block, compute r_W. If r_W < 1e-5 (the certificate's threshold), W is
a no-op for this sweep: W is kept, no CG or line search is run, exit reason `w_converged`. H and b proceed as
frozen. r_W is re-evaluated at the next sweep, after H and b have moved; if it then reaches the threshold, W is
active again.

**Gradient and geometry.** Newton-CG operates in the columnwise tangent space of the unit spheres, restricted to
the free coordinates. The gradient used is the projected (Riemannian) gradient matrix g = R_W restricted to the
free set, the matrix whose Frobenius norm defines the frozen r_W. The operator is the Riemannian Hessian already
used and validated in S1-S5:
Hess[V]_j = P_j((d2Psi V)_j) - (w_j' g_j) V_j,
with P_j the projection onto the tangent space of column j on the free face (zero off the free set, orthogonal
to w_j), g_j here the Euclidean gradient column G_W[, j] as in `wcond.R`, and d2Psi V the analytic Euclidean HVP.
After a trial step, W is retracted to nonnegative unit-norm columns (below).

**Analytic HVP.** d2Psi V is computed directly from every W-dependent term of Psi (reconstruction and survival
contributions, including the W dependence of sigma_j(W)), restricted to the free coordinates. No explicit Hessian
and no finite differences. If an analytic HVP proves infeasible, the implementer stops and reports to the PI; a
finite-difference substitute is not adopted without a ruling.

**Active face.** Free coordinates: W_ij > 0, plus W_ij = 0 with G_W[i,j] < 0 (the gradient points into the
feasible region, so the coefficient may re-enter). These are exactly the two cases in the definition of R_W. All
other coordinates are held at zero for the iteration.

**CG.** Plain truncated CG on Hess d = -g over the free coordinates, starting at d = 0, no preconditioner. Stop at
the first of:
- forcing: ||Hess d + g|| <= min(0.5, sqrt(||g||)) ||g||;
- 100 HVPs in this solve (the direction accumulated so far is used; the start continues);
- nonpositive or numerically flat curvature (below).

**Negative or flat curvature.** For CG search direction p_j, q_j = p_j' Hess p_j / p_j' p_j. Let q_ref be the
first finite positive Rayleigh quotient of this CG solve. Curvature is nonpositive or flat when
q_j <= 1e-12 max(1, |q_ref|); for j = 1, q_ref is q_1 when q_1 is finite and positive, so the j = 1 test is
q_1 <= 1e-12 max(1, |q_1|), which any nonpositive q_1 meets. On detection: if j = 1, terminate CG and use d = -g,
flagged `nc_first_iter = TRUE`; if j > 1, terminate CG and use the direction accumulated so far, flagged
`negative_curvature = TRUE`. The objective is not changed and no damping is added to make Hess positive
definite. A nonfinite q_j at any j is a nonfinite quantity and therefore a numerical failure.

**Retracted search path.** W(t) = N([W + t d]_+), where [.]_+ sets negative entries to zero and N renormalizes each
column to unit norm. Coefficients may hit zero inside the search and the algorithm continues along the boundary;
there is no straight-line feasibility cap. A trial t is rejected (and halved) if [W + t d]_+ annihilates an
entire column, since N is then undefined.

**Line search.** Monotone Armijo backtracking on the actual retracted objective at fixed (H, b): start at t = 1
(the Newton step), halve, and, with Delta W(t) = W(t) - W, accept the first nonzero trial with
<G_W, Delta W(t)> < 0  and  F(W(t)) <= F(W) + 1e-4 <G_W, Delta W(t)>  and  F(W(t)) < F(W).
In exact arithmetic the first two give the third; it is checked explicitly because in floating point
F(W) + 1e-4 <G_W, Delta W(t)> can round to F(W) when the slope is tiny. If <G_W, Delta W(t)> >= 0 (clipping and renormalization
have turned the displacement away from descent), the trial is rejected and t is halved; without this guard the
Armijo inequality could admit an objective increase. The sufficient-decrease term uses the actual displacement of the retracted point (the projected-search form of
Armijo, Bertsekas 1982), so clipping and renormalization are accounted for; for an unclipped step it equals
1e-4 t <g, d> to first order. No fallback to the old first-order step. Face changes (coordinates entering and
leaving the support) are recorded at every accepted step.

**Line-search failure.** The search terminates without acceptance after 60 halvings, or as soon as the retracted
trial W(t) is numerically identical to W (every entry equal in floating point). Either ends the start with status
*nonconverged*, reason `line_search_failed`; the start keeps its last accepted W and remains eligible for the
lowest-F nonconverged fallback. Nonfinite quantities anywhere remain *numerical failure*. The frozen W step floor
t_min = 1e-20 belonged to the replaced step and no longer applies.

**Theory statement.** The PALM convergence theorem (Bolte, Sabach and Teboulle 2014) assumed the frozen
prox-linear W step and its sufficient decrease ||W+ - W||^2/8. The Newton-CG step with Armijo decrease is not
covered by that theorem. The spec makes no convergence claim for the amended algorithm; convergence is certified
per start by the frozen residual r < 1e-5.

## 4. Budget and status

**Per-start budget.** At most 1,000 sweeps, at most 100 HVPs per W solve, at most 25,000 W-block HVPs in total.
All are per start; a 30-start fit has 30 independent budgets. The 25,000 total is the binding work cap; it is not
relaxed because 1,000 x 100 is larger. The 100-HVP ceiling only truncates one CG solve. HVPs used by the Q0 test
and by diagnostics do not count.

**maxit superseded.** For the candidate with the amended W solver, the per-start budget replaces maxit = 5000 in
spec section 4, "Stopping". There is one binding per-start sweep rule, not two counters.

**Status.** Production: a start that reaches the per-start cap (1,000 sweeps or 25,000 total HVPs) is
*nonconverged*, reason `work_cap_reached`; the frozen start-selection rule is unchanged. Qualification: reaching a
per-start cap is a failure.

## 5. Diagnostics recorded per start

Per sweep: F, r, r_W, r_H, r_b, W exit reason (`w_converged` or Newton step), face size, CG iterations, HVPs, CG
exit reason (forcing, HVP ceiling, curvature), `negative_curvature` and `nc_first_iter` flags, accepted t,
halvings, coordinates entering and leaving the support. Per start: wall time, total sweeps, total HVPs, objective
path, final KKT residuals, status and reason. Wall time is recorded, never a criterion.

## 6. Qualification

**Q0. Analytic HVP test (precondition, simulated data).** On each T1 simulated problem, for 10 random directions
V on the free face (seeded, unit Frobenius norm), compare d2Psi V with
(grad_W F(W + eps V) - grad_W F(W - eps V)) / (2 eps) on the fixed grid eps = 1e-2, 1e-3, ..., 1e-8. Relative
error e(eps) = ||fd - analytic||_F / ||analytic||_F. Pass, for every direction:
- min over the grid of e(eps) <= 1e-6; and
- empirical order log10(e(eps_k) / e(eps_k+1)) >= 1.5 over at least two consecutive refinements before roundoff
  (the expected second-order decay of the central difference).
Finite differences test the HVP; they are not used to implement it. The assembled Riemannian operator Hess is
tested the same way, against the central difference composed with the tangent projection and the Weingarten term,
as in `wcond.R`.

**Q-T. Full simulated requalification (first ruling, item 4).** Rerun T1-T12 with the amended W step:
gradients, zero material F increases, the same independent KKT oracle, recovery, timing, edge cases, mutation
tests, fold-locality and transport. Every test passes under its existing criterion; none is relaxed. T7 is the
test whose p = 3000 failure motivated this amendment; its outcome is reported, and Q1, not T7 timing, is the
binding realistic-size criterion.

*Replaced-mechanism assertions.* Assertions that test mechanisms this amendment removes are replaced one for one by
the amended counterpart; every other criterion is unchanged.

| Frozen assertion (test) | Removed mechanism | Amended counterpart |
|---|---|---|
| accepted W steps <= 2 (T4) | step cap t_max = 2 | accepted t <= 1 on the retracted path; every accepted step satisfies the Armijo condition of section 3 |
| W step floor -> numerical failure, W unchanged (T9) | t_min = 1e-20 floor | 60 halvings or no-change guard -> *nonconverged*, `line_search_failed`, W unchanged |
| maxit reached -> nonconverged (T9), maxit 5000 (T4, T7) | maxit = 5000 | per-start cap -> *nonconverged*, `work_cap_reached`; T4 and T7 run to the per-start budget |
| W step cap test (T9) | t_max = 2 | trial annihilating a column is rejected and halved |
| mutant B7, accept at the floor (T10) | floor | mutant accepting the last trial after 60 halvings must be killed |
| mutant B19_cap, W step uncapped (T10) | cap | mutant B19_retraction, Armijo slope computed as t <G_W, d> instead of <G_W, W(t) - W>, must be killed |

Added assertions: W no-op when r_W < 1e-5 at the start of the block; `nc_first_iter` path uses -g.

**Q1. Binding realistic-size qualification.** Warm-start from each of the four stalled p = 3000, n = 218 end
points (starts 1, 5, 9, 16), full sweeps (H, W, b) with the amended W step. All four (the prespecified proportion
is 4 of 4) must:
1. reach the unchanged certificate r < 1e-5;
2. respect monotonicity (no F increase under the frozen scale-aware definition);
3. stay within the per-start budget.
Reaching a per-start cap, a line-search failure or a numerical failure in any of the four fails Q1; a CG solve
truncated at its 100-HVP ceiling does not. F_final <= F_stalled_endpoint follows from criterion 2 (every block is
monotone from the warm start) and is checked as part of it; the size of any improvement is reported, not
required. Loading and prediction agreement with the stalled end point are reported descriptively only.

**Q2. Estimator equivalence and non-regression at p = 150.** On the p = 150 control settings, both solvers run
the frozen 30-start fit from the same initializations. Qualification requires:
- both solvers converge under the same r < 1e-5 certificate;
- the same start is selected;
- relative objective difference |F_new - F_old| / max(1, |F_old|) <= 1e-6;
- after optimal (Hungarian) column matching, every W-column cosine similarity >= 0.999;
- fitted standardized risk scores eta~: Pearson r >= 0.9999 and
  RMSE(eta~_new - eta~_old) / SD(eta~_old) <= 0.01.
H and b differences are reported diagnostically and are not separate qualification gates. Rationale (PI): both
solvers stop at KKT 1e-5 in an ill-conditioned model, so two converged points can share objective and predictions
while differing in parameter space by much more than 1e-6. These criteria are invariant or near-invariant to that
flat-basin freedom; a genuinely different solution would fail them.

**Q3. End-to-end check (not the amendment-3 gate).** Full frozen 30-start p = 3000 fits with the amended W step.
Reported, not pass/fail for amendment 3; must show sensible multi-start behavior (status counts, F spread across
starts, selected start) before the amended solver is used for the definitive recomputation. No numerical
threshold for Q3 is fixed by this amendment.

## 7. Failure rule

If Q0, Q-T, Q1 or Q2 fails, Newton-CG is not tuned against T7, Q1 or Q2 outcomes. Diagnose, and write another
amendment before changing the algorithm. A preconditioner enters only through a later amendment and only with
evidence that CG itself has become the bottleneck.

## 8. Chronology (first ruling, item 5; verbatim for the spec and any write-up)

"The statistical specification was frozen before implementation. Prespecified simulated-data qualification
subsequently showed that the reference W-block optimizer failed to reach the frozen KKT tolerance at realistic
dimensionality despite monotone objective decrease. Before any real-data gate analysis, the specification was
amended once to replace only the W-block numerical solver; the objective, convergence certificate, tuning
procedure, and gate rules were unchanged."

## 9. Record of the PI ruling on the draft's [PI] items (2026-10-05)

| Item | Draft proposal | Ruling |
|---|---|---|
| 1 Sphere | Tangent-space Newton system | Approved; wording fixed: g is the projected gradient matrix whose norm is r_W, not the scalar |
| 2 Face | {W > 0} union {W = 0, G_W < 0} | Approved |
| 3 Jamming | Keep straight-line cap, drop freed coordinates with d < 0 | **Changed**: retracted path N([W + t d]_+), no cap; reject only column annihilation; Armijo on the retracted objective; record face changes |
| 4 Line-search failure | 60 halvings; nonconverged | Approved, plus no-change guard (trial identical to W); nonfinite remains numerical failure |
| 5a HVP test | min relative error <= 1e-6 | Approved, plus deterministic eps grid and order >= 1.5 over two consecutive refinements |
| 5b p = 150 equivalence | Parameter-level tolerances (1e-8 F, 1e-6 cosine, 1e-4 H and b) | **Changed**: invariant criteria in Q2; H and b diagnostic only |
| 6 maxit | Budget supersedes maxit = 5000 | Approved; one sweep counter |
| 7 Curvature | 1e-12 relative to first quotient | **Changed**: scale-aware Rayleigh test relative to the first finite positive quotient, floor 1; first-iteration hit uses -g |
| 8 W after r_W converges | Always attempt W | **Changed**: W no-op when r_W is below threshold at the start of the block |

Also ruled with the items: the p = 3000 qualification is the certificate, the prespecified proportion, monotonicity
and the budget; a lower F is evidence, not the criterion. This replaces "F_final <= F_stalled_endpoint" as a
separate binding criterion from the third ruling; it is retained as an implied check under monotonicity.

Post-approval consistency corrections (PI, 2026-10-05, not a reopening): (a) the T4/T9/T10 remapping in Q-T,
confirmed; (b) the Armijo acceptance requires <G_W, Delta W(t)> < 0 for a nonzero trial, closing a hole in the
approved rule; (c) the B19 replacement mutant is B19_retraction (slope t <G_W, d>), the implementation error
specific to the retracted path. The coordinator choices below were confirmed. Amendment 3 is closed.

Implementation choices made by the coordinator to write the rulings down, open to the PI's correction before code
is written: the strict "<" in the W no-op guard (the certificate uses r < 1e-5); the projected-search Armijo term
<G_W, W(t) - W>; t starting at 1; the Q0 grid (1e-2 to 1e-8), 10 directions per T1 problem and unit-norm V; Q1's
proportion stated as 4 of 4 from the third ruling; Q2 applied to the frozen 30-start fit so that "same start
selected" is defined.
