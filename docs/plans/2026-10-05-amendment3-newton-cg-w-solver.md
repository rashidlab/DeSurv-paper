# Amendment 3: active-face Newton-CG W solver (DRAFT, 2026-10-05, not yet approved)

Amends `2026-10-04-estimator-gate-spec.md` section 4, "Optimizer", step 2 only. Written before any Newton-CG code
exists and before any Newton-CG output has been seen. Sources: PI rulings 2026-10-05 (first, second, third) in
`docs/HANDOFF-2026-10-05-candidate-optimizer.md`, and the PI's mechanics specification of the same day. The
amendment window closed after amendment 2; this amendment is opened on the PI's explicit instruction and is
limited to the problem S1-S5 established.

Items marked **[PI]** are places where the mechanics specification cannot be implemented literally against the
frozen constraint set or the frozen spec. Each carries a proposed resolution. Nothing marked [PI] is final until
ruled. All other text restates rulings already made.

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
F-nonincrease acceptance), sweep order H, W, b, normalization S_+^k = {W >= 0, ||w_j|| = 1}, initialization and
seed schedule, the stopping residual r = max(r_W, r_H, r_b) < 1e-5 under the amendment-2 definitions, the
F-increase definition, start statuses and the start-selection rule (lowest-F converged start, else lowest-F
nonconverged start, flagged), and every element of BO, scoring, Gate A and Gate B.

## 3. The amended W step

Replaces spec section 4, "Optimizer", step 2. Everything is evaluated at fixed (H, b).

**Outer iteration.** One complete W-block update at fixed (H, b): identify the current active face, solve for a
Newton-CG direction on its free coordinates, perform the line search, accept one W, then recompute the frozen
KKT/residual quantities. One outer iteration is one W update. Because the PALM sweep makes exactly one W update,
outer iterations and sweeps coincide one to one (see [PI-6]).

**Geometry [PI-1].** Each column w_j lies on the nonnegative part of the unit sphere. The specification's g_F and
H_W, read literally as Euclidean quantities on the free coordinates, do not vanish or become meaningful at a
stationary point. The Euclidean gradient keeps a radial component (w_j' g_j) w_j that the frozen residual r_W
removes, so a forcing term built on the Euclidean ||g_F|| never tends to zero. Proposed: pose the Newton system
in the tangent space of each sphere on the free face, which is the operator S1-S5 measured:
- g_F = the binding residual R_W restricted to the free coordinates (tangent-projected gradient);
- Hess[V]_j = P_j((d2Psi V)_j) - (w_j' g_j) V_j, where P_j projects onto the tangent space of column j on the
  free face and d2Psi V is the analytic Euclidean HVP;
- retraction = the existing column normalization (the frozen P_S+ applied to a nonnegative matrix).
The Weingarten term -(w_j' g_j) V_j can make Hess indefinite; the negative-curvature rule below handles that.

**Analytic HVP.** d2Psi V is computed directly from every W-dependent term of Psi (reconstruction and survival
contributions, including the W dependence of sigma_j(W)), restricted to the free coordinates. No explicit Hessian
and no finite differences. If an analytic HVP proves infeasible, the implementer stops and reports to the PI; a
finite-difference substitute is not adopted without a ruling.

**Active face [PI-2].** The specification leaves the face rule open, and a face equal to {W > 0} can never release
a zero coordinate. At the p = 3000 end points, S4 found W = 0 entries with negative gradient carrying a
nonzero (under 50%) share of ||R_W||^2. With a face fixed at {W > 0}, r_W cannot reach 1e-5 whenever such entries
exist. Proposed: free set = {W_ij > 0} union {W_ij = 0 and G_W[i,j] < 0}, matching exactly the two cases in the
definition of R_W; all other coordinates are held at zero for the iteration.

**CG.** Plain CG on Hess d = -g_F over the free coordinates, starting at d = 0, with no preconditioner. Stop at
the first of: relative residual ||Hess d + g_F|| <= min(0.5, sqrt(||g_F||)) ||g_F|| (forcing sequence), 100 HVPs
in this solve, or negative curvature.

**Negative curvature [PI-7].** If a CG search direction p has p' Hess p <= 1e-12 s ||p||^2, with
s = |Rayleigh quotient of the first CG direction| (recorded), truncate CG and use the direction accumulated so
far. Record `negative_curvature = TRUE` for that outer iteration. The objective is not changed and no damping is
added to make Hess positive definite. Gap: if negative curvature occurs on the first CG iteration the accumulated
direction is zero. Proposed: use the steepest-descent direction -g_F in that case (the standard Steihaug choice),
flagged separately as `nc_first_iter = TRUE`.

**Feasibility and step [PI-3].** The specification caps the trial step at the largest step that preserves the W
constraints, then applies the existing normalization. A nonzero straight-line step never stays on the sphere, so
the cap can only be taken against nonnegativity: t_max = the largest t with W + t d >= 0 on the free set (capped
at 1, the Newton step). Two consequences to rule on:
- a freed zero coordinate with d_ij < 0, or a free coordinate at 1e-12, gives t_max at or near zero, the
  "jamming" failure of capped active-set methods. At p = 3000 many W entries are near zero, so this is a live
  risk of re-creating a stall by a different route;
- the alternative used by projected Newton methods (Bertsekas 1982) searches along the projected path
  W(t) = P_S+(W + t d), using only the existing projection, and lets coordinates hit zero inside the search.
Proposed: keep the PI's cap (it is the more conservative change) and, in addition, exclude from the free set any
freed zero coordinate whose CG direction is negative (it is returned to the active set before t_max is
computed). Record t_max and the number of coordinates that hit zero at each iteration so jamming is visible. If
the PI prefers the projected path, the line search below applies unchanged to W(t).

**Line search.** Monotone Armijo backtracking on the full objective F at fixed (H, b), starting at t = t_max,
halving, accepting the first t with F(R(W + t d)) <= F(W) + 1e-4 t <g_F, d>, where R is the column normalization.
Actual decrease is required: F(W_new) < F(W). No fallback to the old first-order step. If g_F = 0 exactly
(the W block is stationary on its face), the W update is a no-op: W is kept, no line search is run, and the
outer iteration is recorded with exit reason `w_stationary`. The H and b blocks proceed as frozen.

**Line-search failure [PI-4].** The specification makes a line-search failure "nonconverged". Two gaps:
- a floor is needed. Proposed: 60 halvings (t below about 1e-18 t_max), matching the scale of the frozen
  t_min = 1e-20;
- in the frozen spec, the step floor is a *numerical failure*, not nonconverged, and the third ruling says
  production status logic is unchanged. Proposed: a W line-search failure ends the start with status
  *nonconverged*, reason `line_search_failed`; it is not a numerical failure because all values are finite and F
  did not increase. The start keeps its last accepted W and remains eligible for the lowest-F nonconverged
  fallback. In the qualification benchmark it is a qualification failure.
A related risk, not a decision: near stationarity, predicted decreases can fall below the F-increase tolerance
1e-13 max(1, |F|), so Armijo can fail on rounding before r < 1e-5. If this happens, the qualification fails and
is reported under the failure rule (section 6); the amendment does not prespecify a remedy. The same applies
once r_W < 1e-5 while r_H or r_b has not yet converged: W steps then predict decreases near rounding, and a
line-search failure there would end the start before H and b finish. Skipping the W update when r_W < 1e-5 would
avoid this but changes the sweep, so it is listed as [PI-8] rather than adopted.

**Budget (per start).** At most 1,000 outer iterations, at most 100 HVPs per W solve, at most 25,000 W-block
HVPs in total. All three are per start; a 30-start fit has 30 independent budgets. The 25,000 total binds; it is
not relaxed because 1,000 x 100 is larger. HVPs used by the directional-derivative test and by diagnostics do
not count.

**Status [PI-6].** The 100-HVP ceiling only truncates one CG solve: the direction accumulated so far goes to
the line search and the start continues. Production: a start that reaches the per-start cap (1,000 outer
iterations or 25,000 total HVPs) is *nonconverged*, reason `work_cap_reached`; the frozen start-selection rule is unchanged. Qualification: reaching a per-start cap is a failure.
Consequence to rule on: the frozen maxit = 5000 sweeps cannot bind while the 1,000-outer-iteration cap exists,
because each sweep makes one W update. In production the effective sweep limit becomes 1,000. Proposed: state
this explicitly in the spec ("maxit = 5000 is superseded by the per-start Newton-CG budget"), since it changes
which starts are nonconverged relative to the frozen spec.

**Theory statement.** The PALM convergence theorem (Bolte, Sabach and Teboulle 2014) assumed the frozen
prox-linear W step and its sufficient decrease ||W+ - W||^2/8. The Newton-CG step with Armijo decrease is not
covered by that theorem. The spec makes no convergence claim for the amended algorithm; convergence is certified
per start by the frozen residual r < 1e-5.

## 4. Diagnostics recorded per start

Per outer iteration: F, r, r_W, r_H, r_b, face size, CG iterations, HVPs, CG exit reason (forcing, cap,
negative curvature), negative_curvature flag, t_max, accepted t, halvings, coordinates hitting zero. Per start:
wall time, total outer iterations, total HVPs, objective path, final KKT residuals, status and reason. Wall time
is recorded, never a criterion.

## 5. Qualification

**Q-T. Full simulated requalification (first ruling, item 4).** Rerun T1-T12 with the amended W step:
gradients, zero material F increases, the same independent KKT oracle, recovery, timing, edge cases, mutation
tests, fold-locality and transport. Every test passes under its existing criterion; none is relaxed. T7 is the
test whose p = 3000 failure motivated this amendment; its outcome is reported, and Q1 below, not T7 timing, is the
binding realistic-size criterion.

**Q0. Analytic HVP test (precondition, simulated data, small problems).** For random directions V on the free
face, compare d2Psi V with (grad_W F(W + eps V) - grad_W F(W - eps V)) / (2 eps). Finite differences test the
HVP; they are not used to implement it. Run at several eps and require the error to show the expected O(eps^2)
decay before rounding takes over. Separately test the assembled Riemannian operator Hess against the same
construction composed with the tangent projection and the Weingarten term, as in `wcond.R`. Pass threshold
**[PI-5a]**: proposed minimum over eps of the relative error <= 1e-6, at all tested directions.

**Q1. Binding realistic-size qualification.** Warm-start from each of the four stalled p = 3000, n = 218 end
points (starts 1, 5, 9, 16), full sweeps (H, W, b) with the amended W step. Each must, within the per-start
budget:
- reach the frozen criterion r < 1e-5;
- satisfy F_final <= F_stalled_endpoint (no minimum improvement prespecified).
Loading and prediction agreement with the stalled end point are reported descriptively only.
Any cap, line-search failure or numerical failure in any of the four fails Q1.

**Q2. Estimator equivalence and non-regression at p = 150.** On the converged p = 150 control settings, from the
same initializations, the amended and frozen solvers must agree on: convergence status; F within tolerance;
fitted W, H, b within tolerance; predictions within tolerance; and, for the 30-start control fits, the same
selected start (first ruling, item 4(ii)). **[PI-5b]** The third ruling refers to
"already-defined equivalence tolerances". None are defined in the gate spec, the ruling, the handoff, or the
candidate tests. Proposed, fixed now and before either solver is run for Q2:
- F: |F_new - F_old| <= 1e-8 max(1, |F_old|);
- W, H: after Hungarian column matching (identity expected), max column-wise 1 - cosine <= 1e-6 for W, and
  ||H_new - H_old||_F / ||H_old||_F <= 1e-4;
- b: ||b_new - b_old||_inf <= 1e-4 max(1, ||b_old||_inf);
- predictions (linear predictor on the control's held-out samples): Pearson r >= 0.9999 and
  max |difference| <= 1e-4 times the predictor's SD.
Rationale: both solvers stop at r < 1e-5, so the solutions agree only to the accuracy that tolerance implies
along well-conditioned directions; the F tolerance is tighter because F is second-order in the residual.

**Q3. End-to-end check (not the amendment-3 gate).** Full frozen 30-start p = 3000 fits with the amended W step.
Reported, not pass/fail for amendment 3; must show sensible multi-start behavior (status counts, F spread
across starts, selected start) before the amended solver is used for the definitive recomputation.

## 6. Failure rule

If Q0, Q1 or Q2 fails, Newton-CG is not tuned against T7 or Q1 outcomes. Diagnose, and write another amendment
before changing the algorithm. A preconditioner enters only through a later amendment and only with evidence
that CG itself has become the bottleneck.

## 7. Chronology (first ruling, item 5; verbatim for the spec and any write-up)

"The statistical specification was frozen before implementation. Prespecified simulated-data qualification
subsequently showed that the reference W-block optimizer failed to reach the frozen KKT tolerance at realistic
dimensionality despite monotone objective decrease. Before any real-data gate analysis, the specification was
amended once to replace only the W-block numerical solver; the objective, convergence certificate, tuning
procedure, and gate rules were unchanged."

On approval, the amendment window closes again. Amendment 3 is appended to the spec's Amendments list with its
approval date.

## 8. Open items for the PI

| Item | Question | Proposed |
|---|---|---|
| PI-1 | Euclidean or tangent-space Newton system? | Tangent space (the operator S1-S5 measured); g_F = R_W on the free face |
| PI-2 | Face rule | {W > 0} union {W = 0, G_W < 0}, the two cases of R_W |
| PI-3 | Step cap vs projected path | Keep the cap; drop freed zero coordinates with d < 0; record jamming metrics |
| PI-4 | Line-search floor and status | 60 halvings; status nonconverged, reason `line_search_failed` |
| PI-5a | HVP test threshold | min-over-eps relative error <= 1e-6 |
| PI-5b | p = 150 equivalence tolerances (not defined anywhere yet) | as in Q2 |
| PI-6 | maxit = 5000 vs 1,000-iteration cap | state that the per-start budget supersedes maxit |
| PI-8 | W update once r_W < 1e-5 but r_H or r_b is not | Keep the frozen sweep (always attempt W); report any rounding-driven line-search failure |
| PI-7 | Negative-curvature tolerance; first-iteration case | 1e-12 relative to first Rayleigh quotient; -g_F, flagged |
