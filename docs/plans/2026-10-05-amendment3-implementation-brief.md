# Amendment 3 implementation brief: analytic W-block HVP (Stage A, Q0) and Newton-CG W step (Stage B)

Coordinator brief, 2026-10-05. Implements `docs/plans/2026-10-05-amendment3-newton-cg-w-solver.md` (APPROVED,
FROZEN, window closed). That document is the authority. Where this brief and the frozen text disagree, the frozen
text wins: **report the conflict, do not improvise.** Where this brief and the code you find disagree (function
names, argument order, data layout), report it before working around it.

Working tree: `/home/naimrashid/Downloads/DeSurv-paper-recompute`, branch `recompute/consensus-desurv`. Code lives in
`analysis/estimator-gate/candidate/`. Simulated data only; nothing reads real data.

Standing rules for every stage:
- No git commits, no pushes, no edits to anything under `docs/` (the coordinator commits).
- No finite differences inside any implementation function. Finite differences appear only in test code.
- No change to existing functions in Stage A (additions only). Stage B changes are limited to the inventory below.
- Lab style: base R, no tidyverse. Comment density matches `palm_lib.R`.
- Run R as `Rscript` (the tests set `.libPaths` themselves via `tests/common.R`).
- Your final reply is at most 150 words: verdicts, counts, files written, and any conflict or deviation. The full
  account goes in the report file named in each stage.

---

## Stage A: analytic HVP and Q0 (dispatch now)

### A1. Mathematics (evaluate at fixed H, b; V is a p x k direction)

Notation follows `palm_lib.R`: X is p x n, C centers over samples, P centers over genes (`sweep(W, 2, colMeans(W))`).

- Zc = C X'W, Yc = C X'V (`Y <- crossprod(D$X, V); Yc <- sweep(Y, 2, colMeans(Y))`).
- sigma_j = sqrt(||Zc_j||^2/(n-1) + eps^2); dsigma_j = <Zc_j, Yc_j> / ((n-1) sigma_j).
- beta = b / sigma; dbeta_j = -beta_j dsigma_j / sigma_j.
- Every survival term has the gradient form already used in `grad_W_terms`:
  G_term = q beta' - (X Zc) diag(s),  s_j = c_j beta_j / ((n-1) sigma_j^2),  c = W'q.
  Its directional derivative along V is
  dG = dq beta' + q dbeta' - (X Yc) diag(s) - (X Zc) diag(ds),
  ds_j = (dc_j beta_j + c_j dbeta_j) / ((n-1) sigma_j^2) - 2 s_j dsigma_j / sigma_j,
  dc_j = <V_j, q> + <W_j, dq>.
- Cox: q = X r, r = -(2/ne) g(eta), eta = Zc beta. deta = Yc beta + Zc dbeta; dr = -(2/ne) Hll[deta]; dq = X dr.
- Breslow Hessian-vector product (exact ties, same tie groups `D$gs`, `D$ge` and sort `D$o` as `breslow()`):
  with e_l = exp(eta_l), S_i = sum_{m: y_m >= y_i} e_m, h_l = sum_{events i: y_i <= y_l} d_i / S_i,
  T_i = sum_{m: y_m >= y_i} e_m v_m,
  Hll[v]_l = -e_l v_l h_l + e_l sum_{events i: y_i <= y_l} d_i T_i / S_i^2.
  Fast path: the same single shift as `breslow()` (the products e_l h_l and e_l T_i/S_i^2 are shift invariant).
  Underflow path: mirror `breslow()`'s log-sum-exp fallback; T_i is signed, so split v into positive and negative
  parts and log-sum-exp each. Use the fast path under exactly the condition `breslow()` uses.
- rho-ridge: q = rho sb2 u, u = Wc beta (Wc = P W); du = Vc beta + Wc dbeta (Vc = P V); dq = rho sb2 du.
- Reconstruction: dG_recon = V H H' / ||X||^2.
- Weighted Euclidean HVP: d2Psi V = (1 - alpha) dG_recon + alpha (dG_cox + dG_ridge).
- Free face (amendment section 3): F = (W > 0) | (W == 0 & G_W < 0), G_W = `palm_grad_W` at the current point.
  Pface(V) = V * F, then per column subtract w_j <w_j, (V*F)_j>.
  Projected gradient g = Pface(G_W); check that ||g||_F equals `PALM_DEFS$resid_W` (it must, to rounding).
- Riemannian operator on the tangent face: Hess[V] = Pface(d2Psi V) - V diag(mu), mu_j = <w_j, G_W[, j]>,
  for V already in the tangent face. Same assembly as `make_op` in `tests/out/coordinator/wcond.R`, with the face
  widened to include W == 0 & G_W < 0.

If any step of this derivation is wrong or an analytic HVP proves infeasible, **stop and report**; do not
substitute finite differences.

### A2. Code (additions to `palm_lib.R` only; no existing line changed or deleted)

Place a new section after the gradients. Suggested names (report if you choose differently and why):
- `breslow_hvp(eta, D, v)`: Hll[v], both paths; an internal argument may force the log-sum-exp path for testing.
- `hvp_W_terms(W, H, b, D, hp, V)`: unweighted list(recon, cox, ridge), mirroring `grad_W_terms`.
- `palm_hvp_W(W, H, b, D, hp, V)`: weighted Euclidean HVP.
- `w_face(W, GW)`: logical free mask. `face_tangent(V, W, free)`: Pface.
- `riem_hess_W(W, H, b, D, hp, V, free, GW)`: Hess[V].

### A3. Test `tests/Q0_hvp.R` (sources `tests/common.R`; writes `tests/out/Q0_hvp.txt` and `.rds`; ends with `verdict()`)

Binding Q0, exactly as frozen:
- Problems: the 270 T1 Part A points (54 grid corners x 5 points), regenerated with T1's seeds and code so each
  (W, H, b, hp) is identical to T1's. Data: `sim_data(11, 80, 50)`.
- Directions: 10 per problem, seeded (`set.seed(700000 + 100 * gi + pt)` then draw the 10), restricted to the free
  face, unit Frobenius norm.
- Euclidean test: fd(eps) = (palm_grad_W(W + eps V) - palm_grad_W(W - eps V)) / (2 eps) on
  eps = 1e-2, 1e-3, ..., 1e-8; e(eps) = ||fd - analytic||_F / ||analytic||_F.
- Riemannian test: at the same problem with W replaced by its column-normalized version N(W) (the operator is
  defined on the sphere; coordinator choice, see A5), V drawn in the tangent face at N(W), compare `riem_hess_W`
  against Pface(fd(eps)) - V diag(mu), same grid, same e(eps).
- Pass, for every direction: min_eps e(eps) <= 1e-6, AND there exist two consecutive refinements k, k+1 with
  log10(e_k / e_{k+1}) >= 1.5 and log10(e_{k+1} / e_{k+2}) >= 1.5 (coordinator's operational reading of "over at
  least two consecutive refinements before roundoff").
- Q0 PASS iff both tests pass on all 2,700 directions each.

Report-only diagnostics (not criteria; never used to excuse a binding failure):
- Per-term unweighted checks (recon, cox, ridge) against FD of the matching `grad_W_terms` component, same grid.
- Symmetry: |<U, Hess V> - <V, Hess U>| / (||U|| ||Hess V||) for one pair per problem.
- A boundary-face point: one p = 150 problem with W containing zeros and some (W == 0, G_W < 0) coordinates, so the
  widened face is exercised.
- `breslow_hvp` fast path vs forced log-sum-exp path on points where both are valid; one point with predictor
  contrast large enough to trigger the fallback (checked against FD of `breslow()$g`).
- g identity: max |‖Pface(G_W)‖_F - resid_W| over the 270 points.

Predicted risk, stated before any output: the reconstruction gradient is linear in W, so its central difference has
no truncation error. Where reconstruction dominates (alpha = 1e-4), e(eps) may sit at roundoff over the whole grid
and the order clause may fail although the HVP is correct. Apply the criterion literally. Classify every failing
direction as (a) min e > 1e-6, (b) order clause failed with e(eps) <= 1e-9 at every eps (flat, roundoff-limited),
or (c) other; give counts by class and by alpha x rho. Do not change the grid, the directions or the criterion.

Bugs found by any check (binding or diagnostic) may be fixed so the code matches A1; list every fix in the report.

### A4. Verification-by-counting (put each count in the report)

1. `git diff --numstat analysis/estimator-gate/candidate/palm_lib.R` shows 0 deleted lines.
2. Q0 counts: 2,700 Euclidean and 2,700 Riemannian directions evaluated; the number passing each; the
   classification counts above.
3. No finite difference in implementation code: `grep -n 'eps \* V\|h \* V\|/ (2 \*' palm_lib.R` hits only lines
   that existed before (list them).
4. Regression: rerun `tests/T1_gradients.R` and `tests/T2_objective.R`; both still `VERDICT: PASS`.

Report file: `analysis/estimator-gate/candidate/AMENDMENT3_REPORT.md`, section "Stage A". Include wall time.

### A5. Coordinator choices in Stage A (open to PI correction)

- Riemannian Q0 at N(W): T1 points are ambient (off the sphere); the Weingarten term is only meaningful on it.
- The seed scheme for directions and the operational reading of the order clause (A3).
- T1 points have W > 0 everywhere, so the binding Q0 never exercises W == 0 free coordinates; the boundary-face
  point is a diagnostic, not a criterion, because the frozen text names T1 problems.

**Stop after Stage A.** Do not start Stage B.

---

## Stage B: Newton-CG W step and requalification (dispatch only after the coordinator verifies Q0 PASS)

### B1. Algorithm (amendment section 3; implement literally)

New `w_step_ncg(W, H, b, D, hp, hvp_left)`, at fixed (H, b), G_W at the current W:
1. r_W = `PALM_DEFS$resid_W(W, G_W)`. If r_W < 1e-5: return W unchanged, exit `w_converged`, 0 HVPs.
2. free = `w_face(W, G_W)`; g = Pface(G_W).
3. CG on Hess d = -g, d0 = 0, residual r0 = g, p0 = -g. Per iteration: Hp = Hess[p] (1 HVP); q = <p, Hp> / <p, p>;
   q_ref = first finite positive q of this solve; if q nonfinite -> numerical failure; if
   q <= 1e-12 max(1, |q_ref|) (for j = 1 the test is q_1 <= 1e-12 max(1, |q_1|)): j = 1 -> d = -g,
   `nc_first_iter`; j > 1 -> keep accumulated d, `negative_curvature`. Else standard CG update; stop on forcing
   ||r|| <= min(0.5, sqrt(||g||)) ||g|| (r the recursively updated residual, coordinator choice B5), or 100 HVPs in
   this solve, or the per-start HVP total reaching its cap.
4. Retracted path W(t) = N([W + t d]_+); t = 1, halve. Reject a trial that annihilates a column. Accept the first
   nonzero trial with <G_W, W(t) - W> < 0 and F(W(t)) <= F(W) + 1e-4 <G_W, W(t) - W> and F(W(t)) < F(W).
   Fail (`line_search_failed`) after 60 halvings or when W(t) is identical to W in every entry.
5. Record per sweep every field in amendment section 5 (face size, CG iterations, HVPs, CG exit reason, flags,
   accepted t, halvings, coordinates entering and leaving the support).

### B2. Breakage inventory (every existing object Stage B touches; anything else touched is a deviation to report)

| Object | Change |
|---|---|
| `palm_fit_one` | new arg `w_solver = c("ncg", "pg")`; under `ncg`, `maxit` is the sweep cap (default 1000L), add HVP total cap 25,000; statuses: cap -> nonconverged `work_cap_reached`, line-search failure -> nonconverged `line_search_failed`; nonfinite -> failure; per-sweep diagnostics into `trace` |
| `palm_fit` | pass `w_solver` through; default `ncg`; summary gains total HVPs and W exit counts |
| `w_step`, `PALM_DEFS$step_init` | unchanged; reachable only with `w_solver = "pg"` (needed by Q2 and by `wcond.R`) |
| `tests/mutants.R` | B7 replaced by a mutant accepting the last trial after 60 halvings; B19_cap replaced by B19_retraction (slope t <G_W, d>); others unchanged |
| `tests/T4_monotone.R` | "accepted W steps <= 2" -> "accepted t <= 1 and every accepted step satisfies Armijo"; maxit 5000 -> per-start budget |
| `tests/T7_timing.R` | maxit 5000 -> per-start budget; outcome reported, Q1 binding |
| `tests/T9_edge_cases.R` | floor assertion -> 60 halvings / no-change guard -> nonconverged `line_search_failed`, W unchanged; cap test -> column-annihilation rejection; maxit -> `work_cap_reached`; `step_init` assertion moved under `w_solver = "pg"`; added: W no-op at r_W < 1e-5, `nc_first_iter` uses -g |
| `tests/T10_mutation.R` | mutant table rows for B7 and B19_cap replaced per the amendment Q-T table |

Every other criterion in T1-T12 is unchanged. Then Q1 (`tests/Q1_warm.R`: warm-start from `T7_fit.rds`
`fit$fits[[i]]` for i in 1, 5, 9, 16, carrying W, H, b, eps, hp), Q2 (`tests/Q2_p150.R`), Q3 (report only), each as
frozen in amendment section 6.

### B3. Coordinator choices in Stage B (raise with the PI before dispatch)

- B5 forcing uses the recursively updated CG residual (equal to ||Hess d + g|| in exact arithmetic; saves one HVP
  per solve).
- When the per-start HVP total reaches 25,000 inside a solve: CG stops, the line search on the accumulated direction
  runs (it needs no HVPs), the sweep completes, the certificate is checked; converged if r < 1e-5, else
  nonconverged `work_cap_reached`.
- `maxit` keeps its name as the single sweep counter under `ncg` (tests that pass small `maxit` values are harness
  shortcuts).

---

## PI ruling 2026-10-05 on A5 and B3

All six coordinator choices approved as written: A5.1 (Riemannian Q0 at N(W)), A5.2 (order clause as two
consecutive 10-fold refinements each >= 1.5), A5.3 (zero-coordinate point diagnostic only; binding Q0 stays the
frozen T1 set); B3.1 (recursive CG residual), B3.2 (HVP cap inside CG: finish the direction and line search
without further HVPs, finish the sweep, then apply the certificate), B3.3 (`maxit` is the NCG sweep counter; the
report and code comments must state that under `w_solver = "ncg"` it is the amendment-3 sweep budget, not the old
PG iteration count).

**Q0 discipline.** A Q0 failure confined to the predicted roundoff-limited class (b) is still a Q0 failure under the
frozen criterion. It is not waived because the derivative looks right. Counts go to the PI first; the PI then
decides whether the test itself was malformed.

**Added Stage B verification invariant (diagnostic, not a gate).** In the qualification harness, at CG exit
(every solve in Q1, a seeded sample in Q-T), compute ||Hess d + g|| directly and compare with the recursive
residual; report max relative disagreement and where it occurs. This checks finite-precision drift of the
recursive residual on the ill-conditioned system. It binds nothing unless the frozen amendment already requires it.

**Sequence.** Stage A, then coordinator inspection of Q0, then Stage B dispatch. No Stage B work before that.

## PI ruling 2026-10-05, Q0 result

Q0 failed at 045f593 on the order clause only. PI ruling: defect of the test, not of the HVP. Test-only Q0
amendment recorded in the frozen amendment, section 10 (half-decade grid 1e-1 to 1e-8, order by eps ratio, two
consecutive p_k >= 1.5, descriptive roundoff classification). HVP code untouched. Original result retained;
amended result appended. Stage B blocked until amended Q0 passes on all 5,400 directions.

## PI ruling 2026-10-05, Stage B interpretation (recorded before any Stage B result exists)

Recorded and committed while the Stage B subagent was running, before it returned any result.

- **Residual drift, tier (c), binding.** At a CG exit whose reason is forcing, the recursive residual met
  min(0.5, sqrt(||g||)) ||g|| but the directly computed ||Hess d + g|| did not. Any such event in an
  instrumented solve (every Q1 solve, the seeded Q-T sample) is a qualification failure.
- **All other drift is descriptive only.** There is no negligible/appreciable threshold, now or later, and it has
  no effect on qualification. Report the raw quantities: maximum relative disagreement, quantiles (median, IQR and
  tails), broken down by CG exit reason, plus a statement of whether any tier-(c) event occurred. This supersedes
  the (a)/(b) split in the instruction sent to the subagent.
- **Q2 settings.** If the frozen documents do not uniquely fix the dataset, seed and hyperparameters for "the
  p = 150 control settings", Q2 is not run, not even exploratory. The subagent returns the candidate definitions
  with file:line, and the PI selects one prospectively.
- **Stop on failure.** A failure of any Q-T test, Q1 below 4/4, a tier-(c) event, or Q2 failing any frozen
  criterion is recorded and the run stops. No tolerance changes, CG changes, budget increases, or reruns under
  altered settings (amendment section 7).
- **Q3 not run** in Stage B. It is descriptive and not part of amendment-3 qualification.

Qualification outcome: amendment 3 qualifies if Q-T passes, Q1 is 4/4, no tier-(c) event occurs, and Q2 meets
the frozen equivalence criteria.
