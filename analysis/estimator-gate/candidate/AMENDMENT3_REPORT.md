# Amendment 3 implementation report

Authority: `docs/plans/2026-10-05-amendment3-newton-cg-w-solver.md` (frozen), sections 3 and 6 Q0. Brief:
`docs/plans/2026-10-05-amendment3-implementation-brief.md`, Stage A. Branch `recompute/consensus-desurv`, no commits.

## Stage A

### Verdict

**Q0 FAIL under the literal criterion.** Euclidean 1,120 / 2,700 directions pass; Riemannian 2,624 / 2,700 pass.
Every failure is class (c). The accuracy clause (min e <= 1e-6) holds on all 5,400 directions (worst 1.5e-7
Euclidean, 9.2e-8 Riemannian). Every failure is on the order clause: the central difference reaches its roundoff
floor after one clean second-order refinement (or none), so two consecutive refinements of order >= 1.5 are not
available on the frozen grid. Per the failure rule (amendment section 7) nothing was tuned; Stage B not started.

### Derivation notes (A1 checked against the code)

- A1 matches `grad_W_terms` term for term. Checks made before coding:
  - Cox: `grad_W_terms` uses q = X r without the sample-centering C. The Breslow score sums to zero, so C r = r,
    and the derivative of the code's own gradient is X dr exactly. Hll[v] also sums to zero.
  - Ridge: u = Wc beta already has zero gene-mean, so P u = u and dq = rho sb2 du as written.
  - sigma chain: d(c_j beta_j / ((n-1) sigma_j^2)) gives the stated ds_j, including -2 s_j dsigma_j / sigma_j.
  - Breslow HVP: differentiating g_l = d_l - e_l h_l gives the stated Hll (dh_l = -sum d_i T_i / S_i^2).
- Fast path: same shift and same condition as `breslow()` (min S over events > 1e-200, ll and g finite).
  T_i / S_i^2 is computed as (T_i / S_i) / S_i, because S_i^2 can underflow while S_i > 1e-200 does not. This is
  the same quantity and leaves the path condition unchanged.
- Log-sum-exp path: suffix log-sums of e and of e v^+ and e v^- (T split by sign), then prefix log-sums of
  d_i T_i^{+/-} / S_i^2 and d_i / S_i, with the same `gs` / `ge` tie groups. A helper `lse2` handles -Inf
  operands (v^+ or v^- can be all zero on a suffix).
- Face and tangent projection: `face_tangent` masks with F, then subtracts w_j <w_j, (V*F)_j>. The brief's form
  assumes unit columns, which holds at N(W) where the Riemannian operator is used. The identity
  ||Pface(G_W)||_F = resid_W holds algebraically at any W, because w_j is zero off F.

### Code (A2): additions to `palm_lib.R` only

New section "W-block Hessian-vector products", inserted after `kkt_b` (the end of the gradients block). Names are
as suggested, plus one helper:
- `lse2(a, b)`: log-sum-exp of two terms with -Inf handling. This is the one extra helper.
- `breslow_hvp(eta, D, v, force_lse = FALSE)`
- `hvp_W_terms(W, H, b, D, hp, V)` returns list(recon, cox, ridge), unweighted.
- `palm_hvp_W(W, H, b, D, hp, V)`
- `w_face(W, GW)`, `face_tangent(V, W, free)`
- `riem_hess_W(W, H, b, D, hp, V, free, GW)`

No finite differences in any of them.

### Test (A3): `tests/Q0_hvp.R`

Outputs: `tests/out/Q0_hvp.txt` and `tests/out/Q0_hvp.rds`.
- Problems: the 270 T1 Part A points. The test reuses T1's seed (`1000 + 100 gi + pt`) and T1's draw order for
  W, H and b, and stops drawing before T1's dW, dH and db, which come after W, H, b and eps. Data come from
  `sim_data(11, 80, 50)`.
- Directions: `set.seed(700000 + 100 gi + pt)`, then 10 raw p x k normal matrices. The Euclidean test uses
  V = raw * F(W) at the T1 point, normalized to unit Frobenius norm. The Riemannian test uses
  V = face_tangent(raw, N(W), F(N(W))), normalized, and compares against Pface(fd) - V diag(mu).

Two readings of the brief are mine. Both are open to coordinator correction:
1. For the ambient Euclidean test, "free face" is read as the coordinate mask V * F. W is off the sphere there,
   so a tangent projection is not defined.
2. The same 10 raw draws feed both tests.

### Results

Binding Q0 counts:

| | evaluated | pass | fail | class a | class b | class c |
|---|---|---|---|---|---|---|
| Euclidean | 2,700 | 1,120 | 1,580 | 0 | 0 | 1,580 |
| Riemannian | 2,700 | 2,624 | 76 | 0 | 0 | 76 |

Failures by alpha x rho, Euclidean:

| alpha \ rho | 0.01 | 1 | 1000 |
|---|---|---|---|
| 1e-4 | 164 | 177 | 190 |
| 0.3 | 174 | 192 | 154 |
| 0.9 | 126 | 220 | 183 |

Failures by alpha x rho, Riemannian:

| alpha \ rho | 0.01 | 1 | 1000 |
|---|---|---|---|
| 1e-4 | 14 | 25 | 0 |
| 0.3 | 19 | 1 | 1 |
| 0.9 | 6 | 10 | 0 |

Class (c) failures broken down by the longest run of consecutive refinements with order >= 1.5:

| run length | Euclidean | Riemannian | what it means |
|---|---|---|---|
| 1 | 1,479 | 76 | one clean order-2 refinement, then the roundoff floor |
| 0 | 101 | 0 | no clean refinement; the HVP error already sits near the floor at eps = 1e-2 |

- Among failing directions, the first-ratio order has median 2.000.
- Failing directions are those with small truncation error. Median e(1e-2) is 2.2e-6 for failing and 5.9e-5 for
  passing Euclidean directions. All 76 Riemannian failures are at pt = 5, the large-eps-floor point, where the
  survival terms are nearly quadratic.

Median e(eps) by alpha shows clean second-order decay to a roundoff floor near eps = 1e-4, with error growing
roughly as 1/eps below that:

| test, alpha | 1e-2 | 1e-3 | 1e-4 | 1e-5 | 1e-6 | 1e-7 | 1e-8 |
|---|---|---|---|---|---|---|---|
| Euclidean, 1e-4 | 1.27e-6 | 1.27e-8 | 1.21e-9 | 5.70e-9 | 5.67e-8 | 4.99e-7 | 5.60e-6 |
| Euclidean, 0.3 | 1.16e-5 | 1.16e-7 | 1.24e-8 | 6.72e-8 | 6.97e-7 | 7.23e-6 | 7.15e-5 |
| Riemannian, 0.3 | 2.79e-4 | 2.79e-6 | 2.83e-8 | 9.27e-9 | 9.36e-8 | 9.15e-7 | 9.65e-6 |

**Conflict to raise: the class (b) definition is unreachable.** Class (b) requires "e <= 1e-9 at every eps". On
these problems central-difference roundoff grows to about 1e-6 to 1e-5 at eps = 1e-8, so a roundoff-limited
direction can never meet that condition, including the alpha = 1e-4 case the brief predicted. Every
roundoff-limited failure therefore falls in (c), with the run-length breakdown above as the description. I did not
reclassify anything.

Assessment, offered for the coordinator and not as a criterion: the evidence points to a correct HVP failing an
order clause that the frozen grid cannot resolve.
- The accuracy clause holds everywhere.
- The truncation ratio is exactly 2.0 where resolved.
- Symmetry error is 2.5e-13.
- Every diagnostic below agrees.

Under amendment section 7 this is a Q0 failure. It needs a PI ruling (diagnose and amend). The criterion has not
been changed.

### Diagnostics (report only)

- **Per-term unweighted checks**, same grid: worst min e is 3.9e-12 for recon, 1.5e-7 for cox and 1.5e-7 for
  ridge. Under the literal order clause, recon passes 0 / 2,700, as predicted: the term is linear in W, so its
  difference quotient has no truncation error. Cox passes 1,160 / 2,700 and ridge 1,119 / 2,700, with the same
  roundoff-floor pattern as the binding test.
- **Symmetry** |<U, Hess V> - <V, Hess U>| / (||U|| ||Hess V||): 270 pairs, max 2.5e-13, median 8.4e-15.
- **g identity**, max |‖Pface(G_W)‖_F - resid_W|: 0 at the 270 ambient T1 points and 0 at N(W).
- **`breslow_hvp`, fast path vs forced log-sum-exp**: 270 points, max relative difference 2.4e-14.
- **Fallback points.** The single-shift min S is 0 at both, so the fast path is refused in `breslow()` and in
  `breslow_hvp()`.
  - Point 1 (my first design): a strictly decreasing sorted predictor over a contrast of 900. ||Hll v|| is only
    2.4e-8 there (each risk set is dominated by its first member), so FD roundoff of the O(1) gradient swamps the
    relative error. Absolute errors stay between 1e-11 and 5e-8. Kept in the output, reported in absolute terms.
  - Point 2, added: one early outlier at 900 with the rest N(0,1), so ||Hll v|| = 0.68. e(eps) runs 6.7e-7,
    6.7e-9, 6.9e-11, then the floor. The point passes the Q0-style criterion.
- **Boundary-face point**: `sim_data(11, 150, 50)` with W = proj_S(normal). W has 228 zeros, of which 133 have
  G_W < 0, so they are free. Face size 355; g identity 0. Euclidean 10 / 10 and Riemannian 10 / 10 pass the Q0
  criterion (worst min e 1.4e-8 and 1.7e-8).

### Fixes

No bug was found by any check, and no implementation line changed after the first run. One change to test code
only: the first fallback point was uninformative for relative error (vanishing HVP), so I added point 2 and
reported point 1 in absolute terms.

### A4 verification by counting

1. `git diff --numstat palm_lib.R`: `72 0`, so 72 added lines and 0 deleted.
2. Q0: 2,700 Euclidean evaluated, 1,120 pass. 2,700 Riemannian evaluated, 2,624 pass. Classes: Euclidean a 0,
   b 0, c 1,580; Riemannian a 0, b 0, c 76.
3. `grep -n 'eps \* V\|h \* V\|/ (2 \*' palm_lib.R` hits only lines 73 (`nmf`) and 74 (`hreg`). Both existed
   before; neither is a finite difference.
4. Regression: `tests/T1_gradients.R` `VERDICT: PASS`, `tests/T2_objective.R` `VERDICT: PASS`. Rerunning them
   rewrote `tests/out/T1_gradients.txt` and `tests/out/T2_objective.txt`; the timestamp is the expected change.

### Wall time

Q0 test 27.1 s end to end (270 problems, 5,400 binding directions, all diagnostics). Implementation and debugging
took one session.

### Files

- `palm_lib.R`: additions only.
- `tests/Q0_hvp.R`: new.
- `tests/out/Q0_hvp.txt`, `tests/out/Q0_hvp.rds`: new.
- `tests/out/T1_gradients.txt`, `tests/out/T2_objective.txt`: rerun outputs.
- This report.

## Stage A, amended Q0 (amendment 3 section 10; appended, the original result above is retained)

"The original prespecified Q0 failed solely because its epsilon grid did not contain enough pre-roundoff points to
evaluate the requested convergence order. The analytic HVP and algorithm were unchanged; only the finite-difference
verification grid/order calculation was amended before Stage B."

- Amendment text and script committed at 68dc7e9 before the amended run; `palm_lib.R` unchanged since 045f593.
  The default script path still reproduces the original failed run (rerun output identical apart from the timestamp).
- Run: `Q0_AMENDED=1 Rscript tests/Q0_hvp.R` -> `tests/out/Q0a_hvp.{txt,rds}`, wall time 51.2 s.
- **Amended Q0 PASS.** Euclidean 2,700/2,700, Riemannian 2,700/2,700; classes a/b/c all 0. Worst min e 6.9e-8
  (Euclidean), 5.2e-8 (Riemannian).
- Coordinator recomputation from the saved 2,700 x 15 error matrices with independent criterion code: 2,700/2,700
  on both; median empirical order over the first four refinements 2.00, 2.00, 2.00, 2.00 (Euclidean) and
  1.99, 2.00, 2.00, 2.00 (Riemannian).
- Diagnostics unchanged: symmetry max 2.5e-13; g identity exact; boundary-face point 10/10 on both tests.

Stage B is unblocked.

## Stage B

Authority: the frozen amendment, sections 3-7 and 10. Brief: Stage B (B1-B3) and the PI rulings, including the
Stage B rulings sent during this run: the drift tier (c) rule, Q2 not to be run until its settings are pinned, stop
on failure, Q3 out of scope. No commits; nothing under `docs/` edited; simulated data only.

### Verdict

**Stopped on failure (PI ruling 3).** Two stop conditions were met:
- **Q-T T9 FAIL**, 55 of 56: the unchanged sample-permutation invariance check.
- **Q1 FAIL**, 1 of 4: starts 1, 5 and 16 reached the 1,000-sweep cap.

The algorithm, tolerances and budgets were not changed after either result. After the stop, no new test was
started: T7, T8 and T10 were **not run**. T4 and the Q-T drift sample were already running when the Q1 result was
read, and they ran unaltered to completion. **Q2 was not run**, because its settings are not pinned (below).
**Residual drift: 0 tier (c) events** in Q1 (every solve) and in the Q-T sample.

### Code changes (`palm_lib.R`)

`git diff --numstat palm_lib.R`: 130 added lines, 14 deleted. `git diff -U0` hunks touch only these objects:
- **Added (3):**
  - `w_retract(W, d, t)`: the retraction N([W + t d]_+), returning NULL when a column is annihilated.
  - `ncg_line_search(W, d, GW, H, b, D, hp, F0, max_halv = 60L)`: the section 3 search.
  - `w_step_ncg(W, H, b, D, hp, hvp_left, tol = 1e-5, defs = PALM_DEFS, max_hvp = 100L)`: B1.
- **Changed (2), both in the B2 table:**
  - `palm_fit_one`: new arguments `w_solver = c("ncg", "pg")` and `hvp_cap = 25000L`; `maxit = NULL` resolves to
    1,000 (ncg) or 5,000 (pg). The new status reasons are `work_cap_reached` and `line_search_failed`. Per-sweep
    section 5 diagnostics go into `trace` under ncg, and the return value gains `hvp_total`, `w_exits`, `cg_exits`
    and the nc counts.
  - `palm_fit`: passes `w_solver` and `hvp_cap` through, default `ncg`. The summary gains `hvp_total`, `w_newton`,
    `w_noop` and `w_ls_failed`. The comment block above `pick_start` also changed; `pick_start` code did not.
- **Unchanged:** `w_step`, `PALM_DEFS` (including `step_init`), all Stage A HVP functions, and every other function.
  The pg path executes the old code verbatim inside `if (!ncg)`; `trace` keeps its old columns under pg.

Code comments in `palm_fit_one`, `palm_fit`, T4, T7 and T9 state that under `w_solver = "ncg"`, `maxit` is the
amendment-3 sweep budget, not the old PG iteration count (ruling B3.3).

The finite-difference grep (`eps \* V|h \* V|/ (2 \*`) hits only lines 73 and 74, both pre-existing.

Implementation readings, each open to correction:
1. **Trial count.** "60 halvings" is read as trials at t = 2^0 ... 2^-60 (61 trials). The amendment's mutant
   "accepting the last trial after 60 halvings" implies the trial at 2^-60 is evaluated.
2. **Guard threshold.** The W no-op guard uses palm_fit_one's `tol` (default 1e-5), which is "the certificate's
   threshold". The brief writes the literal 1e-5; the two are identical at the default.
3. **Line search as a helper.** The line search is a separate function so that T9 can test column-annihilation
   rejection directly. Annihilation is unreachable with a tangent d in exact arithmetic, because
   <w_j, w_j + t d_j> = 1, so T9 injects a non-tangent d. This is a finding about the frozen rule: the rejection
   branch is a floating-point safeguard only.
4. **Line-search failure.** It ends the start immediately after the W block (b skipped). `r` is the last completed
   sweep's value, as for the existing failure exits.
5. **Armijo slope in the trace.** The trace slope is recomputed in palm_fit_one as <G_W, W_new - W_old> from the
   returned G_W, independently of the line search's own slope. This lets T4's Armijo recount detect B19_retraction.
6. **Nonfinite trials.** A nonfinite trial point or trial F inside the search is a numerical failure ("nonfinite
   quantities anywhere").

### Test changes (B2 table; nothing else)

Before any rerun, `tests/out/T*.txt/.rds/.log`, `tests/out/mutants/` and `tests/out/quick/` were copied to
`tests/out/superseded/pre-amendment3/`. The copy includes `T7_fit.rds` (md5 a3a5ede5...), which Q1 reads.

- **`mutants.R`:**
  - B7 is replaced by `B7_ls60`: `ncg_line_search` returns its last trial as accepted after 60 halvings.
  - B19_cap is replaced by `B19_retraction`: slope `t * sum(GW * d)`.
  - All other mutants are unchanged.
- **`T10_mutation.R`:** the B7 rows (T4, T9) and B19_cap rows (T4, T9) are replaced by the new ids. Row count
  unchanged (10 controls, 28 mutant runs).
- **`T4_monotone.R`:**
  - "max t <= 2" becomes "max accepted t <= 1 and every accepted step satisfies the Armijo recount" (slope < 0,
    F_W <= F_H + 1e-4 slope, F_W < F_H).
  - MAXIT 5000 becomes the sweep budget 1000. QUICK mode (150) is unchanged.
  - Added: reporting of W exits, HVPs and nonconverged reasons.
- **`T7_timing.R`:** MAXIT 5000 becomes 1000; reporting adds reasons, HVPs and W exits. Not run (stopped).
- **`T9_edge_cases.R`:** 54 checks become 56. The four one-for-one replacements are:
  1. The floor test becomes line-search failure on a sparse start, with F mutated to rise off the start. It
     reaches the 60-halving path.
  2. The cap test becomes column-annihilation rejection.
  3. The maxit test becomes the sweep cap (3 sweeps) plus the HVP cap (`hvp_cap = 5`), both giving
     `work_cap_reached`. The HVP case checks that the sweep completes (ruling B3.2).
  4. The step_init assertion now runs alongside a `w_solver = "pg"` fit (accepted t <= 2).

  Two assertions are added: the W no-op at r_W < 1e-5 (H = 0, b = 0 gives G_W = 0), and `nc_first_iter` (Hess
  mutated to -I) giving d = -g exactly with 1 HVP. One observation is added: a dense start under the same F
  mutation also ended on the 60-halving path, never the no-change guard.

### Q-T results

| Test | Passed / total | Verdict | Note |
|---|---|---|---|
| T1 gradients | 1,554 / 1,554 | PASS | unchanged code path |
| T2 objective | 810 / 810 | PASS | |
| T3 b block | 360 / 360 | PASS | |
| T4 monotone | 135 / 135 runs pass every check | PASS | p = 150, ncg; see below |
| T5 projection/draws | 9 / 9 | PASS | |
| T6 HALS | 36 / 36 | PASS | |
| T7 timing | not run | n/a | stopped (ruling 3); also the only 30-start p = 3000 fit in Q-T |
| T8 recovery | not run | n/a | descriptive; reads T7's output |
| T9 edge cases | 55 / 56 | **FAIL** | sample permutation, 40 sweeps: W rel diff 7.3e-3 (criterion 1e-8) |
| T10 mutation | not run | n/a | stopped; mutant kill counts therefore not available |
| T11 fold-local | 7 / 7 | PASS | |
| T12 transport | 7 / 7 | PASS | |

**T4 detail.** 135 runs and 63,621 sweeps (190,863 block steps checked). Wall time 1.33 h on 4 cores, against
1,200 s per start (harness).
- **Monotonicity:** F increases 0 by the counter and 0 by the independent recount.
- **Accepted steps:** 62,967 W steps accepted, max t 1, and all 62,967 satisfy the Armijo recount.
- **Statuses:** 80 converged (59.3%, pre-amendment pg 49.6%), 54 nonconverged, 0 failures, 1 stopped by the time
  budget. The nonconverged runs split into 53 `work_cap_reached` and 1 `line_search_failed`.
- **W and CG events:** 653 W no-op sweeps, 341 CG curvature exits (1 nc_first_iter), 557 HVP-ceiling exits. HVPs per
  run: median 1,171, max 17,739.
- **By corner.** The weak-lambda corners (lambda 1e-3, nu 0) all run to the 1,000-sweep cap: medians r 1.8e-5,
  4.0e-2, 2.7e-5 and 3.7e-4. Every lambda = 1 corner, including the box centre, converges in a median of 19 to 95
  sweeps; the box centre converged 0 of 15 under pg.
- **The one incomplete run** is at corner 3 (alpha 1e-4, rho 1e3, lambda 1e-3). Its b block took 20 million FISTA
  iterations in total, so it was stopped by the time budget, not by the W solver.

**T9 failure diagnosis** (scratch script, not a test; nothing changed afterwards). Same data, W0, hp and 40 sweeps;
X columns permuted.
- **The pg solver is still invariant:** W rel diff 3.9e-14.
- **Under ncg the two runs agree to about 1e-15 in F through sweep 11.** At sweep 12, the entry states differ by
  1.7e-13 (W), 7.6e-14 (b) and 7.1e-14 (G_W).
- **One CG solve amplifies that difference.** The free masks are identical (0 differing coordinates), and so are
  the CG iteration count (15), the forcing exit and t = 0.5. Yet d differs by 7.0e-6 on |d| about 1, an
  amplification of about 5e7 in that single solve.
- **This is not residual drift.** For that solve the direct ||Hess d + g|| equals the recursive residual
  (4.563e-3 vs threshold 4.615e-3).
- **Curvature at that point.** The Hessian's lambda_max is about 7.8 and lambda_min about 3.6e-4.
- **After sweep 12 the paths separate.** Face changes and CG counts diverge, and by sweep 40 W differs by 7.3e-3
  (unconverged, r 8.5e-3).

Reading: plain unpreconditioned CG on an ill-conditioned system loses floating-point reproducibility under data
permutation, which projected gradient does not. I found no departure from the frozen text. Per section 7 this needs
a PI diagnosis and, if anything changes, a new amendment. The criterion was not relaxed.

### Q1 (binding realistic size): FAIL, 1 of 4

Source: the pre-amendment `T7_fit.rds` (superseded copy). The precondition holds for all four starts: 5,000 PG
iterations, nonconverged, fit eps, stored F reproduced to 1e-14. hp alpha 0.0291276, lambda 1, nu 0.475,
rho 3.16228, eps 1.000768e-5; p = 3000, n = 218. Wall time 1.3 min on 4 cores.

| Start | Status (reason) | Sweeps | HVPs | Final r (r_W / r_H / r_b) | F_stalled | F_final | Monotone (counter / recount) |
|---|---|---|---|---|---|---|---|
| 1 | nonconverged (work_cap_reached) | 1,000 | 3,668 | 2.97e-3 (2.97e-3 / 6.8e-11 / 9.5e-11) | 0.731083 | 0.730911 | 0 / 0 |
| 5 | nonconverged (work_cap_reached) | 1,000 | 3,241 | 1.96e-3 (1.96e-3 / 7.7e-11 / 9.8e-11) | 0.731084 | 0.730909 | 0 / 0 |
| 9 | **converged** | 9 | 119 | 5.70e-6 (5.70e-6 / 1.2e-7 / 2.7e-11) | 0.737414 | 0.737263 | 0 / 0 |
| 16 | nonconverged (work_cap_reached) | 1,000 | 3,081 | 1.53e-3 (1.53e-3 / 5.5e-12 / 7.5e-11) | 0.732985 | 0.732810 | 0 / 0 |

CG exit counts and flags per start:

| Start | forcing | hvp_ceiling | hvp_total_cap | curvature | negative_curvature | nc_first_iter | line-search failures | no-op sweeps | max t | halvings | entered / left support |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 990 | 4 | 0 | 6 | 6 | 0 | 0 | 0 | 1 | 130 | 240 / 265 |
| 5 | 993 | 2 | 0 | 5 | 5 | 0 | 0 | 0 | 1 | 83 | 66 / 69 |
| 9 | 9 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 4 | 25 / 7 |
| 16 | 991 | 6 | 0 | 3 | 3 | 0 | 0 | 0 | 1 | 92 | 75 / 76 |

Criteria per start (certificate / monotone including F_final <= F_stalled / budget):
- start 1: F / T / F
- start 5: F / T / F
- start 9: T / T / T
- start 16: F / T / F

Monotonicity holds everywhere. F improved by 1.5e-4 to 1.75e-4 in every start, about 10 times what 15,000 extra
PG iterations achieved. Descriptively, the minimum matched-column cosine to the stalled W is 0.78 to 0.83, and the
Pearson correlation of fitted eta is >= 0.999997.

**Diagnosis of the three failures** (start 1 end state; scratch script; no change made):
1. **Late sweeps take tiny Newton steps.** After about sweep 20, CG exits on forcing after 1 iteration, the step is
   t = 1, and the slope <G_W, dW> is about 4e-11. r_W sits at 2.5e-3 to 4.5e-3 for 980 sweeps. Only about 3.7 HVPs
   per sweep were used, so the 1,000-sweep cap binds long before the 25,000-HVP cap.
2. **The W block works; the coupling undoes it.** At the end state, g lies entirely on W > 0 coordinates, and 84%
   of it is in column 3. g is nearly an eigenvector of Hess: cos(g, Hg) = 0.9993 with Rayleigh quotient 2.50e5,
   all of it from the Euclidean HVP (Weingarten term 3e-4). A test-side central difference confirms Hess g to
   5.5e-6. One W step at fixed (H, b) moves W by 1.5e-9 and cuts r_W from 2.97e-3 to 1.08e-4. The next H and b
   update restores r_W to about 3e-3.
3. **The cause of the stiffness.** The run has moved to a region where the score sds sigma_j(W) are about 10 times
   smaller than at the stalled points (start 1: 2.3, 2.2, 1.3 vs 19.8, 19.8, 19.9; starts 5 and 16 similar). b is
   essentially unchanged. Curvature through sigma_j grows as 1/sigma^2, giving 2.5e5 against the lambda_max of 47
   to 75 found by S1-S5 at the stalled points.

So the binding obstacle is now W-H(-b) coupling along a stiff sigma direction, not CG accuracy. This bears on the
2026-10-04 implicit-regularization finding (exact optimization of Eq. 1 moves away from the stalled basin).
Section 7 applies: diagnose and amend, no tuning.

### Residual-drift diagnostic (PI rulings; tier (c) binding, all other drift raw only)

Implementation is in `tests/drift_hook.R`. At every CG exit, ||Hess d + g|| is computed directly with one extra HVP
that is not counted toward any budget. The quantity reported is |direct - recursive| / direct. nc_first_iter exits
(d = -g) are counted and not compared; there were none.

| Set | Exits | min | q25 | median | q75 | q90 | q99 | max |
|---|---|---|---|---|---|---|---|---|
| Q1, all (every solve, 4 starts) | 3,010 | 0 | 2.7e-15 | 6.4e-15 | 1.4e-14 | 3.8e-14 | 6.2e-10 | 3.45e-7 |
| Q1, forcing | 2,980 | 0 | 2.6e-15 | 6.3e-15 | 1.4e-14 | 3.3e-14 | 5.7e-10 | 3.45e-7 |
| Q1, hvp_ceiling | 12 | 2.7e-12 | 1.4e-10 | 2.7e-10 | 5.3e-10 | 2.5e-9 | 1.0e-8 | 1.1e-8 |
| Q1, curvature | 14 | 1.0e-12 | 5.2e-12 | 8.6e-12 | 1.8e-11 | 3.4e-11 | 1.3e-10 | 1.5e-10 |
| Q-T sample, all | 8,030 | 0 | 1.3e-14 | 5.2e-14 | 1.6e-13 | 3.9e-13 | 2.1e-12 | 4.8e-9 |
| Q-T sample, forcing | 7,960 | 0 | 1.3e-14 | 5.1e-14 | 1.6e-13 | 3.7e-13 | 1.6e-12 | 4.8e-9 |
| Q-T sample, hvp_ceiling | 42 | 2.3e-14 | 7.2e-13 | 3.2e-12 | 8.5e-12 | 2.8e-11 | 4.3e-11 | 4.9e-11 |
| Q-T sample, curvature | 32 | 0 | 2.3e-15 | 1.6e-14 | 4.7e-14 | 9.8e-14 | 1.0e-12 | 1.3e-12 |

- **Maxima.**
  - Q1: 3.45e-7, at start 16 sweep 1 (forcing, 36 CG iterations, recursive 3.643e-6, direct 3.643e-6, threshold
    4.394e-6).
  - Q-T sample: 4.8e-9, at T4 job 32 sweep 6 (forcing, residual 1.2e-13).
- **Margins.** The smallest margin (threshold - direct) / threshold over forcing exits is 1.6e-4 in Q1 and 1.2e-5
  in the Q-T sample.
- **Tier (c) events: 0 in Q1, 0 in the Q-T sample.**
- **The Q-T sample.** `set.seed(20261005)`, 12 of the 135 T4 runs (jobs 1, 2, 3, 8, 32, 65, 68, 85, 89, 90, 97,
  100), rebuilt exactly as in T4 and run to the T4 stopping rule (`tests/drift_QT.R`). Six converged and six reached
  the sweep cap. The capped runs are at T4 corners 1, 2 and 4, with r from 2.2e-5 to 4.0e-2. These hp corners
  differ from Q1's.

### Q2: blocked, settings undefined (not run, per PI ruling 2)

`tests/Q2_p150.R` is written: both solvers, the frozen 30-start fit from the same draws, and all five criteria.
Column matching is exhaustive over the k! permutations, which equals the Hungarian optimum. eta~ is the training eta
of the selected start divided by its sd. The settings block is NULL, so the script stops with an explanatory error.
The frozen documents do not uniquely pin dataset, seed and hyperparameters. Candidate definitions:
1. **The wcond.R control.** `tests/out/coordinator/wcond.R:95-98`: `sim_data(401, 150, 60)`, the T7 box-midpoint hp
   (alpha 0.0291276, lambda 1, nu 0.475, rho 3.16228, lambdaH 1e-6), fit seed 123, but **5 starts**, not 30.
   - The frozen amendment, section 1 item 1 (`docs/plans/2026-10-05-amendment3-newton-cg-w-solver.md:20`), refers
     to "the p = 150 control" via S1-S5.
   - `docs/HANDOFF-2026-10-05-candidate-optimizer.md:75` describes it as "a converged p = 150 control (dataset 401,
     n = 60, same hyperparameters; all 5 control starts converge)".
   - Open: the start count (5 vs the "frozen 30-start fit") and therefore the pooled eps.
2. **The p = 150 settings where the PG solver converges.** `docs/HANDOFF-2026-10-05-candidate-optimizer.md:59`
   reads "solver-equivalence test on the p = 150 settings where the current solver converges". That points to the
   T4 corners that converged:
   - `REPORT.md:28-32` lists (1e-4, 1e-2, 1, 0.95), (1e-4, 1e3, 1, 0.95) and (0.9, 1e3, 1, 0.95), each 15 of 15
     converged, and (1e-4, 1e3, 1e-3, 0) at 14 of 15.
   - Datasets `sim_data(401..403, 150, 60)`, fit seed 7, 5 starts (`tests/T4_monotone.R`).
3. **REPORT.md:95.** "At p = 150 (dataset 401, start 1) the logit midpoint (alpha 0.0291) converges in 417
   iterations." This is the same point as (1) with a single start.

The choice is the PI's, prospectively.

### Deviations and conflicts

- **T7 overlaps the scope instruction.** T7 is part of Q-T under the frozen text, but the brief said not to run
  30-start p = 3000 fits. Moot here (stopped), but it needs a ruling before a Q-T rerun.
- **Q1 source path.** Q1 reads the superseded copy of `T7_fit.rds`, not `tests/out/T7_fit.rds`. A T7 rerun
  overwrites the latter with ncg end points.
- **Pre-existing T10 status.** The pre-amendment T10 was already "FAIL as defined" (24 of 28; four T4 misses whose
  trigger never fires on correct code). A literal Q-T T10 can only pass if those misses disappear.
- **Possible T10 misses.** B7_ls60 is likely to be missed by T4 for the same reason as before: 60 halvings never
  occur on correct code. T9 should kill it, because its line-search test reaches the 60-halving path.
- **The T9 step_init assertion.** It now also checks a pg fit's accepted t <= 2. This merges the reachability check
  of the frozen step into the moved assertion, keeping the count one for one.
- **Duplicate Q1 run.** Q1 ran twice by accident (a launch I believed had died). Both runs are identical; the final
  `Q1_warm.txt` is the 09:25 run.

### Files

- Changed: `palm_lib.R`, `tests/mutants.R`, `tests/T4_monotone.R`, `tests/T7_timing.R`, `tests/T9_edge_cases.R`,
  `tests/T10_mutation.R`.
- New: `tests/Q1_warm.R`, `tests/Q2_p150.R`, `tests/drift_hook.R`, `tests/drift_QT.R`.
- Outputs:
  - `tests/out/Q1_warm.{txt,rds,log}`
  - `tests/out/drift_QT.{txt,rds,log}`
  - `tests/out/T4_monotone.{txt,rds,log}`
  - rerun T1, T2, T3, T5, T6, T9, T11 and T12 `.txt`
- Preserved: `tests/out/superseded/pre-amendment3/`.

## Post-failure diagnostics (D-Q1, D-T9a, D-T9b)

Run 2026-10-05, per "PI approval 2026-10-05: post-failure diagnostics" in
`docs/plans/2026-10-05-amendment3-implementation-brief.md`. Descriptive only: nothing here changes the amendment-3
FAIL recorded at 0c44cca, and no threshold or setting was altered after a result was seen. `palm_lib.R` is
byte-identical to HEAD throughout (`git diff --quiet -- palm_lib.R` confirmed before and after every run).
Instrumentation is by wrapper only, in the global environment, following the pattern in `tests/drift_hook.R`;
every wrapper calls the unmodified original and returns its result unchanged. New files:
`tests/DQ1_block_replay.R`, `tests/DT9_perm.R`; outputs `tests/out/DQ1_block_replay.{txt,rds}`,
`tests/out/DT9_perm.{txt,rds}`. No existing test file was changed, no commits made.

### D-Q1: instrumented block replay of Q1 starts 1, 5, 16

`tests/DQ1_block_replay.R` wraps `hals_sweep`, `w_step_ncg` and `PALM_DEFS$b_accept`/`solve_b` so the state after
each H, W and b block can be read off, while `palm_fit_one` runs its real, unmodified control flow (the wrappers
call the originals and pass their return values through unchanged; the added computations are `palm_resid` and
`palm_scores`, both pure and RNG-free). Inputs exactly as Q1: the four stalled end points at
`tests/out/superseded/pre-amendment3/T7_fit.rds`, `sim_data(700, 3000L, 218L)`, `hp <- fit$hp`, `palm_fit_one`
defaults (`w_solver = "ncg"`, `maxit = 1000`, `tol = 1e-5`, `hvp_cap = 25000`). Wall time 2.2 min (3 starts, 3
cores).

**GATE: PASS, 3 of 3 starts.** The replay's own trace (produced by the unwrapped control flow) is an exact match
to the committed `tests/out/Q1_warm.rds` trace for `Fstart`, `FH`, `FW`, `Fb`, `rW_start` and `rW` over all 1,000
sweeps, and the final status (`nonconverged`, `work_cap_reached`) for starts 1, 5 and 16. (This is expected: the
wrappers add a side channel only and never touch the arguments or return values on the path `palm_fit_one` itself
uses, so the gate is a sanity check on the instrumentation, not a reproducibility test of the optimizer.)

Block pattern, medians [range] over the last 100 sweeps, against the pre-stated hypothesis (r_W ~3e-3 -(W)->
~1e-4 -(b)-> ~3e-3, r_H negligible):

| Start | r_W after H | r_W after W | r_W after b | r_H (max over the three blocks) |
|---|---|---|---|---|
| 1 | 3.124e-3 [2.555e-3, 3.702e-3] | 1.288e-4 [8.033e-5, 3.179e-1]* | 3.124e-3 [2.555e-3, 3.702e-3] | 5.359e-4 |
| 5 | 2.256e-3 [1.977e-3, 2.653e-3] | 8.661e-5 [5.944e-5, 1.224e-4] | 2.251e-3 [1.961e-3, 2.653e-3] | 7.885e-7 |
| 16 | 1.531e-3 [1.310e-3, 1.681e-3] | 4.888e-5 [3.386e-5, 1.613e-1]* | 1.531e-3 [1.310e-3, 1.681e-3] | 8.551e-5 |

\* the upper end of the post-W range at starts 1 and 16 comes from a small number of sweeps with a face change
(`n_enter`/`n_leave` != 0) in the last 100; the median is representative of the steady-state sweep.

**The observed pattern matches the pre-stated hypothesis closely** at all three starts: r_W drops by a median
factor of 24.3x (start 1), 26.0x (start 5) and 31.3x (start 16) on the W block, i.e. to the 1e-4/1e-5
range, and is restored to its pre-H-block value by the b block in the same sweep; r_H stays at or below 1e-3
throughout and is typically 1e-7 to 1e-14 (negligible by any reasonable reading). sigma(W) and beta = b/sigma are
essentially static across the window (medians stable to 3-4 significant figures start to start, consistent with
the 2026-10-04 Stage B diagnosis that the run sits in a stiff-sigma regime). The zig-zag is therefore a real,
repeating W-then-b cycle, not an artifact of the residual definition or of the instrumentation.

### D-T9a: seeded replay of the T9 permutation pair

`tests/DT9_perm.R` reproduces `tests/T9_edge_cases.R` sections 1-3 verbatim (output suppressed via
`capture.output`, which does not touch the RNG stream) to put the RNG in the same state T9 reaches at
`pm <- sample(D$n)`, then builds the pair exactly as T9 does (`W0`, `palm_init_H`, `b = 0`, `hp`, `maxit = 40`,
default `w_solver`).

**GATE: PASS.** Printed to T9's own precision (`%.1e`): replay W 7.3e-03, b 5.4e-05, F 5.3e-07, H 5.2e-03,
identical to the committed `tests/out/T9_edge_cases.txt` line 45.

Comparison (40-iteration, truncated pair; no D1 claim):

| Quantity | Member a | Member b |
|---|---|---|
| F | 1.679994938 | 1.679994054 |
| r / r_W / r_H / r_b | 8.493e-3 / 1.827e-4 / 8.493e-3 / 6.060e-11 | 6.203e-3 / 1.088e-3 / 6.203e-3 / 5.960e-11 |

Hungarian-matched (k = 3, exhaustive over the 3! permutations) W-column cosines: 1.0000, 1.0000, 1.0000. Risk-score
(eta) Pearson correlation, b matched into a's order via `pm`: 1.000000; SD ratio (b/a): 1.0000. Top-5-loading
overlap per matched column: 1.00, 1.00, 1.00 (N = 5 chosen in advance: it is `markers_per_factor` for p = 60 in
`sim.R`'s generator, `max(5, round(150*60/3000))`, i.e. the true marker-set size, not a post hoc choice). Recovery
of the simulated truth (T8_recovery.R's `match_cos` convention, cosine to `W_true` after exhaustive column
matching): a 0.3944, 0.4307, 0.4794; b 0.3945, 0.4308, 0.4795.

So at the parameter level the permutation sensitivity reported in Stage B (W rel diff 7.3e-3) sits on top of a
pair that is, by every descriptive measure above, the same fit: matched cosines 1.0000, risk scores
r 1.000000 with SD ratio 1.0000, recovery equal to within 1e-4 per column (approximately equivalent, not identical). The KKT residuals differ by about 2x (both well above tol), consistent with two
different not-yet-converged iterates on the same basin rather than two different basins.

### D-T9b: the same pair, run to convergence (nonbinding)

Settings frozen in the brief: same D, Dp, W0, H init, b = 0, hp as D-T9a; `w_solver = "ncg"`, `tol = 1e-5`,
amendment-3 per-start budget (1,000 sweeps, 25,000 HVPs) -- i.e. `palm_fit_one` defaults.

| Member | Status | Sweeps | HVPs | Final F | r |
|---|---|---|---|---|---|
| a | converged | 188 | 5,157 | 1.679845161 | 9.731e-6 |
| b | converged | 190 | 5,295 | 1.679845162 | 9.946e-6 |

Both members converged inside the budget (outcome, not a tuning choice). Comparison at convergence: F rel diff
4.543e-11; r/r_W/r_H/r_b a = 9.731e-6/9.731e-6/8.390e-6/6.804e-11, b = 9.946e-6/9.946e-6/9.010e-6/7.777e-11;
Hungarian-matched W-column cosines 1.0000, 1.0000, 1.0000; eta Pearson 1.000000, SD ratio 1.0000; top-5 overlap
1.00/1.00/1.00; recovery a 0.3988, 0.4346, 0.4828, b 0.3988, 0.4346, 0.4828 (matches a to 4 decimal places).
Converged-path W rel diff 7.152e-5, two orders of magnitude below the truncated-path 7.293e-3.

**Interpretation: (i) fits** -- the truncated fits differ (W rel diff 7.3e-3) but the converged fits agree (W rel
diff 7.2e-5, identical cosines, Pearson, overlap and recovery to 4 decimal places). (ii) is false (the converged
fits do not differ materially) and (iii) does not apply (there is no persistent parameter disagreement at
convergence to begin with). Reading, consistent with the Stage B diagnosis: unpreconditioned CG on the
ill-conditioned W-H-b system loses floating-point reproducibility under data permutation during the transient,
early part of the trajectory, and this washes out once both members reach the certificate. The raw-W T9 criterion,
unchanged in this run, is probably the wrong invariant to hold a truncated NCG path to; a future test should
target a model invariant (prediction, recovery) rather than raw W at a fixed, small iteration count.

### Conflicts and deviations

None. Both gates passed on the first run; no settings were changed; no existing test was touched; `palm_lib.R`
is unmodified (confirmed by `git diff --quiet`).
