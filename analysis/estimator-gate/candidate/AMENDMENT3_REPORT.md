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
