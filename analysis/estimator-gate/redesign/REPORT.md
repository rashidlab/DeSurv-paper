# Redesign estimator (spec DRAFT revision 4): implementation and p <= 300 smoke tests (2026-10-05)

Written under PI ruling 13 (implementation and small-p smoke tests before freeze; no p = 3,000, timing,
qualification or real data). Brief: `docs/plans/2026-10-05-redesign-implementation-brief.md`. Implementation by a
subagent; this report and the independent verification are the coordinator's. Smoke outputs are mathematical checks
only and changed nothing in the objective, ranges, thresholds or algorithm.

## Results (counts from `tests/out/S*.txt`)

| Test | Check | Checks run / passed | Max error | Verdict |
|---|---|---|---|---|
| S1 | FD gradients, H, W, theta subproblem | 108 / 108 | 1.4e-9 | PASS |
| S2 | F vs independent `ref_objective_theta` (coxph) | 40 / 40 | 3.6e-16 | PASS |
| S3 | survival + rho terms equal the frozen candidate's under theta = s0 b / sigma(W), eps = 0 | 20 / 20 | 4.7e-15 | PASS |
| S4 | HALS row KKT, F non-increasing over sweeps | 477 / 477 | 4.5e-15 | PASS |
| S5 | theta solve, independent oracle KKT <= 1e-10 | 60 / 60 | 9.1e-11 | PASS |
| S6 | W-block HVP vs FD; dense Hessian symmetric, lambda_max <= L_W | 20 / 20 | 7.7e-11 | PASS |
| S7 | invalid inputs error (n < 2, n_e = 0, nonfinite X, y, d separately, ||X|| = 0, s0 = 0) | 7 / 7 | 0 | PASS |
| S8 | 36 corner configurations x 3 starts, p = 80, n = 50: no F increase, no numerical failure | 216 / 216 | 0 | PASS |
| S9 | end to end p = 300, train 120 / test 200 from one draw, 5 starts; rerun bit-identical | 4 / 4 | 0 | PASS |
| S10 | start draws identical to the candidate's stream | 270 / 270 | 0 | PASS |

**Reported only (never acted on).** S8: 46 of 108 starts converged, 62 reached the 5,000-sweep work limit. The
nonconverged starts sit mainly at alpha = 1e-4 (residual about 3e-5, reconstruction-dominated) and at lambda = 1e-3
with alpha >= 0.3 (residual up to 0.15, weak theta penalty). 3,205 theta steps were rejected. S9: all 5 starts
converged; held-out C 0.851, correlation with the true linear predictor 0.974; |G| = 147 (ntop 50 per factor).

## Independent coordinator verification (scratch scripts, not committed)

- F recomputed from the spec text with my own code (survival via coxph offset): relative difference <= 1.7e-16.
- Central FD of that F vs `palm_grad_theta` (plus the L1 subgradient), `palm_grad_W`, `palm_grad_H` at two
  hyperparameter settings: worst relative error 8.8e-9. The H check uses step 1e-2 (F is quadratic in H; at 1e-6 the
  roundoff floor is about 1e-5 relative because the H gradient is about 1e-5).
- `hals_sweep` equals the closed form max(0, w_j^T R_j / (1 + lambdaH)) to 2.8e-13.
- **Rejected theta steps are roundoff, not defects.** Re-running S8 configuration 18 (alpha 0.9, lambda 1e3,
  rho 1e3, nu 0; 108 rejections) with a logging acceptance rule: every rejected step raised F by 1.6e-16 to
  3.1e-16 relative (1 to 2 ulp). The rule is the spec's "accepted only if F does not increase" (section 4). The F
  increase guard (1e-13 relative) recorded zero increases.

## Breakage inventory

Changed: palm_data (s0, s0 = 0 error, oracle), palm_scores (no sigma or eps), palm_parts / palm_F / palm_Psi
(theta, hreg / X2), palm_beta (theta / s0), palm_grad_H, grad_W_terms / palm_grad_W (no sigma chain rule),
bsub_* -> tsub_* (M without sigma), palm_grad_b -> palm_grad_theta, kkt_b -> kkt_theta (same form), hvp_W_terms /
palm_hvp_W (new), palm_LW (new), resid_H (a_H per ruling 12), hals_sweep (1 + lambdaH), w_step (new F and
gradient), solve_b -> solve_theta (independent oracle), palm_init_H (new HALS), palm_fit_one / palm_fit (no eps, no
NCG, optional recording), palm_model / palm_score (v = W theta / s0). Unchanged: breslow, breslow_hvp, w_face,
face_tangent, proj_S, palm_draw, pick_start, palm_rank, palm_fold_data. Dropped: palm_eps, riem_hess_W, w_retract,
ncg_line_search, w_step_ncg.

## Subagent-reported mismatches, resolved

1. Sweep-3,750 checkpoint "not in the spec": it is required by the qualification prespec section 6 (Q3 window),
   not by Q6. No change.
2. Per-sweep recording set: the prespec section 6 set (sigma_j, r_j, g_j, theta_j) governs. No change.
3. Seven invalid-input cases vs six: the spec's "nonfinite X, y or d" was tested as three cases. No change.
