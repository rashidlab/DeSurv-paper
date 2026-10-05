# Implementation brief: theta-parameterization estimator (spec revision 4) + small-p smoke tests (2026-10-05)

Authority: `docs/plans/2026-10-05-estimator-redesign-spec.md` DRAFT revision 4, PI ruling 13 (implementation and
smoke tests at p <= 300 allowed before freeze; NO p = 3,000 run, NO timing run, NO qualification, NO real data).
The spec is the source of truth. Where this brief and the spec disagree, or where the spec cannot be implemented as
written, STOP on that item and report the mismatch. Do not improvise. Do not commit (the coordinator commits).

## Files

New folder `analysis/estimator-gate/redesign/`. The candidate folder `analysis/estimator-gate/candidate/` is the
failure record: read it, copy from it, never edit it.

| File | Content |
|---|---|
| `theta_lib.R` | the estimator: data object, objective, gradients, W-block HVP, H/W/theta blocks, residuals, fit, model, scoring |
| `kkt_oracle_theta.R` | independent theta-block KKT evaluator, written WITHOUT theta_lib.R (own scores, own dense Breslow gradient), pattern of `candidate/kkt_oracle.R` |
| `ref_objective_theta.R` | independent evaluator of F (survival via `survival::coxph(..., init = theta / s0, iter.max = 0)`), pattern of `candidate/ref_objective.R` |
| `tests/common.R`, `tests/S*.R` | smoke tests below; outputs in `tests/out/` |
| `REPORT.md` | results table + mismatches |

Base R (+ `survival` in the reference evaluator only). Style of `candidate/palm_lib.R`.

## Mathematics to implement (spec sections 1, 2, 4, 6)

Notation: X p x n ranks, C sample centering, P gene centering, k = 3, W p x k on S_+^k (W >= 0, unit columns),
H k x n >= 0, theta in R^k. s0 = sqrt(sb2), sb2 = mean over genes of the row variance (denominator n - 1), computed
once in the data object from X only. Zc = C X^T W (no sigma anywhere). eta = Zc theta / s0. X2 = ||X||_F^2.

F = (1 - a) [ ||X - W H||^2 / (2 X2) + lambdaH ||H||^2 / (2 X2) ]
  + a [ -(2/ne) l(eta) + lambda (nu ||theta||_1 + (1 - nu)/2 ||theta||^2) + (rho/2) ||P W theta||^2 ]

Gradients of the smooth part Psi (F without a lambda nu ||theta||_1 and the indicators); g = dl/deta from `breslow`:
- grad_H = (1 - a) ( -W^T (X - W H) + lambdaH H ) / X2
- grad_W = (1 - a) ( -(X - W H) H^T / X2 ) + a ( -(2/ne) X (C g) theta^T / s0 + rho (P W theta) theta^T )
- alpha-free theta subproblem: f(theta) = -(2/ne) l(Zc theta / s0) + lambda (1 - nu)/2 ||theta||^2 + (rho/2) theta^T M theta,
  M = (P W)^T (P W); grad = -(2/ne) Zc^T g / s0 + lambda (1 - nu) theta + rho M theta; L1 weight lambda nu.
- W-block Hessian (Euclidean, fixed H, theta) applied to V:
  (1 - a) V H H^T / X2 + a [ (2/ne) X C Hl[C X^T V theta / s0] theta^T / s0 + rho P V theta theta^T ],
  Hl[.] the NEGATIVE log-likelihood Hessian-vector product in eta (reuse `breslow_hvp` from the candidate; check its sign
  convention against a finite difference of g).
- L_W bound (spec section 2): (1 - a) ||H||_2^2 / X2 + a ||theta||^2 ( lambda_1(X C X^T) / s0^2 + rho ).

Blocks (spec section 4; mechanics of the candidate unchanged except as stated):
1. H: one HALS sweep, rows j = 1..k: h_j = max(0, w_j^T R_j / (1 + lambdaH)), R_j = X - sum_{l != j} w_l h_l.
   (Candidate `hals_sweep` used 1 + lambdaH X2/(n k): ruling 12 changes this.)
2. W: the candidate's projected-gradient `w_step` (t = min(2 t_prev, t_max = 2), halving, t_min = 1e-20 failure,
   sufficient decrease, `proj_S` with the argmax rule), with the new gradient and F. NOT Newton-CG: drop
   `w_step_ncg`, `ncg_line_search`, `riem_hess_W`, `w_retract` and the HVP budget logic.
3. theta: the candidate's FISTA `solve_b`, warm-started, with the new smooth function, KKT certificate
   ||KKT||_inf <= 1e-10 by the INDEPENDENT oracle, 1e5-iteration cap = numerical failure, accept only if F does not
   increase (else keep theta, record a rejected step).
Order H, W, theta. No extrapolation. Residuals (section 6): r_W as the candidate; r_H = ||min(H, G_H / a_H)||_F /
max(1, ||H||_F), a_H = (1 - a)(1 + lambdaH) / X2; r_theta = oracle KKT. Converged: r < 1e-5 with no F increase.
F increase: F_new > F_old + 1e-13 max(1, |F_old|), after each block. maxit 5,000 sweeps = work limit
(nonconverged, not failure). Statuses, fit-level rule (lowest-final-F converged start, else flagged nonconverged,
else failed) as in the candidate.

Initialization (spec section 5): exactly the candidate's `palm_draw` stream (W0 = runif(p k, 0, max X), then H0
drawn and discarded), projection, 20 HALS sweeps (with the NEW HALS), theta = 0. 30 starts by default, start s uses
set.seed(S + s - 1). No eps anywhere: remove `palm_eps` and every eps argument.

Invalid input (spec section 2): fail deterministically for n < 2, n_e = 0, nonfinite X/y/d, ||X|| = 0, s0 = 0.

Recording (qualification prespec section 6): every 10th sweep per column sigma_j = sd(Zc_j), r_j = ||P w_j||,
g_j = sigma_j / r_j, theta_j; (W, theta) at sweeps 10, 50, 100, 500, 3750 and the endpoint (a start that stops
earlier stores its endpoint for later checkpoints); every F increase; per fit s0 and lambda_1(X C X^T).
Make recording optional by argument so smoke tests stay light.

Model and scoring (spec section 2, "Prediction vector"): beta = theta / s0, v = W beta, v~ = P v over the training
gene universe, G = `DeSurv::desurv_get_top_genes` on W with ntop (as candidate `palm_model`), score = X_G^T v~_G,
v~ never recentered. Store genes, v~, G, theta, beta, s0.

## Smoke tests (simulated data via `candidate/sim.R`; p <= 300; all deterministic; each prints VERDICT)

| ID | Check | Pass |
|---|---|---|
| S1 | central finite differences of F (smooth part) vs grad_H, grad_W, theta-subproblem grad, at 3 random interior points, p = 80, n = 50 with ties, hp corners alpha {1e-3, 0.5}, rho {1e-2, 1e3} | max rel err <= 1e-6 |
| S2 | theta_lib F vs `ref_objective_theta` F at 20 random points | rel diff <= 1e-10 |
| S3 | identity with the frozen candidate: with eps = 0 and theta = s0 b / sigma(W), the survival term and rho term of theta_lib equal those of `candidate/ref_objective.R` (`ref_F` with eps = 0) | rel diff <= 1e-12 |
| S4 | HALS: after one sweep, each row h_j satisfies its KKT (min(h_j, grad_hj) = 0 to 1e-12 relative); F non-increasing over 50 sweeps at fixed W | as stated |
| S5 | theta solve: oracle KKT <= 1e-10 after `solve_theta` at 10 random (W, hp); oracle and theta_lib gradients agree to 1e-10 | as stated |
| S6 | W-block HVP vs finite difference of grad_W at 5 points; dense Hessian (p = 30, k = 3, by HVP on the basis) symmetric to 1e-8 and lambda_max <= L_W (1 + 1e-8) | as stated |
| S7 | invalid inputs each error with a clear message (the six cases) | all error |
| S8 | full fits, p = 80, n = 50: corner grid alpha {1e-4, 0.3, 0.9} x lambda {1e-3, 1, 1e3} x rho {1e-2, 1e3} x nu {0, 0.95}, 3 starts, maxit 5,000 | zero F increases, zero numerical failures; convergence counts REPORTED only |
| S9 | end to end, p = 300, n = 120 (train) + 200 (test) from ONE `sim_data` draw split by column, centre hp (alpha 0.05, lambda 1, nu 0.5, rho 1, ntop 50), 5 starts | runs; statuses reported; held-out C reported (descriptive only); rerun with the same seed is bit-identical |
| S10 | `palm_draw`-equivalent W0 for (seed, s) equals the candidate's | identical |

Smoke-test outputs are mathematical checks only (ruling 13). A failing S1-S7 or S10 is an implementation defect:
fix it and rerun. S8 failures (F increase or numerical failure) are reported with the configuration; do not change
the algorithm, tolerances or caps to make them pass. Convergence rates in S8/S9 are reported, never acted on.

## Breakage inventory (what changes vs the candidate, all must be accounted for in REPORT.md)

palm_data (adds s0, s0 = 0 error), palm_scores (no sigma/eps), palm_parts / palm_F / palm_Psi (theta, new hreg),
palm_beta (beta = theta/s0), palm_grad_H (new hreg), grad_W_terms / palm_grad_W (no sigma chain rule),
bsub_* / palm_grad_b (theta subproblem, M without sigma), kkt_b (unchanged form), hvp_W_terms / palm_hvp_W (new),
w_face / face_tangent / proj_S (unchanged), riem_hess_W / w_retract / ncg_line_search / w_step_ncg (dropped),
resid_original / palm_resid (a_H), hals_sweep (denominator), w_step (new F/gradient), solve_b -> solve_theta
(oracle), palm_draw (unchanged), palm_init_H (new HALS), palm_eps (dropped), palm_fit_one / palm_fit (no eps,
no NCG, recording), pick_start (unchanged), palm_model / palm_score (v = W theta / s0), palm_rank / palm_fold_data
(unchanged). PALM_DEFS: keep only what still applies; list what was removed.
