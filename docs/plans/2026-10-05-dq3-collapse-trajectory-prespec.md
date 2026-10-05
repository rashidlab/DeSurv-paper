# D-Q3 pre-specification: collapse trajectory and one-dimensional collapse path (frozen 2026-10-05)

Descriptive diagnostics. Frozen and committed before anything is run. The estimator, the Amendment-3 solver and every
gate element are unchanged. Simulated fixture only: `sim_data(700, 3000, 218)`, with hp from
`tests/out/superseded/pre-amendment3/T7_fit.rds`. **No real-data diagnostic** (PI ruling 2026-10-05). Amendment 4 is
not drafted until these results are in.

Background: `docs/plans/2026-10-05-dq1-outliers-and-score-collapse.md` (commits 70bb55e, 9964e70).

## D-Q3a: extended continuation

- **Starts and warm start.** Starts 1, 5 and 16, warm-started from their D-Q1 end points (W, H, b after 1,000
  sweeps, `tests/out/DQ1_block_replay.rds`). Before continuing, gate on r_W, F and b of those saved end points: they
  must equal the committed Q1 end points (identical, or max abs diff 0). If they differ, stop and report.
- **Solver.** `palm_fit_one`, unchanged: w_solver ncg, tol 1e-5, frozen certificate, `PALM_DEFS` unchanged.
  maxit = 19,000, which brings each start to 20,000 sweeps in total.
- **Work cap.** `hvp_cap` is raised to 1e6 so that the HVP cap does not end the run. This is a diagnostic setting
  only. It is not an amendment, and it says nothing about the Q1 budget. If a start meets the certificate before
  19,000 sweeps, it stops there and that is reported.
- **Per-sweep record**, taken at b-block exit through the same wrapping pattern as D-Q1:
  - F, and r_W, r_H, r_b
  - sigma_j and beta_j = b_j / sigma_j, j = 1..3
  - b
  - max |cor| among the centered score columns
  - cos(w_j, u), where u is the uniform unit vector
  - from the solver trace: CG iterations, t, face changes, negative curvature
- **Run.** One process per start, 3 in parallel.

**Reported** (no thresholds):
- sigma_j, beta_j, F, max |cor| and median r_W at sweeps 1,000, 2,000, 5,000, 10,000 and 20,000 (or the last sweep)
- the log-log slope of the coefficient-carrying sigma_j over sweeps 1,000-5,000 and 5,000-20,000
- F drop per 1,000 sweeps
- the fraction of the log-sigma drop carried by sweeps with 10 or more CG iterations
- whether the b sign pattern changes

**Pre-stated interpretation** (PI ruling 2026-10-05, verbatim in substance):
- **Collapse continues.** sigma keeps moving systematically downward over the extended run, the score columns stay
  collapsed/collinear, and F keeps decreasing while beta compensates. This is evidence of an
  estimator/parameterization degeneracy, not an outer-solver problem. Do not proceed with Option A.
- **Interior stabilization.** sigma clearly stabilizes away from the numerical floor (eps = 1.0e-5), and the frozen
  KKT residual then contracts. The remaining problem is optimization, and Option A remains plausible.
- **Mixed or ambiguous.** Reported as observed, without forcing it into either category. The 100 x eps figure from
  the earlier draft is not a cutoff.

## D-Q3b: one-dimensional collapse path

**Path.** At each start's D-Q1 end point (W, H, b), let J = {j : b_j != 0} (starts 1 and 5: {1, 2, 3}; start 16:
{1, 2}). For j in J, write w_j = c_j u + d_j, where c_j = u'w_j and d_j is orthogonal to u. Then for each s:

  w_j(s) = (c_j u + s d_j) / ||c_j u + s d_j||;  w_j(s) = w_j for j not in J.

- Since c_j u + s d_j = (1 - s) c_j u + s w_j, the path stays nonnegative and unit-norm for s in (0, 1].
- s > 1 is infeasible wherever w_ij = 0, and every end point has zero entries, so it is not used.
- Grid: s in {1, 0.9, 0.75, 0.5, 0.25, 0.1, 0.03, 0.01, 0.003, 0.001, 3e-4, 1e-4}.

**Why the score direction is held fixed.** Each column of X holds the within-sample ranks, so X'u is constant across
samples. The centered score Zc_j(s) is then exactly proportional to Zc_j(1), and the standardized score is unchanged
up to the eps floor. The rho term uses P W diag(1/sigma); with P u = 0, its numerator and its denominator both scale
with s. **Prediction stated before running:** along the path, surv, l1, l2 and ridge (at fixed b) stay constant up to
the eps effect, so any change in F comes from the reconstruction term (nmf + hreg).

**At each s, two variants:**
- (A) H fixed at the end point H, b re-solved by the frozen `solve_b` (tol 1e-10), warm-started from the previous s.
- (B) H profiled by HALS sweeps until the max relative change in H is below 1e-12 (cap 5,000 sweeps; hitting the
  cap is reported), then b re-solved as in (A).

**Reported per s and variant:**
- F, and each part of F: nmf, hreg, surv, l1, l2, ridge
- sigma_j and beta_j
- b
- cor(eta(s), eta(1)), and the correlation of each standardized score column with its s = 1 value
- r_W at (W(s), H, b(s)), for information only

**Pre-stated interpretation:**
- **Degeneracy demonstrated directly.** In variant (B), F decreases as s goes to 0 while cor(eta(s), eta(1)) stays
  essentially 1. The objective improves by shrinking the factor deviation toward the uniform limit with the
  predictive score unchanged.
- **Not demonstrated along this path.** F increases as s moves below 1 (a minimum along the path at or near s = 1).
  The drift seen in D-Q1 then runs along some other direction, and that is reported as such.
- **Mixed.** Non-monotone F, a disagreement between variants A and B, or a score correlation that drops well below 1
  is reported as observed.
- If the stated prediction (survival parts invariant) fails beyond the eps effect, report it first, because it
  would mean the path construction does not hold the score fixed.

## Outputs

- `tests/DQ3a_continuation.R` writes `tests/out/DQ3a_continuation.{rds,txt}`.
- `tests/DQ3b_collapse_path.R` writes `tests/out/DQ3b_collapse_path.{rds,txt}`.
- Results are committed as observed, together with a short results section appended to the D-Q1 follow-up record.
- Nothing is pushed.
