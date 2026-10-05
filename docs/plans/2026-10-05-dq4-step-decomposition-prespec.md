# D-Q4 pre-specification: decomposition of the sigma-losing W steps (drafted 2026-10-05)

Descriptive and non-gating (PI rulings 2026-10-05, fifth and sixth). Its purpose is to say what the estimator
redesign has to change, not to rescue the current candidate. The estimator, the Amendment-3 solver, `palm_lib.R`
and every gate element are unchanged. Simulated fixture only: `sim_data(700, 3000, 218)`, hp from
`tests/out/superseded/pre-amendment3/T7_fit.rds`. No real data. No new optimization. No Amendment 4.

**Status: drafted, NOT frozen.** Two open items (Q-A, Q-B) need PI confirmation. Nothing is replayed until they are
settled and this file is committed as frozen.

Background: `docs/plans/2026-10-05-dq1-outliers-and-score-collapse.md`, section 6 (D-Q3 results, 7b15de6) and
section 7 (PI classification of D-Q3a).

## 1. Replay (PI ruling, sixth)

Instrumentation-only, gated exactly like D-Q1 (`tests/DQ1_block_replay.R`).

- **Runs.** Starts 1, 5, 16 from the exact D-Q3a warm starts (the D-Q1 end points in
  `tests/out/DQ1_block_replay.rds`): `palm_fit_one`, w_solver ncg, tol 1e-5, maxit 19,000, `hvp_cap` 1e6,
  `PALM_DEFS` unchanged. Same solver, budgets, stopping rules and arithmetic order; no changed settings.
- **Wrapper.** `w_step_ncg` only. It calls the original unchanged and passes its result through untouched. The
  decomposition is computed online from (W_in, ws$W, H, b) by pure, read-only, RNG-free code, and only component
  summaries are saved. No W snapshots.
- **Reproduction gate (binding before any new quantity is read).** Per start, the replay must reproduce the saved
  D-Q3a run exactly (`identical()`, or max abs diff 0):
  - every column of `f$trace`, which includes F (Fstart, FH, FW, Fb), r_W, CG iterations, W exit (accepted step),
    face changes and negative curvature;
  - every column of the per-sweep `rec` (sigma, beta, b, max |cor|, cos(w, u));
  - the number of sweeps, the terminal `line_search_failed` sweep, status and reason;
  - end-point W, H, b, F.
  `palm_lib.R` is confirmed unmodified by `git diff --quiet`. If any start fails, D-Q4 stops for that start: the
  first differing sweep and column are reported and its decomposition is not interpreted.

## 2. Primary decomposition: centered score space (PI ruling, sixth)

For factor j: z_j = C X' w_j, Delta z_j = C X' Delta w_j, sigma~_j = ||z_j|| / sqrt(n - 1)
(sigma_j = sqrt(sigma~_j^2 + eps^2); the eps gap is reported). Write

  Delta z_j = a_j z_j + e_j,   a_j = <Delta z_j, z_j> / ||z_j||^2,   e_j orthogonal to z_j.

a_j z_j is the shrink part; e_j is the rotation part (change of score direction).

**Conflict with the interpretation, found while specifying (Q-A).** Because e_j is orthogonal to z_j,

  ||z_j + Delta z_j||^2 = (1 + a_j)^2 ||z_j||^2 + ||e_j||^2,

so, exactly,

  Delta log sigma~_j = log|1 + a_j|  +  (1/2) log(1 + ||e_j||^2 / ((1 + a_j)^2 ||z_j||^2)).
                       (shrink)          (rotation, always >= 0)

The rotation term can only raise sigma. In this decomposition every unit of sigma loss is assigned to shrink by
construction, so the shrink share of the decline is at least 100% whatever the data say, and "mostly rotation" can
never be observed. As a split of *why sigma falls*, the score-space form is a tautology: sigma falling is the score
shrinking.

What score space does measure well is the size of the concurrent direction change, which is what the
"mostly rotation" reading is about (the optimizer finding new score directions while their variance drops). So it
stays, as the primary description of the score trajectory, with these summaries:
- per step: a_j, ||e_j|| / ||z_j||, and the angle between z_j,in and z_j,out;
- per start: cumulative score angle (sum of per-step angles) and net angle (z_j at sweep 1,001 vs the end point),
  set against total Delta log sigma~_j.

**Proposed attribution of the sigma loss (Q-A, coordinator proposal).** Use the loading-side identity, which can
assign sigma loss to either part. X has within-sample ranks in its columns, so X'1 is constant and z_j depends on
w_j only through P w_j, P = I - 11'/p. Write P w_j = r_j v_j (||v_j|| = 1) and g(v) = ||C X' v|| / sqrt(n - 1),
the score size per unit loading. Then sigma~_j = r_j g(v_j), and exactly

  Delta log sigma~_j = Delta log r_j   +   Delta log g(v_j)
                        (loading shrink)     (rotation of the loading toward lower-variance directions of X)

Here a rotation lowers sigma whenever it moves v_j toward directions in which X varies less across samples, so both
categories are observable. This is the PW-space form the ruling lists as secondary; the proposal is to make it the
attribution and keep score space as the description of direction change.

**Structural fact, stated before running.** The standardized score C X' v_j / g(v_j) and the rho term's vector
P W diag(1/sigma) b = sum_j b_j v_j / g(v_j) (up to eps) depend on W only through the directions v_j. At fixed b,
the radius r_j enters only reconstruction. D-Q3b confirmed this invariance along one path; here it is a check.

## 3. Per-step record (every accepted W step; exit not `line_search_failed` or `numerical_failure`)

For each j = 1..3, with H and b as passed to the W block:
1. Score space: a_j, ||e_j|| / ||z_j||, angle, and the two terms of the score identity.
2. Loading space: Delta log r_j, Delta log g(v_j), and 1 - cor(standardized score column j, in vs out).
3. Term changes by partial steps. W_rad has columns m_j,out 1 + r_j,out v_j,in; W_rot has columns
   m_j,in 1 + r_j,in v_j,out (m_j the column mean; both unit norm; evaluated only, never accepted, so nonnegativity
   is not required). Report Delta nmf, Delta surv, Delta ridge, Delta F from `palm_parts` (H and b fixed; hreg, l1,
   l2 do not depend on W) for in -> rad, in -> rot and in -> out, with the interaction residual.
4. Stratifiers already in the trace: CG iterations, W exit, face changes, negative curvature.

**Check before reading anything else.** Along in -> rad, Delta surv and Delta ridge must be zero up to the eps
effect. If not, report that first.

## 4. Reported (per start)

All accepted W steps are decomposed; nothing is selected after the fact.
- Coefficient-carrying columns J = {j : b_j != 0 at the D-Q3a end point} (starts 1 and 5 {1, 2, 3}, start 16 {1, 2}).
- Totals over all accepted steps, per j in J and pooled over J: total Delta log sigma~; its loading-shrink and
  loading-rotation parts; cumulative and net score angle.
- The same totals stratified afterwards by recorded CG iterations (>= 10 vs < 10, the D-Q3a definition).
- Summed Delta nmf, Delta surv, Delta ridge for the rad and rot partial steps, and the summed interaction residual.
- Fraction of accepted steps with Delta log r < 0 and with Delta log g < 0.

## 5. Pre-stated interpretation (PI ruling, fifth; cutoffs Q-B)

Headline: S = (sum Delta log r_j) / (sum Delta log sigma~_j), pooled over j in J and all accepted steps, per start.
S can fall outside [0, 1] if one part moves sigma up.

- **Mostly shrink.** The normalization by sigma(W) is itself the primary structural problem.
- **Mostly rotation.** The optimizer finds changing score directions while driving their variance down; the
  interaction of standardization with the rho term is more central than simple scale collapse.
- **Both contribute materially.** Treat the parameterization as jointly problematic; do not patch one term.

**Q-B (coordinator proposal).** Mostly shrink: S >= 0.75. Mostly rotation: S <= 0.25. Both: otherwise. Assigned
per start; if the starts disagree, report per start with no pooled category. The CG-stratified S is reported
alongside and does not override the all-steps category.

## 6. Coverage note (for the PI; not an open item unless the PI makes it one)

The D-Q3a warm starts are the sweep-1,000 points, so this replay covers the later part of the decline only. Most of
the log-sigma drop happened in Q1 sweeps 1-1,000 (min coefficient-carrying sigma, start 1: 16.6 -> 1.31 in
sweeps 1-1,000, then 1.31 -> 0.643; start 16: 10.3 -> 1.10, then 1.10 -> 0.98). That is roughly 78% (start 1) and
95% (start 16) of the total log drop before the D-Q3a window opens. If the redesign question is "what drives the
collapse", the first 1,000 sweeps carry most of it. Adding them is the same replay-and-gate pattern, against
`tests/out/Q1_warm.rds` from the committed Q1 starts, which D-Q1 already passed once. Not included unless the PI says
so.

## Outputs

- `tests/DQ4_step_decomposition.R` writes `tests/out/DQ4_step_decomposition.{rds,txt}`.
- Smoke test at 30 sweeps, outputs discarded, before the full run.
- Results committed as observed, with a short section appended to the D-Q1 follow-up record.
- Nothing is pushed.
