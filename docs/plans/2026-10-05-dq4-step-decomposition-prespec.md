# D-Q4 pre-specification: decomposition of the sigma-losing W steps (frozen 2026-10-05)

Descriptive and non-gating (PI rulings 2026-10-05, fifth to seventh). Its purpose is to say what the estimator
redesign has to change, not to rescue the current candidate. The estimator, the Amendment-3 solver, `palm_lib.R`
and every gate element are unchanged. Simulated fixture only: `sim_data(700, 3000, 218)`, hp from
`tests/out/superseded/pre-amendment3/T7_fit.rds`. No real data. No new optimization. No Amendment 4. No categorical
cutoffs.

**Status: FROZEN** (PI ruling, seventh, settled Q-A, Q-B and coverage). Committed before anything is replayed.

Background: `docs/plans/2026-10-05-dq1-outliers-and-score-collapse.md`, section 6 (D-Q3 results, 7b15de6) and
section 7 (PI classification of D-Q3a).

## 1. Replay: two segments, each behind an exact-reproduction gate

Instrumentation-only, gated exactly like D-Q1 (`tests/DQ1_block_replay.R`). Starts 1, 5, 16. Same solver, budgets,
stopping rules and arithmetic order; no changed settings.

- **Segment Q1 (sweeps 1-1,000).** From the committed Q1 initial points with the committed Q1 settings, exactly as
  D-Q1 reconstructed them. This segment carries most of the collapse (about 78% of the log-sigma drop for start 1,
  95% for start 16).
- **Segment D-Q3a (sweeps 1,001-end).** From the D-Q3a warm starts (the D-Q1 end points in
  `tests/out/DQ1_block_replay.rds`): `palm_fit_one`, w_solver ncg, tol 1e-5, maxit 19,000, `hvp_cap` 1e6,
  `PALM_DEFS` unchanged.
- **Wrapper.** `w_step_ncg` only. It calls the original unchanged and passes its result through untouched. The
  decomposition is computed online from (W_in, ws$W, H, b) by pure, read-only, RNG-free code; only per-step
  summaries are saved. No W snapshots.
- **Gate, segment Q1 (binding before any new column is read).** Per start, exact equality (`identical()`, or max
  abs diff 0) with `tests/out/Q1_warm.rds`: every column of `f$trace` (this includes Fstart, FH, FW, Fb, r_W,
  rW_start, cg_iter, w_exit, face changes, negative curvature), the number of sweeps, status, reason and end-point
  W, H, b, F; and per-sweep sigma at W-block exit against the `sigma1..3` columns of the D-Q1 `blk` record
  (block "W") in `tests/out/DQ1_block_replay.rds`.
- **Gate, segment D-Q3a.** Per start, exact equality with `tests/out/DQ3a_continuation.rds`: every column of
  `f$trace`, every column of the per-sweep `rec` (sigma, beta, b, max |cor|, cos(w, u)), the number of sweeps, the
  terminal `line_search_failed` sweep, status, reason and end-point W, H, b, F.
- `palm_lib.R` is confirmed unmodified by `git diff --quiet`. If a gate fails for a start and segment, that start
  and segment stop: the first differing sweep and column are reported and its decomposition is not interpreted.

## 2. Primary: exact decomposition of Delta log sigma (PI ruling, seventh, Q-A)

X has within-sample ranks in its columns, so X'1 is constant and the centered score z_j = C X' w_j depends on w_j
only through P w_j, P = I - 11'/p (P W is `sweep(W, 2, colMeans(W))` in `palm_lib.R`). Write P w_j = r_j v_j,
r_j = ||P w_j||, ||v_j|| = 1, and g(v) = ||C X' v|| / sqrt(n - 1). Then sigma~_j = ||z_j|| / sqrt(n - 1) = r_j g(v_j)
and, exactly, for each accepted W step,

  Delta log sigma~_j = Delta log r_j       +   Delta log g(v_j)
                        (loading-norm part)     (orientation-to-X part)

sigma_j = sqrt(sigma~_j^2 + eps^2); the gap between Delta log sigma and Delta log sigma~ is reported. The H and b
blocks do not change sigma, so the accepted W steps carry all of the sigma movement.

Why not score space for this: with Delta z_j = a_j z_j + e_j and e_j orthogonal to z_j,
||z_j + Delta z_j||^2 = (1 + a_j)^2 ||z_j||^2 + ||e_j||^2, so the orthogonal part can only raise sigma and every
unit of sigma loss is assigned to shrink by construction.

**Caution on naming (PI ruling, seventh).** p = 3,000 > n = 218, so C X' has a large null space. A change in v_j
toward it lowers g(v_j) while the standardized score direction can stay essentially unchanged. That is a real
change of loading orientation that behaves scientifically like a rescaling of the same score. The orientation part
is therefore never called "rotation" on its own; section 3 says whether the score direction moved.

**Structural fact, stated before running.** The standardized score C X' v_j / g(v_j) and the rho term's vector
P W diag(1/sigma) b = sum_j b_j v_j / g(v_j) (up to eps) depend on W only through the directions v_j. At fixed b,
the radius r_j enters only reconstruction.

## 3. Secondary: standardized-score direction change

Per accepted step and factor: the cosine and angle between the standardized score columns before and after the
step. Per start and segment: the cumulative angle (sum over steps) and the net angle (first vs last point of the
segment, and sweep 1 vs the D-Q3a end point for the combined trajectory).

## 4. Per-step record (every accepted W step; exit not `line_search_failed` or `numerical_failure`)

For each j = 1..3, with H and b as passed to the W block:
1. Delta log sigma~_j, Delta log r_j, Delta log g(v_j), and the eps gap.
2. Standardized-score cosine and angle (section 3).
3. Term changes by partial steps. W_norm has columns m_j,out 1 + r_j,out v_j,in; W_orient has columns
   m_j,in 1 + r_j,in v_j,out (m_j the column mean; both unit norm; evaluated only, never accepted, so nonnegativity
   is not required). Report Delta nmf, Delta surv, Delta ridge, Delta F from `palm_parts` (H and b fixed; hreg, l1,
   l2 do not depend on W) for in -> norm, in -> orient and in -> out, with the interaction residual
   (out - in) - (norm - in) - (orient - in).
4. Already-recorded stratifiers from the trace: CG iterations, W exit, face changes, negative curvature.

**Check before reading anything else.** Along in -> norm, Delta surv and Delta ridge must be zero up to the eps
effect. If not, this is reported first: the decomposition then does not separate what it claims to separate.

## 5. Reported (per start; segment Q1, segment D-Q3a, and combined)

Every accepted W step is decomposed; nothing is selected after the fact.
- Coefficient-carrying columns J = {j : b_j != 0 at the D-Q3a end point} (starts 1 and 5 {1, 2, 3}, start 16
  {1, 2}). Non-carrying columns are reported too, labeled as such.
- Cumulative sums per j and pooled over J: Delta log sigma~, Delta log r, Delta log g, as continuous quantities.
  No share, no ratio and no cutoff is the headline, since the two parts can have opposite signs.
- The same sums stratified afterwards by recorded CG iterations (>= 10 vs < 10, the D-Q3a definition).
- Cumulative and net standardized-score angles (section 3), and their distribution over steps in each CG stratum.
- Summed Delta nmf, Delta surv, Delta ridge for the norm and orient partial steps, and the summed interaction
  residual.
- Fraction of accepted steps with Delta log r < 0, and with Delta log g < 0.

## 6. Reading (descriptive; PI ruling, seventh, Q-B)

No numerical threshold. Language, used only when one contribution is clearly larger in absolute magnitude:
- **Predominantly norm-driven** (|sum Delta log r| clearly larger): literal loading shrinkage. The normalization by
  sigma(W) is itself the primary structural problem.
- **Predominantly g(v)-driven** (|sum Delta log g| clearly larger), then by the section 3 angles:
  - score angle essentially zero: orientation into a low-variance or null-space direction of X with the score
    direction preserved, a rescaling of the same score;
  - score angle material: the optimizer finds changing score directions while driving their variance down; the
    interaction of standardization with the rho term is more central than simple scale collapse.
- **Neither clearly larger, or opposite signs of comparable size:** reported as a combination, and the
  parameterization is treated as jointly problematic rather than patched one term at a time.

Segments are read separately first. The Q1 segment describes the mechanism that generated the collapse; the D-Q3a
segment describes the late stiff regime. If they differ, both are reported and neither overrides the other.

## Outputs

- `tests/DQ4_step_decomposition.R` writes `tests/out/DQ4_step_decomposition.{rds,txt}`.
- Smoke test at 30 sweeps per segment, outputs discarded, before the full run.
- Results committed as observed, with a short section appended to the D-Q1 follow-up record.
- Nothing is pushed.
