# D-Q1 follow-up: outlier sweeps, b active set, and score collapse (2026-10-05)

Descriptive only. These numbers come from the saved fits (`tests/out/DQ1_block_replay.rds`, `tests/out/Q1_warm.rds`,
`tests/out/superseded/pre-amendment3/T7_fit.rds`) and the fixture `sim_data(700, 3000, 218)`, which has true
beta (2, 1, 0). No new optimization was run. Scripts: `tests/DQ1b_outlier_sweeps.R`, `tests/DQ1c_score_collapse.R`.
Hyperparameters: alpha 0.0291, lambda 1, nu 0.475, rho 3.16, eps 1.0e-5.

## 1. b active/sign set (PI protection 2)

- Over all 1,000 sweeps of starts 1, 5 and 16, the sign pattern of b never changes.
- Starts 1 and 5: all three b_j are nonzero, at about +/-0.423.
- Start 16: b_3 = 0 exactly, with |g_3| / (lambda nu) = 0.62. That is strictly inside the subdifferential, so the
  zero is stable.
- At these end points, the active-set Schur complement over the free b coordinates is well defined.

## 2. The outlier sweeps

| start | sweep | r_W after W | r_W after b | CG iters | t | W face changes | neg. curvature | dW | db |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 935 | 3.2e-1 | 3.3e-3 | 93 | 1/16 (4 halvings) | 0 | 0 | 8.2e-3 | 2.1e-4 |
| 16 | 949 | 1.6e-1 | 1.6e-3 | 77 | 1/16 (4 halvings) | 2 leave | 0 | 3.9e-3 | 3.3e-5 |

- Start 5 has no sweep in its last 100 with r_W after W above 10x the median.
- Both outliers exit CG on the forcing test, with no negative-curvature exit and no b sign change.
- In the last 100 sweeps, 83 to 93 W steps use a single CG iteration and t = 1.
- W face changes come in periodic enter/leave pairs about every 10 sweeps in starts 1 and 16, and do not occur at
  all in start 5, which stalls in the same way. So W bound chatter is not needed for the stall.

The outliers are a long fixed-b Newton step that badly mismatches b. The next b solve repairs most of that mismatch.
This is the coupled mechanism at large step size, not a qualitatively different event.

## 3. Score collapse at every end point (not anticipated by the D-Q1 interpretation)

Standardized score columns Zs = Zc / sigma(W) at the end points:

| class | F | score columns | b | sd(eta) | T7 starts (of 30) |
|---|---|---|---|---|---|
| three collinear | 0.73108-0.73110 | all \|cor\| = 1.0000 | about +/-0.423 each | 1.27 | 1, 5, 13, 15, 21, 22, 30 (7) |
| two collinear | 0.73298-0.73300 | one pair \|cor\| = 1.0000 | two at +/-0.589, one 0 | 1.18 | 2, 12, 16, 17, 19, 20, 23-28 (12) |
| single | 0.73741 | max \|cor\| 0.16-0.65 | one at +/-0.982, two 0 | 0.98 | 3, 4, 6-11, 14, 18, 29 (11) |

- Every start's Cox predictor is a single score direction. When columns are collinear, b is split equally across
  them, which is the ridge grouping effect.
- Lower F goes with more collinear columns. The lowest-F class (three collinear) is the most degenerate.
- Q1 start 9 (the third saved Q1 fit) is the one start that converged (9 sweeps). It is in the single class.
- The three stalled Q1 starts (1, 5, 16) are in the collinear classes. Their W columns have moved from cosine about
  0.8 (T7) to 0.993-0.999, and their H rows have |cor| 1.000.

## 4. sigma keeps shrinking along the stalled trajectories

For the coefficient-carrying factors, sigma_j falls throughout the 1,000 sweeps and beta = b / sigma rises by the
same factor:

| start | sigma at sweep 1 | 100 | 1000 | log-log slope, sweeps 100-1000 | F drop, last 100 sweeps |
|---|---|---|---|---|---|
| 1 | 19.4 / 19.3 / 16.6 | 4.6 / 4.3 / 3.0 | 2.3 / 2.2 / 1.3 | -0.33 | 4.5e-7 |
| 5 | 19.1 / 13.0 / 19.2 | 5.7 / 3.6 / 5.8 | 2.8 / 1.5 / 2.8 | -0.40 | 4.7e-9 |
| 16 | 10.9 / 10.3 / (5.8) | 1.6 / 1.6 / (6.4) | 1.1 / 1.1 / (6.6) | -0.16 | 8.2e-8 |

Parentheses mark the factor with b_3 = 0. Its sigma does not shrink.

- The decay follows a power law with no plateau. The current sigma values are 5 orders of magnitude above the
  eps floor (1e-5).
- 59-71% of the total log-sigma drop comes from 10 sweeps, and those sweeps ran 35-100 CG iterations, against a
  median of 1 CG iteration elsewhere. The long-CG sweeps, including the outliers in section 2, are the steps that
  move the iterate along this path.
- All W columns have cosine 0.88-0.91 with the uniform vector. For a uniform w, X'w is constant across samples,
  because each column of the fixture's rank-transformed X is a permutation of 1..p. So the score variance comes
  from a deviation of W around a near-uniform profile, and that deviation is shrinking.

## 5. What this does and does not show

Shown:
- The stalled starts are not sitting near a stationary point and zig-zagging. F is still decreasing.
- The productive steps drive sigma toward zero along a collapsed-factor path.
- The W-b zig-zag is the local symptom of moving along that path.

Not shown:
- Whether F reaches a minimum at a finite sigma above eps, or keeps decreasing until the eps floor.
- Whether the same thing happens on real data. The frozen spec also uses within-sample ranks there (spec lines 25, 56),
  so X'u is constant across samples on real data as well; the structural precondition is not specific to the fixture.
- Whether scale invariance (standardization by sigma(W), with the rho term built from P W diag(1/sigma) b) makes the
  supervised part of the objective nearly free to move away from reconstruction. That is a hypothesis, not a result.

Consequence for Amendment 4:
- If the infimum is approached as sigma goes to eps, Option A would speed the solver down this path instead of
  reaching a usable KKT point.
- The frozen r < 1e-5 certificate might then be met only at a degenerate end point, or never.
- That would be a statistical-design question (estimator specification), not a solver question.
- Recommend a pre-registered descriptive check before any Amendment 4 design: does sigma plateau? A suggested shape
  is in the handoff reply. It is not frozen and not run.

## 6. D-Q3 results (prespec 50a79be; scripts aa3515f, d01d2b0)

### D-Q3b, the one-dimensional collapse path (`tests/out/DQ3b_collapse_path.txt`)

- Premise holds: sd of X'u over samples = 1.3e-10.
- The stated prediction holds. Along the path, cor(eta(s), eta(1)) = 1.00000000. The survival and ridge terms at
  fixed b change by at most 2e-9 for s >= 0.03, and drift only for s <= 1e-3, where sigma approaches eps.
- In all three starts and both variants, F **increases** monotonically as s falls below 1. The rise comes from the
  reconstruction term: +1.1e-4 to +1.2e-4 at s = 0.001 (starts 1 and 5) and +8e-5 (start 16).
- Pre-stated reading: **not demonstrated along this path.** Shrinking the deviation from uniform with the score held
  fixed is not a free descent direction. The supervised terms are exactly invariant along it, and reconstruction
  pays.

### D-Q3a, the extended continuation (`tests/out/DQ3a_continuation.txt`; the report comes from `DQ3a_report.R`)

Gate PASS for all three starts. **No start reached 20,000 sweeps.** Each stopped on `line_search_failed`, which the
prespec did not anticipate.

| start | sweeps | min coefficient-carrying sigma, 1000 -> end | log-log slope | max \|cor\| | F, 1000 -> end | median r_W, 1000 -> end |
|---|---|---|---|---|---|---|
| 1 | 4,604 | 1.31 -> 0.643 | -0.48 | 1.00000 | 0.7309114204 -> 0.7309086246 | 3.1e-3 -> 2.0e-3 |
| 5 | 2,859 | 1.51 -> 0.859 | -0.50 | 1.00000 | 0.7309089263 -> 0.7309073543 | 2.3e-3 -> 1.2e-3 |
| 16 | 1,581 | 1.10 -> 0.980 | -0.31 | 1.00000 | 0.7328099851 -> 0.7328098467 | 1.5e-3 -> 1.3e-3 |

- beta rises in proportion as sigma falls. The b sign pattern never changes.
- Essentially all of the log-sigma drop (share 1.00) comes from the 4-33 sweeps that ran 10 or more CG iterations.
- Terminal sweeps:
  - The Newton-CG directional slope falls to between -1e-17 and -2e-20, and the per-sweep F change reaches 1e-15 to
    1e-12, around the double-precision resolution of F (about 0.73).
  - The backtracking search then fails after 39-49 halvings, while r_W stays at 1.2e-3 to 1.9e-3.
  - The remaining residual lies entirely on free entries of ordinary size (W about 0.019). The zero-bound part is 0.

### Post-hoc, outside the prespec: curvature along the stuck residual (`tests/DQ3c_curvature_posthoc.R`)

Rayleigh quotient of the W-block Hessian along the face-tangent residual g:

| start | T7 (sigma about 20) | D-Q1, sweep 1000 | D-Q3a end |
|---|---|---|---|
| 1 | 1.8 | 2.5e5 | 4.3e6 |
| 5 | 3.8 | 1.4e5 | 1.3e6 |
| 16 | 0.38 | 1.5e6 | 2.4e6 |

- At the D-Q1 and D-Q3a end points, almost all of this curvature comes from the rho ridge term
  (rho sigma_bar^2/2) ||P W diag(1/sigma) b||^2. Reconstruction contributes about 2e-7, and the Cox term is negative.
- In start 1, curvature rises about 17x while sigma halves, consistent with sigma^-4 scaling.
- The largest F decrease available along g, ||g||^2 / (2 g'Hg/g'g), falls from about 1e-8 (T7) to 3e-13 to 5e-13 at
  the end points. That is close to the resolution of F.

### Reading against the pre-stated categories (for the PI; not forced)

- Over the whole window the solver could run, every "collapse continues" feature is present: sigma moves
  systematically downward and the slope steepens, the score columns stay collinear, and F decreases while beta
  compensates.
- "Interior stabilization" is not observed: sigma never stabilizes.
- The window ended early on a numerical-resolution failure that the prespec did not anticipate. The post-hoc
  curvature check suggests a cause: along the collapse trajectory, the rho term makes the residual direction
  increasingly stiff, so the frozen certificate (r_W < 1e-5) is not reachable in double precision from these end
  points.
- D-Q3b shows the mechanism is not a free shrink toward uniform. The descent that lowers sigma runs along a
  different direction, which these diagnostics have not identified.
