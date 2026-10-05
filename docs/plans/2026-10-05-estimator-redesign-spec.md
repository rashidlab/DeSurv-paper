# Estimator-redesign specification: rank-3 DeSurv in the beta parameterization (2026-10-05)

**Status: DRAFT, not frozen.** Written under PI ruling 8 (`docs/HANDOFF-2026-10-05-candidate-optimizer.md`, "PI ruling
2026-10-05, eighth"). No implementation, no real data and no compute until the PI freezes this file. Once frozen, every
later output records this file's commit hash. Amendments then follow the rule in
`2026-10-04-estimator-gate-spec.md`: dated, with a reason, and only before any Gate A output is read.

**Relation to the frozen gate spec.** This is a new estimator specification, not Amendment 4. It replaces the
candidate arm (section 4 of `2026-10-04-estimator-gate-spec.md`, including amendments 2 and 3 as they apply to that
section). Sections 1-3 and 5-10 of the gate spec are adopted by reference and are unchanged unless a line below says
otherwise: the question, data and outer/inner resampling (partition seeds 3001-3004, 20 outer folds), comparator arm,
scoring, Gate A, Gate B, the decision table and the safeguards. Amendment 3 stands as a failure record. Its Newton-CG
solver and per-start budget are not carried forward.

**Why the estimator changes (D-Q4, `2026-10-05-dq1-outliers-and-score-collapse.md` section 8).** In the frozen
candidate the survival predictor is eta = Z_c(W) diag(1/sigma(W)) b. Moving a loading column toward a low-variance or
null-space direction of X shrinks z_j and sigma_j together. The standardized score z_j / sigma_j barely moves (net
angle at most 0.13 degrees), while log sigma falls by 6 to 9, and more than 90% of that loss comes from g(v), not
from ||Pw||. The 1/sigma factor in eta and in the rho term produces the 1e6-level curvature. Solving faster only moves
faster along this direction. The fix is to remove 1/sigma(W) from the objective.

## 1. Model, rank and normalization

- **Rank:** k = 3, fixed in both arms (gate spec section 1).
- **Data:** X is p x n within-sample ranks (average ties) of the training columns, p = 3,000 genes (the frozen
  `ngene = 3000` filter), as in `palm_rank` / `palm_data`. C is the n x n sample-centering matrix. P is the p x p
  gene-centering matrix. Because X has constant column sums (ranks), C X^T 1 = 0.
- **Constants:** sigma_bar^2 = mean gene variance of X (`sb2` in `palm_lib.R`), and s0 = sigma_bar. For a loading
  vector w drawn uniformly from the unit sphere, E Var(X^T w) = sigma_bar^2, so s0 is the score SD of a typical
  unit-norm loading. s0 is a data constant: it does not depend on W.
- **W normalization (kept):** W lies in S_+^k = {W >= 0, ||w_j||_2 = 1 for every j}. This is a constraint enforced
  by projection, not a post-hoc rescaling. It fixes the W-H and W-beta gauges: H and beta carry all scale.
- **H normalization:** H >= 0, with no norm constraint and a ridge lambdaH = 1e-6 (unchanged). H is not rescaled.
- **Scores:** Z_c(W) = C X^T W (n x k), the centered scores of the unit-norm loadings. No division by the score SD
  anywhere in the objective.

## 2. Objective (beta parameterization)

F(W, H, beta) = (1 - alpha) [ ||X - W H||_F^2 / (2 ||X||_F^2) + lambdaH ||H||_F^2 / (2 n k) ]
              + alpha [ -(2/n_e) l(Z_c(W) beta)
                        + lambda ( nu ||s0 beta||_1 + (1 - nu)/2 ||s0 beta||_2^2 )
                        + (rho sigma_bar^2 / 2) ||P W beta||_2^2 ]
              + iota_{S+}(W) + iota_{>=0}(H),

where l is the Breslow partial log-likelihood, n_e is the number of events and beta is in R^k. No zeta rescaling.
beta is optimized and penalized directly: the identity beta = b / sigma(W) appears nowhere.

Equivalent form, for intuition only: with theta = s0 beta, the penalty is lambda(nu ||theta||_1 +
(1 - nu)/2 ||theta||^2) + (rho/2) ||P W theta||^2 and eta = Z_c theta / s0. The scaling of the penalty by s0
(draft choice D1, section 10) sets the units. It is a fixed constant and does not depend on W.

**What the redesign does to the pathology.** Rotating w_j toward a direction with small g keeps eta fixed only if
beta_j grows by sigma_j(old) / sigma_j(new). The lambda term and the rho term both grow with beta_j^2, so the
preserved-score direction is now penalized instead of free. A column with beta_j = 0 can still sit in a
low-variance direction. That column is a reconstruction-only factor, which is legitimate.

**Structure.** For fixed (H, beta), the survival term is convex in W because eta is linear in W. The rho and NMF
terms are convex quadratics in W, so the only nonconvexity in the W block is the sphere constraint. For fixed W,
the beta block is the convex elastic-net Cox problem plus a positive semidefinite quadratic, strongly convex when
nu < 1. F is continuous and coercive in (H, beta) and W ranges over a compact set, so a minimizer exists and every
sublevel set is compact. The Cox log-likelihood is real-analytic in beta and W, and S_+^k and the nonnegative
orthant are semialgebraic, so F has the Kurdyka-Lojasiewicz property. Bolte, Sabach and Teboulle (2014) then give
whole-sequence convergence of PALM to a critical point, subject to the inexact b-block caveat already stated in
gate spec amendment 2.

**Curvature bound (σ-free).** The Hessian of -(2/n_e) l in eta has largest eigenvalue at most 1: each event
contributes a softmax covariance with eigenvalue at most 1/2, and there are n_e events, with Breslow ties included.
On the W block (vectorized),

L_W <= (1 - alpha) ||H||_2^2 / ||X||_F^2 + alpha ||beta||_2^2 ( lambda_1(X C X^T) + rho sigma_bar^2 ).

Every term of F is nonnegative, so on the sublevel set {F <= F_0}, where F_0 is F at the start (beta = 0),
||beta||_2^2 <= 2 F_0 / (alpha lambda (1 - nu) s0^2), which is finite because nu <= 0.95. Both bounds are
polynomial in the iterates and contain no sigma. The old objective's W-block curvature grew like sigma^-4. This
bound is a numerical check in qualification Q4.

**Prediction vector.** v = W beta and v~ = P v over the frozen training gene universe. The scored predictor is
X_G^T v~_G, where G is the support chosen by the unchanged ntop rule (`desurv_get_top_genes` on W). Centering
happens before the support restriction, and v~ is never recentered on another dataset's genes. This is the
amendment-2 rule with W diag(1/sigma) b replaced by W beta.

**Reported effects (post-fit only).** beta as fitted, plus standardized Cox effects beta_j sd(Z_c[, j]) for
interpretation. sigma enters reporting only, never the objective, tuning or start selection.

**Invalid input** fails deterministically, with no fallback constant, when any of the following holds: n < 2,
n_e = 0, any nonfinite entry of X, y or d, ||X|| = 0, or sigma_bar = 0. The eps constant of the frozen candidate is
removed because nothing in F divides by a data-dependent quantity.

**Fallback (not co-primary).** If the beta parameterization fails qualification, the ruling's fallback is a
denominator-floor variant (eta = Z_c diag(1/max(sigma_j, tau)) b). It is not specified here. Specifying it requires
a new PI ruling and its own frozen text.

## 3. Penalties and tuning ranges

| Hyperparameter | Role | Range (inner BO) | Scale |
|---|---|---|---|
| alpha | survival vs reconstruction weight | [1e-4, 0.9] | logit |
| lambda | elastic-net strength on s0 beta | [1e-3, 1e3] | log10 |
| nu | L1 share of the elastic net | [0, 0.95] | linear (nu < 1 keeps the beta block strongly convex) |
| rho | gene-space ridge on P W beta | [1e-2, 1e3] | log10 |
| ntop | genes in the scored support | [50, 300] | integer |
| lambdaH | H ridge | 1e-6 fixed | |
| k | rank | 3 fixed | |

**Why these ranges carry over.** With the s0 scaling, the lambda and rho terms are numerically equal to the frozen
candidate's terms at any W where sigma_j(W) = s0. They differ only by charging low-variance factors more per unit of
effect, which is the intended change. The frozen box is therefore on a comparable scale. Bounds are not changed
after any outcome is seen, and boundary hits are reported. If the PI rules raw beta instead (open item 2), the
lambda range must be restated, because the elastic-net mix then changes meaning (section 10, item 2).

**Rho term.** (rho sigma_bar^2 / 2) ||P W beta||^2 is a ridge on the gene-level coefficient vector v = W beta, with
the uniform gene direction removed. With beta = b / sigma(W) it is identical to the frozen term
(rho sigma_bar^2 / 2) ||P W diag(1/sigma) b||^2, so the scaling constant sigma_bar^2 / 2 keeps its meaning. It
penalizes the components of v that do not change eta (the null space of C X^T), and it makes the null-space rotation
costly.

## 4. Optimizer

The frozen PALM mechanics (gate spec section 4, steps 1-3, amendment 2) are reused unchanged except for the
variable. This reuses existing block steps and is not new solver work. Newton-CG, Schur-complement and outer
acceleration are excluded (ruling 8.1).

1. **H:** one HALS sweep, rows j = 1..k in order.
2. **W:** projected gradient step onto S_+^k (column-wise closed form, including the argmax rule when no entry is
   positive), with step t = min(2 t_prev, t_max), t_max = 2 and t_prev = 1 at the first iteration, halved until the sufficient-decrease condition holds,
   and t_min = 1e-20 counted as a numerical failure.
3. **beta:** convex elastic-net Cox plus the rho quadratic, solved by FISTA warm-started from the current beta, with
   an independent KKT certificate ||KKT_beta||_inf <= 1e-10 and a 1e5-iteration cap (reaching the cap is a
   numerical failure). beta+ is accepted only if F does not increase. Otherwise beta is kept and a rejected step is
   recorded.

The block order is H, W, beta. No outer inertial extrapolation.

## 5. Initialization

Unchanged from the frozen candidate except that beta replaces b. Each fit uses 30 starts. Start s of a fit with seed
S runs `set.seed(S + s - 1)`, draws W0 = runif(p k, 0, max(X)), then draws H0 = runif(n k, 0, max(X)) and discards
it so the random stream matches the comparator. W0 is projected onto S_+^k. H is then set by 20 HALS sweeps at fixed
W, and beta = 0. The seed namespaces of gate spec section 4 are unchanged: inner-fold 300101-300834, candidate outer
final fit 10,001,100-10,004,529, Gate B 3000-3429 and 20,000,000-20,000,029. Before freeze, check that none of these
seeds has been used on real data. Training C is never used to choose among starts.

At beta = 0 the first W step has no survival gradient, so the survival term enters from the first beta block on.
This was also true of the frozen candidate with b = 0.

## 6. Convergence and KKT definition

r = max(r_W, r_H, r_beta), with the amendment-2 definitions and beta in place of b:

- r_W = ||R_W||_F, the KKT residual on {w_j >= 0, ||w_j|| = 1}. R_W[i,j] = G_W[i,j] - (w_j^T G_W[,j]) W[i,j] if
  W[i,j] > 0, and min(G_W[i,j], 0) if W[i,j] = 0.
- r_H = ||min(H, G_H / a_H)||_F / max(1, ||H||_F), with a_H = (1 - alpha)(1/||X||^2 + lambdaH/(n k)).
- r_beta = ||KKT_beta||_inf of the alpha-free beta subproblem.

Stop when r < 1e-5. maxit = 5,000 outer iterations, which restores the pre-amendment-3 budget. An F increase is
F_new > F_old + 1e-13 max(1, |F_old|), checked after each block. Start statuses (converged, nonconverged, numerical
failure), the fit-level rule (lowest-final-F converged start, otherwise flagged nonconverged, otherwise failed),
diagnostics and failure propagation are those of gate spec section 4 and amendment 2. In addition, each start
records sigma_j(W) and ||beta|| at every 10th sweep. These values are for qualification diagnostics only and bind
nothing.

## 7. Nested tuning and the comparator

- **Candidate:** inner BO exactly as in gate spec section 4, over the section 3 box above. It uses the same BO
  engine, candidate pool, exploration weight and seed rule, a budget of n_init = 60 and n_iter = 120, 5-fold inner CV
  with 30 starts per fold, and the inner-fold seeds above. A configuration is feasible only if all five inner-fold
  fits converge. Infeasible configurations are NA with a flag and are excluded from the surrogate and from
  selection. The selected configuration is the feasible one with the highest mean inner-CV C, kept at full
  precision. The final outer-training fit uses 30 starts and is scored with v~ on G.
- **Comparator:** gate spec section 3, unchanged. It is the production DeSurv 1.0.1 procedure (sha `afb00d5`) with
  k = 3, tuned by a fresh inner BO in every outer fold, followed by the `code/04` consensus final fit. It is not the
  PALM candidate and not the published hyperparameters. This reading awaits PI confirmation (open item 4).

## 8. Gates

- **Gate A (unchanged):** the candidate is eligible if mean_f dC_f >= -0.02 over the 20 outer folds, with the
  Nadeau-Bengio corrected interval reported in full and no noninferiority language. Convergence precondition: at
  least 90% of outer-fold final fits are converged under section 6, no F increase is recorded, and no fit failed.
  The Gate A outcome is computed by a script that reads only C values and convergence records, and it is committed
  together with this file's hash.
- **Gate B (hidden until Gate A is committed):** gate spec section 7 and Appendix A are carried over by reference
  without change. No Gate B script is written or run, and no real-data biology object (the published W, D1 loadings,
  reference programs or external cohorts) is loaded by any candidate code, including qualification code, until the
  Gate A outcome file exists in a commit. The Gate B script refuses to run otherwise (gate spec section 10).
- **Decision table:** gate spec section 8, unchanged.

## 9. Qualification on simulated data (before any Gate A compute)

Simulated data only, from the repo generator through `analysis/estimator-gate/candidate/sim.R` (scenario R0, within-sample ranks, integer-rounded
tied times), at p = 3,000. A qualification prespec fixes the replicate counts, grid and generator knobs below and is
committed before any qualification run. Thresholds marked (P) are draft proposals for the PI to rule on. Nothing is
tuned on qualification outcomes. A failure is reported and returned to the PI.

**Scenarios.**
- S1, signal: R0 with generator beta = (2, 1, 0).
- S2, null: R0 with beta = (0, 0, 0).
- S3, trap: R0 with the carrying factor's score variance reduced so its true score SD is at most 0.25 times that of
  the other factors (P). The generator knob is to be fixed in the prespec. This scenario reproduces the condition in
  which a low-variance direction carries survival signal.

**Sizes.** n = 175 (inner-training size) and n = 218 (outer-training size). Each training set gets an independent
test set of n = 1,000 from the same generator.

**Configurations.** A fixed grid inside the section 3 box: alpha in {1e-3, 0.05, 0.5} x log10 lambda in {-2, 0, 2}
at nu = 0.5, log10 rho = 0 and ntop = 150, plus log10 rho in {-2, 3} at the centre, giving 11 configurations (P).
The grid runs on 5 replicate datasets per scenario at n = 175 (P), and the centre configuration also runs at
n = 218. One reduced-budget nested tuning run per scenario (inner BO n_init = 20, n_iter = 40, P) provides the tuned
configuration for Q5.

**Criteria.** "Carrying column" means |s0 beta_j| > 1e-3 at the endpoint.

| ID | Tests | Pass if |
|---|---|---|
| Q1 | Ordinary convergence | At every grid configuration, at least 90% of starts are converged (r < 1e-5 within 5,000 iterations), with no F increase and no numerical failure (P: the 90% mirrors the Gate A precondition) |
| Q2 | No sigma -> 0 | For every carrying column of every converged start, log(sigma_j(end) / sigma_j(sweep 50)) >= -1, and over the last 25% of sweeps the absolute change in log sigma_j is at most 0.01 (P) |
| Q3 | No preserved-score collapse | Reusing the D-Q4 decomposition (`analysis/estimator-gate/candidate/tests/DQ4_step_decomposition.R`; log sigma = log r + log g and the score angle depend on W only, so the adaptation changes only the objective-term accounting, b to beta), no carrying column has cumulative change in log g <= -1 while its net standardized-score angle is <= 0.01 rad, over the whole trajectory or over the last 25% of sweeps (P) |
| Q4 | Bounded curvature | The largest eigenvalue of the W-block Hessian on the active face, by power iteration on HVPs at sweeps 1, 10, 100 and the end, is (a) always at most the section 2 bound L_W (a violation is an implementation defect, not a qualification outcome) and (b) within a factor of 100 between its smallest and largest recorded values per start (P). The slope of log lambda_max on log min_j sigma_j is reported |
| Q5 | Multi-start stability | At the centre and tuned configurations, among converged starts, the median Hungarian-matched loading Pearson r of carrying columns with the best-F start is >= 0.90, and the range of test C across converged starts is <= 0.02 (P). The share of starts within 1e-4 relative F of the best is reported |
| Q6 | Predictive recovery at convergence | At the tuned configuration: test C >= oracle C (true linear predictor) - 0.03 in S1 and S3; mean test C in [0.47, 0.53] in S2. The early-stopping check also applies: mean test C at convergence >= max over checkpoints (sweeps 10, 50, 100, 500) - 0.01 (P). Matched cosine of fitted vs true W is reported |

Q6 carries the F2 lesson (open item 1). The F2 audit found that exact optimization of the published objective
overfit, with CV falling from 0.638 to 0.571, because the stall had acted as implicit regularization. A clean
convergence that loses discrimination to its own early iterates means the stall was still doing the regularizing,
and the estimator should fail.

**Risk to watch (descriptive).** Raw-beta penalization makes high-variance score directions cheap. Factors may
therefore drift toward the leading principal directions of X, at the expense of reconstruction structure. The matched
W cosine in Q6 and the NMF term at the endpoint are reported to show whether this happens.

**Outcome.** If every criterion passes, the PI freezes this spec and Gate A compute is authorized. Otherwise the
failure is reported, mechanics are not adjusted to the outcome, and the PI decides whether to open the fallback.

## 10. Open items for PI decision

From the handoff (coordinator notes, not decisions):

1. **F2 link.** Exact optimization of the published objective overfit (CV 0.638 -> 0.571), because the stall acted
   as implicit regularization. The draft addresses this through Q6: held-out C at convergence is compared with the
   oracle and with the run's own early iterates, in addition to clean convergence. Decide whether Q6 as written is
   sufficient, and whether its margins (0.03 to the oracle, 0.01 to the early iterates) are right.
2. **Scale of the beta penalty.** Penalizing raw beta is scale-dependent: a factor with larger score variance pays
   less penalty per unit of effect. That dependence is what blocks the collapse, because a per-factor rescale by
   sigma_j(W) would reintroduce the singularity. What remains open is the units. The draft scales by the fixed data
   constant s0 = sigma_bar (draft choice D1), which keeps the frozen lambda box comparable. The alternative is raw
   beta with a restated lambda range. The two are not a pure relabelling of lambda when 0 < nu < 1. Decide s0 or
   raw.
3. **Rho term scaling.** The term moves from P W diag(1/sigma) b to P W beta. Under beta = b / sigma(W) it is the
   same quantity, so the draft keeps the constant sigma_bar^2 / 2 and the frozen rho range [1e-2, 1e3]
   (section 3). Confirm, or restate.
4. **Comparator.** Confirm that the "current implemented estimator" is the production DeSurv fitting procedure, as
   in gate spec section 3: k = 3, a fresh inner BO in every outer fold, and the `code/04` consensus final fit. It is
   not the PALM candidate and not the published hyperparameters.

Draft choices made here that the PI may override:

- **D1.** The s0 penalty scaling (item 2).
- **D2.** The PALM block steps are reused with maxit = 5,000. The amendment-3 budget is dropped together with
  Newton-CG.
- **D3.** W keeps the unit-sphere constraint, as ruled. A relaxation to the unit ball would make the W block fully
  convex, and every term of F favors ||w_j|| = 1 whenever H_j or beta_j is nonzero. It is noted here, not adopted.
- **D4.** All (P) thresholds, the S3 generator definition and the qualification grid size (section 9). Compute for
  the grid is unknown until a timing run, which itself needs authorization.
