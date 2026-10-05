# Estimator-redesign specification: rank-3 DeSurv in the theta parameterization (2026-10-05)

**Status: DRAFT (revision 4: revision 3 with the H-ridge rescaling of PI ruling 12), awaiting independent Codex review, not frozen.** Written under PI ruling 8
(`docs/HANDOFF-2026-10-05-candidate-optimizer.md`, "PI ruling 2026-10-05, eighth"). Revision 2 applied PI ruling 9
revision 3 applies PI ruling 10, and revision 4 applies PI ruling 12 (all in section 10). Next: an independent Codex review, then the PI freezes
this file. No implementation, no real data and no compute until then. Once frozen, every later output records this
file's commit hash. Amendments then follow the rule in `2026-10-04-estimator-gate-spec.md`: dated, with a reason,
and only before any Gate A output is read.

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

**Scope of the D-Q evidence (ruling 12).** D-Q1 to D-Q4 diagnosed the old candidate under the unscaled H penalty
(lambdaH ||H||^2 / (2 n k), a HALS denominator of about 1 + 3,000 at p = 3,000). Those diagnostics remain valid as
evidence about that failed estimator. They are not assumed to describe the redesigned estimator, which also changes
the H penalty scale (section 1).

## 1. Model, rank and normalization

- **Rank:** k = 3, fixed in both arms (gate spec section 1).
- **Data:** X is the p x n matrix of within-sample ranks (average ties) of the fit's training columns, p = 3,000 genes
  (the frozen `ngene = 3000` filter), as in `palm_fold_data` (`palm_rank` on the training columns, then `palm_data`). C is the n x n sample-centering matrix.
  P is the p x p gene-centering matrix. Because X has constant column sums (ranks), C X^T 1 = 0.
- **Unit constant s0 (fixed, not adaptive).** s0^2 = sigma_bar^2 = (1/p) sum_{i=1..p} var_i(X), the mean over genes
  of the sample variance (denominator n - 1) of each gene's row of X. Equivalently,
  s0^2 = tr(X C X^T) / (p (n - 1)). This is `sb2 = mean(apply(X, 1, var))` in `palm_data` (`palm_lib.R:31`).
  - It is computed from X only, once per fit, on that fit's training columns: the inner-training set for an
    inner-fold fit, the outer-training set for an outer final fit, and all 273 samples for the Gate B full-data fit.
  - It does not depend on W, y, d, any hyperparameter or any start. Every start and every BO evaluation of a fit
    uses the same s0.
  - It is stored with the fit and never recomputed on validation, test or external data. Scoring new data does not
    use s0 except through the stored beta (section 2).
  - Interpretation: if w is uniform on the unit sphere, E[w w^T] = I/p, so
    E Var(X^T w) = tr(X C X^T) / (p (n - 1)) = s0^2. s0 is the score SD of a typical unit-norm loading. This is a
    property of the definition, not part of it.
- **W normalization (kept, ruling 9 D3):** W lies in S_+^k = {W >= 0, ||w_j||_2 = 1 for every j}. This is a
  constraint enforced by projection, not a post-hoc rescaling. It fixes the W-H and W-theta gauges: H and theta carry
  all scale. The unit ball is not used, because it would change the estimator and add a scale direction.
- **H normalization (ruling 12):** H >= 0, with no norm constraint and the ridge lambdaH ||H||_F^2 / (2 ||X||_F^2),
  lambdaH = 1e-6. H is not rescaled. The ridge is divided by ||X||_F^2, like the reconstruction term, so lambdaH is
  dimensionless on the scale of the normalized reconstruction loss. The frozen candidate divided it by n k instead.
  For ranks at p = 3,000 that made the HALS denominator 1 + lambdaH ||X||_F^2 / (n k), about 1 + 3,000, which
  crushed H and left W governed almost entirely by the survival and penalty terms. The numerical value 1e-6 is
  kept and not retuned.
- **Scores:** Z_c(W) = C X^T W (n x k), the centered scores of the unit-norm loadings. No division by any
  W-dependent score SD anywhere in the objective.

## 2. Objective (theta parameterization)

The optimized survival coordinate is theta in R^k, with eta = Z_c(W) theta / s0. The reported raw-score Cox
coefficient is beta = theta / s0.

F(W, H, theta) = (1 - alpha) [ ||X - W H||_F^2 / (2 ||X||_F^2) + lambdaH ||H||_F^2 / (2 ||X||_F^2) ]
               + alpha [ -(2/n_e) l(Z_c(W) theta / s0)
                         + lambda ( nu ||theta||_1 + (1 - nu)/2 ||theta||_2^2 )
                         + (rho / 2) ||P W theta||_2^2 ]
               + iota_{S+}(W) + iota_{>=0}(H),

where l is the Breslow partial log-likelihood and n_e is the number of events. There is no zeta rescaling. Because
s0^2 = sigma_bar^2, the rho term equals the frozen form (rho sigma_bar^2 / 2) ||P W beta||^2, and the identity
beta = b / sigma(W) appears nowhere. s0 is a fixed unit conversion (section 1), not another adaptive normalization.

**What the redesign does to the pathology.** Rotating w_j toward a direction with small g keeps eta fixed only if
theta_j grows by sigma_j(old) / sigma_j(new). The lambda term and the rho term both grow with theta_j^2, so the
preserved-score direction is now penalized instead of free. A column with theta_j = 0 can still sit in a low-variance
direction. That column is a reconstruction-only factor, which is legitimate. The opposite bias is now possible:
high-variance score directions are cheap per unit of effect. Scenario S3 (section 9) tests this.

**Structure and convergence status (ruling 9, theorem language).** The redesigned objective removes the
1/sigma(W) singularity. Under the stated constraints and a positive penalty (lambda > 0, nu <= 0.95), it is bounded
below by 0 (every term is nonnegative), coercive in (H, theta), and defined on a compact W set, so its sublevel sets
are bounded. For fixed (H, theta), the survival, rho and NMF terms are convex in W, because eta is linear in W, and
the only nonconvexity in the W block is the sphere constraint. For fixed W, the theta block is the convex elastic-net
Cox problem plus a positive semidefinite quadratic, strongly convex when nu < 1.

**Convergence claim (ruling 10, 9d).** The block updates are not claimed to instantiate the exact Bolte, Sabach
and Teboulle (2014) PALM iteration: H is an exact row minimization (HALS), W a projected gradient step with
backtracking, and theta an inexact convex solve to a KKT certificate with an F-acceptance check. Qualification
therefore relies on monotone objective descent, bounded iterates, independently checked block and KKT residuals,
and convergence of the full residual r (section 6). Formal convergence under the mixed exact and inexact block
updates is not claimed here. The H block is not changed to fit a theorem. A broader block-coordinate convergence
result may be investigated later; it is not a precondition for implementation.

**Curvature bound (sigma-free).** The Hessian of -(2/n_e) l in eta has largest eigenvalue at most 1: each event
contributes a softmax covariance with eigenvalue at most 1/2, and there are n_e events, with Breslow ties included.
On the W block (vectorized),

L_W <= (1 - alpha) ||H||_2^2 / ||X||_F^2 + alpha ||theta||_2^2 ( lambda_1(X C X^T) / s0^2 + rho ).

On the sublevel set {F <= F_0}, where F_0 is F at the start (theta = 0), ||theta||_2^2 <= 2 F_0 / (alpha lambda
(1 - nu)). Both bounds are polynomial in the iterates and contain no sigma, whereas the old objective's W-block
curvature grew like sigma^-4. The first bound is a floating-point check in qualification Q4. It is an upper bound,
not a prediction of the observed curvature, which is expected to be much smaller.

**Prediction vector.** v = W theta / s0 = W beta, and v~ = P v over the frozen training gene universe. The scored
predictor is X_G^T v~_G, where G is the support chosen by the unchanged ntop rule (`desurv_get_top_genes` on W).
Centering happens before the support restriction, and v~ is never recentered on another dataset's genes. This is the
amendment-2 rule with W diag(1/sigma) b replaced by W beta. The fitted model stores the training gene universe, v~,
G, theta, beta and s0.

**Reported effects (post-fit only).** beta = theta / s0 (raw-score Cox coefficient), plus standardized effects
beta_j sd(Z_c[, j]) for interpretation. sigma_j(W) enters reporting and qualification diagnostics only, never the
objective, tuning or start selection.

**Invalid input** fails deterministically, with no fallback constant, when any of the following holds: n < 2,
n_e = 0, any nonfinite entry of X, y or d, ||X|| = 0, or s0 = 0. The eps constant of the frozen candidate is removed
because nothing in F divides by a W-dependent quantity.

**Fallback (not co-primary).** If the theta parameterization fails qualification, the ruling's fallback is a
denominator-floor variant (eta = Z_c diag(1/max(sigma_j, tau)) b). It is not specified here. Specifying it requires
a new PI ruling and its own frozen text.

## 3. Penalties and tuning ranges

| Hyperparameter | Role | Range (inner BO) | Scale |
|---|---|---|---|
| alpha | survival vs reconstruction weight | [1e-4, 0.9] | logit |
| lambda | elastic-net strength on theta | [1e-3, 1e3] | log10 |
| nu | L1 share of the elastic net | [0, 0.95] | linear (nu < 1 keeps the theta block strongly convex) |
| rho | gene-space ridge on P W theta | [1e-2, 1e3] | log10 |
| ntop | genes in the scored support | [50, 300] | integer |
| lambdaH | H ridge, divided by ||X||_F^2 (ruling 12) | 1e-6 fixed | |
| k | rank | 3 fixed | |

**Why these ranges carry over (ruling 9, items 2 and 3).** At any W where sigma_j(W) = s0, the lambda and rho terms
are numerically equal to the frozen candidate's terms. They differ only by charging low-variance factors more per
unit of effect, which is the intended change. The frozen box is therefore on the same scale, and it is not altered
in advance because of what the new model is expected to do. Bounds are not changed after any outcome is seen.
Boundary hits are reported. If qualification or tuning shows persistent boundary selection, it is recorded and
decided separately.

**Rho term.** (rho / 2) ||P W theta||^2 = (rho sigma_bar^2 / 2) ||P W beta||^2 is a ridge on the gene-level
coefficient vector, with the uniform gene direction removed. Under beta = b / sigma(W) it is identical to the frozen
term (rho sigma_bar^2 / 2) ||P W diag(1/sigma) b||^2. It penalizes the components of v that do not change eta (the
null space of C X^T), and it makes the null-space rotation costly.

## 4. Optimizer (ruling 9 D2: approved)

The frozen PALM mechanics (gate spec section 4, steps 1-3, amendment 2) are reused unchanged except for the
variable and the H-ridge scale (ruling 12, step 1). This reuses existing block steps and is not new solver work. Newton-CG, Schur-complement and outer
acceleration are excluded (ruling 8.1).

1. **H:** one HALS sweep, rows j = 1..k in order. Because ||w_j|| = 1 and the H terms share the factor
   (1 - alpha) / ||X||_F^2, the exact row update is h_j = max(0, w_j^T R_j / (1 + lambdaH)), with
   R_j = X - sum_{l != j} w_l h_l (ruling 12; the frozen candidate's denominator was 1 + lambdaH ||X||_F^2 / (n k)).
2. **W:** projected gradient step onto S_+^k (column-wise closed form, including the argmax rule when no entry is
   positive), with step t = min(2 t_prev, t_max), t_max = 2 and t_prev = 1 at the first iteration, halved until the
   sufficient-decrease condition holds, and t_min = 1e-20 counted as a numerical failure.
3. **theta:** convex elastic-net Cox plus the rho quadratic, solved by FISTA warm-started from the current theta, with
   an independent KKT certificate ||KKT_theta||_inf <= 1e-10 and a 1e5-iteration cap (reaching the cap is a
   numerical failure). theta+ is accepted only if F does not increase. Otherwise theta is kept and a rejected step is
   recorded.

The block order is H, W, theta. No outer inertial extrapolation.

## 5. Initialization

Unchanged from the frozen candidate except that theta replaces b. Each fit uses 30 starts. Start s of a fit with seed
S runs `set.seed(S + s - 1)`, draws W0 = runif(p k, 0, max(X)), then draws H0 = runif(n k, 0, max(X)) and discards
it so the random stream matches the comparator. W0 is projected onto S_+^k. H is then set by 20 HALS sweeps at fixed
W, and theta = 0. The seed namespaces of gate spec section 4 are unchanged: inner-fold 300101-300834, candidate outer
final fit 10,001,100-10,004,529, Gate B 3000-3429 and 20,000,000-20,000,029. Before freeze, check that none of these
seeds has been used on real data. Training C is never used to choose among starts.

At theta = 0 the first W step has no survival gradient, so the survival term enters from the first theta block on.
This was also true of the frozen candidate with b = 0.

## 6. Convergence and KKT definition

r = max(r_W, r_H, r_theta), with the amendment-2 definitions and theta in place of b:

- r_W = ||R_W||_F, the KKT residual on {w_j >= 0, ||w_j|| = 1}. R_W[i,j] = G_W[i,j] - (w_j^T G_W[,j]) W[i,j] if
  W[i,j] > 0, and min(G_W[i,j], 0) if W[i,j] = 0.
- r_H = ||min(H, G_H / a_H)||_F / max(1, ||H||_F), with a_H = (1 - alpha)(1 + lambdaH) / ||X||^2 (ruling 12).
- r_theta = ||KKT_theta||_inf of the alpha-free theta subproblem.

The estimator is defined by the certificate: a start is *converged* when r < 1e-5 with no recorded F increase. The
5,000-sweep cap is a work limit, not part of the mathematical definition (ruling 9 D2). A start that reaches it
without a numerical failure is *nonconverged* and is recorded as such. If ordinary PALM cannot reach the certificate
within the cap, that is reported. It is not a reason to rebuild a second-order solver. An F increase is
F_new > F_old + 1e-13 max(1, |F_old|), checked after each block. Start statuses (converged, nonconverged, numerical
failure), the fit-level rule (lowest-final-F converged start, otherwise flagged nonconverged, otherwise failed),
diagnostics and failure propagation are those of gate spec section 4 and amendment 2. In addition, each start
records sigma_j(W) and ||theta|| at every 10th sweep, and keeps (W, theta) at the section 9 checkpoints. These are for
qualification diagnostics and bind nothing in Gate A.

## 7. Nested tuning and the comparator

- **Candidate:** inner BO exactly as in gate spec section 4, over the section 3 box above. It uses the same BO
  engine, candidate pool, exploration weight and seed rule, a budget of n_init = 60 and n_iter = 120, 5-fold inner CV
  with 30 starts per fold, and the inner-fold seeds above. Each inner-fold fit computes its own s0 from its
  inner-training columns. A configuration is feasible only if all five inner-fold fits converge. Infeasible
  configurations are NA with a flag and are excluded from the surrogate and from selection. The selected
  configuration is the feasible one with the highest mean inner-CV C, kept at full precision. The final
  outer-training fit uses 30 starts, computes s0 from the outer-training columns, and is scored with v~ on G.
- **Comparator (ruling 9, item 4: confirmed):** gate spec section 3, unchanged. It is the production DeSurv 1.0.1
  procedure (sha `afb00d5`) with k = 3, tuned by a fresh inner BO in every outer training fold, followed by the
  production consensus final fit of `code/04`, on the same outer folds as the candidate. It is not the PALM candidate,
  not the published hyperparameters and not DeSurv 1.2.0. The corrected tie derivative in 1.2.0 may be reported as a
  sensitivity analysis, since its measured effect was negligible. It does not replace the comparator.

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

Simulated data only, at p = 3,000. A qualification prespec fixes the replicate counts, grid and the S3 wrapper below
and is committed before any qualification run. Nothing is tuned on qualification outcomes. A failure is reported and
returned to the PI. After this spec is frozen, a timing-only simulated run may be made to plan compute and
parallelism (ruling 9). Its outputs are wall times only and are never used to change the objective, tuning ranges,
thresholds or algorithm.

**Three kinds of criterion (ruling 9 D4).**
- *Mathematical checks* use tolerances set by floating-point accuracy. A failure is an implementation defect, not a
  qualification outcome.
- *Pathology criteria* (binding) test whether a failure mechanism recurs: the D-Q4 collapse (cost-free loss of
  scale with the same score preserved), its opposite bias (the effect attached to a high-variance nonpredictive
  direction), or the F2 overfit. They are not scalar ranges copied from the failed run (ruling 10, 9b).
- *Validity measures* (reported, not binding) show whether the simulation was informative and how the fit
  represents the truth: the oracle C gap, score and loading recovery, realized curvature and sigma trajectories.

**Scenarios** (all built on the repo generator through `analysis/estimator-gate/candidate/sim.R`, scenario R0, with
within-sample ranks and integer-rounded tied times):
- **S1, signal:** R0 with generator beta = (2, 1, 0).
- **S2, null:** R0 with beta = (0, 0, 0).
- **S3, competing directions (ruling 10, 9c):** a genuinely predictive low-variance factor A competes with a
  higher-variance nonpredictive factor B. The generator computes the survival linear predictor on standardized
  scores (`simulate_survival.R:189-200`), so scaling the rows of the true H changes a factor's variance in X without
  changing its predictive strength. S3 is R0 with generator beta = (2, 0, 0), where factor 1 is A and factor 2 is B,
  and true H rows multiplied by c = (c_A, c_B, 1) before X is formed, with c_A and c_B set so that A's true score SD
  is 0.25 times B's (variance ratio 1/16). This is the primary and only binding stress case. A and B are kept
  uncorrelated by construction: R0's H rows are independent (`correlated_H = FALSE`) and markers do not overlap
  (`marker_overlap = 0`). The realized test-set correlation of the true A and B scores is reported. The generator has
  no per-factor variance argument, so the qualification prespec specifies the wrapper and checks the realized SD
  ratio. S3 tests the bias opposite to the old estimator's: direct theta penalization could favor B, because a
  smaller coefficient on a high-variance score buys the same linear-predictor magnitude.
- **S3b, sensitivity (reported only):** as S3 with an SD ratio of 0.50, run if compute allows. No binding
  criterion applies to S3b.

**Sizes.** n = 175 (inner-training size) and n = 218 (outer-training size). Each training set gets an independent
test set of n = 1,000 from the same generator.

**Configurations.** A fixed grid inside the section 3 box: alpha in {1e-3, 0.05, 0.5} x log10 lambda in {-2, 0, 2}
at nu = 0.5, log10 rho = 0 and ntop = 150, plus log10 rho in {-2, 3} at the centre, giving 11 configurations. The grid
runs on 5 replicate datasets per scenario at n = 175, and the centre configuration also runs at n = 218. One nested
tuning run per scenario and replicate, with a reduced inner BO (n_init = 20, n_iter = 40), gives the tuned
configuration. The grid and the reduced BO budget set compute, not decisions. The final counts are fixed in the
qualification prespec.

**Checkpoints (frozen here).** Sweeps T = {10, 50, 100, 500}. For a fit, the checkpoint iterate at sweep t is the
returned start's (W, theta) at sweep t. If that start converged before sweep t, it is the converged iterate (the
difference is then 0). No checkpoint is chosen after held-out C is seen.

A *carrying column* is one with |theta_j| > 1e-3 at the endpoint, a numerical zero test on the theta scale.

| ID | Kind | Tests | Criterion and reason |
|---|---|---|---|
| Q0 | Math | Implementation | Gradient checks against finite differences on every block, F monotone, theta-block KKT verified by an independent evaluator: floating-point tolerances as in gate spec section 4 preconditions |
| Q1 | Pathology | Ordinary convergence | At every grid configuration and every tuned configuration, at least 90% of starts are converged, with no F increase and no numerical failure. Reason: this is the Gate A convergence precondition applied to the same procedure; a design that fails it on simulations cannot pass it on data |
| Q2 | Validity (demoted, ruling 10) | Scale trajectory | For every carrying column of every converged start, log(sigma_j(end) / sigma_j(sweep 50)) and the sigma_j trace are reported. Not binding: S3 contains a genuinely predictive low-variance factor, and a proper estimator must be allowed to find low variance. The pathology is cost-free collapse with the score preserved, which Q3 tests |
| Q3 | Pathology (binding on mechanism, ruling 11) | No recurrence of the score-preserving collapse | The D-Q4 mechanism is scale lost along a low-variance or null-space direction while the same standardized predictor is preserved, without the run converging. Using the D-Q4 decomposition (`analysis/estimator-gate/candidate/tests/DQ4_step_decomposition.R`; log sigma = log r + log g and the score angle depend on W only, so the adaptation changes only the objective-term accounting, b to theta), a start shows recurrence if it is not converged (section 6) and, for some carrying column over its last 25% of sweeps, all of the following hold: (i) the summed change in log sigma_j is negative; (ii) the loss is g-driven, abs(sum Delta log g) > abs(sum Delta log r); (iii) the score is preserved relative to the scale change, net standardized-score angle (rad) < abs(sum Delta log sigma_j). Each condition is a comparison, not a magnitude. Qualification fails if any start at any grid or tuned configuration shows recurrence. A converged start cannot show recurrence: a large legitimate change in score variance (as S3 invites) passes. Together with Q1 (KKT convergence) and Q4 (the sigma-free analytic curvature bound), this is the binding test of the pathology. Reported descriptively against the failed estimator, never as pass/fail limits: cumulative log g and log sigma changes and net angles per carrying column, compared with the D-Q4 reference (log sigma loss 6 to 9, net angle at most 0.0023 rad), and realized curvature compared with the old 1e6 level |
| Q4 | Math (binding); realized values reported | Bounded curvature | The largest eigenvalue of the W-block Hessian on the active face, by power iteration on HVPs at sweeps 1, 10, 100 and the end, never exceeds the section 2 analytic bound L_W beyond floating-point tolerance. The realized curvature, its ratio across recorded sweeps, the condition number and the slope of log lambda_max on log min_j sigma_j are reported, with no fixed limit (ruling 10: the old factor-100 contrast is historical, not a property of the new estimator) |
| Q5 | Pathology (predictive part binding) | Multi-start stability | At the centre and tuned configurations, among converged starts, the median across starts of the test-C difference from the returned start is within 0.01, the F2 margin, so that the choice among converged starts cannot cost more than one C point. Reported, not binding: the median Hungarian-matched loading Pearson r of carrying columns with the returned start (0.80 is the Gate B reference value, shown for context), and the share of starts within 1e-4 relative F of the best |
| Q6 | Pathology (binding F2 safeguard) | Convergence does not destroy generalization | In S1 and S3, at the tuned configuration, for every checkpoint t in T: mean over replicates of the paired difference C_conv - C_t (test C of the converged returned start minus test C of its checkpoint-t iterate) >= -0.01, with each iterate scored by v~ on its own ntop support. Reason (ruling 9): a properly converged estimator may not give up more than one C point to its own finite trajectory. This is deliberately stricter than Gate A's -0.02 against the published estimator. S2 is reported |
| Q7 | Pathology (binding, predictor level, ruling 10) | The effect lands on the predictive direction, not the high-variance one | In S3, at the tuned configuration, per replicate on the held-out test set: standardize the fitted predictor eta_hat = X_G^T v~_G and the true scores A and B, regress eta_hat on (A, B) jointly, and set D_AB = abs(coef_A) - abs(coef_B). Pass if the one-sided Monte Carlo 95% lower confidence bound for mean D_AB over replicates (mean - t_{0.95, R-1} SD / sqrt(R), R the number of S3 replicates) is above 0. The joint-regression form equals the marginal form abs(cor(eta_hat, A)) - abs(cor(eta_hat, B)) when A and B are uncorrelated and stays unambiguous if they are not; the marginal form is reported alongside |
| V1 | Validity | Oracle gap | Test C of the fit vs test C of the true linear predictor, per scenario and size. Reported only, because it depends on signal strength and censoring (ruling 9, item 1) |
| V2 | Validity | Predictive-score recovery | Pearson and Spearman correlation, on the test set, of the fitted predictor X_G^T v~_G with the true linear predictor, in S1, S3 and S3b; in S3 and S3b also the marginal D_AB |
| V3 | Validity | Loading and program recovery | In S1, S3 and S3b: Hungarian-matched cosine of fitted vs true W columns, whether the column with the largest abs(theta_j) is matched to A (S3, S3b) or to the true predictive factor (S1), and top-ntop marker overlap with the true marker set. Reported, not binding: factorizations can split or mix a predictive direction while producing the correct predictor (ruling 10, 9a). Reported because Gate B concerns loadings, and D-Q4 showed that loadings can move while prediction barely does |

**Outcome.** If every binding criterion passes, Gate A compute is authorized. Otherwise the failure is reported,
mechanics are not adjusted to the outcome, and the PI decides whether to open the fallback.

## 10. PI rulings and open items

**PI ruling 2026-10-05, ninth (on draft revision 1, 8772b86).** The theta redesign is approved in principle, not yet
frozen.
1. **Q6 / F2:** the structure is approved. The binding F2 test is converged versus prespecified early checkpoints, a
   mean paired difference >= -0.01 at every checkpoint, with the checkpoint set frozen in advance and never chosen
   after the fact. The 0.03 margin to the oracle is dropped as a binding criterion and is reported (V1) with score
   and W recovery.
2. **Penalty scale:** keep the fixed s0 = sigma_bar and optimize theta = s0 beta. s0 comes from X only, is computed
   once per training fold, and is independent of W, outcomes, hyperparameters and starts, and is never recomputed on
   validation or test data. Its formula is stated before freezing (section 1).
3. **Rho:** approved. (rho sigma_bar^2 / 2) ||P W beta||^2 = (rho/2) ||P W theta||^2, range [1e-2, 1e3], not altered
   prospectively; persistent boundary selection is recorded and decided separately.
4. **Comparator:** confirmed. DeSurv 1.0.1, k = 3, fresh inner BO in every outer training fold, the production
   `code/04` consensus final fit, same outer folds. 1.2.0 is not the primary comparator; the tie-derivative
   correction may be a sensitivity analysis.
5. **D2:** PALM with maxit 5,000 approved. 5,000 is a work limit, not part of the estimator's definition. Failure to
   reach the certificate is recorded, not answered with Newton-CG.
6. **D3:** keep the unit sphere. No unit ball.
7. **D4:** no arbitrary qualification numbers for completeness. Mathematical checks use floating-point tolerances;
   performance thresholds correspond to the pathology tested. S3 must contain a genuinely predictive low-variance
   direction competing with a higher-variance nonpredictive one. A timing-only simulated run is allowed after
   freeze, for compute planning only.
8. **Theorem language:** no claim that the Bolte-Sabach-Teboulle theorem applies. State coercivity and boundedness,
   and check the PALM assumptions explicitly before implementation (superseded by ruling 10, 9d: no theorem claim, gap recorded).
9. **Recovery:** qualification assesses both predictive-score recovery and loading/program recovery in S1 and S3.

**PI ruling 2026-10-05, tenth (on draft revision 2, 593a4d1).** Pre-freeze changes, after which the spec goes to
an independent Codex review and, if no substantive problem is found, is frozen before any implementation code.
- **9a:** the S3 bias test is binding at the predictor level only (Q7): the Monte Carlo lower 95% confidence bound
  for mean D_AB must exceed 0, with standardized regression coefficients on (A, B) used so the comparison is
  unambiguous if A and B correlate. "Largest abs(theta) matches A" is not binding; it is reported after optimal
  factor matching as a representation-validity measure (V3).
- **9b:** qualification fails on the same mechanism recurring, not on a scalar entering a range reminiscent of the
  old run. The KKT/convergence criterion stays binding (Q1). The raw log sigma limit is reported, not binding (Q2).
  The fixed factor-100 curvature limit is dropped; the derived analytic bound is the binding mathematical check and
  realized curvature and condition number are reported (Q4). In Q5 the 0.01 predictive component is binding as
  part of the F2 safeguard, and the 0.80 loading threshold is descriptive, because rotational ambiguity should not
  sink an otherwise correct estimator.
- **9c:** SD ratio 0.25 is the primary and only binding S3 stress case. 0.50 is an optional reported sensitivity
  (S3b). No 0.10 case. A and B must not be strongly correlated by construction, or the regression form is used.
- **9d:** HALS and the FISTA theta solve are kept. The spec declines a classical PALM theorem claim and relies on
  monotone descent, bounded iterates, independently checked residuals and convergence of the full residual
  (section 2).
- The Q6 checkpoint set {10, 50, 100, 500} stays exactly as frozen; no checkpoints are added after implementation.

**PI ruling 2026-10-05, eleventh (on revision 3, 4c0fa56).** (1) Review route: wait for the independent Codex review
(available 2026-10-09); revision 3 stays unfrozen and nothing is implemented until it is in. (2) Q3 is reworded to bind
only on the mechanism: KKT convergence (Q1), the derived curvature bound (Q4), and no recurrence of the D-Q4
score-preserving collapse (Q3). Historical magnitudes (one sixth of the old log sigma loss, curvature below 1e2, the
old log g and angle limits) are reported against the failed estimator, not used as pass/fail limits. Principle:
qualification fails if the pathology recurs, not because a scalar crosses a number chosen from the previous
failure. No further edits to substantive criteria before the Codex review.

**PI ruling 2026-10-05, twelfth (on review item B1 of `2026-10-05-estimator-redesign-claude-review.md`).**
(1) B1 = option (a). The H ridge is rescaled to lambdaH ||H||_F^2 / (2 ||X||_F^2), with HALS and r_H changed
consistently (sections 1, 2, 4, 6). The old scaling was not strong regularization but inconsistent with the
normalized reconstruction term: at p = 3,000 it crushed H and made the reconstruction term almost irrelevant, an
artifact not carried into the redesign. lambdaH stays at 1e-6, not retuned prospectively. (2) The spec states that
D-Q1 to D-Q4 describe the failed estimator under the unscaled penalty and are not assumed to describe the
redesign. (3) No implementation before freeze. Allowed before freeze: the implementation brief, exact algebra and
gradient derivations, the caller and breakage inventory, unit-test design, and file and folder scaffolding with no
estimator logic. Not allowed: new objective or gradient code, small-p runs, timing or qualification, until the
independent Codex review is complete and the spec is frozen. Review item B3 (Q3 persistence and compensation)
remains open for a PI ruling.

**Open items:** B3 (above). One coordinator note for the reviewer and the PI:
- Q7's power depends on the S3 replicate count R, which the qualification prespec fixes. With R = 5 the one-sided
  bound needs mean D_AB above 2.13 SD / sqrt(5), about 0.95 SD of D_AB; a larger R for S3 alone may be
  worth its compute.
