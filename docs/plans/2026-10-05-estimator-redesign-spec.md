# Estimator-redesign specification: rank-3 DeSurv in the theta parameterization (2026-10-05)

**Status: DRAFT (revision 2), not frozen.** Written under PI ruling 8 (`docs/HANDOFF-2026-10-05-candidate-optimizer.md`,
"PI ruling 2026-10-05, eighth"). Revision 2 applies PI ruling 9 (section 10), which approved the redesign in principle
and ruled on the four open items and draft choices D2-D4. No implementation, no real data and no compute until the PI
freezes this file. Once frozen, every later output records this file's commit hash. Amendments then follow the rule in
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
- **Data:** X is the p x n matrix of within-sample ranks (average ties) of the fit's training columns, p = 3,000 genes
  (the frozen `ngene = 3000` filter), as in `palm_rank` / `palm_fold_data`. C is the n x n sample-centering matrix.
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
- **H normalization:** H >= 0, with no norm constraint and a ridge lambdaH = 1e-6 (unchanged). H is not rescaled.
- **Scores:** Z_c(W) = C X^T W (n x k), the centered scores of the unit-norm loadings. No division by any
  W-dependent score SD anywhere in the objective.

## 2. Objective (theta parameterization)

The optimized survival coordinate is theta in R^k, with eta = Z_c(W) theta / s0. The reported raw-score Cox
coefficient is beta = theta / s0.

F(W, H, theta) = (1 - alpha) [ ||X - W H||_F^2 / (2 ||X||_F^2) + lambdaH ||H||_F^2 / (2 n k) ]
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

This spec makes no convergence-theorem claim. Before implementation, the PALM convergence assumptions are checked
explicitly, one by one, for the block updates of section 4. The results go in a written assumption memo, committed
before any qualification run. The memo either proves each condition or states which one fails. In the latter case
the operational definition in section 6 (monotone decrease plus the r < 1e-5 KKT certificate) is the only claim
made. The checklist:

| # | Condition (Bolte, Sabach and Teboulle 2014 framework) | To check for this objective |
|---|---|---|
| A1 | Nonsmooth parts are proper and lower semicontinuous, with computable proximal maps | iota_{S+} (nonconvex set, projection possibly multivalued, the argmax tie rule), iota_{>=0}, the L1 term on theta |
| A2 | Smooth coupling term has block gradients that are Lipschitz, with moduli bounded above and below on the iterates' bounded set | The Cox term is C^2, so it is locally Lipschitz. Use the sublevel-set bound on ||theta|| and the curvature bound below |
| A3 | Bounded iterates | Monotone decrease plus bounded sublevel sets; holds if every accepted block step decreases F |
| A4 | KL property of F | The Breslow log-likelihood involves log-sum-exp, so definability is in the log-exp o-minimal structure, not the semialgebraic one; cite the precise result |
| A5 | Each block update is the prox-linear step of the theorem (or is covered by a stated extension) | H uses HALS (exact row minimization, not a prox-linear step). W uses projected gradient with backtracking (prox-linear; needs a step lower bound). theta uses an inexact FISTA solve with a KKT certificate and an F-acceptance check (inexact block minimization). Each needs a cited extension or a proof, or the gap is recorded |

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
| lambdaH | H ridge | 1e-6 fixed | |
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
variable. This reuses existing block steps and is not new solver work. Newton-CG, Schur-complement and outer
acceleration are excluded (ruling 8.1).

1. **H:** one HALS sweep, rows j = 1..k in order.
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
- r_H = ||min(H, G_H / a_H)||_F / max(1, ||H||_F), with a_H = (1 - alpha)(1/||X||^2 + lambdaH/(n k)).
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
- *Pathology criteria* (binding) are tied to a specific failure mode: the D-Q4 collapse, its opposite bias, or the F2
  overfit. Each one states the reason for its number.
- *Validity measures* (reported, not binding) show whether the simulation was informative: the oracle C gap, and
  score and loading recovery.

**Scenarios** (all built on the repo generator through `analysis/estimator-gate/candidate/sim.R`, scenario R0, with
within-sample ranks and integer-rounded tied times):
- **S1, signal:** R0 with generator beta = (2, 1, 0).
- **S2, null:** R0 with beta = (0, 0, 0).
- **S3, competing directions:** a genuinely predictive low-variance factor A competes with a higher-variance
  nonpredictive factor B. The generator's survival linear predictor is computed on standardized scores
  (`simulate_survival.R:189-200`), so scaling the rows of the true H changes a factor's variance in X without changing
  its predictive strength. S3 is R0 with generator beta = (2, 0, 0), where factor 1 is A, and true H rows multiplied
  by c = (c_A, c_B, 1) before X is formed. c_A and c_B are chosen so that A's true score SD is at most 0.25 times B's.
  The generator has no per-factor variance argument, so the qualification prespec specifies a wrapper and its own
  check that the realized SD ratio is at most 0.25. S3 is where direct penalization of theta could attach the
  survival effect to B, the bias opposite to the old estimator's.

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
| Q2 | Pathology | No sigma -> 0 | For every carrying column of every converged start, log(sigma_j(end) / sigma_j(sweep 50)) >= -1. Reason: the D-Q4 collapse lost 6 to 9 in log sigma per carrying column, so a loss of 1 (a factor of e) is under one sixth of the smallest observed collapse |
| Q3 | Pathology | No preserved-score collapse | Using the D-Q4 decomposition (`analysis/estimator-gate/candidate/tests/DQ4_step_decomposition.R`; log sigma = log r + log g and the score angle depend on W only, so the adaptation changes only the objective-term accounting, b to theta), no carrying column has cumulative change in log g <= -1 while its net standardized-score angle is <= 0.01 rad, over the whole trajectory. Reason: the D-Q4 signature was a log g loss of 6 or more with a net angle of at most 0.0023 rad, so both limits sit well inside the pathology |
| Q4 | Math (a), Pathology (b) | Bounded curvature | (a) The largest eigenvalue of the W-block Hessian on the active face, by power iteration on HVPs at sweeps 1, 10, 100 and the end, never exceeds the section 2 bound L_W beyond floating-point tolerance. (b) Its largest recorded value per start is within a factor of 100 of its smallest. Reason: the pathology was 1e6-level growth tracking sigma^-4, so 1e2 is four orders below it. The slope of log lambda_max on log min_j sigma_j is reported |
| Q5 | Pathology | Multi-start stability | At the centre and tuned configurations, among converged starts: (i) the median Hungarian-matched loading Pearson r of carrying columns with the returned start is >= 0.80, the Gate B "recovers" threshold, so that a Gate B classification cannot depend on which start was returned; (ii) the median across starts of the test-C difference from the returned start is within 0.01, the F2 margin below. The share of starts within 1e-4 relative F of the best is reported |
| Q6 | Pathology (binding F2 safeguard) | Convergence does not destroy generalization | In S1 and S3, at the tuned configuration, for every checkpoint t in T: mean over replicates of the paired difference C_conv - C_t (test C of the converged returned start minus test C of its checkpoint-t iterate) >= -0.01, with each iterate scored by v~ on its own ntop support. Reason (ruling 9): a properly converged estimator may not give up more than one C point to its own finite trajectory. This is deliberately stricter than Gate A's -0.02 against the published estimator. S2 is reported |
| Q7 | Pathology (opposite bias) | Effect lands on the predictive direction | In S3, see open item 9a |
| V1 | Validity | Oracle gap | Test C of the fit vs test C of the true linear predictor, per scenario and size. Reported only, because it depends on signal strength and censoring (ruling 9, item 1) |
| V2 | Validity | Predictive-score recovery | Pearson and Spearman correlation, on the test set, of the fitted predictor X_G^T v~_G with the true linear predictor, in S1 and S3; in S3 also its correlation with factor B's true score |
| V3 | Validity | Loading and program recovery | In S1 and S3, where the generative W is identifiable: Hungarian-matched cosine of fitted vs true W columns, the matched cosine of the carrying column with the true predictive factor, and top-ntop marker overlap with the true marker set. Reported because Gate B concerns loadings, and D-Q4 showed that loadings can move while prediction barely does |

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
   and check the PALM assumptions explicitly before implementation (section 2 checklist).
9. **Recovery:** qualification assesses both predictive-score recovery and loading/program recovery in S1 and S3.

**Open items for the PI (revision 2):**

- **9a. Q7, opposite-bias criterion in S3: binding or reported?** Proposed binding form, which uses comparisons
  rather than a numeric cutoff: averaged over replicates, the fitted predictor's test correlation with the true
  linear predictor (A) exceeds its correlation with B's true score, and the column with the largest |theta_j| is
  Hungarian-matched to A rather than B in a majority of replicates. The alternative is to report V2 and V3 in S3
  without a binding rule. Ruling 9 says S3 must test this bias but does not say whether it gates.
- **9b. Pathology limits in Q2, Q3, Q4b and Q5.** Each now carries a stated derivation: the D-Q4 magnitudes for Q2
  to Q4b, and the Gate B 0.80 and F2 0.01 margins for Q5. Confirm that these count as tied to the pathology, or
  demote any of them to reported.
- **9c. S3 variance ratio.** A true score SD of at most 0.25 times B's (variance ratio at most 1/16) is a DGP
  choice, not a threshold. Confirm, or set a second ratio to run as well.
- **9d. Assumption memo timing.** The section 2 checklist is completed before implementation. If A5 cannot be met
  as written (HALS and the inexact theta solve are not prox-linear steps), choose between recording the gap and
  relying on the operational certificate, or replacing the H step with a prox-linear step. The second changes the
  optimizer and needs a ruling before freeze.
