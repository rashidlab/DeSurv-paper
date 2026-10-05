# Frozen specification: regularized DeSurv estimator, selection rule and gates (2026-10-04)

Step 2 (and the step 3 selection rule) of `2026-10-04-estimator-gate-ruling.md`. Committed before any step 1
triage output is read. Every later output of this experiment records the commit hash of this file. Changes after
the commit are allowed only before step 1 results are read, each as a dated amendment below with its reason;
after that, none.

PI decisions recorded 2026-10-04: margin 0.02 used as a decision threshold (Gate A), not a noninferiority test;
k fixed at 3; full independent tuning of every other hyperparameter in both arms.

## 1. Question

If the estimation procedure for the rank-3 DeSurv model is replaced by a convergent optimizer of an explicitly
regularized objective, does prediction remain within a practical threshold of the implemented procedure (Gate A),
and do the prespecified D1 biological findings reproduce (Gate B)?

**Why k = 3 is fixed.** k = 3 is the latent rank underlying the manuscript's D1 biological analysis. Rank is held
constant in both arms to isolate the effect of the estimator and optimizer specification. Letting k vary would
change the estimand as well as the estimator, and a D1 counterpart might not exist. A variable-k comparison is an
optional secondary analysis (section 9) and enters no decision.

## 2. Data and resampling

- Training data: the pooled TCGA + CPTAC training set exactly as built by `code/02_load_data.R` (n = 273), gene
  filter `ngene = 3000` as in the published BO, within-sample ranks (`method_trans_train = "rank"`).
- **Outer resampling:** 5-fold cross-validation, stratified by cohort x event, repeated over **4 partitions** with
  partition seeds **3001, 3002, 3003, 3004** (20 outer folds). The seeds were drawn for this experiment and have
  not been used before. Seeds 123-126 (item 4), 2001-2003 (feasibility check) and 1,000,001-1,000,400 (Track 3a)
  are excluded from every use here.
- Why 4 partitions: the decision uses the mean paired difference; its sampling variance has a floor from the
  test/train ratio (Nadeau and Bengio 2003: variance factor 1/J + 1/4, i.e. 0.30 at J = 20 and 0.27 at J = 50; SE
  multiplier sqrt(1/J + 1/4) = 0.55 and 0.52) that further repeats barely reduce, while compute scales with J.
- **Inner resampling:** inside each outer training set, 5-fold CV, stratified by cohort x event, fold seed
  derived as `outer_seed * 100 + outer_fold`. Identical inner folds for both arms.
- External cohorts (Dijk, Moffitt, PACA-AU x2, Puleo) and all previously examined partitions are development data.
  They enter no choice in this experiment. External cohorts are reported descriptively only.

## 3. Comparator arm (implemented estimator)

DeSurv 1.0.1 (`rashidlab/DeSurv@submission/pnas-2026`, sha `afb00d5`), the published procedure with k fixed:

- Inner BO via `DeSurv::desurv_cv_bayesopt_refine` with the published `BO_COMMON` (`engine = "warmstart"`,
  `tol = 1e-5`, `maxit = 4000`, `n_init = 50`, `n_iter = 100`, `candidate_pool = 4000`,
  `exploration_weight = 0.01`, `max_refinements = 0`, `lambdaW = lambdaH = 0`, `n_starts = 30` per evaluation),
  BO seed `outer_seed * 100 + outer_fold`.
- Search space: alpha in [0, 1] (linear), lambda in [1e-3, 1e3] (log10), nu in [0, 1], ntop in [50, 300]
  (integer). k = 3 fixed (the k bound removed).
- Selection: `overall_best$params` of the BO (the published extraction path, `published_rule.R`; with k fixed,
  `select_bo_k_by_cv_se` is a no-op). Full precision, never rounded.
- Final fit on the outer training set: the published `code/04` procedure (100 seed fits,
  `desurv_consensus_seed`, final fit).

## 4. Candidate arm (frozen estimator)

Primary design of `analysis/audit-verification/regularization_design_memo.md` section 2, frozen here. Notation as
Methods Eq. 1; X is p x n within-sample ranks.

**Objective.** W in S_+^k = {W >= 0, ||w_j||_2 = 1}, H >= 0, b in R^k (standardized coefficients);
sigma_j(W) = sqrt(||C X^T w_j||^2/(n-1) + eps^2), Zs(W) the centered, standardized scores, C the sample-centering
matrix, P the gene-centering matrix, sigma_bar^2 the mean gene variance of X:

F(W,H,b) = (1-alpha)[ ||X-WH||^2/(2||X||^2) + lambdaH ||H||^2/(2nk) ]
         + alpha[ -(2/n_e) l(Zs(W) b) + lambda( nu||b||_1 + (1-nu)/2 ||b||^2 ) + (rho sigma_bar^2/2) ||P W diag(1/sigma(W)) b||^2 ]
         + iota_{S+}(W) + iota_{>=0}(H),

l the Breslow partial log-likelihood, n_e the number of events. No zeta rescaling anywhere.

- Fixed constants: lambdaH = 1e-6. eps is one value per fit: 1e-6 x the median of sigma_j(W) (computed with
  eps = 0) over j = 1..k and over all 30 starts' projected initial W. Every start of a fit therefore minimizes the
  same F.
- Normalization: the unit-norm gauge on W is a constraint, enforced by the projection; no post-hoc rescaling.
  Reported beta_j = b_j / sigma_j(W) (the floored scale in the objective, never the raw sd).
- **Prediction vector (amendment 2).** v = W diag(1/sigma(W)) b and v~ = P v = v - mean(v) 1, both over the frozen
  training gene universe. The scored predictor is X_G^T v~_G, where G is the gene support chosen by the unchanged
  ntop rule (`desurv_get_top_genes` on W). Centering happens before support restriction. The fitted model stores
  the training gene universe, v~ and G; test and external cohorts are scored with these stored values and v~ is
  never recentered on the genes available in another dataset.
- **Invalid input** fails deterministically, with no fallback constant: n < 2, n_e = 0, any nonfinite entry of X,
  y or d, ||X|| = 0, or a pooled median sigma (above) <= 0.

**Optimizer.** Block prox-linear (PALM) with a convex b block solved to a certified tolerance:
1. H: one HALS sweep (Cichocki and Phan 2009), rows j = 1..k in order. No multiplicative updates.
2. W: W+ = P_S+(W - t grad_W Psi), column-wise closed form (v_+/||v_+||; if no positive entry, e_i for
   i = argmax of the unclipped v, first index on ties). Step t = min(2 t_prev, t_max) with t_max = 2 and
   t_prev = 1 at the first iteration, halved until
   Psi(W+) <= Psi(W) + <grad, W+-W> + ||W+-W||^2/(4t), so every accepted step decreases F by at least
   ||W+-W||^2/8. If t falls below t_min = 1e-20 the start is a numerical failure; no step is accepted at the floor.
3. b: convex block (elastic-net Cox plus the rho quadratic) solved by FISTA, warm-started from the current b
   (FISTA's internal acceleration is allowed). Stop when the KKT residual of the alpha-free b subproblem, evaluated
   by a routine separate from the solver, satisfies ||KKT_b||_inf <= 1e-10 (FISTA cap 1e5 iterations; reaching the
   cap is a numerical failure). Acceptance: b+ replaces b only if F(W, H, b+) <= F(W, H, b); otherwise b is kept
   and a rejected b step is recorded. This is an inexact block solve. The convergence theorem (Bolte, Sabach and
   Teboulle 2014) is stated for exact blocks; the implementation certifies the block's subgradient error at
   1e-10, five orders below the stopping tolerance, and makes no claim of mathematical exactness.

No outer inertial extrapolation in the frozen version.

**Stopping (amendment 2).** r = max(r_W, r_H, r_b), with G_W = grad_W Psi and G_H = grad_H Psi:
- r_W = ||R_W||_F, the KKT residual on {w_j >= 0, ||w_j|| = 1}: R_W[i,j] = G_W[i,j] - (w_j^T G_W[,j]) W[i,j] if
  W[i,j] > 0, and min(G_W[i,j], 0) if W[i,j] = 0 (the projection and HALS produce exact zeros).
- r_H = ||min(H, G_H / a_H)||_F / max(1, ||H||_F), elementwise min, a_H = (1-alpha)(1/||X||^2 + lambdaH/(n k))
  (the HALS row curvature): the relative size of a curvature-scaled projected step, defined at H = 0.
- r_b = ||KKT_b||_inf of the alpha-free b subproblem (the stricter reading; it bounds the F-scaled residual
  because alpha < 1).
Stop when r < 1e-5; maxit = 5000 outer iterations. The original residual
r_original = max( ||W - P_S+(W - G_W)||_F, ||min(H, G_H)||_F / ||H||_F, r_b ) is recorded alongside and binds
nothing. An F increase is F_new > F_old + 1e-13 max(1, |F_old|), checked after each block; the largest raw
increase is also recorded.

**Start and fit status (amendment 2).** Each start ends in one status: *converged* (r < 1e-5 within maxit, no F
increase), *nonconverged* (maxit reached, all values finite, no numerical failure), or *numerical failure*
(nonfinite value, step floor, FISTA cap, or any recorded F increase under the scale-aware definition above,
which the block acceptance rules exclude in exact arithmetic; differences within that floating-point tolerance
are not increases). A start ends at the first of: r < 1e-5, maxit, or a numerical failure. Every start keeps its status and full diagnostics (r, r_W, r_H, r_b,
r_original, final F, iterations, increases, rejected b steps, backtracking counts). A fit returns the
lowest-final-F *converged* start. If no start converged, the fit is nonconverged and carries the lowest-F
nonconverged start, flagged; if every start failed numerically, the fit failed. Numerical failures are never
scored.

**Initialization.** 30 starts per fit, drawn as the comparator's raw random draw (DeSurv 1.0.1
`desurv_alpha_warmstart.R`): for start s = 1..30 of a fit with seed S, `set.seed(S + s - 1)`, then
W0 = matrix(runif(p k, 0, max(X)), p, k) and H0 = matrix(runif(n k, 0, max(X)), k, n) in that order (H0 is drawn
so the stream matches and is then discarded). W0 is projected onto S_+^k; H by 20 HALS sweeps at fixed W; b = 0.
The comparator's `init()` short runs are not used (they run the zeta objective and choose by training C). Seeds:
inner fold f = 1..5 of every BO evaluation in outer fold (outer_seed, outer_fold) uses
S = outer_seed * 100 + outer_fold + 100 (f - 1) (the comparator's `seed_fold` rule with one hyperparameter set per
call). Start s always uses S + s - 1. Seed namespaces, disjoint by construction (nothing run with them before):
- inner-fold fits (both arms' comparable randomness): S = outer_seed * 100 + outer_fold + 100 (f - 1), the
  comparator's rule, range 300101-300834;
- candidate outer final fit: S = 10,000,000 + 1000 (outer_seed - 3000) + 100 outer_fold, range
  10,001,100-10,004,529;
- Gate B full-data inner-fold fits: S = 3000 + 100 (f - 1), range 3000-3429; Gate B full-data final fit:
  S = 20,000,000, range 20,000,000-20,000,029.
None overlaps another or the excluded seeds of section 2. Training C is never used to choose among starts.

**Tuning space (inner BO, k = 3 fixed).** logit(alpha) in [logit(1e-4), logit(0.9)]; log10 lambda in [-3, 3];
nu in [0, 0.95] (nu < 1 keeps the b block strongly convex); log10 rho in [-2, 3]; ntop in [50, 300] (integer).
alpha is searched on the logit scale because its meaning changes without zeta (the published alpha = 0.334
corresponds to a stated-objective weight near 0.007). Bounds are not changed after any outcome is seen;
boundary hits are reported.

**Inner BO.** Same BO engine, candidate pool, exploration weight and seed rule as the comparator. Budget
n_init = 60, n_iter = 120 (20% above the comparator's 150 evaluations, for the extra dimension). Each evaluation:
inner 5-fold CV C-index with 30 starts per fold, scored with the same ntop truncation and scoring code as the
comparator, applied to the stored centered vector v~. A configuration is *feasible* only if all five inner-fold
fits are converged. An infeasible configuration has no performance value: it is recorded with an explicit
infeasible flag and its score is NA, never a number. The engine fits the surrogate and takes the best only over
rows with status "ok" and a non-missing score (DeSurv 1.0.1 `desurv_cv_bayesopt.R:361`), so infeasible
evaluations enter neither the GP nor selection. They count toward the n_init + n_iter budget. No feasibility
surrogate is added (same engine as the comparator).

**Selection.** The feasible evaluated configuration with the highest mean inner-CV C (the same rule as the
comparator), full precision. If no configuration is feasible, the outer-fold final fit is failed.

**Final fit on the outer training set.** 30 starts at the selected configuration; the fit returned under the
status rule above. A converged final fit is scored. A nonconverged final fit is scored only as a flagged
diagnostic prediction and counts against the 90% convergence precondition (section 6). If any outer fold has a
failed final fit (no finite valid start, or no feasible configuration), Gate A is ineligible; no C-index is
imputed.

**Implementation preconditions** (before any outer-fold fit): gradient check against finite differences on every
block; F monotone on test problems; b block KKT verified by an independent evaluator; the three pilot-code limits
fixed (projection argmax, KKT stop, all residuals recorded). Debugging uses simulated data only; no held-out C is
computed on real data during development.

## 5. Scoring

For each outer fold and each arm: held-out C-index (Harrell, same code as the comparator's CV scoring) of the arm's
final fit applied to the outer test fold, with the arm's selected ntop truncation (the candidate scores with its
stored centered vector v~, section 4). Paired difference
dC_f = C_candidate,f - C_comparator,f over the 20 outer folds.

## 6. Gate A: estimator selection (statistical criteria only)

**Predictive-performance gate.** The candidate is eligible to replace the implemented estimator if the mean paired
out-of-sample C-index difference, mean_f dC_f, is at least -0.02 over the 20 prespecified fresh outer folds. The
corrected 95% confidence interval for the paired difference (Nadeau and Bengio 2003 corrected variance, factor
1/J + n_test/n_train, J = 20) is reported in full. This is a decision threshold for estimator selection and is not
interpreted as demonstrating statistical noninferiority. Results are described as "the point estimate met the
prespecified performance criterion" with the interval, never as "noninferior" or "equivalent".

**Convergence precondition.** The candidate exists to be a convergent estimator of a stated objective. It is
eligible only if at least 90% of its 20 outer-fold final fits are converged (section 4: r < 1e-5 under the
amendment-2 residual, within maxit, no F increase), no fit records an F increase, and no outer-fold final fit failed. If this fails, the candidate is ineligible regardless of dC, and the result is reported as an
implementation or specification failure, not as evidence about D1.

**Reported, not decision-binding:** training-test C gap per arm; start-to-start stability of W (median pairwise
loading correlation across the 30 starts, after Hungarian matching); stability of selected hyperparameters across
outer folds; inner-BO best-value traces; wall time.

**Gate A outcome** is computed by a script that reads only the outer-fold C values and convergence records, and is
committed (with this file's hash) before any Gate B script is run.

## 7. Gate B: biological replication (downstream of Gate A)

Runs only after the Gate A outcome is committed. It cannot change the Gate A outcome: if the candidate clears Gate
A and the biology differs, that is a scientific result to confront and report, not a reason to keep the implemented
estimator.

**Fit evaluated.** The candidate's full-data model: inner-BO procedure of section 4 run on all 273 training
samples (5-fold CV, fold and BO seed 3000; start seeds by the section 4 rule with outer_seed * 100 + outer_fold
replaced by 3000), selected configuration, 30 starts, lowest-F converged start. Gate B requires this fit to be
converged; otherwise the biology is reported as not evaluable.

**D1 counterpart.** Hungarian matching of the candidate's three W columns to the published fit's three
(`tar_fit_desurv_tcgacptac`) by Pearson correlation of loadings on the common gene set; D1's match is the D1
counterpart.

**Prespecified findings and metrics.** The D1 claims in the manuscript at DeSurv-paper-clean 999fe7e, listed in
Appendix A, each assessed for direction and nominal significance as in the paper. Plus, for the D1 counterpart:
loading Pearson r with published D1; overlap of the top-ntop gene sets; Spearman correlation of sample-level D1
scores on the training samples; reference-program correspondence, basal/classical relationship, PurIST and DeCAF
association, pathway enrichment, as in ruling step 1.

**Classification.**
- Recovers: loading r >= 0.80, sample-score Spearman >= 0.80, and every Appendix A claim reproduces in direction.
- Materially different: loading r < 0.50, or any Appendix A claim reverses direction.
- Otherwise partial. Each Appendix A claim is reported individually.

## 8. Decision table (replaces the ruling's pre-declared outcomes)

| Gate A | Gate B | Consequence |
|---|---|---|
| Clears | Recovers | Candidate becomes the primary estimator |
| Clears | Partial or materially different | Candidate becomes primary; narrow the D1 interpretation to what reproduces; hold the paper if the central conclusion disappears |
| Does not clear (dC < -0.02) | Recovers | Implemented estimator stays primary; candidate reported as a major robustness analysis |
| Does not clear | Partial or materially different | Implemented estimator stays primary; report which D1 findings depend on the estimator and narrow claims accordingly |
| Clears | Not evaluable (the frozen procedure yields no converged full-data fit) | Candidate is not eligible to become primary. Report the full-data convergence failure; do not interpret biology. No refitting outside the frozen procedure |
| Ineligible (convergence precondition fails) | not run | Report as specification failure; another principled route needed (simpler supervised factor model or different regularization structure) |

## 9. Optional secondary analysis (no decision weight)

The same nested design with k tuned in both arms (k in [2, 12], rev-3 1-SE rule), if compute allows, answering
whether the candidate is competitive as a general DeSurv estimator. Run only after Gate A is committed.

## 10. Safeguards

- Selection and Gate A scripts compute no biology outputs. The Gate B script refuses to run unless the Gate A
  outcome file exists in a commit.
- No script reads step 1 triage outputs. Triage outputs are not opened until this file is committed.
- All seeds, package shas and this file's hash are written into every result object.

## Appendix A: D1 claims (replication targets)

Extracted from the manuscript at DeSurv-paper-clean 999fe7e (results, discussion, SI). Every analysis is rerun with
the candidate's D1 counterpart in place of D1, using the same code, cohorts and covariates as the paper. "Reproduces"
means the direction criterion below holds; nominal significance is reported for each claim but is not part of the
classification. "Weak" means |rho| < 0.5 throughout.

| ID | Claim (source) | Rerun | Reproduces if |
|---|---|---|---|
| A1 | D1 loads on classical-tumor and restCAF-like reference programs (results:437, Fig. 2a) | Spearman of top-50 loadings vs the 45-program reference universe | Classical-tumor and restCAF-like programs are the two highest-rho programs, each rho > 0.2 |
| A2 | Classical component supported by GATA6 ISH, COMPASS n = 33 (results:437) | Spearman of D1-counterpart score vs GATA6 ISH | rho > 0 |
| A3 | D1 has the smallest reconstruction share but the largest survival contribution among the three factors (results:443) | Shapley reconstruction share and drop-one delta-loglik, as `fig_variation_explained` | D1 counterpart has the smallest share and the largest delta-loglik |
| A4 | D1 not captured by matched-rank standard NMF (results:446) | Max Spearman with `fit_std_desurvk` factors; projection R^2 onto the NMF basis | Max rho weak and R^2 < 0.5 |
| A5 | D1 not recovered by the matched-budget alpha = 0 control (results:449) | Best-match Spearman of `tar_fit_desurv_a0k3` factors with the D1 counterpart | Weak |
| A6 | D1 is the only factor with a significant pooled external association, protective in all five cohorts (results:463) | Per-factor dataset-stratified Cox on the five validation cohorts | D1 counterpart HR < 1 pooled and in all five cohorts, and it has the smallest pooled P of the three factors |
| A8 | D1 retains prognostic value after PurIST + DeCAF adjustment (results:566; SI 1730-1767) | `d1_variance_partition_stats` procedure | Adjusted HR < 1 and LR chi-square for adding D1 > 0 with concordance not decreased |
| A9 | D1 weakly correlated with purity in PACA-AU array and seq (results:569) | Spearman with purity, `desurv_clinical_adjustment` procedure | Weak in both |
| A10 | D1 prognostic in treated metastatic O'Kane/COMPASS, protective in all three arms, no arm interaction (results:575; SI 2472-2478). The paper frames this as a test of D1's tumor-associated direction only (epithelial-enriched tissue); kept as a claim about the projected D1 score | Per-arm and arm-stratified Cox, interaction test | Stratified HR < 1 and per-arm HR < 1 in all three arms |
| A11 | Independent supervised methods (superpc, Cox-PLS, lasso Cox, tuned versions) recover D1, while the leading unsupervised PC does not (results:612, 615; SI 1500-1560) | `supervised_recovery_stats`, `desurv_vs_supervised_tuned` with the D1 counterpart | Every supervised comparator's absolute r with the D1 counterpart > 0.5 and larger than with either other factor; leading PC weak |

Not targets: the composite D1 + D3 predictor HR (A7; reported, not classified, since it is not a D1-only claim), the
top-270 gene list (A12; covered by the top-ntop overlap metric in section 7), and the bulk-not-cellular scope caveat
(A13; not testable).

## Amendments

- 2026-10-04 (before any step 1 output was read; codex review of bca989d): reported beta uses sigma_j(W), not
  sd(Z_j), to match the objective; Nadeau-Bengio variance factor and SE multiplier stated separately (numbers
  unchanged); A11 criterion reworded to remove pipe characters that split the table cell. No decision changed.
- 2026-10-04, amendment 2 (before any step 1 output was read, on the PI's instruction; codex scientific review of
  the candidate plan, `2026-10-04-candidate-optimizer-codex-review.md`, findings 1-6). (1) The b block is an
  inexact convex solve with an independent KKT certificate at 1e-10 and an F-nonincrease acceptance check; "exact"
  removed. (2) W step capped at t_max = 2, giving decrease >= ||W+-W||^2/8; no step accepted at the floor.
  (3) The binding stopping residual is replaced: r_W is the sphere/nonnegativity KKT residual (the unit-step
  projection residual need not vanish at a stationary point) and r_H is curvature-scaled with a zero convention;
  the original r is recorded, not binding. (4) The candidate's scored predictor uses the gene-centered vector
  v~ = P W diag(1/sigma) b, centered before ntop truncation over the frozen training gene universe and never
  recentered on another dataset's genes, because truncation otherwise reintroduces the uniform direction that
  the objective leaves unpenalized; the support rule G is unchanged. (5) Initializer (raw comparator draw, full
  seed schedule) and eps (pooled median over 30 starts) written into the spec; deterministic invalid-input
  failure. (6) Start statuses and failure propagation through starts, BO evaluations (infeasible = 0.5, never
  selected), outer final fits (Gate A ineligible on any failed fit) and Gate B (requires a converged fit).
- 2026-10-04, amendment 2 rulings (PI, same day, before any step 1 output was read): decision-table row for Gate A
  clears with Gate B not evaluable (candidate not eligible; technical failure, no biology); infeasible BO
  configurations carry NA plus an infeasible flag and are excluded from the surrogate and selection, replacing
  the 0.5 value; seeds given as disjoint namespaces (the earlier +500 final-fit offset overlapped inner-fold seeds
  of other outer folds in 102 cases; no seed had been used); the F-increase failure uses the existing
  scale-aware tolerance.
- **Amendment window closed** after amendment 2 (PI, 2026-10-04). No further amendments.
