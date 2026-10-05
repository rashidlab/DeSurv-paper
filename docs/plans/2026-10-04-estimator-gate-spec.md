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

- Fixed constants: lambdaH = 1e-6; eps = 1e-6 x median_j sigma_j at initialization.
- Normalization: the unit-norm gauge on W is a constraint, enforced by the projection; no post-hoc rescaling.
  Reported beta_j = b_j / sigma_j(W) (the floored scale in the objective, never the raw sd).

**Optimizer.** Block prox-linear (PALM) with an exact convex b block:
1. H: one HALS sweep (Cichocki and Phan 2009). No multiplicative updates.
2. W: W+ = P_S+(W - t grad_W Psi), column-wise closed form (v_+/||v_+||; if no positive entry, e_i for
   i = argmax of the unclipped v). Backtracking from 2 t_prev by halving until
   Psi(W+) <= Psi(W) + <grad, W+-W> + ||W+-W||^2/(4t).
3. b: exact minimization (elastic-net Cox plus the rho quadratic) by FISTA, warm-started, stopped on the KKT
   residual ||KKT_b||_inf < 1e-10 (not on successive coefficient change).

No inertial extrapolation in the frozen version.

**Stopping.** r = max( ||W - P_S+(W - grad_W Psi)||_F, ||min(H, grad_H Psi)||_F / ||H||_F, ||KKT_b||_inf ).
Stop when r < 1e-5; maxit = 5000 outer iterations. Every fit records r, final F, iteration count, whether r
reached tolerance, and the number of F increases (must be 0).

**Initialization.** 30 starts per fit: W from the same `DeSurv` initializer as the comparator, projected onto
S_+^k; H by 20 HALS sweeps at fixed W; b = 0. The fit returned is the start with the **lowest final F**. Training C
is never used to choose among starts.

**Tuning space (inner BO, k = 3 fixed).** logit(alpha) in [logit(1e-4), logit(0.9)]; log10 lambda in [-3, 3];
nu in [0, 0.95] (nu < 1 keeps the b block strongly convex); log10 rho in [-2, 3]; ntop in [50, 300] (integer).
alpha is searched on the logit scale because its meaning changes without zeta (the published alpha = 0.334
corresponds to a stated-objective weight near 0.007).

**Inner BO.** Same BO engine, candidate pool, exploration weight and seed rule as the comparator. Budget
n_init = 60, n_iter = 120 (20% above the comparator's 150 evaluations, for the extra dimension). Each evaluation:
inner 5-fold CV C-index with 30 starts per fold, scored with the same ntop truncation and scoring code as the
comparator.

**Selection.** The evaluated configuration with the highest mean inner-CV C (the same rule as the comparator),
full precision.

**Final fit on the outer training set.** 30 starts at the selected configuration, lowest F returned.

**Implementation preconditions** (before any outer-fold fit): gradient check against finite differences on every
block; F monotone on test problems; b block KKT verified; the three pilot-code limits fixed (projection argmax,
KKT stop, all three residuals recorded). Debugging uses simulated data only; no held-out C is computed during
development.

## 5. Scoring

For each outer fold and each arm: held-out C-index (Harrell, same code as the comparator's CV scoring) of the arm's
final fit applied to the outer test fold, with the arm's selected ntop truncation. Paired difference
dC_f = C_candidate,f - C_comparator,f over the 20 outer folds.

## 6. Gate A: estimator selection (statistical criteria only)

**Predictive-performance gate.** The candidate is eligible to replace the implemented estimator if the mean paired
out-of-sample C-index difference, mean_f dC_f, is at least -0.02 over the 20 prespecified fresh outer folds. The
corrected 95% confidence interval for the paired difference (Nadeau and Bengio 2003 corrected variance, factor
1/J + n_test/n_train, J = 20) is reported in full. This is a decision threshold for estimator selection and is not
interpreted as demonstrating statistical noninferiority. Results are described as "the point estimate met the
prespecified performance criterion" with the interval, never as "noninferior" or "equivalent".

**Convergence precondition.** The candidate exists to be a convergent estimator of a stated objective. It is
eligible only if at least 90% of its 20 outer-fold final fits reach r < 1e-5 within maxit and no fit records an F
increase. If this fails, the candidate is ineligible regardless of dC, and the result is reported as an
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
samples (5-fold CV, fold and BO seed 3000), selected configuration, 30 starts, lowest F.

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
