# Qualification prespecification for the estimator-redesign spec (2026-10-05)

**Status: FROZEN 2026-10-05 with spec revision 4 (PI ruling 14).** It operationalizes section 9 of `docs/plans/2026-10-05-estimator-redesign-spec.md`
(revision 4, frozen under ruling 14; the Codex review runs in parallel as quality control). It freezes together with that spec, or immediately
after it. If the review or a PI ruling changes a section 9 criterion, this file follows the change before freezing.
It never changes a criterion on its own authority. Once frozen, it must be committed before the timing-only run
(section 8) and before any qualification fit, and every qualification output records the hashes of both files.

**What this file fixes (spec section 9, "A qualification prespec fixes ...").** The replicate counts, including R for
S3, which sets Q7's power. The configuration grid and the reduced tuning run. The S3 data wrapper, which sets
per-factor H scaling, and the check of the realized SD ratio. The seed namespaces. What each fit records. For every
binding criterion, which fits it reads and how they are aggregated. The verifier controls. The order of execution and
the stopping rules. Where the spec text is ambiguous, this file states its reading and lists it in section 10 for
the PI.

**Work already done for this file (generator only).** `analysis/estimator-gate/qualification/sim_qual.R` (wrapper)
and `s3_calibrate.R` (calibration). Output: `s3_calibration.txt` and `.rds`. No estimator code was sourced and no fit
was run. Calibration is deterministic: rerunning it after freeze gives the same c. Seeds 49,000,001-49,000,020
(calibration) and 49,000,101-49,000,200 (validation) are now used, and only by the generator.

## 1. Data construction (all scenarios)

- **One generator draw per replicate.** `sim_qual(seed, scenario, cvec)` calls the repo generator once
  (`simulate_desurv_scenario`, scenario R0, G = 3,000, K = 3, 150 markers per factor, 500 background genes, as in
  `candidate/sim.R`) with N = 218 + 1,000. Columns 1-218 are the training set. Columns 219-1,218 are the test set.
  Training and test therefore share W, marker sets, gene baselines and score standardization. A separate generator
  call would redraw all of these (`simulate_data.R:91-111`, `simulate_X.R:51`), so the test set would have a
  different truth. The generator standardizes the true scores over all N columns. This changes only the scale of
  the true predictor.
- **Ranks and ties as in the pipeline.** Within-sample ranks (average ties) per column, and survival times rounded
  up to integers (tie unit 1).
- **Sizes.** n = 218 is the full training set (outer-training size). n = 175 is columns 1-175 of the same draw
  (inner-training size). Both sizes are scored on the same 1,000 test columns.
- **True scores.** "True A" and "true B" (Q7, V2, V3) are columns 1 and 2 of the generator's standardized scores
  `scores = std(X_latent^T Wtilde)`, taken on the test columns, where X_latent = W H^T is noise-free and Wtilde is W
  restricted to each factor's markers (`simulate_survival.R:175-197`). The true linear predictor is the generator's
  `linpred` on the test columns.
- **One code path.** S1, S2, S3 and S3b all go through `sim_qual`. For S1 and S2, cvec = (1, 1, 1).

## 2. Scenarios and replicate counts

| Scenario | Generator beta | cvec | Replicates R | Used for |
|---|---|---|---|---|
| S1, signal | (2, 1, 0) | (1, 1, 1) | 5 for the grid; **3** for tuning (ruling 15) | grid, tuning, Q1-Q6 |
| S2, null | (0, 0, 0) | (1, 1, 1) | 5 for the grid; **3** for tuning (ruling 15) | grid, tuning, Q1-Q5 (Q6 reports S2 only) |
| S3, competing directions (binding) | (2, 0, 0) | (0.1911, 1, 1) | **10** for tuning (ruling 15; was 20); replicates 1-5 also run the grid | grid, tuning, Q1-Q7 |
| S3b, sensitivity | (2, 0, 0) | (0.474, 1, 1) | **dropped from qualification (ruling 15)**; may be run later for the paper | none |

**Why R = 20 for S3 (superseded by ruling 15: R = 10; the table below gives its power).** Q7 passes if mean D_AB - t_{0.95, R-1} SD / sqrt(R) > 0. That is, the replicate-level
standardized effect delta = mean / SD must exceed t_{0.95, R-1} / sqrt(R). The probability of passing when the true
effect is delta (noncentral t) is:

| R | delta needed | Pass probability at delta = 0.5 | 0.75 | 1.0 |
|---|---|---|---|---|
| 5 | 0.95 | 0.24 | 0.40 | 0.58 |
| 10 | 0.58 | 0.43 | 0.71 | 0.90 |
| 20 | 0.39 | 0.70 | 0.94 | 1.00 |
| 30 | 0.31 | 0.85 | 0.99 | 1.00 |

When the true mean D_AB is 0, the false-pass rate is 5% at every R by construction. At R = 5, an estimator whose
predictor lands on A in most replicates still fails Q7 about half the time when delta = 1. So the 5-replicate count
of the spec draft is kept for S1 and S2 and raised to 20 for S3 alone, where the cost is the 15 extra tuning runs
(section 9). R = 20 is fixed before any timing or fit output is seen.

## 3. The S3 wrapper

### 3.1 Mechanism

The generator's H is N x K (`simulate_H.R:17-20`), so "H rows" in the spec means H columns here: one per factor.
`gen_scaled` temporarily rebinds `simulate_H` in the global environment to `sweep(simulate_H(...), 2, cvec, "*")`,
calls `simulate_desurv_scenario` once, and restores the original on exit. Because the multiplication draws no
random numbers, the RNG order is unchanged: W and marker sets, then H, then survival, then X noise. The scaled H
enters both X_latent for survival and `simulate_X`. Scaling H column j changes factor j's variance in X without
changing its predictive strength, because the generator's linear predictor uses standardized scores.

**Wrapper self-check (every S3 and S3b replicate, before any fit).** The override was called exactly once. The draw
with cvec and the draw with (1, 1, 1) at the same seed have identical W and marker sets. H(cvec) equals
H(1, 1, 1) times cvec column-wise to 1e-15. Any failure stops the run as a harness defect.

### 3.2 Which SD ratio is the target (reading of spec section 9, item P1 in section 10)

The spec asks for "A's true score SD 0.25 times B's". Three spaces are candidates, and they disagree.
Median over 100 validation draws at the S3 c:

| Space | SD ratio A/B | corr(A, B) |
|---|---|---|
| Ranked training X, unit-norm marker-masked true loadings (**binding target**) | 0.252 | -0.41 (test set) |
| Generator scale: noise-free X_latent, marker-masked loadings | 0.193 | 0.05 (generator scores, test set) |
| Ranked training X, unit-norm full true loadings | 0.370 | not computed in validation; 0.72 at c_A = 0.25 in the probe |

The binding target is the first: the ranked data the estimator fits. S3 tests whether the theta penalty makes a
high-variance score direction cheaper per unit of effect, and that cost is set by score variance in the space the
penalty acts on. This target is also the stricter choice. In generator units, A ends up weaker (0.19) than the
literal reading would give (c = (0.25, 1, 1), ranked ratio about 0.30).

Full loadings cannot carry the target. Every true column loads the same 500 background genes. At c_B = 1, the
ranked full-loading ratio does not fall below about 0.31 even at c_A = 0.01. Reaching 0.25 needs c_B near 4, where
the full-loading A and B scores correlate at 0.96 (generator-only probe, five seeds). The full-loading ratio is
reported as a validity measure.

**"Uncorrelated by construction" holds on the generator scale only.** True A and B (generator scores) correlate at
median 0.05 (range -0.03 to 0.19). On ranked X, the marker-masked scores correlate at about -0.41. Within-sample
ranking is zero-sum: when B's 150 markers rise, every other gene's rank falls. Q7 uses the generator-scale scores
and the joint regression, so it remains unambiguous. Both correlations are reported per replicate.

### 3.3 Calibration (done, deterministic)

c_B = c_3 = 1. c_A was set by bisection on log c_A over [0.01, 1] so that the median, over calibration seeds
49,000,001-49,000,020, of the realized ranked masked SD ratio on the 218 training columns equals the target within
0.002. Fixed values:

- **S3:** c = (0.1911, 1, 1), calibration median 0.2517 after 6 iterations. Validation over seeds
  49,000,101-49,000,200: median 0.252, 5th-95th percentile 0.207-0.318, range 0.198-0.351.
- **S3b:** c = (0.474, 1, 1), calibration median 0.5008. Validation: median 0.504, 5th-95th 0.400-0.669, range
  0.363-0.721.
- Training event rate in validation: median 0.61 in both scenarios.

c is not recalibrated after any fit.

**Frozen mechanics (PI, 2026-10-05).** The wrapper and calibration are fixed as committed. Changes require a PI
ruling. Git blob hashes (first 12 characters): `sim_qual.R` 2379bb00de3f, `s3_calibrate.R` 0f2017a0ac58,
`s3_calibration.rds` d8e7d1f9f534. Generator files: `simulate_data.R` 301920a23600, `simulate_W.R` f21b1af5d47a,
`simulate_H.R` 77a66ef5fb8e, `simulate_X.R` d0e55aede746, `simulate_survival.R` ac52c17b1e01, `scenario_defaults.R`
dca30349342d. Constants: c = (0.1911, 1, 1) for S3 and (0.474, 1, 1) for S3b, N = 218 + 1,000, the section 4 seeds.
Every qualification run checks these hashes before it generates data, and stops on a mismatch.

### 3.4 Realized-ratio check (binding harness check, before any fit)

For S3, every replicate's realized ranked masked ratio must lie in [0.18, 0.37], which contains the full validation
range. The median over the 10 S3 tuning replicates (ruling 15) must lie in [0.22, 0.29]: the validation SD is about 0.034,
so the median of 10 has an SD near 0.013, and the window is about plus or minus 2.7 SD. A failure means the wrapper or the generator
does not behave as calibrated. It stops the run before any fit and is reported. It is never fixed by replacing a
replicate. S3b is not part of qualification (ruling 15); if run later, its ratios are reported without a check.

## 4. Seeds (disjoint namespaces, none used before)

sc = 1, 2, 3, 4 for S1, S2, S3, S3b. r = replicate. g = fit index 1-12 (section 5 table). f = inner fold. Start s of a fit with
seed S uses `set.seed(S + s - 1)`, s = 1..30, as in spec section 5.

| Use | Seed |
|---|---|
| Dataset | 40,000,000 + 100,000 sc + r |
| Grid fit | 41,000,000 + 100,000 sc + 1,000 r + 50 (g - 1) |
| Tuning run: inner-fold partition (5 folds, stratified by event status) and BO seed | 42,000,000 + 100,000 sc + r |
| Tuning run: inner-fold fit, shared by every BO evaluation (as in the gate spec) | 43,000,000 + 100,000 sc + 1,000 r + 100 (f - 1) |
| Tuning run: final fit at the selected configuration, n = 218 | 44,000,000 + 100,000 sc + 1,000 r |
| S3 calibration / validation (used, generator only) | 49,000,001-49,000,020 / 49,000,101-49,000,200 |
| Timing-only run (section 8) | 49,100,001 |

Disjointness: within a namespace, start seeds advance by at most 29 under spacings of 50, 100 and 1,000, and the
namespaces occupy separate millions. None overlaps the gate-spec real-data namespaces (300,101-300,834,
10,001,100-10,004,529, 3,000-3,429, 20,000,000-20,000,029) or the earlier simulation seeds (11-23, 123, 600, 700,
1,000-5,000 ranges, 31,337, 424,242-424,243, 700,000 ranges, 20,261,005). A repo-wide search on 2026-10-05 found no
use of 40,000,000-49,999,999.

## 5. Configurations

**Grid (fixed, inside the spec section 3 box): 11 grid configurations at n = 175, plus 1 separate reference fit at
n = 218, for 12 fits per replicate.** The *centre configuration* is defined by its values: alpha = 0.05,
lambda = 1, rho = 1, nu = 0.5, ntop = 150. It is one of the 11 grid configurations (at n = 175). The reference fit
is the same centre configuration on all 218 training columns. g is the fit index used by the harness and the seed
rule (section 4). Every row below is one fit. nu = 0.5 and ntop = 150 in every row.

| g | alpha | lambda | rho | n | Role |
|---|---|---|---|---|---|
| 1 | 1e-3 | 1e-2 | 1 | 175 | grid |
| 2 | 1e-3 | 1 | 1 | 175 | grid |
| 3 | 1e-3 | 1e2 | 1 | 175 | grid |
| 4 | 0.05 | 1e-2 | 1 | 175 | grid |
| 5 | 0.05 | 1 | 1 | 175 | grid; centre configuration |
| 6 | 0.05 | 1e2 | 1 | 175 | grid |
| 7 | 0.5 | 1e-2 | 1 | 175 | grid |
| 8 | 0.5 | 1 | 1 | 175 | grid |
| 9 | 0.5 | 1e2 | 1 | 175 | grid |
| 10 | 0.05 | 1 | 1e-2 | 175 | grid |
| 11 | 0.05 | 1 | 1e3 | 175 | grid |
| 12 | 0.05 | 1 | 1 | 218 | reference fit: centre configuration at n = 218 |

Grid fits and reference fits run on replicates 1-5 of S1, S2 and S3: 3 scenarios x 5 replicates x 12 fits = 180
fits of 30 starts.

**Tuning run (one per scenario and replicate: 3 + 3 + 10 = 16 under ruling 15).** This is the candidate procedure of
spec section 7 on the 218 training columns, with three changes: the reduced budget n_init = 20, n_iter = 40; inner
folds stratified by event status only (there is no cohort); and the seeds of section 4. Everything else is
unchanged: the same BO engine, candidate pool and exploration weight; the five-dimensional section 3 box (alpha,
lambda, nu, rho, ntop); 5-fold inner CV with 30 starts per fold; per-inner-fold s0; feasibility only when all five
inner-fold fits converge; infeasible = NA plus flag, excluded from the surrogate and from selection; selection = the
feasible configuration with the highest mean inner-CV C, at full precision. The final fit at the selected
configuration uses all 218 training columns, 30 starts and its own s0, and is scored with v~ on G. If no
configuration is feasible, the tuning run has no tuned configuration. That is recorded as a Q1 failure, because the
procedure cannot produce an estimator.

## 6. What every fit records (enables Q2-Q6; binds nothing in Gate A)

Per start: status and reason; sweeps run; final F, r_W, r_H, r_theta; every F increase (there should be none). Every
10th sweep, per column j: sigma_j, r_j = ||P w_j||, g_j = sigma_j / r_j, and theta_j. (W, theta) at sweeps 10, 50,
100, 500, 3,750, 4,375 and the endpoint. For a start that converges before one of these sweeps, the endpoint stands in
(spec section 9). Per fit: s0, lambda_1(X C X^T), the returned start, and the Q4 power-iteration records (section 7,
Q4). Test C of every converged start and of every checkpoint iterate of the returned start, each scored by its own v~
on its own ntop support. The test-C scorer is the gate-spec section 5 scorer.

## 7. Binding criteria: scope and aggregation

Unit throughout: a *fit* is one 30-start fit, and a *cell* is (scenario, configuration, size).

- **Q0 (math, before any qualification fit).** These checks run on small simulated sets (seeds from the 49,200,000
  range), with the floating-point tolerances of the gate spec section 4 preconditions:
  - finite-difference gradients for the H, W and theta blocks;
  - F monotone on the T4-style corner grid;
  - the theta-block KKT residual confirmed by an evaluator written independently of the solver;
  - the analytic L_W of spec section 2 evaluated against a dense Hessian at p = 80;
  - under the map theta = s0 diag(1/sigma(W)) b, the survival term and the rho term equal the frozen candidate's
    terms (`ref_objective.R`) to 1e-12 relative. This last check is the identity on which the spec's rho and range
    carry-over rests.

  A failure is an implementation defect. It is fixed as a conformance bug, Q0 is rerun, and no other output exists
  yet.
- **Q1.** Scope: every grid fit, every tuned final fit, and the five inner-fold fits at each selected configuration.
  Pass, per fit: at least 27 of 30 starts converged, and no start with an F increase or a numerical failure. An F
  increase in any other qualification fit (non-selected BO evaluations) is a Q0-class implementation defect and is
  reported as such. Q1 is stricter than the fit-level Gate A precondition and is meant to be: BO feasibility needs
  all five inner-fold fits converged, and a design whose starts do not converge ordinarily cannot supply that.
- **Q2 (reported).** As in the spec, from the section 6 records.
- **Q3.** As frozen in spec section 9 (rulings 11 and 14): window sweeps 3,751-5,000, halves split at 4,375; (W, theta) stored at 3,750, 4,375 and 5,000; conditions (i)-(v); scope, not-evaluable and numerical-failure rules and active-set transitions as stated there. Pass: no evaluable start shows recurrence.
- **Q4.** Scope: the returned start of every grid fit and every tuned final fit, at sweeps 1, 10, 100 and the
  endpoint. Object: the Euclidean Hessian, in vec(W), of the smooth part of F at fixed (H, theta), restricted to the
  free coordinates W_ij > 0. A principal submatrix of a PSD matrix has lambda_max no larger than the full matrix's,
  so the restriction is still bounded by L_W at the same (H, theta). Power iteration on analytic HVPs: start vector
  ones on the free coordinates, normalized; stop at relative change < 1e-8 or 500 iterations. Power iteration
  approaches lambda_max from below, so the check is necessary, not sufficient: an estimate above the bound is a
  defect, and an estimate below it does not prove the bound. Pass: every estimate <= L_W (1 + 1e-8).
  Iterations used and nonconvergence of the power iteration are reported. Realized values are reported as in the
  spec.
- **Q5.** Cells: the centre configuration at n = 175 and the reference fit at n = 218 (replicates 1-5) and the tuned final fits, in S1,
  S2 and S3. Per fit: m = median over converged non-returned starts of (C_s - C_returned). Pass, per cell: the mean
  of m over replicates <= 0.01. This is a one-sided reading (section 10, P2). Loading r and the share of starts
  within 1e-4 relative F are reported.
- **Q6.** Cells: tuned final fits in S1 (R = 3) and S3 (R = 10) (ruling 15). Checkpoints t in {10, 50, 100, 500}, as frozen.
  Pass: for every t and both scenarios, the mean over replicates of (C_conv - C_t) >= -0.01. S2 is reported.
- **Q7.** Cell: S3 tuned final fits, R = 10 (ruling 15). On the 1,000 test columns, standardize eta_hat = X_G^T v~_G, true A
  and true B (section 1). Fit the joint least-squares regression eta_hat ~ A + B and set
  D_AB = |coef_A| - |coef_B|. Pass: mean D_AB - t_{0.95, 9} SD(D_AB) / sqrt(10) > 0. No top-up replicates (ruling 15). The marginal form and the
  realized A-B correlation on both scales are reported. A constant fitted predictor (theta = 0) gives D_AB = 0 for that replicate (PI
  ruling 16, before any S3 output); a nonfinite predictor is a numerical failure (Q1).

**Verifier controls (run after Q0, before any qualification fit; a failure stops the run, is reported, and is never
answered by adjusting a criterion).**
- **Q3 positive control.** The detector, unchanged, is applied to the failed estimator's four stalled T7 starts
  (1, 5, 9, 16; `sim_data(700, 3000, 218)`). It uses W at the first and last sweep of the final 25% of the
  exact-gated D-Q4 replay trajectory, with the theta-equivalent s0 b_j / sigma_j(W) for the carrying test. It must
  flag recurrence in all four. If it misses a start that D-Q4 showed to be g-driven, the detector is not validated.
- **Q3 negative control.** On a converged start, the detector returns "no recurrence" by definition. This case is
  checked so that the scope filter is exercised.
- **Q7 controls.** On S3 replicate 1's test set: eta_hat = true A gives D_AB > 0.9; eta_hat = true B gives
  D_AB < -0.9; eta_hat = A + B gives |D_AB| < 0.1. Each threshold is set from the arithmetic of exact inputs, not
  from any fit.
- **Q6 control.** For a start converged before sweep 10, every checkpoint difference is exactly 0.

## 8. Timing-only run (step 5, after this file is frozen)

The timing-only run uses one dataset, `sim_qual(49,100,001, "S1")`, outside every qualification namespace. It fits
the centre configuration at n = 175 (g = 5) and the reference fit at n = 218 (g = 12), with 30 starts each, and one inner-fold BO evaluation. The script
writes wall time per start, sweeps per second, and cores used. It discards start status, F, residuals and C before
writing, and its output file has no field for them. Its only use is compute planning, through section 9.

## 9. Compute and the contingency order (fixed now)

Work in starts: grid 180 fits x 30 = 5,400. Tuning runs: (20 + 40) evaluations x 5 inner folds x 30 starts, plus 30
for the final fit = 9,030 per run, x 35 runs = 316,050. The old estimator ran 0.113 s per sweep at n = 175 and 0.128
s at n = 218 (T7, `candidate/REPORT.md`). If every start hits the 5,000-sweep cap, that is about 870 core-hours for
the grid and about 50,000 for tuning. Starts that converge stop early, so the real figure scales with the median
sweeps to convergence, which only the timing run's wall times show.

If the timing-based projection exceeds the compute the PI allots, the PI rules before any qualification fit. The
only reductions available, in this order, are:

1. drop S3b;
2. reduce the S1 and S2 tuning runs from 5 to 3 replicates (grid unchanged);
3. reduce the S3 tuning replicates from 20 to 10 (Q7 then passes at delta = 1 with probability 0.90).

Nothing else changes: not the grid, the BO budget, the starts, the cap or any criterion.

### 9.1 PI ruling 15 (2026-10-05): ceiling, reduced design, fail-fast order

Written after the first three timing fits (wall times only: about 700 s per 30-start fit on 15 cores, roughly
3 core-hours per fit, on an oversubscribed machine) and before any qualification output exists.

- **P6 = 18,000 projected core-hours.** A planning ceiling, not a mid-job kill switch. The frozen reduction order is
  applied before submission until the projection is below it: (1) S3b dropped; (2) S1 and S2 tuning at 3 replicates
  each; (3) S3 tuning at 10 replicates. Projection: 16 tuning runs x about 870 core-hours = about 13,900, plus about
  500 for the 180-fit grid, leaving headroom for the verifier controls, timing variation and restarts.
- **Unchanged:** the 180-fit grid (replicates 1-5 of S1, S2 and S3), 30 starts, five inner folds, the BO budget
  (20 + 40), and every criterion. Reducing starts, folds or BO budget would test a different estimator.
- **No top-up.** If Q7's prespecified lower bound at R = 10 does not clear zero, qualification does not pass; no
  further replicates are added. Q7 passes at delta = 1 with probability 0.90 and at delta = 0.75 with 0.71.
- **Fail-fast execution order (frozen).** Each stage is evaluated by `qual_eval.R` on that stage's complete job set
  (counts from the manifest) before the next stage is submitted. A binding failure at any stage stops all remaining
  stages, and the failure is diagnosed on that failure only.
  1. **Stage 1, cheap and structural:** Q0 (already passed as smoke tests S1-S7 at p <= 300; rerun at the start of
     the stage), the verifier controls, the S3 realized-ratio check on the replicates generated for the stage, and
     the 180 grid and reference fits. Binding at this stage: Q1, Q3 and Q4 on grid and reference fits, Q5 at the
     centre configuration (n = 175) and the reference fit (n = 218).
  2. **Stage 2, S3 tuning (10 runs).** Binding: Q1, Q3, Q4, Q5 on the tuned fits; Q6 S3 part; Q7; the realized-ratio
     median check.
  3. **Stage 3, S1 tuning (3 runs).** Binding: Q1, Q3, Q4, Q5; Q6 S1 part.
  4. **Stage 4, S2 tuning (3 runs).** Binding: Q1, Q3, Q4, Q5.
  Qualification passes only if every stage passes. Q6 needs both S1 and S3; it is evaluated per scenario at its stage.
- **Longleaf:** parallelize across independent scenario x replicate jobs, not by adding cores per fit beyond useful
  scaling.

### 9.2 PI rulings 16-17 (2026-10-05, before any qualification output)

- **Ruling 16 (Q7, constant predictor).** If a tuned S3 fit has theta = 0, its predictor is constant and puts no
  effect on A or B: D_AB = 0 for that replicate. Dropping it would preferentially remove an unfavourable result.
- **Ruling 17 (execution).** The design phase is over. The harness is frozen after the Q3 positive control passes;
  afterwards only execution-blocking defects are fixed. No refactors, diagnostics, power analyses, resource
  optimizations, documentation passes or robustness checks; no inspection of timing-run convergence; no change to the
  ceiling, replicate counts, Q7 top-up or S3b. Fold-level parallelism inside a tuning run is enabled by a scheduler
  switch only (QUAL_FOLD_PAR = 5; results verified identical to sequential folds at p = 300).
- **Q3 positive-control input (implementation correction, recorded here).** Section 7 named the D-Q4 replay as the
  control's input. Those Newton-CG trajectories store no W inside the Q3 window, and on them start 9 converged, so the
  control could not be applied as written. The control uses the failed candidate's original projected-gradient
  trajectories for T7 starts 1, 5, 9 and 16 (sim_data(700, 3000, 218), T7 hyperparameters, fit seed 123, maxit
  5,000; all four stall at the cap), replayed with a recording-only hook (`q3pos/q3pos_replay.R`), and the unchanged
  detector with theta-equivalent s0 b_j / sigma_j(W) (`q3pos/q3pos_eval.R`). Pass rule unchanged: all four flagged.

## 10. Readings and points for the PI

- **P1. S3 SD-ratio space** (section 3.2). Reading: the ranked training X with unit-norm marker-masked loadings,
  giving c = (0.1911, 1, 1). The alternative literal reading (generator scale, c = (0.25, 1, 1)) gives a ranked ratio
  near 0.30, so it is a weaker stress. The fresh-context Claude review
  (`2026-10-05-estimator-redesign-claude-review.md`, S4) proposed the generator-scale reading. This file departs from
  it for the reason in section 3.2.
- **P2. Q5 direction and aggregation.** The spec says "within 0.01". This file reads it one-sided, as its stated
  reason ("cannot cost more than one C point") implies, with the per-fit median averaged over replicates per cell.
  S2 is kept binding, because the spec does not exclude it. In the null scenario the start-to-start C differences
  are pure noise. Decide whether S2 should be reported only.
- **P3. RESOLVED by PI ruling 14 (review B3 approved with refinements; spec section 9, Q3).**
- **P4. Q1 unit.** Per fit: at least 27 of 30 starts converged (review S2). Section 7 reworded the reason for Q1.
- **P5. RESOLVED by PI ruling 12 (review B1 = option (a)).** The H ridge is now lambdaH ||H||_F^2 / (2 ||X||_F^2),
  lambdaH = 1e-6 kept; HALS denominator 1 + lambdaH; a_H = (1 - alpha)(1 + lambdaH) / ||X||^2 (spec revision 4).
  Review B4 (test-set construction) is operational and is adopted in section 1.
- **P6. RESOLVED by PI ruling 15:** 18,000 projected core-hours; reduced design and fail-fast order in section 9.1.
