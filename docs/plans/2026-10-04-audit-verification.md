# Verification of the independent scientific audit (2026-10-04)

Audit verified: `docs/HANDOFF-independent-scientific-audit-2026-10-04.md` in `~/Downloads/DeSurv-paper-clean`
(branch revision/desurv-step6, HEAD fb731e7, working tree clean before and after).
Scope per owner rulings of 2026-10-04: read-only verification; F2 verified against DeSurv 1.0.1 (manuscript chain,
site library, RemoteSha afb00d5) and DeSurv 1.2.0 (rev 3, private library lib-rev3, source DeSurv-clean
feat/reported-estimator 1f3d3f0). The owner additionally authorized, during the session, a scratch-only
fold-level CV experiment (item 4 below). No code, cached result, manuscript, or pipeline script was changed.
F2 was verified by the coordinator; F1 and F3-F7 by a subagent whose central claims (F5, F6) the coordinator
re-checked against source.

## 1. Bottom line

| Finding | Audit | Verdict | Revised severity | Disclosed in manuscript |
|---|---|---|---|---|
| F1 tie derivative | HIGH | CONFIRMED (defect); CV impact negligible | MEDIUM | No (likelihood ties yes, derivative defect no) |
| **F2 optimizer vs objective** | HIGH | **CONFIRMED in 1.0.1 and 1.2.0; consequence reversed** | **HIGH, different remedy** | Mechanisms yes (SI L414-440, L1816; Methods L81) |
| F3 final config outside 1-SE | MEDIUM | CONFIRMED | LOW-MEDIUM | Yes (Methods L86, SI 811-879) |
| F4 grid filter outside folds | MEDIUM | CONFIRMED | LOW | No |
| F5 supervised-PCA surrogate | MEDIUM | CONFIRMED | LOW-MEDIUM | No; Methods L119 inaccurate |
| F6 zero-fill before ranking | MEDIUM | CONFIRMED mechanism; impact negligible | LOW | No; contradicts Methods L163 wording |
| F7 provenance gaps | MEDIUM | CONFIRMED | MEDIUM (reproducibility) | Partly |

The gating answer: **F2 holds in 1.2.0 exactly as in 1.0.1** (the loss, beta and W code is byte-identical between
the versions; only the tie fix differs), so every rev 3 BO and rank-selection result inherits it. But the audit's
implied remedy (correct the optimizer, rerun) is wrong. Optimizing the stated objective exactly **overfits and
loses held-out C** (20-fold paired difference -0.067, SE about 0.01, 19 of 20 folds worse). The stall of the
implemented algorithm is acting as implicit early-stopping regularization. The validated estimator is the
algorithm, not the minimizer of Methods Eq. 1.

## 2. F2: coefficient and W updates inconsistent with the stated objective

### 2.1 Source (both versions)

`git diff afb00d5 1f3d3f0 -- src/functions.cpp` touches only `cox_suffix_stats_full` (tie-group end
accumulation) and adds the `cox_stats_cpp` test hook. Everything below is identical in 1.0.1 and 1.2.0
(line numbers from 1.2.0 `src/functions.cpp`):

- `calc_loss_cpp` (L140): survival term `cs.loglik * 2.0 / n_event - lambda*penalty_beta`, penalty on
  `beta % sdZ` (standardized coefficients b).
- `cdfit_cox_dh_one_lambda` (L341): coordinate step uses `xwr / n_event` and curvature `xwx / n_event` with the
  same lambda, so the sweeps target `-(1/n_e) l + lambda P(b)`, equivalently penalty 2*lambda relative to the
  objective's likelihood weight.
- `update_beta_cpp` (L426): raw `beta0` passed as the standardized starting point (`beta0_keep`, no
  multiplication by sd); the candidate is accepted only if it does not worsen the stated penalized survival loss,
  otherwise a damped step toward it, otherwise beta0 is retained. With `alpha == 0 && !flag_beta`, beta is
  returned unchanged (frozen at its start) inside the loop.
- `update_W_damped_backtrack` (L190): multiplicative proposal with the Cox gradient rescaled so its Frobenius
  norm equals `alpha * ||grad_NMF||` (heuristic zeta), clamp [0.2, 5], backtracking up to 30 times on the full
  loss; if all trials are rejected W is unchanged.
- `optimize_loss_cpp` (L540) stopping rule: `while ((epsW > tol || epsH > tol) && it < maxit)` with
  `eps = 1 - cos(previous, current)` for W and H as whole matrices. The objective does not enter the rule. A fully
  rejected W step gives epsW = 0 exactly, so the rule can be met by H converging with W frozen.
- Additional observation not in the audit: column normalization of W each iteration makes
  `||W||_F^2 = k` constant, so the lambdaW penalty is the constant 1/(2p) (stored `penalty_W` 0.000254 =
  1/(2*1970)). lambdaW is vacuous as implemented and cannot serve as a regularizer.

The SI discloses the 2*lambda effective penalty, the unconverted warm start, the acceptance rule, zeta, and the
termination-by-H behavior. The issue is not disclosure.

### 2.2 Independent numerical checks

Independent base-R implementation of Eq. 1 (Breslow ties exact), analytic W gradient including the dependence of
sd(Z) in the penalty, finite-difference checks, exact beta subproblem solved by FISTA to KKT residual below 1e-9.

| | Manuscript k3 consensus (1.0.1) | Rev 3 k2 consensus, reported recipe (1.2.0) | Inline k12 at rank-selection eval-40 config (1.2.0) | Same, nu = 0 (1.2.0) |
|---|---|---|---|---|
| Independent vs stored objective | 2.941507360 = 2.941507360 | 5.150104262 = 5.150104262 | 5.569184445 = 5.569184445 | 5.443865851 = 5.443865851 |
| Last W change at termination | 3.3e-16 | 1.5e-8 | 0 | 1.1e-16 |
| ||grad_W L|| (FD agrees to 7 digits) | 0.4985 | 0.7063 | 2.504 | |
| dL for feasible step pmax(W - 0.1 G) | -0.0111 | -0.0249 | -0.0200 (step 0.01) | |
| survival+penalty / reconstruction gradient norm | 23x | 80x | about 1,500x | |
| beta KKT residual at stated lambda | (0.292, 0, -0.066) | (0.335, -0.172) | about 0.025 per active coordinate | |
| beta KKT residual at 2*lambda | (-0.0017, 0, -0.0008) | (-0.0075, -0.0091) | | |
| exact b* at stated lambda | (-1.010032, 0, 0.172344) | (-0.7926, 0.3473) | | |
| training C, stored vs exact beta | 0.7860 vs 0.7858 | 0.7835 vs 0.7835 | 0.8676 vs 0.8672 | |

All audit numbers for the k3 fit reproduce to the printed digits (objective 2.941507, gradient norm 0.4985242,
descents 0.0002197453 and 0.01111326, KKT 0.2916716 and -0.06562946). Corrections to the audit:

- The audit's candidate beta (-1.010101, 0.000454, 0.1732) is not the optimum; the exact solution has the middle
  coefficient exactly 0.
- At k3 the stored b is the exact optimum of the **2*lambda** problem (KKT about 1e-3). The normalization mismatch
  fully explains the beta discrepancy; the raw warm start is immaterial once the sweeps converge (package
  `update_beta_cpp` from start 0 or from stored beta returns the same point). At k12 the stored beta lies between
  the lambda and 2*lambda optima: the acceptance rule pins beta where the last improving step left it, so beta is
  path-dependent.
- Replacing beta by the exact optimum (W fixed) changes the five external C-indices by at most 0.0011 and the
  20-fold held-out C by +0.0002. **The beta mechanism is real and scientifically negligible.**
- The W mechanism is the substantive one. DeSurv's own algorithm cannot move W further: restarting 1.0.1 from the
  archived fit with tol 1e-14 ran 803 iterations with all 803 W steps rejected; W unchanged (column r = 1),
  objective 2.941507 to 2.941495 from H and beta only. A stricter tolerance or more iterations does not change the
  estimate.

### 2.3 Impact: what exact optimization of the stated objective does

Probe: from each fit, block descent on Eq. 1 with no package code (H multiplicative, W projected gradient with
Armijo on the full loss, beta exact, columns renormalized; objective invariant under the renormalization because
lambdaH = lambdaW = 0). Probes are local and capped (not certified optima). External scoring uses the production
rule of `code/05` `compute_val_cindex()` (union of top-ntop genes, frozen W and beta); the archived column
reproduces the cached `val_latent_desurv_tcgacptac` C-indices to four decimals.

| | Published k3, 600 iter | Inline k12, 400 iter | Inline k12 nu = 0, 400 iter |
|---|---|---|---|
| Objective | 2.942 to 2.273 | 5.569 to 4.944 | 5.444 to 4.587 |
| Log partial likelihood | -598.5 to -392.6 | -549.7 to -441.9 | |
| Training C | 0.786 to **0.997** | 0.868 to **0.991** | 0.876 to **0.992** |
| External C, n-weighted mean | 0.618 to **0.597** | 0.592 to **0.565** | 0.586 to **0.559** |
| Per cohort (Dijk/Moffitt/PACA-arr/PACA-seq/Puleo) | .607/.562/.649/.634/.636 to .569/.538/.620/.653/.616 | .586/.553/.616/.580/.607 to .566/.497/.609/.557/.585 | .568/.557/.615/.587/.596 to .547/.499/.616/.566/.575 |
| Programs | D1 loading r 0.69, top-270 overlap 0.64, score r 0.64; D2 r 0.9997, D3 r 0.996 (overlap 0.99, 0.96); union Jaccard 0.77 | only nonzero-beta columns move (r 0.85, 0.99; overlap 0.83, 0.97); zero-beta columns r = 1 | two columns move (r 0.81, 0.90; overlap 0.80, 0.79); others overlap at least 0.93 |

The change concentrates in the program carrying the largest beta (the survival gradient for column j is
proportional to beta_j, so zero-beta columns receive no supervision at all). At k3 that is D1, the paper's
central program.

### 2.4 Decisive check: fold-level held-out C (item 4, owner-authorized, scratch only)

Published k3 configuration at full precision (alpha 0.333935556177224, lambda 0.349448737941329,
nu 0.0556512168045738), cached 1,970-gene training matrix, 4 partitions x 5 folds stratified by cohort x event.
In each training fold: DeSurv 1.0.1 `desurv_fit` (30 inits, tol 1e-7, maxit 4000), then the same fit continued
by exact descent of Eq. 1 on the training fold; held-out C on the test fold with the top-270 union rule.

| Descent iterations | Objective | Training C | Held-out C | Paired diff vs stalled | Folds worse (of 20) | Lead-program loading r |
|---|---|---|---|---|---|---|
| 0 (stalled, as implemented) | 2.932 | 0.726 | **0.638** | | | 1.00 |
| 25 | 2.526 | 0.891 | 0.611 | -0.027 (SE 0.009) | 13 | 0.89 |
| 50 | 2.382 | 0.931 | 0.590 | -0.048 | 16 | 0.85 |
| 100 | 2.242 | 0.965 | 0.580 | -0.058 (SE 0.012) | 18 | 0.80 |
| 400 | 2.079 | 0.996 | 0.569 | -0.069 | 17 | 0.74 |
| 1000 | 2.033 | 0.999 | **0.571** | -0.067 (SE 0.009 fold-level, 0.011 partition-level) | 19 | 0.71 |

Stalled W with exact beta: held-out 0.6378 vs 0.6376. All 20 stalled fits terminated on the stopping rule; in 12
the last W change was exactly 0. Limitations: fixed full-data gene filter (not fold-local), random-init
`desurv_fit` rather than the 100-seed consensus, one configuration. The direction and size of the effect are far
outside these caveats.

### 2.5 Verdict and consequence

- **CONFIRMED** in 1.0.1 and 1.2.0: the beta and W updates are inconsistent with the stated objective, and the
  converged fits are not stationary for it.
- **The audit's consequence is reversed.** The stated objective is under-regularized: nothing in Eq. 1 limits
  how far the supervised columns of W tilt toward the 139 training events (the penalty on standardized b is
  invariant to that tilt, and lambdaW is vacuous under column normalization). The zeta-scaled proposal stalls
  before that happens, and the stall is what generalizes. "Fixing the optimizer" alone produces a worse
  estimator.
- **D1 is determined by the algorithm's stopping point, not by Eq. 1.** The program-stability results
  (`code/20`, `code/21`, consensus recovers the basis at r 0.97-0.99) were all produced by the same algorithm;
  they establish reproducibility of the procedure, not that D1 solves the written objective. Defensible:
  "this fitting procedure identifies D1 reproducibly and it validates externally." Not defensible: "the DeSurv
  objective identifies D1" or any language implying the reported fits minimize Eq. 1.
- **Rev 3 rank selection** is internally valid as a comparison of the implemented procedure across ranks. At
  k12 the fit is 9 unsupervised columns (beta = 0, cosine 1 under exact descent) plus 3 lightly supervised ones,
  so "twelve survival-shaped programs" is not a supportable reading. The ruling against reinterpreting the
  manuscript around k12 should stand.
- Remediation options: (a) manuscript: define the estimator as the algorithm (zeta-scaled block descent from
  consensus initialization, terminated by the stated rule), with Eq. 1 as its guiding criterion and the stall
  stated as implicit regularization, plus a sensitivity analysis; or (b) method: add an explicit,
  scale-invariant regularizer tuned by BO, optimize it with a convergent algorithm, and revalidate. A design memo
  for (b) is in preparation separately.
- Smallest next decisive check for (a): repeat section 2.4 with the production recipe (fold-local filter,
  consensus init) at the published configuration; for (b): the same 20-fold design with the regularized
  estimator, target held-out C at least 0.638 with stable D1 top genes.

## 3. alpha = 0 comparators and the second stage (question raised during verification)

Code facts:

- BO CV fold scoring (`R/cv_cold.R:307-327`, 1.0.1) is identical for every alpha: `desurv_fit` with beta0 = 0,
  then `fit$beta` applied to top-gene-truncated held-out scores.
- At alpha = 0, W and H are plain NMF (no Cox term in the W step or the loss) and beta stays at 0 in the loop;
  the reported beta comes from the same post-loop penalized-Cox refit (`flag_beta = TRUE`) used at alpha > 0, with
  the arm's own BO-tuned lambda and nu. The alpha = 0 BO arm and the fair comparator (`code/22`, a0k3) therefore
  already share DeSurv's second stage. Only the `NMF::nmf` comparator (`fit_std_desurvk`, `code/04` L217) uses a
  different one (`R/fit_cox_model.R`: cv.glmnet over a nu grid, 1-SE lambda).

Matched second stage (item 3): each program set scored externally with its own beta, a common unpenalized Cox, and
a common exact elastic net at DeSurv's (lambda, nu), all refit on the same truncated-support training scores.

| External C, n-weighted mean | own beta | common unpenalized | common elastic net |
|---|---|---|---|
| DeSurv k3 (alpha 0.334), ntop 270 | 0.6181 | 0.6180 | 0.6182 |
| alpha = 0 fair comparator k3, ntop 270 | 0.5095 | 0.5132 | 0.5085 |
| alpha = 0 fair comparator k3, own ntop 299 | 0.5073 | 0.5032 | 0.5087 |

The second stage moves C by less than 0.004; the gap is in the programs. Supervision of the (stalled) programs
matters. The fair comparator's own-beta column reproduces cached `val_cindex_desurv_a0k3` exactly.

## 4. F1, F3-F7

### F1, tied-time Cox derivative: CONFIRMED, MEDIUM (was HIGH)

1.0.1 `cox_suffix_stats_full` accumulates `h` and `B` in storage order, so subjects early in a tie group miss
tied events stored after them; the log-likelihood itself is correct. Finite differences on random tied data:
max gradient error 0.508 (1.0.1) vs 7e-9 (1.2.0); gradient sum (should be 0) 3.16 vs 1e-14; untied data identical.
Training data: 60 subjects and 30 events at tied times. At the archived linear predictor the 1.0.1 gradient
differs for 18 subjects (max 0.030). Downstream CV C shift -0.0003 to +0.0003 (k3 ladder report section 8, S1 to
S2). Effect on archived W and beta not measured. Remediation: done in 1.0.2/1.2.0; optional deterministic refit
comparison.

### F3, final hyperparameters outside the recomputed one-SE set: CONFIRMED, LOW-MEDIUM

`code/03:112-123` overwrites only k in the global (k7) optimum. Full-precision values: threshold 0.6435500975,
final configuration 0.6419248411 (SE 0.0069), best k3 0.6466307793. Matches amendment 1c section 1. Fully
disclosed in Methods L86 and SI 811-879. No further action beyond what is disclosed.

### F4, sensitivity grid filters genes outside folds: CONFIRMED, LOW

`code/06:80-101` filters once on pooled training data; fold-local filters differ by 73-95 genes per fold (exact
audit counts). Unsupervised filter, no outcome leakage; primary BO unaffected. Not disclosed. Remediation: one SI
sentence, or rerun the k3 and k9 grid rows fold-locally.

### F5, tuned supervised PCA scored through a correlation surrogate: CONFIRMED, LOW-MEDIUM

`code/19:83-87` projects external data through gene weights `cor(gene, training score)`, not the native
`superpc.predict`. Native predictor: pooled external C 0.596 vs reported surrogate 0.603 (per cohort exact to the
audit). Correcting it strengthens DeSurv's relative position. Methods L119 ("the fitted model was frozen and
applied") is inaccurate. Remediation: native prediction, update the pooled value and the D2 |r| in Results.
Coordinator re-checked the source.

### F6, absent assay genes zero-filled before ranking: CONFIRMED mechanism, LOW

`code/05:55` `zero_fill_missing = TRUE` (`R/preprocess_helpers.R:69-81`). Dropping missing genes instead changes C
by +0.00026 (Dijk) and 0 elsewhere (constant inserted ranks shift scores by a constant). The Methods
"no imputation" sentence and SI 936-946 omit the zero fill. Remediation: manuscript-only sentence. Coordinator
re-checked the source.

### F7, provenance: CONFIRMED, MEDIUM for reproducibility

`code/01_install.R` installs unpinned `../DeSurv` when present; promised session-info file absent
(`paper/paper.Rmd:150`); versioning chunk `eval=FALSE` reading a targets store not used here; subtype object
fallback to `../DeSurv-paper/data/derv/` with no generator; no raw-to-curated survival derivation. Remediation:
manifest with hashes, sessionInfo, pin enforcement, fix the session-info promise.

### Positive findings spot-checked

Per-cohort external C, PACA-AU dedup (570 patients, 388 events, pooled 0.6135), stored beta, union support 810,
full-data preprocessing equal to cache: all match the audit.

## 5. Answers to the audit's closing questions

- **Does the current manuscript chain differ from the audited chain?** No. HEAD fb731e7, DeSurv 1.0.1
  (afb00d5), cached results unchanged.
- **Are the two HIGH findings mathematically correct?** Yes, both. F1 is downgraded to MEDIUM on measured impact.
  F2 stands as HIGH, but its scientific consequence is the opposite of the audit's framing (section 2.5).
- **Are newer corrections complete and tied to regenerated outputs?** The tie fix (F1) is complete in 1.2.0 and
  measured. Nothing in 1.2.0 addresses F2; the rev 3 BO and rank selection inherit it unchanged.
- **Which principal results remain defensible?** External associations, cohort counts, HRs and C-indices of the
  reported estimator (reproduced); supervision of the programs beats an alpha = 0 comparator under a matched
  second stage (section 3); program reproducibility under the implemented procedure. Not defensible: any claim
  that the reported fits minimize Eq. 1 or that D1 is the objective's solution.
- **What must be resolved before submission?** The F2 framing in Methods (estimator defined as the algorithm,
  with the implicit-regularization finding and a sensitivity analysis), or a regularized, truly optimized
  estimator revalidated. F5 and F6 Methods sentences. F7 provenance items. F4 disclosure.

## 6. Provenance of this verification

Scripts, outputs and the subagent's F1/F3-F7 report are in `analysis/audit-verification/` (see its README).
Scripts run from the root of `~/Downloads/DeSurv-paper-clean` against its cached `results/`; 1.2.0 runs read
`DESURV_120_LIB`. Two relocated scripts (`f2_check.R` on the k3 fit, `item3_matched_beta.R`) were rerun from the
repository copy and reproduced the reported values.
