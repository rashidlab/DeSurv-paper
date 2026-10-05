# Ruling: decision-gated evaluation of a regularized DeSurv estimator (PI, 2026-10-04)

Context: the audit verification (`2026-10-04-audit-verification.md`, section 2) showed the reported fits do not
minimize Eq. 1 and that continued optimization of Eq. 1 changes D1 (full-data probe: loading r 0.69, top-270
overlap 0.64) while lowering held-out C. The manuscript remediation (DeSurv-paper-clean c92ae83 / 999fe7e)
documents the implemented estimator accurately but does not show that D1's biology is independent of the
optimizer. The PI's concern: the current algorithm may be perceived as flawed, and the biology would be better
justified by a properly specified model.

## Ruling (PI's words)

> Proceed with the gated investigation. First characterize whether the existing continued-optimization
> sensitivity preserves the major D1 biological conclusions. Then develop and freeze a principled no-zeta,
> convergently optimized, explicitly regularized estimator. Hyperparameters must be selected solely from
> prespecified statistical criteria using fresh nested resampling; neither D1 similarity nor the current C-index
> may enter tuning or estimator selection. Freeze the estimator-selection rule before fitting the final candidate.
> After the selected model is frozen, evaluate D1 biological concordance as an independent scientific validation.
> Do not modify the manuscript's primary estimator until this experiment is complete.

## Structure

1. **Triage (not validation).** On the existing continued-optimization fit: reference-program correspondence,
   basal/classical relationship, PurIST, DeCAF, pathway/enrichment interpretation, and sample-level D1 scores
   (not only gene lists).
2. **Freeze the candidate estimator** before evaluating it: scale-invariant regularizer, the exact objective,
   optimizer (no zeta), line search and convergence criteria, normalization convention, initialization, tuning
   space, inner-loop tuning procedure, number and definition of resamples (fresh partitions; the item-4
   partitions and the external cohorts are development data and may not be used for any choice), primary
   performance criterion, hyperparameter-selection rule.
3. **Primary-estimator decision on statistical properties only**, rule frozen before the final candidate is fit:
   nested out-of-sample discrimination against a prespecified noninferiority margin (fixed on scientific grounds,
   not derived from observed values such as 0.638 vs 0.621, with outer-resample uncertainty), calibration if
   relevant, convergence, numerical stability, stability to initialization and resampling, complexity.
   D1 biology is not an input.
4. **Biology as the independent validation**, run only after step 3 is frozen and committed.

## Pre-declared outcomes

| Prediction (step 3) | Biology (step 4) | Consequence |
|---|---|---|
| Passes margin | Recovers D1 biology | Candidate may become primary |
| Fails margin | Recovers D1 biology | Original algorithm stays primary; regularized fit is a major robustness analysis |
| Any | Materially different | Determine which parts of the D1 interpretation survive and narrow claims; hold or rework if the central conclusion disappears |
| Fails badly | Different | Not evidence against D1 (poor estimator, poor control); needs another principled route (simpler supervised factor model or different regularization structure) |

Prior evidence, stated so it is not forgotten: on fresh partitions (seeds 2001-2003) every regularized pilot arm
lost held-out C against the implemented procedure (best -0.017). One arm gained +0.023 on the already-examined
item-4 folds, and that gain did not replicate on fresh folds (-0.030); see
`analysis/audit-verification/regularization_design_memo.md` section 5. The candidate is a competing method, not a
presumed improvement.

## Process safeguards (coordinator)

- The specification (step 2) and the selection rule (step 3) are committed before step 1 results are read, so
  the triage cannot steer the regularizer's design. The commit hash of the frozen specification is recorded in
  every later output.
- Partition seeds 123-126 (item 4) and 2001-2003 (feasibility check) have been examined and may not be reused for
  tuning, selection or evaluation.
- Selection scripts compute no biology outputs; the step-4 script runs only against a committed selection.
- Open PI decision required before step 2 is frozen: the noninferiority margin for held-out C.
- The closure items (F1 disclosure, F4 sentence, 09a fallback, manifest) proceed independently. The F2 wording
  fixes ("continued gradient-based optimization", "the finite fitting procedure provides substantial implicit
  regularization") apply regardless of the outcome.
