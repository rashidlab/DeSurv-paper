# Jen Jen Revision Plan

Skipping comment #4 (exocrine citation/claim) per author decision.

---

## 1. Abstract rewrite — defuse circular reasoning (Comment #1)
**Status:** Complete

**Files:** `paper/paper.Rmd`, `paper/02_introduction_REVISED.Rmd`

- **Abstract simulation sentence:** Replace vague "more reliably recovers" with quantified ground-truth precision (0.50 vs. 0.07 against ground-truth marker genes).
- **Abstract closing sentence:** Replace "clearer survival separation than their unsupervised counterparts" with out-of-sample framing — "generalize to five independent PDAC cohorts without retraining, with stronger out-of-sample survival stratification."
- **Intro para 3:** Add sentence explicitly pre-empting circularity — two non-circular tests: ground-truth simulations and out-of-sample generalization to five independent cohorts.
- **Intro para 4:** Match abstract's new framing for the summary sentence.

---

## 2. CAF nomenclature audit (Comment #2)
**Status:** Complete

**Files:** `paper/02_introduction_REVISED.Rmd`, `paper/04_results_REVISED.Rmd`, `paper/05_discussion_REVISED.Rmd`

D1 correlates with both Elyada iCAF and DeCAF restCAF per Fig 3A — describe precisely throughout. Do not write anything that reads as "iCAF and restCAF are the same thing."

- Audit every mention of iCAF, restCAF, myCAF, deCAF across all Rmd files.
- Settle on consistent phrasing, e.g. "iCAF/restCAF-associated stromal signatures."
- Add note in Discussion that these two CAF taxonomies overlap but are not equivalent.
- Check intro line 23: "iCAF-associated stroma" — qualify or expand.
- Check results lines 259–266, 359 for any conflation.
- Check discussion line 8: "restCAF-associated stromal signature" — ensure iCAF is also named.

---

## 3. Nomenclature consistency — basal/classical terms (Comments #3 and page 10)
**Status:** Complete

**Files:** `paper/paper.Rmd`, `paper/02_introduction_REVISED.Rmd`, `paper/04_results_REVISED.Rmd`, `paper/05_discussion_REVISED.Rmd`, `paper/si_appendix.Rmd`, `code/08_figures.R`

- Pick one convention: "basal-like" (not "basal"), "classical" with consistent capitalization.
- Verify Fig 3A heatmap labels match text conventions (in `code/08_figures.R`).
- Verify Puleo et al. terminology ("Pure Classical," "Tumor Classical," etc.) is correctly cited in SI.

---

## 4. ~~Exocrine citation or own finding (Comment #4)~~
**Status:** Skipped

---

## 5. Treatment-naive cohort — acknowledge prominently (Comment #5)
**Status:** Complete

**Files:** `paper/04_results_REVISED.Rmd`, `paper/05_discussion_REVISED.Rmd`

- Verify factually whether all training and validation cohorts are treatment-naive (check cohort metadata / Table S1).
- Add explicit statement in Results cohort description.
- Sharpen limitation 3 in Discussion to note applying DeSurv to treated/contemporary cohorts as a future direction.

---

## 6. Remove "rather than assessing them retrospectively" (Comment #6)
**Status:** Complete

**File:** `paper/05_discussion_REVISED.Rmd`, line 6

Simple deletion: "Incorporating survival outcomes directly into NMF factorization rather than assessing them retrospectively reorganizes..." → "Incorporating survival outcomes directly into NMF factorization reorganizes..."

---

## 7. SCISSORS labels — add "-like" if that's what Leary uses (Comment #7)
**Status:** Complete

**Files:** `paper/si_appendix.Rmd`, `code/08_figures.R`

- Check Leary SCISSORS paper (ref 25) for exact label wording (iCAF-like / myCAF-like vs iCAF / myCAF).
- If they use "-like," update heatmap labels in `code/08_figures.R` and SI caption accordingly.

---

## 8. Make k=3 vs NMF-k=7 comparison more visible in main text (Comment #9)
**Status:** To do

**File:** `paper/05_discussion_REVISED.Rmd`, paragraph 3 (around line 10)

Add one sentence to the concordance plateau paragraph, e.g.: "NMF achieves comparable concordance only at k=7, and this signal is largely captured by existing classifiers (Tables S5–S6), whereas DeSurv k=3 retains independent prognostic value after adjustment for PurIST and DeCAF."

---

## 9. Follow up with Jen Jen on published HR comparison (Comment #10)
**Status:** To do (not a manuscript edit)

Email asking which published HRs she has in mind (PurIST? DeCAF? Moffitt?). Point to Tables S3/S6 (PurIST/DeCAF-adjusted analysis) as the stronger version of this test already in the paper.

---

## Execution order
4 (done) → 6 → 1 (done) → 3 → 2 → 5 → 8 → 7 → 9
