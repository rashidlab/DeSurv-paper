# Response to Jen Jen Yeh — Summary of Changes

---

## Comment 1 — Circular reasoning: "more reliably" / "clearer survival separation"

**Concern:** Claims in the abstract read as tautological if the model is trained on survival.

**Changes made:**

- *Abstract simulation sentence* (`paper.Rmd`): "more reliably recovers" replaced with a quantified ground-truth claim: precision against known marker genes is 0.50 vs. 0.07. This is measured against a known answer, not training survival.

- *Abstract closing sentence* (`paper.Rmd`): "clearer survival separation than their unsupervised counterparts" replaced with out-of-sample framing: the survival-aligned factors "generalize to five independent PDAC cohorts without retraining, with stronger out-of-sample survival stratification than their unsupervised counterparts (pooled validation HR per SD = 1.50, 95% CI 1.31–1.72, P < 0.001)." The validation cohorts had no access to survival during training.

- *Introduction paragraph 3* (`02_introduction_REVISED.Rmd`): Added sentence pre-empting the circularity concern, naming the two non-circular tests: ground-truth simulations (where the correct programs are known) and out-of-sample generalization (where the validation cohorts were withheld entirely during factorization).

---

## Comment 2 — CAF subtype conflation (iCAF vs. restCAF vs. deCAF proCAF)

**Concern:** "restCAF-associated (iCAF)" and "iCAF-associated stroma" used interchangeably, conflating subtypes.

**Changes made:**

Audited all CAF terminology across all Rmd files. Adopted iCAF framing throughout (not a combined restCAF/iCAF formulation):

- Results: D1 described as coupling "Classical tumor programs with an iCAF-associated stromal signature" and as a "Classical--iCAF-associated co-program" (replacing prior "classical--restCAF" language).
- Discussion: "iCAF-associated stromal signature" (replacing prior "restCAF-associated stromal signature").
- No location in the revised text equates iCAF and restCAF directly; the Discussion notes that D1 correlates with both Elyada iCAF and DeCAF restCAF categories per Fig. 3A.

---

## Comment 3 — Nomenclature consistency (basal/classical terminology)

**Concern:** Inconsistent "basal" vs. "basal-like," inconsistent capitalization.

**Changes made:**

Adopted uniform convention across all Rmd files: "Basal-like" (not "basal" or "basal-like") and "Classical" (capitalized consistently):

- Introduction: "Basal-like/Classical dichotomy."
- Results: factor labels "Classical tumor" and "Basal-like tumor" throughout; Fig. 3A heatmap context revised to match.
- SI: "iCAF-like and myCAF-like CAF subtypes from SCISSORS"; "Classical and Basal-like tumor subtypes" in heatmap captions; Puleo et al. subtype terminology verified ("Pure Classical," "Tumor Classical," etc.) and cited consistently.

---

## Comment 4 — "Exocrine-compositional variation that dominates standard NMF" — citation or own finding?

**Decision:** Skipped per author decision. No change made.

---

## Comment 5 — Treatment-naive cohorts: acknowledge as both strength and limitation

**Concern:** All training cohorts are treatment-naive. This is a strength but limits applicability to treated populations. Validation cohorts should also be confirmed as treatment-naive and nonmetastatic.

**Changes made:**

- *Results cohort description* (`04_results_REVISED.Rmd`, line 199): Added explicit statement that training cohorts consist of non-metastatic, treatment-naive specimens; added parallel statement (line 353) confirming all five external validation cohorts likewise consist of non-metastatic, treatment-naive specimens.

- *Discussion limitation 3* (`05_discussion_REVISED.Rmd`): Limitation 3 completely rewritten. The previous text discussed proportional hazards assumptions and algorithmic convergence only. The new text explicitly frames treatment-naive homogeneity as reducing confounding during factorization but limiting direct applicability to treated patient populations, names neoadjuvant chemotherapy as the specific gap, and identifies extending DeSurv to neoadjuvant/multimodal therapy cohorts as a natural and clinically important next step.

---

## Comment 6 — Remove "rather than assessing them retrospectively"

**Concern:** The phrase is unnecessary and awkward.

**Changes made:**

- Deleted "rather than assessing them retrospectively" from the Discussion opening sentence (`05_discussion_REVISED.Rmd`).

---

## Comment 7 — SCISSORS heatmap labels: add "-like"?

**Concern:** Heatmap labels "SCISSORS: iCAF" and "SCISSORS: myCAF" — should these be "iCAF-like" / "myCAF-like"?

**Changes made:**

- Updated all SCISSORS references to "iCAF-like and myCAF-like CAF subtypes from SCISSORS" in SI text and figure captions, matching the terminology used in the Leary et al. paper.

---

## Comments 8 & 9 — NMF k=7 comparison; comparison to published hazard ratios

**Concern 8:** Is the k=3 vs. k=7 comparison fair? Standard NMF at its optimal rank should be shown.

**Concern 9:** Comparison to published HRs from existing classifiers requested.

**Changes made:**

Several locations across the paper contained factually incorrect claims about Tables S5 and S6. The actual Table S6 values (computed from production RDS files) are:

| Model | Unadj HR (95% CI) | Adj HR (95% CI) | Adj P |
|---|---|---|---|
| DeSurv k=3 | 1.50 (1.31–1.72) | 1.33 (1.13–1.56) | 0.0005 |
| NMF k=7 (BO, α=0) | 1.48 (1.29–1.70) | 1.47 (1.28–1.69) | <0.001 |

Both methods are independently significant after adjustment. Prior text at multiple locations incorrectly claimed NMF k=7's effect "attenuates and is no longer significant" and that its prognostic content is "largely redundant with existing classifiers." These claims were corrected or deleted. The conceptual point — that the NMF k=7 BO model is not standard NMF, because it uses survival information to select k and ntop via BO (α=0 means survival does not enter the factorization, but does enter hyperparameter selection) — was added to make the comparison accurate.

*Specific edits:*

1. *Abstract* (`paper.Rmd`): Updated hardcoded HR from 1.45 (1.29–1.63) to 1.50 (1.31–1.72) to match the dynamically computed production run value.

2. *Results §"Survival-aligned programs generalize"* (`04_results_REVISED.Rmd`): Replaced "produces a fragmented structure whose prognostic signal largely overlaps with existing molecular classifiers" with: "However, this concordance gain requires more than twice as many factors. DeSurv at k=3 retains a significant hazard ratio after adjustment for PurIST and DeCAF (SI Appendix, Table S6), achieving equivalent independent prognostic value more parsimoniously."

3. *Results §"Factorization rank k=3 is robust"* (`04_results_REVISED.Rmd`): Deleted the sentence asserting that k=7's "prognostic content is largely recoverable from existing classifiers applied to the same data" — directly contradicted by Table S6, where NMF k=7 BO retains HR 1.47, P < 0.001 after adjustment.

4. *Discussion paragraph 3* (`05_discussion_REVISED.Rmd`): Removed "whose prognostic content overlaps with existing molecular classifiers (SI Appendix, Tables S2--S3)" from the sentence describing NMF's fragmented structure.

5. *SI — Table S5 context* (`si_appendix.Rmd`): Replaced "k=7's prognostic signal largely overlaps with existing molecular classifiers" with the accurate k-sensitivity breadth result: k=3 retains adjusted significance across more supervision strengths than k=7 (5 of 7 versus 2 of 7; Tables S2--S4).

6. *SI — Table S6 introduction* (`si_appendix.Rmd`): Rewrote the paragraph introducing Table S6. The prior text falsely stated DeSurv's unadjusted effect was "substantially larger," that NMF k=7's effect "attenuates and is no longer significant," and that NMF k=7 is "largely redundant with existing classifiers." The new text states both methods retain significant independent prognostic value after adjustment, and clarifies the key distinction: the NMF k=7 model uses survival to select hyperparameters via BO but not to guide the factorization (α=0); DeSurv uses survival at both levels, concentrating the prognostic signal into three factors rather than seven.

7. *Discussion paragraph 3* (`05_discussion_REVISED.Rmd`): Added sentence making the k=3 vs. k=7 comparison explicit in the main text: "Even so, DeSurv at k=3 achieves equivalent independent prognostic value with fewer factors: both methods retain significance after classifier adjustment (SI Appendix, Table S6), and k=3 is more robust across supervision strengths (5 of 7 versus 2 of 7 for k=7; SI Appendix, Table S3)."

8. *Discussion paragraph 4* (`05_discussion_REVISED.Rmd`): Added two sentences addressing Comment 9 — that DeSurv is not proposed as a competing molecular classifier (discrete subtype labels would be a downstream application), and that Table S6 provides the relevant comparison with existing classifiers: the DeSurv k=3 linear predictor retains a significant adjusted HR after conditioning on PurIST and DeCAF (HR = 1.33, 95% CI 1.13–1.56, P < 0.001).

---

## Comment 11 — "All samples nonmetastatic?"

Addressed as part of Comment 5. Both training and validation cohorts confirmed as non-metastatic, treatment-naive, stated explicitly in Results.
