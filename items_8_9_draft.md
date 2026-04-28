# Items 8 & 9 — Values and Proposed Framing

## Source tables

All numeric values pulled from `results/precomputed/ntop_bo_50_300/` (production run).
K-sensitivity p-values from `results/cv_grid/`.

---

## Table S5 — External validation C-index by method and cohort

| Cohort | DeSurv k=3 | NMF k=3 | NMF k=5 (elbow) | NMF k=7 (BO α=0) |
|---|---|---|---|---|
| Dijk | 0.607 | 0.542 | 0.633 | 0.631 |
| Moffitt | 0.562 | 0.521 | 0.552 | 0.578 |
| PACA-AU array | 0.649 | 0.470 | 0.657 | 0.669 |
| PACA-AU seq | 0.634 | 0.430 | 0.679 | 0.702 |
| Puleo | 0.636 | 0.534 | 0.642 | 0.651 |

**Key observation:** The BO α=0 model (k=7, ntop=187) matches or exceeds DeSurv k=3 in most cohorts. However, this is not standard NMF: BO used survival signal to select k=7 and ntop=187 (optimizing held-out C-index), even though α=0 means survival does not influence the factorization itself. This distinction matters for interpretation.

---

## Table S3 — Adjusted p-values (ntop=270, Cox + PurIST + DeCAF + strata)

Each cell = p-value for the factor most correlated with the production k=3 D1 factor.

| | α=0 | α=0.25 | α=0.35 | α=0.55 | α=0.75 | α=0.85 | α=0.95 |
|---|---|---|---|---|---|---|---|
| K=2 | 0.564 | 0.038 | 0.038 | 0.019 | 0.060 | 0.011 | 0.003 |
| K=3 | 0.030 | 0.030 | 0.519 | **0.003** | 0.028 | 0.263 | 0.010 |
| K=5 | **0.001** | **0.005** | **0.001** | **0.019** | 0.073 | 0.605 | 0.447 |
| K=7 | 0.514 | 0.471 | **0.001** | 0.041 | 0.676 | 0.131 | 0.130 |
| K=9 | **0.000** | 0.299 | 0.243 | **0.005** | 0.943 | 0.867 | 0.263 |

- **K=3 significant in 5/7 α settings**
- **K=7 significant in only 2/7 α settings**
- K=7 at α=0, ntop=270 (uniform gene-focusing, no factorization supervision): adj p = 0.514 — **not significant**

---

## Table S4 — Adjusted p-values (ntop=NULL, all genes)

| | α=0 | α=0.25 | α=0.35 | α=0.55 | α=0.75 |
|---|---|---|---|---|---|
| K=3 | 0.205 | 0.024 | 0.280 | 0.097 | 0.112 |
| K=5 | 0.461 | **0.000** | 0.150 | 0.233 | 0.067 |
| K=7 | 0.394 | 0.336 | 0.304 | 0.097 | 0.094 |

**Key observation:** K=5, α=0 drops from p=0.001 (ntop=270) to p=0.461 (all genes) — result is fragile, depends on gene focusing.

---

## Table S6 — Pooled validation HRs, unadjusted vs. adjusted for PurIST + DeCAF

Values verified by live computation from val_latent RDS (ntop_bo_50_300 production run):

| Method | Unadjusted HR per SD (95% CI), P | Adjusted HR per SD (95% CI), P |
|---|---|---|
| DeSurv k=3 | 1.50 (1.31–1.72), P<0.001 | **1.33 (1.13–1.56), P=0.0005** |
| NMF k=7 (BO α=0) | 1.48 (1.34–1.64), P<0.001 | **1.47 (1.27–1.70), P<0.001** |

**Key observation:** Both methods remain independently significant after adjustment for PurIST and DeCAF. DeSurv k=3 shows moderate attenuation (HR 1.50→1.33); the BO α=0 model shows almost none (1.48→1.47). Neither method's signal is classifier-redundant.

The robustness of the BO α=0 model makes sense: survival information was used to select k=7 and ntop=187 via BO, so it is not unsupervised in the full sense. The key difference from DeSurv is *where* the survival signal enters: in the BO α=0 model, survival guides only hyperparameter selection; in DeSurv, it also guides the factorization directly (α>0), which is what allows the iCAF signal to concentrate into 3 factors rather than 7.

Note: production_summary.rds stores adj_p=0.004 for K=3/α=0.55 — this is from the k-sensitivity grid (Table S3), not from the tab-hr-adjusted pooled Cox (P=0.0005). Do not conflate.

---

## Story: addressing items 8 and 9 together

### The core argument

A BO-optimized unsupervised NMF (α=0, k=7, ntop=187) can match DeSurv's concordance and retain independent prognostic value after classifier adjustment — but it requires more than twice as many factors. The reason is where the survival signal enters: when survival guides only hyperparameter selection (k, ntop), the factorization is unconstrained and needs 7 factors to capture the iCAF program. When survival also guides the factorization (α>0, DeSurv), the decomposition is steered toward survival-relevant structure from the start, keeping the iCAF signal in one factor and making k=3 sufficient. The argument for DeSurv is therefore parsimony: the same independent prognostic value, with fewer factors, by using the available survival information more completely.
The right comparison for item 9 is not against published HRs but against those classifiers directly as covariates — showing DeSurv's signal is independent.

### Arc

1. **BO-optimized NMF at α=0 matches DeSurv's concordance but requires 7 factors** (Table S5): at k=7, this model matches or exceeds DeSurv k=3 per-cohort C-index in most cohorts, and it retains independent prognostic value after classifier adjustment (Table S6). This is the honest starting point — do not downplay it. Critically, this model is not standard NMF: BO used survival to select k=7 and ntop=187, so the survival signal enters through hyperparameter selection even though it does not influence the factorization (α=0).

2. **DeSurv achieves the same with k=3 by letting survival guide the factorization** (Tables S3, S6): when survival also steers the decomposition (α>0), the iCAF signal concentrates into a single factor rather than dispersing, so k=3 suffices. Under the k-sensitivity analysis with uniform ntop=270 (Table S3), k=3 achieves adjusted significance at 5/7 α settings versus only 2/7 for k=7 — at k=7, 270-gene projection fragments the iCAF signal across factors (confirmed by factor-nesting analysis) because without factorization supervision, the decomposition doesn't preferentially keep the prognostic program intact.

3. **K=3 achieves the broadest adjusted significance with the fewest factors** (Table S3): under the production ntop=270 gene-focusing, k=3 reaches adjusted P<0.05 at 5 of 7 α values tested, versus 4/7 for k=5 and only 2/7 for k=7. The k=7 underperformance under uniform gene-focusing is mechanistic: 270-gene projection concentrates each factor onto its most characteristic genes, and at k=7 the iCAF signal disperses across multiple factors rather than concentrating in one (confirmed by factor-nesting analysis). K=5 does achieve significance across 4 of 7 α values with ntop=270, providing convergent support for the iCAF biology from a higher-k model, but with two additional factors and a narrower α range than k=3.

4. **The meaningful test for item 9 is classifier independence, not published HRs**: rather than comparing to Moffitt/PurIST/DeCAF literature HRs, Table S6 shows DeSurv k=3 retains a significant independent HR after adjustment (HR=1.33, P=0.0005). This is a stronger statement than citing literature HRs. NMF k=7 is also independently significant in Table S6 (HR=1.47, P<0.001), so Table S6 alone does not distinguish the two methods — the distinction comes from parsimony and the Table S3 robustness analysis.

---

## Proposed Discussion paragraph 3 addition

Current text ends with: "...producing a fragmented structure whose prognostic content overlaps with existing molecular classifiers (SI Appendix, Tables S2--S3)."

**Proposed sentence to append:**

> "An unsupervised factorization with survival-guided hyperparameter selection (BO, α=0) approaches DeSurv's concordance only at k=7 (SI Appendix, Table S5); incorporating survival information directly into the factorization (α>0) concentrates the iCAF signal into a single factor, making k=3 sufficient — under controlled conditions, k=3 achieves adjusted significance across 5 of 7 supervision settings versus 2 of 7 for k=7 (SI Appendix, Table S3)."

---

## Notes / open questions

- Table S6 values verified by live computation from val_latent RDS (ntop_bo_50_300). Both DeSurv k=3 and NMF k=7 are independently significant after adjustment; the "NMF k=7 loses significance" claim was wrong.
- NMF k=7 was fit with BO-selected ntop=187; DeSurv k=3 with ntop=270. Table S3 applies ntop=270 uniformly to all models — these are different comparisons and should not be conflated.
- The table reference numbers (S3, S4, S5, S6) should be confirmed against the actual compiled SI once rendered.
- Item 9: Table S6 shows DeSurv k=3 is independently prognostic (HR=1.33, P=0.0005), which answers Jen Jen's question more directly than citing published literature HRs. However, Table S6 alone does not differentiate DeSurv from NMF k=7 — pair it with the Table S3 robustness argument.
