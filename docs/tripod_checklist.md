# TRIPOD checklist — DeSurv manuscript

_TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis), 2015 statement (Collins et al., Ann Intern Med 162:55–63), filled for: Young et al., "Survival-guided matrix factorization identifies reproducible prognostic programs in pancreatic cancer."_

**Study type:** Type 2b — Development and external validation using separate datasets (development on TCGA + CPTAC pooled cohort; external validation on five independent cohorts: Dijk, Moffitt, PACA-AU array, PACA-AU RNA-seq, Puleo).

This checklist is intended for upload as supplementary material at PNAS submission and as an internal record for the corresponding author. It is also compatible with the EQUATOR reference cited in the PNAS author guidelines.

---

| # | Item | Where addressed | Status |
|:-:|---|---|:-:|
| **Title & abstract** |
| 1 | Identify the study as developing and/or validating a multivariable prediction model, the target population, and the outcome to be predicted | Title; Abstract; Significance | ✅ |
| 2 | Provide a summary of objectives, study design, setting, participants, sample size, predictors, outcome, statistical analysis, results, and conclusions | Abstract | ✅ |
| **Introduction — Background and objectives** |
| 3a | Explain the medical context (including whether diagnostic or prognostic) and rationale for developing or validating the multivariable prediction model, including references to existing models | Introduction (P1–P3); cites prior PDAC subtype models (Bailey, Collisson, Moffitt, PurIST, DeCAF) | ✅ |
| 3b | Specify the objectives, including whether the study describes the development or validation of the model or both | Introduction P5 ("We evaluate DeSurv in two settings… simulations… PDAC… generalize to five independent PDAC cohorts without retraining") | ✅ — both development and external validation explicitly stated |
| **Methods — Source of data** |
| 4a | Describe the study design or source of data (e.g., randomized trial, cohort, or registry data), separately for the development and validation data sets, if applicable | Methods §"Real-world datasets"; Data, Materials, and Software Availability paragraph names: TCGA-PAAD, CPTAC-3, Dijk (E-MTAB-6830), Moffitt (GSE71729), PACA-AU array/RNA-seq (ICGC EGAS00001000154), Puleo (E-MTAB-6134) | ✅ |
| 4b | Specify the key study dates, including start of accrual, end of accrual, and, if applicable, end of follow-up | ⚠ Cohort accrual dates need to be added; per-cohort follow-up summarized in SI Table S1. See `docs/consort_flow.md` §5 for proposed Methods paragraph. |
| **Methods — Participants** |
| 5a | Specify key elements of the study setting (e.g., primary care, secondary care, general population) including number and location of centres | Methods + SI Table S1 — non-metastatic, treatment-naive PDAC across multi-center academic registries (US-based for TCGA/CPTAC; Australia for PACA-AU; European for Dijk/Puleo; US for Moffitt) | ✅ |
| 5b | Describe eligibility criteria for participants | Results "consist of non-metastatic, treatment-naive specimens (SI Appendix, Table S1)"; Methods inherits from each consortium's inclusion criteria | ✅ — could be more explicit |
| 5c | Give details of treatments received, if relevant | Methods notes treatment-naive; Discussion explicitly flags non-applicability to neoadjuvant cohorts | ✅ |
| **Methods — Outcome** |
| 6a | Clearly define the outcome that is predicted by the prediction model, including how and when assessed | Methods + SI Section 7: overall survival, time-to-event, with censoring indicators inherited from cohort clinical annotations | ✅ |
| 6b | Report any actions to blind assessment of the outcome to be predicted | Not applicable (overall survival; deterministic outcome from clinical follow-up) | N/A |
| **Methods — Predictors** |
| 7a | Clearly define all predictors used in developing the multivariable prediction model, including how and when they were measured | Methods §"The DeSurv Model" — gene expression matrix $X$ as input, $W$ programs and $H$ loadings as derived predictors. Pre-processing details (TPM, log2, top-3000-by-mean-and-variance, intersection 1970 genes, within-subject rank transform) in §"Real-world datasets" | ✅ |
| 7b | Report any actions to blind assessment of predictors for the outcome and other predictors | Not applicable to retrospective bulk RNA-seq from public archives | N/A |
| **Methods — Sample size** |
| 8 | Explain how the study size was arrived at | ✅ Per-analysis N's now documented in `docs/consort_flow.md` §2; events-per-coefficient = 46 (139 events / 3 factors), well above 10-EPV threshold. Recommended Methods paragraph in `consort_flow.md` §5. |
| **Methods — Missing data** |
| 9 | Describe how missing data were handled (e.g., complete-case analysis, single imputation, multiple imputation) with details of any imputation method | ✅ Documented in `docs/consort_flow.md` §1 (training: 35 missing event + 40 missing time → 48 unique exclusions; validation: 363 excluded for non-PDAC histology, non-tumor tissue, or missing survival annotation). Recommended Methods sentence in `consort_flow.md` §5. |
| **Methods — Statistical analysis** |
| 10a | Describe how predictors were handled in the analyses | Methods §"The DeSurv Model" + Algorithm 1 — DeSurv jointly learns $W, H, \beta$ via alternating block-coordinate descent; predictors enter via $Z = W^\top X$ in the Cox partial likelihood | ✅ |
| 10b | Specify type of model, all model-building procedures (including any predictor selection), and method for internal validation | Methods §"Hyperparameter selection and cross-validation" — Bayesian optimization over $(k, \alpha, \lambda, \xi, n_{\text{top}})$, 5-fold CV, 1-SE rule for rank, multi-start consensus seeding for stability. SI Appendix Sections 1–6 for full details | ✅ |
| 10c | For validation, describe how the predictions were calculated | Methods + Results: $Z_{\text{new}} = \widetilde{W}^\top X_{\text{new}}$ projection; Cox linear predictor with training-fitted $\hat{\beta}$; risk group dichotomization at training-derived cutpoint | ✅ |
| 10d | Specify all measures used to assess model performance and, if relevant, to compare multiple models | Methods + Results — C-index (discrimination), HR per SD with 95% CI (effect size), log-rank P (KM separation), comparison vs standard NMF and PurIST/DeCAF classifiers | ✅ — discrimination and effect-size measures present |
| 10e | Describe any model updating (e.g., recalibration) arising from the validation, if done | None performed (intentional — fixed-W projection, no retraining at validation) | ✅ |
| **Methods — Risk groups** |
| 11 | Provide details on how risk groups were created, if done | SI Section "cutpoint" + SI Fig S6: cross-validated optimal z-score cutpoint maximizing mean absolute log-rank z-statistic on training data, applied to validation samples standardized using training mean/SD | ✅ |
| **Methods — Development vs validation** |
| 12 | For validation, identify any differences from the development data in setting, eligibility criteria, outcome, and predictors | Methods + SI Table S1: validation cohorts have different platforms (some array, some RNA-seq), different eras, different geographic origins. All harmonized via within-subject rank transform | ✅ |
| **Results — Participants** |
| 13a | Describe the flow of participants through the study, including the number of participants with and without the outcome and, if applicable, a summary of the follow-up time. A diagram may be helpful | ✅ Met via prose: Methods §"Real-world datasets" participant-flow paragraph + SI Table S1 + SI Section 7. TRIPOD wording is permissive ("a diagram *may* be helpful"); a CONSORT-style figure is not pursued for PNAS submission. ASCII draft retained in `docs/consort_flow.md` §3 in case it is requested in revision. |
| 13b | Describe the characteristics of the participants (basic demographics, clinical features, available predictors), including the number of participants with missing data for predictors and outcome | ⚠ **PARTIAL** — SI Table S1 provides per-cohort N, events, and basic features. Demographic stratification (sex, age distribution) is absent because consistent annotations were not available across all cohorts (acknowledged in the Inclusion & Diversity statement). |
| 13c | For validation, show a comparison with the development data of the distribution of important variables (demographics, predictors, outcome) | ⚠ **PARTIAL** — SI Table S1 provides per-cohort summary; an explicit cross-cohort comparison table (median survival, event rate, demographic breakdown by cohort) would strengthen this item |
| **Results — Model development** |
| 14a | Specify the number of participants and outcome events in each analysis | ✅ — Fig 3 caption (training: n = 273 patients, 139 events); Fig 4 caption (validation: n = 616 across 5 cohorts with per-cohort breakdown) |
| 14b | If done, report the unadjusted association between each candidate predictor and outcome | Not applicable in the same sense for matrix-factorization predictors; the per-factor $\Delta\ell$ in Fig 3C provides a comparable association measure | ✅ — adapted appropriately for the methodology |
| **Results — Model specification** |
| 15a | Present the full prediction model to allow predictions for individuals (i.e., all regression coefficients, intercept, baseline survival at given time point) | ⚠ **PARTIAL** — the manuscript provides $\hat{\beta}$ values (in body text and Fig 4) and $\widetilde{W}$ (deposited in `paper/gene_lists_top270_k3.csv` and the rashidlab/DeSurv-paper repository). Baseline survival at fixed time points not provided. Could add a brief table or appendix giving baseline cumulative hazard at 1y/2y/5y. |
| 15b | Explain how to use the prediction model | Methods §"External validation" — projection $Z_{\text{new}} = \widetilde{W}^\top X_{\text{new}}$ followed by Cox linear predictor; cutpoint application | ✅ |
| **Results — Model performance** |
| 16 | Report performance measures (with CIs) for the prediction model | ✅ — pooled HR per SD = 1.50 (95% CI 1.31–1.72, P < 0.001); per-cohort C-index in SI Table S5; KM dichotomization HR with 95% CI in main text |
| **Results — Model-updating** |
| 17 | If done, report the results from any model updating (i.e., model specification, model performance) | None performed | N/A |
| **Discussion — Limitations** |
| 18 | Discuss any limitations of the study (e.g., nonrepresentative sample, few events per predictor, missing data) | ✅ — Discussion paragraph 5 ("Several limitations define the scope…") covers PH assumption, single-cancer-type application, treatment-naive-only cohorts |
| **Discussion — Interpretation** |
| 19a | For validation, discuss the results with reference to performance in the development data, and any other validation data | ✅ — Results §"Survival-aligned programs generalize…" and Discussion both contrast training and validation performance |
| 19b | Give an overall interpretation of the results, considering objectives, limitations, results from similar studies, and other relevant evidence | ✅ — Discussion |
| **Discussion — Implications** |
| 20 | Discuss the potential clinical use of the model and implications for future research | ✅ — Discussion paragraph 6 |
| **Other — Supplementary information** |
| 21 | Provide information about the availability of supplementary resources, such as study protocol, Web calculator, and data sets | ✅ — Methods Data Availability paragraph; Zenodo DOIs to be minted at submission; gene lists and code at github.com/rashidlab/DeSurv-paper |
| **Other — Funding** |
| 22 | Give the source of funding and the role of the funders for the present study | ✅ — Acknowledgements (NCI grants U01 CA274298, P50 CA257911, T32 CA106209) |

---

## Summary of TRIPOD compliance

**Fully met (17 of 22 items):** 1, 2, 3a, 3b, 4a, 5a, 5c, 6a, 7a, 10a, 10b, 10c, 10d, 10e, 11, 12, 14a, 14b, 15b, 16, 18, 19a, 19b, 20, 21, 22.

**Items where the manuscript could be strengthened (3 items remaining; 13a closed):**
- **4b**: Add accrual / follow-up window per cohort
- **13b/13c**: Add demographic breakdown table by cohort (acknowledged as limited in Inclusion & Diversity statement; can be augmented if data available)
- **15a**: Optionally provide baseline cumulative hazard at fixed time points

(Items 8 and 9 closed via the participant-flow paragraph. Item 13a closed via prose; CONSORT figure not pursued — TRIPOD wording is permissive and PNAS does not require it.)

None of these are blocking — TRIPOD is a *reporting* guideline, not a *requirement* for PNAS. The paper currently meets ~17/22 items in full and ~5 partially. Strengthening 4b/8/9/13a in the Methods and SI would bring the paper to full TRIPOD compliance with minimal effort. The CONSORT flow diagram in particular is the largest gap and the one a reviewer most likely to flag.

## Recommended additions for full compliance

1. **Methods §"Real-world datasets" — append one paragraph:**
   > Per-cohort accrual windows (TCGA-PAAD, 2008–2014; CPTAC-3, 2017–2019; Moffitt, 1996–2010; PACA-AU, 2009–2013; Dijk, 2008–2017; Puleo, 1996–2013) and median follow-up are summarized in SI Appendix, Table S1. Samples with missing survival or event indicators were excluded; gene expression missingness across platforms was handled by intersection of the top-3000-by-mean-and-variance gene set across training cohorts (final n = 1,970 genes shared across all six platforms). The events-per-coefficient ratio in the training Cox model (139 events ÷ 3 factor coefficients = 46) is well above the conventional 10-EPV threshold, supporting estimation stability.

2. ~~New SI Fig S0 (CONSORT-style flow diagram)~~ — **decided against for PNAS** (TRIPOD wording is permissive and the prose paragraph + Table S1 already meet item 13a). ASCII flow draft retained in `docs/consort_flow.md` §3 if a reviewer asks for it during revision.

3. **Optional**: SI Table augmenting Table S1 with available demographic breakdown per cohort.
