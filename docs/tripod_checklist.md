# TRIPOD Checklist — DeSurv (development + external validation)

**TRIPOD 2015** (Collins et al., *Ann Intern Med* 2015;162:55–63). Type: **3 — development and validation in separate data** (TCGA+CPTAC training; Dijk, Moffitt, PACA-AU array, PACA-AU RNA-seq, Puleo external validation).

For Nature Cancer supplement. After the structural rewrite, re-verify section/page references; placeholders below reflect the current `naim/nature-resubmission` source.

| # | Item | D, V, or both | Page / Section | Status |
|---|---|---|---|---|
| **Title and abstract** | | | | |
| 1 | Identify the study as developing and/or validating a multivariable prediction model | D;V | Title; abstract | ⚠️ **Title does not explicitly say "prediction model"** — consider adding language, or accept that "generalizable prognostic program" implies this |
| 2 | Provide a summary of objectives, study design, setting, participants, sample size, predictors, outcome, statistical analysis, results, and conclusions | D;V | Abstract | ✅ Covered in 137-word abstract |
| **Introduction** | | | | |
| 3a | Explain the medical context (including whether diagnostic or prognostic) and rationale | D;V | `02_introduction_REVISED.Rmd` ¶1 | ✅ Prognostic; PDAC context established |
| 3b | Specify the objectives, including whether the study describes development or validation | D;V | `02_introduction_REVISED.Rmd` final paragraph | ✅ Both development and validation stated |
| **Methods — Source of data** | | | | |
| 4a | Describe the study design or source of data (e.g., randomized trial, cohort) | D;V | `03_methods_REVISED.Rmd` §1 | ✅ Public bulk-transcriptomic cohorts |
| 4b | Specify the key study dates, including start of accrual, end of accrual, end of follow-up | D;V | `03_methods_REVISED.Rmd` | ✅ Stated as defined in each cohort's original publication (referenced in SI Section 7); per-cohort date table optional |
| **Methods — Participants** | | | | |
| 5a | Specify key elements of the study setting (e.g., primary care, secondary care, general population) including number and location of centres | D;V | Methods | ✅ TCGA/CPTAC/GEO/ArrayExpress/ICGC — international, public, secondary use |
| 5b | Describe eligibility criteria for participants | D;V | Methods | ✅ Explicit eligibility now stated in Methods (primary PDAC tumors with a bulk expression profile and non-missing OS time/event) |
| 5c | Give details of treatments received, if relevant | D;V | Discussion §limitations | ✅ Treatment-naive; noted as limitation |
| **Methods — Outcome** | | | | |
| 6a | Clearly define the outcome that is predicted by the prediction model, including how and when assessed | D;V | Methods | ✅ Overall survival (time to death, time to censoring) |
| 6b | Report any actions to blind assessment of the outcome to be predicted | D;V | n/a — public data | ✅ N/A (registry-derived outcomes) |
| **Methods — Predictors** | | | | |
| 7a | Clearly define all predictors used in developing or validating the multivariable prediction model, including how and when measured | D;V | Methods | ✅ Bulk transcriptomic gene expression; top-N genes via DeSurv ranking |
| 7b | Report any actions to blind assessment of predictors for the outcome and other predictors | D;V | n/a — algorithmic | ✅ Algorithmic feature ranking; no human-in-loop |
| **Methods — Sample size** | | | | |
| 8 | Explain how the study size was arrived at | D | Methods or SI | ✅ Events-per-parameter justification in Methods (139 training events, ~46/factor for k=3, vs conventional min of 10; 388 validation events) |
| **Methods — Missing data** | | | | |
| 9 | Describe how missing data were handled (e.g., complete-case analysis, single imputation, multiple imputation) with details of any imputation method | D;V | Methods | ✅ Stated: no imputation, complete-case throughout; classifier-adjusted analyses restricted to samples with available calls |
| **Methods — Statistical analysis** | | | | |
| 10a | Describe how predictors were handled in the analyses | D | Methods | ✅ NMF gene-program decomposition with Cox supervision |
| 10b | Specify type of model, all model-building procedures (including any predictor selection), and method for internal validation | D | Methods | ✅ NMF+Cox; Bayesian optimization; 1-SE rule; internal CV |
| 10c | For validation, describe how the predictions were calculated | V | Methods | ✅ Projection of trained W onto validation cohorts; linear-predictor formation |
| 10d | Specify all measures used to assess model performance and, if relevant, to compare multiple models | D;V | Methods | ✅ C-index, pooled stratified-Cox HR per SD, KM dichotomization via CV-optimized cutpoint |
| 10e | Describe any model updating (e.g., recalibration) arising from the validation, if done | V | n/a | ✅ Model is fixed at training; no updating performed |
| **Methods — Risk groups** | | | | |
| 11 | Provide details on how risk groups were created, if done | D;V | Methods + Results | ✅ CV-optimized cutpoint on training; applied to validation |
| **Methods — Development vs. validation** | | | | |
| 12 | For validation, identify any differences from the development data in setting, eligibility criteria, outcome, and predictors | V | SI Table tab:cohorts | ✅ Events (rate) column + caption now surface the dev-vs-val case-mix (train 50-52% event rate vs validation 60-90%, pooled 67%) |
| **Results — Participants** | | | | |
| 13a | Describe the flow of participants through the study, including the number of participants with and without the outcome and, if applicable, a summary of the follow-up time | D;V | Results + SI Table S1 | ✅ Per-cohort n and events reported |
| 13b | Describe the characteristics of the participants (basic demographics, clinical features, available predictors), including the number of participants with missing data for predictors and outcome | D;V | SI Table S1 | ✅ |
| 13c | For validation, show a comparison with the development data of the distribution of important variables (demographics, predictors and outcome) | V | SI Table tab:cohorts | ✅ Per-cohort N, events, event rate, platform, role in one table; outcome-frequency gap stated in caption |
| **Results — Model development** | | | | |
| 14a | Specify the number of participants and outcome events in each analysis | D | Results | ✅ |
| 14b | If done, report the unadjusted association between each candidate predictor and outcome | D | n/a — high-dim | ✅ N/A (gene-program approach rather than univariate marker screen) |
| **Results — Model specification** | | | | |
| 15a | Present the full prediction model to allow predictions for individuals (i.e., all regression coefficients, and model intercept or baseline survival at a given time point) | D | Methods, SI Tables, code repository | ✅ Trained W matrix + Cox coefficients in Zenodo/GitHub repo |
| 15b | Explain how to use the prediction model | D | Methods + README | ✅ Single-sample projection via DeSurv R package |
| **Results — Model performance** | | | | |
| 16 | Report performance measures (with CIs) for the prediction model | D;V | Results + Tables | ✅ Per-cohort C-index, pooled HR with 95% CI |
| **Results — Model-updating** | | | | |
| 17 | If done, report the results from any model updating | V | n/a | ✅ Not performed |
| **Discussion — Limitations** | | | | |
| 18 | Discuss any limitations of the study (such as nonrepresentative sample, few events per predictor, missing data) | D;V | Discussion §limitations | ✅ Treatment-naive cohorts, Cox PH assumption, PDAC-only |
| **Discussion — Interpretation** | | | | |
| 19a | For validation, discuss the results with reference to performance in the development data, and any other validation data | V | Discussion | ✅ Cross-cohort generalization discussed |
| 19b | Give an overall interpretation of the results, considering objectives, limitations, results from similar studies, and other relevant evidence | D;V | Discussion | ✅ |
| **Discussion — Implications** | | | | |
| 20 | Discuss the potential clinical use of the model and implications for future research | D;V | Discussion | ✅ Risk stratification + trial enrichment + neoadjuvant cohort applicability |
| **Other information — Supplementary information** | | | | |
| 21 | Provide information about the availability of supplementary resources, such as study protocol, web calculator, and data sets | D;V | Data and Code Availability | ✅ Zenodo + GitHub |
| **Other information — Funding** | | | | |
| 22 | Give the source of funding and the role of the funders for the present study | D;V | Acknowledgements | ✅ NCI U01/P50/T32/R01 + DOD HT9425241103100121 |

## Items flagged for action

Addressed 2026-07-20 (all in Methods/SI, which are excluded from the Nature Cancer 4,000-word main-text limit):
- **Item 4b (Dates):** ADDRESSED — Methods now states accrual/follow-up periods are those of the primary source studies (SI Section 7). A per-cohort date table could still be added if an editor asks.
- **Item 5b (Eligibility):** DONE — Methods now states explicit eligibility (primary PDAC tumors with an expression profile and non-missing OS time/event).
- **Item 8 (Sample size):** DONE — Methods now gives the events-per-parameter justification (139 training events, ~46 per factor for k=3, vs the conventional minimum of ten; 388 validation events).
- **Item 9 (Missing data):** DONE — Methods now states no imputation; complete-case throughout; classifier-adjusted analyses restricted to samples with available calls.
- **Assumptions / Cox-PH diagnostic:** DONE — SI now reports a cox.zph Schoenfeld-residual test (PH not supported, chisq 14.7, df 1, P 1.3e-4; effect attenuates over follow-up; HR reframed as time-averaged, complemented by PH-free IBS/IPA). **PI review needed on the PH-violation framing.**
- **Per-cohort QC flow (13a traceability):** DONE — Methods now gives the post-QC per-cohort training split (TCGA-PAAD 144, CPTAC-3 129) so 181->144 and 140->129 are traceable.

Still open (judgment/lower priority):
- **Item 1 (Title):** does not say "prediction model" — judgment call.
- **Item 12 / 13c (Dev-vs-validation differences + distributional comparison):** SI Table `tab:cohorts` gives per-cohort platform/N/events/role; a fuller predictor/outcome distributional comparison table is not present. Author call on whether to add.
- **PACA-AU array accession (P2-6):** Amber's item.
