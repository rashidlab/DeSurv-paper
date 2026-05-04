# CONSORT-style flow diagram and analysis map

_Per-cohort participant flow and analysis-to-cohort mapping for the DeSurv manuscript. All N's extracted directly from the cached `.rds` objects in `results/precomputed/` and `data/original/`. **Status:** the recommended Methods paragraph in §5 has been added to `03_methods_REVISED.Rmd`; a CONSORT-style flow figure (originally proposed as SI Fig S0) was decided against for PNAS submission (TRIPOD is permissive, prose paragraph + Table S1 are sufficient). This document is retained as the reference draft if a reviewer asks for the figure during revision._

## 1. Participant flow per cohort

### Training cohorts (development)

| Cohort | Source | Raw N (in dataset) | Excluded for missing survival or QC fail | Analytic N | Events |
|---|---|---:|---:|---:|---:|
| TCGA-PAAD | NCI Genomic Data Commons | 181 | 37 | 144 | 75 |
| CPTAC-3 | NCI Proteogenomic Translational Research Centers | 140 | 11 | 129 | 64 |
| **Pooled training** | TCGA + CPTAC | **321** | **48** | **273** | **139** |

Of the 48 excluded training samples, 35 had missing event indicators and 40 had missing follow-up time (with overlap; net 48 unique exclusions). Inclusion was implemented via the `keep` flag in the cohort metadata, set during data loading (`R/load_data_internal.R`).

### Validation cohorts (external)

| Cohort | Source / accession | Raw N (in dataset) | Excluded (non-PDAC, non-tumor, missing survival) | Analytic N | Events |
|---|---|---:|---:|---:|---:|
| Dijk | ArrayExpress E-MTAB-6830 | 90 | 0 | 90 | 81 |
| Moffitt GEO array | NCBI GEO GSE71729 | 357 | 234 | 123 | 83 |
| PACA-AU array | ICGC EGAS00001000154 | 131 | 68 | 63 | 38 |
| PACA-AU RNA-seq | ICGC EGAS00001000154 | 92 | 40 | 52 | 31 |
| Puleo array | ArrayExpress E-MTAB-6134 | 309 | 21 | 288 | 181 |
| **Pooled validation** | 5 cohorts | **979** | **363** | **616** | **414** |

The large Moffitt exclusion (234 of 357) reflects the original dataset's inclusion of normal pancreas, chronic pancreatitis, and metastatic samples — only primary, treatment-naive PDAC tumors with survival annotation are retained. The PACA-AU exclusions reflect the same kind of filter plus duplicate samples between array and RNA-seq submissions.

### Gene/feature flow

| Step | Genes |
|---|---:|
| Raw expression matrices (per cohort, varies) | varies (~20,000–60,000 transcripts) |
| Filtered to top 3,000 by mean+variance ranking, **per training cohort** | 3,000 each |
| Intersection across TCGA + CPTAC training cohorts (used for all analyses) | **1,970** |
| Validation cohorts restricted to this same 1,970-gene set | 1,970 |
| `ntop` per-factor projection (Bayesian-optimized) | **270 per factor** |

## 2. Analysis-to-cohort mapping

Each panel/result in the manuscript uses a specific subset; this table shows which cohorts contributed to which analysis and the N actually entered.

| Analysis | Figure / table | Cohorts entering | N | Events |
|---|---|---|---:|---:|
| Simulation studies (3 scenarios × 100 reps) | Fig 2A,B,C; SI Fig S2 | synthetic (no human data) | per scenario p=3,000 genes, n=200 samples, k=3 factors, 100 replicates | per scenario |
| Bayesian optimization (5-fold CV) for $(k, \alpha, \lambda, \xi, n_{\text{top}})$ | Fig 2D | TCGA + CPTAC pooled training | 273 | 139 |
| Final multi-start consensus fit at BO-selected hyperparameters | Fig 3A,B,C,D | TCGA + CPTAC pooled training | 273 | 139 |
| External validation, per-cohort C-index + log-rank | Fig 4A; SI Table S5 | Dijk, Moffitt, PACA-AU array, PACA-AU seq, Puleo | per cohort: 90, 123, 63, 52, 288 | per cohort: 81, 83, 38, 31, 181 |
| Pooled validation, stratified Cox per SD of linear predictor | Fig 4 caption + Discussion | 5 validation cohorts pooled | 616 | 414 |
| Dichotomized risk-group KM (cutpoint cross-validated on training) | Fig 4B,C | training cutpoint selection: 273; validation application: 616 | — | — |
| Adjustment for PurIST + DeCAF molecular classifiers | SI Table S6 | validation cohorts where PurIST/DeCAF annotations are available (subset of 616) | subset of 616 | subset of 414 |
| $k \times \alpha$ sensitivity analysis | SI Tables S2–S4 | 5 validation cohorts pooled (with covariate-adjusted Cox) | 616 (subset where PurIST/DeCAF available, for adjusted) | varies by combination |
| Subtype composition of risk groups | SI Fig S7 | validation cohorts pooled | 616 | — |
| C-index by rank $k$ analysis | SI Fig S3 | training (CV) + validation (pooled) | 273 / 616 | 139 / 414 |
| NMF rank-selection diagnostics (cophenetic, silhouette, residual) | SI Fig S4 | TCGA + CPTAC pooled training | 273 | 139 |
| $k=3$ vs $k=7$ factor correspondence | SI Fig S5 | training | 273 | 139 |
| Convergence trajectories | SI Fig S1 | training | 273 | 139 |

## 3. Draft flow-diagram structure (for SI Fig S0)

The recommended SI figure should depict the participant flow as a vertical waterfall, with a left panel for training and a right panel for validation. Suggested boxes (top → bottom):

```
TRAINING                                      VALIDATION
─────────                                     ──────────

TCGA-PAAD          CPTAC-3                    Dijk    Moffitt   PACA-AU   Puleo
n = 181            n = 140                    n = 90  n = 357   n = 223*  n = 309
   │                  │                          │       │         │        │
   │ excluded:        │ excluded:                │       │ exclude │        │
   │  37 (missing     │  11 (missing             │       │  234    │ 108    │  21
   │   survival/QC)   │   survival/QC)           │       │ (non-   │ (non-  │
   ▼                  ▼                          ▼       ▼  PDAC  ▼  PDAC, ▼
n = 144            n = 129                    n = 90  n = 123   PACA  n = 288
events: 75         events: 64                 evt=81  evt=83    arr+seq evt=181
   │                  │                          │       │      n=63+52    │
   └────────┬─────────┘                          │       │      evt=38+31  │
            ▼                                    └───────┴─────────┴───────┘
   POOLED TRAINING                                          │
   n = 273, events = 139                                    ▼
   genes after filtering: 1,970                       POOLED VALIDATION
            │                                       n = 616, events = 414
            ▼
   Bayesian optimization (5-fold CV)
   k, α, λ, ξ, ntop selected
            │
            ▼
   Final consensus fit at k = 3, α = 0.334
   ntop = 270 per factor
            │
            └──────────► W projected to validation: Z = W^T X
                                  │
                                  ▼
                        Per-cohort C-index, HR
                        Pooled HR per SD = 1.50
                        (95% CI 1.31–1.72)
```

*PACA-AU is split into array (n=131) and RNA-seq (n=92) submissions in the original ICGC archive; combined raw is 223. Both are reported separately as analytic cohorts.

## 4. Where this should land

**Status (2026-05-03):**
- The recommended Methods paragraph in §5 below has been added to `03_methods_REVISED.Rmd` and references `Table S1` (which exists). The prose paragraph contains all per-cohort participant-flow N's.
- The CONSORT-style figure (originally proposed as SI Fig S0) is **not yet drawn**. The Methods prose paragraph closes TRIPOD item 13a textually, but a visual flow diagram is still recommended for full visual TRIPOD compliance and as defensive against potential reviewer requests.
- If you decide to draw the figure: render as a TikZ flowchart or Inkscape diagram and insert in SI. Because main text no longer references "Fig. S0", you can place it at the end of the SI as the next available number (would currently be S10 since SI has S1–S9), or insert at the beginning and renumber existing S1–S9 to S2–S10. The SI numbering uses `\renewcommand{\thefigure}{S\arabic{figure}}` so insertion is automatic.

## 5. Recommended Methods §"Real-world datasets" addition

> Per-cohort participant flow is summarized in SI Fig S0 (also see SI Appendix, Section 7 and Table S1). Of 321 pooled training-cohort samples (TCGA-PAAD, n = 181; CPTAC-3, n = 140), 48 were excluded for missing survival time, missing event indicator, or quality-control failure, yielding 273 analytic training samples (139 events). Of 979 pooled validation-cohort samples, 363 were excluded for non-PDAC histology, non-tumor tissue, or missing survival annotation, yielding 616 analytic validation samples (414 events) across five cohorts: Dijk (n = 90, 81 events), Moffitt GEO array (n = 123, 83), PACA-AU array (n = 63, 38), PACA-AU RNA-seq (n = 52, 31), and Puleo array (n = 288, 181). The events-per-coefficient ratio in the training Cox model (139 ÷ 3 = 46) is well above the conventional 10-EPV threshold, supporting estimation stability.

This addition would close TRIPOD items 4b, 8, 9, 13a, and 14a in a single Methods paragraph, plus the SI figure.
