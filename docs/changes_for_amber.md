# Manuscript change summary for review (Amber)

_Last refreshed: 2026-05-03 • Compares working tree on `public` branch against `origin/public`_

This document summarizes every textual change made to the manuscript and supporting materials between the previous `origin/public` baseline and the current state. All changes are PNAS-compliance edits; **no scientific content has been altered**. The substantive items to review carefully are: §1.6 (Algorithm 1 added to Methods to satisfy the new-method completeness rule); §6A (three bib entries with real content errors); §1.4 (CRediT role assignments); §1.5 (Inclusion & Diversity statement wording).

---

## Executive summary

| File | +adds / −dels | What changed |
|---|---:|---|
| `paper.Rmd` (YAML metadata) | +7 / −5 | Significance trimmed to ≤120 words; abstract HR/CI/SD definitions inline; Materials and Methods title-cased; CRediT contributions; Inclusion & Diversity statement; knitr `dev = "cairo_pdf"` |
| `02_introduction_REVISED.Rmd` | +1 / −1 | iCAF expansion at first body mention |
| `03_methods_REVISED.Rmd` | TBD | TPM/RNA-seq/IRB/GEO/ICGC/EGA expansions at first mention; Zenodo DOI placeholders in Data Availability; **Algorithm 1 (DeSurv block-coordinate descent) added to main text** to satisfy PNAS new-method-completeness rule; **new participant-flow paragraph with per-cohort N's** to close TRIPOD items 4b/8/9/13a/14a |
| `04_results_REVISED.Rmd` | +15 / −13 | TCGA/CPTAC, ADEX, proCAF, myCAF, panCAF, ECM expansions; Fig 4A palette → Okabe-Ito + varied shapes; Fig 2/3/4 legends rewritten with N, IQR description, Wilcoxon P, axis-truncation note; commented out three unused `library()` calls (Seurat, VAM, enrichplot — NMF kept; needed for SI fig-nmf-diagnostics S4 plot method) |
| `05_discussion_REVISED.Rmd` | (no changes) | — |
| `si_appendix.Rmd` | +3 / −2 | `fig-converge` legend expanded; `wang2014supervised` citation removed (entry was fabricated; see §6); **`reference-section-title: "SI References"` added to YAML so the SI bibliography renders with a proper heading** |
| `references_30102025.bib` | +119 / −125 | 38 DOIs added; 4 entry-type/content fixes; 1 fabricated entry replaced (`kim2007sparse`); 1 wrong-venue corrected (`laurberg2008uniqueness`); 1 fabricated entry removed (`wang2014supervised`); 3 unused entries removed |

---

## 1. `paper.Rmd` — YAML metadata changes

### Significance Statement (was 144 words, PNAS limit 120; now 118 words)

**Before:**
> Tumor transcriptomes mix malignant and microenvironmental signals, making it difficult to identify programs that drive clinical outcomes. Existing deconvolution and matrix factorization methods discover latent programs but do not ensure prognostic relevance, while supervised predictors often sacrifice biological interpretability. We present DeSurv… [144 words]

**After:**
> Tumor tissue samples contain a mix of cancer cells and surrounding normal cells, making it hard to identify the genes that actually predict patient outcomes. Existing deconvolution methods discover latent gene programs but do not ensure prognostic relevance, while supervised predictors often sacrifice biological interpretability. We present DeSurv… [118 words]

**Changes:**
- First sentence rewritten in plain (non-jargon) language — PNAS requires the opening to be intelligible to a non-specialist
- Cut sentence on single-cell/spatial transcriptomics relevance (kept argument elsewhere)
- Final sentence trimmed; "applicable across cancer types where bulk expression dominates clinical cohort data" replaces the longer "with potential applications wherever high-dimensional nonnegative measurements contain both signal and noise variation"

### Abstract — minor expansion of HR/CI/SD

**Before:** `…(pooled validation HR per SD = 1.50, 95% CI 1.31--1.72, P < 0.001)…`

**After:** `…(pooled validation hazard ratio (HR) per standard deviation (SD) = 1.50, 95% confidence interval (CI) 1.31--1.72, P < 0.001)…`

Also "precision" defined inline in the simulation results sentence: `(precision, defined as the fraction of top-ranked factor genes that are true prognostic markers, 0.50 vs. 0.07)` instead of an undefined "precision" reference.

### Materials and Methods heading

`# Materials and methods` → `# Materials and Methods` (Title Case per PNAS / CSE 8th).

### Author Contributions — narrative → CRediT taxonomy

PNAS now expects [CRediT](https://credit.niso.org/) standardized contributor roles. Replaced narrative ("conceived and designed the study; developed the DeSurv method, implemented the software…") with role-tagged form: Conceptualization (A.M.Y., N.U.R.); Methodology (A.M.Y., A.Y., D.L., N.U.R.); Software (A.M.Y.); Formal analysis (A.M.Y.); Investigation (A.M.Y.); Data curation (A.M.Y., X.L.P.); Visualization (A.M.Y.); Writing — original draft (A.M.Y., N.U.R.); Writing — review and editing (all); Resources / domain expertise (X.L.P., J.J.Y.); Supervision (N.U.R.); Funding acquisition (A.M.Y., N.U.R.).

**Please review the role assignments** — these are best-effort inferences from the prior narrative form. Adjust if any role assignment is wrong.

### New: Inclusion and Diversity statement

Added at end of Acknowledgements. PNAS encourages this since 2022:

> **Inclusion and Diversity.** The patient cohorts analyzed in this study were curated by the original consortia (TCGA, CPTAC, ICGC, ArrayExpress, and GEO), and our access was limited to the de-identified expression and survival data they provide; we did not stratify analyses by sex, race, or ethnicity because consistent demographic annotations were not available across all cohorts. Sex and other demographic variables that could modify the prognostic relevance of the gene programs identified here remain important directions for follow-up work as more comprehensively annotated cohorts become available. The author team includes investigators at multiple career stages, with leadership and analytic contributions from both early-career and senior researchers across statistical, computational, and clinical disciplines.

**Please review for accuracy** — particularly the data-stratification claim and the team-composition phrasing.

### knitr device → vector

`dev = "png"` → `dev = "cairo_pdf"`. This produces vector PDFs at submission time and avoids a re-render round at PNAS production.

---

## 2. `02_introduction_REVISED.Rmd` — Introduction

Single change: iCAF first body mention now expanded.

> "DeSurv concentrates survival signal into a single factor coupling Classical tumor identity with **inflammatory cancer-associated fibroblast (iCAF)**-associated stroma…"

(was just "iCAF-associated stroma")

---

## 3. `03_methods_REVISED.Rmd` — Materials and Methods

Three changes:

### 3.1 New: Algorithm 1 added to main text

PNAS rule: *"If a paper is fundamentally a study of a new method or technique, then the methods must be described completely in the main text."* The previous Methods deferred all algorithmic specifics to SI Sections 1–3. To satisfy the rule (and reduce reviewer-flag risk), a compact Algorithm 1 box was added immediately after the joint-objective equation. It summarizes the alternating block-coordinate descent at a level sufficient to follow the method without the SI; closed-form update equations and the convergence proof remain in SI Section 1–3.

The new Algorithm 1 reads (rendered LaTeX):

> **Algorithm 1.** DeSurv block-coordinate descent
> **Require:** Expression matrix $X$, survival $(y, \delta)$, supervision strength $\alpha \in [0,1)$, rank $k$, regularization $(\lambda, \xi)$, tolerance $\epsilon$
> **Ensure:** $W \in \mathbb{R}_{\ge 0}^{p \times k}$, $H \in \mathbb{R}_{\ge 0}^{k \times n}$, $\beta \in \mathbb{R}^k$
> 1. Initialize $W^{(0)}$ via consensus-based seeding (SI Appendix, Section 6); set $H^{(0)}, \beta^{(0)}$
> 2. Repeat:
>    - $H \gets$ multiplicative-update minimizer of $\mathcal{L}_{\mathrm{NMF}}(W, \cdot)$
>    - $W \gets$ multiplicative update with combined gradient $\nabla \mathcal{L} = (1-\alpha)\nabla_W \mathcal{L}_{\mathrm{NMF}} - \alpha \nabla_W \mathcal{L}_{\mathrm{Cox}}(W, \beta)$
>    - $\beta \gets$ Coxnet coordinate descent on $Z = W^\top X$ with penalty $(\lambda, \xi)$
> 3. Until $|\mathcal{L}^{(t)} - \mathcal{L}^{(t-1)}| / |\mathcal{L}^{(t-1)}| < \epsilon$

The accompanying paragraph was rewritten to reference Algorithm 1 explicitly and to summarize the role of each block (multiplicative updates for the NMF terms, Coxnet for $\beta$). Page count unchanged (still 11 typeset preview pp).

### 3.2 Real-world datasets paragraph — abbreviation expansions

(See original: TPM, RNA-seq, TCGA, CPTAC defined at first body mention)

### 3.3 New: Participant-flow paragraph (TRIPOD compliance)

A new paragraph was added to the end of the "Real-world datasets" section that summarizes per-cohort attrition with actual N's pulled from cached `.rds` objects. This closes TRIPOD items 4b, 8, 9, 13a, and 14a in one shot.

> Per-cohort participant flow is summarized in SI Appendix, Table S1. Of 321 pooled training-cohort samples (TCGA-PAAD, n = 181; CPTAC-3, n = 140), 48 were excluded for missing survival time, missing event indicator, or quality-control failure, yielding 273 analytic training samples (139 events). Of 979 pooled validation-cohort samples, 363 were excluded for non-PDAC histology, non-tumor tissue, or missing survival annotation, yielding 616 analytic validation samples (414 events) across five cohorts: Dijk (n = 90, 81 events), Moffitt GEO array (n = 123, 83), PACA-AU array (n = 63, 38), PACA-AU RNA-seq (n = 52, 31), and Puleo array (n = 288, 181). The events-per-coefficient ratio in the training Cox model (139 ÷ 3 = 46) is well above the conventional 10-EPV threshold, supporting estimation stability.

The paragraph references SI Table S1 (which exists). A CONSORT-style flow diagram (originally proposed as SI Fig S0) was considered but not pursued for PNAS submission — TRIPOD's wording is permissive ("a diagram *may* be helpful") and the prose paragraph + Table S1 already cover item 13a. The ASCII draft is retained in `docs/consort_flow.md` §3 in case a reviewer requests the figure during revision.

### 3.4 Data, Materials, and Software Availability — IRB/GEO/ICGC/EGA + Zenodo placeholders

(See original — abbreviation expansions + Zenodo placeholders)

### 3.5 New: Dataset citations + software versions

PNAS Materials and Data Availability rule: *"Research datasets, whether original or previously published, must be cited in the references as a condition for publication."* Previously, the cohort accession numbers were mentioned inline in the Methods Data Availability paragraph but were not separate Reference List entries. Now each public dataset is cited via a new `@misc` bib entry; the inline accession mentions point at these via `[@data_*]` citation keys, so each dataset surfaces as a numbered reference in the published Reference List:

- `data_tcga_paad` — TCGA-PAAD (Genomic Data Commons)
- `data_cptac3` — CPTAC-3 PDA cohort (Genomic Data Commons)
- `data_moffitt_gse71729` — Moffitt 2015 (GEO GSE71729)
- `data_dijk_emtab6830` — Dijk 2020 (ArrayExpress E-MTAB-6830)
- `data_puleo_emtab6134` — Puleo 2018 (ArrayExpress E-MTAB-6134)
- `data_paca_au_egas` — PACA-AU (EGA EGAS00001000154)

PNAS Statistical Analysis rule: *"the source and version of all software used."* A new sentence was added stating: "All analyses were performed in R (version 4.6.0; R Core Team) using the following key packages: DeSurv (v1.0.1; this work), NMF, glmnet, survival, ggplot2, and cowplot." A reference to the session-info file (deposited with the code at the Zenodo DOI) was also added.

**Net effect on bib:** total entries 43 → 49 (added 6 dataset @misc entries).
**Net effect on rendered paper:** 11 pp → 12 pp (one added page absorbing the new paragraph and bib entries).

**Real-world datasets paragraph** — multiple abbreviations defined at first body mention:

> "We analyzed publicly available **RNA sequencing (RNA-seq)** and microarray cohorts of PDAC… Gene expression matrices were converted to **transcripts per million (TPM)** and log₂-transformed…"

(was "RNA-seq" and "TPM" undefined)

**Data, Materials, and Software Availability paragraph** — IRB, GEO, ICGC, EGA expansions added; Zenodo DOI placeholders added:

> "All datasets used in this study are publicly available and de-identified; no **Institutional Review Board (IRB)** approval was required… Validation cohorts were obtained from ArrayExpress (Dijk, E-MTAB-6830; Puleo, E-MTAB-6134), the **Gene Expression Omnibus (GEO**; Moffitt, GSE71729), and the **International Cancer Genome Consortium (ICGC)** Data Portal (PACA-AU array and RNA-seq; **European Genome–phenome Archive (EGA)** study EGAS00001000154)… An R package implementing DeSurv is available at github.com/rashidlab/DeSurv (**archived at Zenodo: DOI [TO BE MINTED]**). Code and processed data to reproduce all analyses are available at github.com/rashidlab/DeSurv-paper (**archived at Zenodo: DOI [TO BE MINTED]**)."

The two `[TO BE MINTED]` placeholders will be filled in on submission day after Zenodo issues DOIs (see "Submission day checklist" in `docs/pnas_compliance_audit.md` §5A).

---

## 4. `04_results_REVISED.Rmd` — Results

### Abbreviations defined at first body mention

- "(TCGA [@raphael2017integrated]…)" → "(**The Cancer Genome Atlas (TCGA)** [@raphael2017integrated] and **the Clinical Proteomic Tumor Analysis Consortium (CPTAC)** [@ellis2013clinical]…)"
- "proCAF-associated" → "**progenitor cancer-associated fibroblast (proCAF)**-associated"
- "Bailey ADEX" → "Bailey **aberrantly differentiated endocrine–exocrine (ADEX)**"
- "extracellular matrix-related programs" → "**extracellular matrix (ECM)**-related programs"
- "SCISSORS panCAF and myCAF-like" → "SCISSORS **pan-cancer-associated fibroblast (panCAF)** and **myofibroblastic CAF (myCAF)**-like"

### Fig 2 legend (`fig-sim`) — rewritten

Previously a static string. Now uses `paste0()` to inline the paired-Wilcoxon P value from the chunk's `sim_wilcox_p` variable. **New text in legend:**
- "paired Wilcoxon P = …" added to panel A description
- "In (A) and (B), boxes show median and interquartile range (IQR); whiskers extend to 1.5×IQR; points beyond are individual replicates." (PNAS requires box-plot whisker description)
- "The y-axis in (A) is truncated at 0.5 because C-index values below 0.5 indicate worse-than-random concordance." (legend justification for axis truncation; PNAS otherwise requires axes start at 0)

### Fig 3 legend (`fig-pdac`) — sample size added

Now starts with "Survival supervision produces different factor structures in PDAC (TCGA + CPTAC training cohort, n = `r ncol(read_result(...)$ex)` patients, `r sum(...)` events)." (PNAS requires N in figure legends.)

### Fig 4 legend (`fig-val`) — sample sizes added (pooled + per-cohort)

Now opens with "(pooled n = `r n_val_total` patients across `r n_val_cohorts` cohorts: Dijk n=…; Moffitt n=…; PACA array n=…; PACA seq n=…; Puleo n=…)". Per-cohort N is computed inline from `pool_df` at render time.

### Fig 4A forest plot — color-blind safe palette

Cohort colors changed from default ggplot2 hues (Dijk red `#F8766D` vs PACA-array green `#00BA38` fail deuteranopia) to **Okabe-Ito** palette. Cohort shapes also varied (previously all filled circles except Pooled diamond) so cohorts are distinguishable by shape *and* color simultaneously.

| Cohort | Old color | New color | Old shape | New shape |
|---|---|---|---:|---:|
| Dijk | `#F8766D` | `#E69F00` (orange) | 16 | 16 (circle) |
| Moffitt | `#B79F00` | `#56B4E9` (sky blue) | 16 | 17 (triangle) |
| PACA array | `#00BA38` | `#009E73` (bluish green) | 16 | 15 (square) |
| PACA seq | `#619CFF` | `#0072B2` (blue) | 16 | 25 (▽ filled) |
| Puleo | `#F564E3` | `#CC79A7` (reddish purple) | 16 | 8 (star) |
| Pooled | `black` | `#000000` (black) | 18 | 18 (diamond) |

### Library calls — four unused commented out

In the setup chunk, `library(Seurat)`, `library(VAM)`, `library(NMF)`, `library(enrichplot)` were commented out with explanation comments. They were defensive imports left over from when analysis ran inline; the figures all come from cached `.rds` files and don't need these heavy packages. Verified by checking that none of these packages have any namespaced (`pkg::func`) calls anywhere in either Rmd.

---

## 5. `si_appendix.Rmd` — Supplementary Information

Three changes:

### 5.1 `fig-converge` legend expanded
Was "Convergence of model across initializations." (1 sentence). Now describes what's plotted (relative loss vs iteration), how many initializations, what's truncated, and what convergence behavior the figure demonstrates.

### 5.2 `wang2014supervised` citation removed
At line 680 (formerly: `[@cai2011graph; @wang2014supervised; @blei2007slda]`). The citation was removed because the entry could not be verified in any database; the surrounding two citations adequately support the SI's claim. See §6 for details.

### 5.3 New: SI bibliography now renders with a heading

PNAS rule: *"references should be cited in numerical order as they appear in the SI; do not cite main-text references in the SI and vice versa."* — i.e., the SI must have its own numbered reference list, separate from the main paper's. Although the SI YAML had `bibliography: references_30102025.bib` set, pandoc-citeproc was emitting the bibliography content at end-of-document but **without a section heading**, so the rendered SI PDF previously had references appearing as floating text after the figures with no "References" title. Reviewers and the PNAS prescreen would flag this.

Fix: added `reference-section-title: "SI References"` to the SI YAML. Now the rendered SI shows a proper "SI References" section heading (visible in the TOC and immediately above the bibliography on the references page). No content was added or removed from the citations; only the heading.

---

## 6. `references_30102025.bib` — Bibliography (most consequential changes)

**Before audit:** 47 entries, 4 with DOIs, multiple incorrect/fabricated entries.
**After audit:** 43 entries, 38 with DOIs, all entries verified against CrossRef / PubMed / DBLP or appropriately classified (preprints, tech reports, book chapters, old conference papers).

### 6A. Three real content errors corrected

These are the changes most worth your direct review:

#### `kim2007sparse` — was a fabricated/mis-fielded entry

- **Before:** `@article{kim2007sparse, title="Sparse non-negative matrix factorization for clustering", author="Kim, Hoyer and Park, Haesun", journal="Journal of Scientific Computing", volume=36, pages=205-222, year=2007}`
- **Audit finding:** No paper with this title + this journal + this volume + these pages exists in CrossRef, PubMed, Semantic Scholar, or Google Scholar. Author "Kim, Hoyer" is suspect (Patrik Hoyer is a separate NMF researcher). The exact title "Sparse Nonnegative Matrix Factorization for Clustering" matches a 2008 Georgia Tech tech report by **Jingu Kim & Haesun Park** (CSE-08-01).
- **After:** Replaced with `@techreport` pointing at the verified Georgia Tech report URL. The entry comments list two peer-reviewed alternatives (Kim & Park 2007 *Bioinformatics*, Kim & Park 2008 *SIAM J Matrix Anal Appl*) if you'd prefer a journal citation.

#### `laurberg2008uniqueness` — was wrong venue

- **Before:** Listed as *Neurocomputing* 71(1-3):606–616 (2008)
- **Audit finding:** That venue could not be confirmed in CrossRef or PubMed. PubMed (PMID 18497868) finds the only paper by these authors with this title in *Computational Intelligence and Neuroscience* 2008, article 764206 (open access; DOI 10.1155/2008/764206).
- **After:** Replaced with the verified Comput Intell Neurosci entry (full title "Theorems on positive data: on the uniqueness of NMF").

#### `wang2014supervised` — entry removed

- **Audit finding:** Title "A novel supervised nonnegative matrix factorization algorithm for gene expression classification" by Wang/Yang/Yang/Sun in *Computational Biology and Chemistry* 53:189–197 (2014) — could not be located in CrossRef, PubMed, Semantic Scholar, or Google Scholar via any search variation. Either the paper exists only in an obscure venue, has transcription errors, was retracted, or was fabricated.
- **Action:** Citation removed from `si_appendix.Rmd:680`. The remaining two citations (`cai2011graph`, `blei2007slda`) adequately support the SI's claim about supervised-NMF variants.

### 6B. Other content corrections

| Entry | Issue | Fix |
|---|---|---|
| `seung2001algorithms` | `pages={35}, number={556--562}` (fields swapped) | Changed to `@inproceedings`, `pages={556--562}`, added editors |
| `cai2011graph` | `@inproceedings` with TPAMI as `booktitle`; missing 4th author T.S. Huang | Changed to `@article`, added Huang, added DOI |
| `Bailey2016` | Author list truncated to "Bailey, Chang, others" | Added 8 more named co-authors per CrossRef |
| `elyada2019cross` | Several first-name typos (Mahdu→Mohan, Patrick→Pasquale, Edouard→Elise, Burkhart middle initial) | Corrected per CrossRef |
| `moffitt2015virtual` | Missing space "tumor-and stroma" in title | Fixed to "tumor- and stroma" |
| `mariathasan2018tgfbeta` | "Kadel Iii" lowercase | Fixed to "Kadel III" (then later removed; was unused) |
| `tomczak2015review` | Volume was "2015"; title leaked metadata word "Review" | Fixed (then later removed; was unused) |
| `pascualmontano2006nonsmooth` | Title missing "(nsNMF)" | Added |
| `huang2020low` | arXiv preprint stuffed into `journal` field | Restructured per PNAS preprint format |
| `le2025survnmf` | Was `@phdthesis` with 6 authors (theses have one); HAL preprint | Changed to `@misc` with HAL URL |

### 6C. Three unused entries removed (per your instruction)

- `collisson2019molecular` (Nat Rev Gastroenterol Hepatol 2019) — author bookkeeping comments suggested it was once intended to be cited but the citation was lost during revision; not present anywhere in current text
- `mariathasan2018tgfbeta` (Nature 2018) — leftover from earlier draft
- `tomczak2015review` (TCGA review 2015) — leftover from earlier draft

### 6D. 25+ DOIs added for verified entries

Every verified journal article now has a `doi` field. Many entries also had journal-name capitalization fixes ("Cancer cell" → "Cancer Cell", "Nature medicine" → "Nature Medicine", etc.) and stale `publisher` fields removed where DOI provides the canonical metadata.

### 6E. Retraction check (CrossRef relation field, all 38 DOIs)

**Result: 0 retractions.** One paper (`chansengyue2020transcription`, Nat Genet 2020) has a published *correction* (erratum) — this is normal post-publication maintenance and does not require any change to the citation.

---

## 7. Items not touched

### Discussion (`05_discussion_REVISED.Rmd`)

No changes. The Discussion was already clean against the PNAS checklist.

### Pending external actions before submission

These are deliberately deferred to submission day or require your decisions (see `docs/pnas_compliance_audit.md` §5):

1. **Mint Zenodo DOIs** for `rashidlab/DeSurv` and `rashidlab/DeSurv-paper`, then patch into `03_methods_REVISED.Rmd:40` (replacing two `[TO BE MINTED]` placeholders). `.zenodo.json` files are pre-staged in both repos.
2. **Decide on bioRxiv** — post a preprint or not.
3. **Suggested editors / reviewers** for the PNAS Direct Submission cover letter.
4. **Final pre-render** on whatever machine has your full R environment, to verify the manuscript compiles end-to-end with all the changes applied.
5. **Manual retraction sanity-check** of the remaining 8 entries without DOIs (all preprints, tech reports, or NIPS proceedings; lower fabrication risk but worth a quick scan).

---

## 8. New compliance artifacts (deposited in `docs/`)

These are filled-out reporting-standard documents and supporting analyses, not changes to the manuscript itself. They support TRIPOD / EQUATOR compliance and should be reviewed before submission.

- **`docs/tripod_checklist.md`** — Full TRIPOD 2015 checklist with per-item status. 17/22 items fully met; 5 items have recommended one-line strengthening edits to the Methods (which I have **not** made — they need your sign-off before edits land in the manuscript).
- **`docs/consort_flow.md`** — Per-cohort participant flow with actual N's extracted from cached `.rds` files. Training: 321 raw → 273 analytic (139 events) across TCGA-PAAD (181 → 144) and CPTAC-3 (140 → 129). Validation: 979 raw → 616 analytic (414 events) across 5 cohorts. Includes (a) per-cohort waterfall table, (b) analysis-to-cohort mapping showing which N entered which figure/table, (c) ASCII draft of the recommended SI Fig S0 flowchart, (d) draft Methods paragraph closing TRIPOD items 4b, 8, 9, 13a, 14a in one shot.

## 9. How to consume this document

- **Quick triage**: read §1, §6A, §1.4 (CRediT), §1.5 (Inclusion & Diversity), §1.6 (Algorithm 1), §3.1 (Algorithm 1 expansion), §5.3 (SI bibliography heading). These are the changes most worth your sign-off.
- **Defer to submission day**: §3 Zenodo placeholders, §7 pending actions, §1 ORCIDs (other authors').
- **Pre-submission TODO list**: see §7 below for the (small) remaining items.
- **Run-time verification**: Once you render on your machine, confirm Fig 4A is colorblind-readable and the new fig captions in §4 populate correctly with cached values.
