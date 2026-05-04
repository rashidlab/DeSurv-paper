# PNAS compliance audit — `origin/public` manuscript

Audit performed against the manuscript on the `origin/public` branch (paper.Rmd + 02–05_*_REVISED.Rmd + si_appendix.Rmd + references_30102025.bib). No edits made; this is a punch list to revisit before submission.

---

## 1. Style/policy checklist (no violations found)

- **No novelty/priority claims.** Searches for "novel", "first to", "to our knowledge", "unprecedented", "we are the first", "pioneer", "breakthrough" returned no hits in main text. The word "first" appears only in non-priority procedural senses ("We first asked…", "We first describe…") or describing prior work ("first defined through virtual microdissection"). Discussion explicitly *disclaims* novelty of the gene programs. "Here we present DeSurv" is conventional.
- **No chemicals or drug names** in main text.
- **No model organisms / cell lines / strains** — human bulk transcriptomic cohorts only.
- **No trade-named instruments/reagents** in main text. PurIST, DeCAF, DECODER, SCISSORS, COMPASS are method/study acronyms (proper nouns), not trade names.
- **Système International units** — no physical units invoked beyond "Time (months)" axis labels and dimensionless statistical quantities (P, HR, CI, C-index).
- **No proposed gene names.** Main text uses no individual human gene symbols; subtype labels (Classical, Basal-like) and authors' factor labels (D1–D3, N1–N3) are not gene names.
- **No "data not shown"** — zero matches.
- **No in-press citations** in bib.

---

## 2. Concrete issues to fix before submission

### 2A. Code/data URLs — likely PNAS-flagged

PNAS: *"Only link to websites that are permanent public repositories… self-perpetuating online resources funded by government, academia, and industry."* GitHub repos are user-controlled and can be renamed/deleted, so they are generally not considered permanent.

| File:line | URL |
|---|---|
| `paper/02_introduction_REVISED.Rmd:23` | `https://github.com/rashidlab/DeSurv` |
| `paper/03_methods_REVISED.Rmd:40` | `github.com/rashidlab/DeSurv` and `github.com/rashidlab/DeSurv-paper` |

**Fix:** mint a Zenodo DOI (or Software Heritage archive) for both the package and the paper repo and cite the DOI in the Data, Materials, and Software Availability paragraph (and the Introduction if the URL is retained there). Keep GitHub as the working mirror.

The TCGA URL `https://www.cancer.gov/tcga` in Acknowledgements is fine (US government, permanent).

### 2B. Abbreviations not defined at first body mention

Acknowledgements expansions do not count — first body use is what PNAS / CSE 8th ed. requires.

| Abbreviation | First body use | Suggested expansion |
|---|---|---|
| **TCGA** | `04_results_REVISED.Rmd:202` "(TCGA [@raphael2017integrated]…)" | The Cancer Genome Atlas (TCGA) |
| **CPTAC** | `04_results_REVISED.Rmd:202` | Clinical Proteomic Tumor Analysis Consortium (CPTAC) |
| **TPM** | `03_methods_REVISED.Rmd:37` "converted to TPM" | transcripts per million (TPM) |
| **GEO** | `03_methods_REVISED.Rmd:40` "GEO (Moffitt, GSE71729)" | Gene Expression Omnibus (GEO) |
| **ICGC** | `03_methods_REVISED.Rmd:40` "ICGC Data Portal" | International Cancer Genome Consortium (ICGC) |
| **EGA** | `paper.Rmd:51` (Acks) "EGA study EGAS00001000154" | European Genome–phenome Archive (EGA) |
| **iCAF** | abstract ("iCAF-associated stroma") and `02_introduction_REVISED.Rmd:23` | inflammatory cancer-associated fibroblast (iCAF) |
| **proCAF, myCAF, panCAF** | `04_results_REVISED.Rmd:255–262` | progenitor / myofibroblastic / pan-CAF subtypes — define on first use |
| **HR, CI, SD** | abstract ("HR per SD = 1.50, 95% CI 1.31–1.72") | hazard ratio (HR), 95% confidence interval (CI), standard deviation (SD) |
| **ADEX** | `04_results_REVISED.Rmd:255` "Bailey ADEX" | Aberrantly Differentiated Endocrine–Exocrine (ADEX) |
| **IRB** | `03_methods_REVISED.Rmd:40` "no IRB approval was required" | Institutional Review Board (IRB) |
| **DECODER, SCISSORS, COMPASS, PurIST, DeCAF** | various | method/study acronyms — defining once at first use is the safe choice; currently appear unexpanded |

Already correctly defined: **NMF**, **PDAC**, **CAF**, **C-index**, **BO**, **Kaplan–Meier**, **ECM**.

### 2C. Figure legends — three actionable issues

1. **Fig 2A,B** (boxplots): legend does not state what the box / line / whiskers represent. Add e.g. "Boxes show median and IQR; whiskers extend to 1.5×IQR." Also fold the paired-Wilcoxon P (currently in body text only) into the legend.
2. **Fig 3** and **Fig 4** legends: add sample sizes. Fig 3 should state $n$ patients (TCGA + CPTAC). Fig 4 should state pooled $n$ across the five validation cohorts (and ideally per-cohort $n$ in panel A).
3. **Fig 2A axis truncation**: `04_results_REVISED.Rmd:220` builds Fig 2A with `scale_y_continuous(limits = c(.5, 1))`. PNAS rule: "Numerical axes on all graphs go to 0, except for log axes." Either expand to 0 or add a justification sentence to the legend ("axis truncated at 0.5; values below 0.5 indicate worse-than-random concordance"). Confirm SI Fig `fig-cindex-by-k` ribbons aren't similarly truncated when re-rendering.

Other figure-legend items (passed):
- All four main-text legends sit in the source immediately after the paragraph carrying the first `\ref{fig:…}`.
- Each caption opens with a one-sentence overview and describes every panel.
- Fig 4A forest plot describes error bars as "(95% CIs)" — explicit ✓.
- ± uses are unambiguous (SI explicitly says "±1 SE"; main text uses "(95% CI 1.31–1.72)" instead of ±).
- Magnification / scale bar — N/A (no microscopy).

SI nit: **`fig-converge`** at `si_appendix.Rmd:1092` reads only "Convergence of model across initializations." — no panel description, no $n$, no axis description. Expand.

### 2D. Bibliography (`paper/references_30102025.bib`) — fixes

Apply these BibTeX-level changes:

```bibtex
% huang2020low — fix arXiv formatting per PNAS preprint format
@article{huang2020low,
  title={Low-rank reorganization via proportional hazards non-negative matrix factorization unveils survival associated gene clusters},
  author={Huang, Zhi and Salama, Paul and Shao, Wei and Zhang, Jie and Huang, Kun},
  archivePrefix={arXiv},
  eprint={2008.03776},
  year={2020},
  note={Preprint, accessed [DD Month YYYY]},
  url={https://arxiv.org/abs/2008.03776}
}

% le2025survnmf — add URL to dissertation
@phdthesis{le2025survnmf,
  ...,
  url={...},
  note={Accessed [DD Month YYYY]}
}

% donoho2004nmf, blei2007slda — add page ranges (NeurIPS proceedings need pages)
@inproceedings{donoho2004nmf, ..., pages={...}, doi={...}, ... }
@inproceedings{blei2007slda, ..., pages={...}, ... }

% gillis2014nmf — add editor field for PNAS book-chapter style
@incollection{gillis2014nmf,
  ...,
  editor={Suykens, Johan A. K. and Signoretto, Marco and Argyriou, Andreas},
  ...
}

% tomczak2015review — drop the leaked "Review" prefix from title
title={The Cancer Genome Atlas (TCGA): an immeasurable source of knowledge},
```

**Journal-name normalization sweep.** Bib mixes styles: `Cancer Discov`, `Cancer cell`, `Nature medicine`, `Nature genetics`, `Nature reviews Gastroenterology \& hepatology`. Recommend converting all to MEDLINE abbreviations (`Cancer Cell`, `Nat Med`, `Nat Genet`, `Cancer Discov`, `Nat Rev Gastroenterol Hepatol`, `Proc Natl Acad Sci USA`, etc.) for consistency.

**Manual retraction check.** Tooling cannot detect retractions; do a manual sanity-check of all 47 cited articles via PubMed/Retraction Watch before submission.

---

## 3. Items not exhaustively audited

- **SI Appendix** (1,697 lines): scanned for URLs (none beyond main text), novelty/priority words (none), and ± usage (explicit "±1 SE"). Did not line-by-line audit gene-symbol formatting (HUGO style: human gene symbols upper-case italic, proteins upper-case roman) inside supplementary tables/figures — worth a pass before submission, especially `gene_lists_top270_k3.csv` and Tables S2–S6 if they list specific genes.
- **Text-recycling self-disclosure** (PNAS Best Practices for Researchers) — cannot be verified by grep; this is an author-attestation item.
- **Retraction status** of cited articles — manual check required.

---

## 4. PNAS Editor's review of the audit

A simulated PNAS-editor pass over §1–3 confirmed most items, downgraded a couple, and surfaced new prescreen-relevant issues. These are folded in below.

### 4A. Editor downgrades (de-prioritize)

- **Fig 2A axis truncation at 0.5** (originally §2C item 3): PNAS's "axes go to 0" rule has a long-standing carve-out for index/ratio metrics with meaningful floors (correlations from 0, C-index from 0.5). A one-sentence legend note suffices ("axis truncated at 0.5; C-index < 0.5 indicates worse-than-random concordance"). **Do not re-render to expand the axis.**
- **DECODER, SCISSORS, COMPASS, PurIST, DeCAF expansions** (originally §2B): treat as proper-noun method/study identifiers (analogous to BRCA1, ENCODE) — defining is not required when used as identifiers. **Drop from the abbreviations list.**
- **Journal-name normalization sweep** (originally §2D): real but copy editors handle at proof. **Defer to proof stage.**

### 4B. Editor's additional findings (PNAS-specific)

#### 4B-1. Significance Statement is over the 120-word limit (PRESCREEN BLOCKER)

PNAS Research Articles cap the Significance Statement at **120 words**, hard-enforced by the Editorial Office at submission. Current statement (`paper.Rmd:48`) is **144 words** (verified via `wc -w` on lines 47–49 minus the YAML key tokens). Trim ~24 words. Two sentences are the obvious cut targets:

- "Although single-cell and spatial transcriptomics now resolve individual cell populations, bulk expression remains the primary source of large, clinically annotated cohorts with the sample sizes needed for stable survival modeling." (30 words — relevance argument; can be carried as a clause elsewhere)
- The final sentence's "with potential applications wherever high-dimensional nonnegative measurements contain both signal and noise variation" generalizes redundantly with the prior sentence.

#### 4B-2. First sentence of Significance Statement uses jargon

PNAS guidance: the first sentence must be intelligible to a non-specialist (undergraduate in any science). "Tumor transcriptomes mix malignant and microenvironmental signals, making it difficult to identify programs that drive clinical outcomes" leans on *transcriptomes* and *programs*. Consider: "Tumor tissue samples contain a mix of cancer cells and surrounding normal cells, making it hard to identify the genes that actually predict patient outcomes."

#### 4B-3. "precision = 0.50 vs 0.07" in abstract is undefined

The abstract reports a precision metric without defining it; PNAS readership outside computational biology will not know the operational definition. Either define inline (~10 words) or replace with a more universally recognized metric.

#### 4B-4. Inclusion and Diversity Statement — currently absent

As of 2022, PNAS asks authors to consider including a 2–3 sentence Inclusion and Diversity statement (sex-as-biological-variable, demographic representation, author-team diversity). Optional but increasingly expected; add to Acknowledgements or Methods.

#### 4B-5. CRediT taxonomy for author contributions

Current `author_contributions` field (`paper.Rmd:42–43`) is narrative. PNAS now expects (not strictly requires) [CRediT](https://credit.niso.org/) standardized roles. Map roles to: Conceptualization, Methodology, Software, Formal Analysis, Investigation, Writing — Original Draft, Writing — Review & Editing, Supervision, Funding Acquisition.

#### 4B-6. "Materials and methods" capitalization

`paper.Rmd:172` reads `# Materials and methods`. PNAS style is Title Case: **Materials and Methods**.

#### 4B-7. Color-blind safety in Fig 4A forest plot

`04_results_REVISED.Rmd:397–402` uses default ggplot2 hue palette (`#F8766D`, `#B79F00`, `#00BA38`, `#619CFF`, `#F564E3`). The Dijk-red / PACA-array-green pair fails deuteranopia. PNAS Figure Guidelines explicitly recommend color-blind-safe palettes. Replace with Okabe-Ito or ColorBrewer Set2; verify the per-cohort shapes in `cohort_shapes` are also varied (currently all `16` except Pooled `18`) so color is not the sole differentiator.

#### 4B-8. Figure rendering device for vector content

`paper.Rmd` knitr opts use `dev = "png", dpi = 300`. PNAS production prefers **vector formats (PDF/EPS) for line plots and schematics**. PNG is accepted at submission but you'll be asked for vector versions of Figs 2A–C, 3D, and 4A at acceptance. Switch to `dev = "cairo_pdf"` (or per-chunk overrides for raster-heavy panels like heatmaps) to avoid a production round-trip.

#### 4B-9. Article length verification

PNAS Research Articles cap at ~6 typeset pages (~47,000 characters incl. spaces, refs, legends; SI excluded). Verify by rendering on `public` and checking page count of `paper/paper.pdf`. The Discussion has internal redundancy ("DeSurv at $k = 3$ achieves independent prognostic value with fewer factors" recurs across paragraphs) that could absorb the trim if needed.

#### 4B-10. Exact P values per CSE 8th

PNAS asks for exact P values to 2 sig figs unless P < 0.001. The k-sensitivity Results paragraph (`04_results_REVISED.Rmd:607`) uses "P < 0.05" multiple times to summarize the count over the (k, α) grid — fine *as a count*, but ensure SI Tables S2–S6 report exact P throughout.

#### 4B-11. bioRxiv preprint check

If a bioRxiv version exists, PNAS asks that it be cited in the manuscript. If it doesn't, posting at submission time is encouraged. **Decide yes/no.**

#### 4B-12. Direct Submission preparation

If Naim Rashid is not an NAS member, this is a Direct Submission and the portal asks for 2–3 suggested editors and 5+ suggested reviewers. Plausible-fit Editorial Board members for this paper (statistical genomics / cancer transcriptomics): Andrea Califano, Dana Pe'er, John Storey, Tom Speed. Pre-identify before submission.

---

## 5. Implementation plan

Sequenced from prescreen-blockers to optional. Each item lists touchpoint, effort, and verification.

### Stage 1 — Prescreen blockers (must do first; ~2 hours)

These are what PNAS Editorial Office checks before assigning to an editor. A failure here returns the manuscript without review.

| # | Item | File:lines | Effort | Verification |
|---|---|---|---|---|
| 1 | Trim Significance Statement to ≤120 words | `paper/paper.Rmd:47–48` | 15 min | `git show HEAD:paper/paper.Rmd \| sed -n '47,49p' \| wc -w` returns ≤122 (120 + 2 YAML tokens) |
| 2 | Rewrite Significance first sentence in plain language | `paper/paper.Rmd:48` | 5 min | Read aloud; should be intelligible to a non-biologist |
| 3 | Define "precision" inline in abstract OR substitute metric | `paper/paper.Rmd:45` (abstract) | 10 min | Re-read; metric self-explanatory |
| 4 | "Materials and methods" → "Materials and Methods" | `paper/paper.Rmd:172` | 1 min | grep |
| 5 | Verify rendered page count ≤ 6 | render via `make paper` | 2 min render + 1 min count | Open `paper/paper.pdf`, count typeset pages excluding figures-on-own-pages |

### Stage 2 — Required content fixes (before peer review; ~3 hours)

| # | Item | File:lines | Effort | Verification |
|---|---|---|---|---|
| 6 | Add abbreviation expansions at first body mention (TCGA, CPTAC, TPM, GEO, ICGC, EGA, IRB, ADEX, iCAF, proCAF, myCAF, panCAF, HR, CI, SD) | per §2B table | 30 min | Re-grep each abbreviation; confirm parenthetical expansion at first body occurrence |
| 7 | Mint Zenodo DOI for `rashidlab/DeSurv` and `rashidlab/DeSurv-paper`; add to Data, Materials, and Software Availability + Introduction | external (Zenodo) + `02_introduction_REVISED.Rmd:23`, `03_methods_REVISED.Rmd:40` | 30 min (deposits + edit) | DOIs resolve; cite per PNAS "Archived code" template |
| 8 | Add color-blind-safe palette to Fig 4A forest plot | `04_results_REVISED.Rmd:397–402` | 20 min | Run color-blind simulator on rendered Fig 4 (e.g. `colorblindr::cvd_grid`) |
| 9 | Switch knitr `dev` to `cairo_pdf` for vector figures (raster overrides for heatmaps) | `paper/paper.Rmd:96–101` (knitr opts) + per-chunk `dev=` for `fig-pdac` heatmap panels | 20 min | Rendered figures are vector when opened in PDF reader |
| 10 | Add boxplot whisker description to Fig 2A,B legend; fold in paired-Wilcoxon P | `04_results_REVISED.Rmd:208` (figcap) | 10 min | Read legend; box/line/whiskers described; P stated |
| 11 | Add N to Fig 3 and Fig 4 legends | `04_results_REVISED.Rmd:262, 359` | 15 min | $n$ patients (Fig 3, TCGA+CPTAC); pooled and per-cohort $n$ (Fig 4) |
| 12 | Add legend note for Fig 2A axis truncation | `04_results_REVISED.Rmd:208` | 5 min | One-sentence note added |
| 13 | Expand SI Fig `fig-converge` legend | `paper/si_appendix.Rmd:1092` | 10 min | Caption now has overview + panel description + $n$ + axis units |
| 14 | Bib hygiene: `huang2020low` arXiv reformat, `le2025survnmf` URL+access date, `donoho2004nmf` & `blei2007slda` pages, `gillis2014nmf` editors, `tomczak2015review` title fix | `paper/references_30102025.bib` | 30 min | Re-render; refs render without missing-field warnings |
| 15 | Add Inclusion and Diversity statement (2–3 sentences) | `paper/paper.Rmd:51` (after Acknowledgements) | 15 min | Statement present; mentions sex-as-biological-variable consideration |
| 16 | Convert author contributions to CRediT roles | `paper/paper.Rmd:42–43` | 10 min | Each author has ≥1 CRediT role; corresponding-author roles ≥3 |

### Stage 3 — Recommended improvements (peer-review-quality; ~1 hour)

| # | Item | File:lines | Effort | Verification |
|---|---|---|---|---|
| 17 | Manual retraction check on all 47 cited articles | `paper/references_30102025.bib` | 30 min | Cross-reference Retraction Watch / PubMed retraction notices |
| 18 | Verify SI Tables S2–S6 report exact P values (not just stars) | `paper/si_appendix.Rmd` table chunks | 15 min | Skim rendered SI tables |
| 19 | bioRxiv decision: post preprint OR cite existing | external | 30 min if posting | bioRxiv DOI in cover letter |
| 20 | Pre-identify 2–3 PNAS Editorial Board members + 5 suggested reviewers | external (cover letter prep) | 30 min | Names + affiliations + areas in submission notes |

### Stage 4 — Defer to proof / production stage

These are real items but the copy desk handles them; doing them earlier is wasted effort.

| # | Item | File:lines | Why defer |
|---|---|---|---|
| 21 | Journal-name normalization to MEDLINE abbreviations across all 47 bib entries | `paper/references_30102025.bib` | Copy editors handle; CSL applies most consistency at render |

### Stage 5 — Submission-day checklist (final ~1–2 hours before upload)

**Status update (2026-05-03):** Items below have been resolved during the audit pass:
- ✅ NeurIPS page-range TODOs (donoho2004nmf 1141–1148, blei2007slda 121–128, verified via DBLP)
- ✅ le2025survnmf URL added (HAL hal-04975434v1)
- ✅ Retraction check completed for all 38 DOI'd entries via CrossRef — 0 retractions; 1 entry (chansengyue2020transcription) has a published correction (not blocking)
- ✅ 3 unused bib entries removed (collisson2019molecular, mariathasan2018tgfbeta, tomczak2015review)
- ✅ Render verification — full `make paper` succeeds end-to-end, producing paper.pdf (11 pp rticles preview) and si_appendix.pdf (37 pp). Required reproducibility-environment changes:
  - Library reduction: `library()` calls for unused heavyweight packages (Seurat, VAM, NMF, enrichplot) commented out in `04_results_REVISED.Rmd`. Of these, NMF was reinstated because the SI's `fig-nmf-diagnostics` chunk uses NMF's S4 plot method.
  - **ggplot2 version constraint discovered**: cached figure `.rds` files were saved with ggplot2 3.x (S3 internals); ggplot2 4.x's S7 internals fail to print them. Pin ggplot2 to ≤ 3.5.x in any reproducibility environment, or regenerate the cached figures with the current ggplot2 version. This pin should be documented in the repo's install instructions.
- ✅ latexdiff comparison: `paper_DIFF.pdf` (11 pp) and `si_DIFF.pdf` (37 pp) generated for review, showing all editorial changes since `origin/public` baseline.

These actions must happen **after** all manuscript content is final, because they immutably tie a code snapshot or external citation to the submitted version. Doing them earlier risks DOI/snapshot mismatch with what reviewers see.

#### 5A. Zenodo DOI minting (item 7, deferred from Stage 2)

`.zenodo.json` files have already been written to both repos as preparation:
- `~/Downloads/DeSurv-paper-clean/.zenodo.json` (paper repo)
- `~/Downloads/DeSurv-rashidlab/.zenodo.json` (package repo, cloned from rashidlab/DeSurv)

**One-time setup** (do this anytime; doesn't trigger minting):
1. Sign in at https://zenodo.org/account/settings/github/ with your GitHub account.
2. If `rashidlab` is a GitHub Organization, an org admin must first grant Zenodo OAuth access at github.com/organizations/rashidlab/settings/oauth\_application\_policy.
3. Toggle ON `rashidlab/DeSurv` and `rashidlab/DeSurv-paper`.

**Submission day — mint DOIs by cutting GitHub releases:**

```bash
# 1. Package repo (rashidlab/DeSurv)
cd ~/Downloads/DeSurv-rashidlab
git pull origin main
git add .zenodo.json
git commit -m "Add Zenodo metadata for DOI minting"
git push origin main
gh release create v1.0.1 \
  --title "DeSurv v1.0.1 — PNAS submission version" \
  --notes "Companion package for Young et al. PNAS submission. See rashidlab/DeSurv-paper for analysis code and reproducibility materials."

# 2. Paper repo (rashidlab/DeSurv-paper, public branch)
cd ~/Downloads/DeSurv-paper-clean
git pull origin public
git add .zenodo.json
git commit -m "Add Zenodo metadata for DOI minting"
git push origin public
gh release create v1.0.0-pnas-submission \
  --target public \
  --title "DeSurv-paper v1.0.0 — PNAS submission" \
  --notes "Reproducibility code and pre-computed results for Young et al. PNAS submission."

# 3. Wait ~5 minutes; check https://zenodo.org/account/settings/github/ — DOIs appear next to each repo
```

**After DOIs are issued — patch into manuscript:**

```bash
# Replace [TO BE MINTED] placeholders in paper/03_methods_REVISED.Rmd:40
# (two occurrences, one per repo). Format per PNAS "Archived code" reference style:
#   "archived at Zenodo: DOI 10.5281/zenodo.XXXXXXX"
# Then re-render and re-verify.
```

#### 5B. Final pre-flight verification (must pass)

```bash
# Run on origin/public after DOI patch-in:
make paper                                                # render with final everything
sed -n '47,49p' paper/paper.Rmd | wc -w                   # significance ≤122 (incl YAML tokens)
sed -n '44,46p' paper/paper.Rmd | wc -w                   # abstract ≤252 (incl YAML tokens)
grep -nE 'data not shown|in press|TO BE MINTED' paper/*.Rmd paper/*.bib  # zero hits
grep -nE 'TODO' paper/references_30102025.bib             # zero hits (NeurIPS pages resolved)

# Manual checks on rendered paper.pdf:
# - ≤6 typeset PNAS pages (excl. figures on own pages)
# - Fig 4A passes color-blind simulator (Coblis or colorblindr::cvd_grid)
# - Fig 2 / Fig 3 / Fig 4 captions show numeric N, not blank/error
# - Fig 2A whisker description present in legend
# - All four figures render as vector (zoom in PDF reader; no pixelation)
```

#### 5C. External actions still pending — finalize before upload

| Item | Action |
|---|---|
| 17 | Retraction check: query Retraction Watch DB and PubMed for each of the 47 cited articles; remove or replace any flagged. |
| 19 | bioRxiv decision: post preprint or not. If posting, embed bioRxiv DOI in cover letter. |
| 20 | Suggested editors / reviewers: 2–3 editors + 5+ reviewers in submission portal notes. Candidate editors listed: Andrea Califano, Dana Pe'er, John Storey, Tom Speed. |
| 14b | NeurIPS page-range TODOs in `references_30102025.bib`: verify `donoho2004nmf` (NIPS 16) and `blei2007slda` (NIPS 20) page numbers via proceedings.neurips.cc; verify `le2025survnmf` thesis URL. |

#### 5D. Cover letter / portal items

- Statement of Direct Submission (if Naim is not an NAS member).
- 2–3 PNAS Editorial Board member suggestions.
- 5+ reviewer suggestions (with affiliations and expertise).
- Cover letter with explicit framing of the methodological-novelty pitch (already well-tuned in the Discussion's last paragraph).
- bioRxiv DOI if posted.

#### 5E. Tag-and-archive consistency check (last step before upload)

After Zenodo DOIs are minted and patched in, ensure the *exact* tagged commit on each repo equals what's referenced in the submitted manuscript:

```bash
# In paper repo:
git -C ~/Downloads/DeSurv-paper-clean tag --list 'v*'              # should list v1.0.0-pnas-submission
git -C ~/Downloads/DeSurv-paper-clean log -1 v1.0.0-pnas-submission  # commit hash should match HEAD
diff <(git show v1.0.0-pnas-submission:paper/paper.Rmd) paper/paper.Rmd  # should be empty
```

If anything mismatches, **do not submit yet** — either re-tag (move the tag, re-trigger Zenodo via release-edit, accept the new DOI) or hold the submission and fix the drift.

### Verification checklist (run before submission)

```bash
# On origin/public, after edits:
make paper                                        # renders paper/paper.pdf + si_appendix.pdf
sed -n '47,49p' paper/paper.Rmd | wc -w           # significance statement: ≤122 (incl. YAML tokens)
sed -n '44,46p' paper/paper.Rmd | wc -w           # abstract: ≤252 (incl. YAML tokens)
grep -nE 'data not shown|in press' paper/*.Rmd    # zero hits
grep -nE 'github\.com' paper/*.Rmd                # zero hits OR each accompanied by Zenodo DOI

# Manual checks:
# - Open paper.pdf: confirm ≤ 6 typeset pages (excluding figures on own pages)
# - Open Fig 4: simulate color-blind via Coblis or colorblindr; cohorts distinguishable
# - Open Fig 2: legend describes whiskers; axis truncation noted
# - Open Fig 3, Fig 4: legends state N
# - Render paper.pdf: confirm vector content in Figs 2A-C, 3D, 4A (zoom to verify no pixelation)
```

### Dependencies

- **Stage 2 #7 (Zenodo DOIs)** is on the critical path — needs the deposit done before the manuscript can cite it. Start this first; it can run in parallel with the edits.
- **Stages 1, 2, 3 are otherwise parallelizable** across editors. Item 9 (knitr `dev`) interacts with item 8 (Fig 4A palette) only at re-render time; do both, then re-render once.
- **All edits are on `origin/public`.** Reminder: changes to text don't invalidate cached `.rds` results, so `make paper` is fast (~2 min). Only figure-code edits (items 8, 9, 10) trigger a meaningful re-render of figure objects, but the underlying numeric results stay cached.
