# Pending edits: clinical-significance ("so what") closers

_Drafted 2026-05-04. Status: **not yet applied to manuscript.** Review and edit before committing._

The current closers of the abstract, significance statement, and introduction all end on method/scope framing. The goal of this set of edits is to add an explicit "who benefits, why does this matter" sentence at each of those three exits — grounded in what the paper actually demonstrates (prognostic transportability across five cohorts), not what it doesn't (treatment selection, prospective utility).

---

## Post-review summary (added 2026-05-04)

External review of this document on 2026-05-04 evaluated each of §§1–8 below. **The recommendations and final drafts in this summary supersede the original drafts in §§1–8 wherever they conflict.** The original drafts and rationale are retained below for documentation; consult the post-review final draft in this section when applying.

**Scope note:** this status table covers §§1–8 only. For §§9–29, see the cluster organization at the top of this file (Cluster A–J) and the per-section coordination notes embedded in each section. Section-level streamlining updates from the 2026-05-05 paragraph-level preview pass are folded into individual sections (see "Application priority across clusters" below for sequence).

### Status table

| Section | Status | Key change vs. original draft |
|---|---|---|
| §1 — Abstract closer | Apply with revision | Sharpened to address variance-prognosis misalignment specifically rather than restating "transportable signatures." |
| §2a — Significance closer | Apply with revision | Soften "therapeutic development" phrase to align with §3 (no overreach into target discovery). |
| §2b — Significance opening trim | **Reverse** | Original "Tumor tissue samples contain a mix..." opening is more accessible for non-specialists. Cut "and using automatic model selection" elsewhere instead. |
| §3 — Introduction closer | Apply with revision | Strip "struggled to generalize" clause (contradicts intro's PurIST/Basal-like positioning). Soften "target-discovery hypotheses" to "biological insight into tumor-stroma coupling" (defensible — paper does identify this coupling). |
| §4 — scRNA-seq preemption | Apply compact version + Discussion retrim (revised 2026-05-05; long version retained as fallback) | Verify Discussion paragraph exists via grep before editing (principle is sound; specific text needs confirming). |
| §5a — k-selection literature in Results | Apply | No change. |
| §5b — k-selection clause in Intro P2 | **Skip** | Length grounds: undoes recent 1087→790 trim. §5a alone delivers the defensive payoff. |
| §6 — BO citation + why-clause | Apply | No change. |
| §7 — top-weighted clarification | Apply after SI verification | Verification is **mandatory, not optional**, with both wording-alignment and bug-detection value. |
| §8 — precision/C-index gap | Apply Option B + bridge | Option A is redundant with Option B's mechanism description. |

### Post-review final drafts

#### §1 — Abstract closer (final, post-2026-05-05 streamlining)

Replaces the original draft addition. To append after the existing closer "...standard discover-then-evaluate paradigm.":

> **DeSurv yields reproducible prognostic gene programs for the bulk-transcriptomic cohorts that dominate clinical cancer research, addressing the variance-prognosis misalignment that limits unsupervised factorization.**

**Why this revision (post-paragraph-preview, 2026-05-05):** Two streamlining edits over the prior post-review draft:

1. **Dropped the lead clause** "By prioritizing transcriptional programs that survive cross-cohort generalization, " — paragraph-level preview revealed this paraphrases the immediately preceding sentence's "transportable transcriptional signatures than the standard discover-then-evaluate paradigm." Same beat twice in 30 words.
2. **"biomarkers" → "gene programs"** — "biomarkers" carries CLIA-grade analytically-validated-assay connotations; the paper validates prognostic gene programs, not clinical assays. "Gene programs" stays inside what the paper formally demonstrates.

Net delta vs original abstract: ~+22 words (was ~+32). Total abstract: ~230 words (within 250-word PNAS limit, with margin).

**Earlier post-review draft (now superseded):** *"By prioritizing transcriptional programs that survive cross-cohort generalization, DeSurv yields reproducible prognostic biomarkers for the bulk-transcriptomic cohorts that dominate clinical cancer research, addressing the variance-prognosis misalignment that limits unsupervised factorization."*

#### §2 — Significance statement (final)

**Reverse §2b — keep the original opening sentence as-is:**

> Tumor tissue samples contain a mix of cancer cells and surrounding normal cells, making it hard to identify the genes that actually predict patient outcomes.

**Trim "and using automatic model selection" from the existing fourth sentence** to recover the words instead:

| Original | Final |
|---|---|
| "By embedding outcome information into discovery **and using automatic model selection**, DeSurv yields clinically relevant, reproducible programs across cohorts." | "By embedding outcome information into discovery, DeSurv yields clinically relevant, reproducible programs across cohorts." |

**Apply §2a closing-sentence replacement, with the "therapeutic development" phrase softened to align with §3:**

| Original | Final |
|---|---|
| "Reproducible programs of this kind are a prerequisite for patient risk stratification, clinical-trial enrichment, and **the prioritization of biological hypotheses for therapeutic development** from the bulk-transcriptomic cohorts that dominate clinical cancer research." | "Reproducible programs of this kind are a prerequisite for patient risk stratification, clinical-trial enrichment, and **biological insight into tumor-stroma coupling and other prognostic axes** from the bulk-transcriptomic cohorts that dominate clinical cancer research." |

**Why the revision:** "biological hypotheses for therapeutic development" parallels §3's "target-discovery hypotheses" — both drift past what the paper demonstrates. The paper *does* identify a Classical-tumor + iCAF coupling and presents it as a substantive observation; "biological insight into tumor-stroma coupling" stays inside that scope.

#### §3 — Introduction closer (final)

Replaces both the original draft and the review's suggested replacement (which stripped "target-discovery" without preserving any contribution claim). To replace the existing closing sentence at `paper/02_introduction_REVISED.Rmd:23`:

| Original draft (rejected) | Post-review final |
|---|---|
| "...offering a route to reproducible prognostic biomarkers **and target-discovery hypotheses** for cancers — like PDAC — **where existing expression-based signatures have struggled to generalize across cohorts**." | "...offering a route to reproducible prognostic biomarkers **and biological insight into tumor-stroma coupling** in cancers like PDAC **where compositional variation dominates expression**." |

**Why the revision (combining review's catch with my pushback):**
- "Struggled to generalize" contradicts the intro's own positioning of PurIST and Basal-like/Classical as cohort-confirmed signatures (review's catch — correct).
- "Target-discovery hypotheses" drifts past what the paper tests (review's catch — partially correct), but the suggested blanket strip loses defensible contribution scope. The paper identifies a specific Classical-tumor + iCAF coupling and presents it as a biological observation, so "biological insight into tumor-stroma coupling" preserves the claim within demonstrated bounds.
- "Where compositional variation dominates expression" cites a mechanism the paper actually invokes (`aran2015systematic` already in intro for pan-cancer purity-variance evidence; the §3 closer now points at the same mechanism rather than an unspecified generalization gap).

#### §4 — Pre-edit verification step (added)

Before applying the §4 long-version intro insertion, **verify the Discussion paragraph exists in the form assumed by the review**. Run:

```bash
grep -nE "single.cell|spatial transcriptomics|cohort sizes|stable survival modeling" paper/05_discussion_REVISED.Rmd
```

Expected hit: a paragraph beginning approximately *"Although single-cell and spatial transcriptomics now resolve tumor compartments directly..."*. If found, plan a coordinated retrim of that paragraph concurrently with the intro addition to avoid two near-identical defenses appearing 6 pages apart. If the Discussion text differs from what was assumed, adjust the coordination plan; do not apply §4 in isolation.

#### §7 — Verification step reframed (mandatory, with bug-detection value)

The original §7 description called the SI verification step "optional" / "before committing." Replace with the following framing:

The SI Section 6 verification is **mandatory, not optional**, for two distinct reasons:

1. **Wording alignment.** The inline gloss in 7a was inferred from `R/get_top_genes.R`, not from SI Section 6. The Results inline gloss must match SI Section 6 verbatim where possible.
2. **Bug detection (independent of whether 7a is applied).** If SI Section 6 describes a different procedure than `R/get_top_genes.R` actually implements (e.g., SI says "top genes by absolute loading" but code does column-max-normalized differential ranking), that is a manuscript-vs-code consistency bug that must be resolved in one place or the other — independent of whether the Results inline gloss is added. The verification step is therefore an audit checkpoint for code-vs-manuscript consistency, not just a wording check.

Outcome of verification dictates next action:
- **SI matches code:** apply 7a/7b/7c with wording aligned to SI.
- **SI differs from code:** flag the discrepancy as a separate bug, decide whether to fix the code or fix the SI, then proceed with 7a/7b/7c.

#### §8 — Final preference (Option B + bridge)

Apply **Option B** (the ~+30-word version) plus the bridge sentence. Option A is redundant with Option B's mechanism description. The bridge sentence is recommended (not optional) because it is purely factual and makes the simulation→external-validation relationship explicit without any new claim.

### Items the review missed or under-weighted

For documentation, recording the issues raised in my critique of the review that did not appear in the review itself but are addressed in the final drafts above:

1. **§2 / §3 parallel "therapeutic development" phrasing.** Review caught §3's "target-discovery hypotheses" but not §2's "prioritization of biological hypotheses for therapeutic development" — same overreach, two locations. Both softened to "biological insight" in the final drafts above.
2. **§1 sharpening could create coherence drift with §4.** Review's suggested sharpening ("where unsupervised factorization remains the field default") read as more antagonistic to scRNA-seq-aware bulk methods than is consistent with §4's complementarity framing. Final draft uses "addressing the variance-prognosis misalignment that limits unsupervised factorization" instead — names the specific limitation rather than dismissing the alternative wholesale.
3. **§7 verification has bug-detection value beyond wording alignment.** Review treated verification as wording-only; the final §7 reframing makes the audit-checkpoint use case explicit.
4. **§4 Discussion-redundancy needs grep-verification before editing.** Review asserted the Discussion text without verifying; the final §4 verification step makes this explicit.

---

## Word-budget situation (PNAS limits)

| Section | Current | Limit | Headroom |
|---|---:|---:|---:|
| Abstract | 207 | 250 | 43 |
| Significance | 118 | 120 | 2 |
| Introduction | uncapped | — | — |

Significance statement is at ceiling — adding without trimming will push over. Plan below trims the first sentence to compensate.

---

## Cluster organization (added 2026-05-05)

The 28 numbered sections (§1–§29; §16 was never filed) are grouped below into 10 thematic clusters. Section numbers are preserved for traceability with the post-review summary's status table and existing cross-references throughout the document.

| Cluster | Theme | Sections |
|---|---|---|
| A | Clinical-translation framing ("so what" closers) | §1, §2, §3 |
| B | Matched-rank labeling audit | §9, §14, §15, §18 |
| C | NMF higher-rank transparency + principled-rank-selection contribution claim | §10, §19, §20 |
| D | Type III / joint-vs-individual / "weak not no" framing calibration | §8, §12 |
| E | Methods accuracy and citations | §5, §6, §7, §25, §26, §27 |
| F | Discussion structural cleanup | §17, §21, §22, §24 |
| G | Convergent biological validation | §28 |
| H | scRNA-seq positioning + clinical deployment | §4, §23 |
| I | Mechanical clarity / wording fixes | §11, §13 |
| J | SI consistency | §29 |

Within each cluster, sections appear in their original numerical order. Cross-cluster dependencies (e.g., §20a depends on §10; §29a affects §10's draft) are documented inline within each section's coordination notes.

## Application priority across clusters

The cluster organization above is **thematic**, not priority-driven. The recommended **application sequence** is roughly inverse to the thematic order at the top — factual corrections first, rhetorical strengthening last. A 5-wave structure that maps onto the clusters:

| Wave | Theme | Sections |
|---|---|---|
| 1 | **Factual corrections** (verifiable errors a careful reader will catch) | §7 (top-weighted definition + SI cross-reference, Cluster E); §10 (NMF k=5/k=9 acknowledgment, Cluster C); §27 (TPM provenance, Cluster E); §29 (SI Section 20 internal consistency, Cluster J) |
| 2 | **Matched-rank audit** (single coordinated grep + edit pass) | All of Cluster B (§9, §14, §15, §18) |
| 3 | **Free wins / mechanical** | §11 (Cluster I); §22 (Cluster F); §23 (Cluster H); §25, §26 (Cluster E) |
| 4 | **Contribution-claim arc** (coordinated set; apply together) | §10 + §17 + §19 + §20 + §21 + §24 (across Clusters C, F); plus §28 (Cluster G) and §29 (Cluster J) for the convergent-validation thread |
| 5 | **Closing rhetorical edits** | All of Cluster A (§1, §2, §3); §4 (Cluster H); §5 (Cluster E); §13 (Cluster I) |

**Hard dependencies that constrain ordering:**

- **§20a depends on §10** (apply both in same wave; §20a's conditional framing requires §10's k=5 acknowledgment).
- **§29a precedes §10 application** (§29a fixes SI Section 20 wording so the SI is consistent with main text at apply-time; §10's actual draft cites Section 16 + Fig. S10, not Section 20, so no §10 *draft* revision is needed — only the application sequencing constraint).
- **§4 application requires Discussion line 14 retrim** (specified in §4's checklist; §22+§23 already touch line 14 in Wave 3).
- **§17 + §19 + §20 + §21 form a coordinated Discussion rewrite** in Wave 4; §21's redundancy fixes assume §19/§20 have moved content out.

The waves can deploy independently if word-budget or co-author bandwidth forces sequencing — Waves 1–3 (factual + matched-rank + free wins) are stable as a partial deployment without committing to Waves 4–5 (contribution arc + closers). The reverse is **not** safe: Waves 4–5 without Waves 1–3 leave reviewer trip-hazards in place under stronger rhetorical claims.

---

# Cluster A — Clinical-translation framing ("so what" closers)

The original three "so what" closers (abstract, significance, intro), plus shared notes/dials and an application checklist that originally accompanied this set. **§24 (Discussion close)** is the de facto fourth member of this arc but is filed in Cluster F because of its primary Discussion-cleanup nature.

---

## 1. Abstract

**Current closer (paper.Rmd line 60):**

> ...demonstrating that outcome supervision during factorization yields more transportable transcriptional signatures than the standard discover-then-evaluate paradigm.

**Proposed (append one sentence):**

> ...demonstrating that outcome supervision during factorization yields more transportable transcriptional signatures than the standard discover-then-evaluate paradigm. **By prioritizing transcriptional programs that survive cross-cohort generalization, DeSurv supports the development of reproducible prognostic biomarkers and risk-stratification tools for the bulk-transcriptomic cohorts that dominate clinical cancer research.**

Net: ~233 words, within the 250-word cap.

---

## 2. Significance statement

**Two-part edit** (replace closing sentence + trim opening sentence to stay under 120 words).

### 2a. Replace closing sentence

**Current (paper.Rmd line 63):**

> This provides a general framework for identifying prognostically informative gene programs without requiring post hoc filtering or pre-specified signatures, applicable across cancer types where bulk expression dominates clinical cohort data.

**Proposed:**

> This provides a general framework for identifying prognostically informative gene programs without requiring post hoc filtering or pre-specified signatures. **Reproducible programs of this kind are a prerequisite for patient risk stratification, clinical-trial enrichment, and the prioritization of biological hypotheses for therapeutic development from the bulk-transcriptomic cohorts that dominate clinical cancer research.**

### 2b. Tighten opening sentence (compensates for added length)

**Current:**

> Tumor tissue samples contain a mix of cancer cells and surrounding normal cells, making it hard to identify the genes that actually predict patient outcomes.

**Proposed:**

> Bulk tumor samples mix cancer cells with surrounding normal cells, obscuring the genes that truly predict patient outcomes.

Net: ~122 words. Still ~2 over the 120 cap — may need one more small trim elsewhere (e.g., "and using automatic model selection" → drop, since selection is implementation detail not significance) to fully clear the limit. Recheck word count after editing.

---

## 3. Introduction

**Current closer (02_introduction_REVISED.Rmd line 23, second-to-last sentence):**

> DeSurv provides a general framework for extracting outcome-relevant programs from bulk tumor transcriptomes wherever unsupervised factors conflate biological signal with compositional noise. The method is implemented in an open-source R package (https://github.com/rashidlab/DeSurv).

**Proposed (extend the framework sentence; leave the package sentence as-is):**

> DeSurv provides a general framework for extracting outcome-relevant programs from bulk tumor transcriptomes wherever unsupervised factors conflate biological signal with compositional noise**, offering a route to reproducible prognostic biomarkers and target-discovery hypotheses for cancers — like PDAC — where existing expression-based signatures have struggled to generalize across cohorts**. The method is implemented in an open-source R package (https://github.com/rashidlab/DeSurv).

---

## Notes / dials to consider before committing

- **Audience mix.** Drafts above lean toward translational researchers + clinical biomarker development. If you'd rather emphasize patients more directly (e.g., "patients with PDAC, where five-year survival remains < 12%"), the abstract/intro additions can absorb a single concrete framing clause without breaking the word budget.
- **Overclaim watch.** All three drafts deliberately stop short of "guides treatment decisions" or "enables precision therapy" — the paper validates *prognostic* transportability (HR per SD across five cohorts), not treatment-selection utility. If a co-author or reviewer pushes for stronger language, the demonstrated claim is "supports patient risk stratification" / "trial enrichment," not direct clinical use.
- **Parallel phrasing.** Drafts share the cadence "reproducible … programs/biomarkers + bulk-transcriptomic cohorts that dominate clinical cancer research" across all three exits. Repetition is intentional (theme reinforcement across abstract → significance → intro), but if it reads heavy you can vary the framing in one of the three.
- **PDAC vs. general.** The intro draft names PDAC explicitly ("cancers — like PDAC — where..."); the abstract and significance keep it general. This mirrors the paper's actual scope: PDAC is the demonstration, generality is the claim. If you'd rather keep all three general (or all three PDAC-anchored), adjust the intro accordingly.

---

## When ready to apply

Files to touch:
1. `paper/paper.Rmd` — abstract block (line 60) and significance block (line 63).
2. `paper/02_introduction_REVISED.Rmd` — final paragraph (line 23) for the so-what closer; **paragraph 1 (line 15) for the scRNA-seq preemption (§4 below)**.
3. After editing, recount significance-statement words; trim further if > 120.
4. Re-render with `make paper` and visually check the front-matter renders cleanly in `paper.pdf`.

---

# Cluster B — Matched-rank labeling audit

Four instances of the same audit rule ("when comparing DeSurv to NMF, label the $k$ used for NMF") at four manuscript sites. Apply as a coordinated set.

**Cross-cluster note:** §15 + §18 establish the manuscript-wide audit grep that should also catch any additional unspecified-$k$ instances. **Earlier candidate site dropped after paragraph-level preview (2026-05-05):** Results line 261 was initially flagged as a Cluster B audit candidate ("light coverage") for the matched-rank qualifier. Paragraph-level preview confirmed this is redundant after §9 application — §9 establishes "we fit standard NMF at the same rank ($k = 3$)" at line 253, so all subsequent rank references in the same paragraph zone are unambiguously $k = 3$. **No edit needed at line 261.**

---

## 9. Soften the matched-rank framing and disambiguate "rank" in the PDAC factor-comparison setup

### Why this matters

The matched-rank comparison sentence at `paper/04_results_REVISED.Rmd:253` has two distinct issues that are easier to fix together than separately:

1. **Asymmetric framing reads as handicapping.** The phrase *"we fit standard NMF at DeSurv's BO-selected rank ($k = 3$)"* leads with a possessive attribution to DeSurv. A skeptical reader can land on this and interpret the comparison as "we constrained NMF to use our preferred $k$, ensuring we'd win," before the justification arrives in the next sentence. This is a rhetorical asymmetry, not a methodological flaw — NMF's own rank-selection diagnostics did not converge on a clear alternative (Fig. S4 already establishes this), so matching ranks is the principled apples-to-apples comparison. The fix is to remove the possessive and add a callback to Fig. S4 so the reader sees the justification at the same time as the choice.

2. **"Rank" is overloaded across the paragraph and at one bridging phrase the conflation risk is sharp.** Within the surrounding text:
   - "the same rank ($k = 3$)" — factorization rank, disambiguated by `$k$`
   - "factor-specific gene **rankings**" — gene ordering by loading
   - "Matching **ranks**" — factorization rank, **no disambiguator**, sits adjacent to "gene rankings"
   - "**rank**-selection diagnostics" — compound noun, contextually clear
   - "across a range of **ranks**" — sensitivity context, contextually clear

   The bare "Matching ranks" is the most ambiguous use because it lands in the same line as "gene rankings" with no math symbol or compound-noun anchor. A reader scanning the paragraph can chain them. The fix is to replace the bare phrase with one that uses the math symbol `$k$` directly, breaking the lexical chain without scrubbing the standard NMF terminology elsewhere in the paragraph.

The two fixes coordinate cleanly into one revised sentence pair, so they should be applied together.

### Proposed edit — consolidated final draft

**Current (`paper/04_results_REVISED.Rmd:253`):**

> To compare the factor structures produced with and without survival supervision at the same dimensionality, we fit standard NMF at **DeSurv's BO-selected rank** ($k = 3$) and examined the overlap between each method's factor-specific gene rankings and established PDAC gene programs (Fig. \ref{fig:pdac}A--B). **Matching ranks ensures** that differences in factor content reflect the presence or absence of survival supervision rather than differences in model complexity; whether NMF's performance improves at higher ranks is evaluated in the sensitivity analysis below and in SI Appendix, Table S5.

**Proposed:**

> To compare the factor structures produced with and without survival supervision at the same dimensionality, we fit standard NMF at **the same rank** ($k = 3$) and examined the overlap between each method's factor-specific gene rankings and established PDAC gene programs (Fig. \ref{fig:pdac}A--B). **Matching $k$ across methods isolates** the effect of survival supervision on factor content from differences in model complexity; **standard NMF's own rank-selection diagnostics did not converge on a clear alternative (Fig. S4),** and the sensitivity analysis below (and SI Appendix, Table S5) evaluates NMF performance across a range of ranks and supervision strengths.

### Three coordinated changes

1. **"DeSurv's BO-selected rank" → "the same rank"** (drops possessive attribution; both methods are simply described as being fit at $k = 3$).
2. **"Matching ranks ensures" → "Matching $k$ across methods isolates"** (math symbol breaks lexical chain with "gene rankings"; "across methods" reinforces between-method framing).
3. **Add Fig. S4 callback** (*"standard NMF's own rank-selection diagnostics did not converge on a clear alternative (Fig. S4)"*) — pre-empts the "but what's NMF's natural $k$?" question by reminding readers Fig. S4 already showed there isn't one. No new analysis needed.

The verb swap from "ensures" to "isolates" is incidental but slightly stronger — "isolates the effect" is the standard methodological phrasing for a controlled comparison.

### Why these specific choices

**On the handicapping framing:**
- Don't over-soften. Apologetic framings ("acknowledging that fitting NMF at this rank is conservative...") draw *more* attention to the issue than they resolve. The minimal fix — drop the possessive, add the Fig S4 callback — neutralizes the asymmetry without flagging it as a vulnerability.
- The Fig S4 callback is load-bearing. Without it, "the same rank" is just a passive description; with it, the reader sees the methodological justification (NMF *can't* prefer a different rank because its own diagnostics didn't converge) at the same moment as the choice.
- The sensitivity analysis later in the paragraph empirically backs this up by showing NMF doesn't rescue itself at higher ranks (line 609 region: "$k = 3$ achieved significance ... across [N] of [N] supervision strengths, and $k = 7$ at a comparable [N] of [N]"). The framing softening doesn't need to mention this — the data does the work — but it's worth knowing that the sentence revision is consistent with the empirical backstop already in the paper.

**On the rank/ranking disambiguation:**
- Don't try to scrub every instance of "rank" from the paragraph. "Rank-selection diagnostics" and "across a range of ranks" are standard NMF terminology that readers expect. The targeted fix is the one bare bridging phrase.
- "Matching $k$ across methods" preserves the parallel structure with the original "Matching ranks" rhetorically, while using the math symbol to disambiguate. "Holding $k$ fixed" was considered but loses the explicit between-method framing.

### Cost/benefit

- **Word delta:** ~+8 words (the Fig S4 callback). The two within-sentence fixes (possessive drop + math-symbol swap) are net-neutral.
- **No bib changes.** No new analyses required.
- **Reader benefit:** removes the asymmetric handicapping framing that would catch a skeptical reviewer's eye, and breaks the lexical chain that risks conflating factorization rank with gene rankings.
- **Risk:** very low. The methodological content is unchanged; the matched-rank choice is the same; the figure references are the same.

### Application checklist

1. Apply the consolidated final draft above to `paper/04_results_REVISED.Rmd:253`, replacing the existing sentence pair.
2. Verify the Fig. S4 reference is correct in the final rendered manuscript (the citation key for the unsupervised-NMF-rank-diagnostics figure — likely `\ref{fig:nmf-diagnostics}` per `paper/si_appendix.Rmd:1241–1243`; double-check exact reference syntax for SI figures called from main text).
3. Re-render with `make paper` and confirm the paragraph still reads cleanly with the revised framing.

---

## 14. Symmetric framing of the validation transition sentence; specify $k$ for each method

### Why this matters

The transition sentence to the external-validation Results at `paper/04_results_REVISED.Rmd:353` reads:

> "A prognostic factorization is useful only if it generalizes beyond the training cohort. Having established that survival supervision and standard NMF produce different factor structures in the training data, we next tested whether **DeSurv's** structure transfers to independent cohorts without retraining."

Two issues stack:

1. **Asymmetric framing** ("whether DeSurv's structure transfers") implies only DeSurv is being tested. A reader who scans this sentence and then skips the result paragraph will conclude "DeSurv was tested for transfer; NMF apparently was not" — exactly the handicapping interpretation. The actual analysis at line 355 *is* symmetric: it explicitly tests Standard NMF at three ranks ($k = 3$ matched, $k = 5$ elbow-selected, $k = 7$ BO-selected) and reports that NMF at $k = 7$ "matched or exceeded DeSurv's per-cohort C-index in most cohorts." The framing-vs-analysis mismatch is the issue.

2. **$k$ asymmetry between methods is unstated** in the framing. DeSurv is tested at one rank ($k = 3$); NMF is tested at three. The framing sentence as currently written gives no signal that NMF is being given multiple ranks, so a reader builds expectations that don't match what line 355 then describes. This is more noticeable in the post-§14-fix version that says "each method's structure" — symmetric in language but asymmetric in configuration unless disambiguated.

The fix is one sentence revision that combines both: change the possessive to symmetric framing AND add a parenthetical specifying the per-method configurations.

### Proposed edit (combined fix)

**Current (`paper/04_results_REVISED.Rmd:353`):**

> Having established that survival supervision and standard NMF produce different factor structures in the training data, we next tested whether **DeSurv's** structure transfers to independent cohorts without retraining.

**Proposed (recommended):**

> Having established that survival supervision and standard NMF produce different factor structures in the training data, we next tested whether **each method's** structure transfers to independent cohorts without retraining **(DeSurv at $k = 3$; NMF at multiple ranks, detailed below)**.

Two coordinated changes:
1. **"DeSurv's" → "each method's"** — symmetric framing. The downstream paragraph already tests both, so this aligns the setup with the analysis.
2. **Parenthetical specifying configurations** — transparent about asymmetric $k$ (DeSurv = 1 rank, NMF = multiple ranks). Pre-empts reader surprise when line 355 introduces $k = 5$ and $k = 7$ without warning.

Net word delta: ~+10 words.

### Alternative drafts (if the parenthetical is too heavy)

**Option B — Symmetric framing only, no $k$ specification:**

> Having established that survival supervision and standard NMF produce different factor structures in the training data, we next tested whether **each method's** structure transfers to independent cohorts without retraining.

Two-word swap. Resolves the handicapping framing but leaves the $k$ asymmetry unspecified — the reader still has to wait for line 355 to learn NMF is tested at multiple ranks. Use this if the parenthetical reads as too detail-heavy for a transition sentence.

**Option C — Most explicit:**

> Having established that survival supervision and standard NMF produce different factor structures in the training data, we next tested whether **DeSurv's $k = 3$ gene programs and standard NMF's gene programs (at the matched rank and at NMF's own selected ranks)** transfer to independent cohorts without retraining.

Names everything explicitly. Heavier (~+18 words) but unambiguous. Use only if a reviewer or co-author has flagged the lack of clarity.

### Coordination with §10 (line 355 vs. line 609)

Worth noting that **line 355 already partially implements §10 Path A** for the validation Results: it explicitly acknowledges *"Standard NMF at elbow-selected $k = 5$ and BO-selected $k = 7$ substantially improved over NMF at $k = 3$, with $k = 7$ matching or exceeding DeSurv's per-cohort C-index in most cohorts"* and then reframes DeSurv's value as parsimony + adjusted significance. This is the same move §10 Path A recommends for the sensitivity-analysis paragraph at line 609.

Implication for §10 application: the line 355 framing is the model to follow when applying §10 Path A at line 609. Coordinating the two paragraphs so they handle "NMF at higher $k$ does well" with consistent framing (acknowledge → reframe to parsimony/adjustment) prevents the paper from sounding like it engages with the issue in one place and ducks it in another.

This is metadata for the §10 application step, not a new §14 edit. But it's worth knowing while applying §10 that line 355 has already done much of the work and §10 Path A's drafts can be lightly adapted to that paragraph's existing tone.

### Cost/benefit

- **Word delta:** ~+10 words (recommended draft), or +2 (Option B), or +18 (Option C).
- **No new analyses, no bib changes.**
- **Reader benefit:** removes the framing-vs-analysis asymmetry that creates the handicapping interpretation, AND pre-empts surprise about NMF being tested at multiple ranks.
- **Risk:** very low. The downstream analysis content is unchanged.

### Application checklist

1. Choose the recommended combined draft (default), Option B (minimal), or Option C (maximal).
2. Apply the chosen edit to `paper/04_results_REVISED.Rmd:353`, replacing the existing transition sentence.
3. If applying §10 Path A at line 609, coordinate with the line 355 framing so the "NMF at higher $k$ does well → DeSurv reframed as parsimony/adjustment" move lands consistently across both paragraphs.
4. Re-render with `make paper` and confirm the validation paragraph still reads cleanly with the revision.

---

## 15. Specify matched-rank $k = 3$ at unspecified-$k$ HR sentences (line 355 pooled HR; line 357 dichotomized HR); manuscript-wide $k$-audit recommended

### Update (2026-05-05): paragraph-level preview suggests one of 15a/15b is redundant after §14

Paragraph-level preview of the validation paragraph (post-§14 + §15a + §15b) revealed verbosity: §14's transition sentence already specifies "(DeSurv at $k = 3$; NMF at $k = 3$, 5, and 7)" upfront, then §15a inserts "at the matched rank ($k = 3$)" at the pooled HR comparison, then §15b inserts it again at the dichotomized HR. **Three matched-rank specifications in two paragraphs is too much.**

Streamlining recommendation (revised 2026-05-05): apply **§14 + §15a only**. **Drop §15b.** Reasoning: the pooled-HR sentence (§15a's site) is where the matched-rank ambiguity *first arises* — it sits immediately after the paragraph discusses "NMF at $k = 5$ and $k = 7$" and is where the reader's expectation about $k$ is being formed. Specifying matched-rank at §15a sets the context that the dichotomized comparison (§15b's site, in the next paragraph) inherits. The dichotomized sentence can rely on the matched-rank context once §15a has established it; specifying it twice is redundant.

Earlier recommendation (now superseded) was §14 + §15b. The reversal reflects that §15a is the more dangerous ambiguity site because it directly follows the higher-rank NMF discussion. Either option remains technically defensible; the revised recommendation is the stronger of the two. Do **not apply both 15a and 15b**.

The original §15a and §15b drafts below stand as wording references for whichever of the two is applied.

### Why this matters

This is the fourth occurrence of the unspecified-$k$ pattern (after §9 / line 253; §14 / line 353; and the two new sites here). At each, a comparison sentence reports a result for "the analogous NMF" or "the same procedure applied to NMF" without specifying which $k$ NMF was fit at. From the variable names in the supporting code (`val_latent_std_desurvk_tcgacptac` at line 142 and `fig_median_survival_std_desurvk_tcgacptac` at line 179, both with `std_desurvk` = standard NMF at DeSurv's $k = 3$), the answer is consistently $k = 3$ matched. The choice is methodologically right (apples-to-apples matched-rank comparison), but the implicit $k$ creates two problems:

1. **Reader confusion** — by line 357 the reader has just been walked through "NMF at $k = 7$ matched or exceeded DeSurv's per-cohort C-index" (line 355 first half). When the next sentences report HR results for "NMF" without specifying $k$, the natural inference is "NMF at $k = 7$ — the version that just won on per-cohort C-index" — which is wrong.
2. **Handicapping interpretation** — a skeptical reader can read the unspecified $k$ as cherry-picking NMF's weakest configuration ($k = 3$, where it is matched but underperforms) for the dichotomized comparison while having shown $k = 7$ doing well for per-cohort C-index. The comparison structure is correct; the framing is not transparent about it.

The line 357 case is the most exposed of the four because the priming from the immediately preceding sentences makes the unstated $k$ assumption least likely to be guessed correctly.

### The two specific edits

#### 15a. Line 355 — pooled HR sentence

**Current:**

> ...the DeSurv linear predictor showed a consistent survival association (stratified Cox model: pooled HR per SD 1.50; 95% CI 1.31--1.72; P < 0.001 ...). **The analogous NMF linear predictor** showed a weaker pooled effect (HR per SD 1.11; 95% CI 1--1.23; P = 0.057); a direct comparison using dichotomized risk groups is presented below.

**Proposed (preferred):**

> ...the DeSurv linear predictor showed a consistent survival association.... **At the matched rank ($k = 3$), the analogous NMF linear predictor** showed a weaker pooled effect (HR per SD 1.11; 95% CI 1--1.23; P = 0.057); a direct comparison using dichotomized risk groups is presented below.

Net delta: ~+5 words.

#### 15b. Line 357 — dichotomized HR sentence

**Current:**

> The DeSurv linear predictor separated survival trajectories in the pooled validation cohort (HR = ?, 95% CI ?, P < 0.001; Fig. 4B). **The same procedure applied to NMF** yielded weaker stratification (HR = 1.48, 95% CI 1.04--2.11, P = 0.028; Fig. 4C).

**Proposed:**

> The DeSurv linear predictor separated survival trajectories in the pooled validation cohort.... **The same procedure applied to NMF at the matched rank ($k = 3$)** yielded weaker stratification (HR = 1.48, 95% CI 1.04--2.11, P = 0.028; Fig. 4C).

Net delta: ~+6 words.

### Optional broader transparency move (consider with §10 coordination)

Reporting NMF's pooled HR / dichotomized HR **at $k = 7$** in addition to $k = 3$ would be the most fully transparent move — analogous to line 355's per-cohort C-index discussion which already reports NMF at multiple ranks. This would address the §10 concern at lines 355 and 357 in the same spirit as line 355's first-half acknowledgment of NMF at $k = 5/7$ matching per-cohort. Sub-decision:

- **If applying §10 Path A elsewhere** (e.g., at line 609 sensitivity analysis): consider reporting NMF $k = 7$ pooled HR and dichotomized HR here too, to keep the "transparent across all comparisons" framing consistent.
- **If applying §10 Path B** (minimal acknowledgment): the matched-$k = 3$ specification (15a, 15b above) is sufficient; the line 355 first-half NMF-at-multiple-$k$ discussion already does the §10 Path A move for per-cohort C-index.

This sub-decision affects whether the validation paragraph reports two HRs per metric (matched $k = 3$ + NMF's preferred $k = 7$) or just one (matched $k = 3$). Either is defensible; the dual-reporting version requires looking up or computing the NMF $k = 7$ pooled HR + dichotomized HR if not already in the cached results.

### Pattern observation: manuscript-wide $k$-audit recommended

This is the fourth occurrence of the same issue. Three or four piecemeal fixes risk catching some sites and missing others; a coordinated audit pass is small effort and catches the whole set:

```bash
grep -nE "NMF|standard NMF" paper/04_results_REVISED.Rmd | grep -iE "HR|hazard|cox|c-index|c.index|stratif|kaplan|KM|cutpoint|linear predictor|forest|pooled|across.*cohorts"
```

Plus a similar pass on `paper/05_discussion_REVISED.Rmd` and `paper/si_appendix.Rmd`. For each hit, check whether the implied $k$ for NMF is specified or implicit. The audit rule:

> **When comparing DeSurv to NMF, label the $k$ used for NMF.** If NMF is at the matched rank, write "at the matched rank ($k = 3$)" or "(matched rank)". If NMF is at a different rank, write "(NMF at $k = N$)". Never leave the NMF $k$ implicit when DeSurv's is also implicit.

This is the same rule §9 and §14 implicitly applied; the manuscript-wide audit just makes it consistent.

### Coordination with §12 framing calibration

While re-reading §12 in light of the user's question about HR 1.48 (apparent tension with "NMF has no prognostic association"), one phrase in §12 is too strong and should be softened on a re-read before applying:

- §12 currently uses calibrated phrases like "barely prognostic," "marginal over chance," "weak model-level prognostic content" — these are accurate.
- But §12 also contains the harder phrasing "does not have prognostic content at this rank, jointly or individually" — this is too absolute. NMF $k = 3$ has **weak** prognostic content, not **no** prognostic content. External-validation HR per SD = 1.11 (P = 0.057, borderline) and dichotomized HR = 1.48 (P = 0.028) demonstrate that NMF $k = 3$ is weakly prognostic in absolute terms, just much weaker than DeSurv at the same rank.
- **Recommended §12 framing calibration**: review the §12 drafts and replace any "does not have prognostic content" phrasing with "has weak prognostic content" or similar. The empirical claim (DeSurv $k = 3$ cv-C-index 0.65 vs NMF $k = 3$ cv-C-index 0.55, plus larger external-validation HRs across all metrics) is fully consistent with this calibrated framing and does not require the absolute claim.

The substantive §12 argument (per-factor $\Delta\ell$ contrast reflects overall model-level differences, not multicollinearity) is unchanged. Only the absolute-vs-weak language needs a small calibration.

### Cost/benefit (combined)

- **Word delta:** ~+11 words across both 15a and 15b (or ~+25–35 if also adding NMF $k = 7$ HR per the optional transparency move).
- **No new analyses** for 15a and 15b. NMF $k = 7$ pooled/dichotomized HR may require a quick look at cached results if applying the transparency move.
- **No bib changes.**
- **Reader benefit:** removes ambiguity at the two most-exposed instances of the unspecified-$k$ pattern; pre-empts handicapping interpretation; calibrates §12 framing so external validation results don't contradict it.
- **Risk:** very low. Methodological content unchanged; only labeling and framing.

### Application checklist

1. Apply 15a to `paper/04_results_REVISED.Rmd:355` (pooled HR sentence).
2. Apply 15b to `paper/04_results_REVISED.Rmd:357` (dichotomized HR sentence).
3. Run the manuscript-wide grep audit (above). For any additional unspecified-$k$ instances, apply the same audit rule.
4. Re-read §12 in this parking-lot file and soften any "does not have prognostic content" phrasing to "has weak prognostic content" (or equivalent).
5. Decide on the optional dual-HR transparency move (15c) in coordination with §10 Path A vs. Path B decision.
6. Re-render with `make paper` and confirm the validation paragraphs still read cleanly.

---

## 18. Extend §15's matched-rank audit to the Discussion file (line 8 demonstrated catch)

### Why this matters

§15 flagged a manuscript-wide pattern of unspecified-$k$ comparisons in `paper/04_results_REVISED.Rmd` (lines 253, 353, 355, 357 — the matched-rank choice was correct but the $k$ was implicit, leaving comparisons readable as handicapping). The recommended audit was a grep across the Results file. **The same pattern appears in the Discussion file**, and at least one specific instance is a clear catch:

`paper/05_discussion_REVISED.Rmd:8` (within the opening Discussion paragraph):

> "Standard NMF, by contrast, allocated a dedicated factor to exocrine-compositional variation, a signal that reflects tumor purity rather than cancer biology [@rashid2020purity], **leaving less capacity for the microenvironmental programs that DeSurv separates**."

The "less capacity" claim is **true at matched $k = 3$ and not at higher $k$**. At $k = 3$, NMF spends one of three factors (N2) on exocrine variation, leaving only two factors for tumor + microenvironmental programs. At $k = 5$ or $k = 7$, NMF has slack to separate exocrine *and* keep tumor + microenvironmental programs distinct — and the per-cohort C-index in line 355 of Results confirms this empirically (NMF at $k = 7$ matches DeSurv). Without a matched-rank qualifier, a reader who has already read line 355 can object: "NMF doesn't have *less capacity* in absolute terms — it had less at the matched rank you chose."

### Three drafts for the immediate fix

**Option A — Append "at the same rank" (lightest, ~+4 words):**

> "Standard NMF, by contrast, allocated a dedicated factor to exocrine-compositional variation, a signal that reflects tumor purity rather than cancer biology [@rashid2020purity], leaving less capacity for the microenvironmental programs that DeSurv separates **at the same rank**."

**Option B — Insert "at the matched rank ($k = 3$)" upfront (~+5 words):**

> "Standard NMF **at the matched rank ($k = 3$)**, by contrast, allocated a dedicated factor to exocrine-compositional variation..."

**Option C — Append "at the matched rank ($k = 3$)" (most explicit, consistent with §9/§14/§15 pattern, ~+6 words, preferred):**

> "Standard NMF, by contrast, allocated a dedicated factor to exocrine-compositional variation, a signal that reflects tumor purity rather than cancer biology [@rashid2020purity], leaving less capacity for the microenvironmental programs that DeSurv separates **at the matched rank ($k = 3$)**."

### Recommendation: Option C

Two reasons:

1. **Uniform phrasing across the manuscript.** §9 (line 253), §14 (line 353), §15 (lines 355, 357) all use "at the matched rank ($k = 3$)" or close variants. Using the same phrasing in the Discussion helps readers build a mental model that "matched rank" = $k = 3$ consistently. Option A is acceptable but breaks the phrasing pattern.
2. **The Discussion is sometimes read independently.** A reader who jumps to Discussion without reading the full Results may not have "the matched rank" in working memory. Specifying $k = 3$ explicitly removes ambiguity.

### Extend the §15 audit to include the Discussion file

The §15 manuscript-wide audit checklist currently grep-references `04_results_REVISED.Rmd`. **Update it to also include `05_discussion_REVISED.Rmd` and `si_appendix.Rmd`.** Specifically:

```bash
# §15-style audit, extended:
grep -nE "NMF|standard NMF" paper/04_results_REVISED.Rmd paper/05_discussion_REVISED.Rmd paper/si_appendix.Rmd | grep -iE "HR|hazard|cox|c-index|c.index|stratif|kaplan|KM|cutpoint|linear predictor|forest|pooled|across.*cohorts|capacity|fewer|more.*factor|concentrate|separate|distribute|allocate"
```

The added keyword filters (`capacity`, `fewer`, `more.*factor`, `concentrate`, `separate`, `distribute`, `allocate`) catch the kind of claim made in the line 8 Discussion sentence — "less capacity," "DeSurv separates," "allocated a dedicated factor," "distribute across multiple factors," etc. These are claims that are true at matched rank but not necessarily at higher $k$, and need the qualifier.

The audit rule extends:

> **When comparing DeSurv to NMF — whether on HR, C-index, dichotomized stratification, factor capacity, or any structural property — label the $k$ used for NMF.** If NMF is at the matched rank, write "at the matched rank ($k = 3$)" or "(matched rank)". If NMF is at a different rank, write "(NMF at $k = N$)". Never leave the NMF $k$ implicit when DeSurv's is also implicit.

### Other likely catches in the Discussion (worth verifying)

Without doing the full grep, other candidates in the Discussion paragraph at line 8 to check:

- "**unsupervised methods distribute across multiple factors**" — same kind of claim. True at matched $k = 3$; less obviously true at higher $k$. May want "at the same rank" appended.
- "**The factor with the largest survival contribution (D1)**" — implicit DeSurv $k = 3$. Probably OK in Discussion since DeSurv's $k = 3$ is the model the paper defends, but worth verifying.
- "**This co-occurrence within a single factor**" — implicit matched rank. Probably OK because "within a single factor" is descriptive of the result, not a comparison claim.

A grep pass would catch any others.

### Cost/benefit (for §18 specifically)

- **Word delta:** ~+6 words for the line 8 Discussion fix. Plus possibly ~+10–15 words across other Discussion instances if the audit catches more.
- **No bib changes, no new analyses.**
- **Reader benefit:** removes the "less capacity in absolute terms" misreading; aligns Discussion framing with the matched-rank conventions established in Results.
- **Risk:** very low. Methodological content unchanged.

### Application checklist

1. Apply Option C (or Option A/B per preference) to `paper/05_discussion_REVISED.Rmd:8` — append "at the matched rank ($k = 3$)" to the "DeSurv separates" claim.
2. Run the extended §15 audit grep across `04_results_REVISED.Rmd`, `05_discussion_REVISED.Rmd`, and `si_appendix.Rmd` to catch any additional unspecified-$k$ instances.
3. For each catch, apply the audit rule (label NMF $k$ — matched if implicit, explicit otherwise).
4. Re-render with `make paper` and confirm the Discussion paragraph reads cleanly with the revisions.

---

# Cluster C — NMF higher-rank transparency + principled-rank-selection contribution claim

Distributed argument: §10 acknowledges NMF at $k = 5$/$k = 9$ in the sensitivity analysis with the SI-aligned "convergent evidence" framing; §19 surfaces the plateau argument and "NMF needs BO" finding earlier (Results lines 206 and 355) and foreshadows principled rank selection in abstract/intro; §20 sharpens to the explicit two-part contribution framing (gradient + model selection framework).

**Hard dependency:** §20a depends on §10 — apply both or neither (otherwise §20a's "matches only when DeSurv's tuning framework applied" framing is exposed without the §10 acknowledgment that justifies it).

**Coordination with Cluster J:** §29 fixes SI Section 20's framing to align with §10's main-text reframing.

---

## 10. Address standard NMF performance at $k = 5$ and $k = 9$ in the sensitivity analysis

**Status:** Multiple Path A drafts and a new Path C captured below. No final recommendation — to be decided later. **The original Path A draft below contains factual errors flagged in subsequent review and should NOT be applied without revision.**

### Update (2026-05-05, later): SI verification reveals Section 16 already provides the right framing

**Critical SI finding (verified 2026-05-05).** The SI at `paper/si_appendix.Rmd` line 1353 already contains the exact rhetorical framing the main-text §10 should adopt — verbatim:

> "$k = 5$ at $\alpha=0.55$ achieves the highest factor-program fidelity in the $n_{\text{top}} = 270$ analysis (H-cor = ... vs. $k = 3$'s ...) and also achieves adjusted significance ($P = ...$), **providing convergent evidence for the iCAF biology from a higher-$k$ model**; an independent $k = 5$ all-genes fit at $\alpha=0.25$ likewise achieves adjusted significance..."

This is **substantially stronger than concession-and-reframe**. The SI doesn't say "k=5 also generalizes (oh well)" — it says "k=5 generalizes *and that's convergent evidence for the same biology DeSurv recovers at k=3*." All of the prior Path A drafts (original, review-revised, extended) read as defensive after-the-fact justifications when held against this SI framing.

**SI-aligned Path A draft (supersedes original Path A, review-revised Path A, and extended Path A) — softened 2026-05-05:**

> Standard NMF at $k = 5$ also achieves external significance (PurIST/DeCAF-adjusted $P = 0.0005$), **providing convergent evidence for the iCAF biology that DeSurv concentrates at $k = 3$ (SI Appendix, Section 16)**. The supervised $k = 3$ solution achieves the same biology more parsimoniously: the additional factors at $k = 5$ capture transcriptional variance without new factor-level prognostic contribution (SI Appendix, Fig. S10). DeSurv's contribution is therefore **not unconditional dominance over every NMF rank, but a supervised selection framework that identifies a compact, interpretable, externally validated solution without relying on post hoc inspection of multiple unsupervised ranks**.

**Why this softening (2026-05-05):** The earlier wording said *"DeSurv's principled rank selection avoids the diagnostic-disagreement problem that leaves NMF without a reliable path to $k = 5$ (Fig. S4)"* — but Results line 355 already mentions "elbow-selected $k = 5$" as one of NMF's diagnostic outputs. A skeptical reviewer could object that NMF's elbow does provide a path to $k = 5$, even if it's not consensus across all diagnostics. The softer reframing — "not unconditional dominance over every NMF rank, but a supervised selection framework that identifies a compact, interpretable, externally validated solution without relying on post hoc inspection" — is harder to attack while making the same contribution claim.

**Earlier wording (now superseded):** *"...DeSurv's principled rank selection avoids the diagnostic-disagreement problem that leaves NMF without a reliable path to $k = 5$ (Fig. S4)."*

Word delta vs original draft: ~+50 words. Strengths over the prior Path A drafts:

- **Frames k=5 as corroboration not concession** — directly mirrors SI Section 16's "convergent evidence" language. Main text and SI now tell the same story.
- **Cites SI Section 16 as source of the convergent-evidence framing** — reader can verify that the SI substantively backs the main-text claim, not the reverse.
- **Cites Fig. S10 for the "no new factor-level contribution" claim** — grounds the parsimony argument in a specific figure rather than asserting it.
- **Cites Fig. S4 for the diagnostic-disagreement claim** — same as the prior Path A but more concise.
- **Drops the "fails at p=0.56, 0.03, 0.51" framing entirely** — no more overstating individual ranks.

The "no prognostically relevant structure" language elsewhere in the SI (Section 20 region) is consistent with this framing because it refers specifically to **factor-level Type III log-likelihood contribution at k=5** (Fig S10 shows N4, N5 contribute negligibly), not to overall model performance. The SI distinguishes "no new independent factor-specific survival contribution" (per-factor $\Delta\ell$) from "the full higher-rank NMF linear predictor still validates" (external dichotomized stratification). Both are true simultaneously — different metrics. The SI-aligned main-text draft above respects this distinction by citing Fig S10 for the per-factor claim and SI Section 16 for the convergent-evidence claim.

**This SI-aligned draft replaces the recommended Path A in the prior update subsection below.** Path C (count-metric disambiguation) remains a valid lighter complementary fix.

### Update (2026-05-05): revised drafts and a complementary lighter option

Two subsequent observations have refined the §10 approach beyond the original Path A / Path B framing:

1. **The original Path A draft contained a factual error** flagged by external review of this document. The original draft claimed the cophenetic, silhouette, and residual criteria favor $k \in \{2, 3, 7\}$ and that the unsupervised baseline "fails" at those ranks ($p = 0.56$, $0.03$, $0.51$ adjusted). Two problems: (a) the elbow criterion does converge on $k = 5$ (per line 355 of the Results: *"elbow-selected $k = 5$"*), so saying "diagnostics don't identify NMF-success ranks" is inaccurate — the right framing is "no consensus across diagnostics"; (b) calling $p = 0.03$ "failing" is overstatement — that's borderline-significant, not failing. The **review-revised Path A draft** below corrects both.

2. **The count-metric ambiguity is a separate but related issue.** The sensitivity analysis paragraph reports "$k = 3$ retained significance at 5 of 7 supervision strengths, while $k = 7$ dropped to 2 of 7" — but the denominator of 7 includes one $\alpha = 0$ (NMF) and six $\alpha > 0$ (DeSurv configurations). The count conflates NMF with DeSurv, and a reader can't tell whether NMF specifically is in the significant set at any given $k$. From the cached `adj_p_270_matrix.rds`: at $k = 3$, NMF ($\alpha = 0$) IS in the 5 significant ($p = 0.03$, weakly); at $k = 7$, NMF is NOT in the 2 significant ($p = 0.51$). Disambiguating the count metric without engaging the broader contribution reframe is a lighter complementary fix (**Path C** below).

#### Review-revised Path A (corrects factual errors in original draft)

> Standard NMF at $k = 5$ also achieved external significance (PurIST/DeCAF-adjusted $p = 0.0005$). However, no single unsupervised rank-selection criterion consistently identifies $k = 5$ as preferred (Fig. S4): a researcher using only standard NMF would not reliably arrive at this rank, and the unsupervised baseline at the ranks its own diagnostics do favor yields markedly weaker external validation. The contribution of survival-supervised model selection is therefore principled, generalizable rank choice — combined with the interpretability of the resulting Classical-tumor + iCAF coupling — rather than unconditional improvement over standard NMF at every rank.

Word delta: ~+70 words. Two improvements over the original Path A:
- "No single criterion consistently identifies $k = 5$ as preferred" — accurate even though elbow does favor $k = 5$, because the *other* diagnostics don't agree (no consensus is the right framing).
- "Markedly weaker external validation" — accurate as a description of $p \in \{0.56, 0.03, 0.51\}$ as a set, without overclaiming "failing."
- Adds the interpretability beat (Classical + iCAF coupling) that was missing in the original.

#### Extended Path A (full NMF transparency at every tested $k$)

If you want to address the user's specific question ("what about NMF at $k = 3$?") and pull NMF out of the count-metric aggregate at every rank, the extended draft is:

> Standard NMF retained adjusted external significance at $k = 3$ (weakly, $p = 0.03$), $k = 5$ ($p = 0.0005$), and $k = 9$ ($p = 0.0003$), but failed at $k = 2$ and $k = 7$ ($p = 0.56$ and $p = 0.51$). NMF's own rank-selection diagnostics (Fig. S4) point to different ranks without consensus — the cophenetic, silhouette, residual, and elbow criteria each favor a different value — so a researcher using only standard NMF would not reliably converge on a rank that generalizes across cohorts. The contribution of survival-supervised model selection is therefore principled, generalizable rank choice — combined with the interpretability of the resulting Classical-tumor + iCAF coupling — rather than unconditional improvement over standard NMF at every rank.

Word delta: ~+95 words. Trade-off vs. review-revised Path A:
- **More transparent**: every NMF result is on the page, including the borderline $k = 3$ retention. A reader can see exactly what NMF does at each tested rank.
- **More words**: ~+25 over review-revised Path A.
- **Removes the count-metric ambiguity**: by listing NMF's results per rank, the reader doesn't need to reverse-engineer which $\alpha$ values contribute to "5 of 7" vs "2 of 7."
- **Flag**: the empirical claim that diagnostics "each favor a different value" needs verification against Fig. S4 before applying. From the manuscript text, we know elbow → $k = 5$; the assignments for cophenetic, silhouette, and residual need to come from Fig. S4 directly.

#### Path C — Count-metric disambiguation (lighter complementary fix)

If you want lighter transparency without committing to the full Path A reframe, you can disambiguate the count metric directly. Two variants:

**Path C light (~+15 words):**

> "$k = 3$ retained adjusted significance across 5 of 7 supervision strengths (including standard NMF, $\alpha = 0$, weakly), while $k = 7$ dropped to 2 of 7 (neither at $\alpha = 0$)."

**Path C medium (~+35 words, more disclosure):**

> "$k = 3$ achieved significance ($P < 0.05$) across 6 of 7 supervision strengths (including standard NMF at $\alpha = 0$), and $k = 7$ at 5 of 7 (excluding standard NMF, which failed). After PurIST/DeCAF adjustment, $k = 3$ retained significance at 5 of 7 (including NMF, weakly, $p = 0.03$), while $k = 7$ dropped to 2 of 7 (still excluding NMF). Standard NMF retained adjusted significance specifically at $k = 3$, $k = 5$, and $k = 9$ (SI Appendix, Table S3)."

Path C addresses the count-metric ambiguity but does NOT articulate the principled-rank-selection reframing of the contribution claim. Use Path C alone if you want lighter transparency without engaging the contribution-claim reframe; use it together with Path A (review-revised or extended) for both.

#### Decision matrix: which option(s) to apply

| Combination | Transparency about NMF | Contribution reframe | Word delta | Comments |
|---|:---:|:---:|---:|---|
| Original Path A | Partial ($k = 5$ only) | Yes | +65 | **Has factual errors — DO NOT APPLY** |
| Review-revised Path A | $k = 5$ only | Yes | +70 | Defensible if you only want to pull out the strongest counter-example |
| Extended Path A | Every tested $k$ | Yes | +95 | Most rigorous; requires Fig. S4 verification of diagnostic assignments |
| Path B (count reframe to DeSurv-only) | None | Implicit | +10–20 | Leaves NMF's results undefended; reads as evasive |
| Path C light | Partial (which $\alpha$ in the count) | None | +15 | Honest about counts but doesn't reframe contribution |
| Path C medium | Full (all $k$) | None | +35 | Most disclosure without contribution reframe |
| **Review-revised Path A + Path C light** | $k = 5$ + count disambiguation | Yes | +85 | Reasonable middle ground |
| **Extended Path A alone** | Every tested $k$ | Yes | +95 | Subsumes Path C; cleanest at the cost of more words |
| **Extended Path A + Path C medium** | Heavy redundancy | Yes | +130 | Probably over-disclosure |

**Recommended pairings:**

- If you want **maximum rigor and full transparency**: Extended Path A alone. It subsumes the count-disambiguation and articulates the contribution reframe in one paragraph. The trade-off is ~+95 words in an already-long Results paragraph.
- If you want **lighter transparency without committing to contribution reframe**: Path C light or medium. Honest about the counts; doesn't engage the §10 strategic question. Adequate if §15-style framing fixes are applied at lines 355 and 357 elsewhere.
- If you want the **middle ground**: Review-revised Path A + Path C light. Names $k = 5$ explicitly, reframes contribution, AND disambiguates the count metric. ~+85 words but covers both aspects.
- **Don't use**: original Path A (factual errors), Path B alone (deflective), Path C alone if §15 framing fixes aren't also applied (count transparency without surrounding-paragraph clarity is partial).

#### Coordination with §12 framing calibration

All Path A variants (revised, extended) and the discussion above are consistent with §12's framing-calibration note: NMF $k = 3$ has "weak prognostic content" (HR per SD 1.11; adj $p = 0.03$; cv-C-index 0.55), not "no prognostic content." The extended Path A explicitly says "weakly, $p = 0.03$" at $k = 3$; this language matches the §12 framing fix and the §15-flagged calibration. Apply consistently across §10, §12, and the rest of the manuscript.

#### Verification step before applying any Path A variant

Read Fig. S4 directly to confirm which rank each diagnostic favors. Currently known from the manuscript:
- **Elbow**: $k = 5$ (per line 355: "elbow-selected $k = 5$")
- **DeSurv BO max within tested range**: $k = 7$ (per line 609; this is the joint $k \times \alpha$ surface max)

Need to verify from Fig. S4:
- **Cophenetic correlation**: which $k$?
- **Silhouette width**: which $k$?
- **Reconstruction residual / elbow on residuals**: which $k$? (May be the same as the "elbow" criterion above, or a different elbow on a different metric.)

If the diagnostics actually do converge on a single rank (against my current assumption), the extended Path A draft needs to be calibrated to that data. If they diverge as assumed, the draft stands.

---

### Why this matters

The sensitivity analysis at `paper/04_results_REVISED.Rmd:609` reports significance counts for $k = 3$ vs. $k = 7$ (the BO global maximum), but does not specifically engage with what standard NMF ($\alpha = 0$) does at intermediate ranks $k = 5$ and $k = 9$. Inspection of the underlying matrices (`results/cv_grid/unadj_p_270_matrix.rds` and `adj_p_270_matrix.rds`, verified 2026-05-04) shows that **standard NMF at $k = 5$ achieves a stronger PurIST/DeCAF-adjusted p-value in external validation than DeSurv at $k = 3$ at any supervision strength**:

| $k$ | Unadj. p (α=0) | PurIST/DeCAF-adj. p (α=0) | Significant after adjustment? |
|---:|---:|---:|---|
| 2 | 0.8088 | 0.5639 | No |
| 3 | 0.0009 | 0.0304 | Yes |
| **5** | **0.0000** | **0.0005** | **Yes — most significant α=0 row** |
| 7 | 0.8784 | 0.5135 | No |
| **9** | **0.0000** | **0.0003** | **Yes — comparably strong** |

For comparison, DeSurv at $k = 3$ achieves a best PurIST/DeCAF-adjusted p of **0.0028** (at $\alpha = 0.55$), which is *less* significant than standard NMF at $k = 5$ (0.0005) or $k = 9$ (0.0003).

The data is in SI Appendix Table S3, so it is technically not hidden — but the main-text narrative compares $k = 3$ vs. $k = 7$ (the BO global maximum), and the count-based metric ("$k = 3$ achieved significance across [N] of [N] supervision strengths, $k = 7$ at a comparable [N] of [N]") conflates $\alpha = 0$ (standard NMF) with $\alpha > 0$ (DeSurv) and never specifically addresses $k = 5$. A reviewer who reads SI Table S3 will notice this in minutes and the main text has no ready answer.

The relationship is **non-monotonic** in $k$ for standard NMF (sig at $k = 3$, $5$, $9$; not sig at $k = 2$, $7$), so a defense cannot rest on "we showed standard NMF doesn't improve at higher $k$" — that's true for $k = 7$ specifically, not for $k = 5$ or $k = 9$. The $k = 7$ comparison was selected because it was the BO global maximum, but it has the side effect of selecting NMF's *worst* high-$k$ result and skipping its *best*.

### What the paper's defensible argument actually is

The argument that survives this finding (but is not currently made in the main text):

1. **Principled rank selection is the contribution, not "DeSurv beats NMF at every $k$."** Standard NMF at $k = 5$ succeeds *in retrospect*, but a researcher using only standard NMF would not have known to pick $k = 5$. Fig. S4 shows that standard NMF's own rank-selection diagnostics (cophenetic, silhouette, residuals) variously favor $k \in \{2, 3, 7\}$ — exactly the configurations where the unsupervised baseline fails (adjusted p = 0.56, 0.03, 0.51). DeSurv with cross-validated C-index reliably picks an interpretable, generalizable rank.
2. **Factor interpretability and parsimony are part of the contribution.** A 5-factor unsupervised solution does not isolate the Classical-tumor + iCAF coupling that DeSurv's 3-factor solution delivers. External *p-value* and biological interpretability are different things — the paper's downstream sections lean on the latter.
3. **The cohort-generalization claim is empirical, not structural.** DeSurv at $k = 3$ generalizes; standard NMF at $k = 5$ also generalizes here. The argument cannot rest on "DeSurv generalizes better in absolute terms"; it has to rest on "DeSurv generalizes via a principled, interpretable, parsimonious decomposition that standard NMF does not deliver via its own selection workflow."

### Path A — Explicit acknowledgment in main text

Add a sentence to the sensitivity-analysis paragraph (`paper/04_results_REVISED.Rmd:609`) acknowledging the $k = 5$ result and reframing the contribution as principled rank selection.

**Draft:**

> Standard NMF at $k = 5$ also achieved external significance (PurIST/DeCAF-adjusted $p = 0.0005$), comparable to or stronger than DeSurv on this metric. However, standard NMF's own rank-selection diagnostics (Fig. S4) do not converge on $k = 5$; the cophenetic, silhouette, and residual criteria variously favor $k \in \{2, 3, 7\}$, where the unsupervised baseline fails ($p = 0.56$, $0.03$, $0.51$ adjusted). The contribution of survival-supervised model selection is therefore principled, generalizable rank choice — not unconditional improvement over standard NMF at every rank a researcher could in principle have tried.

**Word delta:** ~+65 words.

**Pros:**
- Honest about the data; converts a potential reviewer ambush into a clear methodological argument
- Aligns the manuscript narrative with what's actually in SI Table S3
- Strengthens the "principled rank selection" framing as the paper's actual contribution
- A reviewer who reads SI Table S3 will see the main text has already addressed it — credibility win

**Cons:**
- Surfaces a finding ($k = 5$ NMF stronger than DeSurv $k = 3$) that some readers might over-weight relative to the paper's broader argument
- ~+65 words in an already-long Results paragraph
- Requires careful framing — "the contribution is principled rank selection" must be clearly articulated or the acknowledgment reads as a concession without a counter

### Path B — Reframe count metric, do not directly engage with $k = 5$

Adjust the existing $k = 3$ vs. $k = 7$ framing to emphasize that the comparison is across *DeSurv tuning configurations* rather than *across all $k$*, sidestepping the $\alpha = 0$ data without explicitly engaging with it.

**Draft:**

> Across DeSurv supervision strengths ($\alpha > 0$), $k = 3$ achieved significance in [N] of [N] configurations, with the parsimonious solution outperforming both larger and smaller ranks after adjustment for known classifiers. Standard NMF ($\alpha = 0$) results across $k$ are reported in SI Appendix, Table S3.

**Word delta:** ~+10–20 words (depending on how the existing sentence is restructured).

**Pros:**
- Minimal word delta
- Doesn't surface the $k = 5$ unsupervised result in the main text
- Frames the existing analysis around what DeSurv tuning achieves

**Cons:**
- Leaves the $k = 5$ NMF result undefended in the main text
- A reviewer who reads SI Table S3 will see the strongest unsupervised result is unaddressed and may take the framing as evasive
- Does not articulate the "principled rank selection" argument that survives the $k = 5$ finding
- "Standard NMF results across $k$ are reported in SI Appendix" is a soft pointer — a reviewer can interpret this as deflection rather than transparency

### Application checklist (whichever path)

1. Decide between Path A and Path B based on co-author preferences and risk tolerance for reviewer pushback.
2. Apply the chosen draft to the sensitivity-analysis paragraph at `paper/04_results_REVISED.Rmd:609`.
3. If Path A: ensure the framing of "principled rank selection" is internally consistent with §5 (which adds the prior-literature connection on inconsistent NMF rank-selection diagnostics) — the two arguments reinforce each other.
4. If Path B: be aware that the "Standard NMF results across $k$ are reported in SI Appendix, Table S3" pointer must be accurate (verify the table reference) and that this path does not preempt the reviewer question — it only minimizes main-text exposure.
5. Re-render with `make paper` and confirm the sensitivity-analysis paragraph still reads cleanly with the chosen revision.

### Underlying data for reference (not for inclusion in manuscript)

This data was extracted on 2026-05-04 from `results/cv_grid/unadj_p_270_matrix.rds` and `adj_p_270_matrix.rds`. Recorded here so future revisions don't need to re-derive it:

```
=== Unadjusted p-values (top-270 genes, by k × α) ===
    a=0.00 a=0.25 a=0.35 a=0.55 a=0.75 a=0.85 a=0.95
K=2 0.8088 0.0253 0.0047 0.0003 0.0001 0.0000 0.0000
K=3 0.0009 0.0009 0.1972 0.0000 0.0000 0.0169 0.0000
K=5 0.0000 0.0000 0.0000 0.0000 0.0004 0.1134 0.1205
K=7 0.8784 0.0060 0.0000 0.0014 0.1112 0.0028 0.0272
K=9 0.0000 0.0392 0.0063 0.0000 0.5887 0.3075 0.0091

=== PurIST/DeCAF-adjusted p-values (top-270 genes, by k × α) ===
    a=0.00 a=0.25 a=0.35 a=0.55 a=0.75 a=0.85 a=0.95
K=2 0.5639 0.0376 0.0379 0.0194 0.0603 0.0113 0.0025
K=3 0.0304 0.0302 0.5186 0.0028 0.0284 0.2629 0.0095
K=5 0.0005 0.0046 0.0008 0.0187 0.0727 0.6045 0.4471
K=7 0.5135 0.4713 0.0006 0.0408 0.6756 0.1314 0.1304
K=9 0.0003 0.2991 0.2430 0.0048 0.9427 0.8672 0.2628
```

---

## 19. Surface plateau argument and "NMF needs BO" argument in Results; foreshadow principled-rank-selection contribution claim in abstract/intro

### Why this matters

The Discussion paragraph at `paper/05_discussion_REVISED.Rmd:10` contains five distinct arguments doing heavy contribution-claim work:

1. **Plateau argument** ($k = 3$ to $12$ cv-C-index plateau on DeSurv's supervised dimension) — purely empirical
2. **"NMF needs DeSurv's BO" argument** (NMF only reaches $k = 7$ via Bayesian optimization with $\alpha = 0$; its own diagnostics don't converge there) — sharp methodological argument
3. **Diagnostic-by-diagnostic mapping** (cophenetic / silhouette / residual don't point to $k = 7$; elbow → $k = 5$) — empirical detail
4. **Robustness comparison** ($k = 3$ at 5 of 7 vs. $k = 7$ at 2 of 7) — already in Results line 609
5. **"Model selection framework is part of the methodological contribution"** — explicit contribution-claim framing

Three of these are stronger when surfaced earlier rather than waiting for Discussion. **Combined effect**: a reader who hits Discussion with the contribution claim already foreshadowed reads the Discussion paragraph as empirical confirmation; a reader who hits it as new framing reads it as defensive after-the-fact justification (especially in light of §10's NMF-at-$k = 5$ acknowledgment).

### Audit table — where each contribution-claim layer currently lives

| Contribution-claim layer | Abstract | Significance | Intro close (§3 post-review) | Results | Discussion |
|---|:---:|:---:|:---:|:---:|:---:|
| Biological insight (tumor-stroma coupling) | ✓ | ✓ | ✓ | ✓ (line 261) | ✓ (implicit) |
| Theoretical extension (Cook + IB + 4 constraints) | — | — | partial (line 19 raises open question) | — | ✓ (post §17) |
| Empirical demonstration (5 cohorts, transportable) | ✓ | ✓ | — | ✓ | ✓ |
| **Plateau argument ($k = 3$–$12$)** | — | — | — | **— (gap)** | ✓ |
| **NMF-needs-BO argument** | — | — | — | **— (gap)** | ✓ |
| **Principled rank selection (model selection framework as contribution)** | — | — | — | partial (§10 Path A acknowledgment) | ✓ |

The bottom three rows are gaps — content that belongs earlier but currently appears only in Discussion (or only at the §10 Path A reframing point).

### Three recommendations

#### 19a. Add plateau range to Results around line 206 (~+15 words)

**Current (`paper/04_results_REVISED.Rmd:206`):**

> "The resulting C-index surface identified a coherent region of high performance (Fig. \ref{fig:sim}D). Model selection followed the one-standard-error rule (Methods)..."

**Proposed:**

> "The resulting C-index surface identified a coherent region of high performance — **a broad plateau spanning $k = 3$ through $k = 12$ along DeSurv's supervised dimension** (Fig. \ref{fig:sim}D), indicating that three supervised factors capture nearly as much prognostic information as twelve. Model selection followed the one-standard-error rule..."

Net effect: gives Results-level empirical grounding for why the 1-SE rule could legitimately select $k = 3$. The plateau range is the load-bearing fact; surfacing it here makes the 1-SE choice in the next sentence land empirically rather than by author preference. Also reduces redundancy with the Discussion's plateau sentence (which can become a brief callback rather than introducing the fact for the first time).

#### 19b. Add "NMF needs BO" argument to Results around line 355 (~+15 words, post-2026-05-05 streamlining)

**Current (`paper/04_results_REVISED.Rmd:355`):**

> "Standard NMF at elbow-selected $k = 5$ and BO-selected $k = 7$ substantially improved over NMF at $k = 3$, with $k = 7$ matching or exceeding DeSurv's per-cohort C-index in most cohorts."

**Proposed (em-dash qualifying clause within existing sentence, post-paragraph-preview):**

> "Standard NMF at elbow-selected $k = 5$ and BO-selected $k = 7$ substantially improved over NMF at $k = 3$, with $k = 7$ matching or exceeding DeSurv's per-cohort C-index in most cohorts **— though $k = 7$ is reached only through Bayesian optimization with $\alpha = 0$, since NMF's own rank-selection heuristics (cophenetic, silhouette, residual; Fig. S4) do not converge on this rank.**"

**Why this revision (post-paragraph-preview, 2026-05-05):** The original draft inserted the "NMF needs BO" content as a standalone sentence between the existing comparison and the "However, this concordance gain requires more than twice as many factors" sentence. Paragraph-level preview revealed the standalone insertion creates referent ambiguity: the "However" was supposed to refer to NMF matching DeSurv per-cohort C-index, but with the inserted sentence in between, "However" reads as connected to the inserted sentence instead. Folding the BO-borrowing point into the existing sentence as a qualifying clause (em-dash) eliminates the structural digression and saves ~15 words.

**Earlier post-review draft (now superseded):** *standalone sentence after the existing comparison: "NMF only reaches $k = 7$ through Bayesian optimization with $\alpha = 0$ — applying DeSurv's tuning infrastructure to an unsupervised model — since NMF's own rank-selection heuristics (cophenetic, silhouette, residual; Fig. S4) do not converge on $k = 7$."*

**Coordination with §20a:** §20a proposed a sharper version of this same content. After this paragraph-preview streamlining, §20a's sharper conditional framing ("matches only when DeSurv's tuning framework applied") can be folded into the same em-dash qualifying clause if desired, OR §20a can be applied at a later sentence in the same paragraph as a separate beat. See §20 entry for the post-preview alignment.

#### 19c. Foreshadow principled-rank-selection contribution claim in abstract or intro close (~+10 words)

This is the §18-precursor recommendation from the prior turn (offered but not filed at that point). The contribution claim **"the model selection framework is part of the methodological contribution, not merely the survival gradient"** currently lives only in Discussion. It should be foreshadowed at one of: abstract closer, significance statement, or intro close. Two options:

**Option 1 — Add to §3 intro close (preferred — Discussion-aligned framing):**

Current §3 close (post-review):
> "DeSurv provides a general framework for extracting outcome-relevant programs from bulk tumor transcriptomes wherever unsupervised factors conflate biological signal with compositional noise, offering a route to reproducible prognostic biomarkers and biological insight into tumor-stroma coupling in cancers like PDAC where compositional variation dominates expression."

Revised:
> "DeSurv provides a general framework for extracting outcome-relevant programs from bulk tumor transcriptomes wherever unsupervised factors conflate biological signal with compositional noise. **Survival-supervised cross-validation simultaneously resolves the rank-selection ambiguity that unsupervised diagnostics leave open**, offering a route to reproducible prognostic biomarkers and biological insight into tumor-stroma coupling in cancers like PDAC where compositional variation dominates expression."

**Option 2 — Add to §1 abstract closer (lighter alternative):**

Current §1 closer (post-review):
> "By prioritizing transcriptional programs that survive cross-cohort generalization, DeSurv yields reproducible prognostic biomarkers..."

Revised:
> "By prioritizing transcriptional programs that survive cross-cohort generalization **and selecting model complexity by cross-validated outcome prediction rather than unsupervised heuristics**, DeSurv yields reproducible prognostic biomarkers..."

Either sets up the Discussion claim ("model selection framework is part of the contribution") as confirmation rather than novel framing.

### What NOT to move

- **Diagnostic-by-diagnostic mapping** (elbow → $k = 5$; others → not $k = 7$): leave in Discussion. Interpretive elaboration that benefits from Discussion-level context. 19b telegraphs the relevant subset.
- **1-SE rule philosophy** ("This parsimony is a property of the 1-SE rule rather than an intrinsic consequence of supervision..."): leave in Discussion. Methodological commentary belongs there.
- **Robustness comparison** ($5/7$ vs $2/7$): already in Results line 609; current Discussion mention is appropriate elaboration.

The principle: **empirical observations and contribution claims foreshadow earlier; interpretive philosophy stays in Discussion.**

### Coordination notes

- **19a + 19b strengthen §10 Path A.** With the plateau argument and NMF-needs-BO argument inline in Results, §10 Path A's reframing of the contribution as principled rank selection lands as the obvious interpretation rather than as a defensive move. The two Results additions essentially do §10 Path A's work inline at the moment each comparison is made — and §10 Path A's main-text-acknowledgment-then-reframe move can be lighter or even folded into 19b.
- **19c is the same as the §18-precursor** I'd flagged in the prior-turn audit. Filing it here consolidates the cross-paper contribution-claim distribution audit into one section.
- **19a's plateau range partly subsumes the Discussion's plateau sentence.** When applying §17 (Discussion reframing), check that the plateau range isn't introduced twice — the Discussion can become a brief callback ("the broad plateau noted in Results...") rather than re-introducing the fact.

### Cost/benefit

- **Word delta:** +15 (19a) + 30 (19b) + 10 (19c) = ~+55 words across three locations.
- **No new analyses, no bib changes.** All facts already in cached results or existing main-text claims.
- **Reader benefit:** the contribution claim's empirical and methodological pieces foreshadow earlier in the paper; the Discussion's heavy contribution-claim paragraph reads as confirmation rather than defensive after-the-fact justification.
- **Risk:** very low. Substantive content unchanged; only the distribution of contribution-claim foreshadowing shifts.

### Application checklist

1. Apply 19a to `paper/04_results_REVISED.Rmd:206` (plateau range in Results, near Fig. 2D introduction).
2. Apply 19b to `paper/04_results_REVISED.Rmd:355` (NMF-needs-BO sentence after the existing $k = 5$ / $k = 7$ statement).
3. Choose 19c Option 1 (intro close) or Option 2 (abstract closer); apply to corresponding location.
4. When applying §17 (Discussion reframing) and the §17-paragraph plateau sentence, check that 19a's plateau range isn't introduced twice — Discussion should become a callback to the Results-level statement.
5. Re-render with `make paper` and confirm all three additions read cleanly with their surrounding text.

---

## 20. Surface the two-part contribution framing earlier — model selection framework as contribution + NMF's matching performance is conditional on borrowing DeSurv's infrastructure

### Update (2026-05-05): paragraph-level preview adds two coordination notes

**Coordination with §19b (em-dash streamlining):** §19b was originally a standalone sentence at line 355; paragraph-level preview revised it to an em-dash qualifying clause within the existing sentence. §20a's "sharper conditional framing" can now either fold into the same em-dash clause (compact), OR appear as a separate beat at a later sentence in the same paragraph (preserves §20a's "matches only when" emphasis). The em-dash-fold variant is more concise; the separate-beat variant is more rhetorically explicit. Editor's choice.

**Discussion sentence 7 becomes callback after §20b/§20c lands:** §20b (intro close two-part contribution framing) or §20c (abstract closer) surfaces the two-part contribution claim earlier in the paper. After either is applied, the Discussion sentence at line 10 *"The model selection framework — joint BO over rank and hyperparameters with the 1-SE rule — is thus itself part of the methodological contribution, not merely the survival gradient"* becomes a callback rather than first-introduction. **Recommended adjustment when §20b/§20c is applied:** light reframing of Discussion sentence 7 as callback. Possible:

> "As foreshadowed in the Introduction, the model selection framework — joint BO with the 1-SE rule — is itself part of the methodological contribution: parsimony is a property of the 1-SE rule, but supervision creates the flat concordance surface under which it can operate effectively."

This wording also folds in §21's merged-sentences-7-8 streamlining (see §21's update). The "As foreshadowed in the Introduction" phrase explicitly bridges to §20b's intro framing. Apply §20b/§20c + §21 together to get the coordinated effect.

### Why this matters

The Discussion sentence at `paper/05_discussion_REVISED.Rmd:10` contains a sharp rhetorical move:

> "**The model selection framework — joint BO over rank and hyperparameters with the 1-SE rule — is thus itself part of the methodological contribution, not merely the survival gradient.**"

This sentence does two distinct rhetorical jobs:

1. **Names the model selection framework specifically** (joint BO + 1-SE rule) — concrete description.
2. **Claims it as a separate contribution component** alongside the survival gradient — the "not merely" construction explicitly preempts the reading that DeSurv = "just supervised NMF."

Combined with the §19b argument (NMF only reaches $k = 7$ via Bayesian optimization with $\alpha = 0$), this exposes the underlying rhetorical structure: **NMF's matching performance at $k = 7$ is conditional on borrowing DeSurv's tuning infrastructure**. NMF in its native usage (cophenetic / silhouette / residual / elbow rank-selection diagnostics) never gets to $k = 7$. We have effectively been giving NMF the benefit of DeSurv's model-selection apparatus when comparing the two methods at NMF's "best" configuration.

This observation is two things at once:

- **Positive claim about DeSurv:** the model selection framework (joint BO + 1-SE) is part of the methodological contribution, not just the survival gradient.
- **Negative claim about NMF:** NMF's matching performance is conditional on borrowing DeSurv's infrastructure.

Both currently live only in the Discussion. **Surfacing them earlier converts what looks like a concession ("NMF at $k = 7$ matches") into a strength ("NMF only matches with DeSurv's infrastructure")** — a rhetorically much stronger framing that makes §10's NMF-at-$k = 5$ acknowledgment less of a vulnerability.

### Distinction from §19c (lighter version of the same idea)

§19c foreshadows "principled rank selection" generally: *"Survival-supervised cross-validation simultaneously resolves the rank-selection ambiguity that unsupervised diagnostics leave open"* (intro close) or *"selecting model complexity by cross-validated outcome prediction rather than unsupervised heuristics"* (abstract). These are lighter touches.

§20 is the more explicit two-part framing — names "joint BO over rank and hyperparameters with the 1-SE rule" concretely, and asserts the framework as a separate contribution component. §19c and §20 are complementary, not redundant: pick one for any given location based on how explicit you want to be.

### Three placement options for surfacing the two-part contribution earlier

#### 20a. Sharpen §19b's Results sentence to foreground the conditionality

Current §19b draft (already filed):

> "NMF only reaches $k = 7$ through Bayesian optimization with $\alpha = 0$ — applying DeSurv's tuning infrastructure to an unsupervised model — since NMF's own rank-selection heuristics (cophenetic, silhouette, residual; Fig. S4) do not converge on $k = 7$."

Sharper variant:

> "**Standard NMF matches DeSurv's per-cohort C-index only at $k = 7$, and reaches $k = 7$ only when DeSurv's tuning framework is applied to it ($\alpha = 0$ in the joint BO over rank and hyperparameters)**. NMF's own rank-selection heuristics (cophenetic, silhouette, residual, elbow; Fig. S4) do not converge on $k = 7$."

Net delta vs §19b: ~+10 words. The sharper variant explicitly names the conditional ("matches only when... is applied") rather than just mentioning the BO usage. Same data, more rhetorically explicit.

#### 20b. Foreshadow the model selection framework as a separate contribution in intro close (~+15 words on top of §19c Option 1)

Replaces §19c Option 1 with a more explicit version:

> "DeSurv provides a general framework for extracting outcome-relevant programs from bulk tumor transcriptomes wherever unsupervised factors conflate biological signal with compositional noise. **The methodological contribution is twofold: a survival-supervised gradient on the gene-program matrix $W$, and a joint Bayesian optimization framework over rank and hyperparameters (with the 1-SE rule for parsimonious selection) that resolves the rank-selection ambiguity unsupervised diagnostics leave open.** This offers a route to reproducible prognostic biomarkers and biological insight into tumor-stroma coupling in cancers like PDAC where compositional variation dominates expression."

Net delta vs §19c Option 1: ~+20 words. Names both contribution components explicitly. Trade-off: heavier intro close, but matches the Discussion's "two-part contribution" framing.

#### 20c. Add to abstract closer with explicit two-part framing (~+15 words on top of §19c Option 2)

Replaces §19c Option 2 with a more explicit version:

> "By prioritizing transcriptional programs that survive cross-cohort generalization **and resolving rank-selection ambiguity through joint Bayesian optimization with one-standard-error model selection**, DeSurv yields reproducible prognostic biomarkers for the bulk-transcriptomic cohorts that dominate clinical cancer research, addressing the variance-prognosis misalignment that limits unsupervised factorization."

Net delta vs §19c Option 2: ~+12 words. Names the model selection framework explicitly without going full two-part-contribution exposition. Lighter than 20b.

### Recommendation

**Apply 20a unconditionally** (replacing §19b's Results draft). The sharper conditional framing ("matches only when DeSurv's tuning framework is applied") is the load-bearing rhetorical move that converts §10's NMF-at-$k = 5$ acknowledgment from concession to expected qualification. Cost is ~+10 words on top of §19b; benefit is significant rhetorical strengthening.

**Choose 20b or 20c (or §19c Option 1 / Option 2) based on word budget tolerance**:

- **20b** if intro close can absorb +20 words and you want maximum explicitness about the two-part contribution.
- **20c** if abstract has room for the model-selection-framework specification.
- **§19c Option 1 or 2** (lighter wording) if 20b/20c feel too heavy for the location.

**Don't apply both 20b and 20c** — that would duplicate the model-selection-framework framing in two places. Pick one location for the explicit version.

### Coordination notes

- **20a directly strengthens §10 Path A.** §10 Path A's central reframing — that DeSurv's contribution is principled rank selection — is much sharper when 20a's "matches only when DeSurv's tuning framework is applied" is in Results. The §10 Discussion-level acknowledgment then reads as confirmation, not justification.
- **20b/20c overlap with §19c.** Pick one path per location. If applying 20b, skip §19c Option 1 (intro close). If applying 20c, skip §19c Option 2 (abstract).
- **§17 Discussion plateau callback.** When applying §17 (Discussion theoretical-extension reframing), the existing Discussion sentence about model selection framework as contribution can become slightly lighter — readers will already have seen the two-part framing earlier. The Discussion's job becomes elaboration, not first-introduction.

### Why this strengthens, doesn't weaken, the paper

The "leg up" observation could superficially read as "we admit we gave NMF an artificial advantage and it still didn't match us." But the actual rhetorical move is the opposite: by making explicit that NMF's matching performance requires DeSurv's infrastructure, the paper foregrounds the value of *both* contribution components (gradient + framework), and turns NMF's borrowed-infrastructure performance from a weakness in DeSurv's claim into a strength in DeSurv's contribution scope.

Symbolically:
- **Without 20a/20b/20c**: "NMF at BO-selected $k = 7$ matches DeSurv per-cohort." Reader: "So NMF does just as well? What's the point?"
- **With 20a (and 20b or 20c)**: "NMF matches DeSurv only when given DeSurv's tuning framework, and even then with twice as many factors." Reader: "So DeSurv's contribution includes the tuning framework that lets *any* method find a generalizable rank — and DeSurv's supervised version does it more parsimoniously."

The second framing is a stronger contribution claim. The data is unchanged.

### Cost/benefit

- **Word delta:** +10 (20a, replaces §19b's +30 → net same word count) + 15-20 (20b or 20c, replaces §19c's +10 → net +5-10 over §19c).
- **No new analyses, no bib changes.**
- **Reader benefit:** strongest available reframing of the comparison structure; converts NMF's matching performance from latent vulnerability into explicit contribution.
- **Risk:** very low. Substantive content unchanged; only the rhetorical framing of contribution and conditional matching shifts.

### Application checklist

1. Apply 20a to `paper/04_results_REVISED.Rmd:355` (replaces §19b's draft with sharper conditional framing).
2. Choose 20b (intro close) or 20c (abstract closer); apply to corresponding location, replacing §19c's lighter version.
3. When applying §17 (Discussion reframing), check that the Discussion's "model selection framework is itself part of the methodological contribution" sentence becomes elaboration rather than first-introduction — a callback to the earlier two-part framing.
4. Re-render with `make paper` and confirm the contribution-claim arc reads consistently across abstract/intro → Results → Discussion.

---

# Cluster D — Type III / joint-vs-individual / "weak not no" framing calibration

Two sections, plus calibration threads embedded in §10, §15, and §29. The recurring rule: replace "no prognostic content" / "negligible contribution" wording with "weak prognostic content" / "negligible per-factor (Type III) contribution" wording where it appears, distinguishing per-factor metrics from joint-model metrics where relevant.

**Additional candidate sites identified during cluster-expansion analysis:**
- Results line 261 ("no NMF factor achieves comparable survival contribution at the same rank") — Type III qualifier could be added if §12 is applied.
- SI Section 10 (formal $\Delta\ell$ definition) — could include a brief sentence noting per-factor metric may underestimate joint contribution. Mirrors §29c's qualifier.
- SI line 1582 ("standard NMF continues to distribute survival contribution negligibly across all factors") — extend §29c-style qualifier here.

---

## 8. Address the precision/C-index gap in the simulation Results

### Why this matters

The simulation Results in `paper/04_results_REVISED.Rmd:244` report DeSurv vs. standard NMF as median C-index 0.839 vs. 0.724 (Δ = 0.115) alongside a much wider precision gap (near-zero vs. 0.50). A reader can reasonably wonder why the C-index gap is so much smaller than the precision gap, and may walk away thinking "DeSurv is marginally better" instead of "standard NMF can predict moderately well without recovering the right biology." The current text does not close this loop. A short clarifying sentence is enough to do so.

### Methodology note (for our reference; no manuscript change needed)

The reported C-index is a **held-out test-set evaluation within each simulation replicate**, not the BO inner-loop cross-validated value. Verified in `code/sim_helpers.R`:

- `split_simulation_samples()` (line 290) creates the train/test split per replicate
- `compute_dataset_cindex()` (line 408) computes C-index from a fitted model applied to a dataset
- Lines 1404–1405 call it on `processed_train` and `processed_test` separately
- Line 1228 reports `test_cindex` as the headline `cindex`

Fig 1A caption already reads *"Test-set concordance index"*, so the manuscript is correct on this point — no edit needed.

### Proposed edit (clarifying sentence in Results)

After the existing C-index sentence at `paper/04_results_REVISED.Rmd:244` (*"...DeSurv achieved median C-index 0.839 versus 0.724 for standard NMF (Δ = 0.115, paired Wilcoxon P = ...).")*, append one sentence. Two drafts at different lengths:

**Option A — Compact (~+50 words):**

> The smaller gap in C-index than in precision reflects a key distinction: standard NMF's factor scores can predict survival moderately well via incidental linear combinations of gene expression, even when the factors themselves do not recover the true prognostic genes. C-index alone is therefore a permissive criterion — DeSurv's advantage is most visible at the gene level, where biological interpretation and cross-cohort transportability are determined.

**Option B — Tighter (~+30 words):**

> The C-index gap is smaller than the precision gap because moderate survival prediction can arise from incidental linear combinations of gene expression in the factor scores, even when factors do not recover the true prognostic genes — a permissive criterion that masks the gene-level failure precision exposes.

Either version is defensible. Option B is preferred if the surrounding paragraph is already long; Option A if the paragraph has room and clearer reader-facing prose is wanted.

### Optional soft bridge to external validation (one sentence, defensible)

The simulation tests within-cohort sample-level generalization (held-out test set drawn from the same generative process as training). Cross-cohort transportability is tested separately in the PDAC validation analyses on five independent cohorts. Adding a single bridging sentence makes the relationship between the two halves of the paper explicit without making any new claim:

> These results characterize within-cohort sample-level generalization; cross-cohort transportability is assessed separately in the PDAC external validation analyses.

Place at the end of the simulation paragraph, after the clarifying sentence above. Purely factual; doesn't claim the simulation predicts external validation outcomes.

### What NOT to do (avoid these overreach traps)

These were considered during drafting and explicitly rejected. Recording here so future revisions don't reintroduce them:

- **Don't cite `aran2015systematic` or `rashid2020purity`** to support a "compositional-prognosis confounding" mechanism here. Those citations are used in the intro to argue the *opposite* — that PDAC compositional factors (exocrine, low purity) are outcome-neutral. Invoking them here would contradict the intro's positioning and is unnecessary: the simulation is constructed with outcome-neutral high-variance backgrounds, so compositional confounding is excluded by design and the result is explained without that mechanism.
- **Don't cite Venet et al. 2011** ("Most random gene expression signatures are significantly associated with breast cancer outcome") as a parallel demonstration. The analogy with random gene panels in breast cancer requires careful defense — different data type, different mechanism, not directly comparable to factor-based methods — and is not needed for the point being made. Self-contained framing is cleaner.
- **Don't expand into deep linear-algebra mechanism exposition** (column space of $W$, projection $W\beta$, subspace overlap). The existing "incidental linear combinations" wording in Options A and B above is sufficient. Going deeper invites reviewers to ask for formal characterizations of when subspace overlap is large/small, which is a can of worms.
- **Don't make strong predictive claims** about how subspace-overlap factors will behave in independent cohorts. The bridge sentence above ("characterize within-cohort sample-level generalization; cross-cohort transportability is assessed separately") is the maximum defensible scope. Stronger framings ("subspace-overlap factors are expected to degrade in independent cohorts") are mechanistically sensible but invite reviewer requests for analysis to back the prediction up.

### Cost/benefit

- **Word delta:** ~+30–50 words for the clarifying sentence + ~+15 words for the optional bridge = ~+45–65 words total. Within budget.
- **No bib changes.** No new analyses required.
- **Reader benefit:** closes a comprehension gap on the simulation's headline numbers and makes the precision/C-index discrepancy land as the paper's central methodological point rather than an apparent contradiction.
- **Risk:** very low if the drafts above are used verbatim. The "What NOT to do" list captures the specific overreach traps that were considered and rejected.

### Application checklist

1. Choose Option A or Option B for the clarifying sentence; insert at end of the paragraph in `paper/04_results_REVISED.Rmd:244` (immediately after the existing C-index numbers sentence).
2. Optionally add the soft bridge sentence after the clarifying sentence.
3. Re-render with `make paper` and confirm the simulation paragraph still reads cleanly.

---

## 12. Report overall NMF vs. DeSurv cv-C-index at matched rank to pre-empt multicollinearity / Type-III concern in the per-factor $\Delta\ell$ analysis

### Why this matters

The per-factor contribution paragraph at `paper/04_results_REVISED.Rmd:259` reports survival contribution as the change in Cox partial log-likelihood ($\Delta\ell$) when each factor is dropped from the $k$-factor Cox model. This is structurally a Type-III-like test — it measures the *unique* contribution of each factor conditional on all others. When covariates are correlated, Type-III tests can show small individual contributions even when the joint contribution is meaningful: dropping one correlated factor transfers its signal to the remaining factors, leaving each factor's $\Delta\ell$ small while the overall model remains prognostic.

A reviewer reading "All three NMF factors had negligible $\Delta\ell$" can reasonably ask: *is the NMF model overall as weak as the per-factor decomposition suggests, or is multicollinearity between NMF factors hiding a meaningful joint contribution?* The current main text does not engage with this question. The reader has to either trust the per-factor $\Delta\ell$ at face value or piece the answer together from the BO surface (Fig. 2D), which shows the cv-C-index landscape but does not isolate the $\alpha = 0$, $k = 3$ point.

### Empirical answer (verified 2026-05-04)

The overall model-level cv-C-index at matched rank empirically resolves the concern. Extracted directly from the cached BO results in `results/precomputed/ntop_bo_50_300/desurv_bo_results_alpha0_tcgacptac.rds` and `desurv_bo_results_tcgacptac.rds`:

| Method (matched $k = 3$) | Best cv-C-index | Note |
|---|---:|---|
| Standard NMF ($\alpha = 0$) | **0.554** | Marginally above chance (0.5) |
| DeSurv (best $\alpha$) | **0.646** | Meaningfully prognostic |

**Interpretation:** the overall NMF $k = 3$ Cox model is barely prognostic. The near-zero per-factor $\Delta\ell$ across all three NMF factors is consistent with this overall weakness, **not** with multicollinearity hiding joint contribution. Multicollinearity can hide *attribution* between correlated covariates; it cannot manufacture joint signal that does not exist. The 0.092 cv-C-index gap (0.646 vs. 0.554) at matched rank is real model-level superiority, and the $\Delta\ell$ decomposition (D1 = 66.4; all NMF factors near zero) reflects this rather than being an artifact of the conditional metric.

For broader context, NMF cv-C-index by $k$ at $\alpha = 0$ (best across other hyperparameters):

| $k$ | NMF ($\alpha = 0$) cv-C-index |
|---:|---:|
| 2 | 0.530 |
| 3 | 0.554 |
| 5 | 0.622 |
| 7 | 0.643 |
| 9 | 0.642 |
| 11 | 0.646 |

For DeSurv (best $\alpha$ at each $k$): 0.573 (k=2), 0.646 (k=3), 0.645 (k=5), 0.655 (k=7), 0.652 (k=9). DeSurv at $k = 3$ matches NMF's best across any $k$ tested (NMF reaches 0.646 only at $k = 11$). This is broader context for §10 — the "NMF at higher $k$ does well" finding (§10) was about external-validation p-value, not training cv-C-index magnitude. In magnitude, DeSurv at $k = 3$ matches NMF's best at any rank.

### Proposed edit (Results, per-factor contribution paragraph)

Insert a single sentence in the per-factor contribution paragraph at `paper/04_results_REVISED.Rmd:259`. Recommended placement: between the qualitative summary ("All three NMF factors had negligible $\Delta\ell$ ... while DeSurv concentrated nearly all survival contribution into D1") and the quantitative breakdown ("Specifically, D1 explained 7% of variance ..."). This places the overall comparison as the framing for the factor-level decomposition rather than a defensive afterthought.

**Three drafts at increasing levels of explicitness:**

**Option A — Compact (~+25 words, preferred):**

> "These per-factor contrasts reflect overall model-level differences: at matched rank, the NMF $k = 3$ Cox model achieved a cross-validated C-index of 0.55 (marginal over chance), while DeSurv achieved 0.65."

**Option B — Names the multicollinearity concern explicitly (~+45 words):**

> "These per-factor contrasts reflect overall model-level differences rather than redistribution of signal across correlated factors: at matched rank, the NMF $k = 3$ Cox model achieved a cross-validated C-index of 0.55 — marginally above chance — versus DeSurv's 0.65, confirming that the near-zero per-factor $\Delta\ell$ for NMF reflects weak model-level prognostic content."

**Option C — Most defensive (~+55 words):**

> "These per-factor contrasts reflect overall model-level differences rather than Type-III artifacts of correlated factors: at matched rank ($k = 3$), the NMF Cox model achieved a cross-validated C-index of 0.55 (marginally above chance), versus DeSurv's 0.65. Multicollinearity between factors can redistribute conditional contribution but cannot manufacture joint signal absent from the overall model — confirming that NMF's near-zero per-factor $\Delta\ell$ reflects weak overall prognostic content."

### Recommended choice

Option A is preferred for most contexts. It addresses the concern empirically without invoking technical jargon ("Type III," "multicollinearity") that may itself raise questions. Option B is the right call if you want the multicollinearity rebuttal to be explicit (e.g., if a co-author has flagged this concern). Option C is most rigorous but reads heavy and may oversell the point. The empirical resolution (cv-C-index 0.55 vs. 0.65) is doing all the work in any of the three; the wording differences are about how loudly to advertise the rebuttal.

### Coordination with §5 and §10

This edit reinforces §5 (k-selection literature) and §10 (NMF at higher $k$):

- **§5 + §12 together**: §5 establishes that NMF rank-selection diagnostics give inconsistent guidance; §12 shows that even at the matched rank, NMF's overall model is weak. Together: NMF can't pick a good $k$ from its own diagnostics, *and* matched-rank NMF is itself weak — DeSurv's principled rank selection is the contribution, not just selection of a different $k$.
- **§10 + §12 together**: §10 acknowledges NMF at $k = 5$ achieves external significance; §12 establishes that even there, NMF's training cv-C-index (0.622) is comparable but not superior to DeSurv at $k = 3$ (0.646). The §10 "NMF wins on adjusted p-value" is about statistical significance with reasonable n, not magnitude. §12 closes the gap in the magnitude argument.

If applying §10 Path A (explicit acknowledgment of $k = 5$), §12 Option B or C is more compatible because both engage the unfavorable comparison directly. If applying §10 Path B (reframe count metric), §12 Option A is sufficient.

### Cost/benefit

- **Word delta:** ~+25 to +55 words depending on Option chosen.
- **No bib changes, no new analyses.** All numbers come from already-cached BO results.
- **Reader benefit:** pre-empts the multicollinearity / Type-III pushback entirely; turns the per-factor $\Delta\ell$ analysis from a vulnerable Type-III argument into a corroborated overall-model-level finding.
- **Risk:** very low. The cited cv-C-index values (0.55, 0.65) are training cross-validated, not external validation, so they're consistent with the BO-tuning context already described. No new claims about external generalization are made.

### Application checklist

1. Decide between Options A, B, and C based on co-author preferences (or skip entirely if §10 Path B is chosen and the multicollinearity concern is not flagged elsewhere).
2. Apply the chosen draft to `paper/04_results_REVISED.Rmd:259`, between the qualitative summary and the quantitative breakdown sentences.
3. Verify the cv-C-index values (0.55 for NMF, 0.65 for DeSurv at $k = 3$) by re-reading from the cached BO results — these are the rounded values from `desurv_bo_results_alpha0_tcgacptac.rds` (best at $k = 3$ = 0.554) and `desurv_bo_results_tcgacptac.rds` (best at $k = 3$ = 0.646). Round consistently with how cv-C-index values are reported elsewhere in the manuscript (typically 2 decimal places in body text, 3 decimal places in figures).
4. Re-render with `make paper` and confirm the per-factor contribution paragraph still reads cleanly.

### Underlying data for reference (not for inclusion in manuscript)

Extracted 2026-05-04 from `results/precomputed/ntop_bo_50_300/`:

```
NMF (alpha=0) cv-C-index by k (best across lambda/nu/ntop):
  k=2:  0.530
  k=3:  0.554   ← matched-rank comparison value
  k=5:  0.622
  k=7:  0.643
  k=9:  0.642
  k=11: 0.646   ← NMF overall best (alpha=0)

DeSurv cv-C-index by k (best across alpha/lambda/nu/ntop):
  k=2:  0.573
  k=3:  0.646   ← matched-rank comparison value
  k=5:  0.645
  k=7:  0.655
  k=9:  0.652

DeSurv overall best: k=7, cv-C-index=0.655
```

---

# Cluster E — Methods accuracy and citations

Six methods-section fixes, each pre-empting a specific reviewer comment about citation completeness or methodological accuracy. Mostly free wins; can be applied independently of each other and of the larger clusters.

**Additional candidate sites identified during cluster-expansion analysis:**
- Verify Methods Algorithm 1 cites NMF multiplicative updates (`lee1999learning`) — likely cited but worth confirming.
- Within-subject rank transformation step in Methods is uncited; check whether a standard reference applies. Probably non-blocking.
- SI Section 6 (Bayesian Optimization) should reference §6's BO citations after they are added to Methods, or stay silent on citations. Coordination check.

---

## 5. Connect the "inconsistent guidance across candidate ranks" observation to prior k-selection literature

### Why this matters

The Results currently contains this sentence (`paper/04_results_REVISED.Rmd:204`):

> Standard NMF diagnostics (reconstruction residuals, cophenetic correlation, mean silhouette width) yielded inconsistent guidance across candidate ranks (SI Appendix, Fig. S4).

As written, this reads as a novel observation we made on this dataset. That framing is rhetorically weaker than it could be. A reviewer can reasonably ask "did you just pick three diagnostics that happened to disagree on this one cohort?" The defensive framing is to ground the observation in the **established literature on NMF rank selection**, where metric disagreement is a documented phenomenon — not a quirk of our PDAC training data. Reframing it as "we observed the documented phenomenon" shifts the burden away from us having to defend why these three particular criteria conflict, and toward acknowledging this is a known limitation of the unsupervised toolkit that motivates the supervised approach.

The intro file's top note records that the recent revision *trimmed* a rank-selection pre-emption from P3 (1087 → 790 word reduction). That means the paper currently introduces the variance-vs-prognosis problem but does not foreground the rank-selection problem in the introduction at all. This subsection proposes (a) a surgical augmentation of the Results sentence as the minimum-risk fix, and (b) optionally a short clause in intro P2 to telegraph the issue upfront.

### Available citations (already in bib — no new entries needed)

Two relevant entries already exist in `references_30102025.bib`:

- **`brunet2004metagenes`** — Brunet, Tamayo, Golub, Mesirov. "Metagenes and molecular pattern discovery using matrix factorization," *PNAS* 101(12):4164–4169 (2004). The canonical introduction of cophenetic correlation as an NMF rank-selection criterion. Already cited in the intro (line 17, P2 of `02_introduction_REVISED.Rmd`).
- **`frigyesi2008nmf`** — Frigyesi & Höglund. "Non-negative matrix factorization for the analysis of complex gene expression data: identification of clinically relevant tumor subtypes," *Cancer Informatics* 6:275–292 (2008). Explicitly addresses NMF rank-selection limitations in cancer genomics and discusses the limitations of cophenetic correlation as a sole selection criterion. Currently in the bib but only cited in Results (per the citation-key inventory at top of `02_introduction_REVISED.Rmd`).

These two together cover both *introduction* of the standard metric and its *known failure modes* — a paired citation pattern that signals to reviewers that the authors know the methodological history.

**If a stronger pan-method case is wanted,** candidates worth adding (none currently in bib — would require new entries):
- Hutchins, Murphy, Singh, Graber. "Position-dependent motif characterization using non-negative matrix factorization," *Bioinformatics* 24(23):2684–2690 (2008). DOI: 10.1093/bioinformatics/btn526. Discusses quality-of-fit metrics for NMF rank.
- Tan & Févotte. "Automatic relevance determination in nonnegative matrix factorization with the β-divergence," *IEEE TPAMI* 35(7):1592–1605 (2013). DOI: 10.1109/TPAMI.2012.240. Tackles automatic rank selection, motivating the need for principled approaches.
- Stein-O'Brien et al. "Enter the matrix: factorization uncovers knowledge from omics," *Trends in Genetics* 34(10):790–805 (2018). DOI: 10.1016/j.tig.2018.07.003. Reviews dimensionality-reduction in genomics including NMF rank-selection challenges.

For a defensive intro hook, the two already-in-bib citations are sufficient and avoid bib churn.

### Proposed edits

#### 5a. Surgical: augment the Results sentence (preferred minimum)

**Current (`paper/04_results_REVISED.Rmd:204`):**

> Standard NMF diagnostics (reconstruction residuals, cophenetic correlation, mean silhouette width) yielded inconsistent guidance across candidate ranks (SI Appendix, Fig. S4).

**Proposed:**

> Standard NMF diagnostics (reconstruction residuals, cophenetic correlation, mean silhouette width) yielded inconsistent guidance across candidate ranks **— a known limitation of unsupervised rank-selection criteria, which can disagree on the same dataset [@brunet2004metagenes; @frigyesi2008nmf]** (SI Appendix, Fig. S4).

Word delta: ~+18 words. No new bib entries.

#### 5b. Optional: add a setup clause in introduction P2

This is the more substantive change — it telegraphs the rank-selection problem in the introduction so the Results observation lands as expected rather than as new information. Trade-off: partially re-inflates the P3 trim that was just done.

**Current sentence at end of P2 (`02_introduction_REVISED.Rmd:17`):**

> Identifying the prognostically relevant programs among the discovered factors requires extensive downstream evaluation: over a decade of PDAC subtyping was needed to converge on a robust Basal-like/Classical dichotomy, first defined through virtual microdissection [@moffitt2015virtual] and subsequently confirmed in independent cohorts [@rashid2020purity], single-cell analyses [@chansengyue2020transcription; @werba2023single], and prospective clinical profiling [@aung2018compass].

**Proposed (insert one sentence between this sentence and the existing one, before "Identifying the prognostically relevant programs..."):**

> ...yet the objective optimized during unsupervised discovery (reconstruction error) differs fundamentally from the criterion used during evaluation (survival association) [@cook2007fisher]. **Even within the unsupervised pipeline, model selection itself is unresolved: standard rank-selection diagnostics — reconstruction residuals, cophenetic correlation, silhouette indices — frequently disagree on the same dataset [@brunet2004metagenes; @frigyesi2008nmf], leaving practitioners without a principled criterion for choosing the number of programs to extract.** Identifying the prognostically relevant programs among the discovered factors requires extensive downstream evaluation: over a decade of PDAC subtyping was needed...

Word delta: ~+35 words. Re-adds a portion of the rank-selection pre-emption that was trimmed from P3 in the recent revision (so this is partially undoing that decision — flag for co-author discussion).

#### 5c. Combined effect if both 5a + 5b are applied

The intro establishes "rank-selection metrics are known to disagree" with two citations. The Results delivers the empirical instantiation: "we observed the documented phenomenon, see Fig. S4." Together they form a clean setup-and-payoff arc, which is rhetorically the strongest configuration. Total word cost across both: ~+53 words, still small relative to the recent intro trim of ~300 words.

### Recommendation

**Implement 5a unconditionally** — it is a small, surgical change with no bib churn, no intro-length impact, and clear defensive value against a methods reviewer.

**Treat 5b as opt-in based on co-author preference about the intro length budget.** If the recent P3 trim is something co-authors fought hard for, leave 5b out. If the intro can absorb +35 words, 5b materially strengthens the rhetorical setup for 5a.

### Dials to consider

- **Pan-method breadth.** Both already-in-bib citations are NMF-specific and cancer-genomics-flavored. If a reviewer pushes back that the criticism is too narrow (e.g., "cophenetic correlation is only one heuristic"), the Hutchins / Tan-Févotte / Stein-O'Brien citations listed above offer broader methodological coverage at the cost of 1–3 new bib entries.
- **Tone.** The drafts above frame the conflict as "a known limitation" — neutral, descriptive. If a co-author prefers a more pointed framing ("…motivating the supervised approach we develop here"), the §5b clause can absorb a connector that explicitly bridges to DeSurv. Risk: this veers from observation into thesis-statement and may feel premature in P2, where the contribution claim arrives in P3.
- **Placement of 5b.** The proposed insertion sits between the variance/prognosis sentence and the "extensive downstream evaluation" sentence. An alternative placement is at the very end of P2, as a standalone closer — this gives the rank-selection problem more rhetorical weight but also reads as a longer detour. Mid-P2 placement (as drafted) ties it more tightly to the variance/prognosis argument, which is preferable.

---

## 6. Cite Bayesian Optimization (BO) and explain briefly why it is used

### Why this matters

BO is named with an acronym in Methods (`paper/03_methods_REVISED.Rmd:23, 40, 48, 58`) and Results (`paper/04_results_REVISED.Rmd:194, 206, 208, 242`) but is **uncited anywhere in the manuscript or bibliography** (verified 2026-05-04: zero hits for `snoek` / `shahriari` / `frazier` / `jones1998` in `references_30102025.bib`, and zero hits for any of those in `paper/*.Rmd` and `paper/*.tex`). For a PNAS general-scientific audience this term is not common knowledge — the recent compliance pass added 38 DOIs to close exactly this kind of gap, so leaving BO uncited is inconsistent with the rest of the audit.

A second reason to fix this: BO's *suitability* for this problem is non-trivial and worth surfacing. Each cross-validated evaluation is an expensive black-box objective (full nested-CV refit of DeSurv), and the hyperparameter space is mixed integer/continuous ($k$, $n_{\text{top}}$ integer; $\alpha, \lambda, \xi$ continuous). Grid search is wasteful, random search is better but still inefficient, and BO's Gaussian-process surrogate plus acquisition-function policy is the textbook fit. One short clause covers the why and signals to readers unfamiliar with BO that this was a principled choice rather than a buzzword.

### Stale-comment cleanup (related)

The introduction file's keyword inventory at `paper/02_introduction_REVISED.Rmd:12` currently lists `snoek2012practical` as *"(in Methods/SI)"* — but that comment is stale; the citation is not actually placed anywhere in the manuscript. After this edit the comment will become accurate again. Worth a quick re-verification of the rest of that comment block during cleanup.

### Proposed edit (Methods, line 40 — first definition of BO)

**Current:**

> Hyperparameters $(k,\alpha,\lambda,\xi,n_{\text{top}})$ were selected by maximizing the cross-validated C-index using Bayesian optimization (BO), with final rank $k$ chosen by the one-standard-error rule...

**Proposed (citation + brief why-clause, ~+22 words):**

> Hyperparameters $(k,\alpha,\lambda,\xi,n_{\text{top}})$ were selected by maximizing the cross-validated C-index using Bayesian optimization (BO) **[@snoek2012practical; @shahriari2016taking], a sample-efficient global-search strategy for expensive black-box objectives over mixed integer–continuous hyperparameter spaces**, with final rank $k$ chosen by the one-standard-error rule...

The other in-text mentions of "Bayesian optimization (BO)" / "BO" can stay as-is once Methods is augmented — Methods is the right home for the citation and justification, and the Results occurrences function as cross-references back to it.

### Bib entries to add (two new entries)

- **`snoek2012practical`** — Snoek J, Larochelle H, Adams RP. "Practical Bayesian optimization of machine learning algorithms." *Advances in Neural Information Processing Systems* 25 (NeurIPS 2012). arXiv:1206.2944. URL: https://proceedings.neurips.cc/paper/2012/hash/05311655a15b75fab86956663e1819cd-Abstract.html
- **`shahriari2016taking`** — Shahriari B, Swersky K, Wang Z, Adams RP, de Freitas N. "Taking the human out of the loop: A review of Bayesian optimization." *Proceedings of the IEEE* 104(1):148–175 (2016). DOI: 10.1109/JPROC.2015.2494218.

### Cost/benefit and dials

- **Two new bib entries** — small relative to the 38 DOIs added in the recent compliance pass.
- **If only one citation can be used,** prefer **`snoek2012practical`** — it is the directly-applied "BO for ML hyperparameters" paper and is the citation the stale intro inventory was already trying to reference. `shahriari2016taking` is the conventional review-paper companion; including both gives "applied + comprehensive review" coverage that mirrors how methods sections in PNAS typically pair foundational and review citations.
- **If word budget is tight in Methods,** the why-clause can shrink to "[@snoek2012practical; @shahriari2016taking], a sample-efficient strategy for expensive black-box objectives" (drops "over mixed integer–continuous hyperparameter spaces" — saves ~6 words). The shortest defensible version is just the bare citation with no why-clause (~+5 words total), but that does not address the original concern about whether the term is common knowledge for the PNAS audience.

### Application checklist

1. Add the two bib entries to `paper/references_30102025.bib`.
2. Edit `paper/03_methods_REVISED.Rmd:40` per the proposed text above.
3. Update the keyword inventory comment at `paper/02_introduction_REVISED.Rmd:12` so the `snoek2012practical (in Methods/SI)` entry now reflects reality (or simply re-verify; the comment will become accurate after step 2).
4. Re-render with `make paper` and confirm the citations resolve cleanly in `paper.pdf`.

---

## 7. Clarify "top-weighted genes" on first use in Results, and unify terminology with abstract

### Update (2026-05-05): SI structure verified; apply with reviewer's wording

**SI verification is complete.** The SI table of contents confirms:

- **Section 6** = "Bayesian Optimization" (covers BO procedure and `ntop` selection)
- **Section 9** = "Gene Program Correlation Analysis" — contains the formal factor specificity score definition: $s_{ij} = W^*_{ij} - \max_{j' \neq j} W^*_{ij'}$ (verified verbatim at SI Appendix line 893–897)
- **Section 21** = "Factor-Specific Gene Lists" with Tables S7–S9 listing top 270 genes per factor

The current Methods pointer to "SI Appendix, Section 6" for top-gene truncation is technically defensible because Section 6 covers ntop selection, but the formal definition of "top-weighted" lives in Section 9. The right inline gloss separates the two distinct operations (which-genes vs how-many-genes):

**Revised §7a draft (supersedes the earlier proposed wording):**

> Here precision denotes the fraction of a factor's **top-weighted genes — those with the largest factor specificity score $s_{ij}$ (SI Appendix, Section 9), with $n_{\text{top}}$ chosen by Bayesian optimization (SI Appendix, Section 6)** — that overlap with a true prognostic program's gene set.

The dual-pointer structure (Section 9 for "which genes," Section 6 for "how many genes") matches the SI's actual structure exactly. Adopt SI Section 9's terminology verbatim — "**factor specificity score $s_{ij}$**" — so the main text and SI use the same term.

§7b (Fig 1B caption) and §7c (abstract terminology unification, two `top-ranked` → `top-weighted` swaps) stand as previously drafted; only §7a's wording changes.

### Why this matters

The Results body sentence at `paper/04_results_REVISED.Rmd:244` defines _precision_ in terms of "top-weighted genes" without telling the reader what "top-weighted" means or where to find a formal definition:

> "Here precision denotes the fraction of top-weighted genes in a learned factor that overlap with a true prognostic program's gene set."

Three problems stack:

1. **"Top-weighted" is used as a primitive.** A reader hitting this term has to either (a) guess from context — likely landing on "highest values in column of $W$" — or (b) flip to Methods + SI to find the formal definition. The Methods sentence at line 40 mentions $n_{\text{top}}$ in the PDAC validation context but does not define "top-weighted" for the simulation-precision context, and the formal procedure lives in SI Appendix, Section 6 with no inline pointer from Results.

2. **The naive reading of "top-weighted" is wrong.** Inspection of `R/get_top_genes.R` (verified 2026-05-04) shows the actual procedure is **factor-specific differential ranking after column-max normalization**, not raw column ranking. Concretely: each column of $W$ is normalized to its column maximum (lines 6–10); for each factor $i$, `diff_vector = current_col - max_other` is computed with `max_other` the elementwise max across all _other_ factor columns (lines 32–35); the top $n_{\text{top}}$ genes are those with the largest values of this difference (line 36). The result is genes that are most distinctive to factor $i$ relative to all other factors — a markedly different set from "highest absolute $W$-loadings". This was a deliberate design choice (it makes top-gene sets factor-specific markers rather than bulk high-loaders), but the term gives readers no signal that this is happening.

3. **Terminology mismatch between abstract and Results.** The abstract (`paper/paper.Rmd:60`) uses **"top-ranked"** for the same concept (twice: "top-ranked genes of each factor" and "top-ranked factor genes"), while Results uses **"top-weighted"**. Cancer-genomics readers are sensitive to terminology drift between abstract and body — these should be the same word. "Top-weighted" is the more accurate label (the metric is loading-based) and the more common term in current Results, so the unification should standardize on **top-weighted**.

### Proposed edits

#### 7a. Inline gloss + pointer in Results (preferred)

**Current (`paper/04_results_REVISED.Rmd:244`):**

> Here precision denotes the fraction of top-weighted genes in a learned factor that overlap with a true prognostic program's gene set.

**Proposed (~+15 words):**

> Here precision denotes the fraction of a factor's **top-weighted genes — those with the largest factor-specific loadings in $W$ after column-max normalization (Materials and Methods; SI Appendix, Section 6)** — that overlap with a true prognostic program's gene set.

This names the matrix ($W$), names the operation (factor-specific, normalized), and points to the formal definition. The phrase "factor-specific" is the load-bearing word — it signals that the ranking is differential vs. other factors, not absolute, which is the non-obvious part.

#### 7b. Update the figure caption to match

**Current (`paper/04_results_REVISED.Rmd:208`, Fig 1B caption):**

> (B) Precision (fraction of top-weighted genes in a learned factor overlapping a true prognostic program) across replicates.

**Proposed:**

> (B) Precision (fraction of a factor's top-weighted genes — defined in main text — overlapping a true prognostic program) across replicates.

Minimal change; avoids re-stating the full definition in the caption while pointing readers who skim figures back to the body.

#### 7c. Standardize terminology in abstract

**Current (`paper/paper.Rmd:60`, two occurrences):**

> ...DeSurv places more true prognostic genes among the **top-ranked** genes of each factor than standard NMF (precision, defined as the fraction of **top-ranked** factor genes that are true prognostic markers, 0.50 vs. 0.07)...

**Proposed (replace both occurrences):**

> ...DeSurv places more true prognostic genes among the **top-weighted** genes of each factor than standard NMF (precision, defined as the fraction of **top-weighted** factor genes that are true prognostic markers, 0.50 vs. 0.07)...

Two-word swap. No net word count change. Aligns abstract with the Results terminology and the SI Section 6 formal definition.

#### 7d. Optional: alternative drafts for the Results sentence

**Pointer only (~+8 words, lighter):**

> Here precision denotes the fraction of top-weighted genes in a learned factor that overlap with a true prognostic program's gene set **(top-weighted genes defined in Materials and Methods; SI Appendix, Section 6)**.

Resolves the "where do I find a definition" half of the gap but not the "what does it mean" half. Choose this if the inline gloss in 7a feels heavy.

**Full inline definition (~+30 words, most rigorous):**

> Here precision denotes the fraction of a factor's **top-weighted genes — defined as those whose column-max-normalized loading in that factor exceeds their loading in every other factor by the largest margin (Materials and Methods; SI Appendix, Section 6) — that overlap with a true prognostic program's gene set**.

Use this only if SI Section 6 itself is terse and you want self-contained Results. For most readers 7a is the better balance.

### Verify the SI definition before committing

Before applying 7a/7b/7c, **read SI Appendix Section 6** to confirm two things:

1. The formal definition there matches the inline gloss ("largest factor-specific loadings in $W$ after column-max normalization"). If SI Section 6 phrases it differently (e.g., uses "discriminative loading" or "factor-distinctive" language), align the inline gloss to that wording so the two locations match exactly.
2. SI Section 6 is the right pointer. If the definition actually lives in a different SI section, update the pointer.

The wording in 7a above was inferred from `R/get_top_genes.R`, not from SI Section 6 directly. The code is authoritative for what the procedure does, but the manuscript-side terminology should match SI Section 6 for consistency.

### Cost/benefit

- **Word delta:** ~+15 in Results body, ~+5 in figure caption, 0 net in abstract. Total ~+20 words across three locations. Well within budget.
- **No bib changes.**
- **Reader benefit:** removes a concrete clarity gap on a term that's central to the simulation Results (it appears in Fig 1B, the body precision claim, and the abstract's headline 0.50 vs. 0.07 number). A reviewer asking "what does top-weighted mean?" should not happen after this edit.
- **Risk:** very low. The unified terminology and inline gloss don't change any results or claims; they only make existing definitions self-locating.

### Application checklist

1. Read SI Appendix Section 6 and confirm the formal definition matches the proposed inline gloss in 7a (adjust wording if needed).
2. Apply 7a (Results body sentence at line 244) and 7b (Fig 1B caption at line 208) in `paper/04_results_REVISED.Rmd`.
3. Apply 7c (abstract terminology unification, two `top-ranked` → `top-weighted` swaps) in `paper/paper.Rmd:60`.
4. Search the rest of the manuscript and SI for any remaining occurrences of "top-ranked" referring to this concept (`grep -nE "top.ranked|top-ranked|top weighted|top-weighted" paper/*.Rmd paper/si_appendix.Rmd`) and unify if any are found.
5. Re-render with `make paper` and visually confirm the abstract still reads cleanly with the swap.

---

## 25. Disambiguate or rethink the events-per-variable (EPV) sentence in Methods

### Update (2026-05-05): two layers of question — wording and necessity

The original §25 entry below addresses the **wording** issue Amber flagged (the bare "$\div 3$" is ambiguous about whether 3 is the count or a value). A subsequent question raised whether the EPV sentence is **necessary at all**, given that DeSurv's $\beta$ is jointly estimated with $W$ under elastic-net regularization rather than fit by standard MLE Cox regression. The 10-EPV rule (Peduzzi 1996; Vittinghoff & McCulloch 2007) was derived for standard MLE Cox regression where $\beta$ is the only parameter being estimated.

This update reframes §25 around three paths, of which the original entry below now serves as **Path A**.

### Three paths

#### Path A — Keep and clarify (recommended, original §25 entry's territory)

Add a "Cox-stage" qualifier so the claim's scope is explicit. Best draft (refines Option 1 from the original entry below):

> "With 139 events and 3 factor-score covariates (one per DeSurv factor at $k = 3$), the **Cox-stage events-per-variable ratio** is 46, well above the conventional 10-EPV threshold."

The "Cox-stage" qualifier scopes the claim to where the 10-EPV rule actually applies, pre-empting the W-parameters follow-up question without committing to a longer methodological exposition. ~+5 words on the original. Keeps a TRIPOD-checkbox-compatible EPV statement.

#### Path B — Drop the sentence entirely

The participant-flow paragraph (line 53) already conveys 273 samples / 139 events. A reader who cares about EPV can compute 46 from those numbers themselves. Removing the sentence:

- **Saves words** (~−25).
- **Avoids the W-parameters follow-up question** entirely — no claim about estimation stability is made.
- **Loses the TRIPOD-checkbox compliance** if the lab's TRIPOD checklist requires EPV reporting (verify against `docs/tripod_checklist.md`).

#### Path C — Replace with a substantive estimation-stability note

> "Estimation stability was supported by elastic-net regularization on $\beta$ (BO-selected $\lambda$ and $\xi$) and cross-validated hyperparameter selection (see Methods, BO procedure)."

Replaces the EPV statement with a description of what *actually* stabilizes the estimates in DeSurv. More methodologically accurate but substitutes a non-conventional framing for the conventional one — which means it doesn't satisfy a TRIPOD-style box if there is one. ~+15 words net (drops EPV sentence, adds the substantive note).

### Decision criteria

| Condition | Recommended path |
|---|---|
| Lab's TRIPOD checklist requires EPV reporting | **Path A** (clarify, satisfy the box) |
| TRIPOD checklist does not specifically require EPV | **Path A** (default) or **Path B** (if word budget tight) |
| Co-author or reviewer has flagged "rule applies?" concern | **Path C** (substantive replacement) |

Path A is the default unless one of the other conditions is met. Verify against `docs/tripod_checklist.md` before deciding between A and B.

### Why Path A is the default

1. **TRIPOD/methods compliance** routinely includes EPV reporting; dropping (Path B) creates a checkbox gap.
2. **The "Cox-stage" qualifier is a small addition** that pre-empts the W-parameters follow-up without committing to longer exposition.
3. **The numbers genuinely pass with margin** (46 vs 10) — there's no reason to hide a comfortable result.
4. **Path C is more accurate but heavier** and substitutes the conventional EPV framing for a non-conventional one — likely fails a TRIPOD checkbox.

### Original §25 content (Path A wording variants)

The wording variants below are all Path A variants — they keep the EPV statement but resolve the ambiguity Amber flagged. Path A's recommended draft (with "Cox-stage" qualifier) is captured above; the variants below offer additional wording options if a different formulation is preferred.

### Why this matters

The Methods sentence at `paper/03_methods_REVISED.Rmd:53` reads:

> "The events-per-coefficient ratio in the training Cox model ($139 \div 3 = 46$) is well above the conventional 10-EPV threshold, supporting estimation stability."

Amber flagged that the bare "$\div 3$" is ambiguous: a reader can parse it as "the coefficient is 3" rather than "there are 3 coefficients." The "3" is actually the **number of fitted Cox coefficients**, equivalent to $k = 3$ DeSurv factors — but the sentence doesn't make that link explicit.

Two issues compound the confusion:

1. **The "3" is unexplained.** Without "$k = 3$ factor coefficients" or similar disambiguation, the division by 3 looks unmotivated.
2. **Terminology mismatch within the sentence.** "Events-per-coefficient ratio" → "10-EPV threshold." EPV = events-per-variable (Peduzzi 1996 is the canonical citation). The sentence uses the non-standard "events-per-coefficient" then references the standard "10-EPV" — minor inconsistency that adds friction.

### Three rewordings

**Option 1 — Spell out the connection to $k = 3$ (preferred):**

> "With 139 events and **3 factor-score covariates (one per DeSurv factor at $k = 3$)**, the events-per-variable ratio of 46 is well above the conventional 10-EPV threshold, supporting estimation stability."

Restructures the sentence: leads with inputs (139 events, 3 covariates with explicit $k = 3$ link), then states the ratio. Removes the "$\div 3$" ambiguity by naming what's being divided. Switches to "events-per-variable" for terminology consistency with "10-EPV."

**Option 2 — Minimal change, keep inline arithmetic:**

> "The events-per-variable ratio in the training Cox model — 139 events ÷ **3 factor coefficients ($k = 3$)** = 46 — is well above the conventional 10-EPV threshold, supporting estimation stability."

Keeps the existing arithmetic structure but disambiguates "3" with "factor coefficients ($k = 3$)" and switches "events-per-coefficient" → "events-per-variable" for consistency.

**Option 3 — Most pedagogical (~+15 words):**

> "The trained Cox model fits one coefficient per DeSurv factor (3 coefficients at $k = 3$). With 139 events in the training cohort, the events-per-variable ratio is $139/3 \approx 46$, well above the conventional 10-EPV threshold, supporting estimation stability."

Names what's being counted, the ratio, and the rule. Heavier but unambiguous. Optional citation to Peduzzi 1996 (`peduzzi1996simulation`) — not strictly required if the rule is treated as common knowledge.

### Recommendation: Option 1

Clearest reading for someone parsing "3" as potentially ambiguous (Amber's catch). Most defensible against a reviewer reading: "factor-score covariates (one per DeSurv factor at $k = 3$)" makes the methodological claim airtight. Option 2 is fine if you want to preserve the inline arithmetic; Option 3 is heaviest and only needed if a reviewer specifically asks for the EPV rule citation.

### Methodological subtlety (probably doesn't need to surface in the manuscript)

The 10-EPV rule was developed for standard MLE Cox regression where the only learned parameters are the regression coefficients. DeSurv's joint optimization also learns the gene-program matrix $W$ via the supervision gradient, which has many more entries than 3. A strict-reading reviewer could ask "what about the effective parameters in $W$?"

The defensible answer:

- EPV applies to the Cox stage. The prognostic prediction comes from the 3 fitted coefficients applied to projected scores $W^\top X$.
- The W-matrix updates are constrained by reconstruction loss and elastic-net regularization on $\beta$, not by free-parameter MLE estimation.
- At inference time, the projection-based scoring uses the same 3 coefficients applied to new samples — the relevant statistical inference target.

This subtlety probably doesn't need to be addressed in the Methods unless a reviewer raises it. Option 1's wording stays inside the standard Cox-stage interpretation, which is conventional usage.

### Cost/benefit

- **Word delta:** ~+5 (Options 1–2) or +15 (Option 3).
- **No bib changes** unless adding Peduzzi 1996 in Option 3 (and only if not already in bib — needs verification).
- **Reader benefit:** removes the "is the coefficient 3?" misreading that Amber flagged; aligns terminology with the standard EPV rule.
- **Risk:** zero. Methodological content unchanged; only labeling/wording shifts.

### Application checklist

1. Apply Option 1, 2, or 3 to `paper/03_methods_REVISED.Rmd:53`, replacing the existing EPV sentence.
2. If using Option 3 with citation, verify whether `peduzzi1996simulation` is already in the bib; add if needed.
3. Re-render with `make paper` and confirm the Methods paragraph reads cleanly.

---

## 26. Optionally clarify the 1-SE rule's reduction from multi-dimensional BO surface to a 1-D rank choice

### Why this matters

The Methods sentence at `paper/03_methods_REVISED.Rmd:40` reads:

> "...with final rank $k$ chosen by the one-standard-error rule (the smallest $k$ whose predicted performance lay within one standard error of the maximum)."

The phrase **"predicted performance"** at a given $k$ is conventional BO/CV shorthand for "the best mean cross-validated C-index achievable at that $k$ over the remaining hyperparameters ($\alpha, \lambda, \xi, n_\text{top}$)." A reader familiar with BO + 1-SE will infer this, but the text doesn't make it explicit.

Code-verified algorithm (in `R/bo_helpers.R`, `select_bo_k_by_cv_se` lines 13–150 and `select_bo_params_1se` lines 152–208):

1. Group the BO history by $k$. At each unique $k$, take the single best configuration (highest `mean_cindex`, tie-break: smallest $\alpha$). → `best_per_k` (one row per $k$).
2. Find the global best across `best_per_k` (highest `mean_cindex`, tie-break: smallest $k$).
3. Compute threshold = `best_mean − best_se`.
4. Filter `best_per_k` to rows with `mean_cindex ≥ threshold`.
5. Select $k$ = smallest such $k$.
6. At $k = k_\text{selected}$, pick the row with the highest `mean_cindex` and read off the full hyperparameter vector from that row.

Key nuance: **"within 1 SE" is evaluated on per-$k$ bests, not on all BO evaluations.** Suppose at $k = 3$ you have two configs (mean 0.65, 0.62) and at $k = 5$ one config (0.66, the global max). Threshold = 0.66 − SE. The code asks "is $k = 3$'s *best* (0.65) within threshold?" — not "are any of $k = 3$'s configs within threshold?" In practice these usually agree; formally the granularity is per-$k$-best.

After selecting $k$, the chosen hyperparameter vector is the best config at that $k$, which is the within-1-SE candidate from Step 4–5.

### Two options

**Option 1 — Leave as-is (no edit, default for most contexts).** The current text is accurate at a high level; "predicted performance at rank $k$" is conventional shorthand. A BO-aware reader will infer the per-$k$-best aggregation. Most readers don't need the explicit elaboration.

**Option 2 — Add a brief disambiguating clause (~+25 words):**

> "...with final rank $k$ chosen by the one-standard-error rule: the smallest $k$ whose **per-rank best cross-validated C-index (maximum over $\alpha$, $\lambda$, $\xi$, $n_\text{top}$ at that $k$)** lay within one standard error of the global best. The remaining hyperparameters were then taken from the configuration achieving that per-rank best."

Names the per-$k$-best aggregation explicitly and clarifies how the full hyperparameter vector is then recovered. Removes the implicit ambiguity at the cost of length.

### Recommendation: Option 1

The current text is technically correct and conventional in the BO literature. Option 2 is worth using only if a co-author or reviewer specifically asks how the multi-dimensional BO surface is reduced to the 1-D rank-selection rule. For a PNAS general-audience reader, the brief shorthand is sufficient — the detailed Methods description in the SI (Section 6, BO procedure) and the code-level definitive answer cover the rigorous case.

### Cost/benefit

- **Word delta:** 0 (Option 1) or +25 (Option 2).
- **No bib changes, no SI changes.**
- **Reader benefit:** Option 2 removes a rare-edge-case ambiguity. Option 1 keeps the Methods compact.
- **Risk:** zero.

### Application checklist

1. Decide between Option 1 (leave) or Option 2 (apply) based on co-author or reviewer feedback.
2. If applying Option 2, edit `paper/03_methods_REVISED.Rmd:40` to insert the disambiguating clause.
3. Re-render with `make paper` and confirm the Methods paragraph reads cleanly.

---

## 27. Correct the inaccurate "converted to TPM" data-provenance claim in Methods and SI

### Update (2026-05-05, post-Amber-input): drop long version; use minimal correction

Amber's response: *"I think this is more detail than we need. Whether we did it ourselves or not, at the end of the day all of the rnaseq cohorts were analyzed on log tpm scale. If you want more details, Laura is the person to ask."*

**Recommended approach (preferred): minimal correction at both Methods and SI.** Drop the per-cohort breakdown (TCGA-PAAD, CPTAC, Dijk, etc.) entirely. The factual error is the formula `log₂(TPM + 1)` applied universally — TPM doesn't apply to microarray cohorts. The minimum-change correction:

**Methods (line 51) — minimal correction:**

> "Gene expression matrices were analyzed on a log scale: RNA-seq cohorts as log₂(TPM + 1); microarray cohorts in platform-normalized log-scale intensities."

**SI Section 7 (line 822) — minimal correction (same wording):**

> "Expression matrices were analyzed on a log scale: RNA-seq cohorts as log₂(TPM + 1); microarray cohorts in platform-normalized log-scale intensities."

Net delta: net 0 words at Methods (replaces existing sentence with same length); net 0 words at SI. Removes the factual error without overspecifying.

**If reviewers ask for per-cohort detail**, defer to Xianlu Laura Peng (DeCAF first author) for the most precise per-cohort scaling description — she's closest to the DeCAF compilation and would know exact provenance.

The original 27a long version (~+45 words) and compact alternative below remain documented as wording references. **The minimal correction above supersedes both as the recommended apply-version.**

### Why this matters

The Methods sentence at `paper/03_methods_REVISED.Rmd:51` and SI Section 7 sentence at `paper/si_appendix.Rmd:822` both claim:

> "Gene expression matrices were converted to transcripts per million (TPM) and log₂-transformed (log₂(TPM + 1))." [Methods]
> "All expression matrices were log₂-transformed (log₂(TPM + 1)) prior to gene selection." [SI]

**Both claims are factually wrong** on multiple counts:

1. **No TPM conversion is performed by DeSurv.** Code inspection (`R/load_data_internal.R`) confirms expression matrices are loaded from `data/original/<cohort>.rds` and used in whatever format the source provided. DeSurv does not compute TPM from raw counts.
2. **`log₂(x + 1)` is applied only to a subset of cohorts.** Specifically: TCGA-PAAD (line 21), Dijk (line 13), and PACA-AU RNA-seq (line 64). It is **not** applied to CPTAC, Moffitt, PACA-AU array, or Puleo — those cohorts arrive pre-transformed and are used as-is.
3. **TPM is RNA-seq-specific.** Microarray cohorts (Moffitt, PACA-AU array, Puleo) cannot be in TPM by definition — they're normalized intensities (RMA or platform-equivalent). Saying "log₂(TPM + 1)" was applied to microarray data is a category error.

A methods-aware reviewer who knows TPM is an RNA-seq quantification will flag this immediately. The fix is small but needed in two locations.

### Actual data provenance (verified from `R/load_data_internal.R` and DeCAF Methods)

| Cohort | Platform | Source format on disk | Scaling applied by DeSurv |
|---|---|---|---|
| TCGA-PAAD | RNA-seq | Linear TPM (GDC) | `log₂(x + 1)` |
| CPTAC | RNA-seq | Pre-log₂ (GDC) | Used as-is |
| Dijk | RNA-seq | Linear (DeCAF compilation) | `log₂(x + 1)` |
| PACA-AU RNA-seq | RNA-seq | Linear (DeCAF compilation) | `log₂(x + 1)` |
| Moffitt | Microarray | RMA-normalized log-scale | Used as-is |
| PACA-AU array | Microarray | Platform-normalized log-scale | Used as-is + gene-symbol collapse |
| Puleo | Microarray | RMA-normalized log-scale | Used as-is |

Validation data was compiled by the DeCAF group [@peng2024determination], which applied log₂ transformation + column-wise quantile normalization across cohorts (per DeCAF Methods). DeSurv inherited the DeCAF compilation for validation data and obtained training data directly from GDC.

After scaling, all cohorts pass through DeSurv's **within-subject rank transformation** to harmonize across platforms — that is the actual cross-platform alignment step.

### Proposed corrections

#### 27a. Methods (line 51)

**Current:**

> "We analyzed publicly available RNA sequencing (RNA-seq) and microarray cohorts of PDAC with corresponding overall survival outcomes. **Gene expression matrices were converted to transcripts per million (TPM) and log₂-transformed (log₂(TPM + 1)).** For each training cohort separately..."

**Proposed (long version, ~+45 words):**

> "We analyzed publicly available RNA sequencing (RNA-seq) and microarray cohorts of PDAC with corresponding overall survival outcomes. **Expression matrices were used in the form provided by each source: RNA-seq cohorts contributed TPM-quantified or pre-log-transformed expression (TCGA-PAAD and PACA-AU RNA-seq in linear-scale TPM; CPTAC pre-log-transformed by GDC; Dijk in linear scale from the DeCAF compilation [@peng2024determination]); microarray cohorts (Moffitt, PACA-AU array, Puleo) contributed RMA- or platform-normalized log-scale intensities. Linear-scale RNA-seq cohorts were placed on a comparable log-scale via $\log_2(x + 1)$; pre-log-transformed sources were used as-is.** For each training cohort separately..."

**Compact alternative (~+15 words):**

> "...Expression matrices were used in the form provided by each source — RNA-seq cohorts in TPM or pre-log-transformed scale, microarray cohorts in platform-normalized log-scale intensities — with linear-scale RNA-seq data log₂-transformed (log₂(x + 1)) for comparability. For each training cohort separately..."

Use compact version if word budget is tight. Long version is more transparent for reviewers asking about cross-platform handling.

#### 27b. SI Section 7 (line 822)

**Current:**

> "All expression matrices were log₂-transformed (log₂(TPM + 1)) prior to gene selection."

**Proposed:**

> "Expression matrices were placed on a comparable log-scale prior to gene selection: linear-scale RNA-seq cohorts (TCGA-PAAD, Dijk, PACA-AU RNA-seq) were log₂-transformed (log₂(x + 1)); CPTAC was used as-is (pre-log-transformed by GDC); microarray cohorts (Moffitt, PACA-AU array, Puleo) were used as published in their RMA-normalized log-scale form."

Net delta: ~+30 words. Cohort-by-cohort accuracy at the SI level.

### Coordination with the existing dataset-specific paragraph

SI lines 810–818 already describe dataset-specific inclusion and gene-symbol handling. The §27 fix doesn't conflict with that paragraph; it corrects the high-level statement that follows. After applying §27b, the two paragraphs together (lines 810–818 + corrected line 822) give a complete, accurate provenance story. The SI lines 837–840 ("Gene expression data and molecular classifications for the validation cohorts were originally compiled for the DeCAF study") already attribute the compilation correctly and reinforce 27a's `peng2024determination` reference.

### Cost/benefit

- **Word delta:** +15 to +45 in Methods (depending on version), +30 in SI.
- **No bib changes.** `peng2024determination` already cited.
- **Reader benefit:** removes a factual error that any methods reviewer who knows TPM ≠ microarray will catch immediately. Aligns the manuscript's data-provenance description with what the code actually does.
- **Risk:** zero — substantive content unchanged; only the description of preprocessing is corrected to match reality.

### Application checklist

1. Apply 27a (Methods line 51) using either the long or compact version per word-budget preference.
2. Apply 27b (SI Section 7 line 822).
3. Verify against `R/load_data_internal.R` that the per-cohort scaling claim matches exactly (PACA-AU array gene-symbol collapse step is already documented in SI lines 814–817 as it should be).
4. Re-render with `make paper` and confirm both the Methods paragraph and SI Section 7 read cleanly.

### What this fix does NOT change

- The within-subject rank transformation step is the actual cross-platform harmonization mechanism. That description in the existing Methods paragraph stands as-is and should be retained verbatim.
- Gene-set selection (top 3,000 by mean+variance ranks per training cohort, intersected to 1,970 genes) is unchanged.
- Validation gene-set restriction to the 1,970 genes is unchanged.
- The DeCAF compilation attribution (already in SI Section 7 lines 837–840) stands and reinforces 27a's reference.

---

# Cluster F — Discussion structural cleanup

Four Discussion-section edits: theoretical reframing (line 6), redundancy removal (lines 10 and 14), closing paragraph strengthening (line 16). Coordinate carefully — §21 depends on §19 application (which moves plateau/NMF-needs-BO content earlier, making line 10 paragraph more clearly redundant); §24 depends on §17 + §20 (which establish framing the closing should mirror).

**§24 functions as the fourth member of Cluster A's "so what" arc** (closing rhetorical hook for Discussion exit); filed here because of its primary Discussion-cleanup nature.

**Additional candidate site identified during cluster-expansion analysis:** Discussion line 12 (deployment guidance + Table S6 reference) could be lightly enhanced to mention concordance + independence together — currently emphasizes only independence. Coordinate with Cluster G (§28).

---

## 17. Reframe Discussion's connection to sufficient dimension reduction / information bottleneck theory to foreground DeSurv's contribution rather than make it sound like translation

### Why this matters

The opening paragraph of the Discussion at `paper/05_discussion_REVISED.Rmd:6` connects DeSurv to foundational theory:

> "More broadly, the principle that outcome-guided dimensionality reduction targets different subspaces than variance-driven reduction extends beyond cancer genomics: sufficient dimension reduction theory [@cook2007fisher] and the information bottleneck framework [@tishby1999information] both **predict** that supervised compression retains outcome-relevant structure while discarding noise variation. DeSurv **implements this principle within the specific constraints of NMF deconvolution**, where nonnegativity preserves biological interpretability and the factorization structure enables single-sample scoring by projection."

Two related issues with the current framing:

1. **The verb "predict" implies prior theory anticipated DeSurv's specific finding.** A reader who parses "predicts" reads the sentence as "this idea isn't new — these theories already said this." That subtly lowers the perceived novelty of the contribution.
2. **"Implements this principle within the specific constraints" reads passively** — like translation rather than innovation. Sufficient dimension reduction theory was developed for unconstrained, continuous-outcome regression. The information bottleneck is general classification theory. **Neither directly addresses the four constraints DeSurv resolves**: (a) NMF's nonnegativity, (b) censored survival, (c) mixture-coefficient interpretability of $H$, and (d) single-sample transportability via projection. Saying DeSurv "implements" the principle obscures that the empirical extension to this constrained setting is the paper's actual contribution.

The intro at `paper/02_introduction_REVISED.Rmd:19` already frames this correctly:

> "Whether this advantage extends to NMF-based deconvolution, where nonnegative constraints restrict the factor space geometry, censored survival outcomes provide a noisier supervisory signal, and overly aggressive survival gradients could compromise the mixture-coefficient interpretation of sample loadings, remains an open empirical question."

The Discussion can simply close that loop — the open empirical question framed in the intro is resolved here in the affirmative. That framing makes prior theory the *motivation* for DeSurv's empirical question, not a *prediction* that pre-empts the contribution.

### Three options

#### Option 1 — Light touch (~+20 words)

> "...sufficient dimension reduction theory [@cook2007fisher] and the information bottleneck framework [@tishby1999information] **provide formal grounding for the same intuition. DeSurv extends this principle to a setting prior theory does not directly address** — bulk-tumor deconvolution under nonnegativity, censored survival, and mixture-coefficient interpretability constraints — **demonstrating empirically that supervised factorization can recover outcome-relevant structure where the standard discover-then-evaluate pipeline cannot.** The architecture (supervision on $W$, mixture interpretation preserved on $H$) enables single-sample scoring by projection, a property not shared by methods that route supervision through $H$ [@huang2020low; @le2025survnmf]."

Replaces "predict" with "provide formal grounding for the same intuition" (theory motivates rather than pre-empts). Replaces "implements this principle" with "extends this principle to a setting prior theory does not directly address." Adds contrast with $H$-routed methods.

#### Option 2 — Reference back to the intro's open empirical question (~+30 words, **preferred**)

> "...sufficient dimension reduction theory [@cook2007fisher] and the information bottleneck framework [@tishby1999information] both formalize this intuition. **The empirical question framed in the Introduction — whether this principle survives the constraints of nonnegativity, censored survival outcomes, and mixture-coefficient interpretability that bulk-tumor deconvolution imposes — is resolved here in the affirmative.** DeSurv's architecture (supervision on $W$, mixture interpretation preserved on $H$) enables single-sample scoring by projection, a transportability property not shared by methods that route supervision through $H$ [@huang2020low; @le2025survnmf]."

Closes the empirical-question loop the intro opens. Tightest version if intro and Discussion are meant to frame contribution consistently across the paper.

#### Option 3 — Aggressive reframing (~+40 words)

> "...sufficient dimension reduction theory [@cook2007fisher] and the information bottleneck framework [@tishby1999information] **formalized this intuition decades ago for unconstrained, continuous-outcome settings. The empirical question this paper resolves is whether the principle survives the specific constraints of bulk-tumor deconvolution**: NMF's nonnegative cone restricts the factor space geometry, censored survival provides a noisier supervisory signal than the continuous outcomes traditionally studied, and the mixture-coefficient interpretation of sample loadings can be compromised by overly aggressive survival gradients. DeSurv's architecture (supervision on $W$, mixture preserved on $H$) demonstrates that the principle does survive — yielding factors that are transportable across independent cohorts without retraining, a property not shared by methods that route supervision through $H$ [@huang2020low; @le2025survnmf]."

Most explicit about what prior theory doesn't cover and what DeSurv resolves. Heaviest; use only if a co-author or reviewer has specifically flagged this concern.

### Recommendation: Option 2

Three reasons:

1. **It closes the empirical-question loop the intro explicitly opens.** The intro frames DeSurv's value proposition as resolving an open empirical question (does the supervised-dim-reduction principle extend to NMF deconvolution?). The Discussion currently doesn't reference this framing back. Option 2 does, making the paper's framing arc consistent (intro raises question; Discussion confirms answer).
2. **It names the four constraints** (nonnegativity, censoring, mixture interpretability, transportability) that make extension non-trivial. Prior theory doesn't address these jointly; DeSurv's contribution is the joint resolution.
3. **It adds the architectural distinction** ($W$-routing vs. $H$-routing) with citations (`huang2020low`, `le2025survnmf`) already in the bib. This concretely contrasts DeSurv with the closest related methods.

Option 1 is acceptable if Option 2's intro-callback feels too heavy for a Discussion opening. Option 3 is overkill unless the concern has been independently flagged.

### Why this doesn't actually weaken the paper

The user's underlying worry is whether citing prior theory lessens the paper's impact. Counter-argument:

- **Methods without theoretical lineage read as ad-hoc heuristics.** Methods that extend foundational theory read as principled. PNAS reviewers in particular value the latter framing — it signals that the authors understand the broader scientific context.
- **The risk isn't the citations; it's the verb.** "Predicts" is the load-bearing problem word. A reframed sentence using "motivates" / "formalizes the same intuition" preserves the credit to Cook 2007 and Tishby 1999 while foregrounding what the paper actually demonstrates.
- **The contribution claim is empirically substantive even with the theoretical citations.** What DeSurv shows: the principle *does* survive the four constraints, and the resulting factors generalize across five independent cohorts (HR per SD = 1.50). Prior theory does not show this; it predicts only that supervision will help in the unconstrained continuous-outcome setting. The empirical resolution at the constrained censored-survival NMF setting is the paper's real contribution, and Option 2 makes that explicit.

### Cost/benefit

- **Word delta:** ~+20 (Option 1), +30 (Option 2), or +40 (Option 3).
- **No bib changes.** All citations (`cook2007fisher`, `tishby1999information`, `huang2020low`, `le2025survnmf`) already in bib.
- **Reader benefit:** sharpens contribution claim, prevents the "implementation of an existing idea" misreading, closes the intro-discussion framing loop.
- **Risk:** very low. Substantive content unchanged; only the framing of contribution-vs-prior-theory shifts.

### Application checklist

1. Choose Option 1, 2, or 3 based on co-author preference and word budget tolerance.
2. Apply the chosen draft to `paper/05_discussion_REVISED.Rmd:6`, replacing the existing two sentences ("...sufficient dimension reduction theory..." through "...single-sample scoring by projection.").
3. Verify the citations resolve correctly (no new bib entries; only existing keys reused).
4. Re-render with `make paper` and confirm the Discussion opening reads cleanly with the revised framing.

---

## 21. Remove redundant "fewer factors" framing in the Discussion paragraph at line 10

### Update (2026-05-05): paragraph-level preview adds two more streamlining options for the same paragraph

After §19a moves the plateau range to Results line 206 and §19b/§20a moves the "NMF needs BO" content to Results line 355, the Discussion line 10 paragraph compresses substantially (from ~270 words / 8 sentences to ~150 words / 5 sentences). Paragraph-level preview revealed two additional streamlining opportunities beyond the original §21 "fewer factors" fix:

**Streamlining add-on 1 — Merge sentences 7 and 8** (the model-selection-framework-as-contribution claim and the 1-SE parsimony explanation are both about the 1-SE rule's role and can be merged):

> "The model selection framework — joint BO over rank and hyperparameters with the 1-SE rule — is itself part of the methodological contribution, not merely the survival gradient: parsimony is a property of the 1-SE rule, but supervision creates the flat concordance surface under which the rule can operate effectively."

Net delta: ~−15 words, eliminates conceptual repetition between sentences 7 and 8.

**Streamlining add-on 2 — Compress sentences 1–2** (after §19a moves the plateau range out, the introductory framing is disproportionately weighty relative to the remaining empirical content):

> "Survival supervision concentrates prognostic signal into fewer factors, creating a broad concordance plateau that the one-standard-error rule leverages to select parsimonious ranks — a deliberate balance between biological completeness and clinical relevance."

Net delta: ~−15 words. Sentence 1's "balance between biological completeness and clinical relevance" becomes a closing clause rather than a standalone topic sentence.

**Combined effect of all three §21 streamlining steps** (original "fewer factors" fix + add-on 1 + add-on 2): paragraph compresses from ~150 words (post-§19/§20) to ~110 words. That feels right — paragraph stays substantive but stops feeling like a half-empty container.

**Coordination with §20b/§20c:** if §20b/§20c lands earlier in the paper (intro/abstract two-part contribution framing), the merged sentence 7-8 above can be reframed as a callback: *"As foreshadowed in the Introduction, the model selection framework — joint BO with the 1-SE rule — is itself part of the methodological contribution: parsimony is a property of the 1-SE rule, but supervision creates the flat concordance surface under which it can operate effectively."* This is the same content; the "As foreshadowed in the Introduction" prefix explicitly bridges. Apply §21 add-on 1 + §20b/§20c callback prefix together for full coordination.

### Why this matters

The Discussion paragraph at `paper/05_discussion_REVISED.Rmd:10` makes the parsimony argument twice within the same paragraph:

- **Sentence 4** (already in paragraph): *"...for standard NMF, concordance increased steadily with $k$, approaching DeSurv's level by $k = 7$ (SI Appendix, Table S5) but **requiring more than twice as many factors** and producing a fragmented structure."*
- **Sentence 7** (the redundant restatement): *"Even so, DeSurv at $k = 3$ achieves independent prognostic value **with fewer factors**: both methods retain significance after classifier adjustment (SI Appendix, Table S6), and $k = 3$ is more robust across supervision strengths (5 of 7 versus 2 of 7 for $k = 7$; SI Appendix, Table S3)."*

Sentence 7's "fewer factors" framing repeats sentence 4's "twice as many factors" / "fragmented structure" point. The novel content of sentence 7 — adjustment retention (Table S6) and robustness across supervision strengths (Table S3) — is genuinely new, but it sits behind the redundant parsimony framing and reads as supporting evidence for parsimony rather than as separate points.

### Subtle nuance: the two claims in sentence 7 come from two different tables

The sentence merges two distinct claims with a comma:

- *"Both methods retain significance after classifier adjustment (Table S6)"* — a per-configuration claim about the BO-selected DeSurv $k = 3$ and NMF $k = 7$.
- *"$k = 3$ is more robust across supervision strengths (5 of 7 vs 2 of 7 for $k = 7$; Table S3)"* — a robustness-across-supervision claim aggregating multiple $\alpha$ values.

These are two different tables answering two different questions. The current sentence's comma-separated structure makes the relationship slightly fuzzy. A clean rewrite separates them.

### Three drafts

**Option A — Minimal: drop "with fewer factors" (~−3 words):**

Proposed:
> "Even so, DeSurv at $k = 3$ achieves independent prognostic value: both methods retain significance after classifier adjustment (SI Appendix, Table S6), and $k = 3$ is more robust across supervision strengths (5 of 7 versus 2 of 7 for $k = 7$; SI Appendix, Table S3)."

Two-word delete. The novel content is now primary; parsimony stays implied from sentence 4.

**Option B — Acknowledge prior parsimony explicitly (~net 0 words, preferred):**

Proposed:
> "**Beyond this parsimony**, DeSurv at $k = 3$ achieves independent prognostic value: both methods retain significance after classifier adjustment (SI Appendix, Table S6), and $k = 3$ is more robust across supervision strengths (5 of 7 versus 2 of 7 for $k = 7$; SI Appendix, Table S3)."

"Beyond this parsimony" references the earlier "twice as many factors" framing, marking the rest of the sentence as additional content rather than restating the parsimony argument. Net-zero word delta with cleaner paragraph flow.

**Option C — Restructure to separate the two distinct claims (~+5 words, most rigorous):**

Proposed:
> "**Beyond this parsimony, DeSurv at $k = 3$ retains a significant adjusted hazard ratio after conditioning on PurIST and DeCAF (SI Appendix, Table S6) and remains the more robust solution across supervision strengths — significant at 5 of 7 supervision strengths vs. 2 of 7 for $k = 7$ (SI Appendix, Table S3).**"

Drops the "both methods retain significance" framing (technically right but rhetorically muted — by ending with NMF's weaker robustness, the "both retain" beat is undercut anyway). Names DeSurv's retained adjusted HR explicitly and treats the robustness as a separate claim. Two-table claims are clearly separated, no merger ambiguity.

### Recommendation: Option B

Net-zero word delta. Explicitly acknowledges the prior parsimony argument ("Beyond this parsimony") so the reader knows why it's not being restated. Preserves the existing rhetorical structure of "both methods retain + DeSurv more robust." Option A is the lightest touch (just dropping two words); Option C is more rigorous but rewrites the sentence's structure more aggressively.

### Cost/benefit

- **Word delta:** −3 (Option A), 0 (Option B), or +5 (Option C).
- **No bib changes, no new analyses.**
- **Reader benefit:** removes within-paragraph parsimony repetition, foregrounds the genuinely novel content (adjustment retention + supervision-strength robustness).
- **Risk:** zero. Substantive content unchanged; only sentence framing shifts.

### Coordination with §17, §19, §20

When applying §17 (Discussion reframing of theoretical extension) and the §19/§20 contribution-claim foreshadowing in earlier sections, the Discussion paragraph at line 10 may need a coordinated lighter touch:

- §19a moves the plateau range to Results (line 206). When that's applied, the Discussion's plateau sentence becomes a callback rather than first-introduction.
- §19b/§20a move the "NMF needs BO" argument to Results (line 355). When that's applied, the Discussion's "NMF only reaches $k = 7$ through BO with $\alpha = 0$" sentence becomes a callback or can be lighter.
- §20b/§20c move the model-selection-framework contribution claim to abstract/intro. When that's applied, the Discussion's "model selection framework is itself part of the methodological contribution" sentence becomes confirmation rather than first-introduction.
- §21 (this section) handles the "fewer factors" within-paragraph redundancy independently of these moves.

If §19a + §19b + §20 are all applied, the line 10 Discussion paragraph could potentially shrink considerably — much of its empirical content will already be in Results. Worth a once-over of the entire paragraph after applying §19/§20 to verify it reads as elaboration rather than redundant first-introduction.

### Application checklist

1. Choose Option A (minimal), B (preferred), or C (most rigorous).
2. Apply the chosen draft to `paper/05_discussion_REVISED.Rmd:10`, replacing the existing sentence 7.
3. After applying §17/§19/§20, do a once-over of the entire line 10 Discussion paragraph to check whether other sentences have become redundant given the earlier-section additions.
4. Re-render with `make paper` and confirm the paragraph still reads cohesively with the revision.

---

## 22. Remove redundant "Extending DeSurv to..." sentence in Discussion limitations paragraph

### Why this matters

The Discussion limitations paragraph at `paper/05_discussion_REVISED.Rmd:14` contains two consecutive sentences with verbatim parallel structure:

> "Extending DeSurv to cohorts with neoadjuvant or multimodal therapy is therefore a natural and clinically important next step. **Extending DeSurv to multi-omics integration, spatial transcriptomics, and alternative outcome models is a natural next step.**"

Both start with "Extending DeSurv to..." and end with "natural ... next step." Pure structural repetition of a copy-paste/edit residue.

Worse: the second sentence's "spatial transcriptomics" mention is partly redundant with the very next sentence in the paragraph (*"Although single-cell and spatial transcriptomics now resolve tumor compartments directly..."*), so it does double-redundant duty.

The two sentences address different things despite the parallel structure:

- **Sentence A** (kept): logically follows from the preceding limitation discussion (treatment-naive cohorts → extension to neoadjuvant cohorts). Earns its place.
- **Sentence B** (proposed removal): generic future-work filler about multi-omics / spatial / alternative outcome models that doesn't connect to any specific limitation discussed earlier.

### Three drafts

**Option A — Drop sentence B entirely (preferred, ~−15 words):**

Just remove the second sentence. The next paragraph already discusses spatial transcriptomics. Multi-omics and alternative outcome models become unmentioned, but they're generic future-work directions that don't materially strengthen the Discussion. Net delta around −15 words.

**Option B — Combine with semicolon, different structure for second clause (~−5 words):**

> "Extending DeSurv to cohorts with neoadjuvant or multimodal therapy is therefore a natural and clinically important next step; multi-omics integration, spatial transcriptomics, and alternative outcome models represent further directions."

Keeps all content. Uses semicolon transition + "represent further directions" instead of parallel "is a natural next step."

**Option C — Single-sentence list (~−5 words):**

> "Natural next steps include extending DeSurv to cohorts with neoadjuvant or multimodal therapy — clinically important given current standard of care — alongside multi-omics integration, spatial transcriptomics, and alternative outcome models."

Loses the explicit "clinically important next step" emphasis on neoadjuvant.

### Recommendation: Option A

Drop sentence B entirely. Reasons:

1. Sentence B's content is partially handled by the follow-on paragraph (spatial transcriptomics); multi-omics and alternative-outcome-models additions are generic future-work filler.
2. Sentence A's neoadjuvant point earns its place because it follows from the actual limitation discussed earlier (treatment-naive cohorts). Sentence B doesn't have a parallel anchor.
3. Limitations sections benefit from focus. A reader reading "Third, treatment-naive cohorts → extension to treated cohorts is the next step" gets a coherent argument. Adding "...also we should do multi-omics, spatial, alternative outcome models" feels like the limitation paragraph wandered into generic future-work.

### Cost/benefit

- **Word delta:** −15 (Option A), −5 (Option B), −5 (Option C).
- **No bib changes.**
- **Reader benefit:** removes structural repetition; tightens the limitations paragraph.
- **Risk:** zero. Substantive content unchanged or arguably improved.

### Application checklist

1. Apply Option A (delete sentence B), or B/C if you want to preserve the multi-omics / spatial / alternative-outcome-models mentions.
2. Re-render with `make paper` and confirm the limitations paragraph reads cleanly.

---

## 24. Strengthen the Discussion's closing paragraph: compress repetition with line 6 and add specific PDAC achievement + methodological contribution

### Update (2026-05-05): paragraph-level preview reveals scRNA-seq complementarity redundancy with line 14

Paragraph-level preview of the Discussion exit (line 14 post-§22+§23 + line 16 post-§24 Option B) revealed mild redundancy on the scRNA-seq complementarity message. Line 14 (post-§22+§23) already establishes complementarity twice: *"...patient-level survival modeling may complement deconvolution-based approaches; DeSurv's gene programs... can serve as a bridge between bulk-derived prognostic models and emerging cellular-resolution datasets."* Line 16's §24 Option B closing then adds: *"...is more likely to expand alongside, rather than be displaced by, single-cell technologies."* — third or fourth restatement of the same beat in close succession.

**Two coordination options:**

1. **Trim line 14's "may complement deconvolution-based approaches" phrase** while preserving the unique "bridge between bulk-derived... and emerging cellular-resolution datasets" claim. Line 16's "alongside vs displaced by" then carries the complementarity beat without competing. This is the recommended option because it preserves §24's closing momentum (which the §24 entry's "Why this matters" explicitly argues was needed to avoid the original closing's circular generality).
2. **Drop the "alongside, rather than be displaced by, single-cell technologies" tail from §24 Option B.** Line 16 then ends on *"...wherever clinical outcomes are available."* — slightly abrupt but acceptable since line 14 carries complementarity. Net delta: −20 words on the §24 closing.

**Recommendation: Option 1 (line 14 trim).** Preserves §24's forward-looking momentum at the closing while removing the redundancy. Specifically: in the line 14 sentence *"As single-cell cohorts with clinical annotation grow, patient-level survival modeling may complement deconvolution-based approaches; DeSurv's gene programs..."*, drop the *"may complement deconvolution-based approaches; "* clause and make it read *"As single-cell cohorts with clinical annotation grow, DeSurv's gene programs..."* — preserves the bridge claim without restating complementarity.

This coordination is in addition to the §22 (Cluster F) and §23 (Cluster H) edits that already touch line 14.

### Why this matters

The Discussion's closing paragraph at `paper/05_discussion_REVISED.Rmd:16` is structured as a high-altitude recap, but it has two issues:

1. **Cross-paragraph repetition with Discussion line 6.** Sentence 1 says: *"the standard workflow... incurs a structural cost whenever the highest-variance signals are not the most prognostic."* Line 6 already says: *"the principle that outcome-guided dimensionality reduction targets different subspaces than variance-driven reduction extends beyond cancer genomics."* These are the same variance-prognosis-misalignment frame phrased differently. Closings typically recap, but recapping the *exact theoretical frame* opened in the Discussion's first paragraph is unusual.
2. **No specific accomplishments named.** The closing is missing: (a) the Classical-tumor + iCAF-coupling finding, (b) 5-cohort transportability with pooled HR = 1.50, (c) the two-part methodological contribution (supervised gradient + BO/1-SE rule), (d) the parsimonious $k = 3$ + adjustment-retention story. Sentence 3's generic three-part list ("recovers programs directly, quantifies survival contributions, produces compact signatures") could describe almost any survival-supervised factorization method. A reader walks out of the paper without anything concrete to remember.

Sentence 4 ("The framework is general to any setting where NMF-based deconvolution is applied and clinical outcomes are available") is essentially circular ("works wherever it works"). A stronger forward-looking claim would close with momentum.

### Audit table — what each closing-paragraph sentence currently does

| Sentence | Content | Repeats earlier section? |
|---|---|---|
| 1 | Variance-prognosis misalignment is a structural cost | Yes — Discussion line 6 (sufficient dim reduction frame); intro multiple paragraphs; abstract |
| 2 | Pan-cancer generalization via aran2015systematic | Yes — intro line 17 already cites Aran 2015 for the same generalization |
| 3 | DeSurv recovers programs / quantifies $\Delta\ell$ / transfers by projection | Yes — abstract + intro + Discussion line 6 all state these |
| 4 | "Framework is general wherever NMF + outcomes apply" | Partial — circular generality claim |

Conclusion: the entire paragraph is recap. That's not necessarily wrong for a closing, but it's high-altitude recap that doesn't bring the specific accomplishments forward.

### Three options

**Option A — Light tightening (combine sentences 1-2, leave the rest, ~−10 words):**

> "**The standard workflow in cancer transcriptomics — discovering expression programs unsupervised and evaluating their clinical relevance retrospectively — incurs a structural cost whenever the highest-variance signals are not the most prognostic, likely the norm rather than the exception in cancer types where tumor purity and tissue composition dominate expression variation [@aran2015systematic].** DeSurv provides a principled alternative: by aligning the factorization objective with the clinical question from the outset, it recovers prognostic programs directly, quantifies their survival contributions during discovery, and produces compact signatures that transfer to new patients by projection. The framework is general to any setting where NMF-based deconvolution is applied and clinical outcomes are available."

Combines sentences 1-2 to remove cross-paragraph repetition with line 6. Keeps everything else.

**Option B — Substantive revision: name PDAC achievement + specific methodological contribution (~net 0 words, preferred):**

> "The standard workflow in cancer transcriptomics — unsupervised discovery followed by retrospective evaluation — incurs a structural cost whenever the highest-variance signals are not the most prognostic, likely the norm in cancer types where tumor purity and tissue composition dominate expression variation [@aran2015systematic]. DeSurv resolves this by aligning factorization with the clinical question from the outset: **in PDAC, this recovered a Classical-tumor + iCAF-coupling program that transferred to five independent cohorts (pooled HR per SD = 1.50) at a parsimonious $k = 3$ that survives PurIST/DeCAF adjustment**. The combination of **a survival-supervised gradient and joint Bayesian optimization with the one-standard-error rule** provides a principled framework for extracting outcome-relevant programs from bulk tumor transcriptomes wherever clinical outcomes are available — a setting that, as clinically annotated cohorts accumulate, is more likely to expand alongside than be displaced by single-cell technologies.

Net effect: same word count, but the closing now contains:

- A specific PDAC empirical achievement (Classical+iCAF, 5 cohorts, HR 1.50, $k = 3$, adjustment-retention)
- The two-part contribution framing (gradient + BO/1-SE rule), per §20
- A forward-looking claim that pre-empts the "but what about scRNA-seq?" question (per §4 / §23)

**Option C — Aggressive rewrite (forward-looking, deployment-focused, ~+15 words):**

> "DeSurv addresses a structural cost in cancer transcriptomics: the standard discover-then-evaluate workflow misallocates factor capacity whenever the highest-variance signals are not the most prognostic — likely the norm in tumor types where purity and tissue composition dominate expression variation [@aran2015systematic]. **In PDAC, survival-supervised factorization recovered a Classical-tumor + iCAF-coupling program that transferred to five independent cohorts (pooled HR per SD = 1.50) at a parsimonious $k = 3$, retaining significance after adjustment for established classifiers**. The combination of supervised gradient and joint Bayesian optimization with the one-standard-error rule is the methodological contribution; the resulting transportable, projection-scorable factors are the deliverable. **As clinically annotated bulk-transcriptomic cohorts continue to underpin CLIA-validated diagnostics and multi-cohort trial designs, this approach generalizes to any cancer type where the variance-prognosis misalignment is severe enough to bias unsupervised discovery.**"

Most aggressive. Compresses theoretical frame, foregrounds PDAC achievement, names methodological contribution explicitly, ends with deployment-scenario claim that ties to §23 (CLIA-validated diagnostics).

### Recommendation: Option B

Does the most work for the smallest change:

1. Compresses the cross-paragraph repetition with line 6 by combining sentences 1-2.
2. Names specific PDAC achievement (Classical+iCAF, 5 cohorts, HR 1.50, $k = 3$) — gives the reader something concrete to remember.
3. States the two-part methodological contribution (gradient + BO/1-SE) — picks up §20's framing.
4. Forward-looking but defensible about scRNA-seq complementarity rather than displacement — picks up §4 / §23.

Net word delta: ~0 (additions balance compressions).

Option A is a lighter touch if you don't want to revise the closing substantively. Option C is more aggressive but commits the closing to deployment-scenario framing that may be heavier than needed.

### Coordination notes

- **§17 (Discussion line 6 reframing).** If §17 is applied (closing the empirical-question loop the intro opens), Option B's sentence 1 becomes a clean callback rather than redundant restating. The two paragraphs (line 6 and line 16) then complement: line 6 = "we resolved the empirical question"; line 16 = "this matters because variance-prognosis misalignment is the norm."
- **§20 (two-part contribution framing).** Option B's sentence 3 explicitly names "supervised gradient and joint Bayesian optimization with the one-standard-error rule" — consistent with §20's framing of the two-part contribution. If §20b or §20c is applied earlier in the paper, the closing's mention becomes confirmation rather than introduction.
- **§23 (CLIA-validated diagnostic assays).** Option C's deployment-scenario sentence references CLIA-validated diagnostics directly. Option B's "more likely to expand alongside than be displaced by single-cell technologies" implies the scope without naming specifics.
- **§4 (scRNA-seq preemption in intro).** The closing's "alongside than be displaced by" phrasing in Options B and C reinforces §4's intro framing. Apply consistently.

### Cost/benefit

- **Word delta:** −10 (Option A), 0 (Option B), +15 (Option C).
- **No bib changes, no new analyses.** Specific empirical numbers (Classical+iCAF, HR 1.50, $k = 3$) come from elsewhere in the paper.
- **Reader benefit:** the closing brings the specific accomplishments forward and ends with momentum rather than a circular generality claim.
- **Risk:** very low. All claims are already established empirically elsewhere; the closing just makes them explicit at the wrap-up point.

### Application checklist

1. Choose Option A, B, or C based on how substantive a closing rewrite you want.
2. Apply the chosen draft to `paper/05_discussion_REVISED.Rmd:16`, replacing the existing closing paragraph.
3. If applying §17, §20, §22, §23 elsewhere, verify language consistency with the closing (especially the "supervised gradient + BO/1-SE rule" framing and the "alongside than be displaced by" scRNA-seq complementarity framing).
4. Re-render with `make paper` and confirm the closing reads cleanly with the revision.

---

# Cluster G — Convergent biological validation

Singleton — surfaces the SI-buried PurIST/DeCAF concordance claim in main-text Results, with built-in bridge to the existing adjusted-HR (independence) claim to prevent the apparent tension between "DeSurv recovers established subtypes" and "DeSurv adds info beyond them" from registering with a careful reader.

**Coordination with Cluster F:** Discussion line 12 could be lightly enhanced to coordinate with §28's concordance-plus-independence framing in Results.

**Coordination with Cluster J:** §29a's revised wording explicitly invokes Fig. S7 (concordance) — apply §28 and §29a together so main-text and SI tell the same story.

**Convergent-validation cross-cluster thread:** The convergent-validation argument runs through four sections in three clusters: §10 (Cluster C, Path A — NMF k=5 acknowledgment), §17 (Cluster F, Discussion theory reframing), §28 (this cluster, main-text Results addition), and §29a (Cluster J, SI Section 20 correction). Apply as a coordinated set; §29a must precede main-text application of §10 + §28. Within §28, the recommended placement is **end of the first paragraph** (immediately after the existing adjusted-HR sentence), not the originally proposed end-of-second-paragraph location — this puts the concordance bridge adjacent to its referent rather than reaching across a paragraph boundary.

---

## 28. Surface the convergent-validation claim (DeSurv risk groups concordantly recover Basal-like + proCAF subtypes) in main-text Results

### Update (2026-05-05): paragraph-level preview revises placement and bridge tense

Paragraph-level preview of the validation paragraph (post-§14 + §15 + §28) revealed two refinements to the original §28 specification:

**Placement update — move from end of second paragraph (line 357) to end of first paragraph (after the existing adjusted-HR sentence, line 355).** The original §28 placement at line 357 (end of dichotomized-validation paragraph) reaches across a paragraph boundary for its bridge clause: *"...drive its retained adjusted hazard ratio (Table S6)"* refers to the adjusted-HR claim in the *previous* paragraph (line 355). Moving §28 to the end of the first paragraph (immediately after *"DeSurv at $k = 3$ retains a significant hazard ratio after adjustment for PurIST and DeCAF... achieving independent prognostic value more parsimoniously"*) puts the concordance bridge adjacent to its referent. The dichotomized paragraph then stays topically focused on cutpoint stratification.

**Bridge tense correction — "drive" → "captured by".** The original draft says *"adding within-subtype gradations that drive its retained adjusted hazard ratio"* — but Table S6 already showed this; the bridge describes what the table captures, not what new evidence drives it. Corrected: *"adding the within-subtype gradations captured by this retained adjusted hazard ratio."*

**Revised §28 Option 2-refined draft (incorporating both updates) — superseded by post-Amber-input revision below.**

> "Although DeSurv was trained without molecular subtype labels, the resulting high-risk group is significantly enriched for the two established poor-prognosis subtypes (Basal-like by PurIST, proCAF by DeCAF; Fisher's exact, SI Appendix, Fig. S7) — recovering the same consensus biology while adding the within-subtype gradations **captured by** this retained adjusted hazard ratio (Table S6)."

### Update (2026-05-05, post-Amber-input): drop "within-subtype gradations captured by Table S6" bridge; drop "significantly"

Amber flagged two issues with the prior draft:

1. **"Within-subtype gradations captured by Table S6" is unsupported.** Table S6 reports the adjusted hazard ratio for the *continuous linear predictor*, not the *dichotomized risk groups*. The "high-risk group" comes from Fig 4B/4C dichotomization. Forcing a bridge between the dichotomized-group enrichment (Fig S7) and the linear-predictor adjusted HR (Table S6) via "within-subtype gradations" overclaims what the table directly shows.
2. **"Significantly" without inline test wording**. Fisher's exact P-values are in Fig S7's caption, but main-text "significantly enriched" reads as forward-leaning. Softer "enriched" lets Fig S7 carry the test where it lives.

**Revised §28 draft (preferred, post-Amber-input):**

Insert at end of first paragraph (`paper/04_results_REVISED.Rmd:355`), immediately after *"...achieving independent prognostic value more parsimoniously"*:

> "Although DeSurv was trained without molecular subtype labels, the resulting high-risk group was enriched for the two established poor-prognosis subtypes (Basal-like by PurIST, proCAF by DeCAF; Fisher's exact, SI Appendix, Fig. S7), and the DeSurv linear predictor retained significant prognostic value after adjustment for these classifiers (Table S6) — capturing established consensus biology while adding independent prognostic information."

Net effect:

- **Treats the dichotomization (Fig S7) and linear-predictor adjusted HR (Table S6) as separate complementary observations** rather than forcing the unsupported "gradations captured by Table S6" bridge.
- **Drops "significantly"** in favor of "enriched" (softer; Fisher's exact in Fig S7 carries the test).
- **Synthesis sentence "capturing established consensus biology while adding independent prognostic information"** explicitly states the dual-property point Amber's "Are you trying to say... [both]?" question identified.

Word delta: net 0–5 vs. prior Option 2-refined (slightly cleaner phrasing). Same length, materially more defensible.

**§29a coordination:** the §29 SI Section 20 wording revision (see below in §29 entry) similarly drops "while concordantly capturing" in favor of explicit "combines biological concordance and methodological complementarity" framing — both edits reflect the same Amber feedback applied at SI and main-text levels.

The original Options 1, 2-refined, 3 below remain valid as wording references for the placement-and-tense issues, but the post-Amber draft above supersedes them on the bridge wording.

### Why this matters

The SI passage at `paper/si_appendix.Rmd:1408` (section "Overlap of DeSurv risk groups with known molecular subtypes") makes a substantive validation claim that is **completely absent from the main text**:

> "The high-risk group was significantly enriched for Basal-like (PurIST) and proCAF (DeCAF) subtypes, consistent with the known poor prognosis of these molecular classes. This concordance confirms that DeSurv's survival-driven factorization recovers clinically relevant biology without requiring subtype labels during training."

The main text mentions PurIST/DeCAF in three places (Results lines 355, 609; Discussion line 12) — all framing them as **adjustment variables** for an *independence* claim ("DeSurv adds prognostic info beyond PurIST/DeCAF"). The SI's **convergent-validation** claim ("DeSurv risk groups concordantly recover the same consensus biology") is rhetorically distinct and currently buried.

The convergent-validation claim is rhetorically powerful and uniquely positioned for PNAS:
- An unsupervised-with-respect-to-labels model independently rediscovers the field's two consensus poor-prognosis subtypes.
- This pre-empts the "but does this recover real biology?" reviewer concern.
- It complements §10's "Classical+iCAF coupling discovered de novo" framing.

### Apparent tension between the two claims (and why it's resolvable)

A careful reader will notice an apparent tension between the two claims:

- **Independence (existing):** DeSurv adds prognostic content beyond PurIST/DeCAF (adjusted HR = 1.33, $P < 0.001$).
- **Concordance (proposed):** DeSurv risk groups concordantly recover Basal-like and proCAF subtypes (Fisher's exact significance).

These would only conflict if DeSurv's risk score were a *deterministic function* of PurIST and DeCAF — in which case conditional on (PurIST, DeCAF), DeSurv would add zero information and the adjusted HR would collapse to 1.0. **They don't conflict** because:

- "Concordance" is a weak association statement (Fisher's exact rejects the null of independence; doesn't say identical).
- "Independence beyond" is a residual-content statement (adjusted HR captures information not explained by the binary covariates).

A continuous score can be statistically associated with binary classifications **and** add residual content beyond them. The thermometer/fever-indicator analogy: a continuous body-temperature reading is concordant with a binary fever indicator (≥38°C) — the high-temp group is enriched for "fever" status — but the continuous reading also distinguishes 38.5°C from 41°C, which the binary collapses. DeSurv's continuous risk score is the thermometer; PurIST and DeCAF are binary classifiers. Concordance: same buckets. Independence: within-bucket gradations.

The two claims **jointly establish** something stronger than either alone:

> "DeSurv recovers the consensus poor-prognosis biology (Basal-like, proCAF) without supervision, AND captures additional within-subtype prognostic gradations that the binary classifiers can't resolve."

This is a **complementarity** argument that frames DeSurv as a continuous, biologically grounded refinement of established subtyping rather than (a) a competing classifier or (b) a pure black-box score.

### Implication for the proposed sentence

A simple "the high-risk group was enriched for Basal-like and proCAF" sentence inserted into the main text **without bridging to the existing adjusted-HR claim** creates the apparent tension above. A careful reader hits "DeSurv adds info beyond PurIST/DeCAF" (line 355) immediately followed by "DeSurv risk groups align with PurIST/DeCAF subtypes" (the new sentence) and wonders. The proposed sentence needs to **bridge** to the independence claim explicitly, not just sit alongside it.

### Three drafts (recommended: Option 2-refined)

Insert at end of the validation paragraph at `paper/04_results_REVISED.Rmd:357`, after the existing "finer prognostic discrimination" sentence and before the closing "Together, these results indicate..." sentence.

**Option 1 — Compact, no bridge (~+25 words). NOT RECOMMENDED.**

> "Furthermore, the DeSurv high-risk group was significantly enriched for established poor-prognosis subtypes — Basal-like (PurIST) and proCAF (DeCAF) — across the pooled validation cohort (Fisher's exact test, SI Appendix, Fig. S7), confirming that the survival-driven factorization recovers clinically relevant biology without requiring subtype labels during training."

This is rhetorically clean but **creates the apparent tension** with line 355's adjusted-HR claim. A careful reader will notice. Avoid unless there's a clear separation between the two claims in the manuscript flow that prevents the tension from registering.

**Option 2-refined — Bridge + concordance (~+30 words, preferred):**

> "Furthermore, although DeSurv was trained without molecular subtype labels, the resulting high-risk group was significantly enriched for the two established poor-prognosis subtypes (Basal-like by PurIST, proCAF by DeCAF; Fisher's exact, SI Appendix, Fig. S7) — **recovering the same consensus biology while adding within-subtype gradations that drive its retained adjusted hazard ratio (Table S6)**."

Single sentence does three things:

1. Surfaces the convergent-validation finding (concordance with consensus subtypes).
2. Names the surprise (trained without subtype labels).
3. **Explicitly bridges to the adjusted-HR claim** via "within-subtype gradations" — exactly the right mechanistic description of how a continuous score is both concordant with and adds beyond binary classifiers.

This is the recommended draft. It pre-empts the tension the user flagged by integrating the bridge into the same sentence.

**Option 3 — Two-sentence explicit complementarity (~+50 words):**

> "Furthermore, the DeSurv high-risk group was significantly enriched for Basal-like (PurIST) and proCAF (DeCAF) subtypes (SI Appendix, Fig. S7), confirming that the survival-driven factorization recovers clinically relevant biology without subtype labels during training. Combined with the adjusted-HR analysis (Table S6) showing that DeSurv's $k = 3$ linear predictor retains independent prognostic value beyond these classifiers, the validation establishes both biological concordance and methodological complementarity with established PDAC subtyping."

Most explicit framing of the two-claim complementarity. Names "biological concordance" and "methodological complementarity" as separate beats. Heavier (~+50 words) but rhetorically the most rigorous if a co-author wants the complementarity argument made unambiguous.

### Recommendation: Option 2-refined

Reasons:

1. **Resolves the apparent tension in a single sentence** rather than requiring a separate complementarity beat.
2. **Names the surprise** ("trained without molecular subtype labels") that makes the convergent-validation finding land for a PNAS general audience.
3. **Bridges to the adjusted-HR claim** via "within-subtype gradations" — mechanistically accurate description of why concordance and independence coexist.
4. **~+30 words is reasonable** for a substantive validation claim that's currently invisible above the SI fold.

Option 1 is acceptable only if you trust the reader to bridge the two claims themselves. Option 3 is overkill unless a co-author specifically wants the complementarity argument made unambiguous — Option 2-refined captures the same content more economically.

### Coordination with existing PurIST/DeCAF mentions

After applying §28 Option 2-refined, the main text will mention PurIST/DeCAF in four places:

| Location | Current/post-§28 framing |
|---|---|
| Results line 355 | Independence: "retains a significant HR after adjustment for PurIST and DeCAF" |
| Results line 357 (new) | **Concordance + bridge: "recovering the same consensus biology while adding within-subtype gradations"** |
| Results line 609 | Independence: "retained significance at [N] of 7 supervision strengths after adjustment" |
| Discussion line 12 | Independence: "DeSurv's factors capture independent prognostic signal beyond what those classifiers already explain" |

The new line 357 sentence is the **only concordance mention** — the other three remain independence claims. This distribution is correct: concordance is a one-time observation; independence shows up wherever adjusted analyses are reported. With Option 2-refined's bridge, the reader doesn't need to reconcile the tension because the bridge is built into the sentence itself.

### Cost/benefit

- **Word delta:** +25 (Option 1), +30 (Option 2-refined), +50 (Option 3).
- **No bib changes, no new analyses.** SI Fig. S7 already exists.
- **Reader benefit:** surfaces a substantive validation claim that's currently invisible in main text; converts a buried-SI finding into a foregrounded contribution argument.
- **Risk:** very low if Option 2-refined or Option 3 used (both pre-empt the tension). Higher with Option 1 (creates apparent contradiction with existing adjusted-HR claim).

### Application checklist

1. Apply Option 2-refined to `paper/04_results_REVISED.Rmd:357`, inserting the new sentence between the "finer prognostic discrimination" sentence and the "Together, these results indicate..." closing sentence.
2. Re-read the surrounding paragraph after applying — verify the flow:
   - DeSurv pooled HR (line 355)
   - NMF pooled HR (line 355)
   - Dichotomized DeSurv KM (line 357)
   - Dichotomized NMF KM (line 357)
   - DeSurv vs NMF group sizes (147 vs 46)
   - **New: DeSurv risk groups recover established subtypes + bridge**
   - "Together, these results indicate..." closer
3. Verify that the closing sentence ("Together, these results indicate that incorporating survival information during factorization yields gene programs with more transportable prognostic associations than the standard discover-then-evaluate approach") still flows naturally given the new content.
4. If §16 (closer wording: "more transportable prognostic associations than standard NMF — the canonical example of the discover-then-evaluate paradigm") is also being applied, coordinate the closer wording with the new convergent-validation content.
5. Re-render with `make paper` and confirm the validation paragraph reads cleanly with the addition.

---

# Cluster H — scRNA-seq positioning + clinical deployment

Two edits making the bulk-vs-single-cell case at intro P1 (with CIBERSORTx/BayesPrism citations) and Discussion (CLIA-validated assays broadening). Apply with coordination: §4 application requires retrim of the corresponding Discussion paragraph (line 14) to avoid duplicate defenses 6 pages apart.

**Coordination with §24** (Discussion close): Option B/C of §24 already adds forward-looking statement about scRNA-seq complementarity; ensures consistent language across intro → Discussion close.

---

## 4. Introduction — preempt the "why bulk in the scRNA-seq era?" question

### Update (2026-05-05): compact version preferred over long; long retained as fallback

External review flagged that the **long version** (with CIBERSORTx + BayesPrism citations, ~+55 words) may distract from the introduction's existing arc (supervised factorization, BO, rank selection, external validation, PDAC biology, NMF sensitivity). The introduction already establishes that bulk profiling remains dominant and that each profile is a composite — a long single-cell-vs-bulk discussion may overload an already-dense intro.

**Revised recommendation: compact version preferred for application.** The compact draft below:

> "Single-cell and spatial atlases clarify the cellular states present in PDAC, but cohort-scale survival annotation remains largely bulk-transcriptomic; supervised decomposition of bulk mixtures therefore remains necessary for population-scale prognostic discovery."

Net delta: ~+30 words. Handles the "why bulk in the scRNA-seq era?" reviewer concern without turning the intro into a methods discussion. The long version (with CIBERSORTx + BayesPrism citations) remains documented below as a fallback if a reviewer specifically asks for the scRNA-seq-aware method citations or if the compact version reads as too thin.

**Trade-off:** Long version cites methods (CIBERSORTx, BayesPrism) that exemplify bulk + scRNA-seq complementarity, which is rhetorically powerful. Compact version makes the same point without specific method citations — defensible, less detailed.

**Recommendation: try compact version first; promote to long only if reviewer feedback specifically requests citations or expanded discussion.**

### Update (2026-05-05): concrete Discussion line 14 retrim specified

The original §4 entry flagged a coordination requirement: "When applied, retrim the corresponding paragraph in the Discussion ('Although single-cell and spatial transcriptomics now resolve tumor compartments directly...') to avoid two near-identical defenses repeated in both Intro and Discussion." Paragraph-level preview of intro P1 (post-§4) and Discussion line 14 (post-§22 + §23) made the specific retrim concrete:

**After §4 lands at intro P1 (carrying the cohort-size-rarity argument), drop the equivalent sentence from Discussion line 14:**

- **Drop from Discussion line 14:** *"Although single-cell and spatial transcriptomics now resolve tumor compartments directly, the cohort sizes required for stable survival modeling (hundreds of patients with mature follow-up) remain available primarily in bulk expression data, and"*
- **Keep:** *"bulk profiling continues to underpin clinical molecular stratification — from trial programs such as COMPASS to CLIA-validated diagnostic assays in routine clinical care — a deployment landscape favored by the cost, FFPE compatibility, and standardization requirements of the clinical workflow. As single-cell cohorts with clinical annotation grow, patient-level survival modeling may complement deconvolution-based approaches; DeSurv's gene programs, which align with cell-type signatures confirmed by independent single-cell studies, can serve as a bridge between bulk-derived prognostic models and emerging cellular-resolution datasets."*

This trim removes the cohort-size-rarity argument (now carried by intro P1 §4) and keeps the §23 CLIA-validated diagnostics + the unique "bridge" claim that intro P1 doesn't replicate.

**Coordination with §24's update:** §24's update separately flags trimming "may complement deconvolution-based approaches" from line 14 (to remove the redundancy with §24's "alongside vs displaced by" closer). If both §4 and §24 retrims are applied, line 14 reads: *"...bulk profiling continues to underpin clinical molecular stratification — from trial programs such as COMPASS to CLIA-validated diagnostic assays in routine clinical care — a deployment landscape favored by the cost, FFPE compatibility, and standardization requirements of the clinical workflow. As single-cell cohorts with clinical annotation grow, DeSurv's gene programs, which align with cell-type signatures confirmed by independent single-cell studies, can serve as a bridge between bulk-derived prognostic models and emerging cellular-resolution datasets."* — substantially shorter than the current line 14 paragraph, but retains the §23 CLIA content and the bridge claim. Both retrims are coordinated and compatible.

### Why this matters

A methods-aware reviewer (or skeptical reader) will reasonably ask: _"bulk deconvolution feels like an older technique now that we have single-cell and spatial transcriptomics. Why is this paper even needed?"_ The current intro doesn't address this head-on. The opening of P1 says bulk "remains the dominant technology" — true, but soft, and easy to dismiss as inertia rather than necessity. Better to name the question and answer it in two sentences before the reviewer has time to form the objection on their own.

### The four-part answer (for our own reference; the inserted text condenses pieces 1 + 3)

1. **Survival-linked cohorts are overwhelmingly bulk.** TCGA, CPTAC, ICGC PACA-AU, and the five validation cohorts in this paper are all bulk. scRNA-seq cohorts with linked long-term survival are rare, small, and underpowered for cohort-scale survival modeling — and the cost / FFPE / logistics constraints make population-scale single-cell survival cohorts unlikely to materialize in the near term.
2. **Clinical biomarker assays are bulk.** Anything deployable in a CLIA lab today (NanoString, FFPE RNA-seq, qPCR panels) operates on bulk material; a signature designed on bulk has a deployment path that a single-cell signature does not.
3. **scRNA-seq and bulk deconvolution are complementary, not substitutive.** Single-cell tells you _what cell states exist_; bulk + survival supervision tells you _which programs predict outcome at population scale_. The leading scRNA-seq-aware bulk methods (CIBERSORTx, BayesPrism) explicitly use single-cell to inform bulk decomposition rather than replacing it — the field consensus is "both."
4. **The variance/prognosis misalignment DeSurv targets is intrinsic to the data type that has the survival outcomes.** scRNA-seq doesn't suffer from compositional dilution because cells are already separated — but it also doesn't have the survival data needed to identify prognostic programs at scale. The problem this paper solves only exists where the outcomes live.

### Proposed edit (P1 of the introduction, line 15)

**Current sentence in P1 (after the Collisson opener, before the "Nonnegative matrix factorization has been instrumental..." sentence):**

> Bulk transcriptomic profiling remains the dominant technology for the large, clinically annotated cohorts required for such survival modeling [@bair2004semi; @aung2018compass], but each expression profile is a composite of malignant cells, cancer-associated fibroblasts (CAFs), immune infiltrates, and other microenvironmental components [@nguyen2024fourteen].

**Proposed (insert two sentences between the existing two clauses):**

> Bulk transcriptomic profiling remains the dominant technology for the large, clinically annotated cohorts required for such survival modeling [@bair2004semi; @aung2018compass]. **While single-cell and spatial transcriptomic atlases have transformed our understanding of tumor cellular composition [@chansengyue2020transcription; @werba2023single], cohort-scale survival annotation remains an almost exclusively bulk-transcriptomic resource — single-cell datasets with linked long-term outcomes are rare and underpowered for population-scale prognostic discovery. Computational decomposition of bulk mixtures therefore remains the practical route to identifying outcome-relevant transcriptional programs from existing clinical cohorts, complementing rather than being replaced by single-cell technologies.** Each expression profile is a composite of malignant cells, cancer-associated fibroblasts (CAFs), immune infiltrates, and other microenvironmental components [@nguyen2024fourteen].

Word delta: ~+55 words. The intro was trimmed from 1087 → 790 in the recent revision; +55 is small relative to that headroom and addresses a substantive defensive concern. No new bib entries needed (`chansengyue2020transcription` and `werba2023single` are already cited later in P2).

### Compact fallback (if word count gets tight)

If the +55 expansion is too much, the same point compresses into a single clause appended to the existing first sentence:

> Bulk transcriptomic profiling remains the dominant technology for the large, clinically annotated cohorts required for survival modeling [@bair2004semi; @aung2018compass] — a setting where single-cell and spatial transcriptomic atlases [@chansengyue2020transcription; @werba2023single] complement rather than replace bulk profiling, since cohort-scale outcome data remains nearly exclusively bulk.

Word delta: ~+25 words. Loses the "intrinsic to the data type that has outcomes" piece but preserves complementarity framing.

### Dials to consider

- **Tone toward scRNA-seq.** Drafts above frame complementarity, not inferiority. If a co-author wants a sharper "and here is what bulk + supervision uniquely solves" beat, the long version can absorb one more clause naming the variance/prognosis misalignment as a bulk-specific problem (point 4 in the four-part answer above).
- **Citation depth.** The two scRNA-seq cites used here (`chansengyue2020transcription`, `werba2023single`) are PDAC-specific. For pan-cancer methods breadth, see the CIBERSORTx / BayesPrism note immediately below.
- **Placement.** The proposed insertion sits at the front of P1, where bulk is introduced. An alternative is a one-sentence acknowledgement at the end of P5 (the contribution paragraph). Front-of-P1 is preferred because it disarms the question before the reviewer has read the methods claim; end-of-P5 leaves the question hanging through the entire intro.

### Note on CIBERSORTx and BayesPrism citations

**Status (verified 2026-05-04):** Neither method is currently cited in `paper/references_30102025.bib`. A grep for `cibersort|newman.*2019|bayesprism|chu.*2022` returned zero hits. The closest existing citations are `nguyen2024fourteen` (deconvolution review) and `peng2019novo` (DECODER) — neither covers the scRNA-seq-informed bulk methodology specifically.

**Why citing them strengthens the §4 preemption.** The argument that "the field consensus is _both_" is rhetorically much stronger when it names the two leading scRNA-seq-informed bulk-deconvolution methods. Both use single-cell references to inform bulk deconvolution rather than replacing it — they are concrete proof that the field already operates on the complementarity model the paper claims. A methods-aware reviewer asking "why bulk in the scRNA-seq era?" is often someone familiar with these tools; citing them signals the authors are too, and shifts the conversation from "is bulk obsolete?" to "given the field already integrates both, here is the supervised contribution we add."

**Suggested bib entries.** These follow the existing key convention (`<lastname><year><firstcontentword>`):

- **CIBERSORTx** — `newman2019determining`
  > Newman AM, Steen CB, Liu CL, Gentles AJ, Chaudhuri AA, Scherer F, Khodadoust MS, Esfahani MS, Luca BA, Steiner D, Diehn M, Alizadeh AA. Determining cell-type abundance and expression from bulk tissues with digital cytometry. *Nature Biotechnology* 37, 773–782 (2019). DOI: 10.1038/s41587-019-0114-2

- **BayesPrism** — `chu2022cell`
  > Chu T, Wang Z, Pe'er D, Danko CG. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. *Nature Cancer* 3, 505–517 (2022). DOI: 10.1038/s43018-022-00356-3

**Where they fit in the §4 long draft.** Add to the second sentence of the proposed insertion, in a new clause that names the scRNA-seq-informed bulk methodology lineage explicitly. Updated draft (additions in **bold**, original §4 long-draft text kept intact otherwise):

> While single-cell and spatial transcriptomic atlases have transformed our understanding of tumor cellular composition [@chansengyue2020transcription; @werba2023single], cohort-scale survival annotation remains an almost exclusively bulk-transcriptomic resource — single-cell datasets with linked long-term outcomes are rare and underpowered for population-scale prognostic discovery. Computational decomposition of bulk mixtures therefore remains the practical route to identifying outcome-relevant transcriptional programs from existing clinical cohorts, **and recent scRNA-seq-informed methods [@newman2019determining; @chu2022cell] reinforce this complementarity by using single-cell references to refine bulk decomposition rather than supplant it.**

Word delta vs. the existing §4 long draft: +21 words (so the long draft grows from ~+55 to ~+76 over the current intro). Still small in the context of the recent 1087 → 790 trim.

**Compact-fallback variant.** If using the short fallback in §4 instead of the long draft, append a single citation to the existing complementarity clause:

> ...a setting where single-cell and spatial transcriptomic atlases [@chansengyue2020transcription; @werba2023single] complement rather than replace bulk profiling **[@newman2019determining; @chu2022cell]**, since cohort-scale outcome data remains nearly exclusively bulk.

**Cost/benefit.** Two new bib entries with DOIs — small relative to the 38 DOIs added in the recent PNAS-compliance pass (per `docs/changes_for_amber.md`). The benefit is real: it preempts a methods-aware reviewer comment with concrete referents instead of a hand-wave at "the field." Recommend including both unless there's a co-author objection to specific tool naming.

**If only one can be cited,** prefer **BayesPrism (`chu2022cell`)** over CIBERSORTx — it is more recent, was published in *Nature Cancer* (closer to the oncology framing of this paper), and explicitly motivates Bayesian single-cell-bulk integration in cancer contexts. CIBERSORTx is the better-known name in the field, but BayesPrism is the more recent and specific methodological match.

---

## 23. Broaden COMPASS-only mention to CLIA-validated diagnostic assays in the Discussion

### Why this matters

The Discussion paragraph at `paper/05_discussion_REVISED.Rmd:14` (just after the "extending DeSurv" sentences) makes the case for bulk transcriptomics by naming a single research effort:

> "...bulk profiling continues to underpin clinical molecular stratification efforts such as **COMPASS** [@aung2018compass]."

COMPASS is one specific exemplar of clinical bulk-RNA-seq use, but the actual landscape is much broader. CLIA-validated bulk-RNA-seq diagnostic assays span Oncotype DX, MammaPrint, PAM50/ProSigna, Decipher, Prolaris, NanoString panels, FFPE-RNA-seq workflows, and an expanding set of pharmacogenomic / companion diagnostic assays. All operate on bulk material because of cost, FFPE compatibility, sample-handling logistics, and regulatory standardization. None of these constraints apply favorably to single-cell technologies in the foreseeable near term.

This is the **same defensive argument as §4 (scRNA-seq preemption in intro)**, applied to the Discussion's bulk-vs-single-cell paragraph. Naming COMPASS alone narrows what should be a broader scope claim. A reviewer who knows the clinical-assay landscape will recognize the broader framing as accurate; a reviewer who doesn't will read it as an authoritative scope statement.

### Three drafts

**Option A — Light broadening (~+5 words):**

> "...bulk profiling continues to underpin clinical molecular stratification efforts such as COMPASS [@aung2018compass] **and CLIA-validated diagnostic assays more broadly**."

Lightest. Adds the category without specifics or rationale.

**Option B — Broaden with rationale (~+20 words, preferred):**

> "...bulk profiling continues to underpin clinical molecular stratification — **from research efforts such as COMPASS [@aung2018compass] to CLIA-validated diagnostic assays in routine clinical care — a deployment landscape favored by the cost, FFPE compatibility, and standardization requirements of the clinical workflow**."

Names CLIA explicitly and gives the why (cost, FFPE, standardization). Doesn't require specific commercial product citations. Strongest version that stays within the existing paragraph's tone.

**Option C — Most explicit (~+30 words):**

> "...bulk profiling continues to underpin clinical molecular stratification — from research efforts such as COMPASS [@aung2018compass] to CLIA-validated diagnostic assays in routine clinical care. **This deployment landscape is unlikely to shift in the near term: the cost, FFPE-compatibility, sample-handling logistics, and regulatory standardization requirements of clinical workflows favor bulk technologies, and single-cell cohorts with clinical annotation remain rare and small.**"

Most explicit. Forward-looking statement about why this isn't going to change soon.

### Recommendation: Option B

Names CLIA explicitly, points at the structural reasons (cost, FFPE, standardization) that explain why bulk dominates clinical use, doesn't name specific commercial products, and stays within the existing paragraph's tone. ~+20 words to convert a one-exemplar defensive claim into a genuinely persuasive scope argument.

**On not naming specific commercial products:** Oncotype DX, MammaPrint, PAM50, Decipher, etc. are all CLIA-validated bulk assays that would exemplify the broader category. But naming specific products introduces complications (citation requirements, possible perception of advertising, ambiguity about which products to name and which to omit). Better to use the generic category.

### Coordination with §4

§4 (scRNA-seq preemption in intro) addresses the same bulk-vs-single-cell argument at intro P1. If §4 is applied, this Discussion sentence becomes a callback to the intro framing rather than a first-introduction. Phrasing should stay consistent — if intro says "single-cell datasets with linked long-term outcomes are rare and underpowered for population-scale prognostic discovery," the Discussion can implicitly reference back via consistent language ("CLIA-validated diagnostic assays in routine clinical care" + "single-cell cohorts with clinical annotation remain rare and small" if going with Option C).

If §4 is not applied, this §23 broadening serves as the paper's primary scope statement on bulk-vs-single-cell deployment.

### Cost/benefit

- **Word delta:** +5 (Option A), +20 (Option B), +30 (Option C).
- **No bib changes.** No new analyses.
- **Reader benefit:** converts a single-exemplar defensive claim into a genuinely persuasive scope argument; pre-empts "you only named COMPASS" reviewer comments.
- **Risk:** very low. Empirical scope claim is accurate; no new factual claims that need separate verification.

### Application checklist

1. Apply Option A, B, or C to `paper/05_discussion_REVISED.Rmd:14` (the COMPASS sentence).
2. If §4 is also being applied, verify language consistency between intro P1 and this Discussion sentence (use "CLIA-validated diagnostic assays" and "single-cell cohorts with clinical annotation remain rare" in both, for example).
3. Re-render with `make paper` and confirm the paragraph reads cleanly.

---

# Cluster I — Mechanical clarity / wording fixes

Two small surface fixes (Fig 3 panel order; Classical/Basal-like W-correspondence label paradox). Each is independent; no coordination required.

**Note on §16:** Section 16 was discussed during drafting but never filed; the current parking-lot file's section numbering jumps from §15 to §17.

**Additional candidate site identified during cluster-expansion analysis:** quick figure-panel ordering audit across Fig 1, 2, 4 to confirm no other §11-style issues. Likely OK based on what's been observed but worth a verification pass.

---

## 11. Fix out-of-order figure-panel citation (Fig. 3C cited before Fig. 3B in standalone references)

### Why this matters

In `paper/04_results_REVISED.Rmd`, the standalone references to Fig. 3 panels appear out of alphabetical order:

| Line | Reference | Order encountered |
|---:|---|---|
| 253 | `Fig. 3A--B` (range, both panels in one citation) | OK |
| 255 (DeSurv paragraph, multiple) | `Fig. 3A` | A |
| **255 (last sentence)** | **`Fig. 3C`** *(forward reference: "quantified below")* | **C — out of order** |
| 257 (NMF paragraph, multiple) | `Fig. 3B` | B (first standalone) |
| 259 | `Fig. 3C` | C (in order) |
| 261 | `Fig. 3D` | D |

The DeSurv (3A) and NMF (3B) paragraphs are structurally parallel — both discuss factor identity, both end with a forward-pointing transition. But only the DeSurv paragraph attaches an explicit `(Fig. 3C)` to its transition. That asymmetry creates the cite-order issue. Forward references like "quantified below" are technically allowed, but PNAS reviewers reading strictly will flag the out-of-alphabetical-order panel citation, especially because the NMF paragraph (line 257) already ends with its own forward-pointing transition ("quantified next") that *implicitly* points to 3C without the explicit tag.

### Proposed edit (preferred — minimum-intervention)

**Drop** the parenthetical Fig. 3C reference from line 255's last sentence. Line 259 cites Fig. 3C at first actual need, and line 257's "quantified next" already signposts the upcoming panel. No information is lost.

**Current (`paper/04_results_REVISED.Rmd:255`, last sentence):**

> The survival contributions of these factors and their relationship to expression variance are quantified below **(Fig. \ref{fig:pdac}C)**.

**Proposed:**

> The survival contributions of these factors and their relationship to expression variance are quantified below.

Net effect: standalone-reference order becomes A → B → C → D, alphabetical and clean.

### Alternative (preserves explicit signpost in alphabetical position)

If the explicit Fig. 3C tag is wanted as a signpost to the upcoming panel, **shift it from the DeSurv paragraph (line 255) to the NMF paragraph's transition (line 257)**:

**Current (`paper/04_results_REVISED.Rmd:257`, last sentence):**

> Whether these structural differences correspond to differences in prognostic information is quantified next.

**Proposed:**

> Whether these structural differences correspond to differences in prognostic information is quantified next **(Fig. \ref{fig:pdac}C)**.

Combined with the line 255 deletion above, this preserves the forward-pointing signpost but moves it to the natural alphabetical position (after both A and B are introduced, before the C analysis).

### Cost/benefit

- **Word delta:** −3 words (preferred) or net 0 (alternative)
- **No new analyses, no bib changes**
- **Reader benefit:** removes a strict-reading panel-order issue that PNAS reviewers can flag
- **Risk:** zero — purely mechanical fix; no scientific content changes

### Application checklist

1. Choose preferred (drop) or alternative (move) variant.
2. Apply the chosen edit at the specified line in `paper/04_results_REVISED.Rmd`.
3. Re-render with `make paper` and visually confirm panel order on first standalone reference is now A → B → C → D.

---

## 13. Resolve the apparent "Classical → Basal-like" paradox in the W-matrix correspondence sentence

### Update (2026-05-05): paragraph-level preview re-recommends Option 3 over Option 2

Paragraph-level preview of the PDAC factor structure paragraph zone (Results lines 253–263) revealed a notation-style consistency issue: §13 Option 2 strips subtype labels in the first sentence (*"NMF's N1 mapped strongly to DeSurv's D3 and NMF's N3 mapped primarily to DeSurv's D2"*), but the rest of the paragraph reverts to functional descriptors (*"DeSurv's most prognostic factor (D1)"*, *"NMF's exocrine-associated factor (N2)"*, etc.). The reader has to mentally translate bare N1/D3/N3/D2 in the first sentence, then read functional labels in the rest of the paragraph — a notation switch within the same paragraph.

**Updated recommendation: Option 3 instead of Option 2.** Option 3 keeps subtype labels but explicitly explains the apparent paradox, which is consistent with the rest of the paragraph's labeling style:

> "Two clear correspondences emerged: NMF's tumor-identity factor (N1, top-correlated with Classical programs) mapped strongly to DeSurv's tumor factor (D3, top-correlated with Basal-like programs), reflecting shared pan-tumor gene loadings rather than subtype concordance, and NMF's microenvironmental factor (N3) mapped primarily to DeSurv's D2."

Net delta vs current sentence: ~+15 words. The original §13 entry rated Option 3 as "heavier" and "only worth doing if a co-author has flagged the paradox." Paragraph-level preview shows Option 2's notation switch is its own friction; Option 3's slightly more words buy paragraph-style consistency.

The original Options 1 and 2 below remain as wording references; Option 3 is now recommended.

### Why this matters

The current sentence at `paper/04_results_REVISED.Rmd:261` describing the pairwise W-matrix correspondences reads:

> "Two clear correspondences emerged: NMF's Classical tumor factor (N1) mapped strongly to DeSurv's Basal-like tumor factor (D3), and NMF's microenvironmental factor (N3) mapped primarily to DeSurv's D2."

Verbally, this juxtaposes "Classical" and "Basal-like" — typically described as **opposing PDAC tumor subtypes** (Collisson, Moffitt, Bailey) — as if they map together. A reader's first parse is "wait, those are opposite subtypes — how do they map?" This is a clarity issue, not a scientific one: the empirical finding (N1 and D3 have high Spearman correlation in W-matrix gene loadings) is real, but the sentence puts the apparently-contradictory subtype labels front-and-center where they invite confusion rather than the loading-overlap finding the figure actually shows.

### Why the finding is real and not a labeling error

Both factor labelings (N1 = Classical, D3 = Basal-like) were assigned based on each factor's *top correlations* with established gene programs (lines 255 and 257). But each factor's full $W$ column also carries loadings on many *shared* tumor-identity genes — pan-tumor markers that load high in any factor capturing tumor cells, regardless of subtype. The Fig. 3D Spearman correlation reflects total gene-space overlap across all genes, which can be high between two factors whose top subtype-specific correlates differ, because both factors share substantial pan-tumor loadings. So N1 and D3 share many high-loading genes (tumor-identity background) even though their *top* correlated programs differ (Classical for N1, Basal-like for D3). The label-vs-loading-correlation distinction is what the current sentence obscures.

### Three drafts

**Option 1 — Drop the subtype labels in this sentence (preferred, ~−10 words):**

> Two clear correspondences emerged in $W$-matrix gene loadings: N1 mapped strongly to D3, and N3 mapped primarily to D2 (Fig. \ref{fig:pdac}D).

The factor identities (Classical for N1, Basal-like for D3, microenvironmental for N3) are already established two paragraphs earlier. This sentence is about the *correspondence* between factorizations, not the subtype identity, so the labels add confusion without information.

**Option 2 — Name what the correlation reflects (~+8 words):**

> Two clear correspondences emerged: NMF's N1 mapped strongly to DeSurv's D3 and NMF's N3 mapped primarily to DeSurv's D2 — reflecting shared tumor-identity and microenvironmental gene loadings across the two factorizations, rather than subtype-specific concordance.

More explicit about what the Spearman correlation captures (gene-space overlap including shared pan-tumor loadings) without invoking the apparent label paradox directly.

**Option 3 — Acknowledge the apparent tension and explain (~+15 words, most defensive):**

> Two clear correspondences emerged: NMF's tumor-identity factor (N1, top-correlated with Classical programs) mapped strongly to DeSurv's tumor factor (D3, top-correlated with Basal-like programs), reflecting shared pan-tumor gene loadings rather than subtype concordance, and NMF's microenvironmental factor (N3) mapped primarily to DeSurv's D2.

Most explicit; best if a reviewer or co-author has already flagged the apparent contradiction. Risk: actively draws attention to the paradox by naming both subtype labels in the same sentence.

### Recommended choice

**Option 1** is preferred for most contexts. The subtype labels are already in the reader's mind from the prior paragraphs (lines 255 and 257), the figure is right there for visual confirmation, and the Spearman-correlation finding is about gene-loading overlap regardless of subtype label. Stripping the labels in this one sentence removes the paradox without losing information.

**Option 2** is the right call if you want to actively name what the correlation reflects (gene-loading overlap including shared pan-tumor signal) without invoking the contradictory subtype labels.

**Option 3** is the right call only if the apparent contradiction has already been flagged by a co-author or reviewer and a direct rebuttal is needed.

### Cost/benefit

- **Word delta:** −10 (Option 1), +8 (Option 2), or +15 (Option 3).
- **No new analyses, no bib changes.** The Fig. 3D data and factor identities are unchanged; only the wording of the correspondence sentence changes.
- **Reader benefit:** removes a verbally confusing juxtaposition that triggers a double-take and forces the reader to reconcile the apparent paradox before continuing.
- **Risk:** very low. The subsequent sentences (about D1 lacking an NMF counterpart, N2's exocrine factor's weak D3 correlation) are unaffected and continue to use subtype labels appropriately in their own context.

### Coordination with rest of paragraph

After applying any of the three options, the rest of the paragraph at line 261 continues to use subtype labels in places where they are doing real work (e.g., "Because D1 couples Classical tumor and iCAF-associated stromal signatures into a single survival-aligned program..."). These uses are *interpretively meaningful* — the Classical+iCAF coupling is the paper's central biological observation — and should be left alone. The §13 fix is targeted at one specific sentence where the labels create confusion rather than clarity.

### Application checklist

1. Choose Option 1, 2, or 3 based on co-author input on whether the label paradox has been noticed.
2. Apply the chosen draft to `paper/04_results_REVISED.Rmd:261`, replacing the existing first sentence-and-a-half ("Two clear correspondences emerged: NMF's Classical tumor factor..." through "...DeSurv's D2.").
3. Verify the rest of the paragraph still flows cleanly into the next sentence about D1's lack of NMF counterpart.
4. Re-render with `make paper` and visually confirm the paragraph reads naturally with the revision.

---

# Cluster J — SI consistency

Singleton — three internal inconsistency fixes in SI Section 20 (the "independent of known classifiers" claim that contradicts Table S6; the "To quantify this directly" lead-in mismatch; the per-factor "no prognostically relevant structure" claim missing a Type III qualifier).

**Coordination with Cluster C:** §29a fixes SI Section 20's framing for apply-time consistency with §10's main-text draft. §10's actual draft cites Section 16 + Fig. S10 (not Section 20), so no §10 *draft* revision is needed; the constraint is sequencing only (§29a applies to SI before §10 applies to main text).

**Coordination with Cluster G:** §29a's revised wording explicitly invokes Fig. S7 (concordance), aligning with §28's main-text addition.

**Additional candidate sites identified during cluster-expansion analysis:**
- SI Section 10 (formal $\Delta\ell$ definition) could include a brief Type III caveat sentence noting per-factor metric may underestimate joint contribution.
- SI line 1582 wording ("continues to distribute survival contribution negligibly across all factors") could use the same Type III qualifier as §29c. Worth extending §29c.

---

## 29. Fix three internal inconsistencies in SI Section 20 ("Standard NMF factor structure at independently selected rank")

### Why this matters

SI Section 20 (`paper/si_appendix.Rmd:1425–1487`) is the SI's substantive treatment of NMF at independently selected ranks ($k = 5$ elbow, $k = 7$ BO). Its framing — "DeSurv's advantage thus lies not in maximal concordance but in recovering a compact, interpretable factorization" — was identified in §10 as a load-bearing piece of the contribution argument that the main text should mirror. Reviewing the section in detail reveals three specific issues that conflict with the data the same section reports:

1. **"Prognostic content independent of known classifiers" contradicts Table S6.** The section's closing sentence (line 1485) claims DeSurv's advantage is "prognostic content independent of known classifiers." But Table S6 in the same section shows DeSurv's adjusted HR drops 1.50 → 1.33 (~11% drop) after PurIST/DeCAF adjustment, while NMF $k = 7$ drops only 1.48 → 1.47 (essentially flat). So DeSurv is **more** dependent on (concordant with) known classifiers than NMF $k = 7$, not less. The §28 framing established the right interpretation: DeSurv's HR drop is biological *concordance* with established subtypes (good — convergent validation per Fig. S7), and NMF $k = 7$'s smaller drop reflects residual content from less interpretable factors. The current SI wording reverses this.

2. **"To quantify this directly" lead-in (line 1487) doesn't match what Table S6 quantifies.** Para 4 ends with the "DeSurv's advantage lies in ... independent of known classifiers" claim. Para 5 immediately says "To quantify this directly, Table S6 reports..." But Table S6 shows BOTH methods retain adjusted significance, with NMF $k = 7$ retaining slightly more independence than DeSurv. The lead-in sets up Table S6 as substantiating DeSurv's advantage; the actual table shows comparable adjusted retention with NMF $k = 7$ marginally ahead on this metric. The mismatch between lead-in expectation and table content is a logical gap.

3. **"Additional rank adds no prognostically relevant structure" at $k = 5$ (line 1429) needs a Type III qualifier.** This claim collides with the same section's Table S5, which shows NMF $k = 5$ substantially outperforms NMF $k = 3$ at per-cohort C-index across all five validation cohorts. Both can be true simultaneously: per-factor Type III $\Delta\ell$ negligible at $k = 5$ (Fig. S10 supports this) AND joint $k = 5$ model validates externally (Table S5 supports this). But the SI doesn't currently distinguish these levels, leaving a reader who notices the apparent contradiction without resolution.

A methods-aware reviewer comparing Table S6 to the section's claims, or comparing Table S5 to the "no prognostically relevant structure" claim at $k = 5$, will spot these gaps.

### Three proposed fixes

#### 29a. Fix the "independent of known classifiers" framing (line 1485)

**Current:**

> "DeSurv's advantage thus lies not in maximal concordance but in recovering a compact, interpretable factorization whose prognostic content is independent of known classifiers."

**Proposed (revised 2026-05-05 post-Amber-input):**

> "DeSurv's advantage thus combines biological concordance and methodological complementarity: its risk groups identify the same poor-prognosis subtypes as established classifiers (Basal-like, proCAF; Fig. S7), while its linear predictor retains significant prognostic value after adjustment for those classifiers (Table S6)."

**Why this revision (post-Amber-input):** Amber agreed with the direction (drop "independent" → "retains significant after adjustment") but flagged that the prior wording's *"...while concordantly capturing the established Basal-like and proCAF poor-prognosis biology"* was unclear, asking *"Are you trying to say that we have prognostic relevance beyond the classifiers while also capturing classifier information?"* The revised sentence explicitly states the dual-property point — biological concordance + methodological complementarity — rather than relying on "while concordantly capturing" which Amber found unclear. It also distinguishes the two metrics: dichotomized risk groups (Fig S7) vs. linear predictor (Table S6).

**Earlier proposed wording (now superseded):** *"DeSurv's advantage thus lies not in maximal concordance but in recovering a compact, interpretable factorization that retains significant prognostic value after adjustment for known classifiers (Table S6) while concordantly capturing the established Basal-like and proCAF poor-prognosis biology (Fig. S7)."*

Net delta: ~+12 words vs current SI text. Coordinates with §28 (which uses the same dual-property framing in main text).

#### 29b. Fix the "To quantify this directly" lead-in (line 1487)

**Current:**

> "To quantify this directly, Table S6 reports pooled validation hazard ratios (per SD of the linear predictor) for DeSurv at $k = 3$ and the BO-selected NMF model at $k = 7$ ($\alpha = 0$), both unadjusted and after adjustment for PurIST and DeCAF molecular classifiers. Both methods retain significant independent prognostic value after adjustment."

**Proposed:**

> "Table S6 reports pooled validation hazard ratios (per SD of the linear predictor) for DeSurv at $k = 3$ and the BO-selected NMF model at $k = 7$ ($\alpha = 0$), both unadjusted and after adjustment for PurIST and DeCAF molecular classifiers. Both methods retain significant prognostic value after adjustment, with DeSurv's larger HR drop reflecting biological concordance with the established subtypes (per Fig. S7) and NMF $k = 7$'s smaller drop reflecting residual content from its additional, less interpretable factors."

Drops "To quantify this directly" (which promised something Table S6 doesn't deliver). Engages with the actual finding (DeSurv drops more, NMF $k = 7$ drops less) and explains both via the parsimony + concordance framing established in §28. Net delta: ~+30 words.

#### 29c. Add Type III qualifier to "no prognostically relevant structure" at $k = 5$ (line 1429)

**Current:**

> "These biologically redundant factors explain why supervised BO selects $k = 3$ rather than $k = 5$: the additional rank adds no prognostically relevant structure."

**Proposed:**

> "These biologically redundant factors explain why supervised BO selects $k = 3$: the additional factors at $k = 5$ contribute no independent prognostic content beyond the three core programs (Fig. S10), even though the overall higher-rank model still validates externally on per-cohort C-index (Table S5)."

Distinguishes per-factor Type III contribution (Fig. S10 — additional factors negligible) from joint-model performance (Table S5 — $k = 5$ validates well). Resolves the apparent contradiction without changing the substantive claim. Net delta: ~+15 words.

### Optional minor polish (non-required)

#### "Maximal concordance" jargon (line 1485)

"Maximal concordance" is conventional survival-analysis terminology but slightly jargony for a general SI reader. Could simplify to *"DeSurv's advantage thus lies not in absolute discrimination performance but in..."* — defensible but stylistic. Skip if you prefer the conventional term.

#### Sharpen the "uses survival information to select hyperparameters via BO" sentence (line 1487)

**Current:**

> "The NMF $k = 7$ model uses survival information to select hyperparameters via BO but not to guide the factorization itself ($\alpha = 0$); DeSurv uses survival at both levels, concentrating the prognostic signal into three factors rather than seven."

**Proposed (slightly sharper, makes the borrowed-infrastructure point explicit):**

> "The NMF $k = 7$ model selects its rank via DeSurv's BO infrastructure (with $\alpha = 0$, no factorization-level supervision); DeSurv uses survival at both levels, concentrating the prognostic signal into three factors rather than seven."

Net 0 words. Mirrors §20a's framing in main text (NMF reaches $k = 7$ only when DeSurv's tuning framework is applied).

### Recommendation

**Apply 29a + 29b together as a coordinated fix** — they're intertwined (29b's "To quantify this directly" stops making sense once 29a's "independent of known classifiers" is corrected). Apply 29c as a smaller separate fix.

**Optional polish (minor sharpening + jargon simplification)** can be applied per stylistic preference; not strictly required.

### Coordination with main-text edits

| Fix | Coordinates with |
|---|---|
| 29a (drop "independent of known classifiers") | §28 (convergent validation in main text); §10 SI-aligned Path A (sequencing only — see correction below) |
| 29b (Table S6 lead-in) | §10 (which doesn't engage with Table S6 but should remain consistent); §15 (matched-rank framing) |
| 29c (Type III qualifier at $k = 5$) | §12 (framing-calibration: "weak" not "no" prognostic content); main-text Results line 259 if Type III caveat is added there too |

**§10 sequencing constraint (corrected 2026-05-05):** Earlier wording in this entry incorrectly claimed §10's recommended draft cites Section 20's "independent of known classifiers" framing. **It does not.** §10's actual SI-aligned Path A draft cites Section 16 ("convergent evidence") and Fig. S10 ("no new factor-level prognostic contribution") — the "independent of known classifiers" phrasing exists only in SI Section 20 itself, which §29a is correcting. **No §10 *draft* revision is needed post-§29a.** What stands is the **sequencing constraint at apply-time**: §29a must land in the SI before §10 is applied to main text, so SI and main text are consistent at apply-time. The application-priority structure (Wave 1 §29 before Wave 4 §10) handles this automatically.

### Cost/benefit

- **Word delta:** +10 (29a) + 30 (29b) + 15 (29c) = ~+55 words across SI Section 20.
- **No bib changes, no new analyses.** All claims trace to existing tables/figures.
- **Reader benefit:** removes three internal inconsistencies that any methods-aware reviewer comparing the section's claims to its tables will catch. Aligns the SI's framing with Table S6's actual values and with the §28 complementarity framing in main text.
- **Risk:** very low — substantive content unchanged; only the framing of three specific claims is corrected to match the underlying data.

### Application checklist

1. Apply 29a + 29b together to `paper/si_appendix.Rmd` (lines 1485 and 1487 respectively).
2. Apply 29c to `paper/si_appendix.Rmd:1429`.
3. After applying, re-read the full Section 20 to confirm the paragraphs flow coherently with the updated framing.
4. Update §10's recommended SI-aligned Path A draft (in this parking-lot file) to reflect 29a's corrected SI wording.
5. Consider applying optional minor polish (jargon simplification + borrowed-infrastructure sharpening) per stylistic preference.
6. Re-render with `make paper` and confirm SI Section 20 still reads cleanly with the revisions.
