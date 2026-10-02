# What-if: is D1 present in an unsupervised fit at its own selected rank? (2026-10-02)

Isolated analysis. Not in the pipeline, the claims inventory, or the paper. Uses the PUBLISHED DeSurv fit (D1-D3);
repeat on the Step 5 fit when it exists. `Rscript analysis/whatif-d1-own-rank/run.R` (DeSurv >= 1.2.0).
k = 6 alpha = 0 fit made here with the reported recipe (seeds 1..100; params from the step4-dd54827-r2 alpha0 BO;
3 min, 8 cores); it is a preview of Step 5's alpha0 fit. Sanity: the k = 3 rows reproduce the paper's inline values
(D1 projection R^2 0.091; max |Spearman| 0.357 / 0.873 / 0.845).

| basis | D1 gene R^2 | D2 | D3 | D1 best single-factor rho | D1 score R^2 | D2 | D3 | D1 best top-gene Jaccard |
|---|---|---|---|---|---|---|---|---|
| alpha0 k = 6 (new BO) | 0.152 | 0.827 | 0.776 | 0.397 | 0.736 | 0.998 | 0.992 | 0.166 |
| alpha0 k = 7 (paper control) | 0.151 | 0.825 | 0.778 | 0.388 | 0.734 | 0.998 | 0.992 | 0.159 |
| std NMF k = 3 (paper Fig 2) | 0.091 | 0.814 | 0.755 | 0.357 | 0.295 | 0.997 | 0.980 | 0.093 |
| alpha0 k = 3 matched budget | 0.093 | 0.824 | 0.745 | 0.245 | 0.278 | 0.997 | 0.979 | 0.098 |

Gene R^2 = OLS R^2 of the DeSurv loading vector on the basis W (paper's projection_r2). Score R^2 = R^2 of the
DeSurv H row (training patients) on the basis H rows. Scale check: total factor mass per patient is nearly
constant (CV 0.007); R^2 with H as within-patient proportions is 0.733 (k = 6), 0.732 (k = 7), 0.278 (NMF k = 3),
so the score result is not a scale artifact. D1 score on k = 6 factors: single-factor r +0.59, +0.56, -0.50, -0.51
(F1, F2, F5, F6); no single factor carries it.

Reading:
1. k = 6 and k = 7 give the same answer; the rank choice does not matter for this question.
2. Gene program: D1 is still not recovered at the unsupervised control's own rank. Its loadings are 15% explained
   by the six or seven factors (vs 9% at k = 3; D2 and D3 83% and 78%), best single-factor rho 0.40, top-gene
   overlap 0.17.
3. Patient scores: at its own rank the unsupervised fit spans about 73% of the variation in D1's patient scores
   (vs about 29% at k = 3), through a signed combination of several factors, none dominant.
4. Together: unsupervised NMF given more factors captures most of the patient-level variation D1 summarizes, spread
   across several programs, but does not produce D1's gene program. Consistent with D1 as a supervised combination
   ("classical tumor + restCAF-like") that NMF keeps in separate factors. "Did not recover D1" holds at the program
   level beyond matched rank; a statement about patient-level information would need the score result's caveat.

## Follow-up: CAF states in the own-rank alpha = 0 factors (stroma.R, stroma_per_factor.csv)

Same reference lists and views as code/23_stroma_resolution.R (SCISSORS iCAF/myCAF = the paper's restCAF/proCAF).
Top-50 view, k = 6 (k = 7 is the same: its F5 and F4 match k = 6's F2 and F5):

- F5 = activated stroma / myCAF-like. Moffitt Activated 8/24 (hypergeometric p 5e-8), SCISSORS proCAF 5/16
  (3e-5), deCAF proCAF 3/4 (6e-5), Elyada myCAF 2/15 (0.05), DECODER Activated 39/108. Top genes COL11A1, ITGA11,
  ADAM12, CTHRC1, COL10A1, LRRC15, GREM1, FAP, LOX, PRRX1. A clean activated-CAF program.
- F2 = "iCAF-enriched" by list overlap (Elyada iCAF 8/25, p 8e-8; SCISSORS restCAF 7/17, 7e-8; deCAF restCAF 3/8),
  but the hits are quiescent/adventitial fibroblast genes (PI16, DPT, TNXB, ADH1B, GPX3, SFRP1, FBLN5, PLA2G2A,
  FOSB). The inflammatory iCAF genes present in the universe (CXCL1, CXCL2, CXCL12, LIF, C3, CFD, CCL2) are not in
  its top 50. Its top genes are a normal-tissue mix: lymphoid (CCL19, CCL21, CR2), islet/neuroendocrine (CHGA, SCG2,
  PCSK2), adipose (PLIN4, FABP4), plus PTGDS, NGFR.
- Full-gene view is weak for every fit (|rho| <= 0.14 for CAF-state lists; DECODER Activated 0.33 for F5), as
  code/23 already warns.
- DeSurv D1 (k = 3) top 50: 2 Elyada iCAF, 2 SCISSORS restCAF, 2 Elyada myCAF, 3 DECODER Normal: mixed.
  NMF k = 3: no CAF-state hits outside its activated factor.

Reading: the own-rank unsupervised fit separates activated stroma (myCAF-like) from a normal / quiescent stromal
compartment. It does not recover an inflammatory iCAF program; the iCAF-list overlap comes from quiescent
fibroblast genes shared between single-cell iCAF signatures and normal stroma. This agrees with the prior bulk
experience: activated vs normal stroma resolves, iCAF vs myCAF does not. Caveat: one fit, top-50 view, no
stability assessment.

## Competitor check: how D1 got its restCAF / iCAF label (published fit; read-only)

Mechanism (R/figure_plot_helpers.R make_gene_overlap_heatmap, Fig. 2a): take the union of the three factors' top-50
genes (150), keep those in any of the 45 reference lists (75 genes), and for each factor correlate its loadings over
those 75 genes with 0/1 list membership (Spearman); display rows with r > 0.2. SCISSORS iCAF/myCAF are renamed
restCAF/proCAF in the display.

What D1's label rests on:
- Elyada iCAF: r 0.24, FDR 0.19, carried by 2 genes (S100A4, CXCL14).
- "SCISSORS restCAF" (= SCISSORS iCAF): r 0.24, FDR 0.19, carried by 2 genes (LIF, CXCL14). CXCL14 is shared.
- Elyada myCAF: 2 genes in D1's top 50 too (CST1, HOPX), r 0.14, FDR 0.43.
- deCAF restCAF (the DeCAF classifier the "restCAF" name comes from): all 8 genes are in the 1,970-gene data and in
  the reference set, but none is in any DeSurv factor's top 50, so the list contributes nothing to Fig. 2a.
- Full-gene view (code/23): D1 |rho| <= 0.10 for every CAF-state list.

Points a competitor would make:
1. The label rests on 2-3 genes and is not significant after FDR (0.19); it is displayed because r clears 0.2.
2. D1 carries myCAF genes as often as iCAF genes in its top 50.
3. "restCAF" is a renamed single-cell iCAF list, not the DeCAF restCAF genes, which are absent from D1.
4. "restCAF and iCAF corresponded to a DeSurv factor but to no standard NMF factor" (Results L440) holds only at
   matched rank. Fig. 2a's own function (make_gene_overlap_heatmap, same reference set, same r > 0.2 display rule)
   applied to each fit:

   | fit | Elyada iCAF: best factor r, FDR, list genes in support | SCISSORS iCAF ("restCAF"): r, FDR, genes | deCAF restCAF: r, FDR |
   |---|---|---|---|
   | DeSurv k = 3 (D1) | 0.24, 0.19, 2 (S100A4, CXCL14) | 0.24, 0.19, 2 (LIF, CXCL14) | none in top 50 |
   | NMF k = 3 | none in top 50 | none in top 50 | none in top 50 |
   | alpha0 k = 6 (F2) | 0.27, 1.2e-4, 8 | 0.24, 9.2e-4, 8 | 0.15, 0.057, 3 |
   | alpha0 k = 7 (F5) | 0.27, 6.2e-5, 9 | 0.25, 2.1e-4, 9 | 0.17, 0.016, 4 |

   Caveat: the correlation universe is the union of all factors' top-50 genes and grows with k (75 genes at k = 3,
   240 at k = 6, 251 at k = 7), so FDRs are smaller at k = 6/7 partly for that reason alone; compare r, not FDR.
   The correlation is about the same; the difference is support and significance (2 genes, FDR 0.19 for D1 vs 8-9
   genes, FDR 1e-4 for the own-rank control). At k = 7 the control also shows myCAF (r 0.20, FDR 0.004) and
   SCISSORS myCAF (0.23, 0.0013) on a separate factor (F4). By the paper's own metric, the own-rank unsupervised
   fit matches the iCAF/restCAF lists more strongly than D1 does, through quiescent-fibroblast genes (section above).
5. The text already hedges ("among the displayed top-ranked genes"; "not identifications of CAF cell states"), but
   the abstract's "classical tumor + restCAF-like" label carries none of it.
Claims-inventory note: the frozen label test (best reference overlap, per-program argmax) would pass on 2 genes; it
tests stability across reruns, not strength. Repeat this check on the Step 5 fit.

## Fig. 2a-style heatmaps, side by side (heatmap_*.png; programs with r > 0.2, * = FDR < 0.05)

- DeSurv k = 3: D1 classical tumor (DECODER 0.40*, Puleo 0.32) + iCAF/restCAF (0.24, n.s.); D2 immune (0.42*) +
  activated stroma (0.37*) + proCAF (0.23); D3 basal-like (0.34-0.39*). 16 programs shown.
- NMF k = 3: N1 classical AND basal-like mixed; N2 exocrine (0.81*); N3 immune + activated stroma. 15 shown.
- alpha0 k = 6 and k = 7 (near identical): a clean classical tumor factor (DECODER 0.65-0.67*, 9 classical lists), a
  clean basal-like factor (DECODER 0.58-0.61*), exocrine, activated stroma / ECM (DECODER 0.63-0.64*, Maurer ECM
  0.54-0.55*, proCAF 0.23*; k = 7 adds Elyada myCAF 0.20*), and an immune + endocrine + iCAF/restCAF factor
  (Maurer Immune 0.34-0.40*, DECODER Endocrine 0.31-0.37*, iCAF 0.27*, restCAF 0.24-0.25*). k = 6 has one factor
  with nothing above 0.2; k = 7 splits a second, weaker ECM factor and a weak immune factor. 27 / 32 shown.

Reading: the own-rank unsupervised solution is cleaner and more interpretable than either k = 3 fit. It separates
classical from basal-like (NMF k = 3 mixes them), separates activated stroma from the immune/quiescent compartment,
and its classical factor matches the classical lists more strongly than D1 (0.67 vs 0.40). DeSurv's D1 is the
combination of two things this solution keeps apart: classical tumor (its F1) and the quiescent/immune stromal
factor (its F2), matching the score regression above (D1 score r +0.59 with F1, +0.56 with F2). What k = 3 DeSurv
offers over it is a compact supervised combination, not resolution of programs the unsupervised fit misses.
