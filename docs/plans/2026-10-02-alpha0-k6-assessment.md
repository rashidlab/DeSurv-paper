# alpha = 0 BO selected k = 6 (was k = 7): assessment for the review gate (2026-10-02)

Status: provisional. Written while the main BO (job 3148775) was still running; everything below assumes
DeSurv's own k stays at 3. Revisit at the Step 4 owner gate and again after Step 5 numbers exist.

## 1. The result

Stage `alpha0`, run `out/step4-dd54827-r2/alpha0` (BO 3148776 at dd54827; completion job 3456329 at 5261f60 built
only `tar_data_filtered_alpha0`; cached BO/finalist md5s unchanged). Paired 1-SE rule over 3 partitions x 5 folds:

| k | CV C | gap to best (SE) | within 1 SE |
|---|---|---|---|
| 2 | 0.500 | 0.136 (0.015) | no |
| 3 | 0.526 | 0.109 (0.018) | no |
| 4 | 0.500 | 0.136 (0.015) | no |
| 5 | 0.630 | 0.006 (0.004) | no |
| **6** | **0.636** | best | yes (selected) |
| 7 | 0.635 | 0.001 (0.005) | yes (runner-up) |
| 8-10 | 0.630-0.633 | <= 0.006 | yes |
| 11-12 | 0.628 | 0.007 | no |

Selected params: k 6, lambda 0.0563, nu 0.909, ntop 195. k = 6 and 7 are indistinguishable (0.001, SE 0.005).
Caveat: at k = 2 and 4 every finalist had a heavy Cox penalty (C exactly 0.500), so those ranks are
under-searched, not shown to be signal-free. k = 3 scored 0.526 at moderate penalties; the dedicated a0k3 BO
searches k = 3 properly.

## 2. Where the alpha = 0 rank-optimized control is used (not only a supplementary figure)

- Abstract (paper.Rmd L5): "a rank-optimized seven-factor control carried comparable prognostic information".
- Main Table hr-adjusted row (Results L535) and Results L601 (rank set at k = 7; coverage 16 -> 32; multiplicity;
  "Seven is where the one-standard-error rule landed").
- Discussion L12 and L21 ("seven-factor").
- SI: section L2047-2058 (factor structure + `si_fig_nmf_k7_heatmap`), L2076-2088 (coverage anchors), L2128/L2275
  (ntop), C-index table L2144, L2296 (k = 7 control exceeds DeSurv C in 5/5 datasets).
- NOT this control: `bo_best_k7` in code/29 and the Methods/SI final-config CV block. That is the main DeSurv BO's
  overall best (k = 7 in May); it changes with the main BO.

## 3. Mechanical changes if k = 6 stands

- Recomputed automatically: code/04 consensus fit (100 seeds) at k = 6, code/05 validation HR/C, table rows,
  inline coverage/multiplicity, the SI heatmap (6 rows).
- Breaks or typed: code/24:91 `stopifnot(ncol(W) == 7)`; anchor id "nmf_k7_unsupervised" (code/24:197, Results
  L592-594, SI L2076-2078); figure filename `si_fig_nmf_k7_heatmap` (code/09b:224, SI L2058); ~12 typed "7"/
  "seven" sites listed in section 2. Make them read k from `tar_params_best_alpha0`.
- Claims inventory: no re-freeze needed (tests reference the control object, held fixed across reruns).
  Q-R19c ("comparable"), R19d (5/5 datasets) and the coverage claim are scored on the new numbers.

## 4. Interpretation if k = 6 replaces k = 7 in its current role

Unchanged structurally. "Rank-optimized control carried comparable prognostic information" probably holds (CV C
0.636 vs 0.635 for k = 7), but it is a frozen binary claim scored on new validation numbers. "Exceeds DeSurv in
all five datasets" is fragile. Coverage at k = 6 likely a little below k = 7; "some programs represented by more
than one factor" may weaken. The "where the rule landed" hedge becomes more honest if the text says 6 and 7 tie.

## 5. Why the main comparisons use alpha = 0 at k = 3, not the rank-optimized control

Design intent, stated in SI L2050 and Discussion L12: holding rank (and BO budget) fixed isolates the effect of
survival supervision on the factors. A D1 that is absent from k = 3 NMF cannot be blamed on DeSurv having more
or fewer factors. Comparing DeSurv (k = 3) against NMF at k = 6/7 confounds supervision with rank: any program
difference could come from the extra factors. The rank-optimized control answers a different question (how much
prognostic information unsupervised NMF reaches at its own best rank) and was therefore used only for prognosis
and reference-program coverage.

Gap: the paper never tests whether D1 is present among the rank-optimized control's factors. "Did not recover D1"
is established only at matched rank. A reviewer can ask: does unsupervised NMF find D1 once it is allowed its own
best rank?

## 6. Swapping in the k = 6 control as a main-text comparator

Do not replace the matched-rank control; it is the only clean test of supervision. Add the k = 6 control as a
second D1 test, using the same statistic already used at k = 3 (R^2 of D1 on the NMF basis, and max |r| with any
single factor). Cheap: both fits exist after Step 5. Two outcomes:

- D1 not recovered at k = 6 (low R^2, no factor with high |r|): the claim strengthens to "survival supervision
  finds a program that unsupervised NMF misses at matched rank and at its own selected rank". Abstract and
  Discussion L12 can say so.
- D1 recovered at k = 6 (one factor or a combination closely matches D1): D1 is not supervision-specific; it is
  a program unsupervised NMF reaches with more factors. The claim becomes about representation: supervision
  isolates the prognostic program at k = 3 where NMF needs about six. "Did not recover D1" must then be scoped
  explicitly to matched rank, and Discussion L12's "did not appear to reflect rank" needs revision.

Process: this is a new fit-dependent claim. Add it to the claims inventory as an owner-approved amendment BEFORE
any Step 6 rerun output exists, or report it as a descriptive quantity only (D1 R^2 on the k = 6 basis and max
|r| with any k = 6 factor, primary value plus rerun range). A binary form needs a comparator the owner picks
(e.g. "D1 R^2 on the k = 6 basis is below D2's and D3's"); no numeric cutoff, per the 2026-09-30 ruling.
