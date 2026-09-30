# Claims inventory, rev 3 Step 2 (2026-09-30)

Agent B deliverable for plan rev 3, section 4 and the owner clarification in section 8. Frozen before any
rerun output exists. Line numbers are at HEAD 88cf0f7 (the concurrent Agent A edits in the worktree are
line-preserving; checked with `git diff -U0`). "Current value" is read from the tracked `results/*.rds`
(before side, fb731e7 PDFs) in a throwaway R session that evaluated the purled Rmd chunks; nothing in
`results/` was written.

Scope: every substantive claim in the abstract (`paper/paper.Rmd` L5), `paper/04_results_REVISED.Rmd`,
`paper/05_discussion_REVISED.Rmd` and the Results-type SI sections of `paper/si_appendix.Rmd` whose truth
depends on the fitted DeSurv solution (W, H, beta, factor labels, top genes, validation scores). Claims about
data, cohorts, simulations with fixed truth, the method's definition, and the fixed comparators alone are
excluded (borderline cases listed at the end). Order is manuscript position. Nothing is ranked.

Location codes: MAIN = main text, SI = Supplementary Information, ABS = also in the abstract.

---

## 0. Shared definitions (apply to every test below)

**0.1 Alignment (A).** For rerun j, compute the k x k matrix of Pearson correlations between the columns of
W_j and W_primary over all analyzed genes, and solve the assignment problem (Hungarian, `clue::solve_LSAP`,
maximize) on that matrix. Permute W_j columns, H_j rows and beta_j entries into primary order. After
permutation, "D1/D2/D3 in rerun j" means the factor matched to primary D1/D2/D3. Record the matched r and
top-ntop Jaccard per factor (plan section 4 reporting). [OWNER: plan says "Hungarian on correlation"; this
file fixes Pearson on full W columns. Spearman is the manuscript's loading-correspondence statistic; say so if
Spearman is preferred.]

Why the permutation is required: every downstream script identifies programs by column position, not by
label. Examples: `04_results_REVISED.Rmd` L237 (`D1 = Z[, 1]`), L147-148 (`.corrC[1, ]`), L466 guard
`beta[2] == 0`; `code/14_d1_variance_partition.R:40`; `code/11_treated_cohort_analysis.R:57`;
`code/15_supervised_recovery.R:42`; `code/18_clinical_adjustment.R:34,44`;
`code/19_desurv_vs_supervised_tuned.R:82-128`; and the heatmap labels are hard-coded by position at
`code/09a_figures.R:161-163`. Unpermuted rerun outputs would score the wrong factor.

**0.2 Per-rerun recomputation.** A claim is scored from rerun j's own outputs produced by the same code that
produced the current value, applied to the permuted fit, into a rerun-specific results directory: fit (code/04
recipe) -> `val_latent_desurv_tcgacptac` (code/05) -> the named downstream script. Fixed comparators
(`fit_std_desurvk`, `tar_fit_desurv_a0k3`, `tar_fit_desurv_alpha0`, their `val_latent_*`, the supervised
comparators' own fits) are NOT rerun; they are read from the primary build.

**0.3 Significance.** "Significant" / "associated" = two-sided Wald P < 0.05 from the named Cox model,
equivalently the 95% CI excludes 1. "Not significant" = P >= 0.05.

**0.4 Label-assignment test (L), PROPOSED, for owner review.** Labels are anchored to reference programs, not
to the primary fit, so the same test can first be applied to the new primary fit and then to each rerun.
Statistic: the factor x reference-program Spearman matrix from
`R/figure_plot_helpers.R::make_gene_overlap_heatmap(fit, get_top_genes(W, ntop)$top_genes, top_genes,
return_matrix = TRUE)$cor_mat` (top 50 genes per factor, the Fig. 2a statistic; reference lists from
`cmbSubtypes_formatted.RData`). "Argmax" = the factor with the largest rho for that program among the k
factors. A factor carries its label when every condition holds:

| Label | Conditions | Current rho (D1 / D2 / D3) |
|---|---|---|
| L-D1 "classical tumor" | D1 is argmax for DECODER_ClassicalTumor and rho > 0.2 | 0.40 / -0.40 / -0.37 PASS |
| L-D1 "restCAF-like" | D1 is argmax for SCISSORS_restCAF and rho > 0.2 (fallback: Elyada_iCAF, same rule) | restCAF 0.24 / -0.13 / -0.09; iCAF 0.24 / -0.10 / -0.01 PASS (margin 0.04) |
| L-D2 "proCAF-like" | D2 is argmax for SCISSORS_proCAF and rho > 0.2, AND argmax for DECODER_ActivatedStroma with rho > 0.2 | proCAF -0.19 / 0.23 / 0.04; ActStroma -0.06 / 0.37 / 0.06 PASS (proCAF margin 0.03) |
| L-D3 "basal-like" | D3 is argmax for DECODER_BasalTumor and Puleo_Pure Basal-like, each rho > 0.2 | Basal 0.34 (D3) ; Pure Basal 0.39 (D3) PASS |

A rerun fails L for a factor when any condition fails. Any claim naming a label ("classical", "restCAF-like",
"proCAF-like", "basal-like") holds in rerun j only if L holds for that factor in rerun j AND the claim's own
test holds. Claims about D1/D2/D3 by number only (e.g. R7) use alignment A alone. Two of the four conditions
pass the 0.2 threshold by 0.03-0.04 at the current fit (see Findings). The reference list file is NOT in
this worktree: `code/09a_figures.R:147-151` falls back to `../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData`
(Amber's repo); the label test needs that file pinned (copy or checksum) before scoring. [OWNER]

---

## 1. Abstract (paper/paper.Rmd L5; typed, not computed)

**A1.** ABS, MAIN. "DeSurv recovered a three-program tumor--stroma architecture whose survival-aligned program
(D1) combined classical-tumor and restCAF-like stromal features."
Current: k = 3; L-D1 passes (table 0.4); D1 is argmax of Delta-loglik (R5).
Test: L-D1 (both rows) holds AND R5 holds AND R7b (D1 the only significant program) holds. All three required.

**A2.** ABS, MAIN. "the risk score remained prognostic (HR per pooled SD 1.42, 95% CI 1.26--1.59 ...)".
Current: 1.42 (1.26-1.59), P = 5.4e-9. Test: = R10.

**A3.** ABS, MAIN. "1.28, 95% CI 1.11--1.47 after PurIST and DeCAF adjustment".
Current: 1.28 (1.11-1.47), P = 8.3e-4. Test: = R19a.

**A4.** ABS, MAIN, Discussion L12. "A matched-rank, matched-budget unsupervised control did not recover D1".
Current: max |Spearman rho| of D1 with any a0k3 factor = 0.245. Test: = R9a.

**A5.** ABS, MAIN. "whereas a rank-optimized seven-factor control carried comparable prognostic information."
Current: DeSurv unadjusted 1.42 (1.26-1.59) vs k = 7 control 1.44 (1.30-1.59). Test: = R19c.

(The abstract's final sentence is simulation-only and excluded.)

---

## 2. Main-text Results (paper/04_results_REVISED.Rmd)

**R1.** MAIN, L426. "Individual restarts reached similar training concordance (median C-index 0.70, range
0.68--0.73) yet dispersed gene sets (median pairwise Jaccard ... 0.13 against 0.07 by chance)."
Object: `program_stability_stats` (code/20), from the seed-fit pool.
Current: C median 0.704 (0.678-0.732); median of per-factor median Jaccard 0.13 vs chance 0.074.
Test (per rerun, on its own seed bank): max - min training C over seed fits < 0.10 AND median pairwise
top-ntop Jaccard (factor-matched as code/20) < 0.30. [OWNER: thresholds quantify "similar" and "dispersed".]
Note: `program_stability_stats$pairwise` has 435 rows per factor (= 30 choose 2) while `meta$n_restarts` = 100;
confirm which restarts code/20 summarizes before reusing it per rerun.

**R2.** MAIN, L426; SI L1865-1873; Discussion L21. "Consensus fits from each half of 8 split-half partitions ...
agreed closely with the reported basis (median gene-loading correlations 0.97--0.99, median top-ranked-gene
Jaccard 0.73--0.88)." Object: `consensus_endpoint_stats` (code/21).
This is a procedure-level stability claim that the ten reruns measure directly. Test (per rerun j): all three
alignment-A matched Pearson r (rerun vs primary) >= 0.90. [OWNER: threshold; the rerun-vs-primary comparison
replaces the split-half statistic as the operational test; code/21 need not be rerun.]

**R3.** MAIN, L437. "DeSurv resolved three distinct programs (D1--D3): D1, enriched for classical tumor and
restCAF-like stromal signatures; D2, enriched for proCAF-like stroma; and D3, enriched for basal-like tumor
programs (Fig. 2a)." Current: all L conditions pass. Test: L-D1, L-D2, L-D3 all hold.

**R4.** MAIN, L437; Discussion L33. "In the COMPASS trial cohort, D1 scores correlated with GATA6 RNA in situ
hybridization ... (Spearman rho = 0.68, P < 0.001, n = 33)". Object: `treated_cohort_stats$gata6` (code/11).
Current: rho 0.679, P = 1.4e-5. Test: rho > 0 AND P < 0.05.

**R5.** MAIN, L443 (plan-required). "D1 accounted for only 18.8% of the full model's reconstruction gain
(Shapley shares), the smallest share among the three factors."
Object: `fig_variation_explained_tcgacptac$data` (code/09a L~100-143), accessor `.vp` at L114.
Current: D1 18.8, D2 25.2, D3 56.0. Test: D1's normalized Shapley share < D2's and < D3's.

**R6.** MAIN, L443; SI L1289-1294; Fig. 2c caption L451. "omitting D1 and refitting lowered the maximized
partial log-likelihood by 66.4, compared with 1.0 for D3 and 0.02 for D2."
Object: same, `.dl` at L113. Current: D1 66.37, D3 1.00, D2 0.02.
Test (R6a, carries A1 "survival-aligned"): D1 has the largest Delta-loglik of D1-D3.
Test (R6b, "little ... conditional survival contribution" for D3): Delta-loglik(D3) <= 0.1 x Delta-loglik(D1).
[OWNER: 0.1 quantifies "little".]
Test (R6c, SI L1293 "concentrates in DeSurv D1 with all NMF factors small"): Delta-loglik(D1) > max NMF
Delta-loglik (fixed, 1.16).

**R6d.** MAIN, L443. "D3 showed the opposite pattern, accounting for 56% of the reconstruction gain".
Test: D3 has the largest Shapley share of D1-D3 AND R6b holds.

**R7.** MAIN, L446. "The principal NMF tumor and microenvironmental factors corresponded closely to D3
(rho = 0.84) and D2 (rho = 0.87), whereas D1 correlated only weakly with every matched-rank NMF factor
(rho <= 0.36)". Accessors L146-148 (`cor(W_desurv, W_nmf, "spearman")`, max |rho| per DeSurv row).
Current: D3 0.845, D2 0.873, D1 0.357.
Test (R7a): max |rho|(D2, NMF) >= 0.70 AND max |rho|(D3, NMF) >= 0.70.
Test (R7b "weakly"): max |rho|(D1, NMF) < 0.50 AND < min of the D2 and D3 values. [OWNER: 0.70 / 0.50.]

**R8.** MAIN, L446; SI L1435, Table S projection-r2 caption L1430. "D1's gene weights were poorly
reconstructed from the three NMF factors together (R^2 = 0.09)" / SI "essentially unreconstructable".
Accessor `projection_r2(W_desurv, W_nmf)[1]` (L149; R/reconstruction_helpers.R). Current 0.091.
Test: R^2(D1 | NMF basis) < 0.25. [OWNER]

**R8b.** MAIN, L446; SI L1457. "DeSurv reconstructed the analysis matrix less completely than matched-rank NMF
(R^2 0.75 versus 0.81, a difference of 6.0 percentage points)" / SI "a real but bounded reconstruction cost".
Accessor `.r2c_*` L132-134 (`overall_recon_r2_centered`). Current 0.752 vs 0.812, gap 6.0 pp.
Test: 0 < (R^2_NMF - R^2_DeSurv) x 100 < 10. [OWNER: 10 pp bounds "bounded".]

**R8c.** MAIN, L446 (conclusion sentence). "survival supervision preserved most of the variance-dominant tumor
and stromal structure while resolving D1 as a distinct program not captured by matched-rank NMF."
Test: R7a AND R7b AND R8 all hold.

**R9.** MAIN, L449; SI L2284-2289. Matched-budget control "again recovered D2 (rho = 0.86) and D3 (rho = 0.75)
but similarly failed to recover D1 (rho <= 0.25)". Accessor `.corrF_*` L162-166
(`cor(W_desurv, W_a0k3, "spearman")`). Current: D2 0.861, D3 0.749, D1 0.245.
Test (R9a, carries A4 and Discussion L12): max |rho|(D1, a0k3) < 0.50 AND < min of the D2, D3 values.
Test (R9b): max |rho|(D2, a0k3) >= 0.70 AND max |rho|(D3, a0k3) >= 0.70. [OWNER: D3 passes 0.70 by 0.049.]
Test (R9c, SI L2285 "slightly lower than ... NMF::nmf"): max |rho|(D1, a0k3) < max |rho|(D1, NMF).

**R10.** MAIN, L466 (plan-required, ABS A2); Discussion L9 "transferred ... to four independent external
cohorts". "The linear predictor was associated with overall survival in the dataset-stratified analysis (HR per
pooled SD 1.42, 95% CI 1.26--1.59, P = 5.4e-09)". Accessor `pooled_cox_stats(pool_df)` L207-219, L301.
Current 1.42 (1.26-1.59). Test: 95% CI lower bound > 1.

**R11.** MAIN, L463 (plan-required). "D1 had protective point estimates in all five validation datasets".
Statistic: per-dataset Cox HR per within-dataset SD of the D1 score (Fig. 3a left, platform datasets, not
de-duplicated). Current: Dijk 0.77, Moffitt 0.83, PACA-AU array 0.69, PACA-AU seq 0.60, Puleo 0.73.
Test: HR < 1 in all 5 datasets.

**R12.** MAIN, L463 (plan-required, also A1). "D1 ... was the only individual DeSurv program with a
statistically significant combined association in the dataset-stratified analysis (HR 0.75, 95% CI
0.68--0.83, P < 0.001)". Accessor `.pf` L233-249 (de-duplicated, within-dataset standardized,
`strata(dataset)`). Current: D1 0.75 P = 1.6e-8; D2 1.02 P = 0.76; D3 1.09 P = 0.10.
Test: P(D1) < 0.05 AND P(D2) >= 0.05 AND P(D3) >= 0.05. (Direction: D1 HR < 1 as part of the claim.)

**R13.** MAIN, L463 (plan-required, "D3 per-cohort pattern"). "D3 had hazard ratios above 1 in four of the
five datasets (1.09, 95% CI 0.98--1.21, P = 0.10)". Statistic as R11 for D3.
Current: 1.27, 1.07, 0.95, 1.26, 1.04 (4 of 5 above 1). Test: HR > 1 in >= 4 of 5 datasets.
[OWNER: the sentence reads as "exactly four"; this test uses ">= 4". Say if exact.]

**R14.** MAIN, L463. "while D2 was near null overall (1.02, P = 0.76) and varied more across datasets".
Test (R14a): pooled D2 95% CI includes 1 (from `.pf`).
Test (R14b "varied more"): SD across the five datasets of per-dataset log HR (R11 statistic) for D2 > the
same SD for D1 and for D3. Current SDs: D2 0.128, D3 0.127, D1 0.122 (holds by 0.001). [OWNER: comparator of
"more" is unstated in the text; this test assumes "than D1 and D3".]

**R15.** MAIN, L466; SI L1052-1053, L2395; Discussion L18. "The trained DeSurv linear predictor combines D1 and
D3 ...; the elastic-net penalty set the D2 coefficient to exactly zero". Guard `beta[2] == 0` in the Rmd.
Current beta = (-0.00100, 0, 9.36e-5). Test: beta_D2 == 0 exactly AND beta_D1 != 0 AND beta_D3 != 0 (after
permutation A). Note the Rmd guard at L466 would break on an unpermuted rerun.

**R16.** MAIN, L466. "the estimates were similar with each dataset omitted in turn (HRs 1.34--1.45)".
Accessor `.loco` L284-292. Current 1.344-1.453. Test: every LODO point estimate lies inside the full pooled
95% CI of R10. [OWNER]

**R17.** MAIN, L466; SI L1621; Discussion L27. "Scaled Schoenfeld residuals indicated attenuation over follow-up
(chi^2 = 14.7, 1 df, P = 1.3e-4)". Object `validation_heterogeneity_ci$ph` (code/25).
Current chi^2 14.67, P = 1.3e-4, slope -0.92. Test: P < 0.05 AND `slope_time` < 0.

**R18.** MAIN, L469; SI L1640-1646. "Allowing cohort-specific slopes did not significantly improve fit (LR
chi^2 = 1.5, 3 df, P = 0.69) ... the cohort estimates were descriptively homogeneous (Q = 1.5, I^2 = 0%).
Discrimination across cohorts was modest (event-weighted C-index 0.61, 95% CI 0.58--0.65)."
Object `validation_heterogeneity_ci` (code/25). Current LRT P 0.69, I^2 0, C 0.613 (0.579-0.646).
Test (R18a): LRT P >= 0.05. Test (R18b): I^2 < 0.25. [OWNER] Test (R18c "modest"): pooled C < 0.70 AND its
95% CI lower bound > 0.50. [OWNER]

**R18d.** MAIN, L469; SI L1576-1599. Recalibrated IBS 0.159 vs 0.168 KM reference (IPA 5.4%); SI "modest
proper-score improvement over the null". Object `desurv_brier_stats` (code/17).
Test: 0 < IPA (pooled) < 0.10. [OWNER]

**R19.** MAIN, Table hr-adjusted (L478-547) and L598, L601. Objects: `val_latent_desurv_tcgacptac` with
`pooled_cox_unadj_adj` L484-502; controls fixed.
- R19a (ABS A3; Discussion L9 via R20): DeSurv adjusted (PurIST + DeCAF, strata dataset) HR 95% CI lower bound
  > 1. Current 1.28 (1.11-1.47).
- R19b (SI L2296 "Both the k = 7 control and the DeSurv linear predictor retain significant prognostic value
  after adjustment"): = R19a.
- R19c (ABS A5; L601 "carried comparable external prognostic information"; Discussion L12): DeSurv unadjusted HR
  point estimate lies inside the k = 7 control's unadjusted 95% CI (fixed: 1.30-1.59). Current 1.42 PASS.
  [OWNER: operationalizes "comparable".]
- R19d (SI L2296): "At k = 7, the rank-optimized control exceeds DeSurv's C-index in all five validation
  datasets." Objects `val_cindex_desurv_tcgacptac` vs fixed `val_cindex_desurv_alpha0_tcgacptac`.
  Current DeSurv 0.607/0.562/0.649/0.634/0.636 vs 0.631/0.578/0.669/0.702/0.651. Test: control > DeSurv in
  5 of 5.
- R19e (SI L2296, Discussion L12 "matched-rank ... did not recover" context): "additional unsupervised factors
  capture prognostically relevant variation missed at the matched rank" involves only fixed controls: excluded.

**R20.** MAIN, L566; Discussion L9 "D1 remained prognostic after adjustment for ... PurIST and DeCAF".
"binary PurIST and DeCAF classifications together explained 34% of the within-dataset variance in the continuous
D1 score, and DeCAF alone explained 19%. D1 nevertheless remained associated with survival after adjustment for
both classifiers (HR 0.82, 95% CI 0.73--0.93, P = 0.0014). Adding D1 ... improved fit (chi^2 = 10.1, P = 0.0015)
and increased concordance from 0.592 to 0.622. D1 therefore overlaps substantially with these classifications
without being reducible to them." Object `d1_variance_partition_stats` (code/14).
Test (R20a): `prognosis$d1_adj_hi` < 1 (CI excludes 1, protective).
Test (R20b): `lr_p` < 0.05 AND `c_with_d1` > `c_base`.
Test (R20c "overlaps substantially"): 0.20 <= `variance$explained_total` <= 0.80. [OWNER]
Test (R20d, SI L1753-1762, interaction-model sensitivity): `interaction$lr_p` < 0.05.

**R21.** MAIN, L569; SI L1769; Discussion L33. "D1 correlated only weakly with purity (Spearman rho = 0.06 and
-0.05)" and "adjustment for purity left the DeSurv linear predictor's survival association essentially unchanged
(PACA-AU array 1.44 versus 1.45; RNA-seq 1.66 versus 1.65)". Object `desurv_clinical_adjustment` (code/18).
Test (R21a): |rho(D1, purity)| < 0.30 in PACA_AU_array and PACA_AU_seq. [OWNER]
Test (R21b): |HR_purity / HR_unadj - 1| < 0.10 in PACA_AU_array and PACA_AU_seq (SI caption L1790 extends to
TCGA 3.50 vs 3.60: include TCGA_PAAD). [OWNER]

**R22.** MAIN, L569. "In Puleo ... the DeSurv linear predictor also remained prognostic after adjustment for sex
and resection-margin status (HR 1.33, 95% CI 1.14--1.55)". `desurv_clinical_adjustment$Puleo_array$
clinical_adjusted`. Test: lo > 1.

**R23.** MAIN, L575; Fig. 4 caption L578; SI L2403. O'Kane/COMPASS: "D1 remained associated with overall
survival in a Cox model stratified by treatment arm (HR 0.75, 95% CI 0.62--0.90)"; caption and SI "D1 had
protective point estimates in all three arms, with no evidence of an arm-by-score interaction (P = 0.40)".
Object `treated_cohort_stats$okane` (code/11). Current per-arm 0.80 / 0.57 / 0.83.
Test (R23a): stratified `hi` < 1. Test (R23b): per-arm HR < 1 in all 3 arms. Test (R23c): `interaction_p`
>= 0.05. Discussion L30 "In treated cohorts the frozen scores remained informative" = R23a AND R24.

**R24.** MAIN, L575; caption L578; SI L2409. Linehan: "the projected proCAF-like D2 score increased from pre- to
post-treatment biopsies (mean change +0.34 SD, paired Wilcoxon P = 0.013)". `treated_cohort_stats$paired`.
Test: `mean_delta` > 0 AND `wilcox_p` < 0.05 AND L-D2 holds.

**R25.** MAIN, L601; SI L2088. "Reference-program coverage ... increased from 16 for DeSurv at k = 3 to 32 for
the k = 7 control, while multiplicity ... remained low (1.00 versus 1.125)". Object
`reference_program_multiplicity$anchors`, anchor `desurv_k3` (code/24). Test: DeSurv k3 coverage < 32 (fixed
k = 7 value) AND DeSurv k3 multiplicity <= 1.25. [OWNER]

**R26.** MAIN, L601 closing sentence. "Survival supervision was not required to obtain a prognostic signal, but
it altered which programs were identified and how that signal was decomposed." Test: R19c AND R9a.

**R27.** MAIN, L612; SI L1494-1506; Fig. S suprec caption L1510. "All three produced patient scores strongly
correlated with D1 (Pearson |r| = 0.81, 0.83 and 0.85), whereas the leading unsupervised principal component was
only weakly correlated with D1 (|r| = 0.10) and instead tracked the higher-variance D2/D3 structure".
Object `supervised_recovery_stats$table` (code/15; comparator fits fixed, correlations recomputed against the
rerun's D-scores). Current PC1: |r| D2 0.98, D3 0.93.
Test (R27a): |score_cor_D1| >= 0.70 for Supervised PCA, Cox-PLS and Sparse Cox. [OWNER]
Test (R27b): |r(PC1, D1)| < 0.30 AND |r(PC1, D2)| and |r(PC1, D3)| both > |r(PC1, D1)|. [OWNER]
Test (R27c, SI L1498 "each aligned with D1 rather than with D2 or D3"): for each supervised method |r_D1| >
max(|r_D2|, |r_D3|). (Current max |r| with D2/D3 = 0.19.)

**R28.** MAIN, L615; SI L1540-1546; Fig. 3b caption L472. "Independently tuned supervised PCA and lasso Cox had
pooled external C-indices of 0.60 and 0.60, versus 0.61 for DeSurv"; "Their risk scores correlated only weakly
with the ... D2 program (|r| = 0.06; 0.14)"; "correlating with an external activated-stroma and proCAF
signature at |r| = 0.20 and 0.23, against 0.55 for D2"; "the methods differed less in prognostic discrimination
than in how that information was represented". Object `desurv_vs_supervised_tuned` (code/19).
Test (R28a "differed less in discrimination"): |C_DeSurv full LP - C_method| < 0.03 for both methods
(`prediction[, "POOLED"]`). [OWNER]
Test (R28b): |`axis_decomposition`[method, "D2"]| < 0.30 for both. [OWNER]
Test (R28c): `stromal_carriage`["DeSurv D2 (proCAF factor)"] > both supervised values.
Test (R28d, SI L1540 "both tuned methods aligned with the D1 direction"): |r_D1| > max(|r_D2|, |r_D3|) in
`axis_decomposition` for both. Current Sparse Cox 0.71 vs 0.06/0.10; SPCA 0.46 vs 0.14/0.31.
(The basal-like enrichment of the comparators' selected genes, 12/12 seeds, does not depend on the DeSurv fit:
excluded.)

**R29.** MAIN, L440; SI L2217-2230. "Both methods recovered activated-stroma and proCAF-like reference programs
but differed at the tumor-restraining end: among the displayed top-ranked genes, the restCAF and iCAF reference
programs corresponded to a DeSurv factor but to no standard NMF factor". Statistic: Fig. 2a Spearman matrix
(0.4) and `stroma_resolution_stats$top_n` (code/23).
Test (R29a): DeSurv max rho > 0.2 for SCISSORS_restCAF AND for Elyada_iCAF (current 0.24, 0.24).
Test (R29b): DeSurv max rho > 0.2 for DECODER_ActivatedStroma AND SCISSORS_proCAF (current 0.37, 0.23).
Test (R29c, SI L2224-2227): restCAF genes among DeSurv defining genes >= 1 (current 2 of 25) AND iCAF >= 1
(current 2 of 35); NMF side fixed at 0 and 0.

---

## 3. Discussion (paper/05_discussion_REVISED.Rmd)

**DC1.** MAIN, L9. "the architecture comprised a principal survival-aligned program, D1, characterized by
coordinated cross-patient variation in classical tumor and restCAF-like stromal features, alongside separately
resolved proCAF-like stromal (D2) and basal-like tumor (D3) programs. D1 remained prognostic after adjustment
for ... PurIST and DeCAF ... The frozen programs transferred without refitting to four independent external
cohorts." Test: A1 AND R3 AND R20a AND R10.

**DC2.** MAIN, L12 (plan-required). "An unsupervised control matched to DeSurv on rank and tuning budget did not
recover D1 ... A rank-optimized unsupervised control recovered comparable prognostic information ...
Conventional supervised predictors achieved similar discrimination but combined contributions from underlying
expression programs within single scores." Test: R9a AND R19c AND R28a.

**DC3.** MAIN, L18. "the fitted risk score combines D1 and D3, but D2 remains available as proCAF-like stromal
context." Test: R15 AND L-D2.

**DC4.** MAIN, L21. "Individual runs from different random starts recover appreciably different gene sets.
Split-half consensus fits agreed closely with the reported basis". Test: R1 (Jaccard part) AND R2.

**DC5.** MAIN, L27. PH attenuation (= R17) and "Within-cohort analyses argued against the risk score simply
encoding training-cohort membership" (= S3).

**DC6.** MAIN, L30. "In treated cohorts the frozen scores remained informative in exploratory analyses". = R23a
AND R24.

**DC7.** MAIN, L33. "Purity adjustment in the PACA-AU datasets left the risk score's association essentially
unchanged ... GATA6 RNA in situ hybridization independently supports D1's classical tumor component". = R21b AND
R4 AND L-D1 (classical row).

**DC8.** MAIN, L39. "incorporating survival during factorization recovered a three-program tumor--stroma
architecture that transferred across independent cohorts". = A1 AND R10.

---

## 4. SI Results-type sections (paper/si_appendix.Rmd)

**S1.** SI, L1052-1053. D2 coefficient exactly zero; LP a weighted sum of D1 and D3. = R15.

**S2.** SI, L1092-1100 / L1646. Combined event-weighted C-index and CI. = R18c.

**S3.** SI, L1145 (Cohort handling). "the fitted linear predictor remained strongly associated with survival
within TCGA-PAAD and CPTAC separately (in-sample HR per SD 3.50 and 3.72) and in a model stratified by cohort
(3.35, against 3.12 unstratified), while the between-cohort shift in the score was opposite to the between-cohort
hazard difference: the predictor averaged 0.73 SD higher in TCGA-PAAD, yet CPTAC had the higher hazard given the
predictor (HR 1.90)." Accessors `.ct_*` SI L1128-1140 (note: SI scores here use per-program masks).
Test (S3a): within-cohort unadjusted HR 95% CI lower bound > 1 in both TCGA_PAAD and CPTAC
(`desurv_clinical_adjustment`). Test (S3b "opposite"): sign(`.ct_shift`) != sign(log `.ct_cohhr`), where `.ct_shift` = mean LP(TCGA_PAAD) -
mean LP(CPTAC) and `.ct_cohhr` = HR for TCGA_PAAD vs CPTAC (reference level) given the LP. Current
`.ct_shift` = -0.73, `.ct_cohhr` = 1.90: PASS.
DEFECT FOUND (pre-existing, prose only): verified in R, the TCGA_PAAD mean LP is 0.73 SD LOWER than CPTAC
(means -0.344 vs 0.384) and TCGA_PAAD, not CPTAC, has the higher hazard given the predictor (coef
cohTCGA_PAAD = 0.64, HR 1.90). SI L1145 names both cohorts the wrong way round ("0.73 SD higher in TCGA-PAAD,
yet CPTAC had the higher hazard"); it prints `abs(.ct_shift)`, which hides the sign. The "opposite" conclusion
is unaffected. If the owner wants the named direction scored, test S3c: `.ct_shift` > 0 AND `.ct_cohhr` < 1
(currently FAILS).

**S4.** SI, L1246. "Individual shares may be signed in principle but were all positive in the reported PDAC
models." Test: all three DeSurv Shapley shares > 0.

**S5.** SI, L1289-1294. Delta-loglik concentrates in D1, NMF factors small. = R6a AND R6c.

**S6.** SI, L1435-1442; Table projection-r2 caption L1430. "The biologically dominant tumor (N1) and stromal (N3)
programs ... are largely reconstructable from the DeSurv basis (R^2 = 0.90 and 0.88) ... reconstructing the
exocrine factor (N2) less well (R^2 = 0.38)"; "the DeSurv basis resolves tumor signal into basal- and
classical-aligned programs". Accessor `.r2_N` SI L1421. Test (S6a): R^2(N1 | DeSurv) >= 0.80 AND
R^2(N3 | DeSurv) >= 0.80. [OWNER] Test (S6b): R^2(N2 | DeSurv) < min(R^2 N1, R^2 N3) AND < 0.50. Test (S6c):
L-D1 classical row AND L-D3.

**S7.** SI, L1457. Reconstruction gap "real but bounded". = R8b.

**S8.** SI, L1486-1507. D1 axis shared across supervised methods. = R27a-c. Also "D1 is unreconstructable from
the standard NMF basis" = R8; "standard NMF at matched rank does not recover the restCAF and iCAF marker genes
that appear among DeSurv's defining genes" = R29c.

**S9.** SI, L1540-1546. Tuned comparators. = R28a-d; "DeSurv preserves stromal context as an explicit program
even when that program does not contribute independently to the fitted linear predictor" = R15 AND R28c.

**S10.** SI, L1576-1599. Brier. = R18d. "consistent with the moderate discrimination (C = 0.61)" = R18c.

**S11.** SI, L1621. PH diagnostic. = R17.

**S12.** SI, L1640-1646. Heterogeneity LRT, Q/I^2, C-index CI. = R18a-c. Also within-platform standardization
sensitivity LRT (`sensitivity$lrt_within_platform_standardization$p`) >= 0.05 (S12a).

**S13.** SI, L1649-1715 (LODO table tab:loco). = R16.

**S14.** SI, L1720-1762 (Variance partition). = R20a-d. "D1's prognostic information is therefore not captured
by the saturated four-level joint classification" = R20d.

**S15.** SI, L1769 and Table clinical-adjustment caption L1790. "After clinical adjustment the DeSurv score
remains independently prognostic in Puleo, PACA-AU array and both training cohorts, and is attenuated to
non-significance in PACA-AU RNA-seq, Dijk and Moffitt"; "Moffitt is non-significant both unadjusted and after
adjustment". Object `desurv_clinical_adjustment$<cohort>$clinical_adjusted$p` and `$unadjusted$p`.
Current clinical-adjusted P: Puleo 2.1e-4, PACA-AU array 0.025, TCGA 4.6e-13, CPTAC 1.6e-14; PACA-AU seq 0.077,
Dijk 0.19, Moffitt 0.15; Moffitt unadjusted 0.076.
Test (S15a): P < 0.05 in Puleo_array, PACA_AU_array, TCGA_PAAD, CPTAC. Test (S15b): P >= 0.05 in PACA_AU_seq,
Dijk, Moffitt_GEO_array. Test (S15c): Moffitt unadjusted P >= 0.05. Purity part = R21b.

**S16.** SI, L1816-1819 (Empirical optimization stability). "100 of the 100 runs terminated on the stopping
criterion ..., none failed numerically, and in every run the recorded outer-loop objective decreased
monotonically". Guard in the Rmd: `stopifnot(.n_conv == .R, .n_maxit == 0, !any(.nan), all(.mono))`.
Test (per rerun, on its own seed fits): same four conditions.

**S17.** SI, L1851-1862. Restart dispersion and near-constant training C. = R1 (SI adds "Matched gene-loading
correlations were moderate"; test: median matched loading r per factor in [0.3, 0.8]. [OWNER])

**S18.** SI, L1865-1894. Split-half consensus agreement; "The consensus step is what makes the reported
decomposition stable to restart-pool composition". = R2.

**S19.** SI, L1881-1888. "The fitted basis also moved materially from the consensus initialization W0 ...
gene-loading correlations between W0 and the final W were 0.57--0.81 and top-270 Jaccard 0.44--0.82." Object
`zeta_trace$w0` (code/26). Test (per rerun, from its own W0 and final W): min per-factor r(W0, W) < 0.90.
[OWNER: "materially".]

**S20.** SI, L2088. DeSurv anchor coverage/edges/multiplicity. = R25.

**S21.** SI, L2217-2245 (Stromal resolution). Top-gene view = R29a-c. Full-gene view: "no method separates
restCAF from proCAF: every correspondence lies between 0.04 and 0.10 ... whereas the activated-stroma axis is
recovered by all methods at 0.19". Object `stroma_resolution_stats$full_rho` (DeSurv column recomputed per
rerun). Test (S21a): DeSurv |full_rho| for SCISSORS_restCAF, Elyada_iCAF, SCISSORS_proCAF, Elyada_myCAF all
< DeSurv full_rho for DECODER_ActivatedStroma.

**S22.** SI, L2284-2293 (Matched-budget control). = R9a-c. "Matching the tuning budget therefore does not close
the representational gap" = R9a.

**S23.** SI, L2296. k = 7 exceeds DeSurv C-index in all five datasets = R19d; both retain significance after
adjustment = R19b.

**S24.** SI, L2395 and L2403-2409 (treated cohorts). = R15, R23a-c, R24. SI L2409 "indicating a shift toward
proCAF-associated bulk expression" additionally requires L-D2.

---

## 5. Judged fit-independent but borderline (one line each)

- paper.Rmd L5 / Results L423: k = 3 selected by the one-SE rule, alpha, and the BO surface. Depends on the BO,
  not the fitted solution; the reruns hold hyperparameters fixed, so these are not scored.
- Results L423 / SI L2014-2031: cross-validated C-index across k and the cv-grid points are single grid fits
  independent of the primary fit.
- Results L460: "the union of their top-ranked genes (N genes)" is a count of the projection support; it
  changes with W but is a descriptive number, not a claim.
- Results L463: "None of the three standard NMF factors at the same rank showed a significant combined
  association" involves only the fixed NMF fit.
- Results L598 and Table hr-adjusted: the matched-budget control's own non-significance (unadj P = 0.43) involves
  only the fixed control.
- Results L601 / SI L2102: k = 7 coverage (32) and the supervised-curve collapse at k = 6, 8 are fixed control
  or grid fits.
- Results L615 / SI L1543: the comparators' selected genes enriched for basal-like tumor in 12/12 seeds; the
  comparator gene counts. Comparator-only.
- Results L566, SI L1720-1734: "DeCAF alone explained 19%" and PurIST/DeCAF main-effect numbers are reported
  values under R20c, not separate claims.
- SI L1904-1924: 80% subsample single-run stability (DeSurv D1 Jaccard 0.34, "mid-range"). Computed from
  separate single-run refits in code/16 at fixed hyperparameters, not from the primary fit or a rerun; not
  rescored per rerun.
- SI L1103-1117 (independence of discovery and validation) and L1028-1089 (scaling conventions): method
  definitions.
- SI L2253-2281 (matched-budget control hyperparameters, own-ntop sensitivity): control-only.
- SI L1178-1200 (reference universe) and L2248-2250 (label nomenclature): data definitions.
- All simulation claims (Results L637-646, SI L1927-2009): fixed truth.

---

## 6. Findings and conflicts with the brief

1. **The primary fit does not exist yet.** The inventory is frozen at Step 2, before BO (Step 4) and the
   primary fit (Step 5). All current values are from the fb731e7 fit. Tests are written relative to "primary
   D1/D2/D3" (alignment A), which presumes the new primary fit has k = 3 and that its factors can be labeled.
   The label test L is anchored to reference programs so it can first be applied to the new primary fit. If
   Step 4 selects k != 3, or the new primary fit fails L, most of this inventory is undefined and needs owner
   re-freezing before Step 6.
2. **Factors are identified by column position everywhere** (see 0.1). The scoring harness must permute each
   rerun's W, H and beta into primary order before running code/05, 09a, 11, 14, 15, 18, 19, 23, 24, 25. The
   Rmd guard `beta[2] == 0` (Results L466) and the hard-coded heatmap labels (`code/09a_figures.R:161-163`)
   will silently misreport an unpermuted rerun.
3. **Reference gene lists live outside the worktree** (`../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData`,
   Amber's repo). The label test and R29 depend on them. Pin before scoring.
4. **Claims already marginal at the current fit** (likely to land in the 5-8/10 band): R14b "D2 varied more
   across datasets" (SD of per-dataset log HR 0.128 vs D3 0.127, D1 0.122); L-D1 restCAF and R29a (rho 0.24
   vs threshold 0.2); L-D2 proCAF (0.23); R9b control recovers D3 (0.749 vs 0.70); R27c supervised scores' max
   |r| with D2/D3 is 0.19.
5. Several claims use qualitative words ("closely", "weakly", "modest", "comparable", "little", "materially",
   "essentially unchanged"). Each has a proposed numeric rule marked [OWNER]. These thresholds are not in the
   plan and need owner sign-off before they are frozen.
6. The brief says 88 tracked `.rds`; `git ls-files results | grep -c '\.rds$'` gives 83 at 88cf0f7.
7. `program_stability_stats$pairwise` has 435 rows per factor (30 choose 2) while `meta$n_restarts` = 100 and
   the text says R = 100 (R1, S17). Confirm the restart set before reusing code/20 per rerun.
8. **SI L1145 cohort direction is inverted in the prose** (see S3): the score is higher in CPTAC, and TCGA-PAAD
   has the higher conditional hazard. Text fix belongs to Step 7; no edit made here.
9. The worktree had uncommitted Agent A edits during this pass (paper Rmds, code/*.R, R/paths.R). They are
   line-preserving in the manuscript files, so the line numbers above hold at 88cf0f7.
