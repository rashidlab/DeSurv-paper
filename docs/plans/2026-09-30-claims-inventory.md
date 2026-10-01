# Claims inventory, rev 3 Step 2 (2026-09-30, revised on the owner ruling of 2026-09-30)

Frozen before any rerun output exists. Plan context: `2026-09-30-rerun-plan-rev3.md` sections 4 and 8. Line
numbers are at HEAD 88cf0f7 (concurrent Agent A edits were line-preserving; checked with `git diff -U0`).
"Current value" is from the tracked `results/*.rds` (fb731e7 PDFs), read in a throwaway R session that
evaluated the purled Rmd chunks; nothing in `results/` was written. Values carried over from the first draft
were not re-derived.

Scope: every substantive claim in the abstract (`paper/paper.Rmd` L5), `paper/04_results_REVISED.Rmd`,
`paper/05_discussion_REVISED.Rmd` and the Results-type SI sections of `paper/si_appendix.Rmd` whose truth
depends on the fitted DeSurv solution (W, H, beta, factor labels, top genes, validation scores). Location codes:
ABS = abstract, MAIN = main text, SI = Supplementary Information.

## Rules (owner ruling, 2026-09-30)

1. Ordinary descriptive words ("weakly", "bounded", "overlaps substantially", "closely", "comparable",
   "little", "modest", "materially", "essentially unchanged") are not numerical gates. In the final manuscript
   such wording is replaced by the actual quantity. Here they are DESCRIPTIVE QUANTITIES (Section B).
2. Every entry is one of two classes.
   - BINARY CLAIMS (Section A): categorical scientific claims (direction or sign, ordering, significance as the
     manuscript states it, which factor is largest or smallest, presence or absence of an association, exact
     counts). Each has one mechanical, exact test. Scoring rule (plan section 8): holds in >= 9/10 reruns ->
     ordinary finding (frequency reported in the sensitivity SI); 5-8/10 -> described as solution-sensitive
     with its frequency; < 5/10 -> not asserted.
   - DESCRIPTIVE QUANTITIES (Section B): HR magnitudes, C-indices, correlations, overlap values,
     reconstruction shares, R^2, SDs and similar. Reported as primary-fit value plus range and distribution
     over the 10 reruns. No pass/fail.
   A sentence carrying both is split: the categorical part goes to A, the number to B. Where a word was
   categorical only through a threshold the first draft invented, it is now descriptive only ("manuscript to
   state the quantity").
3. D1/D2/D3 identity in a rerun comes from alignment to the frozen primary fit (Hungarian on W correlation),
   with no threshold. A biological label holds when the aligned factor's greatest top-gene correspondence is
   with the same reference program as in the primary fit. The correspondence values themselves are
   descriptive. The fixed 0.20 cutoff and all margins are dropped.
4. A number word in the text is tested literally ("four" means exactly four).
5. [OWNER] tags remain only on categorical claims these rules do not settle (Section C).

---

## 0. Shared definitions

**0.1 Alignment (A).** For rerun j, compute the k x k matrix of Pearson correlations between the columns of
W_j and W_primary over all analyzed genes and solve the assignment problem (Hungarian, `clue::solve_LSAP`,
maximize). Permute W_j columns, H_j rows and beta_j entries into primary order. "D1/D2/D3 in rerun j" means
the factor matched to primary D1/D2/D3. No threshold on the matched r. The matched Pearson r, matched
Spearman r and top-ntop Jaccard per factor are reported as descriptive quantities (Q-ALN).

Why the permutation is required: every downstream script identifies programs by column position, not by
label. Examples: `04_results_REVISED.Rmd` L237 (`D1 = Z[, 1]`), L147-148 (`.corrC[1, ]`), L466 guard
`beta[2] == 0`; `code/14_d1_variance_partition.R:40`; `code/11_treated_cohort_analysis.R:57`;
`code/15_supervised_recovery.R:42`; `code/18_clinical_adjustment.R:34,44`;
`code/19_desurv_vs_supervised_tuned.R:82-128`; heatmap labels hard-coded by position at
`code/09a_figures.R:161-163`. Unpermuted rerun outputs would score the wrong factor.

**0.2 Per-rerun recomputation.** A claim is scored from rerun j's own outputs, produced by the same code that
produced the current value, applied to the permuted fit, into a rerun-specific results directory: fit (code/04
recipe) -> `val_latent_desurv_tcgacptac` (code/05) -> the named downstream script. Fixed comparators
(`fit_std_desurvk`, `tar_fit_desurv_a0k3`, `tar_fit_desurv_alpha0`, their `val_latent_*`, the supervised
comparators' own fits) are NOT rerun; they are read from the primary build.

**0.3 Significance.** "Significant" / "associated" = two-sided Wald P < 0.05 from the named Cox model,
equivalently the 95% CI excludes 1. "Not significant" = P >= 0.05. This is the manuscript's own convention,
not a threshold added here.

**0.4 Top-gene correspondence matrix (M).** The factor x reference-program Spearman matrix
`R/figure_plot_helpers.R::make_gene_overlap_heatmap(fit, get_top_genes(W, ntop)$top_genes, top_genes,
return_matrix = TRUE)$cor_mat` (W loadings over the union of each factor's top 50 genes, correlated with
reference-program membership; the Fig. 2a statistic). Reference lists from `cmbSubtypes_formatted.RData`,
which lives in `../DeSurv-paper/data/derv/` (see D2). "Row argmax" = the reference program with the largest
entry in a factor's row, over all reference programs in M. "Column argmax" = the factor with the largest entry
in a program's column, over the k factors.

---

## Section A. Binary claims

Format: ID (old ID) | location | file:line | claim | test | primary current value. All tests run on the
alignment-A-permuted rerun.

### A.1 Label tests

Owner rule 3 test (row argmax):

- **L-D1** (old 0.4 L-D1) | ABS, MAIN, SI | Results L437, paper.Rmd L5 | D1 is the classical tumor /
  restCAF-like program. Test: row argmax of M for aligned D1 = row argmax of M for primary D1. Current: the
  primary program is recorded at freeze (the first draft recorded column values only; see Q-LBL).
- **L-D2** (old L-D2) | MAIN, SI | Results L437 | D2 is the proCAF-like stromal program. Test: row argmax of M
  for aligned D2 = that of primary D2. Current: recorded at freeze.
- **L-D3** (old L-D3) | MAIN, SI | Results L437 | D3 is the basal-like tumor program. Test: row argmax of M
  for aligned D3 = that of primary D3. Current: recorded at freeze.

Named-component tests (PROPOSED pending owner question C1; each label names two reference programs and a
single row argmax can carry only one of them). Test for each: the column argmax of M for that program is the
same aligned factor as in the primary fit. No threshold, no margin.

- **L-D1c** DECODER_ClassicalTumor -> D1. Current D1 0.40 / D2 -0.40 / D3 -0.37: D1.
- **L-D1r** SCISSORS_restCAF -> D1. Current 0.24 / -0.13 / -0.09: D1. (The first draft's Elyada_iCAF fallback
  is dropped: an either-or weakens the claim; iCAF values stay in Q-LBL.)
- **L-D2a** DECODER_ActivatedStroma -> D2. Current -0.06 / 0.37 / 0.06: D2.
- **L-D2p** SCISSORS_proCAF -> D2. Current -0.19 / 0.23 / 0.04: D2.
- **L-D3b** DECODER_BasalTumor -> D3. Current D3 0.34: D3.
- **L-D3p** Puleo_Pure Basal-like -> D3. Current D3 0.39: D3.

A label "holds" in rerun j when L-Dk holds and, if C1 is approved, each of its named-component tests holds.
Any claim that names a label holds only if the label holds and the claim's own test holds. Claims about D1/D2/D3
by number only use alignment A alone.

### A.2 Elementary binary claims

- **R4** (R4) | MAIN | Results L437; Discussion L33 | D1 scores correlated with GATA6 RNA ISH in COMPASS.
  Test: `treated_cohort_stats$gata6` rho > 0 AND P < 0.05. Current: rho 0.679, P 1.4e-5: PASS.
- **R5** (R5) | MAIN | Results L443 | D1 had the smallest reconstruction share of the three factors. Test: D1
  normalized Shapley share < D2's AND < D3's (`fig_variation_explained_tcgacptac$data`, `.vp` code/09a L114).
  Current: 18.8 < 25.2, 56.0: PASS.
- **R6a** (R6a) | ABS (via A1), MAIN, SI | Results L443; SI L1289-1294; Fig. 2c caption L451 | Omitting D1
  cost the most partial log-likelihood (the "survival-aligned" program). Test: Delta-loglik(D1) > Delta-loglik(D2)
  AND > Delta-loglik(D3) (`.dl`, code/09a L113). Current: 66.37 vs 0.02, 1.00: PASS.
- **R6c** (R6c) | SI | SI L1293 | Survival contribution concentrates in DeSurv D1 rather than any NMF factor.
  Test: Delta-loglik(D1) > max NMF Delta-loglik (fixed, 1.16). Current 66.37 > 1.16: PASS.
- **R6d** (R6d, R6b) | MAIN | Results L443 | "D3 showed the opposite pattern": largest reconstruction share but
  not the survival-aligned program. Test: D3 Shapley share > D1's AND > D2's, AND Delta-loglik(D3) <
  Delta-loglik(D1). Current: 56.0 largest; 1.00 < 66.37: PASS. ("little conditional survival contribution" is
  descriptive: Q-R6.)
- **R7a** (R7a, split) | MAIN | Results L446 | The principal NMF tumor factor corresponds to D3 and the
  microenvironmental factor to D2. Test: which.max(|`.corrC`[3, ]|) = N1 (tumor, column 1) AND
  which.max(|`.corrC`[2, ]|) = N3 (microenvironment, column 3); `.corrC` = `cor(W_desurv, W_nmf, "spearman")`,
  Results L146; NMF column labels from SI L1424. Current: rho 0.845 and 0.873; argmax indices recorded at
  freeze (not captured in the first draft). "closely" is descriptive (Q-R7).
- **R7b** (R7b, split) | MAIN | Results L446 | D1's correspondence to NMF is weaker than that of D2 and D3.
  Test: max|`.corrC`[1, ]| < max|`.corrC`[2, ]| AND < max|`.corrC`[3, ]|. Current 0.357 < 0.873, 0.845: PASS.
  "only weakly (rho <= 0.36)" is descriptive; manuscript to state the quantity (Q-R7).
- **R8b** (R8b, split) | MAIN, SI | Results L446; SI L1457 | DeSurv reconstructs the analysis matrix less
  completely than matched-rank NMF. Test: `overall_recon_r2_centered` NMF > DeSurv (`.r2c_*`, Results
  L132-134). Current 0.812 > 0.752: PASS. Gap and "bounded" are descriptive (Q-R8b).
- **R9a** (R9a, split) | ABS (A4), MAIN, SI | Results L449; SI L2284-2293; Discussion L12 | The matched-budget
  control did not recover D1 while recovering D2 and D3. Test: max|`.corrF`[1, ]| < max|`.corrF`[2, ]| AND
  < max|`.corrF`[3, ]| (`.corrF` = `cor(W_desurv, W_a0k3, "spearman")`, Results L162-166). Current 0.245 <
  0.861, 0.749: PASS. The rho values and "failed to recover (rho <= 0.25)" are descriptive (Q-R9).
- **R9c** (R9c) | SI | SI L2285 | D1's best match in the matched-budget control is lower than in NMF::nmf.
  Test: max|`.corrF`[1, ]| < max|`.corrC`[1, ]|. Current 0.245 < 0.357: PASS. "slightly" is descriptive.
- **R10** (R10, A2) | ABS, MAIN | Results L466; Discussion L9 | The linear predictor was associated with
  overall survival (dataset-stratified). Test: 95% CI lower bound > 1 (`pooled_cox_stats(pool_df)`, Results
  L207-219, L301). Current 1.42 (1.26-1.59): PASS.
- **R11** (R11) | MAIN | Results L463 | D1 had protective point estimates in all five validation datasets.
  Test: per-dataset HR per within-dataset SD of D1 < 1 in 5 of 5 (Fig. 3a statistic). Current 0.77, 0.83,
  0.69, 0.60, 0.73: PASS.
- **R12** (R12) | ABS (via A1), MAIN | Results L463 | D1 was the only individual program with a significant
  combined association, protective. Test: P(D1) < 0.05 AND HR(D1) < 1 AND P(D2) >= 0.05 AND P(D3) >= 0.05
  (`.pf`, Results L233-249). Current D1 P 1.6e-8 HR 0.75; D2 P 0.76; D3 P 0.10: PASS.
- **R13** (R13) | MAIN | Results L463 | D3 had HRs above 1 in four of the five datasets. Test (literal, rule 4):
  number of datasets with per-dataset D3 HR > 1 equals exactly 4. Current 1.27, 1.07, 0.95, 1.26, 1.04 -> 4:
  PASS. Note: a rerun with 5 of 5 fails this test; if the reruns vary the count, the sentence should state the
  count as computed rather than be loosened.
- **R14a** (R14a) | MAIN | Results L463 | D2 was not significant overall. Test: pooled D2 P >= 0.05 (`.pf`).
  Current P 0.76: PASS. "near null" is descriptive (Q-R12).
- **R14b** (R14b) [RULED C2] | MAIN | Results L463 | D2 "varied more across datasets". Proposed test: SD across
  the five datasets of per-dataset log HR (R11 statistic) for D2 > the same SD for D1 AND for D3. Current 0.128
  vs 0.122, 0.127: PASS.
- **R15** (R15, S1) | MAIN, SI | Results L466; SI L1052-1053, L2395; Discussion L18 | The LP combines D1 and
  D3; D2's coefficient is exactly zero. Test: beta_D2 == 0 AND beta_D1 != 0 AND beta_D3 != 0 (after
  permutation A). Current (-0.00100, 0, 9.36e-5): PASS. The Rmd guard at L466 breaks on an unpermuted rerun.
- **R17** (R17, S11) | MAIN, SI | Results L466; SI L1621; Discussion L27 | Scaled Schoenfeld residuals
  indicated attenuation over follow-up. Test: `validation_heterogeneity_ci$ph` P < 0.05 AND `slope_time` < 0.
  Current P 1.3e-4, slope -0.92: PASS.
- **R18a** (R18a) | MAIN, SI | Results L469; SI L1640-1646 | Cohort-specific slopes did not significantly
  improve fit. Test: LRT P >= 0.05 (`validation_heterogeneity_ci`, code/25). Current P 0.69: PASS.
- **R18d** (R18d, split) | MAIN, SI | Results L469; SI L1576-1599 | The recalibrated score improves on the
  KM reference (proper-score improvement over the null). Test: pooled IPA > 0 (`desurv_brier_stats`,
  code/17). Current IPA 5.4% (IBS 0.159 vs 0.168): PASS. "modest" is descriptive.
- **R19a** (R19a, A3, R19b) | ABS, MAIN, SI | Table hr-adjusted Results L478-547, L598; SI L2296 | DeSurv LP
  remains prognostic after PurIST + DeCAF adjustment. Test: adjusted HR 95% CI lower bound > 1
  (`pooled_cox_unadj_adj`, Results L484-502). Current 1.28 (1.11-1.47): PASS.
- **R19d** (R19d) | SI | SI L2296 | The k = 7 control exceeds DeSurv's C-index in all five validation datasets.
  Test: control C > DeSurv C in 5 of 5 (`val_cindex_desurv_tcgacptac` vs fixed
  `val_cindex_desurv_alpha0_tcgacptac`). Current DeSurv 0.607/0.562/0.649/0.634/0.636 vs
  0.631/0.578/0.669/0.702/0.651: PASS.
- **R20a** (R20a) | MAIN, SI | Results L566; Discussion L9; SI L1720-1762 | D1 remained associated with
  survival, protective, after adjustment for both classifiers. Test: `prognosis$d1_adj_hi` < 1
  (`d1_variance_partition_stats`, code/14). Current HR 0.82 (0.73-0.93): PASS.
- **R20b** (R20b) | MAIN, SI | Results L566 | Adding D1 improved fit and increased concordance. Test: `lr_p` <
  0.05 AND `c_with_d1` > `c_base`. Current P 0.0015; 0.622 > 0.592: PASS.
- **R20d** (R20d) | SI | SI L1753-1762 | D1's information is not captured by the saturated four-level joint
  classification. Test: `interaction$lr_p` < 0.05. Current: value in the object (not transcribed in the first
  draft); recorded at freeze.
- **R22** (R22) | MAIN | Results L569 | In Puleo the LP remained prognostic after sex and margin adjustment.
  Test: `desurv_clinical_adjustment$Puleo_array$clinical_adjusted` lower bound > 1. Current 1.33
  (1.14-1.55): PASS.
- **R23a** (R23a) | MAIN, SI | Results L575; SI L2403 | D1 associated with OS in COMPASS stratified by arm,
  protective. Test: stratified `hi` < 1 (`treated_cohort_stats$okane`). Current 0.75 (0.62-0.90): PASS.
- **R23b** (R23b) | MAIN, SI | Fig. 4 caption L578; SI L2403 | D1 protective point estimates in all three arms.
  Test: per-arm HR < 1 in 3 of 3. Current 0.80 / 0.57 / 0.83: PASS.
- **R23c** (R23c) | MAIN, SI | Fig. 4 caption L578; SI L2403 | No evidence of an arm-by-score interaction.
  Test: `interaction_p` >= 0.05. Current 0.40: PASS.
- **R24** (R24) | MAIN, SI | Results L575; caption L578; SI L2409 | The proCAF-like D2 score increased from pre-
  to post-treatment (Linehan). Test: `treated_cohort_stats$paired` `mean_delta` > 0 AND `wilcox_p` < 0.05 AND
  label D2 holds. Current +0.34 SD, P 0.013: PASS.
- **R25** (R25, split) | MAIN, SI | Results L601; SI L2088 | Reference-program coverage increased from DeSurv
  k = 3 to the k = 7 control. Test: DeSurv k3 coverage (`reference_program_multiplicity$anchors`, anchor
  `desurv_k3`, code/24) < 32 (fixed k = 7 value). Current 16: PASS. Multiplicity "remained low" is
  descriptive (Q-R25).
- **R27b** (R27b, split) | MAIN, SI | Results L612; SI L1494-1506 | PC1 tracked D2/D3 rather than D1. Test:
  |r(PC1, D2)| > |r(PC1, D1)| AND |r(PC1, D3)| > |r(PC1, D1)| (`supervised_recovery_stats$table`, code/15).
  Current 0.98, 0.93 > 0.10: PASS. "only weakly" is descriptive.
- **R27c** (R27c) | SI | SI L1498 | Each supervised method aligned with D1 rather than D2 or D3. Test: for SPCA,
  Cox-PLS and Sparse Cox, |r_D1| > max(|r_D2|, |r_D3|). Current |r_D1| 0.81-0.85 vs max 0.19: PASS.
- **R28c** (R28c) | MAIN, SI | Results L615; SI L1540-1546 | D2 carries the external activated-stroma /
  proCAF signature more than either tuned supervised score. Test: `stromal_carriage` DeSurv D2 > each
  supervised value (`desurv_vs_supervised_tuned`, code/19). Current 0.55 > 0.20, 0.23: PASS.
- **R28d** (R28d) | SI | SI L1540 | Both tuned methods aligned with the D1 direction. Test: in
  `axis_decomposition`, |r_D1| > max(|r_D2|, |r_D3|) for both. Current Sparse Cox 0.71 vs 0.06/0.10; SPCA 0.46
  vs 0.14/0.31: PASS.
- **R29c** (R29a, R29b, R29c merged) [RULED C3] | MAIN, SI | Results L440; SI L2217-2230 | The restCAF and
  iCAF reference programs correspond to a DeSurv factor but to no standard NMF factor. Proposed test (presence):
  restCAF genes among DeSurv defining genes >= 1 AND iCAF genes >= 1 (`stroma_resolution_stats$top_n`,
  code/23); NMF side fixed at 0 and 0. Current 2 of 25, 2 of 35: PASS. The first draft's R29a/R29b rho > 0.2
  tests are dropped (rule 3); the activated-stroma / proCAF half of the sentence is L-D2a and L-D2p.
- **S3a** (S3a) | SI | SI L1145 | The LP remained associated with survival within TCGA-PAAD and CPTAC
  separately. Test: within-cohort unadjusted HR lower bound > 1 in both (`desurv_clinical_adjustment`).
  Current HR 3.50, 3.72: PASS. "strongly" is descriptive.
- **S3b** (S3b) | SI | SI L1145 | The between-cohort score shift is opposite in direction to the between-cohort
  hazard difference. Test: sign(`.ct_shift`) != sign(log `.ct_cohhr`) (`.ct_*`, SI L1128-1140; shift = mean
  LP(TCGA_PAAD) - mean LP(CPTAC), cohort HR = TCGA_PAAD vs CPTAC given the LP). Current -0.73 vs HR 1.90: PASS.
  See D3 for the prose defect.
- **S4** (S4) | SI | SI L1246 | All three Shapley shares were positive. Test: each DeSurv share > 0. Current
  18.8, 25.2, 56.0: PASS.
- **S6b** (S6b, split) | SI | SI L1435-1442 | The exocrine factor N2 is reconstructed less well from the DeSurv
  basis than N1 and N3. Test: R^2(N2 | DeSurv) < R^2(N1 | DeSurv) AND < R^2(N3 | DeSurv) (`.r2_N`, SI L1421).
  Current 0.38 < 0.90, 0.88: PASS. The first draft's "< 0.50" is dropped.
- **S12a** (S12 part) | SI | SI L1640-1646 | Within-platform standardization sensitivity: no significant
  improvement from cohort-specific slopes. Test:
  `sensitivity$lrt_within_platform_standardization$p` >= 0.05. Current: recorded at freeze.
- **S15a** (S15a) | SI | SI L1769, caption L1790 | After clinical adjustment the score remains significant in
  Puleo, PACA-AU array, TCGA-PAAD and CPTAC. Test: `clinical_adjusted$p` < 0.05 in all four. Current 2.1e-4,
  0.025, 4.6e-13, 1.6e-14: PASS.
- **S15b** (S15b) | SI | SI L1769 | Attenuated to non-significance in PACA-AU RNA-seq, Dijk and Moffitt. Test:
  `clinical_adjusted$p` >= 0.05 in all three. Current 0.077, 0.19, 0.15: PASS.
- **S15c** (S15c) | SI | caption L1790 | Moffitt non-significant unadjusted. Test: `unadjusted$p` >= 0.05.
  Current 0.076: PASS.
- **S16** (S16) | SI | SI L1816-1819 | 100 of the 100 runs terminated on the stopping criterion, none failed
  numerically, every run's outer-loop objective decreased monotonically. Test (per rerun, its own seed fits):
  `.n_conv == .R` AND `.n_maxit == 0` AND `!any(.nan)` AND `all(.mono)`, with .R = 100 (NINIT_FULL). Current:
  PASS (Rmd guard).
- **S21a** (S21a) | SI | SI L2217-2245 | In the full-gene view no CAF subtype is resolved as well as the
  activated-stroma axis. Test: DeSurv |`full_rho`| for SCISSORS_restCAF, Elyada_iCAF, SCISSORS_proCAF and
  Elyada_myCAF each < DeSurv `full_rho` for DECODER_ActivatedStroma (`stroma_resolution_stats$full_rho`).
  Current 0.04-0.10 < 0.19: PASS.

### A.3 Composite claims (one claim, test = conjunction of A.1/A.2 tests)

- **A1** (A1) | ABS, MAIN | paper.Rmd L5 | A three-program architecture whose survival-aligned program D1
  combined classical-tumor and restCAF-like features. Test: label D1 AND R6a AND R12. (The first draft cited
  "R5" and "R7b" here; the intended tests were R6a and R12.)
- **R3** (R3) | MAIN | Results L437 | D1 classical + restCAF-like, D2 proCAF-like, D3 basal-like (Fig. 2a).
  Test: labels D1, D2, D3 all hold.
- **R8c** (R8c) | MAIN | Results L446 | Supervision preserved the tumor and stromal structure while resolving D1
  as a program not captured by NMF. Test: R7a AND R7b. ("most of" is descriptive.)
- **R26** (R26) | MAIN | Results L601 | Supervision was not required for a prognostic signal but altered which
  programs were identified. Test: R9a. (The prognostic half depends on the fixed k = 7 control only; the
  first draft's R19c term is now descriptive.)
- **DC1** (DC1) | MAIN | Discussion L9 | Survival-aligned D1 plus separately resolved D2 and D3; D1 prognostic
  after adjustment; programs transferred to external cohorts. Test: A1 AND R3 AND R20a AND R10. ("four
  independent external cohorts" counts datasets, not a fit result.)
- **DC2** (DC2) | MAIN | Discussion L12 | Matched control did not recover D1. Test: R9a. (The "comparable" and
  "similar discrimination" parts are descriptive: Q-R19c, Q-R28a.)
- **DC3** (DC3) | MAIN | Discussion L18 | The risk score combines D1 and D3; D2 remains as proCAF-like context.
  Test: R15 AND label D2.
- **DC5** (DC5) | MAIN | Discussion L27 | PH attenuation; within-cohort analyses argue against the score encoding
  training-cohort membership. Test: R17 AND S3a AND S3b.
- **DC6** (DC6) | MAIN | Discussion L30 | In treated cohorts the frozen scores remained informative. Test: R23a
  AND R24.
- **DC7** (DC7, split) | MAIN | Discussion L33 | GATA6 ISH supports D1's classical tumor component. Test: R4
  AND L-D1 AND L-D1c. (Purity "essentially unchanged" is descriptive: Q-R21b.)
- **DC8** (DC8) | MAIN | Discussion L39 | A three-program architecture that transferred across cohorts. Test:
  A1 AND R10.
- **S5** (S5) | SI | SI L1289-1294 | Delta-loglik concentrates in DeSurv D1. Test: R6a AND R6c.
- **S6c** (S6c) | SI | SI L1435-1442 | DeSurv resolves tumor signal into basal- and classical-aligned programs.
  Test: L-D1 AND L-D1c AND label D3.
- **S9b** (S9 part) | SI | SI L1540-1546 | DeSurv preserves stromal context as an explicit program not used by
  the LP. Test: R15 AND R28c.
- **S24** (S24) | SI | SI L2395, L2403-2409 | Treated-cohort block, "shift toward proCAF-associated bulk
  expression". Test: R15 AND R23a AND R23b AND R23c AND R24 (R24 carries label D2).

---

## Section B. Descriptive quantities

Reported as primary value plus the distribution (min, median, max, all 10 values) over the reruns. No
pass/fail. "Text" marks wording the manuscript should replace with the quantity.

| ID (old) | file:line | Quantity | Object / accessor | Current value |
|---|---|---|---|---|
| Q-ALN (0.1) | plan section 4 | Matched Pearson r, Spearman r, top-ntop Jaccard per factor, rerun vs primary | alignment A | primary = 1 by definition |
| Q-LBL (0.4, R29a, R29b) | Results L437, L440; Fig. 2a | M entries for each factor's row argmax program and for the named programs (ClassicalTumor, restCAF, iCAF, ActivatedStroma, proCAF, BasalTumor, Pure Basal-like) | `make_gene_overlap_heatmap(..., return_matrix = TRUE)$cor_mat` | ClassicalTumor 0.40/-0.40/-0.37; restCAF 0.24/-0.13/-0.09; iCAF 0.24/-0.10/-0.01; ActStroma -0.06/0.37/0.06; proCAF -0.19/0.23/0.04; Basal D3 0.34; Pure Basal D3 0.39 |
| Q-R1a (R1, S17) | Results L426; SI L1851-1862 | Training C across restarts, median and range ("similar"; text) | `program_stability_stats` (code/20), per rerun on its own seed bank | 0.704 (0.678-0.732) |
| Q-R1b (R1, DC4, S17) | Results L426; Discussion L21; SI L1851-1862 | Median pairwise top-ntop Jaccard vs chance ("dispersed", "appreciably different"; text) | same | 0.13 vs 0.074 |
| Q-S17 (S17) | SI L1851-1862 | Median matched gene-loading r among restarts ("moderate"; text) | same | not transcribed in first draft; read at freeze |
| Q-R2 (R2, DC4, S18) | Results L426; SI L1865-1894; Discussion L21 | Split-half consensus agreement ("agreed closely"; text). Per rerun the analog is Q-ALN | `consensus_endpoint_stats` (code/21, primary only) | r 0.97-0.99; Jaccard 0.73-0.88 |
| Q-R4 (R4) | Results L437 | Spearman rho D1 vs GATA6 ISH, n = 33 | `treated_cohort_stats$gata6` | 0.679 |
| Q-R5 (R5, R6d, S4) | Results L443 | Shapley reconstruction shares D1/D2/D3 | `fig_variation_explained_tcgacptac$data`, `.vp` | 18.8 / 25.2 / 56.0 |
| Q-R6 (R6, R6b, S5) | Results L443; SI L1289-1294 | Delta-loglik D1/D2/D3; ratio D3/D1 ("little"; text); max NMF (fixed) | `.dl` | 66.37 / 0.02 / 1.00; NMF 1.16 |
| Q-R7 (R7a, R7b) | Results L446 | max abs Spearman rho with NMF for D1, D2, D3 ("closely", "only weakly"; text) | `.corrC`, Results L146-148 | 0.357 / 0.873 / 0.845 |
| Q-R8 (R8, S8) | Results L446; SI L1435 | R^2(D1 given NMF basis) ("poorly", "essentially unreconstructable"; text) | `projection_r2(W_desurv, W_nmf)[1]`, Results L149 | 0.091 |
| Q-R8b (R8b, S7) | Results L446; SI L1457 | Centered R^2 DeSurv vs NMF and gap in pp ("real but bounded"; text) | `.r2c_*`, Results L132-134 | 0.752 vs 0.812; 6.0 pp |
| Q-R9 (R9a, R9b, R9c, S22) | Results L449; SI L2284-2289 | max abs Spearman rho with matched-budget control for D1, D2, D3 ("recovered", "failed to recover"; text) | `.corrF`, Results L162-166 | 0.245 / 0.861 / 0.749 |
| Q-R10 (R10, A2) | Results L466 | Pooled HR per pooled SD, CI, P | `pooled_cox_stats(pool_df)` | 1.42 (1.26-1.59), P 5.4e-9 |
| Q-R11 (R11, R13) | Results L463 | Per-dataset HR per within-dataset SD, D1/D2/D3, five datasets | Fig. 3a statistic | D1 0.77/0.83/0.69/0.60/0.73; D3 1.27/1.07/0.95/1.26/1.04 |
| Q-R12 (R12, R14a) | Results L463 | Pooled per-program HR, CI, P ("near null" for D2; text) | `.pf`, Results L233-249 | D1 0.75 (0.68-0.83); D2 1.02 (P 0.76); D3 1.09 (0.98-1.21, P 0.10) |
| Q-R14 (R14b) | Results L463 | SD of per-dataset log HR, D1/D2/D3 | from Q-R11 | 0.122 / 0.128 / 0.127 |
| Q-R15 (R15) | Results L466 | beta D1/D2/D3 | fit beta | -0.00100 / 0 / 9.36e-5 |
| Q-R16 (R16, S13) | Results L466; SI L1649-1715 | LODO HR range ("similar"; text) | `.loco`, Results L284-292 | 1.344-1.453 |
| Q-R17 (R17, S11) | Results L466; SI L1621 | PH chi^2, P, slope | `validation_heterogeneity_ci$ph` | 14.67, 1.3e-4, -0.92 |
| Q-R18a (R18a, S12) | Results L469; SI L1640-1646 | Heterogeneity LRT chi^2, P | `validation_heterogeneity_ci` | 1.5 (3 df), 0.69 |
| Q-R18b (R18b, S12) | Results L469 | Q and I^2 ("descriptively homogeneous"; text) | same | Q 1.5, I^2 0% |
| Q-R18c (R18c, S2, S10, S12) | Results L469; SI L1092-1100, L1646 | Event-weighted pooled C and CI ("modest", "moderate"; text) | same | 0.613 (0.579-0.646) |
| Q-R18d (R18d, S10) | Results L469; SI L1576-1599 | Recalibrated IBS vs KM IBS, IPA ("modest"; text) | `desurv_brier_stats` | 0.159 vs 0.168; 5.4% |
| Q-R19a (R19a, A3) | Table hr-adjusted | Adjusted HR, CI, P | `pooled_cox_unadj_adj` | 1.28 (1.11-1.47), P 8.3e-4 |
| Q-R19c (R19c, A5, R26, DC2) | paper.Rmd L5; Results L601; Discussion L12 | DeSurv vs k = 7 control unadjusted HR ("comparable"; text) | same; control fixed | 1.42 (1.26-1.59) vs 1.44 (1.30-1.59) |
| Q-R19d (R19d, S23) | SI L2296 | Per-dataset C, DeSurv and k = 7 control | `val_cindex_*` | see R19d |
| Q-R20a (R20a, S14) | Results L566 | D1 HR after PurIST + DeCAF | `d1_variance_partition_stats$prognosis` | 0.82 (0.73-0.93), P 0.0014 |
| Q-R20b (R20b) | Results L566 | LR chi^2, P; C base vs with D1 | same | 10.1, 0.0015; 0.592 -> 0.622 |
| Q-R20c (R20c, S14) | Results L566; SI L1720-1734 | Variance of D1 explained by PurIST + DeCAF, and by DeCAF alone ("overlaps substantially"; text) | `$variance$explained_total` | 34%; 19% |
| Q-R20d (R20d) | SI L1753-1762 | Interaction-model LR statistic, P | `$interaction` | read at freeze |
| Q-R21a (R21a) | Results L569; SI L1769; Discussion L33 | Spearman rho(D1, purity), PACA-AU array and seq ("only weakly"; text) | `desurv_clinical_adjustment` | 0.06, -0.05 |
| Q-R21b (R21b, DC7, S15) | Results L569; caption L1790 | HR unadjusted vs purity-adjusted, PACA-AU array, seq, TCGA-PAAD ("essentially unchanged"; text) | same | 1.45 vs 1.44; 1.65 vs 1.66; 3.60 vs 3.50 |
| Q-R22 (R22) | Results L569 | Puleo clinical-adjusted HR | `$Puleo_array$clinical_adjusted` | 1.33 (1.14-1.55) |
| Q-R23 (R23a-c, S24) | Results L575; caption L578 | Stratified HR; per-arm HR; interaction P | `treated_cohort_stats$okane` | 0.75 (0.62-0.90); 0.80/0.57/0.83; 0.40 |
| Q-R24 (R24, S24) | Results L575; SI L2409 | Mean paired change in D2, Wilcoxon P | `treated_cohort_stats$paired` | +0.34 SD, 0.013 |
| Q-R25 (R25, S20) | Results L601; SI L2088 | DeSurv k3 coverage and multiplicity ("remained low"; text) | `reference_program_multiplicity$anchors` | 16; 1.00 (vs 32; 1.125 fixed) |
| Q-R27a (R27a, S8) | Results L612; SI L1494-1506 | |score_cor_D1| for SPCA, Cox-PLS, Sparse Cox ("strongly"; text) | `supervised_recovery_stats$table` | 0.81 / 0.83 / 0.85 |
| Q-R27b (R27b) | Results L612 | |r(PC1, D1/D2/D3)| ("only weakly"; text) | same | 0.10 / 0.98 / 0.93 |
| Q-R27c (R27c) | SI L1498 | Supervised scores' max |r| with D2/D3 | same | 0.19 |
| Q-R28a (R28a, DC2, S9) | Results L615; SI L1540-1546; Fig. 3b caption L472 | Pooled external C: DeSurv full LP, SPCA, lasso Cox ("differed less in discrimination"; text) | `desurv_vs_supervised_tuned$prediction[, "POOLED"]` | 0.61 vs 0.60, 0.60 |
| Q-R28b (R28b) | Results L615 | |axis_decomposition D2| for both methods ("only weakly"; text) | `$axis_decomposition` | 0.06; 0.14 |
| Q-R28c (R28c) | Results L615 | Stromal carriage: DeSurv D2 vs supervised | `$stromal_carriage` | 0.55 vs 0.20, 0.23 |
| Q-R29 (R29c) | SI L2224-2227 | restCAF and iCAF gene counts among DeSurv defining genes | `stroma_resolution_stats$top_n` | 2 of 25; 2 of 35 |
| Q-S3 (S3) | SI L1145 | Within-cohort HRs; stratified vs unstratified HR; score shift; cohort HR ("strongly"; text) | `.ct_*`, SI L1128-1140 | 3.50, 3.72; 3.35 vs 3.12; -0.73 SD; 1.90 |
| Q-S6 (S6a, S6b) | SI L1435-1442 | R^2(N1, N2, N3 given DeSurv basis) ("largely reconstructable"; text) | `.r2_N`, SI L1421 | 0.90 / 0.38 / 0.88 |
| Q-S12a (S12a) | SI L1640-1646 | Within-platform standardization LRT P | `sensitivity$lrt_within_platform_standardization$p` | read at freeze |
| Q-S15 (S15) | SI L1769, L1790 | Clinical-adjusted and unadjusted P per cohort | `desurv_clinical_adjustment` | see S15a-c |
| Q-S19 (S19) | SI L1881-1888 | r(W0, W) and top-270 Jaccard per factor ("moved materially"; text) | `zeta_trace$w0` (code/26), per rerun | r 0.57-0.81; Jaccard 0.44-0.82 |
| Q-S21 (S21) | SI L2217-2245 | DeSurv full-gene rho for CAF subtypes and activated stroma | `stroma_resolution_stats$full_rho` | 0.04-0.10; 0.19 |

---

## Section C. Owner rulings on the remaining questions (2026-09-30)

- **C1. Two-program labels: RULED, per-program argmax.** For each named reference program, the aligned factor must
  be the factor with the greatest overlap for that program (column argmax), in addition to L-Dk (row argmax).
  Overlap values reported as descriptive. No cutoff.
- **C2. R14b: RULED, spread of log HR.** SD of per-dataset log HR, D2 greater than both D1 and D3. The manuscript
  text is to name this statistic.
- **C3. R29: RULED, marker presence (R29c).** The named marker genes are in the aligned factor's top-ntop genes;
  the Fig. 2a display filter plays no role.

No owner questions remain; the inventory is FROZEN at the commit that adds this section, before any rerun output.

---

## Section D. Notes

- **D1. Primary k.** The inventory assumes the new primary fit has k = 3. All current values are from the
  fb731e7 fit. If Step 4 selects k != 3, the inventory must be re-frozen by the owner before Step 6. The label
  tests take their reference program from the new primary fit, so the primary label assignments (row argmax
  programs, Q-LBL) are recorded when the primary fit is frozen and before any rerun.
- **D2. Reference gene lists.** `cmbSubtypes_formatted.RData` lives in `../DeSurv-paper/data/derv/` (Amber's
  repo; `code/09a_figures.R:147-151` falls back to it). Pin by checksum (or copy) before any scoring. All label
  tests, R29c, S21a and Q-LBL depend on it.
- **D3. SI L1145 direction wording (text fix, not a claim change).** Verified in R: the TCGA-PAAD mean LP is
  0.73 SD LOWER than CPTAC (means -0.344 vs 0.384), and TCGA-PAAD, not CPTAC, has the higher hazard given the
  predictor (coef cohTCGA_PAAD = 0.64, HR 1.90). The sentence names both cohorts the wrong way round; it prints
  `abs(.ct_shift)`, which hides the sign. The "opposite" claim (S3b) is unaffected. Fix in Step 7.
- **D4. Column position.** Factors are identified by column position everywhere (0.1). The scoring harness
  permutes each rerun's W, H and beta into primary order before running code/05, 09a, 11, 14, 15, 18, 19, 23,
  24, 25. The Rmd guard `beta[2] == 0` (Results L466) and the hard-coded heatmap labels
  (`code/09a_figures.R:161-163`) misreport an unpermuted rerun.
- **D5. Alignment statistic.** The ruling says "Hungarian on W correlation" without naming the coefficient;
  this file keeps Pearson on full W columns (first-draft choice). Matched Spearman r is reported in Q-ALN.
  Change here before the freeze if Spearman is wanted. This is a method setting, not a claim question.
- **D6. program_stability_stats.** `$pairwise` has 435 rows per factor (30 choose 2) while `meta$n_restarts` =
  100 and the text says R = 100 (Q-R1a, Q-R1b, Q-S17). Confirm which restarts code/20 summarizes before reusing
  it per rerun.
- **D7. Claims close at the current fit** (candidates for the 5-8/10 band): R14b (SD 0.128 vs 0.127); L-D2p
  (proCAF 0.23 vs D3 0.04, D1 -0.19); R27c (max |r| with D2/D3 0.19 against 0.81); R6d is not close.
- **D8.** The brief said 88 tracked `.rds`; `git ls-files results | grep -c '\.rds$'` gives 83 at 88cf0f7.
- **D9. Fit-independent, not scored** (unchanged from the first draft): k = 3 selection, alpha and the BO
  surface (paper.Rmd L5, Results L423); cv-grid C-indices (Results L423, SI L2014-2031); the projection-support
  gene count (Results L460); NMF-only non-significance (Results L463); the matched-budget control's own
  non-significance (Results L598); k = 7 coverage and the supervised-curve collapse (Results L601, SI L2102);
  comparators' basal-like gene enrichment, 12/12 seeds (Results L615, SI L1543); 80% subsample single-run
  stability (SI L1904-1924, code/16); method and scaling definitions (SI L1028-1117); control hyperparameters
  (SI L2253-2281); reference universe and nomenclature (SI L1178-1200, L2248-2250); all simulation claims.

### Old-ID crosswalk (every first-draft entry)

A1 -> A1 | A2 -> R10, Q-R10 | A3 -> R19a, Q-R19a | A4 -> R9a, Q-R9 | A5 -> Q-R19c | R1 -> Q-R1a, Q-R1b |
R2 -> Q-R2, Q-ALN | R3 -> R3 | R4 -> R4, Q-R4 | R5 -> R5, Q-R5 | R6 (a, b, c) -> R6a, R6c, Q-R6 | R6d -> R6d,
Q-R5 | R7 (a, b) -> R7a, R7b, Q-R7 | R8 -> Q-R8 | R8b -> R8b, Q-R8b | R8c -> R8c | R9 (a, b, c) -> R9a, R9c,
Q-R9 | R10 -> R10, Q-R10 | R11 -> R11, Q-R11 | R12 -> R12, Q-R12 | R13 -> R13, Q-R11 | R14 (a, b) -> R14a,
R14b, Q-R14 | R15 -> R15, Q-R15 | R16 -> Q-R16 | R17 -> R17, Q-R17 | R18 (a, b, c) -> R18a, Q-R18a, Q-R18b,
Q-R18c | R18d -> R18d, Q-R18d | R19 (a, b) -> R19a | R19c -> Q-R19c | R19d -> R19d, Q-R19d | R19e -> dropped |
R20 (a, b, c, d) -> R20a, R20b, R20d, Q-R20a-d | R21 (a, b) -> Q-R21a, Q-R21b | R22 -> R22, Q-R22 | R23 (a, b,
c) -> R23a-c, Q-R23 | R24 -> R24, Q-R24 | R25 -> R25, Q-R25 | R26 -> R26 | R27 (a, b, c) -> Q-R27a, R27b,
R27c, Q-R27b, Q-R27c | R28 (a, b, c, d) -> Q-R28a, Q-R28b, R28c, R28d, Q-R28c | R29 (a, b, c) -> R29c, L-D1r,
L-D2a, L-D2p, Q-LBL, Q-R29 | DC1, DC2, DC3, DC5, DC6, DC7, DC8 -> same IDs | DC4 -> Q-R1b, Q-R2 | S1 -> R15 |
S2 -> Q-R18c | S3 (a, b) -> S3a, S3b, Q-S3 | S3c -> dropped | S4 -> S4 | S5 -> S5 | S6 (a, b, c) -> Q-S6, S6b,
S6c | S7 -> R8b | S8 -> R27c, R29c, Q-R8, Q-R27a | S9 -> R28c, R28d, S9b | S10 -> R18d, Q-R18c | S11 -> R17 |
S12 (+ a) -> R18a, S12a, Q-R18a-c, Q-S12a | S13 -> Q-R16 | S14 -> R20a-d | S15 (a, b, c) -> S15a-c, Q-S15 |
S16 -> S16 | S17 -> Q-R1a, Q-R1b, Q-S17 | S18 -> Q-R2 | S19 -> Q-S19 | S20 -> R25 | S21 (+ a) -> S21a, Q-S21 |
S22 -> R9a, R9c | S23 -> R19d, R19a | S24 -> S24 | 0.4 L-D1..L-D3 -> L-D1..L-D3, named-component tests

---

## Counts

- Binary claims: 70 tests. Label tests 9 (L-D1, L-D2, L-D3 under rule 3; 6 named-component tests pending C1).
  Elementary binary claims 46. Composite claims 15 (conjunctions of the above).
- Descriptive quantities: 50 (Section B rows).
- Owner questions: 3 (C1 two-program labels, C2 R14b comparator, C3 R29 "corresponded").
- Dropped:
  - R19e: involves only fixed controls (unchanged from the first draft).
  - S3c: the named cohort direction is a prose defect in SI L1145 (D3), not a claim; the claim is S3b.
  - Threshold tests replaced by descriptive quantities under rule 1 (content kept in Section B): R1, R2, R6b,
    R7a threshold part, R7b threshold part, R8, R8b bound, R9b, R16, R18b, R18c, R18d bound, R19c, R20c, R21a,
    R21b, R25 multiplicity, R27a, R28a, R28b, S6a, S6b "< 0.50", S17, S19.
  - Label-test cutoffs and margins (rho > 0.2, iCAF fallback) and R29a/R29b rho tests under rule 3.
