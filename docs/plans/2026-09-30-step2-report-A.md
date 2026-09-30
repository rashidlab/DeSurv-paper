# Rev 3 Step 2, Agent A report (2026-09-30): section 5 pre-run fixes

Worktree `~/Downloads/DeSurv-paper-recompute`, branch `recompute/consensus-desurv`, HEAD 88cf0f7, all changes
uncommitted. No `code/NN_*.R`, `make`, `sbatch` or `cache_or_compute` was run. `results/` and `figures/`:
`git status` clean, `find -newer <marker>` empty. Canonical repo untouched (read only).

## Acceptance

- Parse: all 29 edited or new `.R` files `parse()` OK; the 3 edited child/SI Rmds `knitr::purl` + `parse` OK;
  the two new inline expressions in `paper.Rmd` evaluate (`4.6.0`, `1.0.1` under the site library).
- `bash -n`: all 10 edited or new `.sh` files OK. `rev3_env.sh` exercised locally with a stubbed `module`:
  positive case passes with the local rev 3 lib; it refuses `DESURV_RESULTS_DIR=results`, an unset
  `DESURV_RESULTS_DIR`, and a library whose DeSurv is 1.0.1.
- Paper tests (`tests/`): 2 files / 5 blocks / 68 expectations / 0 failed / 0 errors / 0 skipped, identical to
  step 1. Library: the step 1 private rev 3 lib
  `/tmp/claude-1000/-home-naimrashid-Downloads-DeSurv-paper-clean/da9c6c81-47b4-4511-8cfc-1c71451b090a/scratchpad/lib-rev3`
  (DeSurv 1.2.0) plus `~/R/library` for testthat, run from the repo root. The same run started from `~`
  gives 1 error ("no successful seed fits") at HEAD too (pristine `git archive HEAD` copy), so it is a
  library-path artefact, not these edits. Note: the repo's tracked, empty `.Renviron` means `~/.Renviron` is
  not read when R starts in the repo root.
- Build-graph guard (`scripts/check_build_graph.sh` writer regex, emulated without `make`): the writer map for
  all 38 manuscript cache names is identical before and after.

## Item 1. Output isolation

New `R/paths.R` (side-effect free): `RESULTS_BASE` = `DESURV_RESULTS_DIR` (default `"results"`),
`FIGURE_BASE` = `DESURV_FIGURE_DIR` (default `"figures"`), `res_path()`, `fig_path()`.

- `code/00_helpers.R:28-30` before: `RESULTS_DIR <- if (quick) "results/quick" else "results"`,
  `CV_GRID_DIR <- "results/cv_grid"`, `FIGURE_DIR <- if (quick) "figures/quick" else "figures"`.
  After (now 29-34): `source("R/paths.R")`, then the same three built from `RESULTS_BASE`/`FIGURE_BASE`.
  Defaults are byte-identical strings (checked).
- `paper/load_precomputed.R:8`: `PRECOMPUTED_DIR` takes `DESURV_RESULTS_DIR` when set, else the old
  `results` / `../results` rule; adds `CV_GRID_CSV_DIR` and `paper_fig(name, default)`, where `default`
  is the call site's old prefix (`"../figures"` or `"figures"`, which differ by site).
- Standalone scripts (not sourcing 00_helpers) get `source("R/paths.R")` before first use.

Site counts (grep for `"results/`, `file.path("results"`, `"figures/`, `file.path("..", "figures"` etc.):

| Scope | Before (lines) | Routed | Left (why) |
|---|---|---|---|
| `code/`, `R/` | 86 | 65 path uses in 20 scripts + 3 in 00_helpers | 26 lines: 6 in code/13, 13c, 13d (see below); 7 provenance labels or messages, not paths (24:90, 24:205-207, 24:217, 29:88, 16:167); 13 comments and the R/paths.R defaults |
| `paper/*.Rmd` | 16 | 15 (3 cv_grid CSV reads, 12 figure embeds) | 2 `fig.path = "figures/"` (paper.Rmd:83, si:1382): knitr's own chunk-output dir, left |
| `paper/load_precomputed.R` | 1 | 1 | |

Routed scripts: 09a (2), 11 (2), 14 (3), 15 (3), 16 (8), 17 (2), 18 (4), 19 (5), 20 (4), 21 (5), 22 (2),
23 (1), 25 (1), 26 (1), 27 (3), 28 (5), issue7 (3), util_fig3_standalone_preview (6),
util_fig_supervised_recovery (2), util_fig_treated_context (3). 09a:406-408 comment updated (reads the
base dir via `res_path`, as before it read the fixed base `results/`).

Not routed, for the coordinator:
- **code/13, 13c, 13d (withdrawn spatial)**: 13c:17 and 13d:12 `setwd("/home/naimrashid/Downloads/DeSurv-paper-clean")`,
  13:40-42 `setwd(ROOT)` to the same canonical repo. Run from the worktree they would write into the
  canonical `results/`. Not in any runner; `run_rev3_step.sh` refuses `13*`.
- `Makefile` `clean` (127-135) hardcodes `results`/`figures/quick` (not in the grep scope; left).
- `slurm/run_fair_comparator.sh` preflight now also requires `RESULTS_BASE == "results"`.
- Static figures the Rmds embed (e.g. `model_schematic_final.pdf`, and Fig 1/Fig 4 assets built out of tree)
  must be copied into `DESURV_FIGURE_DIR` before any render that sets it.

## Item 2. code/26 version stop

`code/26_zeta_trace.R:45` (now 46-47) `stopifnot(identical(as.character(packageVersion("DeSurv")), "1.0.1"))`
-> `stopifnot("DeSurv >= 1.2.0 required (rev 3)" = packageVersion("DeSurv") >= "1.2.0")`. Also 133 (now 134):
`meta$package` was the typed `"DeSurv 1.0.1 (afb00d5)"`, now `paste0("DeSurv ", packageVersion("DeSurv"))`.
Left: header 14-22 still describes the v1.0.1 C++; the script's own reconstruction assertions will catch a
changed W update under 1.2.0.

## Item 3. code/01 install source

`code/01_install.R:9-27` replaced: `DESURV_PKG_SOURCE` is required (local directory, or GitHub spec
`owner/repo@ref`); unset, not a dir and not a spec, or resolving to `../DeSurv` -> `stop()`. Local installs
use `force = TRUE` and log the git HEAD. After `library(DeSurv)`, the version must be `>= DESURV_MIN_VERSION`
(default `1.2.0`). Stale comment `results/precomputed/` -> `results/`.
Consequence: `make install`, `make quick`, `make main`, `make all` (all start with 01) now fail loudly unless
`DESURV_PKG_SOURCE` is set. Intended; the Slurm runners no longer call make.

## Item 4. paper.Rmd versions

`paper/paper.Rmd:150` before `R (version 4.6.0; ...) using DeSurv (v1.0.1; this work)`; after
`` `r paste(R.version$major, R.version$minor, sep = ".")` `` and
`` `r as.character(utils::packageVersion("DeSurv"))` ``. Both report the RENDER session: rendering locally gives
R 4.6.0, while the rerun computes on Longleaf R 4.5.0; render with the rev 3 library first on `.libPaths()`.
Same line still types `NMF (v0.28)`, `glmnet (v5.0)`, `survival (v3.8.6)`, `ggplot2 (v4.0.3)`,
`cowplot (v1.2.0)` (not changed).
Report only: `paper.Rmd:15` `tar_store: "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main"` is used nowhere
(no `params$tar_store` reference); the only `targets::tar_read` is in chunk `versioning` (paper.Rmd:174),
`eval=FALSE`, which also reads `../DeSurv` git state. `_targets.yaml` is
`store: ../DeSurv-paper/store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main`, referenced by no script in
`code/ R/ scripts/ slurm/ Makefile`. Both are dead config.

## Item 5. Slurm coverage

New:
- `slurm/rev3_env.sh`: shared preamble sourced by every runner: `module purge; module load r/4.5.0`;
  `REV3_LIB=/work/users/n/u/nur2/desurv-rerun/lib/r4.5-rev3` exported as both `R_LIBS` and `R_LIBS_USER`;
  `mkdir -p slurm/logs`; `DESURV_RESULTS_DIR` and `DESURV_FIGURE_DIR` required, made absolute, refused if
  inside tracked `results/`/`figures/`; `DESURV_RECOMPUTE=TRUE` explicit; `DESURV_NCORES` from the allocation;
  an Rscript check that DeSurv >= 1.2.0 AND loads from `REV3_LIB`, before any work.
- `slurm/run_rev3_step.sh <script>...`: generic runner; refuses 01, 13*, 21 (comment: 21 is replaced by the
  section 4 sensitivity analysis, no runner).
- `slurm/submit_rev3_downstream.sh`: one job per script for 14, 15, 16, 17, 18, 20, 22, 23, 24, 25, 26, 29 with
  per-script resources and `afterok` on `MAIN_JOB` (and `CV_JOB` for 16 and 24, which read the CV grid).

Existing runners fixed (all had `module load r/4.4.0`, no rev 3 lib, no version check):
- `run_main_pipeline.sh`: `make main` -> direct `Rscript` over `DESURV_STEPS` (default `02 03 04 05 08 19`),
  never 01; refuses a non-empty results dir for a stage containing 02/03 (clean run); 09a dropped from the
  default because it reads code/11's `treated_cohort_stats.rds`. Walltime 18 h -> 72 h placeholder.
- `run_cv_grid.sh`: `make cv-grid` -> `Rscript code/06_cv_grid.R`; 6 h -> 72 h (see item 6).
- `run_simulations.sh`, `run_merge_simulations.sh`: preamble; merge renders only with `REV3_RENDER=1`
  (render writes tracked `paper/*.pdf`, needs LaTeX and every static figure in `DESURV_FIGURE_DIR`).
- `submit_all.sh`: rewritten as the post-BO chain main(04 05 08 19) -> 06 and 07 after main -> merge; the old
  version launched 06 in parallel with 03 although 06 reads the BO's `tar_params_best`.
- `run_fair_comparator.sh`: marked pre-rev-3 legacy (local nohup, pinned 1.0.1, tracked `results/`).
- `code/16_desurv_vs_supervised.R:37` unconditionally prepended `~/R/x86_64-pc-linux-gnu-library/4.5` to
  `.libPaths()`, which could shadow the rev 3 lib in-process after the preamble's check. Now skipped when
  `REV3_LIB` is set. Same pattern remains in 21:38, 13c:14, 13d:10 (not run).

Brief mismatch: there are no existing per-step runners for 03, 04 or 05; 03-05 ran only inside `make main`
(and `run_fair_comparator.sh`).

Walltimes. Sized from script headers via the inventory: 03 4-8 h at 30 starts, 04 2-4 h, 06 24-48 h, 07 2 h per
task. **Unsized**: the main stage (BO now scores 100 seed fits per fold; set from the step 3 benchmark), and
every downstream job (14 1 h, 15 2 h, 16 24 h, 17 1 h, 18 1 h, 20 2 h, 22 1 h, 23 1 h, 24 2 h, 25 2 h,
26 24 h, 29 24 h/16 cores). These are guesses with margin; 16, 26 and 29 are the likeliest to be wrong.

code/11 inputs absent on Longleaf (all under `~/research/rif-research/`, restricted):
| File | Size | md5 |
|---|---|---|
| `purist-twist/data/external/response_master_canonical.rds` | 78,686,880 B | 5d07017711198e05e02a9b40814f5b0f |
| `shared/data_paths.R` (defines `COHORT`) | 30,620 B | |
| `data/longleaf/public/Processed_datasets/PDAC_RNA_as_is/Linehan.rds` (`COHORT$Linehan`, the only cohort 11 reads, line 68) | 8,815,474 B | 159d82aa8274d30de1cebeb0ec2c97bd |
Proposal: run 11 locally (restricted data stays local; the script is light) after copying the rerun's
`tar_fit_desurv_tcgacptac.rds` into a local isolated dir, with the same `DESURV_RESULTS_DIR`; copy
`treated_cohort_stats.rds` back before 09a. `Linehan.rds` mirrors `/proj/jjyehlab/data/public/...` on
Longleaf, but `response_master_canonical.rds` and `data_paths.R` would still need copying. Also: 11:46 hardcodes
`get_top_genes(fit$W, 270)` and D1-D3 (k = 3).

## Item 6. CV-grid walltime

code/06 has no chunking hook (no env/args; `params_grid` built at 128-134 and processed whole by `mclapply` in
three stages: CV 149-167, full fits 224-256, validation 276-282; each stage caches one list rds). Not
implemented. Proposed split, consistent with how 06 iterates: every grid point uses the fixed `SEED = 123L`,
so points are independent and can be partitioned without changing results. Add `DESURV_CVGRID_CHUNK=i/n`
(or `SLURM_ARRAY_TASK_ID`) that subsets `params_grid` by k (11 tasks, 40 points each) and writes
`cv_grid_result_list_chunk_<i>.rds` and `cv_grid_fit_list_chunk_<i>.rds`; a merge job reassembles both lists
in `params_grid` order and runs the cutpoint, summary and validation stages unchanged. Per-k cost rises with k,
so array by k with 32 cores each needs roughly 4-8 h per task (from the 24-48 h x 32-core estimate), 12 h
walltime. Until then `run_cv_grid.sh` is 72 h unsplit. Inventory mismatch: 06 no longer types lambda/nu/ntop
(55-56 read `tar_params_best`).

## Item 7. cmbSubtypes file

`R/paths.R` adds `cmb_subtypes_path()`: `DESURV_CMBSUBTYPES` if set, else `data/derv/cmbSubtypes_formatted.RData`;
`stop()` if missing, no fallback. Replaced: `code/09a_figures.R:147-155` (9 lines -> 1),
`code/23_stroma_resolution.R:32-34` (3 -> 1), `code/24_reference_program_multiplicity.R:79-86` (8 -> 1).
`grep "DeSurv-paper/data"` in code/ R/: 0 hits outside 13c/13d (`Elyada_umap.Rds`, withdrawn).
File: `/home/naimrashid/Downloads/DeSurv-paper/data/derv/cmbSubtypes_formatted.RData`, 87,698 B,
md5 `ee4c81468936b0a6b514fe4b4e41fe65`. `data/derv/` does not exist in either repo. Not committed.

## Item 8. gene_lists_top270_k3.csv producer

New `R/gene_lists.R` (`build_gene_lists(W_desurv, W_std, ntop)`, pure) and `code/util_gene_lists.R` (reads
`tar_fit_desurv`, `fit_std_desurvk`, `tar_params_best` via `res_path`; no `cache_or_compute`; writes
`gene_lists_top<ntop>_k<k>.csv` to `DESURV_GENE_LIST_DIR`, default `paper/`). Columns: rank, gene,
weight = `W[gene, j]`, contrast_score = `desurv_get_top_genes()$diffs`, factor D1..Dk then N1..Nk, rounded to 6.
Verification (site DeSurv 1.0.1, current cached objects, output to scratch): 1620 rows x 5 as committed;
**header + all 810 D rows byte-identical**; the 810 N rows differ: N-factor gene overlap with the current
`fit_std_desurvk` is 267/270/267 (N1/N2/N3), so the committed N rows came from an earlier std-NMF fit (the CSV
dates from 143ca80, 2026-04-29). The SI reads only D1-D3. Script and function outputs are byte-identical.

SI lines hardcoding 270 or k = 3 (prose not edited): si_appendix.Rmd 2014 and 2020 (`n_top = 270`; 2020 also
types lambda 0.34), 2085 (`n_top = 270`), 2303 (tables `gene-list-270-D1..D3`, "top 270 genes", "three DeSurv
factors"), 2317 (chunk label), 2323 (`paper/gene_lists_top270_k3.csv`), 2327 and 2339 (`D1, D2, D3`), 2329-2334
(45-row split to 270), 2340-2350 (factor captions and descriptions), 2354-2356, 2359-2363 (labels, col.names,
caption "Top 270", "The 270-gene cutoff").

## Item 9. Seeding

- `code/03_bayesian_optimization.R` (brief says 183-187; the call is at 227-229 in HEAD, 230-233 now): added
  `STD_NMF_SEED <- 20261201L` and `seed = STD_NMF_SEED` to `NMF::nmf(..., nrun = STD_NMF_NRUN, .options = "p<nrun>")`.
  `?NMF::nmf` (0.28): a single numeric `seed` drives L'Ecuyer streams across runs so "parallel computation
  are fully reproducible". Checked on a 60 x 20 matrix, rank 2:3, nrun 4, `.options = "p4"`: seeded runs
  identical (residuals and k = 3 basis); unseeded runs differ.
- `R/fit_cox_model.R:1-17`: `fit_cox_model(XtW, df, nfold = 5, fold_seed = 20261202L)`: `foldid` drawn once
  (`set.seed(fold_seed)`, `sample(rep(seq_len(nfold), length.out = n))`) and passed to every `cv.glmnet` over
  the 21 nu values (replacing `nfolds = nfold`); caller's `.Random.seed` saved and restored. Checked:
  repeated calls identical, caller RNG unchanged, a different `fold_seed` changes the result. Callers
  (04:227, 04:241, cv_grid_helpers:313) unchanged; positional `nfolds` still binds.
- Not in the brief: `R/cv_grid_helpers.R:283` has the same unseeded parallel `NMF::nmf` (std-NMF CV used by
  code/08:127). Not changed.

## Item 10. preprocessCore

Paper repo (`code/ R/ paper/*.R`): 0 calls to `preprocessCore` or `normalize.quantiles*`. Package
(`DeSurv-clean`): `R/preprocess_data.R:229-251`, `R/cv_helpers.R:130-138`, `R/predict_methods.R:114-125`, each
only in the `"quant"` branch. Steps 2/5/6 pass `method_trans_train = "rank"` (03:91, 178, 220, 271; 05:46;
06:96), so the quantile path is not on the pipeline. Setup report: preprocessCore 1.70.0 already in the
Longleaf r/4.5.0 system library. No fix needed.

## Item 11. setup_env.sh test invocation

`analysis/longleaf/setup_env.sh` step "testthat on installed package" ran `library(DeSurv)` + `test_file()` on
`test-cox-ties.R` and `test-consensus-fit.R`. Broken: tests ran outside the namespace (internal
`cox_stats_cpp` not found) and `local_mocked_bindings()` falls back to `pkgload::dev_package()`; also
`test-consensus-fit.R` no longer exists in 1.2.0. Reproduced locally: old invocation on test-cox-ties 5/5
errors. Fix: `NOT_CRAN=true` + `testthat::test_dir("$SRC/tests/testthat", package = "DeSurv", load_package = "installed")`,
whole suite, asserts DeSurv loads from `$LIB`, prints counts, exits non-zero on failure. Locally against rev 3
1.2.0 (1f3d3f0): 18 files / 124 blocks / 458 expectations / 0 / 0 / 0 skipped, identical to step 1's devtools
counts (without `NOT_CRAN`, 14 blocks skip). Heredoc expansion checked (`$SRC` expands, `\$` survives, parse OK).
Left: the rest of that script is pre-rev-3 (`LIB=$ROOT/lib/r4.5`, source `feat/consensus-estimator @ b42d7a9`,
preflight `== "1.1.0"` and `desurv_consensus_fit`/`desurv_consensus_gate`, which 1.2.0 removed). Rerunning it
as is fails at preflight; it needs the rev 3 lib path, source and checks.

## Item 12. CLAUDE.md

Mismatch: `CLAUDE.md` is gitignored (`.gitignore:23`) and absent from the worktree; the only copy is the
canonical repo's untracked file, which I may not edit. I wrote a corrected copy at the worktree root
(`CLAUDE.md`, ignored, so it will not appear in a commit). 8 edits vs canonical: `run_full_pipeline.sh` ->
`run_main_pipeline.sh` (17, and the Slurm paragraph 91 rewritten for the rev 3 runners); `run_pipeline.R`
commands (23-25) replaced by direct `Rscript` examples plus the isolated-rerun form; env-var list (47) now
includes `DESURV_RESULTS_DIR`, `DESURV_FIGURE_DIR`, `DESURV_CMBSUBTYPES`, `DESURV_PKG_SOURCE`; path block (55-57)
shows `R/paths.R`; "60 objects" -> 79 top-level `.rds` (88 tracked files under `results/`: 83 `.rds` + 5 `.csv`);
`make install` needs `DESURV_PKG_SOURCE`; `code/09_render_paper.R` -> `code/10_render_paper.R`. DANGER note kept.
To apply to canonical: `cp` the worktree file over `~/Downloads/DeSurv-paper-clean/CLAUDE.md` (canonical-side
text such as the note that code/13c/13d `setwd()` there also applies). Still stale, not in the brief: "There is
no test suite" (29; `tests/` has 2 files), `make quick` "full pipeline" (15; it runs 01-05).

## Files changed

Modified (37): `R/fit_cox_model.R`, `analysis/longleaf/setup_env.sh`, `code/00_helpers.R`, `code/01_install.R`,
`code/03_bayesian_optimization.R`, `code/09a_figures.R`, `code/11_treated_cohort_analysis.R`,
`code/14_d1_variance_partition.R`, `code/15_supervised_recovery.R`, `code/16_desurv_vs_supervised.R`,
`code/17_brier_metrics.R`, `code/18_clinical_adjustment.R`, `code/19_desurv_vs_supervised_tuned.R`,
`code/20_program_stability.R`, `code/21_consensus_endpoint.R`, `code/22_fair_comparator_gate.R`,
`code/23_stroma_resolution.R`, `code/24_reference_program_multiplicity.R`,
`code/25_validation_heterogeneity_ci.R`, `code/26_zeta_trace.R`, `code/27_program_specificity.R`,
`code/28_tuned_predictor_decomposition.R`, `code/issue7_fig2_reconstruction.R`,
`code/util_fig3_standalone_preview.R`, `code/util_fig_supervised_recovery.R`, `code/util_fig_treated_context.R`,
`paper/03_methods_REVISED.Rmd`, `paper/04_results_REVISED.Rmd`, `paper/load_precomputed.R`, `paper/paper.Rmd`,
`paper/si_appendix.Rmd`, `slurm/run_cv_grid.sh`, `slurm/run_fair_comparator.sh`, `slurm/run_main_pipeline.sh`,
`slurm/run_merge_simulations.sh`, `slurm/run_simulations.sh`, `slurm/submit_all.sh`.

New (6, untracked): `R/paths.R`, `R/gene_lists.R`, `code/util_gene_lists.R`, `slurm/rev3_env.sh`,
`slurm/run_rev3_step.sh`, `slurm/submit_rev3_downstream.sh`.

New, gitignored: `CLAUDE.md` (worktree copy, item 12). This report.
