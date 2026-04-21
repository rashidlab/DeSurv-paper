# Remaining Steps to Generate Main Paper Figures

Status of each figure needed by `paper/04_results_REVISED.Rmd` (main paper).
This document does **not** cover supplement figures (`paper/si_appendix.Rmd`).

## Completed (in code/08_figures.R)

| Figure | Function | Status |
|--------|----------|--------|
| `fig_residuals_tcgacptac` | `make_nmf_metric_plot(fit_std, "residuals")` | Done |
| `fig_cophenetic_tcgacptac` | `make_nmf_metric_plot(fit_std, "cophenetic")` | Done |
| `fig_silhouette_tcgacptac` | `make_nmf_metric_plot(fit_std, "silhouette")` | Done |
| `fig_bo_heat_tcgacptac` | `extract_gp_curve()` + `ggplot geom_tile` | Done |
| `fig_variation_explained_tcgacptac` | `build_variance_survival_df()` + `ggplot` | Done |
| `fig_gene_overlap_heatmap_desurv_tcgacptac` | `make_gene_overlap_heatmap()` | Done |
| `fig_gene_overlap_heatmap_std_desurvk_tcgacptac` | `make_gene_overlap_heatmap()` | Done |
| `fig_gene_overlap_heatmap_std_elbowk_tcgacptac` | `make_gene_overlap_heatmap()` | Done |
| `fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac` | `make_gene_overlap_heatmap()` | Done |
| `fig_hr_forest_tcgacptac` | `compute_hrs()` + `ggplot` | Done |
| `fig_desurv_std_correlation_tcgacptac` | `pheatmap::pheatmap()` on Spearman corr | Done |

## Remaining: KM Survival Plots

### Figures needed

| Figure | Function | Blocking dependency |
|--------|----------|---------------------|
| `fig_median_survival_desurv_tcgacptac` | `splot_cutpoint(data_val_filtered, tar_fit_desurv, desurv_lp_stats, ntop)` | `desurv_lp_stats` |
| `fig_median_survival_std_desurvk_tcgacptac` | `splot_cutpoint(data_val_filtered, fit_std_desurvk, std_desurvk_lp_stats, NULL)` | `std_desurvk_lp_stats` |

### What `lp_stats` contains

Each `lp_stats` object is a list:
```r
list(
  lp_mean = <numeric>,                  # mean of training LP scores
  lp_sd = <numeric>,                    # sd of training LP scores
  optimal_z_cutpoint = <numeric>,       # best z by log-rank
  cutpoint_abs = <numeric>,             # z * sd + mean
  optimal_z_cutpoint_cindex = <numeric>,# best z by c-index
  cutpoint_abs_cindex = <numeric>       # z * sd + mean
)
```

### Computation chain

The chain to produce `lp_stats` has an expensive step that must run on HPC:

#### Step A: CV cutpoint search (EXPENSIVE — needs HPC)

For DeSurv:
```r
desurv_cv_cutpoint_result <- run_cv_grid_point(
  data = tar_data_filtered,
  k = tar_params_best$k,
  alpha = tar_params_best$alpha,
  fixed_params = list(
    lambda = tar_params_best$lambda,
    nu = tar_params_best$nu,
    lambdaW = 0,  # from bo_config defaults
    lambdaH = 0,
    ntop = ntop_value
  ),
  nfolds = 5,
  n_starts = 30,
  seed = 123
)
```

For std NMF at DeSurv k:
```r
std_desurvk_cv_cutpoint_result <- run_cv_grid_point_std_nmf(
  data = tar_data_filtered,
  k = tar_params_best$k,
  fixed_params = list(ntop = NULL),
  nfolds = 5,
  n_starts = 30,
  seed = 123
)
```

This is 5 folds x 30 starts = 150 model fits per method. Runs in ~10-30 min on HPC with 30 cores.

#### Step B: Evaluate cutpoints (cheap — can run locally)

```r
desurv_cutpoint_eval <- evaluate_cutpoint_zscores(
  desurv_cv_cutpoint_result,
  z_grid = seq(-2.0, 2.0, by = 0.2)
)

desurv_cutpoint_summary <- desurv_cutpoint_eval |>
  dplyr::group_by(z_cutpoint) |>
  dplyr::summarise(
    mean_cindex_dichot = mean(cindex_dichot, na.rm = TRUE),
    mean_abs_logrank_z = mean(abs(logrank_z), na.rm = TRUE),
    .groups = "drop"
  )

desurv_optimal_z_cutpoint <- desurv_cutpoint_summary |>
  dplyr::slice_max(mean_abs_logrank_z, n = 1, with_ties = FALSE) |>
  dplyr::pull(z_cutpoint)

desurv_optimal_z_cutpoint_cindex <- desurv_cutpoint_summary |>
  dplyr::slice_max(mean_cindex_dichot, n = 1, with_ties = FALSE) |>
  dplyr::pull(z_cutpoint)
```

Same pattern for `std_desurvk_*`.

#### Step C: Assemble lp_stats (cheap)

```r
lp <- compute_lp(
  tar_fit_desurv$W, tar_fit_desurv$beta,
  tar_data_filtered$ex, ntop_value
)
desurv_lp_stats <- list(
  lp_mean = mean(lp),
  lp_sd = sd(lp),
  optimal_z_cutpoint = desurv_optimal_z_cutpoint,
  cutpoint_abs = desurv_optimal_z_cutpoint * sd(lp) + mean(lp),
  optimal_z_cutpoint_cindex = desurv_optimal_z_cutpoint_cindex,
  cutpoint_abs_cindex = desurv_optimal_z_cutpoint_cindex * sd(lp) + mean(lp)
)
```

For std NMF, `compute_lp` is called with `ntop = NULL`.

### Implementation plan

1. **Create `code/04b_cutpoint_selection.R`** — new pipeline step that:
   - Loads `tar_data_filtered`, `tar_params_best`, `tar_fit_desurv`, `fit_std_desurvk` from precomputed
   - Runs `run_cv_grid_point()` for DeSurv (with ntop) and std NMF (without ntop)
   - Evaluates cutpoints, picks optimal z-cutpoints
   - Assembles and saves `desurv_lp_stats_tcgacptac.rds` and `std_desurvk_lp_stats_tcgacptac.rds`
   - Runtime: ~10-30 min on HPC

2. **Add to `slurm/run_ntop_pipeline.sh`** — insert after step 4, before step 5:
   ```bash
   Rscript code/04b_cutpoint_selection.R
   ```

3. **Add KM figures to `code/08_figures.R`**:
   ```r
   desurv_lp_stats      <- load_precomputed("desurv_lp_stats_tcgacptac")
   std_desurvk_lp_stats <- load_precomputed("std_desurvk_lp_stats_tcgacptac")

   fig_median_survival_desurv <- cache_or_compute("fig_median_survival_desurv_tcgacptac", {
     splot_cutpoint(data_val_filtered, tar_fit_desurv, desurv_lp_stats, ntop_value)
   })

   fig_median_survival_std_desurvk <- cache_or_compute("fig_median_survival_std_desurvk_tcgacptac", {
     splot_cutpoint(data_val_filtered, fit_std_desurvk, std_desurvk_lp_stats, ntop = NULL)
   })
   ```

## Remaining: Simulation Figures

| Figure | Status | Notes |
|--------|--------|-------|
| `sim_figs_by_scenario` | Blocked | 6-12 hour HPC run (step 7). Not ntop-dependent. Existing base-directory `.rds` is an LFS pointer — fix with `git lfs pull`, or rerun step 7. |

## Remaining: Non-figure objects loaded by paper Rmd

These are data objects (not figures) that the paper Rmd loads via `load_result()`. Most are produced by steps 3-5 and should already exist in `results/precomputed/ntop_150/`:

| Object | Produced by | In ntop_150? |
|--------|-------------|--------------|
| `desurv_seed_fits_tcgacptac` | Step 4 | Yes |
| `val_cindex_desurv_tcgacptac` | Step 5 | Yes |
| `val_cindex_std_desurvk_tcgacptac` | Step 5 | Yes |
| `val_cindex_std_elbowk_tcgacptac` | Step 5 | Yes |
| `val_cindex_desurv_alpha0_tcgacptac` | Step 5 | Yes |
| `val_latent_desurv_tcgacptac` | Step 5 | Yes |
| `val_latent_std_desurvk_tcgacptac` | Step 5 | Check |
| `val_latent_desurv_alpha0_tcgacptac` | Step 5 | Yes |
| `data_val_filtered_tcgacptac` | Step 5 | Yes |
| `fit_std_elbowk_tcgacptac` | Step 4 | Yes |
| `tar_data_filtered_elbowk_tcgacptac` | Step 3 | Check (may be base only) |

## Reference data dependency

The gene overlap heatmaps require `data/derv/cmbSubtypes_formatted.RData` (contains `top_genes` reference gene lists). This file does not exist in DeSurv-paper-rl — the script falls back to `../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData`. For standalone use, copy this file into the repo.
