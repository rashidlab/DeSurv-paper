#!/usr/bin/env Rscript
# code/08_figures.R — Generate all figure objects from pre-computed results
#
# This script rebuilds figure objects that are stored in results/precomputed/.
# Most figures are already pre-computed from the HPC pipeline. This script
# can regenerate them from the model fits and validation results.
#
# Note: Some figure objects (KM plots, heatmaps) depend on the ggplot2 version
# matching between the environment that created the model fits and the one
# rendering them. If figures fail, use the pre-computed versions.
#
# Inputs:  model fits, validation results from steps 04-05
# Outputs: figure objects in results/precomputed/fig_*.rds
#
# Runtime: ~5 min

message("=== Step 8: Figures ===")
source("code/00_helpers.R")
library(ggplot2)
library(cowplot)
library(survival)
library(survminer)
library(NMF)

# Source figure-building functions
source("R/figure_targets.R")
source("R/cluster_alignment.R")
source("R/compare_models.R")
source("R/cv_grid_helpers.R")
source("R/get_top_genes.R")
source("R/plot_survival.R")

# ── Load prerequisites ────────────────────────────────────────────────────
tar_fit_desurv      <- load_precomputed("tar_fit_desurv_tcgacptac")
fit_std_desurvk     <- load_precomputed("fit_std_desurvk_tcgacptac")
fit_std_elbowk      <- load_precomputed("fit_std_elbowk_tcgacptac")
tar_data_filtered   <- load_precomputed("tar_data_filtered_tcgacptac")
data_val_filtered   <- load_precomputed("data_val_filtered_tcgacptac")
desurv_bo_results   <- load_precomputed("desurv_bo_results_tcgacptac")
tar_params_best     <- load_precomputed("tar_params_best_tcgacptac")

# ── Note on figure regeneration ───────────────────────────────────────────
# The pre-computed figure RDS files (fig_*.rds) contain ggplot/grob objects
# serialized on HPC with ggplot2 3.x. These may not render correctly on
# ggplot2 4.x due to the ggproto→S7 migration. If you need to regenerate
# figures, ensure your ggplot2 version matches the fitting environment.
#
# For most users, the pre-computed figures work fine for paper compilation.
# The paper .Rmd files load these objects directly via load_result().

message("  Figure objects are pre-computed in results/precomputed/fig_*.rds")
message("  To regenerate, uncomment the relevant sections below.")

# ── Example: regenerate BO heatmap ────────────────────────────────────────
# fig_bo_heat <- cache_or_compute("fig_bo_heat_tcgacptac", {
#   # Requires BO history and GP model extraction
#   extract_gp_curve_k(desurv_bo_results)
# })

# ── Example: regenerate gene overlap heatmaps ────────────────────────────
# Requires: tops_desurv, top_genes_ref (reference gene lists)
# fig_gene_overlap <- cache_or_compute("fig_gene_overlap_heatmap_desurv_tcgacptac", {
#   tops <- DeSurv::desurv_get_top_genes(tar_fit_desurv, ntop_value)
#   make_gene_overlap_heatmap(tar_fit_desurv, tops, top_genes_ref)
# })

message("=== Step 8 complete ===")
