#!/usr/bin/env Rscript
# export_precomputed.R — Extract paper-needed objects from targets store
#
# Run ONCE from the DeSurv-paper repo to export all pre-computed
# results needed for the clean repo.
#
# Usage: Rscript export_precomputed.R

library(targets)
tar_config_set(store = "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main")

out_dir <- "results/precomputed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Targets objects needed by the paper ─────────────────────────────────────
paper_targets <- c(
  # Main text figures (04_results_REVISED.Rmd)
  "sim_figs_by_scenario",
  "fig_bo_heat_tcgacptac",
  "fig_gene_overlap_heatmap_desurv_tcgacptac",
  "fig_gene_overlap_heatmap_std_desurvk_tcgacptac",
  "fig_variation_explained_tcgacptac",
  "fig_desurv_std_correlation_tcgacptac",
  "fig_hr_forest_tcgacptac",
  "fig_median_survival_desurv_tcgacptac",
  "fig_median_survival_std_desurvk_tcgacptac",
  # Main text data (for inline statistics)
  "val_latent_desurv_tcgacptac",
  "val_latent_std_desurvk_tcgacptac",
  "tar_k_selection_tcgacptac",
  "desurv_bo_results_tcgacptac",
  "tar_params_best_tcgacptac",
  "tar_data_filtered_tcgacptac",
  # Supplement figures
  "desurv_seed_fits_tcgacptac",
  "fig_residuals_tcgacptac",
  "fig_cophenetic_tcgacptac",
  "fig_silhouette_tcgacptac",
  "fit_std_tcgacptac",
  "fig_gene_overlap_heatmap_std_elbowk_tcgacptac",
  # Supplement tables
  "val_cindex_desurv_tcgacptac",
  "val_cindex_std_desurvk_tcgacptac",
  "val_cindex_std_elbowk_tcgacptac",
  "val_cindex_desurv_alpha0_tcgacptac",
  # Supplement data (adjusted HR analysis)
  "val_latent_desurv_alpha0_tcgacptac",
  "data_val_filtered_tcgacptac",
  "fit_std_elbowk_tcgacptac",
  "tar_data_filtered_elbowk_tcgacptac"
)

cat("Exporting", length(paper_targets), "targets to", out_dir, "\n\n")

failures <- character(0)
for (tgt in paper_targets) {
  out_path <- file.path(out_dir, paste0(tgt, ".rds"))
  if (file.exists(out_path)) {
    cat("  [skip] ", tgt, " (already exists)\n")
    next
  }
  cat("  [export] ", tgt, " ... ")
  result <- tryCatch({
    obj <- tar_read_raw(tgt)
    saveRDS(obj, out_path)
    sz <- file.size(out_path)
    cat(sprintf("%.1f MB\n", sz / 1e6))
    TRUE
  }, error = function(e) {
    cat("FAILED:", conditionMessage(e), "\n")
    failures <<- c(failures, tgt)
    FALSE
  })
}

# Special case: object with no metadata (must read from store directly)
special <- "fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac"
special_path <- file.path(out_dir, paste0(special, ".rds"))
if (!file.exists(special_path)) {
  cat("  [export] ", special, " (direct from store) ... ")
  tryCatch({
    obj <- readRDS(file.path(tar_config_get("store"), "objects", special))
    saveRDS(obj, special_path)
    sz <- file.size(special_path)
    cat(sprintf("%.1f MB\n", sz / 1e6))
  }, error = function(e) {
    cat("FAILED:", conditionMessage(e), "\n")
    failures <- c(failures, special)
  })
} else {
  cat("  [skip] ", special, " (already exists)\n")
}

cat("\n")
if (length(failures)) {
  cat("FAILURES (", length(failures), "):\n")
  cat(paste("  -", failures), sep = "\n")
} else {
  n_exported <- length(paper_targets) + 1
  total_size <- sum(file.size(list.files(out_dir, full.names = TRUE, pattern = "\\.rds$")))
  cat(sprintf("Success: %d objects exported (%.0f MB total)\n", n_exported, total_size / 1e6))
}
