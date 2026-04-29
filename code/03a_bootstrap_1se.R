#!/usr/bin/env Rscript
# code/03a_bootstrap_1se.R — Seed the _1se subfolder with BO inputs
#
# Copies the cached Bayesian optimization outputs from the parent
# ntop-subfolder (e.g., ntop_bo_50_300/) into the _1se subfolder
# (ntop_bo_50_300_1se/) so that step 03b can read them and the BO
# itself does not need to be re-run.
#
# This script is a no-op when DESURV_PARAM_RULE != "1se".

message("=== Step 3a: Bootstrap _1se subfolder ===")
source("code/00_helpers.R")

if (CONFIG$param_rule != "1se") {
  message("  DESURV_PARAM_RULE != '1se'; nothing to do.")
} else {
  parent_subfolder <- sub("_1se$", "", CONFIG$ntop_subfolder)
  if (parent_subfolder == CONFIG$ntop_subfolder) {
    stop("Expected ntop_subfolder to end with '_1se' but got '",
         CONFIG$ntop_subfolder, "'")
  }
  parent_dir <- file.path(RESULTS_BASE, parent_subfolder)
  if (!dir.exists(parent_dir)) {
    stop("Parent results directory not found: ", parent_dir,
         "\nExpected BO outputs from the 'best' run to live there.")
  }

  bo_inputs <- c(
    "desurv_bo_results_tcgacptac.rds",
    "desurv_bo_results_alpha0_tcgacptac.rds",
    "desurv_bo_results_elbowk_tcgacptac.rds",
    "fit_std_tcgacptac.rds",
    "std_nmf_selected_k_tcgacptac.rds"
  )

  copied <- 0L
  skipped <- 0L
  for (fname in bo_inputs) {
    src <- file.path(parent_dir, fname)
    dst <- file.path(RESULTS_DIR, fname)
    if (!file.exists(src)) {
      stop("Required BO input not found in parent: ", src)
    }
    if (file.exists(dst)) {
      message("  Already present, skipping: ", fname)
      skipped <- skipped + 1L
      next
    }
    ok <- file.copy(src, dst, overwrite = FALSE, copy.date = TRUE)
    if (!ok) stop("file.copy failed for ", fname)
    message("  Copied: ", fname)
    copied <- copied + 1L
  }
  message(sprintf("  Bootstrap complete: %d copied, %d already present.",
                  copied, skipped))
}

message("=== Step 3a complete ===")
