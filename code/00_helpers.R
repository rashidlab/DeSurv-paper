#!/usr/bin/env Rscript
# code/00_helpers.R — Shared utilities for all pipeline scripts
#
# Provides:
#   CONFIG       — list with quick/ncores/recompute flags
#   RESULTS_DIR  — path to results/precomputed/
#   precomputed_path(name) — returns path to a precomputed RDS file
#   cache_or_compute(name, expr) — load from cache or compute and save

# ── Parse environment variables ────────────────────────────────────────────
CONFIG <- list(
  quick     = identical(Sys.getenv("DESURV_QUICK"), "TRUE"),
  ncores    = as.integer(Sys.getenv("DESURV_NCORES", "1")),
  recompute = identical(Sys.getenv("DESURV_RECOMPUTE"), "TRUE")
)

if (CONFIG$quick) {
  message("=== QUICK MODE: reduced data/iterations for smoke testing ===")
}

# ── ntop configuration ────────────────────────────────────────────────────
# DESURV_NTOP=150         → fixed ntop, subfolder ntop_150/
# DESURV_NTOP_LOWER=50 DESURV_NTOP_UPPER=250 → BO-tuned ntop, subfolder ntop_bo_50_250/
# (none)                  → ntop=NULL (all genes), original flat layout
ntop_fixed_env <- Sys.getenv("DESURV_NTOP", "")
ntop_lower_env <- Sys.getenv("DESURV_NTOP_LOWER", "")
ntop_upper_env <- Sys.getenv("DESURV_NTOP_UPPER", "")

if (nzchar(ntop_fixed_env) && nzchar(ntop_lower_env)) {
  warning("Both DESURV_NTOP and DESURV_NTOP_LOWER set. ",
          "Using fixed DESURV_NTOP=", ntop_fixed_env, " (BO bounds ignored).")
}

if (nzchar(ntop_fixed_env)) {
  CONFIG$ntop_mode      <- "fixed"
  CONFIG$ntop_value     <- as.integer(ntop_fixed_env)
  CONFIG$ntop_subfolder <- paste0("ntop_", CONFIG$ntop_value)
} else if (nzchar(ntop_lower_env) && nzchar(ntop_upper_env)) {
  CONFIG$ntop_mode      <- "bo"
  CONFIG$ntop_lower     <- as.integer(ntop_lower_env)
  CONFIG$ntop_upper     <- as.integer(ntop_upper_env)
  CONFIG$ntop_subfolder <- paste0("ntop_bo_", CONFIG$ntop_lower, "_", CONFIG$ntop_upper)
} else {
  CONFIG$ntop_mode      <- "default"
  CONFIG$ntop_value     <- NULL
  CONFIG$ntop_subfolder <- ""
}

# DESURV_PARAM_RULE: "best" (default) or "1se".
# When "1se", the final-model parameter selection in step 03b takes the full
# 1-SE row from BO history (not just k), and the subfolder gets a _1se suffix
# so outputs do not overwrite the "best" run.
param_rule_env <- Sys.getenv("DESURV_PARAM_RULE", "best")
if (!nzchar(param_rule_env)) param_rule_env <- "best"
if (!param_rule_env %in% c("best", "1se")) {
  stop("DESURV_PARAM_RULE must be 'best' or '1se' (got: ", param_rule_env, ")")
}
CONFIG$param_rule <- param_rule_env

if (CONFIG$param_rule == "1se") {
  if (!nzchar(CONFIG$ntop_subfolder)) {
    stop("DESURV_PARAM_RULE=1se requires DESURV_NTOP or DESURV_NTOP_LOWER/UPPER ",
         "to be set so outputs land in a dedicated subfolder.")
  }
  CONFIG$ntop_subfolder <- paste0(CONFIG$ntop_subfolder, "_1se")
}

if (nzchar(CONFIG$ntop_subfolder)) {
  message("=== ntop config: ", CONFIG$ntop_mode,
          ", rule: ", CONFIG$param_rule,
          " -> subfolder: ", CONFIG$ntop_subfolder, " ===")
}

# ── Paths ──────────────────────────────────────────────────────────────────
RESULTS_BASE <- "results/precomputed"
RESULTS_DIR  <- if (nzchar(CONFIG$ntop_subfolder)) {
  file.path(RESULTS_BASE, CONFIG$ntop_subfolder)
} else {
  RESULTS_BASE
}
CV_GRID_DIR <- "results/cv_grid"
FIGURE_DIR  <- if (nzchar(CONFIG$ntop_subfolder)) {
  file.path("figures", CONFIG$ntop_subfolder)
} else {
  "figures"
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CV_GRID_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR,  recursive = TRUE, showWarnings = FALSE)

precomputed_path <- function(name) {
  file.path(RESULTS_DIR, paste0(name, ".rds"))
}

# ── Cache-or-compute pattern ──────────────────────────────────────────────
# If the file exists and we're not forcing recompute, load it.
# Otherwise evaluate expr, save, and return.
cache_or_compute <- function(name, expr) {
  path <- precomputed_path(name)
  if (!CONFIG$recompute && file.exists(path)) {
    message("  Loading cached: ", name)
    return(readRDS(path))
  }
  message("  Computing: ", name, " ...")
  t0 <- Sys.time()
  result <- force(expr)
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  saveRDS(result, path)
  message(sprintf("  Saved: %s (%.1f min)", path, as.numeric(elapsed)))
  result
}

# ── Load a previously computed result (error if missing) ──────────────────
# Falls back to base results directory for shared outputs (e.g., step 02 data)
load_precomputed <- function(name) {
  path <- precomputed_path(name)
  if (!file.exists(path) && nzchar(CONFIG$ntop_subfolder)) {
    base_path <- file.path(RESULTS_BASE, paste0(name, ".rds"))
    if (file.exists(base_path)) {
      # Only allow fallback for shared step-2 outputs (training data).
      # Everything else should exist in the subfolder after step 3 runs.
      shared_base_files <- c("tar_data_tcgacptac", "sim_figs_by_scenario")
      if (!name %in% shared_base_files) {
        stop("Required result not found in subfolder: ", path,
             "\n  Found in base directory but refusing silent fallback.",
             "\n  Run step 3 with correct DESURV_NTOP env vars first.")
      }
      message("  Loading from base directory: ", name)
      return(readRDS(base_path))
    }
  }
  if (!file.exists(path)) {
    stop("Required result not found: ", path,
         "\nRun earlier pipeline steps first, or use pre-computed results.")
  }
  readRDS(path)
}
