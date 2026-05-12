#!/usr/bin/env Rscript
# code/00_helpers.R — Shared utilities for all pipeline scripts
#
# Provides:
#   CONFIG       — list with quick/ncores/recompute flags and hard-coded ntop bounds
#   RESULTS_DIR  — "results" (flat, no subfolders)
#   CV_GRID_DIR  — "results/cv_grid"
#   FIGURE_DIR   — "figures"
#   precomputed_path(name) — returns path to a result RDS file
#   cache_or_compute(name, expr) — load from cache or compute and save
#   load_precomputed(name) — load a previously saved result (error if missing)

# ── Parse environment variables ────────────────────────────────────────────
CONFIG <- list(
  quick      = identical(Sys.getenv("DESURV_QUICK"), "TRUE"),
  ncores     = as.integer(Sys.getenv("DESURV_NCORES", "1")),
  # Default TRUE; set DESURV_RECOMPUTE=FALSE to load from cache instead
  recompute  = !identical(Sys.getenv("DESURV_RECOMPUTE"), "FALSE"),
  ntop_lower = 50L,
  ntop_upper = 300L
)

if (CONFIG$quick) {
  message("=== QUICK MODE: reduced data/iterations for smoke testing ===")
}

# ── Paths ──────────────────────────────────────────────────────────────────
RESULTS_DIR <- "results"
CV_GRID_DIR <- "results/cv_grid"
FIGURE_DIR  <- "figures"

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CV_GRID_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR,  recursive = TRUE, showWarnings = FALSE)

precomputed_path <- function(name) {
  file.path(RESULTS_DIR, paste0(name, ".rds"))
}

# ── Cache-or-compute pattern ──────────────────────────────────────────────
# Default behaviour is to recompute and save. Set DESURV_RECOMPUTE=FALSE to
# skip recomputation and load the cached file instead.
cache_or_compute <- function(name, expr) {
  path <- precomputed_path(name)
  if (!CONFIG$recompute && file.exists(path)) {
    message("  Loading cached: ", name)
    return(readRDS(path))
  }
  message("  Computing: ", name, " ...")
  t0     <- Sys.time()
  result <- force(expr)
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  saveRDS(result, path)
  message(sprintf("  Saved: %s (%.1f min)", path, as.numeric(elapsed)))
  result
}

# ── Load a previously computed result (error if missing) ──────────────────
load_precomputed <- function(name) {
  path <- precomputed_path(name)
  if (!file.exists(path)) {
    stop("Required result not found: ", path,
         "\nRun earlier pipeline steps first, or use pre-computed results.")
  }
  readRDS(path)
}
