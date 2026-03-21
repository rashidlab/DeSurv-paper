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

# ── Paths ──────────────────────────────────────────────────────────────────
RESULTS_DIR <- "results/precomputed"
CV_GRID_DIR <- "results/cv_grid"
FIGURE_DIR  <- "figures"

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
load_precomputed <- function(name) {
  path <- precomputed_path(name)
  if (!file.exists(path)) {
    stop("Required result not found: ", path,
         "\nRun earlier pipeline steps first, or use pre-computed results.")
  }
  readRDS(path)
}
