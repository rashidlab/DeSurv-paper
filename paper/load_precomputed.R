# paper/load_precomputed.R
# Drop-in replacement for tar_load/tar_read that reads from precomputed RDS files.
# Source this file at the top of each .Rmd instead of library(targets).
#
# Respects DESURV_NTOP / DESURV_NTOP_LOWER / DESURV_NTOP_UPPER env vars
# to load from ntop-specific subfolders. Falls back to base directory for
# shared results (e.g., step 02 data, simulations).

# ── Determine base directory ─────────────────────────────────────────────
# When knit_root_dir = getwd() (repo root), paths are relative to root.
# When running from paper/ directly, paths need ../
PRECOMPUTED_BASE <- if (file.exists("results/precomputed")) {
  "results/precomputed"
} else if (file.exists("../results/precomputed")) {
  "../results/precomputed"
} else {
  stop("Cannot find results/precomputed/ directory. ",
       "Run from the repo root with knit_root_dir = getwd().")
}

# ── ntop subfolder detection ─────────────────────────────────────────────
.ntop_fixed_env <- Sys.getenv("DESURV_NTOP", "")
.ntop_lower_env <- Sys.getenv("DESURV_NTOP_LOWER", "")
.ntop_upper_env <- Sys.getenv("DESURV_NTOP_UPPER", "")

.ntop_subfolder <- if (nzchar(.ntop_fixed_env)) {
  paste0("ntop_", .ntop_fixed_env)
} else if (nzchar(.ntop_lower_env) && nzchar(.ntop_upper_env)) {
  paste0("ntop_bo_", .ntop_lower_env, "_", .ntop_upper_env)
} else {
  ""
}

PRECOMPUTED_DIR <- if (nzchar(.ntop_subfolder)) {
  subdir <- file.path(PRECOMPUTED_BASE, .ntop_subfolder)
  if (!dir.exists(subdir)) {
    warning("ntop subfolder not found: ", subdir, ". Falling back to base directory.")
    PRECOMPUTED_BASE
  } else {
    subdir
  }
} else {
  PRECOMPUTED_BASE
}

# ── Helper: resolve path with fallback to base dir ───────────────────────
.resolve_result_path <- function(name) {
  path <- file.path(PRECOMPUTED_DIR, paste0(name, ".rds"))
  if (!file.exists(path) && nzchar(.ntop_subfolder)) {
    base_path <- file.path(PRECOMPUTED_BASE, paste0(name, ".rds"))
    if (file.exists(base_path)) return(base_path)
  }
  path
}

load_result <- function(name, envir = parent.frame()) {
  path <- .resolve_result_path(name)
  if (!file.exists(path)) {
    stop("Pre-computed result not found: ", name,
         "\nRun the pipeline or download pre-computed results first.")
  }
  obj <- readRDS(path)
  assign(name, obj, envir = envir)
  invisible(obj)
}

read_result <- function(name) {
  path <- .resolve_result_path(name)
  if (!file.exists(path)) {
    stop("Pre-computed result not found: ", name)
  }
  readRDS(path)
}
