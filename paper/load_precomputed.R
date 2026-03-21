# paper/load_precomputed.R
# Drop-in replacement for tar_load/tar_read that reads from precomputed RDS files.
# Source this file at the top of each .Rmd instead of library(targets).

# When knit_root_dir = getwd() (repo root), paths are relative to root.
# When running from paper/ directly, paths need ../
PRECOMPUTED_DIR <- if (file.exists("results/precomputed")) {
  "results/precomputed"
} else if (file.exists("../results/precomputed")) {
  "../results/precomputed"
} else {
  stop("Cannot find results/precomputed/ directory. ",
       "Run from the repo root with knit_root_dir = getwd().")
}

load_result <- function(name, envir = parent.frame()) {
  path <- file.path(PRECOMPUTED_DIR, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Pre-computed result not found: ", path,
         "\nRun the pipeline or download pre-computed results first.")
  }
  obj <- readRDS(path)
  assign(name, obj, envir = envir)
  invisible(obj)
}

read_result <- function(name) {
  path <- file.path(PRECOMPUTED_DIR, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Pre-computed result not found: ", path)
  }
  readRDS(path)
}
