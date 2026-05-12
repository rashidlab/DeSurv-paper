# paper/load_precomputed.R
# Reads pre-computed RDS files from results/ for use in .Rmd files.
# Source this file at the top of each .Rmd instead of library(targets).
#
# When knit_root_dir = getwd() (repo root), "results/" resolves from root.
# When running from paper/ directly, falls back to "../results/".

PRECOMPUTED_DIR <- if (file.exists("results")) "results" else "../results"

load_result <- function(name, envir = parent.frame()) {
  path <- file.path(PRECOMPUTED_DIR, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Pre-computed result not found: ", name,
         "\nRun the pipeline or download pre-computed results first.")
  }
  obj <- readRDS(path)
  assign(name, obj, envir = envir)
  invisible(obj)
}

read_result <- function(name) {
  path <- file.path(PRECOMPUTED_DIR, paste0(name, ".rds"))
  if (!file.exists(path)) stop("Pre-computed result not found: ", name)
  readRDS(path)
}
