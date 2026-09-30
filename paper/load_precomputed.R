# paper/load_precomputed.R
# Reads pre-computed RDS files from results/ for use in .Rmd files.
# Source this file at the top of each .Rmd instead of library(targets).
#
# When knit_root_dir = getwd() (repo root), "results/" resolves from root.
# When running from paper/ directly, falls back to "../results/".

# DESURV_RESULTS_DIR (see R/paths.R) overrides; use an absolute path when set,
# since knitting may run from paper/ or from the repo root.
.desurv_res_env <- Sys.getenv("DESURV_RESULTS_DIR", "")
PRECOMPUTED_DIR <- if (nzchar(.desurv_res_env)) .desurv_res_env else
  if (file.exists("results")) "results" else "../results"
CV_GRID_CSV_DIR <- file.path(PRECOMPUTED_DIR, "cv_grid")

# Figure files the Rmds embed. `default` is the historical relative path at the
# call site (either "../figures" or "figures", which differ by call site);
# DESURV_FIGURE_DIR, when set, replaces it.
.desurv_fig_env <- Sys.getenv("DESURV_FIGURE_DIR", "")
paper_fig <- function(name, default = "../figures") {
  file.path(if (nzchar(.desurv_fig_env)) .desurv_fig_env else default, name)
}

# Shared PACA-AU de-duplication helper (single source of truth for combined
# validation analyses; keeps one record per unique patient, preferring RNA-seq).
source(if (file.exists("R/paca_dedup.R")) "R/paca_dedup.R" else "../R/paca_dedup.R")

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
