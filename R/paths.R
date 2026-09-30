# R/paths.R: single source of truth for the results and figure base directories.
#
# Side-effect free (no dir.create, no messages) so standalone revision scripts can
# source it without pulling in code/00_helpers.R.
#
#   DESURV_RESULTS_DIR  base results directory (default "results", relative to repo root)
#   DESURV_FIGURE_DIR   base figure directory  (default "figures", relative to repo root)
#
# Rev 3 reruns set both to an EMPTY directory (absolute path recommended) so that
# tracked results/*.rds and figures/*.pdf are never overwritten. Unset, every path
# below is byte-identical to the historical hardcoded "results/..." / "figures/...".

RESULTS_BASE <- Sys.getenv("DESURV_RESULTS_DIR", "results")
FIGURE_BASE  <- Sys.getenv("DESURV_FIGURE_DIR",  "figures")
if (!nzchar(RESULTS_BASE)) RESULTS_BASE <- "results"
if (!nzchar(FIGURE_BASE))  FIGURE_BASE  <- "figures"

# res_path("x.rds") -> "results/x.rds" by default; fig_path("x.pdf") -> "figures/x.pdf"
res_path <- function(...) file.path(RESULTS_BASE, ...)
fig_path <- function(...) file.path(FIGURE_BASE, ...)

# Reference subtype gene lists (cmbSubtypes_formatted.RData: top_genes, colors, ...),
# read by code/09a, code/23 and code/24. Not tracked. One path, no fallback:
#   DESURV_CMBSUBTYPES  if set; otherwise data/derv/cmbSubtypes_formatted.RData.
# Errors if the file is missing (the old silent fallback to ../DeSurv-paper/ is removed).
cmb_subtypes_path <- function() {
  p <- Sys.getenv("DESURV_CMBSUBTYPES", "data/derv/cmbSubtypes_formatted.RData")
  if (!nzchar(p)) p <- "data/derv/cmbSubtypes_formatted.RData"
  if (!file.exists(p)) {
    stop("Reference gene list not found: ", p,
         "\nSet DESURV_CMBSUBTYPES to the cmbSubtypes_formatted.RData file ",
         "(md5 ee4c81468936b0a6b514fe4b4e41fe65) or place it at data/derv/.", call. = FALSE)
  }
  p
}
