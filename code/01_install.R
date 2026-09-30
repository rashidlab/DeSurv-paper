#!/usr/bin/env Rscript
# code/01_install.R — Install the DeSurv package
#
# Installs from the source named in DESURV_PKG_SOURCE (local dir or GitHub spec); errors if unset.
# Also verifies required dependencies are available.

message("=== Step 1: Install DeSurv ===")

# Rev 3 (2026-09-30): the package source must be given EXPLICITLY. There is no
# default and no silent fallback. In particular the sibling checkout ../DeSurv
# (a fork, not the rev 3 package) is never installed.
#
#   DESURV_PKG_SOURCE  either a local package directory, e.g.
#                        DESURV_PKG_SOURCE=~/Downloads/DeSurv-clean   (branch
#                        feat/reported-estimator, 1f3d3f0, v1.2.0; unpushed)
#                      or a GitHub spec "owner/repo@ref", e.g.
#                        DESURV_PKG_SOURCE=rashidlab/DeSurv@<ref>
#   DESURV_MIN_VERSION minimum acceptable installed version (default "1.2.0")
#
# The pre-rev-3 submission build was rashidlab/DeSurv@submission/pnas-2026
# (afb00d5, v1.0.1); install it only with DESURV_PKG_SOURCE plus
# DESURV_MIN_VERSION=1.0.1. Slurm jobs do not run this script; they use the
# Longleaf library /work/users/n/u/nur2/desurv-rerun/lib/r4.5-rev3.
pkg_source  <- Sys.getenv("DESURV_PKG_SOURCE", "")
min_version <- Sys.getenv("DESURV_MIN_VERSION", "1.2.0")
if (!nzchar(pkg_source)) {
  stop("DESURV_PKG_SOURCE is not set. Set it to a local DeSurv package directory ",
       "or to a GitHub spec 'owner/repo@ref'. Refusing to guess (no ../DeSurv fallback).",
       call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = "https://cloud.r-project.org")

pkg_dir <- path.expand(pkg_source)
if (dir.exists(pkg_dir)) {
  sibling <- file.path("..", "DeSurv")
  if (dir.exists(sibling) &&
      identical(normalizePath(pkg_dir), normalizePath(sibling))) {
    stop("DESURV_PKG_SOURCE points at ../DeSurv, which is not the rev 3 package. ",
         "Point it at the rev 3 checkout explicitly.", call. = FALSE)
  }
  git_head <- tryCatch(system2("git", c("-C", shQuote(pkg_dir), "rev-parse", "--short", "HEAD"),
                               stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
  message("Installing DeSurv from local directory: ", normalizePath(pkg_dir),
          " (git HEAD ", paste(git_head, collapse = ""), ")")
  devtools::install_local(pkg_dir, upgrade = "never", force = TRUE)
} else if (grepl("^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+@.+$", pkg_source)) {
  message("Installing DeSurv from GitHub: ", pkg_source)
  devtools::install_github(pkg_source, upgrade = "never", force = TRUE)
} else {
  stop("DESURV_PKG_SOURCE='", pkg_source, "' is neither an existing directory nor ",
       "a GitHub spec 'owner/repo@ref'.", call. = FALSE)
}

# Verify DeSurv loads
library(DeSurv)
message("DeSurv version: ", packageVersion("DeSurv"))
if (packageVersion("DeSurv") < min_version) {
  stop("Installed DeSurv ", packageVersion("DeSurv"), " < required ", min_version,
       " (DESURV_MIN_VERSION).", call. = FALSE)
}

# NOTE: Cached figure .rds files in results/ are saved using the
# current ggplot2 version's class system (ggplot2 >= 4.0 uses S7 with @-slot
# accessors). If you upgrade or downgrade ggplot2 across a major version, you
# may need to regenerate the cached figures by re-running step 8.
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")

# Check other key dependencies
pkgs <- c("survival", "NMF", "ggplot2", "cowplot", "dplyr", "survminer",
           "pheatmap", "ComplexHeatmap", "DiceKriging", "lhs",
           "glmnet", "riskRegression", "prodlim", "superpc", "ggtext")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

message("=== Step 1 complete ===")
