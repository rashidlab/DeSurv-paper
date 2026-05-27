#!/usr/bin/env Rscript
# code/01_install.R — Install the DeSurv package
#
# Installs from GitHub (rashidlab/DeSurv) or from a local checkout if available.
# Also verifies required dependencies are available.

message("=== Step 1: Install DeSurv ===")

# Pinned to the rashidlab/DeSurv tag that produced the cached .rds files in
# results/precomputed/. Reviewers reproducing the manuscript should install
# this exact version; the rashidlab/DeSurv tip may include post-submission
# fixes (e.g., RNG seeding) that change numerical results.
DESURV_REF <- "submission/pnas-2026"   # rashidlab/DeSurv @ afb00d5

# Try local install first (for development), then GitHub
local_path <- file.path("..", "DeSurv")
if (dir.exists(local_path)) {
  message("Installing DeSurv from local checkout: ", local_path)
  message("  NOTE: local checkout is NOT pinned — confirm it is at tag '", DESURV_REF, "'")
  message("        (rashidlab/DeSurv @ afb00d5) before relying on byte-exact reproduction.")
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = "https://cloud.r-project.org")
  devtools::install_local(local_path, upgrade = "never", force = FALSE)
} else {
  message("Installing DeSurv from GitHub: rashidlab/DeSurv@", DESURV_REF)
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = "https://cloud.r-project.org")
  devtools::install_github(paste0("rashidlab/DeSurv@", DESURV_REF), upgrade = "never", force = FALSE)
}

# Verify DeSurv loads
library(DeSurv)
message("DeSurv version: ", packageVersion("DeSurv"))

# NOTE: Cached figure .rds files in results/precomputed/ are saved using the
# current ggplot2 version's class system (ggplot2 >= 4.0 uses S7 with @-slot
# accessors). If you upgrade or downgrade ggplot2 across a major version, you
# may need to regenerate the cached figures by re-running step 8.
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")

# Check other key dependencies
pkgs <- c("survival", "NMF", "ggplot2", "cowplot", "dplyr", "survminer",
           "pheatmap", "ComplexHeatmap", "DiceKriging", "lhs")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

message("=== Step 1 complete ===")
