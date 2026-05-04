#!/usr/bin/env Rscript
# code/01_install.R — Install the DeSurv package
#
# Installs from GitHub (rashidlab/DeSurv) or from a local checkout if available.
# Also verifies required dependencies are available.

message("=== Step 1: Install DeSurv ===")

# Try local install first (for development), then GitHub
local_path <- file.path("..", "DeSurv")
if (dir.exists(local_path)) {
  message("Installing DeSurv from local checkout: ", local_path)
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_local(local_path, upgrade = "never", force = FALSE)
} else {
  message("Installing DeSurv from GitHub: rashidlab/DeSurv")
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("rashidlab/DeSurv", upgrade = "never", force = FALSE)
}

# Verify DeSurv loads
library(DeSurv)
message("DeSurv version: ", packageVersion("DeSurv"))

# NOTE: Cached figure .rds files in results/precomputed/ are saved using the
# current ggplot2 version's class system (ggplot2 >= 4.0 uses S7 with @-slot
# accessors). If you upgrade or downgrade ggplot2 across a major version, you
# may need to regenerate the cached figures by re-running step 8.
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Check other key dependencies
pkgs <- c("survival", "NMF", "ggplot2", "cowplot", "dplyr", "survminer",
           "pheatmap", "ComplexHeatmap")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing)
}

message("=== Step 1 complete ===")
