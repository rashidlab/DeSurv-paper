#!/bin/bash
# DeSurv robust-rerun environment setup on Longleaf.
# Scope: ONLY writes under /work/users/n/u/nur2/desurv-rerun/. Never touches /proj.
# Installs DeSurv (feat/consensus-estimator @ b42d7a9) into a private library and
# runs the install/consensus-fit preflight checks.
#
# Usage: bash setup_env.sh   (run on a Longleaf login node, or via
#   srun -p interact -n 1 -c 4 -t 120 --mem=8g --pty bash setup_env.sh
#   if compilation is refused on the login node)
set -euo pipefail

ROOT=/work/users/n/u/nur2/desurv-rerun
SRC="$ROOT/src/DeSurv"
LIB="$ROOT/lib/r4.5"
LOGDIR="$ROOT/logs"
mkdir -p "$LIB" "$LOGDIR"

module load r/4.5.0
export R_LIBS="$LIB"

echo "== R version ==" | tee "$LOGDIR/00_env.log"
R --version | tee -a "$LOGDIR/00_env.log"
hostname | tee -a "$LOGDIR/00_env.log"

# --- Non-CRAN dependencies ---
# cvwrapr: pulled from CRAN archive (off current CRAN); falls back to GitHub (Yixuan Qiu / trevorhastie) if archive install fails.
# preprocessCore: Bioconductor, installed via BiocManager.
R_LIBS="$LIB" Rscript - <<'RSCRIPT' 2>&1 | tee "$LOGDIR/01_nonCRAN_deps.log"
lib <- Sys.getenv("R_LIBS")
repos <- "https://cloud.r-project.org"

# BiocManager for preprocessCore
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib, repos = repos)
}
if (!requireNamespace("preprocessCore", quietly = TRUE, lib.loc = lib)) {
  BiocManager::install("preprocessCore", lib = lib, update = FALSE, ask = FALSE,
                        site_repository = repos)
}

# cvwrapr: try CRAN archive tarball first, GitHub as fallback
if (!requireNamespace("cvwrapr", quietly = TRUE, lib.loc = lib)) {
  ok <- tryCatch({
    install.packages("https://cran.r-project.org/src/contrib/Archive/cvwrapr/cvwrapr_1.0.tar.gz",
                      repos = NULL, type = "source", lib = lib)
    requireNamespace("cvwrapr", quietly = TRUE, lib.loc = lib)
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    if (!requireNamespace("remotes", quietly = TRUE, lib.loc = lib)) {
      install.packages("remotes", lib = lib, repos = repos)
    }
    library(remotes, lib.loc = lib)
    remotes::install_github("kjytay/cvwrapr", lib = lib)
  }
}
RSCRIPT

# --- Remaining CRAN dependencies (Imports + LinkingTo + Suggests needed for tests) ---
R_LIBS="$LIB" Rscript - <<RSCRIPT 2>&1 | tee "$LOGDIR/02_cran_deps.log"
lib <- "$LIB"
repos <- "https://cloud.r-project.org"
pkgs <- c("Rcpp", "survival", "cluster",
          "RcppArmadillo", "RcppEigen",
          "testthat", "knitr", "rmarkdown", "ncvreg", "DiceKriging", "lhs", "ggplot2")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE, lib.loc = lib)]
if (length(need)) install.packages(need, lib = lib, repos = repos)
RSCRIPT

# --- Install DeSurv itself ---
echo "== R CMD INSTALL DeSurv ==" | tee "$LOGDIR/03_install_desurv.log"
R_LIBS="$LIB" R CMD INSTALL --preclean --library="$LIB" "$SRC" 2>&1 | tee -a "$LOGDIR/03_install_desurv.log"

# --- Preflight checks ---
R_LIBS="$LIB" Rscript - <<RSCRIPT 2>&1 | tee "$LOGDIR/04_preflight.log"
lib <- "$LIB"
.libPaths(c(lib, .libPaths()))
library(DeSurv, lib.loc = lib)

cat("packageVersion(DeSurv):", as.character(packageVersion("DeSurv")), "\n")
stopifnot(packageVersion("DeSurv") == "1.1.0")

cat("cox_stats_cpp exists:", exists("cox_stats_cpp", asNamespace("DeSurv")), "\n")
stopifnot(exists("cox_stats_cpp", asNamespace("DeSurv")))

cat("desurv_consensus_fit exists:", exists("desurv_consensus_fit", asNamespace("DeSurv")), "\n")
stopifnot(exists("desurv_consensus_fit", asNamespace("DeSurv")))

ns_exports <- getNamespaceExports("DeSurv")
cat("desurv_consensus_gate exported:", "desurv_consensus_gate" %in% ns_exports, "\n")
stopifnot("desurv_consensus_gate" %in% ns_exports)

cat("== sessionInfo() ==\n")
print(sessionInfo())
RSCRIPT

# --- testthat on installed package ---
R_LIBS="$LIB" Rscript - <<RSCRIPT 2>&1 | tee "$LOGDIR/05_testthat.log"
lib <- "$LIB"
.libPaths(c(lib, .libPaths()))
library(DeSurv, lib.loc = lib)
library(testthat)

for (f in c("$SRC/tests/testthat/test-cox-ties.R",
            "$SRC/tests/testthat/test-consensus-fit.R")) {
  cat("== Running", f, "==\n")
  res <- test_file(f, reporter = "summary")
}
RSCRIPT

# --- HPC facts for later Slurm jobs ---
echo "== sinfo (general partition) ==" | tee "$LOGDIR/06_hpc_facts.log"
sinfo -p general -o "%P %C %D %m %c" | tee -a "$LOGDIR/06_hpc_facts.log"
hostname | tee -a "$LOGDIR/06_hpc_facts.log"

echo "Setup complete. Logs in $LOGDIR"
