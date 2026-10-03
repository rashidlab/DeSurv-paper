#!/bin/bash
# slurm/rev3_env.sh: shared preamble for every rev 3 Slurm runner (source it, do not run it).
#
# What it guarantees before any R work starts:
#   1. R 4.5.0 (module r/4.5.0) with the rev 3 private library FIRST on .libPaths().
#   2. DeSurv >= 1.2.0, loaded from that library (fails otherwise).
#   3. Output isolation: DESURV_RESULTS_DIR and DESURV_FIGURE_DIR must be set, are made
#      absolute, and must NOT be the tracked results/ or figures/ of this repo.
#   4. DESURV_RECOMPUTE=TRUE explicitly (the runners write into an isolated, initially
#      empty directory; see code/00_helpers.R for why the default recomputes).
#   5. slurm/logs exists (note: #SBATCH --output is opened before the script runs, so
#      create slurm/logs before the FIRST sbatch as well; slurm/submit_*.sh do).
#
# Usage inside a runner (after the #SBATCH header):
#   source "$SLURM_SUBMIT_DIR/slurm/rev3_env.sh"
# Required exports at submit time (sbatch passes the environment by default):
#   DESURV_RESULTS_DIR=/work/users/n/u/nur2/desurv-rerun/out/<run>/results
#   DESURV_FIGURE_DIR=/work/users/n/u/nur2/desurv-rerun/out/<run>/figures
#   DESURV_CMBSUBTYPES=<path to cmbSubtypes_formatted.RData>   (code/09a, 23, 24 only)

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
[ -f code/00_helpers.R ] || { echo "rev3_env: $REPO_ROOT is not the paper repo root" >&2; exit 1; }
mkdir -p "$REPO_ROOT/slurm/logs"

module purge
module load r/4.5.0

export REV3_LIB="${REV3_LIB:-/work/users/n/u/nur2/desurv-rerun/lib/r4.5-rev3}"
[ -d "$REV3_LIB/DeSurv" ] || { echo "rev3_env: DeSurv not installed in $REV3_LIB" >&2; exit 1; }
# R_LIBS is searched before R_LIBS_USER; set both so ~/.Renviron cannot put another lib first.
export R_LIBS="$REV3_LIB"
export R_LIBS_USER="$REV3_LIB"

: "${DESURV_RESULTS_DIR:?rev3_env: set DESURV_RESULTS_DIR to an isolated output directory}"
: "${DESURV_FIGURE_DIR:?rev3_env: set DESURV_FIGURE_DIR to an isolated figure directory}"
DESURV_RESULTS_DIR="$(realpath -m "$DESURV_RESULTS_DIR")"
DESURV_FIGURE_DIR="$(realpath -m "$DESURV_FIGURE_DIR")"
for d in "$DESURV_RESULTS_DIR" "$DESURV_FIGURE_DIR"; do
  case "$d" in
    "$(realpath -m "$REPO_ROOT/results")"*|"$(realpath -m "$REPO_ROOT/figures")"*)
      echo "rev3_env: $d is inside the tracked results/ or figures/; refusing" >&2; exit 1;;
  esac
done
mkdir -p "$DESURV_RESULTS_DIR" "$DESURV_FIGURE_DIR"
export DESURV_RESULTS_DIR DESURV_FIGURE_DIR

export DESURV_RECOMPUTE=TRUE
export DESURV_QUICK="${DESURV_QUICK:-FALSE}"
export DESURV_NCORES="${SLURM_CPUS_PER_TASK:-1}"

# DeSurv version and library check, before any work.
Rscript -e '
  lib <- normalizePath(Sys.getenv("REV3_LIB"))
  v   <- packageVersion("DeSurv")
  p   <- normalizePath(find.package("DeSurv"))
  cat("rev3_env: R", paste(R.version$major, R.version$minor, sep = "."),
      "| DeSurv", as.character(v), "from", p, "\n")
  cat("rev3_env: .libPaths()[1] =", .libPaths()[1], "\n")
  stopifnot("DeSurv >= 1.2.0 required (rev 3)" = v >= "1.2.0",
            "DeSurv must load from REV3_LIB" = startsWith(p, lib))
' || { echo "rev3_env: DeSurv version/library check FAILED" >&2; exit 1; }

echo "rev3_env: results -> $DESURV_RESULTS_DIR"
echo "rev3_env: figures -> $DESURV_FIGURE_DIR"
echo "rev3_env: DESURV_CMBSUBTYPES=${DESURV_CMBSUBTYPES:-<unset; code/09a, 23, 24 need it>}"
echo "rev3_env: DESURV_RECOMPUTE=$DESURV_RECOMPUTE DESURV_NCORES=$DESURV_NCORES DESURV_QUICK=$DESURV_QUICK"
echo "rev3_env: job ${SLURM_JOB_ID:-local} on ${SLURM_NODELIST:-$(hostname)} at $(date)"
