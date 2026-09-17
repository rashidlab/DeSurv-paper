#!/usr/bin/env bash
# run_fair_comparator.sh
# ---------------------------------------------------------------------------
# Builds the fair matched-rank comparator (alpha = 0 at DeSurv's k, with its own
# BO tuning budget) and then runs the pre-registered decision gate.
#
# SAFETY. code/00_helpers.R:18 defines
#     recompute = !identical(Sys.getenv("DESURV_RECOMPUTE"), "FALSE")
# so the DEFAULT IS TO RECOMPUTE AND OVERWRITE every cached object, including
# desurv_bo_results_tcgacptac and tar_params_best_tcgacptac, which anchor every
# number in the manuscript. DESURV_RECOMPUTE=FALSE below is load-bearing, and
# the preflight aborts if it is not in effect.
#
# RESUMABLE. With DESURV_RECOMPUTE=FALSE, cache_or_compute() skips any object
# already on disk. If this dies partway, just re-run it: completed stages load
# instantly and it continues from the first missing object.
#
# Usage:  nohup setsid bash slurm/run_fair_comparator.sh > LOG 2>&1 &
# ---------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

export DESURV_RECOMPUTE=FALSE
export DESURV_NCORES="${DESURV_NCORES:-14}"

say() { echo "###### $* | $(date '+%F %T') ######"; }

say "PREFLIGHT"
Rscript -e '
  suppressMessages(source("code/00_helpers.R"))
  stopifnot("DESURV_RECOMPUTE=FALSE is NOT in effect" = isFALSE(CONFIG$recompute))
  stopifnot("RESULTS_DIR must be results/ (not quick)" = RESULTS_DIR == "results")
  cat("  recompute =", CONFIG$recompute, "| ncores =", CONFIG$ncores,
      "| dir =", RESULTS_DIR, "\n")
  # The canonical objects must LOAD, never recompute. stop() inside the expr
  # fires only if the cache is being bypassed.
  for (n in c("desurv_bo_results_tcgacptac", "tar_params_best_tcgacptac",
              "tar_fit_desurv_tcgacptac", "tar_data_filtered_tcgacptac")) {
    invisible(cache_or_compute(n, stop("CACHE BYPASSED for ", n)))
  }
  cat("  DeSurv", as.character(packageVersion("DeSurv")), "\n")
  stopifnot("DeSurv must be 1.0.1 (pinned submission/pnas-2026)" =
            as.character(packageVersion("DeSurv")) == "1.0.1")
  cat("  PREFLIGHT OK: canonical cache loads, nothing will be overwritten\n")
' || { say "PREFLIGHT FAILED - ABORTING"; exit 1; }

# Background memory heartbeat so a slow death is diagnosable after the fact.
( while true; do
    echo "  [mem $(date '+%T')] $(free -g | awk '/^Mem:/{print "used="$3"G avail="$7"G"}') load=$(cut -d' ' -f1 /proc/loadavg)"
    sleep 300
  done ) &
HEARTBEAT=$!
trap 'kill $HEARTBEAT 2>/dev/null' EXIT

for s in 03_bayesian_optimization 04_fit_models 05_external_validation; do
  say "START $s"
  Rscript "code/${s}.R" || { say "FAILED $s"; exit 1; }
  say "DONE $s"
done

say "START decision gate"
Rscript code/22_fair_comparator_gate.R || { say "FAILED gate"; exit 1; }

say "FAIR COMPARATOR RUN COMPLETE"
