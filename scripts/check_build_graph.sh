#!/usr/bin/env bash
# LESSONS.md L12: every cached object the manuscript reads must be produced by a
# script that `make all` runs. A green render proves nothing about clean-build
# reproducibility when a required cache lives outside the build graph.
# Exits 1 if any read_result()/load_result() name in paper/*.Rmd is written by
# a code/*.R script that `make -n all` does not invoke, or by no script at all.
set -u; cd "$(dirname "$0")/.." >/dev/null
fail=0
names=$(grep -ohE '(read_result|load_result)\("[A-Za-z0-9_]+"\)' paper/*.Rmd | grep -oE '"[^"]+"' | tr -d '"' | sort -u)
allrun=$(make -n all 2>/dev/null | grep -oE 'code/[0-9a-z]+_[A-Za-z0-9_]+\.R' | sort -u)
for n in $names; do
  # out-of-band artifacts (restricted external data) are exempt by design; see Makefile clean:
  case "$n" in treated_cohort_stats|spatial_cooccurrence_stats|spatial_spots_scored|spatial_adjacency_perspot) continue;; esac
  writers=$(grep -lE "saveRDS\([^,]*,[^)]*\b${n}\.rds|cache_or_compute\(\"${n}\"" code/*.R 2>/dev/null | sort -u)
  if [ -z "$writers" ]; then echo "NO WRITER  : $n (no code/*.R writes results/$n.rds)"; fail=1; continue; fi
  ok=0; for w in $writers; do echo "$allrun" | grep -qx "$w" && ok=1; done
  [ $ok = 1 ] || { echo "NOT IN all : $n  (written by: $writers)"; fail=1; }
done
[ $fail = 0 ] && echo "PASS: every manuscript cache is produced by a script in 'make all'"
exit $fail
