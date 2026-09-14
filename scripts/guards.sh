#!/bin/bash
# render + full guard set; prints one line per check
# Usage: scripts/guards.sh   (from anywhere; renders with DESURV_RECOMPUTE=FALSE, never recomputes)
# Manuscript guard set: render, first-use terms, sentence-extractor selftest, build graph,
# prose integrity, unresolved refs, page counts, main-text and abstract word counts.
cd "$(dirname "$0")/.."
DESURV_RECOMPUTE=FALSE timeout 480 make paper > "${GUARDS_LOG:-/tmp/desurv-render.log}" 2>&1; rc=$?; echo "  render exit=$rc"
[ $rc = 0 ] || { L=$(ls -t paper/*.log 2>/dev/null|head -1); [ -n "$L" ] && grep -n -A6 '^! ' "$L" | head -20; tail -20 "${GUARDS_LOG:-/tmp/desurv-render.log}"; exit 1; }
rm -f paper/paper.aux paper/si_appendix.aux paper/si_appendix.lof paper/si_appendix.lot paper/si_appendix.toc
make check-terms 2>&1 | tail -1; python3 code/util_sentence_extract.py --selftest; scripts/check_build_graph.sh | tail -1; scripts/check_prose_integrity.sh | tail -1
for f in paper si_appendix; do echo "  $f: pages=$(pdfinfo paper/$f.pdf|awk '/Pages/{print $2}') unresolved=$(pdftotext paper/$f.pdf -|tr -s '[:space:]' ' '|grep -o '??'|wc -l)"; done
echo "  refs in sources: $(grep -o '\\ref{' paper/*.Rmd | wc -l)"
echo "  main text: $(t=0; for f in 02_introduction_REVISED 04_results_REVISED 05_discussion_REVISED; do n=$(python3 code/util_sentence_extract.py paper/$f.Rmd --list | grep -oP '^\S+\s+\S+\s+\K\d+(?=w)' | paste -sd+ | bc); t=$((t+n)); done; echo $t) / 4000"
echo "  abstract (rendered): $(pdftotext -f 1 -l 1 paper/paper.pdf - | tr -s '[:space:]' ' ' | python3 -c "
import sys; t=sys.stdin.read(); i=t.find('Pancreatic ductal adenocarcinoma (PDAC) outcome'); j=t.find('prognostically relevant, interpretable programs.')+len('prognostically relevant, interpretable programs.'); print(len(t[i:j].split()))") words"
