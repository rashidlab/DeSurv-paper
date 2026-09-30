#!/bin/bash
# render + full guard set; prints one line per check; exits nonzero if render or any guard fails
# Usage: scripts/guards.sh   (from anywhere; renders with DESURV_RECOMPUTE=FALSE, never recomputes)
# Manuscript guard set: render, first-use terms, sentence-extractor selftest, build graph,
# prose integrity, abstract HRs vs cache, unresolved refs, page counts, main-text and abstract word counts.
cd "$(dirname "$0")/.."
DESURV_RECOMPUTE=FALSE timeout 480 make paper > "${GUARDS_LOG:-/tmp/desurv-render.log}" 2>&1; rc=$?; echo "  render exit=$rc"
[ $rc = 0 ] || { L=$(ls -t paper/*.log 2>/dev/null|head -1); [ -n "$L" ] && grep -n -A6 '^! ' "$L" | head -20; tail -20 "${GUARDS_LOG:-/tmp/desurv-render.log}"; exit 1; }
rm -f paper/paper.aux paper/si_appendix.aux paper/si_appendix.lof paper/si_appendix.lot paper/si_appendix.toc
fail=0
make check-terms 2>&1 | tail -1; [ ${PIPESTATUS[0]} = 0 ] || fail=1
python3 code/util_sentence_extract.py --selftest || fail=1
scripts/check_build_graph.sh | tail -1; [ ${PIPESTATUS[0]} = 0 ] || fail=1
scripts/check_prose_integrity.sh | tail -1; [ ${PIPESTATUS[0]} = 0 ] || fail=1
Rscript scripts/check_abstract_hr.R || fail=1
for f in paper si_appendix; do u=$(pdftotext paper/$f.pdf -|tr -s '[:space:]' ' '|grep -o '??'|wc -l); echo "  $f: pages=$(pdfinfo paper/$f.pdf|awk '/Pages/{print $2}') unresolved=$u"; [ "$u" = 0 ] || fail=1; done
echo "  refs in sources: $(grep -o '\\ref{' paper/*.Rmd | wc -l)"
echo "  main text: $(scripts/main_text_words.sh) / 4000"
echo "  abstract (rendered): $(pdftotext -f 1 -l 1 paper/paper.pdf - | tr -s '[:space:]' ' ' | python3 -c "
import sys; t=sys.stdin.read(); i=t.find('Pancreatic ductal adenocarcinoma (PDAC) outcome'); j=t.find('survival signal contributed little expression variation.')+len('survival signal contributed little expression variation.'); print(len(t[i:j].split()))") words"
[ $fail = 0 ] && echo "  guards: PASS" || echo "  guards: FAIL"
exit $fail
