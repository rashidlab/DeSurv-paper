#!/usr/bin/env bash
# LESSONS.md L11: a reference converted into plain text is invisible to the
# unresolved-reference check. Fails on any carriage return in a manuscript
# source (the signature of a regex replacement that ate a backslash) and on the
# literal "ef{" in rendered text. Prints \ref counts so a net loss across an
# edit batch is visible.
set -u; cd "$(dirname "$0")/.." >/dev/null
fail=0
cr=$(cat paper/*.Rmd | tr -cd '\r' | wc -c); [ "$cr" = 0 ] || { echo "FAIL: $cr carriage return(s) in paper/*.Rmd"; fail=1; }
for p in paper/paper.pdf paper/si_appendix.pdf; do [ -f $p ] || continue; n=$(pdftotext $p - 2>/dev/null | grep -c 'ef{'); [ "$n" = 0 ] || { echo "FAIL: literal 'ef{' in $p ($n)"; fail=1; }; done
echo "ref counts: $(for f in paper/0*.Rmd paper/si_appendix.Rmd; do printf '%s=%s ' "$(basename $f .Rmd | cut -c1-2)" "$(grep -o '\\ref{' $f | wc -l)"; done)"
[ $fail = 0 ] && echo "PASS: no carriage returns, no plain-text refs"
exit $fail
