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
# Stale manuscript title: every tracked source that states the title must match paper.Rmd.
# 2026-09-06: the SI subtitle kept the old title after the title changed, and a context
# regex that required 45 characters before the word missed the short subtitle line.
python3 - <<'PYCHK' || fail=1
import re
norm = lambda t: re.sub(r"\s+", " ", t.replace("--", "-").replace("{", "").replace("}", "")).strip()
t = norm(re.search(r'^title:\s*"(.*)"\s*$', open("paper/paper.Rmd").read(), re.M).group(1))
bad = []
si = re.search(r'^subtitle:\s*"(.*)"\s*$', open("paper/si_appendix.Rmd").read(), re.M)
if si and norm(si.group(1)) != t: bad.append("paper/si_appendix.Rmd subtitle")
rd = open("README.md").read()
h = rd.splitlines()[0]
if h.startswith("# ") and norm(h[2:]) != t: bad.append("README.md heading")
bib = re.search(r"title\s*=\s*\{(.*?)\},[ \t]*\n", rd, re.S)   # to the field END, not the first "}," (which "{DeSurv}," contains)
if bib and norm(bib.group(1)) != t: bad.append("README.md bibtex title")
if bad:
    print("FAIL: stale manuscript title in: " + ", ".join(bad)); raise SystemExit(1)
print("title consistent across paper.Rmd, si_appendix.Rmd, README.md")
PYCHK
echo "ref counts: $(for f in paper/0*.Rmd paper/si_appendix.Rmd; do printf '%s=%s ' "$(basename $f .Rmd | cut -c1-2)" "$(grep -o '\\ref{' $f | wc -l)"; done)"
[ $fail = 0 ] && echo "PASS: no carriage returns, no plain-text refs"
exit $fail
