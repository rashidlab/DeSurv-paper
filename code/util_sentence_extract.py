#!/usr/bin/env python3
"""Extract manuscript prose one sentence per line with word and clause counts.

Supports the line-by-line pass: batches of paragraphs, one sentence per line,
so that over-long sentences, clause pile-ups and repeated frames are visible as
data rather than having to be noticed while reading. Reading the rendered PDF
does not substitute, because reflow hides sentence length.

Operates on the .Rmd SOURCE, not the PDF, so that line numbers are editable
addresses. Inline R (`r ...`) is collapsed to a <VAL> token so that a computed
number does not inflate the word count or split a sentence at its decimal point.

Usage:
  util_sentence_extract.py <file.Rmd> [--start N] [--batch N] [--min-words N]
  util_sentence_extract.py <file.Rmd> --list      # paragraph index only
"""
import argparse
import re
import sys

# Sentence-final punctuation followed by whitespace + capital/quote/paren.
# Guarded against the abbreviations and the "Fig." / "et al." forms that occur
# in this manuscript, plus decimals, which would otherwise split mid-number.
ABBREV = r"(?<!\bFig)(?<!\bSupplementary Fig)(?<!\bet al)(?<!\bcf)(?<!\be\.g)(?<!\bi\.e)(?<!\bvs)(?<!\bapprox)(?<!\bRef)(?<!\bNo)(?<!\bDr)"
# No digit or capital-letter lookbehind here. Both were tried and both
# over-blocked: this manuscript ends sentences on factor names ("...recover
# D1.") and acronyms ("...resolved by NMF."), and either guard silently merged
# two sentences into one, which shows up as a fake over-length finding. A
# decimal needs no guard, because "3.5" has no whitespace after the period.
# Single-letter initials ("R. A. Moffitt") are re-merged after the split.
SENT_SPLIT = re.compile(rf"{ABBREV}([.!?])\s+(?=[A-Z\"'(\\$])")
INITIAL_END = re.compile(r"(?:^|[\s(])[A-Z]\.$")

CLAUSE_MARKERS = re.compile(
    r"\b(although|whereas|while|because|since|if|unless|when|where|which|that|who|"
    r"after|before|so that|such that|rather than|instead of|and|but|yet|however|"
    r"therefore|thus|consequently)\b",
    re.I,
)


# LaTeX environments whose BODY is not prose. The original version skipped only
# the \begin and \end lines, which is why si_appendix.Rmd reported 46 \State
# pseudocode lines and every display-math body as prose paragraphs.
SKIP_ENVS = {
    "equation", "align", "gather", "multline", "eqnarray", "displaymath",
    "algorithmic", "algorithm", "tabular", "tabularx", "table", "figure",
    "verbatim", "lstlisting", "array", "cases", "split", "pmatrix", "bmatrix",
}
ENV_BEGIN = re.compile(r"\\begin\{([a-zA-Z]+)\*?\}")
ENV_END = re.compile(r"\\end\{([a-zA-Z]+)\*?\}")

# Structural lines that are not prose. \part was absent, so
# "\part{Supplementary Methods}" survived as a two-word prose paragraph reading
# "{Supplementary Methods}"; the pseudocode and float-furniture commands were
# absent for the same reason.
STRUCT_PREFIXES = (
    "#", "<!--",
    "\\part", "\\chapter", "\\section", "\\subsection", "\\subsubsection",
    "\\paragraph", "\\begin", "\\end",
    "\\toprule", "\\midrule", "\\bottomrule", "\\addlinespace", "\\hline",
    "\\cmidrule", "\\multicolumn", "\\caption", "\\label", "\\centering",
    "\\listoftables", "\\listoffigures", "\\tableofcontents",
    "\\newpage", "\\clearpage", "\\pagebreak", "\\maketitle", "\\appendix",
    "\\phantomsection", "\\addcontentsline",
    "\\input", "\\include", "\\includegraphics", "\\vspace", "\\hspace",
    "\\State", "\\Statex", "\\Require", "\\Ensure", "\\While", "\\EndWhile",
    "\\For", "\\EndFor", "\\If", "\\ElsIf", "\\Else", "\\EndIf",
    "\\Function", "\\EndFunction", "\\Procedure", "\\EndProcedure",
    "\\Return", "\\Comment", "\\Repeat", "\\Until",
)


def strip_markup(text: str) -> str:
    """Reduce Rmd/LaTeX markup to plain prose without changing sentence count."""
    text = re.sub(r"`r [^`]*`", "<VAL>", text)          # inline R -> single token
    text = re.sub(r"\$[^$]*\$", "<MATH>", text)          # inline math -> single token
    # Citation groups collapse to one token. Left expanded, a five-key group
    # adds five words to the count and can trip the length thresholds this
    # tool exists to audit; the reader sees one superscript run, not five words.
    text = re.sub(r"\[[^\]]*@[^\]]*\]", "<CITE>", text)
    text = re.sub(r"\\ref\{[^}]*\}", "<REF>", text)
    text = re.sub(r"\\label\{[^}]*\}", "", text)
    text = re.sub(r"\\[a-zA-Z]+\s*", " ", text)          # residual LaTeX commands
    text = re.sub(r"[*_]{1,2}([^*_]+)[*_]{1,2}", r"\1", text)  # emphasis
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def paragraphs(path: str):
    """Yield (line_no, raw_text) for prose paragraphs, skipping code chunks,
    YAML front matter, display math, and non-prose LaTeX environment bodies."""
    out, buf, start = [], [], None
    in_chunk = in_comment = in_yaml = in_display = False
    env_depth = 0
    lines = open(path, encoding="utf-8").readlines()
    # YAML front matter: only when the very first line opens it. si_appendix.Rmd
    # and paper.Rmd both have one; the child .Rmd files do not.
    if lines and lines[0].strip() == "---":
        in_yaml = True

    for i, raw in enumerate(lines, 1):
        line = raw.rstrip("\n")
        if in_yaml:
            if i > 1 and line.strip() in ("---", "..."):
                in_yaml = False
            continue
        if line.startswith("```"):
            in_chunk = not in_chunk
            continue
        if in_chunk:
            continue
        # HTML comments span lines here (the SI figure-order block), so track
        # the open/close state rather than testing for a leading marker.
        if in_comment:
            if "-->" in line:
                in_comment = False
            continue
        if line.lstrip().startswith("<!--") and "-->" not in line:
            in_comment = True
            continue

        stripped = line.strip()

        # $$ display math, which may span lines.
        if stripped.startswith("$$"):
            if stripped.count("$$") == 1:
                in_display = not in_display
            continue
        if in_display:
            continue

        # Environment BODIES, not just their delimiters. Tracked by depth so a
        # nested align inside an algorithm does not close the outer one early.
        if env_depth:
            if ENV_END.search(line) and ENV_END.search(line).group(1) in SKIP_ENVS:
                env_depth -= 1
            continue
        m = ENV_BEGIN.search(line)
        if m and m.group(1) in SKIP_ENVS:
            env_depth += 1
            if buf:
                out.append((start, " ".join(buf)))
                buf, start = [], None
            continue

        if not stripped or stripped.startswith(STRUCT_PREFIXES):
            if buf:
                out.append((start, " ".join(buf)))
                buf, start = [], None
            continue
        if not buf:
            start = i
        buf.append(stripped)
    if buf:
        out.append((start, " ".join(buf)))
    return out


def sentences(text: str):
    text = strip_markup(text)
    if not text:
        return []
    parts = SENT_SPLIT.split(text)
    # re.split with one capture group interleaves [body, punct, body, punct, ...]
    merged, i = [], 0
    while i < len(parts):
        body = parts[i]
        punct = parts[i + 1] if i + 1 < len(parts) else ""
        s = (body + punct).strip()
        if s:
            # Re-join an initial ("R.") to the sentence it belongs to.
            if merged and INITIAL_END.search(merged[-1]):
                merged[-1] = merged[-1] + " " + s
            else:
                merged.append(s)
        i += 2
    return merged


# Fixture reproducing the three incidents this tool has actually mis-reported on
# si_appendix.Rmd: YAML front matter counted as a prose paragraph, \part
# surviving as "{Supplementary Methods}", and algorithmic/equation BODIES
# counted as prose because only the \begin and \end lines were skipped.
# Perturbation is written from the incidents, not from whatever breaks easiest.
SELFTEST = '''---
title: "Supplementary Information"
output: pdf_document
---

\\listoftables
\\newpage

\\part{Supplementary Methods}

\\section{Model details}
This is the only real prose paragraph in the fixture. It has two sentences.

\\begin{algorithm}
\\begin{algorithmic}
\\Require Expression matrix and survival outcomes
\\State Initialize W and H nonnegatively
\\Ensure Fitted gene programs
\\end{algorithmic}
\\end{algorithm}

\\begin{equation}
\\mathcal{L} = (1-\\alpha)\\,\\mathcal{L}_{\\mathrm{NMF}} - \\alpha\\,\\mathcal{L}_{\\mathrm{Cox}}
\\end{equation}
'''


def selftest() -> int:
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".Rmd")
    with os.fdopen(fd, "w") as fh:
        fh.write(SELFTEST)
    try:
        paras = paragraphs(path)
    finally:
        os.unlink(path)
    ok = True
    if len(paras) != 1:
        print(f"FAIL: expected 1 prose paragraph, got {len(paras)}")
        for ln, txt in paras:
            print(f"  L{ln}: {strip_markup(txt)[:70]}")
        ok = False
    elif "only real prose paragraph" not in paras[0][1]:
        print(f"FAIL: wrong paragraph captured: {paras[0][1][:70]}")
        ok = False
    else:
        nw = len(strip_markup(paras[0][1]).split())
        if nw != 14:
            print(f"FAIL: expected 14 words, got {nw}")
            ok = False
    print("selftest: PASS" if ok else "selftest: FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    ap.add_argument("file")
    ap.add_argument("--start", type=int, default=1, help="first paragraph (1-based)")
    ap.add_argument("--batch", type=int, default=10, help="paragraphs per batch")
    ap.add_argument("--min-words", type=int, default=0, help="only show sentences >= N words")
    ap.add_argument("--list", action="store_true", help="print paragraph index only")
    args = ap.parse_args()

    paras = paragraphs(args.file)

    if args.list:
        print(f"{len(paras)} prose paragraphs in {args.file}\n")
        for n, (ln, txt) in enumerate(paras, 1):
            flat = strip_markup(txt)
            print(f"P{n:<3} L{ln:<5} {len(flat.split()):>4}w  {flat[:88]}")
        return

    lo = args.start - 1
    hi = min(lo + args.batch, len(paras))
    print(f"{args.file}: paragraphs {args.start}-{hi} of {len(paras)}\n")

    for n in range(lo, hi):
        line_no, txt = paras[n]
        sents = sentences(txt)
        total = sum(len(s.split()) for s in sents)
        print(f"=== P{n+1}  (source line {line_no})  {len(sents)} sentences, {total} words")
        for j, s in enumerate(sents, 1):
            w = len(s.split())
            c = len(CLAUSE_MARKERS.findall(s))
            if w < args.min_words:
                continue
            flag = ""
            if w >= 40:
                flag = "  <<< LONG"
            elif w >= 30:
                flag = "  <<< long"
            if c >= 5:
                flag += f"  <<< {c} clause markers"
            print(f"  {j:>2}. [{w:>3}w {c}c]{flag}")
            print(f"      {s}")
        print()


if __name__ == "__main__":
    main()
