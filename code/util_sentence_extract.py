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


def strip_markup(text: str) -> str:
    """Reduce Rmd/LaTeX markup to plain prose without changing sentence count."""
    text = re.sub(r"`r [^`]*`", "<VAL>", text)          # inline R -> single token
    text = re.sub(r"\$[^$]*\$", "<MATH>", text)          # inline math -> single token
    text = re.sub(r"\\ref\{[^}]*\}", "<REF>", text)
    text = re.sub(r"\\label\{[^}]*\}", "", text)
    text = re.sub(r"\\[a-zA-Z]+\s*", " ", text)          # residual LaTeX commands
    text = re.sub(r"[*_]{1,2}([^*_]+)[*_]{1,2}", r"\1", text)  # emphasis
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def paragraphs(path: str):
    """Yield (line_no, raw_text) for prose paragraphs, skipping code chunks."""
    out, buf, start, in_chunk, in_comment = [], [], None, False, False
    for i, raw in enumerate(open(path, encoding="utf-8"), 1):
        line = raw.rstrip("\n")
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
        # Section headers and pure-LaTeX structural lines are not prose.
        is_struct = stripped.startswith(("#", "\\section", "\\subsection", "\\begin",
                                         "\\end", "\\bottomrule", "\\midrule",
                                         "\\toprule", "\\addlinespace", "<!--"))
        if not stripped or is_struct:
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


def main():
    ap = argparse.ArgumentParser()
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
