#!/usr/bin/env python3
"""Convert a knitr .knit.md (LaTeX-flavoured) into a DOCX via pandoc.

The manuscript targets pdf_document, so knitr emits LaTeX floats
(\\begin{figure*} ... \\includegraphics ... \\caption{...}) and kable emits
LaTeX tabulars. Pandoc's markdown reader treats those as raw LaTeX and drops
them from a DOCX, which silently yields a file with the prose but no figures,
no captions and no tables. This flattens them first.

Figure PDFs are rasterised to PNG because DOCX cannot embed PDF images.

Usage: util_knit_to_docx.py <base>   # reads <base>.knit.md, writes <base>.docx
Run from the directory containing the .knit.md (i.e. paper/).
"""
import os
import re
import subprocess
import sys

MEDIA = ".docx_media"


def match_brace(s, i):
    """Index just past the '}' matching the '{' at s[i]."""
    depth = 0
    while i < len(s):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return len(s)


def take_arg(s, cmd):
    """Yield (start, end, inner) for each `cmd{...}`, brace-balanced."""
    out = []
    for m in re.finditer(re.escape(cmd) + r"\s*\{", s):
        b = s.index("{", m.start())
        e = match_brace(s, b)
        out.append((m.start(), e, s[b + 1:e - 1]))
    return out


def replace_cmd(s, cmd, fmt):
    for start, end, inner in reversed(take_arg(s, cmd)):
        s = s[:start] + fmt(inner) + s[end:]
    return s


def resolve(p):
    for cand in (p, p + ".pdf", p + ".png", p + ".jpg"):
        if os.path.exists(cand):
            return cand
    return None


def to_png(src):
    os.makedirs(MEDIA, exist_ok=True)
    base = re.sub(r"[^A-Za-z0-9]+", "_", os.path.splitext(src)[0]).strip("_")
    out = os.path.join(MEDIA, base + ".png")
    if not os.path.exists(out):
        if src.lower().endswith(".pdf"):
            subprocess.run(["pdftoppm", "-png", "-r", "200", "-singlefile",
                            src, out[:-4]], check=True)
        else:
            subprocess.run(["cp", src, out], check=True)
    return out


def latex_table_to_markdown(block):
    """Flatten a kable LaTeX tabular into a pipe table. Content over fidelity."""
    cap = ""
    caps = take_arg(block, r"\caption")
    if caps:
        cap = caps[0][2]
    # kableExtra emits \begin{tabular}[t]{lll}; the optional placement argument
    # must be allowed or the match fails and the table is silently deleted.
    m = re.search(r"\\begin\{tabular\}(?:\[[^\]]*\])?\{[^}]*\}(.*?)\\end\{tabular\}",
                  block, re.S)
    if not m:
        sys.stderr.write("WARNING: table block matched but no tabular found; "
                         "kept verbatim so nothing is lost\n")
        return "\n```\n" + block.strip() + "\n```\n"
    body = m.group(1)
    for junk in (r"\toprule", r"\midrule", r"\bottomrule", r"\addlinespace"):
        body = body.replace(junk, "")
    rows = [r.strip() for r in body.split(r"\\") if r.strip()]
    cells = [[clean_inline(c).strip() for c in r.split("&")] for r in rows]
    cells = [r for r in cells if any(r)]
    if not cells:
        return ""
    ncol = max(len(r) for r in cells)
    cells = [r + [""] * (ncol - len(r)) for r in cells]
    out = ["", "| " + " | ".join(cells[0]) + " |",
           "|" + "|".join(["---"] * ncol) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in cells[1:]]
    if cap:
        out += ["", "*" + clean_inline(cap).strip() + "*"]
    return "\n".join(out) + "\n"


def clean_inline(s):
    """Strip LaTeX markup OUTSIDE math only. Running the command-stripping
    regex over $...$ turns "$\\alpha \\in [0,1)$" into "$ [0,1)$", which pandoc
    then fails to parse and reports as a math error."""
    parts = re.split(r"(\$[^$]*\$)", s)
    return "".join(p if (p.startswith("$") and p.endswith("$") and len(p) > 1)
                   else _clean_text(p) for p in parts)


def _clean_text(s):
    s = replace_cmd(s, r"\textbf", lambda x: "**%s**" % x)
    s = replace_cmd(s, r"\textit", lambda x: "*%s*" % x)
    s = replace_cmd(s, r"\texttt", lambda x: "`%s`" % x)
    s = replace_cmd(s, r"\label", lambda x: "")
    s = replace_cmd(s, r"\ref", lambda x: "[%s]" % x)
    s = replace_cmd(s, r"\multicolumn", lambda x: "")
    s = re.sub(r"\\%", "%", s)
    s = re.sub(r"\\ ", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?", "", s)
    return re.sub(r"[ \t]{2,}", " ", s)


def convert(base):
    s = open(base + ".knit.md", encoding="utf8").read()

    # Tables first: they contain \caption too, and must not be eaten by the
    # figure pass.
    def table_repl(m):
        return latex_table_to_markdown(m.group(0))
    s = re.sub(r"\\begin\{table\*?\}.*?\\end\{table\*?\}", table_repl, s, flags=re.S)

    # Figures: rasterise, keep the caption as an italic paragraph.
    missing = []

    def img_repl(m):
        src = resolve(m.group(2))
        if src is None:
            missing.append(m.group(2))
            return ""
        return "![](%s)" % to_png(src)
    s = re.sub(r"\\includegraphics(\[[^\]]*\])?\{([^}]*)\}", img_repl, s)
    s = replace_cmd(s, r"\caption", lambda x: "\n\n*" + clean_inline(x).strip() + "*\n")
    for env in ("figure\\*", "figure"):
        s = re.sub(r"\\begin\{%s\}(\[[^\]]*\])?" % env, "", s)
        s = re.sub(r"\\end\{%s\}" % env, "", s)
    s = re.sub(r"\{\\centering\s*", "", s)
    s = re.sub(r"^\}\s*$", "", s, flags=re.M)

    # Cross-references cannot be numbered outside LaTeX; keep the target
    # visible rather than dropping it silently.
    s = replace_cmd(s, r"\ref", lambda x: "[%s]" % x)
    s = replace_cmd(s, r"\label", lambda x: "")
    s = re.sub(r"\\externaldocument\{[^}]*\}", "", s)

    md = base + ".docx.md"
    open(md, "w", encoding="utf8").write(s)

    cmd = ["pandoc", md, "-o", base + ".docx", "--citeproc",
           "--bibliography=references_30102025.bib", "--csl=nature.csl",
           "--resource-path=.:..:../figures:figures:" + MEDIA,
           "--from", "markdown+tex_math_dollars-raw_tex"]
    subprocess.run(cmd, check=True)
    return len(re.findall(r"!\[\]\(", s)), s.count("|---"), missing


if __name__ == "__main__":
    for b in sys.argv[1:]:
        n_img, n_tab, miss = convert(b)
        print("%s.docx: %d images, %d tables%s"
              % (b, n_img, n_tab, ", MISSING %s" % miss if miss else ""))
