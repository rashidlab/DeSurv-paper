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


def replace_cmd_n(s, cmd, nargs, pick):
    """Replace `cmd{a}{b}{c}` consuming nargs brace groups, keeping pick(args).

    Needed because \\multicolumn takes three arguments. Consuming only the first
    left "{l}{Training total}" as literal text in the SI cohort table.
    """
    out, i = [], 0
    pat = re.compile(re.escape(cmd) + r"\s*(?=\{)")
    while True:
        m = pat.search(s, i)
        if not m:
            out.append(s[i:])
            break
        out.append(s[i:m.start()])
        j, args = m.end(), []
        for _ in range(nargs):
            if j >= len(s) or s[j] != "{":
                break
            e = match_brace(s, j)
            args.append(s[j + 1:e - 1])
            j = e
        out.append(pick(args))
        i = j
    return "".join(out)


def flatten_math(s):
    """equation/align/gather bodies become $$...$$ so Word renders them as
    equations instead of showing \\begin{equation} as literal source."""
    def repl(m):
        body = m.group(3)
        body = replace_cmd(body, r"\label", lambda x: "")
        body = re.sub(r"\\notag|\\nonumber", "", body)
        # pandoc's math reader does not know \mathbbm (bbm package).
        body = re.sub(r"\\mathbbm\b", r"\\mathbf", body)
        body = body.strip()
        if not body:
            return ""
        if r"\\" in body:
            # A bare \\ is not legal inside $$...$$; aligned is, and pandoc
            # converts it to a multi-line Word equation. Keep the & alignment
            # points in that case, strip them otherwise.
            return "\n\n$$\n\\begin{aligned}\n%s\n\\end{aligned}\n$$\n\n" % body
        body = re.sub(r"(?<!\\)&", " ", body)
        return "\n\n$$\n%s\n$$\n\n" % body
    return re.sub(r"\\begin\{(equation|align|gather|multline|displaymath)(\*?)\}"
                  r"(.*?)"
                  r"\\end\{\1\2\}",
                  repl, s, flags=re.S)


def flatten_algorithms(s):
    """algorithm/algorithmic bodies become an indented plain-text block.
    Readable prose beats correct-looking LaTeX source in a DOCX."""
    def repl(m):
        body = m.group(0)
        title = ""
        caps = take_arg(body, r"\caption")
        if caps:
            title = "**" + clean_inline(caps[0][2]).strip() + "**\n\n"
        inner = re.search(r"\\begin\{algorithmic\}(?:\[[^\]]*\])?(.*?)\\end\{algorithmic\}",
                          body, re.S)
        body = inner.group(1) if inner else body
        lines = []
        for ln in body.splitlines():
            ln = ln.strip()
            if not ln or ln.startswith(("\\begin", "\\end", "\\caption", "\\label")):
                continue
            ln = re.sub(r"^\\Require\b", "Input:", ln)
            ln = re.sub(r"^\\Ensure\b", "Output:", ln)
            ln = re.sub(r"^\\State(x)?\b", "", ln)
            ln = replace_cmd(ln, r"\Comment", lambda x: "  // " + x)
            # Control words take a braced condition; leaving the braces printed
            # "while{$eps \geq tol$ ...}" in the DOCX.
            for cmd, pre, post in ((r"\While", "while ", " do"),
                                   (r"\For", "for ", " do"),
                                   (r"\If", "if ", " then"),
                                   (r"\ElsIf", "else if ", " then")):
                ln = replace_cmd(ln, cmd, lambda x, p=pre, q=post: p + x + q)
            ln = re.sub(r"^\\End(While|For|If|Function|Procedure)\b",
                        lambda mm: "end " + mm.group(1).lower(), ln)
            ln = re.sub(r"^\\(Else|Repeat|Until)\b", lambda mm: mm.group(1).lower(), ln)
            ln = clean_inline(ln).strip()   # math-protecting; _clean_text
                                            # would gut $\mathbb{R}$ etc.
            if ln:
                # NOT indented: four leading spaces make markdown treat the
                # block as code, which leaves $...$ unrendered and prints
                # literal ** markers. Two trailing spaces is a hard line break.
                lines.append(ln + "  ")
        return "\n\n" + title + "\n".join(lines) + "\n\n" if lines else ""
    return re.sub(r"\\begin\{algorithm\}(?:\[[^\]]*\])?.*?\\end\{algorithm\}",
                  repl, s, flags=re.S)


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
    # \cmidrule takes an optional parenthesised trim plus a brace range; the
    # generic command-stripper removed the name and left "(l{3pt}r{3pt}){1-6}".
    body = re.sub(r"\\cmidrule\s*(?:\([^)]*\))?\s*\{[^}]*\}", "", body)
    # \multicolumn takes THREE arguments; keep only the text.
    body = replace_cmd_n(body, r"\multicolumn", 3,
                         lambda a: a[2] if len(a) == 3 else "")
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
    # \( ... \) is inline math too; without this the protection below misses it
    # and the command stripper eats \beta out of "\(\beta\)-update".
    s = s.replace(r"\(", "$").replace(r"\)", "$")
    # Spacing commands take an argument. Stripping them by name alone left
    # "{0.5cm}" printed in the DOCX, the same defect as \multicolumn.
    s = re.sub(r"\\[hv]space\*?\s*\{[^}]*\}", " ", s)
    # Mask math with placeholders rather than SPLITTING on it. A command can
    # span a math delimiter -- \textbf{($\beta$-update)} -- and splitting first
    # handed _clean_text the unbalanced fragment "\textbf{(", whose brace match
    # ran to end-of-fragment and emitted stray ** markers into the DOCX.
    math = []

    def stash(m):
        math.append(m.group(0))
        return "\x00%d\x00" % (len(math) - 1)

    s = _clean_text(re.sub(r"\$[^$]*\$", stash, s))
    return re.sub(r"\x00(\d+)\x00", lambda m: math[int(m.group(1))], s)


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


def flatten_prose(s):
    head, sep, body = s.partition("\n---\n")  # keep the YAML header intact
    if not sep:
        head, body = "", s
    math = []

    def stash(m):
        math.append(m.group(0))
        return "\x00%d\x00" % (len(math) - 1)
    body = re.sub(r"\$\$.*?\$\$|\$[^$\n]*\$", stash, body, flags=re.S)
    body = re.sub(r"\\(part|section|subsection|subsubsection)\*\{", lambda m: "\\%s{" % m.group(1), body)  # starred forms
    body = replace_cmd(body, r"\part", lambda x: "\n\n# %s\n" % x)
    body = replace_cmd(body, r"\section", lambda x: "\n\n# %s\n" % x)
    body = replace_cmd(body, r"\subsection", lambda x: "\n\n## %s\n" % x)
    body = replace_cmd(body, r"\subsubsection", lambda x: "\n\n### %s\n" % x)
    body = replace_cmd(body, r"\textbf", lambda x: "**%s**" % x)
    body = replace_cmd(body, r"\textit", lambda x: "*%s*" % x)
    body = replace_cmd(body, r"\emph", lambda x: "*%s*" % x)
    body = replace_cmd(body, r"\texttt", lambda x: "`%s`" % x)
    body = replace_cmd(body, r"\url", lambda x: "<%s>" % x)
    body = replace_cmd(body, r"\textsubscript", lambda x: "~%s~" % x)
    body = replace_cmd(body, r"\textsuperscript", lambda x: "^%s^" % x)
    body = re.sub(r"\\addcontentsline\{[^}]*\}\{[^}]*\}\{[^}]*\}", "", body)
    body = re.sub(r"\\(newpage|clearpage|noindent|centering)\b", "", body)
    body = re.sub(r"\x00(\d+)\x00", lambda m: math[int(m.group(1))], body)
    return head + sep + body


def convert(base):
    s = open(base + ".knit.md", encoding="utf8").read()

    # Tables first: they contain \caption too, and must not be eaten by the
    # figure pass.
    def table_repl(m):
        return latex_table_to_markdown(m.group(0))
    s = re.sub(r"\\begin\{table\*?\}.*?\\end\{table\*?\}", table_repl, s, flags=re.S)

    # Algorithms before math: an algorithm body can contain display math.
    s = flatten_algorithms(s)
    s = flatten_math(s)

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

    # Prose written in LaTeX rather than markdown (the SI mixes both): section
    # commands and inline formatting outside math. pandoc is run with -raw_tex,
    # so anything left here would print verbatim in the DOCX. Real incident
    # 2026-09-06: the R5/R6 SI sections used \subsection*{} and \texttt{} and
    # appeared as raw markup in the editable file. Math is masked so that
    # commands inside $...$ are untouched.
    s = flatten_prose(s)

    md = base + ".docx.md"
    open(md, "w", encoding="utf8").write(s)

    cmd = ["pandoc", md, "-o", base + ".docx", "--citeproc",
           "--bibliography=references_30102025.bib", "--csl=nature.csl",
           "--resource-path=.:..:../figures:figures:" + MEDIA,
           "--from", "markdown+tex_math_dollars-raw_tex"]
    subprocess.run(cmd, check=True)
    return len(re.findall(r"!\[\]\(", s)), s.count("|---"), missing


if __name__ == "__main__":
    failed = False
    for b in sys.argv[1:]:
        n_img, n_tab, miss = convert(b)
        print("%s.docx: %d images, %d tables" % (b, n_img, n_tab))
        if miss:
            # A converter that reports a dropped figure and exits 0 reintroduces
            # exactly the silent-loss mode this script exists to prevent.
            sys.stderr.write("ERROR: %s: unresolved figure(s), dropped from the "
                             "DOCX: %s\n" % (b, ", ".join(miss)))
            failed = True
    sys.exit(1 if failed else 0)
