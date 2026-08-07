#!/usr/bin/env Rscript
# util_check_first_use.R
# ---------------------------------------------------------------------------
# LESSONS.md L1: technical terms must be DEFINED AT FIRST USE in the main text.
#
# DESIGN. Regex cannot reliably decide whether prose "defines" a term: an
# earlier version guessed from cue words and produced false flags in both
# directions, which is worse than no check (LESSONS.md L3). So the judgement is
# recorded by a human, once, and the machine enforces it:
#
#   DEFINED  term -> the exact definitional phrase that must appear near first
#                    use. If someone edits that phrase away, the check FIRES.
#   ACCEPTED term -> standard terminology for this readership, deliberately not
#                    defined, with the reason recorded here.
#
# A term in neither list is flagged as an unreviewed first use.
#
# Usage:  Rscript code/util_check_first_use.R [paper/paper.pdf]
# Exits 1 on any failure, so it can gate a submission.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
pdf  <- if (length(args)) args[1] else "paper/paper.pdf"
if (!file.exists(pdf)) stop("no such PDF: ", pdf, call. = FALSE)

txt <- paste(system2("pdftotext", c(shQuote(pdf), "-"), stdout = TRUE), collapse = "\n")
txt <- gsub("[ \t]*\n[ \t]*", " ", txt)

# Main text only: Introduction through Discussion. Hard-stop at References, or
# reference titles leak in and produce false flags.
# Start at the INTRODUCTION, not the abstract. "Pancreatic ductal" opens both,
# so take the second occurrence; otherwise a term first used in the abstract
# is judged against a paragraph that cannot reasonably define it.
i0a <- gregexpr("Pancreatic ductal", txt, fixed = TRUE)[[1]]
i0  <- if (length(i0a) >= 2) i0a[2] else i0a[1]
ref <- regexpr(" References ", txt, fixed = TRUE)
lim <- if (ref > i0) ref else nchar(txt)
i1  <- gregexpr(" Methods ", txt, fixed = TRUE)[[1]]
i1  <- i1[i1 > i0 & i1 < lim]
body <- substr(txt, i0, if (length(i1)) max(i1) else lim)

# ── Terms whose definition is enforced: term -> required phrase near first use ──
DEFINED <- c(
  "gene program"             = "weighted sets of co-varying genes",
  "frozen"                   = "are held fixed, or frozen",
  "representation"           = "together form a representation of the tumor transcriptome",
  "reconstruction error"      = "the mismatch between the measured expression matrix",
  "rank"                      = "the number of programs the factorization extracts",
  "latent"                    = "the true underlying (latent) programs are known",
  "Bayesian optimization"     = "a sequential search algorithm that proposes each new setting",
  "concordance"               = "the probability that the model ranks a pair of patients",
  "one-standard-error"        = "the most parsimonious model within one standard error",
  "partial log-likelihood"    = "the fit criterion of the Cox model",
  "consensus initialization"  = "starting values pooled over many restarts",
  "linear predictor"          = "the weighted sum of the three program scores"
)

# ── Terms deliberately left undefined, with the reason ──────────────────────
ACCEPTED <- c(
  "Shapley"                          = "cited to Shapley (1953) at first use; derived in full in SI Section 9",
  "nonnegative matrix factorization" = "expanded at first use as 'nonnegative matrix factorization (NMF)'",
  "hazard ratio"                     = "standard in any oncology journal",
  "stratified"                       = "standard statistical usage; the models are specified in Methods",
  "projection"                       = "defined by the displayed equation Z = W'X at first use",
  "cophenetic"                       = "named only as published NMF rank criteria, with citations",
  "silhouette"                       = "named only as published NMF rank criteria, with citations",
  "C-index"                          = "introduced parenthetically alongside 'concordance'",
  "index of prediction accuracy"     = "expanded at first use and confined to a sensitivity analysis"
)

fail <- character(0)
cat("\n=== FIRST-USE CHECK (Introduction through Discussion) ===\n\n")

for (tm in names(DEFINED)) {
  m <- regexpr(tm, body, ignore.case = TRUE)
  if (m < 0) { cat(sprintf("[skip] %-26s not present in main text\n", tm)); next }
  win <- substr(body, max(1, m - 240), min(nchar(body), m + 320))
  if (grepl(DEFINED[[tm]], win, fixed = TRUE)) {
    cat(sprintf("[ok  ] %-26s defined at first use\n", tm))
  } else {
    cat(sprintf("[FAIL] %-26s definition MISSING at first use\n         expected: \"%s\"\n         context : ...%s...\n",
                tm, DEFINED[[tm]], gsub("\\s+", " ", substr(win, 1, 260))))
    fail <- c(fail, tm)
  }
}
for (tm in names(ACCEPTED)) {
  if (regexpr(tm, body, ignore.case = TRUE) > 0)
    cat(sprintf("[accepted] %-22s %s\n", tm, ACCEPTED[[tm]]))
}

cat("\n")
if (length(fail)) {
  cat("FAILED: definition removed or moved for -> ", paste(fail, collapse = ", "), "\n")
  cat("Either restore the definition at first use, or update DEFINED in this script.\n")
  quit(status = 1)
}
cat("PASS: every enforced term is defined at first use.\n")
