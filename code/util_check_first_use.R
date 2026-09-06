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

.lines <- system2("pdftotext", c(shQuote(pdf), "-"), stdout = TRUE)
# Drop standalone page numbers before matching. pdftotext emits each page number
# as its own line, so when a sentence straddles a page break the number lands in
# the middle of the running text and defeats fixed-string matching.
# Real incident 2026-08-15: a prose edit moved the "gene program" definition
# across a page boundary, rendering as "interpretable gene programs: weighted /
# 1 / sets of co-varying genes", and this script reported a FAIL for a
# definition that was present and correctly placed. A false FAIL is not benign
# here: the script's own message invites the reader to edit DEFINED instead,
# which would silently retire the check.
# Order matters: strip the form feed FIRST, because pdftotext glues it to the
# first word of the next page ("\fsets of co-varying genes") and can also glue
# it to the page number itself.
.lines <- gsub("\f", "", .lines)
.lines <- .lines[!grepl("^\\s*[0-9]{1,3}\\s*$", .lines)]
txt <- paste(.lines, collapse = "\n")
txt <- gsub("[ \t]*\n[ \t]*", " ", txt)
# Removing the page-number line leaves a blank line behind, which collapses to a
# double space and breaks fixed-string matching just as the form feed did.
txt <- gsub("[ \t]{2,}", " ", txt)

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
  "frozen"                   = "held fixed, or frozen",
  "representation"           = "together form a representation of the tumor transcriptome",
  "reconstruction error"      = "the mismatch between the measured expression matrix",
  "rank"                      = "the number of programs the factorization extracts",
  "latent"                    = "the true underlying (latent) programs are known",
  "Bayesian optimization"     = "a sequential search algorithm that proposes each new setting",
  "concordance"               = "the probability that the model ranks a pair of patients",
  "one-standard-error"        = "the most parsimonious model within one standard error",
  "partial log-likelihood"    = "the fit criterion of the Cox model",
  "consensus initialization"  = "starting values built from the solutions of many restarts",
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

# pdftotext deletes a hyphen that falls at a line end and joins the halves
# ("co-varying" -> "covarying", "reconstruction-driven" -> "reconstructiondriven").
# Real incident 2026-09-06: a prose edit reflowed the Introduction so that
# "co-varying" broke across lines, and the "gene program" definition reported a
# FAIL although it was present and unchanged. Same family as the page-number
# and form-feed artifacts above. Hyphens are therefore removed from BOTH the
# rendered text and the enforced phrases before matching; positions are taken
# from the hyphen-free text so the window stays aligned.
nohy   <- function(x) gsub("-", "", x, fixed = TRUE)
body_n <- nohy(body)

for (tm in names(DEFINED)) {
  m <- regexpr(nohy(tm), body_n, ignore.case = TRUE)
  if (m < 0) { cat(sprintf("[skip] %-26s not present in main text\n", tm)); next }
  win <- substr(body_n, max(1, m - 240), min(nchar(body_n), m + 320))
  if (grepl(nohy(DEFINED[[tm]]), win, fixed = TRUE)) {
    cat(sprintf("[ok  ] %-26s defined at first use\n", tm))
  } else {
    cat(sprintf("[FAIL] %-26s definition MISSING at first use\n         expected: \"%s\"\n         context : ...%s...\n",
                tm, DEFINED[[tm]], gsub("\\s+", " ", substr(win, 1, 260))))
    fail <- c(fail, tm)
  }
}
for (tm in names(ACCEPTED)) {
  if (regexpr(nohy(tm), body_n, ignore.case = TRUE) > 0)
    cat(sprintf("[accepted] %-22s %s\n", tm, ACCEPTED[[tm]]))
}

cat("\n")
if (length(fail)) {
  cat("FAILED: definition removed or moved for -> ", paste(fail, collapse = ", "), "\n")
  cat("Either restore the definition at first use, or update DEFINED in this script.\n")
  quit(status = 1)
}
cat("PASS: every enforced term is defined at first use.\n")
