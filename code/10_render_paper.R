#!/usr/bin/env Rscript
# code/10_render_paper.R — Compile manuscript and supplement
#
# Renders the paper from pre-computed results (no pipeline needed).
# Requires: results/*.rds, results/cv_grid/*.rds, figures/*.pdf
#
# Runtime: ~5 min (3-pass build, see below)
#
# CROSS-DOCUMENT REFERENCES (xr-hyper): the main paper and the SI reference each
# other's tables (paper -> SI Table \ref{tab:projection-r2}; SI -> main paper
# \ref{tab:hr-adjusted}). Because they compile as two separate PDFs, \ref can
# only resolve a label in the *other* document if that document's .aux file
# already exists. We therefore (a) keep the .aux files (clean = FALSE,
# tinytex.clean = FALSE) and (b) render in the order paper -> SI -> paper so each
# pass picks up the other document's .aux from the previous pass:
#   pass 1 (paper): SI .aux absent  -> tab:projection-r2 unresolved (expected)
#   pass 2 (SI)   : paper.aux ready -> tab:hr-adjusted resolves
#   pass 3 (paper): SI .aux ready   -> tab:projection-r2 resolves
# After the final pass we delete the LaTeX auxiliaries (keeping .pdf and .tex).

message("=== Step 10: Render Paper ===")
source("code/00_helpers.R")

# Verify pre-computed results exist
n_rds <- length(list.files(RESULTS_DIR, pattern = "\\.rds$"))
if (n_rds < 25) {
  stop("Only ", n_rds, " pre-computed RDS files found in ", RESULTS_DIR,
       "\nExpected >= 25. Run earlier steps or download pre-computed results.")
}
message("  Found ", n_rds, " pre-computed results in ", RESULTS_DIR)

# Keep LaTeX auxiliaries between passes so xr-hyper can read the other
# document's labels (see header note).
old_opts <- options(tinytex.clean = FALSE)
on.exit(options(old_opts), add = TRUE)

render_pass <- function(rmd, out) {
  # envir = globalenv() preserves the original top-level render semantics: the
  # two documents knit into a shared environment, and some SI chunks rely on
  # objects loaded during an earlier render (the SI is not fully self-contained).
  rmarkdown::render(rmd, knit_root_dir = getwd(), quiet = TRUE,
                    output_file = out, clean = FALSE, envir = globalenv())
}

message("  Pass 1/3: paper.Rmd (seeds paper.aux) ...")
render_pass("paper/paper.Rmd", "paper.pdf")
message("  Pass 2/3: si_appendix.Rmd (resolves SI -> main refs) ...")
render_pass("paper/si_appendix.Rmd", "si_appendix.pdf")
message("  Pass 3/3: paper.Rmd (resolves main -> SI refs) ...")
render_pass("paper/paper.Rmd", "paper.pdf")
message("  Output: paper/paper.pdf, paper/si_appendix.pdf")

# Tidy LaTeX auxiliaries left behind by clean = FALSE (keep .pdf and .tex).
aux <- list.files("paper", full.names = TRUE,
  pattern = "\\.(aux|log|out|fls|fdb_latexmk|toc|lof|lot|nav|snm|bbl|blg|synctex\\.gz)$")
if (length(aux)) file.remove(aux)

message("=== Step 10 complete ===")
