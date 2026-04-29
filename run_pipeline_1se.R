#!/usr/bin/env Rscript
# run_pipeline_1se.R — Render the 1-SE-rule variant of the paper
#
# Reuses the cached BO results in results/precomputed/ntop_bo_50_300/ but
# selects the FULL 1-SE row from BO history (not just k) for the final
# DeSurv and alpha=0 model fits, then re-fits, re-validates, regenerates
# figures, and renders the paper + SI to:
#
#   paper/paper_1se.pdf
#   paper/si_appendix_1se.pdf
#
# Outputs land in:
#   results/precomputed/ntop_bo_50_300_1se/
#   paper/figures/ntop_bo_50_300_1se/
#
# The original "best"-rule outputs are not touched.
#
# Usage:
#   Rscript run_pipeline_1se.R                  # default (1 core)
#   Rscript run_pipeline_1se.R --ncores 8       # parallel seed fits

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library(optparse)

option_list <- list(
  make_option("--ncores", type = "integer", default = 1L,
              help = "Number of cores for parallel steps [default: %default]"),
  make_option("--step", type = "integer", default = 1L,
              help = "Start from this step number (1=bootstrap, 2=03b, 3=04, 4=05, 5=08, 6=09)"),
  make_option("--only", action = "store_true", default = FALSE,
              help = "Run only the specified --step")
)

opts <- parse_args(OptionParser(option_list = option_list))

# Pin the BO run + selection rule.
Sys.setenv(
  DESURV_NTOP_LOWER = "50",
  DESURV_NTOP_UPPER = "300",
  DESURV_PARAM_RULE = "1se",
  DESURV_NCORES     = as.character(opts$ncores)
)

steps <- c(
  "code/03a_bootstrap_1se.R",
  "code/03b_select_k.R",
  "code/04_fit_models.R",
  "code/05_external_validation.R",
  "code/08_figures.R",
  "code/09_render_paper.R"
)

for (i in seq_along(steps)) {
  if (i < opts$step) next
  if (opts$only && i != opts$step) next
  if (!file.exists(steps[i])) {
    message(sprintf("  Skipping step %d (file not found: %s)", i, steps[i]))
    next
  }
  cat(sprintf("\n=== Step %d/%d: %s ===\n", i, length(steps), steps[i]))
  t0 <- Sys.time()
  source(steps[i], local = new.env(parent = globalenv()))
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  cat(sprintf("    Completed in %.1f minutes\n", as.numeric(elapsed)))
}

cat("\n=== 1-SE pipeline complete ===\n")
cat("    Results: results/precomputed/ntop_bo_50_300_1se/\n")
cat("    Paper:   paper/paper_1se.pdf\n")
cat("    SI:      paper/si_appendix_1se.pdf\n")
