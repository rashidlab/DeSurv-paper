#!/usr/bin/env Rscript
# run_pipeline.R — Run the DeSurv analysis pipeline
#
# Usage:
#   Rscript run_pipeline.R                       # Render paper only (ntop BO 50-300 by default)
#   Rscript run_pipeline.R --step 9              # Regenerate figures + render
#   Rscript run_pipeline.R --ntop-lower 50 --ntop-upper 300  # Explicit ntop bounds
#   Rscript run_pipeline.R --quick                                 # Quick smoke test (~10 min)
#   Rscript run_pipeline.R --full --ncores 8                       # Full re-computation
#   Rscript run_pipeline.R --step 8                                # Run from step 8 onward
#   Rscript run_pipeline.R --step 8 --only                         # Run only step 8

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}
library(optparse)

option_list <- list(
  make_option("--quick", action = "store_true", default = FALSE,
              help = "Quick mode: reduced data/iterations for smoke testing"),
  make_option("--full", action = "store_true", default = FALSE,
              help = "Full re-computation from raw data"),
  make_option("--step", type = "integer", default = 1L,
              help = "Start from this step number [default: %default]"),
  make_option("--only", action = "store_true", default = FALSE,
              help = "Run only the specified step"),
  make_option("--ncores", type = "integer", default = 1L,
              help = "Number of cores for parallel steps [default: %default]"),
  make_option("--ntop-lower", type = "integer", default = 50L,
              help = "Lower bound for BO-tuned ntop (sets DESURV_NTOP_LOWER) [default: %default]"),
  make_option("--ntop-upper", type = "integer", default = 300L,
              help = "Upper bound for BO-tuned ntop (sets DESURV_NTOP_UPPER) [default: %default]"),
  make_option("--ntop", type = "integer", default = NULL,
              help = "Fixed ntop value; overrides --ntop-lower/--ntop-upper (sets DESURV_NTOP)")
)

opts <- parse_args(OptionParser(option_list = option_list))

# Set environment variables consumed by individual scripts
if (opts$quick) {
  Sys.setenv(DESURV_QUICK = "TRUE", DESURV_RECOMPUTE = "TRUE", DESURV_NCORES = "1")
} else if (opts$full) {
  Sys.setenv(DESURV_RECOMPUTE = "TRUE", DESURV_NCORES = as.character(opts$ncores))
}

if (!is.null(opts[["ntop"]])) {
  Sys.setenv(DESURV_NTOP = as.character(opts[["ntop"]]))
} else {
  Sys.setenv(DESURV_NTOP_LOWER = as.character(opts[["ntop-lower"]]),
             DESURV_NTOP_UPPER = as.character(opts[["ntop-upper"]]))
}

steps <- c(
  "code/01_install.R",
  "code/02_load_data.R",
  "code/03_bayesian_optimization.R",
  "code/04_fit_models.R",
  "code/05_external_validation.R",
  "code/06_sensitivity_analysis.R",
  "code/07_simulations.R",
  "code/08a_cutpoint_analysis.R",
  "code/08_figures.R",
  "code/08b_si_figures.R",
  "code/08c_sim_figures.R",
  "code/09_render_paper.R"
)

# Default: just render paper from precomputed results
if (!opts$quick && !opts$full && opts$step == 1L) {
  opts$step <- length(steps)
  message("No --quick or --full specified. Running render-paper step only.")
  message("Use --quick for smoke test or --full for re-computation.\n")
}

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

cat("\n=== Pipeline complete ===\n")
