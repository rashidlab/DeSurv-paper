#!/usr/bin/env Rscript
# code/10_render_paper.R — Compile manuscript and supplement
#
# Renders the paper from pre-computed results (no pipeline needed).
# Requires: results/*.rds, results/cv_grid/*.rds, figures/*.pdf
#
# Runtime: ~2 min

message("=== Step 10: Render Paper ===")
source("code/00_helpers.R")

# Verify pre-computed results exist
n_rds <- length(list.files(RESULTS_DIR, pattern = "\\.rds$"))
if (n_rds < 25) {
  stop("Only ", n_rds, " pre-computed RDS files found in ", RESULTS_DIR,
       "\nExpected >= 25. Run earlier steps or download pre-computed results.")
}
message("  Found ", n_rds, " pre-computed results in ", RESULTS_DIR)

# Render main paper
message("  Rendering paper/paper.Rmd ...")
rmarkdown::render("paper/paper.Rmd", knit_root_dir = getwd(), quiet = TRUE,
                  output_file = "paper.pdf")
message("  Output: paper/paper.pdf")

# Render supplement
message("  Rendering paper/si_appendix.Rmd ...")
rmarkdown::render("paper/si_appendix.Rmd", knit_root_dir = getwd(), quiet = TRUE,
                  output_file = "si_appendix.pdf")
message("  Output: paper/si_appendix.pdf")

message("=== Step 10 complete ===")
