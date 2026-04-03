#!/usr/bin/env Rscript
# code/09_render_paper.R — Compile manuscript and supplement
#
# Renders the paper from pre-computed results (no pipeline needed).
# Requires: results/precomputed/*.rds, results/cv_grid/*.rds, figures/*.pdf
#
# Runtime: ~2 min

message("=== Step 9: Render Paper ===")
source("code/00_helpers.R")

# Verify pre-computed results exist
# In a subfolder (ntop-specific run), fewer files are expected (~15-20 vs ~30+)
min_rds <- if (nzchar(CONFIG$ntop_subfolder)) 10 else 25
n_rds <- length(list.files(RESULTS_DIR, pattern = "\\.rds$"))
if (n_rds < min_rds) {
  stop("Only ", n_rds, " pre-computed RDS files found in ", RESULTS_DIR,
       "\nExpected >= ", min_rds, ". Run earlier steps or download pre-computed results.")
}
message("  Found ", n_rds, " pre-computed results in ", RESULTS_DIR)

# Ensure figure output directory exists (knitr fig.path is relative to paper/)
fig_dir <- if (nzchar(CONFIG$ntop_subfolder)) {
  file.path("paper/figures", CONFIG$ntop_subfolder)
} else {
  "paper/figures"
}
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Render main paper
message("  Rendering paper/paper.Rmd ...")
rmarkdown::render("paper/paper.Rmd", knit_root_dir = getwd(), quiet = TRUE)
message("  Output: paper/paper.pdf")

# Render supplement
message("  Rendering paper/si_appendix.Rmd ...")
rmarkdown::render("paper/si_appendix.Rmd", knit_root_dir = getwd(), quiet = TRUE)
message("  Output: paper/si_appendix.pdf")

message("=== Step 9 complete ===")
