#!/usr/bin/env Rscript
# code/09_render_paper.R — Compile manuscript and supplement
#
# Renders the paper from pre-computed results (no pipeline needed).
# Requires: results/precomputed/*.rds, results/cv_grid/*.rds, figures/*.pdf
#
# Runtime: ~2 min

message("=== Step 9: Render Paper ===")

# Verify pre-computed results exist
precomputed_dir <- "results/precomputed"
n_rds <- length(list.files(precomputed_dir, pattern = "\\.rds$"))
if (n_rds < 25) {
  stop("Only ", n_rds, " pre-computed RDS files found in ", precomputed_dir,
       "\nExpected ~30. Run earlier steps or download pre-computed results.")
}
message("  Found ", n_rds, " pre-computed results")

# Render main paper
message("  Rendering paper/paper.Rmd ...")
rmarkdown::render("paper/paper.Rmd", knit_root_dir = getwd(), quiet = TRUE)
message("  Output: paper/paper.pdf")

# Render supplement
message("  Rendering paper/supplement.Rmd ...")
rmarkdown::render("paper/supplement.Rmd", knit_root_dir = getwd(), quiet = TRUE)
message("  Output: paper/supplement.pdf")

message("=== Step 9 complete ===")
