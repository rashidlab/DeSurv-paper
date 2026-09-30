# Paths for the pipeline tests. testthat runs each file with the working
# directory set to tests/, so repo files are one level up. These tests source
# only pure-definition files (R/bo_helpers.R, R/consensus_settings.R) and never
# read or write results/.
repo_file <- function(...) normalizePath(file.path("..", ...), mustWork = TRUE)
