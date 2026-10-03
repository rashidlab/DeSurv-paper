# run_may_validation.R — mode 1 (May): validate the reconstruction against the cached
# results. Read-only on results/; writes only to the output directory given as argv[1].
#
# Usage: Rscript analysis/published-procedure-reconstruction/run_may_validation.R <out_dir>
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "/tmp/recon_may_out"

this_dir <- "analysis/published-procedure-reconstruction"
source(file.path(this_dir, "reconstruct.R"))

result <- run_may_validation("results", out_dir)

cat("\n=== May validation: reconstruction output ===\n")
cat("k_selected:", result$selection$k_selected, "\n")
cat("k_best:", result$selection$k_best, "\n")
cat("best_mean:", result$selection$best_mean, "\n")
cat("best_se:", result$selection$best_se, "\n")
cat("lcb_threshold:", result$selection$lcb_threshold, "\n")
cat("reason:", result$selection$reason, "\n")
cat("\npublished (global-best, k replaced) params:\n")
print(result$published_global_best_params)
cat("\nconditional-k params:\n")
print(result$conditional_k_params)
cat("\nper-k table:\n")
print(result$per_k_table[, c("eval_id", "k_numeric", "mean_cindex", "c_se", "clears_threshold")])

cat("\n=== Compare against cached results/tar_k_selection_tcgacptac.rds ===\n")
cached_k_sel <- readRDS("results/tar_k_selection_tcgacptac.rds")
print(cached_k_sel)
cat("identical(selection, cached):", identical(result$selection, cached_k_sel), "\n")
cat("all.equal(selection, cached):\n")
print(all.equal(result$selection, cached_k_sel))

cat("\n=== Compare against cached results/tar_params_best_tcgacptac.rds ===\n")
cached_params <- readRDS("results/tar_params_best_tcgacptac.rds")
print(cached_params)
cat("identical(published_global_best_params, cached):",
    identical(result$published_global_best_params, cached_params), "\n")
cat("all.equal(published_global_best_params, cached):\n")
print(all.equal(result$published_global_best_params, cached_params))
