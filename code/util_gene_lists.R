#!/usr/bin/env Rscript
# code/util_gene_lists.R: produce paper/gene_lists_top<ntop>_k<k>.csv, the Supplementary
# Data gene tables read by paper/si_appendix.Rmd (chunk tab-gene-list-270).
#
# Reads (never recomputes; no cache_or_compute):
#   tar_fit_desurv_tcgacptac   DeSurv fit       -> factors D1..Dk
#   fit_std_desurvk_tcgacptac  std NMF at k     -> factors N1..Nk
#   tar_params_best_tcgacptac  ntop (rounded) and k
# from the results base dir (R/paths.R: DESURV_RESULTS_DIR, default results/).
#
# Output dir: DESURV_GENE_LIST_DIR (default "paper"). The file name encodes ntop and k,
# so a rerun that moves either writes a NEW file and the SI chunk (which hardcodes
# gene_lists_top270_k3.csv) must be updated deliberately.
#
#   Run as:  Rscript code/util_gene_lists.R
source("R/paths.R")
source("R/gene_lists.R")
suppressMessages(library(DeSurv))

fit_d  <- readRDS(res_path("tar_fit_desurv_tcgacptac.rds"))
fit_s  <- readRDS(res_path("fit_std_desurvk_tcgacptac.rds"))
params <- readRDS(res_path("tar_params_best_tcgacptac.rds"))

ntop <- as.integer(round(params$ntop))
k    <- as.integer(round(params$k))
stopifnot("DeSurv W has k columns" = ncol(fit_d$W) == k,
          "std NMF W has k columns" = ncol(fit_s$W) == k)

out <- build_gene_lists(fit_d$W, fit_s$W, ntop)

out_dir <- Sys.getenv("DESURV_GENE_LIST_DIR", "paper")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, sprintf("gene_lists_top%d_k%d.csv", ntop, k))
write.csv(out, out_path, row.names = FALSE)
message(sprintf("Wrote %s (%d rows; ntop = %d, k = %d, DeSurv %s)",
                out_path, nrow(out), ntop, k, as.character(packageVersion("DeSurv"))))
