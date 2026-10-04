# Usage: Rscript f2_val.R <reopt.rds> [ntop]
# Score archived vs re-optimized fit on the cached external cohorts with the production
# rule of code/05 compute_val_cindex(): union of top-ntop genes per factor, frozen W, beta.
suppressMessages(library(DeSurv))
args <- commandArgs(TRUE); ntop <- if (length(args) >= 2) as.integer(args[2]) else 270L
f0 <- readRDS(if (length(args) >= 3) args[3] else "results/tar_fit_desurv_tcgacptac.rds"); r <- readRDS(args[1])
rownames(r$W) <- rownames(f0$W)
dv_list <- readRDS("results/data_val_filtered_tcgacptac.rds")
score <- function(W, beta, dv) {
  g <- unique(unlist(DeSurv::desurv_get_top_genes(W, ntop)$top_genes))
  lp <- as.numeric(t(dv$ex[g, , drop = FALSE]) %*% W[g, , drop = FALSE] %*% beta)
  keep <- dv$sampInfo$keep == 1
  survival::concordance(survival::Surv(dv$sampInfo$time[keep], dv$sampInfo$event[keep]) ~ lp[keep],
                        reverse = TRUE)$concordance
}
out <- t(sapply(dv_list, function(dv) c(archived = score(f0$W, as.vector(f0$beta), dv),
                                         reoptimized = score(r$W, r$beta, dv))))
rownames(out) <- sapply(dv_list, `[[`, "dataname")
n <- sapply(dv_list, function(dv) sum(dv$sampInfo$keep == 1))
print(round(out, 4)); cat("n-weighted mean:", round(colSums(out * n) / sum(n), 4), "\n")
