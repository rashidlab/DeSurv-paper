# What-if (isolated, not part of the pipeline or the paper, 2026-10-02):
# is the published D1 present among the factors of an unsupervised (alpha = 0) fit at its own BO-selected rank?
#   - k = 6: NEW alpha = 0 BO (step4-dd54827-r2, params below), fitted here with the reported recipe
#            (DeSurv 1.2.0 desurv_reported_fit: 100 seed fits -> consensus seed -> final fit).
#   - k = 7: the paper's existing rank-optimized control (results/tar_fit_desurv_alpha0_tcgacptac.rds).
#   - k = 3 baselines: the paper's matched-rank comparators (fit_std_desurvk, tar_fit_desurv_a0k3).
# D1-D3 = the published DeSurv fit (results/tar_fit_desurv_tcgacptac.rds). All descriptive; no thresholds.
# Run from the repo root: Rscript analysis/whatif-d1-own-rank/run.R
suppressMessages(library(DeSurv)); stopifnot(packageVersion("DeSurv") >= "1.2.0")
source("R/reported_settings.R"); source("R/reconstruction_helpers.R"); source("R/get_top_genes.R")
out <- "analysis/whatif-d1-own-rank"
rd  <- function(n) readRDS(file.path("results", paste0(n, ".rds")))
td  <- rd("tar_data_filtered_tcgacptac")   # identical content to the rerun's (md5 177be7ac on Longleaf)
S   <- desurv_reported_settings(FALSE)

f6_path <- file.path(out, "fit_alpha0_k6.rds")
if (!file.exists(f6_path)) {
  p6 <- list(k = 6L, alpha = 0, lambda = 0.05634645, nu = 0.90931298, ntop = 195L)  # step4-dd54827-r2 alpha0
  t0 <- Sys.time()
  f6 <- desurv_reported_fit(X = td$ex, y = td$sampInfo$time, d = td$sampInfo$event,
          k = p6$k, alpha = 0, lambda = p6$lambda, nu = p6$nu, ntop = p6$ntop,
          seeds = seq_len(S$n_seed_fits), tol = S$tol, tol_init = S$tol_init,
          maxit = S$maxit, imaxit = S$imaxit, ncores = as.integer(Sys.getenv("NCORES", "8")))
  attr(f6, "whatif") <- list(params = p6, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                             desurv = as.character(packageVersion("DeSurv")))
  saveRDS(f6, f6_path)
}
f6 <- readRDS(f6_path)
fits <- list(`alpha0 k=6 (new BO)` = f6, `alpha0 k=7 (paper)` = rd("tar_fit_desurv_alpha0_tcgacptac"),
             `std NMF k=3 (paper Fig 2)` = rd("fit_std_desurvk_tcgacptac"),
             `alpha0 k=3 matched budget (paper)` = rd("tar_fit_desurv_a0k3_tcgacptac"))
d  <- rd("tar_fit_desurv_tcgacptac"); Wd <- d$W; Hd <- d$H
ntop <- as.integer(rd("tar_params_best_tcgacptac")$ntop)
topd <- get_top_genes(Wd, ntop)$top_genes
rows <- list()
for (nm in names(fits)) {
  Wb <- fits[[nm]]$W; Hb <- fits[[nm]]$H
  stopifnot(identical(rownames(Wb), rownames(Wd)))
  r2g <- projection_r2(Wd, Wb)                                  # paper's statistic (loadings on basis, OLS R^2)
  rho <- cor(Wd, Wb, method = "spearman")                         # paper's .corrC
  r2s <- vapply(1:3, function(j) summary(lm(Hd[j, ] ~ t(Hb)))$r.squared, 0)   # sample scores
  topb <- get_top_genes(Wb, ntop)$top_genes
  jac <- sapply(1:3, function(j) max(sapply(topb, function(g) length(intersect(topd[[j]], g)) / length(union(topd[[j]], g)))))
  for (j in 1:3) rows[[length(rows) + 1]] <- data.frame(basis = nm, k_basis = ncol(Wb), factor = paste0("D", j),
    gene_R2 = r2g[j], max_abs_spearman = max(abs(rho[j, ])), best_match = which.max(abs(rho[j, ])),
    score_R2 = r2s[j], best_top_jaccard = jac[j])
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out, "d1_own_rank_results.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
cat("k=6 fit minutes:", round(attr(f6, "whatif")$minutes, 1), " DeSurv", attr(f6, "whatif")$desurv, "\n")
