# What-if follow-up (isolated, 2026-10-02): do the alpha = 0 k = 6 / k = 7 factors resolve CAF states
# (iCAF vs myCAF) or only the activated / normal stroma axis? Same reference lists and the same two views as
# code/23_stroma_resolution.R (whose mk_refs is evaluated verbatim), but per factor and for the own-rank fits.
# Writes only into analysis/whatif-d1-own-rank/. Run from the repo root.
source("R/paths.R"); source("paper/load_precomputed.R"); source("R/get_top_genes.R")
suppressMessages(library(purrr))
L <- readLines("code/23_stroma_resolution.R")
load(cmb_subtypes_path())                      # provides top_genes (DESURV_CMBSUBTYPES)
eval(parse(text = L[grep("^mk_refs <- function", L):grep("^refs <- mk_refs", L)]))   # defines mk_refs and refs
KEY <- c("Elyada_iCAF","SCISSORS_restCAF","Elyada_myCAF","SCISSORS_proCAF","deCAF_restCAF","deCAF_proCAF",
         "Moffitt_Activated","Moffitt_Normal","DECODER_ActivatedStroma","DECODER_NormalStroma")
KEY <- KEY[KEY %in% names(refs)]
ntop <- as.integer(read_result("tar_params_best_tcgacptac")$ntop)
FITS <- list(DeSurv_k3 = read_result("tar_fit_desurv_tcgacptac"), NMF_k3 = read_result("fit_std_desurvk_tcgacptac"),
             a0_k6 = readRDS("analysis/whatif-d1-own-rank/fit_alpha0_k6.rds"),
             a0_k7 = read_result("tar_fit_desurv_alpha0_tcgacptac"))
rows <- list()
for (fn in names(FITS)) {
  W <- FITS[[fn]]$W; g <- rownames(W)
  tops <- get_top_genes(W = W, ntop = ntop)$top_genes[1:50, , drop = FALSE]
  for (j in seq_len(ncol(W))) for (s in KEY) {
    m <- as.numeric(g %in% refs[[s]])
    rows[[length(rows) + 1]] <- data.frame(fit = fn, factor = j, signature = s, sig_in_universe = sum(m),
      top50_hits = sum(tops[[j]] %in% refs[[s]]),
      full_rho = if (sum(m) < 3) NA else suppressWarnings(cor(W[, j], m, method = "spearman")))
  }
}
res <- do.call(rbind, rows)
write.csv(res, "analysis/whatif-d1-own-rank/stroma_per_factor.csv", row.names = FALSE)
# compact views: top-50 hits (factor x signature) and full-gene rho, per fit
for (fn in names(FITS)) {
  r <- res[res$fit == fn, ]
  cat("\n==", fn, "| top-50 hits (rows = factors) ==\n")
  print(xtabs(top50_hits ~ factor + signature, r)[, KEY])
  cat("-- full-gene Spearman --\n")
  print(round(xtabs(full_rho ~ factor + signature, r)[, KEY], 2))
}
cat("\nsignature sizes in the 1,970-gene universe:\n"); print(tapply(res$sig_in_universe, res$signature, max)[KEY])
