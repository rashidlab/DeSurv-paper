# R/gene_lists.R: builder for paper/gene_lists_top<ntop>_k<k>.csv (SI Supplementary Data).
#
# One row per (factor, rank): rank 1..ntop in desurv_get_top_genes() order, gene,
# weight = W[gene, factor], contrast_score = the specificity score s_ij
# (desurv_get_top_genes()$diffs), factor label. DeSurv factors are D1..Dk
# (tar_fit_desurv_tcgacptac), standard NMF at DeSurv's k are N1..Nk
# (fit_std_desurvk_tcgacptac). Numbers are rounded to 6 decimals, as in the
# committed CSV. Pure function: no I/O.
build_gene_lists <- function(W_desurv, W_std, ntop, digits = 6) {
  ntop <- as.integer(round(ntop))
  one <- function(W, prefix) {
    tg <- DeSurv::desurv_get_top_genes(W, ntop)
    do.call(rbind, lapply(seq_len(ncol(W)), function(j) {
      g <- as.character(tg$top_genes[[j]])
      data.frame(rank = seq_along(g), gene = g,
                 weight = round(unname(W[g, j]), digits),
                 contrast_score = round(unname(tg$diffs[g, j]), digits),
                 factor = paste0(prefix, j), stringsAsFactors = FALSE)
    }))
  }
  out <- rbind(one(W_desurv, "D"), one(W_std, "N"))
  rownames(out) <- NULL
  out
}
