#!/usr/bin/env Rscript
# 23_stroma_resolution.R
# ---------------------------------------------------------------------------
# Which stromal reference programs does each k = 3 factorization actually
# resolve from bulk?
#
# Motivation. Bulk deconvolution has historically separated ACTIVATED from
# NORMAL stroma; the finer CAF states (myCAF/iCAF, and the proCAF/restCAF
# classification derived from them) were established in single-cell data. This
# script asks, for each method, (a) whether the activated/normal axis is
# recovered, and (b) whether any restraining-stroma (restCAF/iCAF) marker genes
# appear among the genes that define each factor.
#
# TWO VIEWS, deliberately both reported, because they disagree and the honest
# claim depends on which question is being asked:
#   TOP-GENE view  -- are signature genes enriched among a factor's 50 defining
#                     genes? This is what main-text Fig. 2a,b display.
#   FULL-GENE view -- do loadings across all 1,970 genes track membership? Far
#                     more dilute for a 15-25 gene signature, and correspondingly
#                     low powered.
# The full-gene view shows NO method separating restCAF from proCAF. Do not
# report the top-gene view without it.
#
# In-repo only; no new data. Run from the repository root.
# ---------------------------------------------------------------------------

suppressMessages({ library(purrr) })
source("paper/load_precomputed.R")
source("R/get_top_genes.R")

REF <- "data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(REF)) REF <- "../DeSurv-paper/data/derv/cmbSubtypes_formatted.RData"
if (!file.exists(REF)) stop("reference signature file not found: cmbSubtypes_formatted.RData")
load(REF)   # provides `top_genes`

## Assemble reference signatures with the SAME renaming the figure helper uses
## (R/figure_plot_helpers.R). NOTE: SCISSORS iCAF/myCAF are displayed as
## restCAF/proCAF; that relabeling is recorded here so it is not invisible.
mk_refs <- function(g) {
  g$deCAF <- list(
    proCAF  = c("IGFL2","NOX4","VSNL1","BICD1","NPR3","ETV1","ITGA11","CNIH3","COL11A1"),
    restCAF = c("CHRDL1","OGN","PI16","ANK2","ABCA8","TGFBR3","FBLN5","SCARA5","KIAA1217"))
  if (length(g) >= 3)  names(g)[3]  <- "Moffitt"
  if (length(g) >= 4)  names(g)[4]  <- "Moffitt"
  if (length(g) >= 13) names(g)[13] <- "SCISSORS"
  if (length(g) >= 16) names(g)[16] <- "SCISSORS"
  if (length(g) >= 12) names(g)[12] <- "Elyada"
  t <- purrr::list_flatten(g)
  ren <- c("SCISSORS_iCAF" = "SCISSORS_restCAF", "SCISSORS_myCAF" = "SCISSORS_proCAF",
           "SCISSORS_panCAF_vs_peri_top25_panCAF" = "SCISSORS_panCAF")
  i <- names(t) %in% names(ren); names(t)[i] <- ren[names(t)[i]]
  t
}
refs <- mk_refs(top_genes)

ntop <- as.integer(read_result("tar_params_best_tcgacptac")$ntop)
FITS <- list(
  DeSurv    = read_result("tar_fit_desurv_tcgacptac"),
  NMF_k3    = read_result("fit_std_desurvk_tcgacptac"),
  fair_a0k3 = read_result("tar_fit_desurv_a0k3_tcgacptac")
)
KEY <- c("SCISSORS_restCAF","Elyada_iCAF","SCISSORS_proCAF","Elyada_myCAF",
         "DECODER_ActivatedStroma","DECODER_NormalStroma",
         "Moffitt_Activated","Moffitt_Normal")
KEY <- KEY[KEY %in% names(refs)]

## ---- (A) TOP-GENE view: presence and correlation among 50 defining genes ----
topview <- function(fit) {
  tops <- get_top_genes(W = fit$W, ntop = ntop)$top_genes[1:50, ]
  W  <- fit$W[unlist(tops), , drop = FALSE]
  cg <- Reduce(intersect, list(rownames(W), unique(unlist(refs))))
  W  <- W[cg, , drop = FALSE]
  n_present <- vapply(refs[KEY], function(v) sum(cg %in% v), integer(1))
  rho <- vapply(refs[KEY], function(v) {
    m <- as.numeric(cg %in% v)
    if (length(unique(m)) < 2) NA_real_ else
      max(abs(suppressWarnings(apply(W, 2, cor, y = m, method = "spearman"))), na.rm = TRUE)
  }, numeric(1))
  list(n_present = n_present, rho = rho)
}

## ---- (B) FULL-GENE view: loadings vs membership over all analysed genes ----
fullview <- function(fit) {
  g <- rownames(fit$W)
  vapply(refs[KEY], function(v) {
    m <- as.numeric(g %in% v)
    if (sum(m) < 3) NA_real_ else
      max(abs(suppressWarnings(apply(fit$W, 2, cor, y = m, method = "spearman"))), na.rm = TRUE)
  }, numeric(1))
}

top <- lapply(FITS, topview)
ful <- lapply(FITS, fullview)
sig_size <- vapply(refs[KEY], length, integer(1))

res <- list(
  programs   = KEY,
  sig_size   = sig_size,
  top_n      = sapply(top, `[[`, "n_present"),
  top_rho    = sapply(top, `[[`, "rho"),
  full_rho   = do.call(cbind, ful),
  ntop       = ntop,
  note = paste("SCISSORS iCAF/myCAF are displayed as restCAF/proCAF, matching",
               "R/figure_plot_helpers.R. The full-gene view shows no method",
               "separating restCAF from proCAF; only the activated/normal axis",
               "is recovered, and comparably by all methods.")
)
saveRDS(res, "results/stroma_resolution_stats.rds")

cat("\n=== TOP-GENE view: signature genes among each factor's 50 defining genes ===\n")
tn <- res$top_n; colnames(tn) <- names(FITS)
print(cbind(`signature size` = sig_size, tn))
cat("\n=== FULL-GENE view: |Spearman| over all analysed genes ===\n")
print(round(res$full_rho, 3))
cat("\nSaved -> results/stroma_resolution_stats.rds\n")
