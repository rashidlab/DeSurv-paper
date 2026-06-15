#!/usr/bin/env Rscript
# 13_spatial_cooccurrence.R
# ---------------------------------------------------------------------------
# SPATIAL validation of the D1 tumour-stroma coupling (Nature Cancer).
#
# The coupling is a tissue-architecture phenomenon, validated at spot resolution.
# SPATIAL ADJACENCY: do Classical-malignant regions sit next to restraining-CAF
# (restCAF) stroma, while Basal-like regions sit next to promoting-CAF (proCAF)
# stroma? Spots (not single cells) are the unit -- multicellular, far denser than
# scRNA, and the key adjacency test uses marker/DeCAF signature calls -- so it is
# well powered (tens of thousands of spots across 7 sections) and dropout-robust.
#
# Data: UNC Visium CytAssist spatial transcriptomics, GSE311783 (now public via
# Peng et al. Cell Rep Med 2026). Read WITHOUT Seurat -- 10x triplets + Visium
# tissue_positions are plain Matrix / data.frame. Frozen DeSurv W projected per
# spot with the manuscript convention (rank over W namespace -> top-270 W-tilde).
# ---------------------------------------------------------------------------

suppressMessages({ library(Matrix) })
set.seed(1)
ROOT  <- "/home/naimrashid/Downloads/DeSurv-paper-clean"
STDIR <- "/tmp/unc_st"
setwd(ROOT)

## --- trained fit + top-270 basis (identical convention to code/11) -------------
fit <- readRDS(file.path("results","tar_fit_desurv_tcgacptac.rds"))
W <- fit$W; Wg <- rownames(W)
source("R/get_top_genes.R")
PGENES <- unique(unlist(lapply(get_top_genes(fit$W, 270)$top_genes, as.character)))

## --- reference signatures -----------------------------------------------------
proCAF_sig  <- c("IGFL2","NOX4","VSNL1","BICD1","NPR3","ETV1","ITGA11","CNIH3","COL11A1")
restCAF_sig <- c("CHRDL1","OGN","PI16","ANK2","ABCA8","TGFBR3","FBLN5","SCARA5","KIAA1217")
epi_sig     <- c("EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1")
fib_sig     <- c("COL1A1","COL1A2","LUM","DCN","PDGFRB","PDGFRA","FN1","THY1")
classical_sig <- c("GATA6","TFF1","TFF2","TFF3","LGALS4","CLDN18","AGR2","CEACAM6","REG4","FOXA2")
basal_sig     <- c("KRT5","KRT6A","KRT14","KRT17","S100A2","TP63","SPRR3","SPRR1B","FAM83A","LY6D")

## --- per-spot helpers ---------------------------------------------------------
# log1p-CPM normalisation (per spot), then signature = mean of z-scored genes
lognorm <- function(M){            # M: genes x spots counts
  cs <- Matrix::colSums(M); cs[cs == 0] <- 1
  M <- M %*% Diagonal(x = 1e4 / cs); M@x <- log1p(M@x); M
}
sigZ <- function(L, sig){          # L: log-norm genes x spots ; z across spots
  g <- intersect(sig, rownames(L)); if (length(g) < 3) return(rep(NA_real_, ncol(L)))
  s <- Matrix::colMeans(L[g, , drop = FALSE]); as.numeric(scale(s)) }
rankProj <- function(counts, j){   # rank-then-project (manuscript convention), per spot
  sh <- intersect(Wg, rownames(counts))
  xr <- apply(as.matrix(counts[sh, , drop = FALSE]), 2, rank, ties.method = "average")
  rownames(xr) <- sh; g <- intersect(PGENES, sh)
  as.numeric(drop(t(xr[g, , drop = FALSE]) %*% W[g, j, drop = FALSE])) }

samples <- sub("_matrix.mtx.gz$","", list.files(STDIR, pattern="matrix.mtx.gz$"))
cat("== spatial samples:", length(samples), "==\n")

all_spots <- list()
for (s in samples) {
  pf <- file.path(STDIR, s)
  M  <- as(Matrix::readMM(gzfile(paste0(pf,"_matrix.mtx.gz"))), "CsparseMatrix")
  bc <- readLines(gzfile(paste0(pf,"_barcodes.tsv.gz")))
  ft <- read.delim(gzfile(paste0(pf,"_features.tsv.gz")), header=FALSE)
  rownames(M) <- make.unique(ft$V2); colnames(M) <- bc
  pos <- read.csv(gzfile(paste0(pf,"_tissue_positions.csv.gz")))
  pos <- pos[pos$in_tissue == 1, ]
  keep <- intersect(colnames(M), pos$barcode)
  M <- M[, keep, drop=FALSE]; pos <- pos[match(keep, pos$barcode), ]
  if (ncol(M) < 50) next
  L <- lognorm(M)
  df <- data.frame(
    sample = sub("^GSM[0-9]+_","", s), barcode = keep,
    row = pos$array_row, col = pos$array_col,
    D1 = scale(rankProj(M,1))[,1], D2 = scale(rankProj(M,2))[,1], D3 = scale(rankProj(M,3))[,1],
    epi = sigZ(L, epi_sig), fib = sigZ(L, fib_sig),
    classical = sigZ(L, classical_sig), basal = sigZ(L, basal_sig),
    proCAF = sigZ(L, proCAF_sig), restCAF = sigZ(L, restCAF_sig),
    stringsAsFactors = FALSE)
  all_spots[[s]] <- df
  cat(sprintf("  %-12s %5d spots\n", df$sample[1], nrow(df)))
}
spots <- do.call(rbind, all_spots); rownames(spots) <- NULL
cat("total in-tissue spots:", nrow(spots), "\n\n")

## ===========================================================================
## Compartment assignment (per spot): epithelial- vs fibroblast-dominant
## ===========================================================================
spots$compartment <- ifelse(spots$epi > spots$fib, "epithelial", "fibroblast")
# epithelial spots: Classical vs Basal ; fibroblast spots: restCAF vs proCAF
spots$epi_type <- NA_character_
ei <- spots$compartment == "epithelial"
spots$epi_type[ei] <- ifelse(spots$classical[ei] > spots$basal[ei], "Classical", "Basal")
spots$caf_type <- NA_character_
fi <- spots$compartment == "fibroblast"
spots$caf_type[fi] <- ifelse(spots$restCAF[fi] > spots$proCAF[fi], "restCAF", "proCAF")

## ===========================================================================
## (1) Spot-level co-occurrence of programs (within-sample Spearman, combined)
## ===========================================================================
wcor <- function(x, y, grp){           # mean within-sample Spearman (Fisher-z combined)
  zs <- c(); ns <- c()
  for (g in unique(grp)) { k <- grp==g & is.finite(x) & is.finite(y)
    if (sum(k) > 30) { r <- suppressWarnings(cor(x[k], y[k], method="spearman"))
      if (is.finite(r)) { zs <- c(zs, atanh(pmin(pmax(r,-.999),.999))); ns <- c(ns, sum(k)-3) } } }
  if (!length(zs)) return(c(rho=NA, p=NA, k=0))
  zbar <- sum(zs*ns)/sum(ns); se <- sqrt(1/sum(ns))
  c(rho = tanh(zbar), p = 2*pnorm(-abs(zbar/se)), k = length(zs)) }

co_class_rest <- wcor(spots$classical, spots$restCAF, spots$sample)  # coupling: Classical~restCAF
co_class_pro  <- wcor(spots$classical, spots$proCAF,  spots$sample)
co_basal_pro  <- wcor(spots$basal,     spots$proCAF,  spots$sample)
co_basal_rest <- wcor(spots$basal,     spots$restCAF, spots$sample)
co_D1_rest    <- wcor(spots$D1,        spots$restCAF, spots$sample)
co_D1_pro     <- wcor(spots$D1,        spots$proCAF,  spots$sample)

## ===========================================================================
## (2) SPATIAL ADJACENCY: restCAF enriched next to Classical (vs Basal)?
##     Visium hex grid: neighbours are array (row +/-1, col +/-1) and (row, col +/-2).
##     For each epithelial spot, mean restCAF/proCAF over its fibroblast neighbours.
## ===========================================================================
neigh_off <- list(c(-1,-1),c(-1,1),c(1,-1),c(1,1),c(0,-2),c(0,2))
adj_rows <- list()
for (g in unique(spots$sample)) {
  sp <- spots[spots$sample==g, ]
  key <- paste(sp$row, sp$col, sep="_"); idx <- setNames(seq_len(nrow(sp)), key)
  for (i in which(sp$compartment=="epithelial")) {
    nb <- c()
    for (o in neigh_off) { kk <- paste(sp$row[i]+o[1], sp$col[i]+o[2], sep="_")
      if (!is.na(idx[kk])) nb <- c(nb, idx[kk]) }
    nb <- nb[!is.na(nb)]; nbf <- nb[sp$compartment[nb]=="fibroblast"]
    if (length(nbf) >= 1)
      adj_rows[[length(adj_rows)+1]] <- data.frame(
        sample=g, epi_type=sp$epi_type[i],
        nb_restCAF=mean(sp$restCAF[nbf], na.rm=TRUE),
        nb_proCAF =mean(sp$proCAF[nbf],  na.rm=TRUE),
        n_fib_nb=length(nbf)) }
}
adj <- do.call(rbind, adj_rows)
# restCAF-vs-proCAF neighbour signal, Classical vs Basal epithelial spots
adj$nb_diff <- adj$nb_restCAF - adj$nb_proCAF      # >0 => neighbourhood leans restCAF
class_diff <- adj$nb_diff[adj$epi_type=="Classical"]
basal_diff <- adj$nb_diff[adj$epi_type=="Basal"]
adj_test <- wilcox.test(class_diff, basal_diff)
# within-sample permutation null: shuffle epi_type within each sample, recompute mean gap
obs_gap <- mean(class_diff, na.rm=TRUE) - mean(basal_diff, na.rm=TRUE)
nperm <- 2000; perm_gaps <- numeric(nperm)
for (p in 1:nperm) {
  lab <- ave(adj$epi_type, adj$sample, FUN=function(z) sample(z))
  perm_gaps[p] <- mean(adj$nb_diff[lab=="Classical"], na.rm=TRUE) -
                  mean(adj$nb_diff[lab=="Basal"], na.rm=TRUE) }
perm_p <- (1 + sum(abs(perm_gaps) >= abs(obs_gap))) / (nperm + 1)

## ===========================================================================
## Save + report
## ===========================================================================
res <- list(
  meta = list(n_samples=length(unique(spots$sample)), n_spots=nrow(spots),
              n_epi=sum(ei), n_fib=sum(fi),
              n_classical=sum(spots$epi_type=="Classical",na.rm=TRUE),
              n_basal=sum(spots$epi_type=="Basal",na.rm=TRUE),
              n_restCAF=sum(spots$caf_type=="restCAF",na.rm=TRUE),
              n_proCAF=sum(spots$caf_type=="proCAF",na.rm=TRUE), ntop=270),
  spotlevel = list(Classical_restCAF=co_class_rest, Classical_proCAF=co_class_pro,
                   Basal_proCAF=co_basal_pro, Basal_restCAF=co_basal_rest,
                   D1_restCAF=co_D1_rest, D1_proCAF=co_D1_pro),
  adjacency = list(class_nb_diff_mean=mean(class_diff,na.rm=TRUE),
                   basal_nb_diff_mean=mean(basal_diff,na.rm=TRUE),
                   gap=obs_gap, wilcox_p=adj_test$p.value, perm_p=perm_p,
                   n_classical_epi=length(class_diff), n_basal_epi=length(basal_diff)),
  spots_summary = aggregate(cbind(D1,D2,D3,classical,basal,restCAF,proCAF)~sample, spots, mean))
saveRDS(res, file.path("results","spatial_cooccurrence_stats.rds"))
saveRDS(spots[,c("sample","row","col","D1","D2","D3","compartment","epi_type","caf_type",
                 "classical","basal","restCAF","proCAF")],
        file.path("results","spatial_spots_scored.rds"))
saveRDS(adj, file.path("results","spatial_adjacency_perspot.rds"))

fmtp <- function(p) ifelse(is.na(p),"NA", ifelse(p<1e-4,"< 0.0001", sprintf("%.4f",p)))
cat("\n================= SPATIAL CO-OCCURRENCE (UNC Visium, GSE311783) =================\n")
cat(sprintf("samples: %d | in-tissue spots: %d (epi %d, fib %d)\n",
            res$meta$n_samples, res$meta$n_spots, res$meta$n_epi, res$meta$n_fib))
cat(sprintf("epithelial: Classical %d / Basal %d | fibroblast: restCAF %d / proCAF %d\n\n",
            res$meta$n_classical, res$meta$n_basal, res$meta$n_restCAF, res$meta$n_proCAF))
cat("(1) Spot-level co-occurrence (within-sample Spearman, Fisher-combined):\n")
pr <- function(nm,v) cat(sprintf("    %-20s rho = % .3f  p %s  (%d sections)\n", nm, v["rho"], fmtp(v["p"]), v["k"]))
pr("Classical ~ restCAF", co_class_rest); pr("Classical ~ proCAF ", co_class_pro)
pr("Basal ~ proCAF     ", co_basal_pro);  pr("Basal ~ restCAF    ", co_basal_rest)
pr("D1 ~ restCAF       ", co_D1_rest);    pr("D1 ~ proCAF        ", co_D1_pro)
cat("\n(2) Spatial adjacency -- restCAF-minus-proCAF signal in fibroblast NEIGHBOURS:\n")
cat(sprintf("    next to Classical epithelium: % .3f (n=%d spots)\n", mean(class_diff,na.rm=TRUE), length(class_diff)))
cat(sprintf("    next to Basal epithelium    : % .3f (n=%d spots)\n", mean(basal_diff,na.rm=TRUE), length(basal_diff)))
cat(sprintf("    gap (Classical - Basal) = %.3f | Wilcoxon p %s | within-sample permutation p %s\n",
            obs_gap, fmtp(adj_test$p.value), fmtp(perm_p)))
cat("\nSaved -> results/spatial_cooccurrence_stats.rds , results/spatial_spots_scored.rds\n")
