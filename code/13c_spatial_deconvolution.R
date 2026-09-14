#!/usr/bin/env Rscript
# 13c_spatial_deconvolution.R
# ---------------------------------------------------------------------------
# Robust spatial validation of the tumour-stroma coupling via REFERENCE-BASED
# DECONVOLUTION -- the field-standard way to handle Visium multicellularity +
# dropout (cf. RCTD/cell2location/SPOTlight). Instead of per-spot marker calls or
# rank-projection (both dropout-fragile), we deconvolve each Visium spot into
# cell-type PROPORTIONS by NNLS against an Elyada PDAC scRNA reference (additive
# mixture model; borrows strength jointly across all marker genes). Then we test
# the coupling by Squidpy-style neighbourhood enrichment with a within-section
# permutation null: are Classical-malignant-rich spots spatially adjacent to
# restCAF-rich stroma, and Basal-rich spots to proCAF-rich stroma?
# ---------------------------------------------------------------------------
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(Matrix); library(quadprog) })
set.seed(1)
setwd("/home/naimrashid/Downloads/DeSurv-paper-clean")
STDIR <- "/tmp/unc_st"
REF   <- "/home/naimrashid/Downloads/DeSurv-paper/data/derv/Elyada_umap.Rds"

## --- 1. Build reference cell-type profiles from Elyada -------------------------
cat("building Elyada reference profiles...\n")
o   <- readRDS(REF); md <- attr(o,"meta.data"); cnt <- attr(attr(o,"assays")$RNA,"counts")
lf  <- as.character(md$label_fine)
ct  <- ifelse(grepl("Classical|Ductal", lf), "Classical",
       ifelse(grepl("Basal", lf),            "Basal",
       ifelse(lf=="iCAF",                     "restCAF",
       ifelse(lf=="myCAF",                    "proCAF",
       ifelse(lf=="apCAF",                    "apCAF",
       ifelse(grepl("Acinar|Normal Stroma|Fibroblast", lf), "Acinar_Normal",
              "Immune"))))))
types <- c("Classical","Basal","restCAF","proCAF","apCAF","Acinar_Normal","Immune")
# mean CPM profile per cell type (CPM = counts scaled to 1e4 per cell, then averaged)
cs <- Matrix::colSums(cnt); cs[cs==0] <- 1; cpm <- cnt %*% Diagonal(x = 1e4/cs)
ref <- sapply(types, function(t) Matrix::rowMeans(cpm[, ct==t, drop=FALSE]))
rownames(ref) <- rownames(cnt)
# informative marker genes: top 60 by fold-change of each type vs the rest
markers <- unique(unlist(lapply(types, function(t){
  fc <- log2((ref[,t]+1) / (rowMeans(ref[, setdiff(types,t), drop=FALSE])+1))
  names(sort(fc, decreasing=TRUE))[1:60] })))
cat("reference:", length(types), "types,", length(markers), "marker genes\n")

## --- 2. NNLS deconvolution of each Visium spot --------------------------------
nnls_qp <- function(R, y){                      # min ||y - R p||^2 s.t. p >= 0
  D <- crossprod(R) + diag(1e-6, ncol(R)); dv <- crossprod(R, y)
  s <- tryCatch(quadprog::solve.QP(D, dv, diag(ncol(R)), rep(0, ncol(R)))$solution,
                error=function(e) rep(0, ncol(R)))
  s <- pmax(s, 0); if (sum(s) > 0) s/sum(s) else s }

samples <- sub("_matrix.mtx.gz$","", list.files(STDIR, pattern="matrix.mtx.gz$"))
all_spots <- list()
for (s in samples) {
  M  <- as(Matrix::readMM(gzfile(file.path(STDIR, paste0(s,"_matrix.mtx.gz")))), "CsparseMatrix")
  ft <- read.delim(gzfile(file.path(STDIR, paste0(s,"_features.tsv.gz"))), header=FALSE)
  rownames(M) <- make.unique(as.character(ft$V2))
  bc <- readLines(gzfile(file.path(STDIR, paste0(s,"_barcodes.tsv.gz")))); colnames(M) <- bc
  pos <- read.csv(gzfile(file.path(STDIR, paste0(s,"_tissue_positions.csv.gz")))); pos <- pos[pos$in_tissue==1, ]
  keep <- intersect(colnames(M), pos$barcode); M <- M[, keep, drop=FALSE]; pos <- pos[match(keep, pos$barcode), ]
  if (ncol(M) < 50) next
  g <- intersect(markers, rownames(M)); R <- as.matrix(ref[g, , drop=FALSE]); Y <- as.matrix(M[g, , drop=FALSE])
  P <- t(apply(Y, 2, function(y) nnls_qp(R, y)))         # spots x types proportions
  colnames(P) <- types
  df <- data.frame(sample=sub("^GSM[0-9]+_","",s), row=pos$array_row, col=pos$array_col, P,
                   depth=Matrix::colSums(M>0), stringsAsFactors=FALSE)
  all_spots[[s]] <- df
  cat(sprintf("  %-10s %4d spots deconvolved\n", df$sample[1], nrow(df)))
}
spots <- do.call(rbind, all_spots); rownames(spots) <- NULL
cat("total spots:", nrow(spots), "\n\n")

## --- 3. Coupling via neighbourhood enrichment (deconvolution proportions) -----
# epithelial-dominant spot = Classical+Basal proportion exceeds stroma; classify by which tumour type;
# for its fibroblast neighbours, measure restCAF - proCAF proportion.
spots$epi <- spots$Classical + spots$Basal; spots$fibp <- spots$restCAF + spots$proCAF
spots$is_epi <- spots$epi > spots$fibp & spots$epi > 0.2
spots$is_fib <- spots$fibp > spots$epi & spots$fibp > 0.2
spots$epi_type <- ifelse(spots$is_epi, ifelse(spots$Classical > spots$Basal, "Classical", "Basal"), NA)

neigh_off <- list(c(-1,-1),c(-1,1),c(1,-1),c(1,1),c(0,-2),c(0,2))
adj <- list()
for (g in unique(spots$sample)) {
  sp <- spots[spots$sample==g, ]; idx <- setNames(seq_len(nrow(sp)), paste(sp$row, sp$col, sep="_"))
  for (i in which(sp$is_epi)) {
    nb <- na.omit(sapply(neigh_off, function(o) idx[paste(sp$row[i]+o[1], sp$col[i]+o[2], sep="_")]))
    nbf <- nb[sp$is_fib[nb]]
    if (length(nbf) >= 1) adj[[length(adj)+1]] <- data.frame(sample=g, epi_type=sp$epi_type[i],
        nb_rest=mean(sp$restCAF[nbf]), nb_pro=mean(sp$proCAF[nbf]), depth=sp$depth[i]) }
}
adj <- do.call(rbind, adj); adj$nb_diff <- adj$nb_rest - adj$nb_pro
cd <- adj$nb_diff[adj$epi_type=="Classical"]; bd <- adj$nb_diff[adj$epi_type=="Basal"]
obs <- mean(cd,na.rm=T) - mean(bd,na.rm=T)
nperm <- 2000; pg <- numeric(nperm)
for (p in 1:nperm) { lab <- ave(adj$epi_type, adj$sample, FUN=function(z) sample(z))
  pg[p] <- mean(adj$nb_diff[lab=="Classical"],na.rm=T) - mean(adj$nb_diff[lab=="Basal"],na.rm=T) }
perm_p <- (1 + sum(abs(pg) >= abs(obs)))/(nperm+1)

res <- list(meta=list(n_spots=nrow(spots), n_epi=sum(spots$is_epi), n_fib=sum(spots$is_fib),
                      types=types, n_markers=length(markers)),
            mean_prop=round(colMeans(spots[,types]),3),
            adjacency=list(class_nb_diff=mean(cd,na.rm=T), basal_nb_diff=mean(bd,na.rm=T),
                           gap=obs, wilcox_p=wilcox.test(cd,bd)$p.value, perm_p=perm_p,
                           n_class=length(cd), n_basal=length(bd)),
            spots=spots)
saveRDS(res, "results/spatial_deconv_stats.rds")

cat("=== DECONVOLUTION-BASED COUPLING (Elyada reference, NNLS) ===\n")
cat("mean spot composition:\n"); print(res$mean_prop)
cat(sprintf("\nNeighbourhood enrichment (restCAF - proCAF proportion in fibroblast neighbours):\n"))
cat(sprintf("  next to Classical epithelium: %+.3f (n=%d)\n", mean(cd,na.rm=T), length(cd)))
cat(sprintf("  next to Basal epithelium    : %+.3f (n=%d)\n", mean(bd,na.rm=T), length(bd)))
cat(sprintf("  gap = %+.3f | Wilcoxon p = %.3g | within-section permutation p = %.4g\n",
            obs, wilcox.test(cd,bd)$p.value, perm_p))
cat("Saved -> results/spatial_deconv_stats.rds\n")
