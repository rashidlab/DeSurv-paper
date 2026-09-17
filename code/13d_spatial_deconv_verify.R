#!/usr/bin/env Rscript
# 13d_spatial_deconv_verify.R
# Robustness check on the (negative) deconvolution coupling result (13c). Two changes:
#  (1) SECOND REFERENCE DEFINITION: re-label Elyada CAFs as restCAF/proCAF using the
#      actual DeCAF signatures (not the iCAF/myCAF Elyada-clustering proxy), then
#      rebuild the reference and re-deconvolve. Addresses "the proxy was wrong."
#  (2) DEEP-SECTION restriction: repeat the adjacency + spot-level co-occurrence in
#      the 3 highest-depth sections only (best data).
# If the coupling still does not appear, the negative result is robust.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(Matrix); library(quadprog) })
set.seed(1); setwd("/home/naimrashid/Downloads/DeSurv-paper-clean")
STDIR <- "/tmp/unc_st"; REF <- "/home/naimrashid/Downloads/DeSurv-paper/data/derv/Elyada_umap.Rds"
proCAF_sig  <- c("IGFL2","NOX4","VSNL1","BICD1","NPR3","ETV1","ITGA11","CNIH3","COL11A1")
restCAF_sig <- c("CHRDL1","OGN","PI16","ANK2","ABCA8","TGFBR3","FBLN5","SCARA5","KIAA1217")

cat("re-labelling Elyada CAFs by DeCAF restCAF/proCAF signatures...\n")
o<-readRDS(REF); md<-attr(o,"meta.data"); cnt<-attr(attr(o,"assays")$RNA,"counts"); dat<-attr(attr(o,"assays")$RNA,"data")
lf<-as.character(md$label_fine)
is_caf <- lf %in% c("iCAF","myCAF","apCAF") | as.character(md$label_broad)=="Fibroblast"
sig <- function(g){ gg<-intersect(g,rownames(dat)); colMeans(as.matrix(dat[gg,,drop=FALSE])) }
caf_lab <- rep(NA_character_, ncol(cnt))
caf_lab[is_caf] <- ifelse(sig(restCAF_sig)[is_caf] > sig(proCAF_sig)[is_caf], "restCAF", "proCAF")
cat("  re-labelled CAFs: "); print(table(caf_lab))

ct <- ifelse(grepl("Classical|Ductal", lf), "Classical",
      ifelse(grepl("Basal", lf), "Basal",
      ifelse(!is.na(caf_lab), caf_lab,
      ifelse(grepl("Acinar|Normal Stroma", lf), "Acinar_Normal", "Immune"))))
types <- c("Classical","Basal","restCAF","proCAF","Acinar_Normal","Immune")
cs<-Matrix::colSums(cnt); cs[cs==0]<-1; cpm<-cnt %*% Diagonal(x=1e4/cs)
ref <- sapply(types, function(t) Matrix::rowMeans(cpm[, ct==t, drop=FALSE])); rownames(ref)<-rownames(cnt)
markers <- unique(unlist(lapply(types, function(t){ fc<-log2((ref[,t]+1)/(rowMeans(ref[,setdiff(types,t),drop=FALSE])+1))
  names(sort(fc,decreasing=TRUE))[1:80] })))
cat("reference:", length(types), "types,", length(markers), "markers\n")

nnls_qp <- function(R,y){ D<-crossprod(R)+diag(1e-6,ncol(R)); s<-tryCatch(solve.QP(D,crossprod(R,y),diag(ncol(R)),rep(0,ncol(R)))$solution,error=function(e)rep(0,ncol(R))); s<-pmax(s,0); if(sum(s)>0)s/sum(s) else s }
samples <- sub("_matrix.mtx.gz$","", list.files(STDIR, pattern="matrix.mtx.gz$"))
alls<-list()
for (s in samples){ M<-as(Matrix::readMM(gzfile(file.path(STDIR,paste0(s,"_matrix.mtx.gz")))),"CsparseMatrix")
  ft<-read.delim(gzfile(file.path(STDIR,paste0(s,"_features.tsv.gz"))),header=FALSE); rownames(M)<-make.unique(as.character(ft$V2))
  bc<-readLines(gzfile(file.path(STDIR,paste0(s,"_barcodes.tsv.gz")))); colnames(M)<-bc
  pos<-read.csv(gzfile(file.path(STDIR,paste0(s,"_tissue_positions.csv.gz")))); pos<-pos[pos$in_tissue==1,]
  keep<-intersect(colnames(M),pos$barcode); M<-M[,keep,drop=FALSE]; pos<-pos[match(keep,pos$barcode),]; if(ncol(M)<50)next
  g<-intersect(markers,rownames(M)); R<-as.matrix(ref[g,,drop=FALSE]); Y<-as.matrix(M[g,,drop=FALSE])
  P<-t(apply(Y,2,function(y)nnls_qp(R,y))); colnames(P)<-types
  alls[[s]]<-data.frame(sample=sub("^GSM[0-9]+_","",s),row=pos$array_row,col=pos$array_col,P,depth=Matrix::colSums(M>0)) }
spots<-do.call(rbind,alls); rownames(spots)<-NULL

wcor<-function(x,y,grp){zs<-c();ns<-c();for(g in unique(grp)){k<-grp==g&is.finite(x)&is.finite(y);if(sum(k)>30){rr<-suppressWarnings(cor(x[k],y[k],method="spearman"));if(is.finite(rr)){zs<-c(zs,atanh(pmin(pmax(rr,-.999),.999)));ns<-c(ns,sum(k)-3)}}};zb<-sum(zs*ns)/sum(ns);c(rho=tanh(zb),p=2*pnorm(-abs(zb*sqrt(sum(ns)))))}
adj_test<-function(sp){ sp$is_epi<-(sp$Classical+sp$Basal)>(sp$restCAF+sp$proCAF)&(sp$Classical+sp$Basal)>0.2
  sp$is_fib<-(sp$restCAF+sp$proCAF)>(sp$Classical+sp$Basal)&(sp$restCAF+sp$proCAF)>0.2
  sp$epi_type<-ifelse(sp$is_epi,ifelse(sp$Classical>sp$Basal,"Classical","Basal"),NA)
  off<-list(c(-1,-1),c(-1,1),c(1,-1),c(1,1),c(0,-2),c(0,2)); A<-list()
  for(g in unique(sp$sample)){q<-sp[sp$sample==g,];idx<-setNames(seq_len(nrow(q)),paste(q$row,q$col,sep="_"))
    for(i in which(q$is_epi)){nb<-na.omit(sapply(off,function(o)idx[paste(q$row[i]+o[1],q$col[i]+o[2],sep="_")]));nbf<-nb[q$is_fib[nb]]
      if(length(nbf)>=1)A[[length(A)+1]]<-data.frame(epi_type=q$epi_type[i],nb_diff=mean(q$restCAF[nbf])-mean(q$proCAF[nbf]))}}
  A<-do.call(rbind,A); cd<-A$nb_diff[A$epi_type=="Classical"];bd<-A$nb_diff[A$epi_type=="Basal"]
  c(class=mean(cd,na.rm=T),basal=mean(bd,na.rm=T),gap=mean(cd,na.rm=T)-mean(bd,na.rm=T),p=wilcox.test(cd,bd)$p.value) }

cat("\n=== ALL sections (DeCAF-relabelled reference) ===\n")
cat("spot-level: Classical~restCAF rho=",round(wcor(spots$Classical,spots$restCAF,spots$sample)["rho"],3),
    "| Classical~proCAF rho=",round(wcor(spots$Classical,spots$proCAF,spots$sample)["rho"],3),"\n")
a<-adj_test(spots); cat(sprintf("adjacency: Classical-nbhd %+.3f, Basal-nbhd %+.3f, gap %+.3f, Wilcoxon p=%.3g\n",a["class"],a["basal"],a["gap"],a["p"]))
deep<-spots[spots$sample %in% c("ST5473","ST11588","ST5425"),]
cat("\n=== DEEP sections only (ST5473/ST11588/ST5425) ===\n")
cat("spot-level: Classical~restCAF rho=",round(wcor(deep$Classical,deep$restCAF,deep$sample)["rho"],3),
    "| Classical~proCAF rho=",round(wcor(deep$Classical,deep$proCAF,deep$sample)["rho"],3),"\n")
ad<-adj_test(deep); cat(sprintf("adjacency: Classical-nbhd %+.3f, Basal-nbhd %+.3f, gap %+.3f, Wilcoxon p=%.3g\n",ad["class"],ad["basal"],ad["gap"],ad["p"]))
saveRDS(list(spots=spots, all=a, deep=ad), "results/spatial_deconv_verify.rds")
cat("\nSaved -> results/spatial_deconv_verify.rds\n")
