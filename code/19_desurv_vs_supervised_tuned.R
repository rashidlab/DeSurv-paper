#!/usr/bin/env Rscript
# code/19_desurv_vs_supervised_tuned.R
# ---------------------------------------------------------------------------
# Hardened DeSurv-vs-supervised comparison for the Nature Cancer revision.
#
# Question (deliberately narrow): do PROPERLY TUNED, canonical survival-supervised
# methods recover a similar prognostic direction and predict similarly, and if so,
# what biological organization does DeSurv add?
#
# Design choices (see docs / reviewer thread):
#  - Shared 1,970-gene training-only expression/variance prefilter for EVERY method
#    (generic preprocessing, survival- and DeSurv-blind; not circular).
#  - Data is already within-sample rank-transformed (train + validation); project
#    DIRECTLY (no re-ranking, no z-scoring of validation genes).
#  - DeSurv scored exactly as published: BO-selected top-270-per-factor truncation.
#  - Comparators fit from the FULL 1,970 with their OWN training-only tuning:
#      * Sparse Cox  = glmnet family="cox", one lambda.min model (its own genes).
#      * Supervised PCA = superpc (Bair-Tibshirani), CV-tuned threshold + n.components.
#    (Cox-PLS omitted: the plsRcox toolchain is unavailable in this environment and a
#     hand-rolled version is not trustworthy; two canonical methods suffice.)
#  - Frozen external transfer: weights fixed on training; concordance on validation.
#  - DeSurv D1 (dominant program) and full LP (D1-D3 risk combination) kept SEPARATE.
#
# Endpoints: (1) external prediction (frozen transfer C-index); (2) axis/score
# correlation, incl. decomposition across D1/D2/D3; (3) biological compartment
# attribution of each method's genes vs established PDAC programs; (4) superpc
# feature-selection stability across CV seeds (supportive).
#
# Output: results/desurv_vs_supervised_tuned.rds   (read by the manuscript)
# ---------------------------------------------------------------------------
suppressMessages({ library(survival); library(glmnet); library(superpc) })
options(warn = -1)
source("R/get_top_genes.R")

QUICK  <- toupper(Sys.getenv("DESURV_QUICK")) %in% c("TRUE","1","YES")
NSEED  <- if (QUICK) 3L else 12L      # superpc/sparse-Cox stability + attribution seeds
NTHR   <- if (QUICK) 10L else 20L     # superpc CV thresholds
message(sprintf("[19] tuned supervised comparison (NSEED=%d, NTHR=%d)%s", NSEED, NTHR, if (QUICK) " [QUICK]" else ""))

## ---- data (shared training-only 1,970-gene universe; already rank-transformed) ----
dat <- readRDS("results/tar_data_filtered_tcgacptac.rds")
ex  <- dat$ex; si <- dat$sampInfo; y <- si$time; d <- si$event
fit <- readRDS("results/tar_fit_desurv_tcgacptac.rds"); W <- fit$W
dv  <- readRDS("results/data_val_filtered_tcgacptac.rds")
genes <- rownames(W); Xtr <- t(ex[genes,,drop=FALSE])
gsd <- apply(Xtr,2,sd); Xtr <- Xtr[,gsd>0,drop=FALSE]; genes <- colnames(Xtr)

proj <- function(E,w){ g<-intersect(names(w),rownames(E)); as.numeric(drop(t(E[g,,drop=FALSE])%*%w[g])) }
cidx <- function(s,t,e) unname(survival::concordance(Surv(t,e)~s)$concordance)
transfer <- function(w){ ci<-nev<-c(); nm<-c()
  for(co in dv){ keep<-if(!is.null(co$samp_keeps))co$samp_keeps else which(co$sampInfo$keep==1)
    E<-co$ex[,keep,drop=FALSE];tt<-co$sampInfo$time[keep];ee<-co$sampInfo$event[keep];ok<-is.finite(tt)&tt>0&is.finite(ee)
    ci<-c(ci,cidx(proj(E[,ok,drop=FALSE],w),tt[ok],ee[ok]));nev<-c(nev,sum(ee[ok]));nm<-c(nm,co$sampInfo$dataset[1]) }
  c(setNames(ci,nm), POOLED=sum(ci*nev)/sum(nev)) }

## ---- DeSurv: top-270-per-factor truncated (published scoring) ----
tg <- get_top_genes(W,270)$top_genes
Wtr <- matrix(0,length(genes),3,dimnames=list(genes,colnames(W)))
for(j in 1:3){ gj<-intersect(tg[[j]],genes); Wtr[gj,j]<-W[gj,j] }
w_D1 <- Wtr[,1]
Zt <- sapply(1:3,function(j) proj(ex,Wtr[,j])); colnames(Zt)<-c("D1","D2","D3")
b_lp <- coef(coxph(Surv(y,d)~Zt[,1]+Zt[,2]+Zt[,3])); w_LP <- setNames(as.numeric(Wtr%*%b_lp),genes)
s_LP <- proj(ex,w_LP)

## ---- superpc setup (feature scores are deterministic; CV is seed-dependent) ----
dtr <- list(x=t(Xtr), y=y, censoring.status=d, featurenames=genes)
sp  <- superpc.train(dtr, type="survival")
superpc_pick <- function(seed){ set.seed(seed)
  cv <- superpc.cv(sp,dtr,n.threshold=NTHR,n.fold=5,min.features=5,max.features=length(genes))
  lr<-cv$scor; lr[is.na(lr)]<--Inf; bi<-which(lr==max(lr),arr.ind=TRUE)[1,]
  list(ncomp=as.integer(bi[1]), thr=cv$thresholds[bi[2]], scor_by_ncomp=apply(cv$scor,1,max,na.rm=TRUE)) }
superpc_dir <- function(pk){ pr<-superpc.predict(sp,dtr,dtr,threshold=pk$thr,n.components=pk$ncomp,prediction.type="continuous")
  scr<-as.numeric(pr$v.pred[,1]); keepg<-abs(sp$feature.scores)>=pk$thr
  w<-setNames(rep(0,length(genes)),genes); w[keepg]<-apply(Xtr[,keepg,drop=FALSE],2,function(g)cor(g,scr))
  list(w=w, genes=genes[keepg], score=scr) }

## ---- (1)+(2) prediction + axis decomposition (seed 1 reference fits) ----
set.seed(1); cvg <- cv.glmnet(Xtr,Surv(y,d),family="cox",alpha=1,nfolds=10)
w_SC <- setNames(as.numeric(coef(cvg,s="lambda.min")),genes)
pk1 <- superpc_pick(1); sd1 <- superpc_dir(pk1); w_SP <- sd1$w
# Orient to DeSurv D1 (protective): survival::concordance(Surv~x) is >0.5 for a
# protective score, and this reproduces the published DeSurv C-index convention.
al <- function(w) if (cor(proj(ex,w),proj(ex,w_D1))<0) -w else w
M <- lapply(list("DeSurv full LP"=w_LP,"DeSurv D1"=w_D1,"Sparse Cox"=w_SC,"Supervised PCA"=w_SP), al)
val <- t(sapply(M, transfer))
comp_scores <- lapply(M[c("Sparse Cox","Supervised PCA")], function(w) proj(ex,w))
factor_dist <- t(sapply(comp_scores, function(s) c(D1=cor(s,Zt[,1]),D2=cor(s,Zt[,2]),D3=cor(s,Zt[,3]),fullLP=cor(s,s_LP))))

## ---- (3) compartment attribution: signatures ----
load("data/original/cmbSubtypes.RData")
comp <- subtypeGeneList[[3]]
sig <- lapply(comp, function(v) intersect(unique(as.character(unlist(v))),genes))
sig$proCAF  <- intersect(c("IGFL2","NOX4","VSNL1","BICD1","NPR3","ETV1","ITGA11","CNIH3","COL11A1"),genes)
sig$restCAF <- intersect(c("CHRDL1","OGN","PI16","ANK2","ABCA8","TGFBR3","FBLN5","SCARA5","KIAA1217"),genes)
sig <- sig[sapply(sig,length)>=3]
enrich <- function(gs,cs,N=length(genes)) -log10(max(phyper(length(intersect(gs,cs))-1,length(cs),N-length(cs),length(gs),lower.tail=FALSE),1e-300))
hits <- function(gs) names(sig)[sapply(sig,function(cs) enrich(gs,cs))>2]
desurv_attr <- lapply(1:3, function(j) hits(intersect(tg[[j]],genes))); names(desurv_attr)<-c("D1","D2","D3")

## ---- (3)+(4) attribution robustness + superpc stability across seeds ----
stab <- data.frame(); comp_hit_counts <- list("Sparse Cox"=setNames(integer(length(sig)),names(sig)),
                                               "Supervised PCA"=setNames(integer(length(sig)),names(sig)))
ng <- list("Sparse Cox"=integer(),"Supervised PCA"=integer())
for (s in 1:NSEED){ set.seed(s)
  g_sc <- genes[as.numeric(coef(cv.glmnet(Xtr,Surv(y,d),family="cox",alpha=1,nfolds=10),s="lambda.min"))!=0]
  pk <- superpc_pick(s); sdd <- superpc_dir(pk); g_sp <- sdd$genes
  for (c in hits(g_sc)) comp_hit_counts[["Sparse Cox"]][c] <- comp_hit_counts[["Sparse Cox"]][c]+1L
  for (c in hits(g_sp)) comp_hit_counts[["Supervised PCA"]][c] <- comp_hit_counts[["Supervised PCA"]][c]+1L
  ng[["Sparse Cox"]]<-c(ng[["Sparse Cox"]],length(g_sc)); ng[["Supervised PCA"]]<-c(ng[["Supervised PCA"]],length(g_sp))
  stab <- rbind(stab, data.frame(seed=s, sp_ncomp=pk$ncomp, sp_thr=round(pk$thr,2), sp_ngene=length(g_sp),
                sp_pooledC=round(transfer(al(sdd$w))["POOLED"],3), scox_ngene=length(g_sc))) }

## ---- per-gene annotation of the seed-1 reference sets ----
annot_set <- function(gs) sapply(gs, function(g){ a<-names(sig)[sapply(sig,function(cs) g%in%cs)]; if(length(a)==0)"(none)" else paste(a,collapse="/") })

res <- list(
  meta = list(n_genes=length(genes), n_train=nrow(Xtr), n_events=sum(d),
              val_events=sapply(dv,function(z) z$sampInfo$dataset[1]), nseed=NSEED,
              methods=c("DeSurv (top-270 x3)","Sparse Cox (glmnet)","Supervised PCA (superpc)")),
  prediction = round(val,3),
  axis_decomposition = round(abs(factor_dist),2),   # report |r| (orientation-invariant)
  superpc_ncomp_scores = round(pk1$scor_by_ncomp,2),
  desurv_attribution = desurv_attr,
  supervised_attribution_across_seeds = comp_hit_counts,
  genes_per_seed = lapply(ng, function(v) c(median=median(v),min=min(v),max=max(v))),
  stability = stab,
  ref_seed1 = list("Sparse Cox"=annot_set(genes[w_SC!=0]), "Supervised PCA"=annot_set(sd1$genes)),
  signature_sizes = sapply(sig,length))
saveRDS(res, "results/desurv_vs_supervised_tuned.rds")
message("[19] saved results/desurv_vs_supervised_tuned.rds")

cat("\n== prediction (pooled transfer C) ==\n"); print(res$prediction[,"POOLED",drop=FALSE])
cat("\n== axis decomposition (comparator vs D1/D2/D3/fullLP) ==\n"); print(res$axis_decomposition)
cat("\n== attribution across", NSEED, "seeds (# seeds each compartment hit) ==\n")
for (m in names(comp_hit_counts)) cat(sprintf("  %-16s %s | genes/seed med=%d [%d,%d]\n", m,
   paste(sprintf("%s:%d", names(comp_hit_counts[[m]])[comp_hit_counts[[m]]>0], comp_hit_counts[[m]][comp_hit_counts[[m]]>0]),collapse=" "),
   res$genes_per_seed[[m]]["median"], res$genes_per_seed[[m]]["min"], res$genes_per_seed[[m]]["max"]))
cat("  DeSurv:", paste(sprintf("%s={%s}",names(desurv_attr),sapply(desurv_attr,paste,collapse=",")),collapse="  "),"\n")
