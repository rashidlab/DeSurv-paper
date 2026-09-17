#!/usr/bin/env Rscript
# 28_tuned_predictor_decomposition.R
# ---------------------------------------------------------------------------
# Source-separation decomposition of the EXACT tuned predictors whose pooled
# external C-indices are reported in the Results (code/19, seed-1 reference
# fits: DeSurv full LP, sparse Cox lambda.min, superpc CV-tuned).  Predictors
# are gene-weight vectors w applied to raw (rank-transformed) expression,
# score = w'X (code/19 proj()).  Under the fit X ~ WH,  score ~ a'H with
# a = W'w, so the score's program-space component is a'H and:
#   * R2_program_space     = R2 of score on (H1,H2,H3)         [training fit]
#   * semipartial R2 D2D3  = R2(H1,H2,H3) - R2(H1)  = variance moved by D2/D3
#                            at fixed D1, as a fraction of total score variance
#   * within-program shares: |a_j| sd(H_j) / sum_j |a_j| sd(H_j) (marginal
#     sensitivities; H rows are correlated so no exact additive split exists),
#     plus each program's own semipartial R2 given the other two.
# Identity check: transfer()["POOLED"] must equal desurv_vs_supervised_tuned
# $prediction[, "POOLED"] for every predictor, else stop.
# D3 carries a positive Cox coefficient in DeSurv: dependence on D3 is a genuine
# secondary prognostic source, not noise.  Output:
# results/tuned_predictor_decomposition.rds
# ---------------------------------------------------------------------------
suppressMessages({ library(survival); library(glmnet); library(superpc) })
source("R/get_top_genes.R")
QUICK <- identical(Sys.getenv("DESURV_QUICK"), "TRUE")
NTHR  <- if (QUICK) 10L else 20L     # must match code/19 for the seed-1 fit to be identical

## ---- identical setup to code/19 ------------------------------------------
dat <- readRDS("results/tar_data_filtered_tcgacptac.rds")
ex  <- dat$ex; si <- dat$sampInfo; y <- si$time; d <- si$event
fit <- readRDS("results/tar_fit_desurv_tcgacptac.rds"); W <- fit$W; H <- fit$H
dv  <- readRDS("results/data_val_filtered_tcgacptac.rds")
genes <- rownames(W); Xtr <- t(ex[genes,,drop=FALSE])
gsd <- apply(Xtr,2,sd); Xtr <- Xtr[,gsd>0,drop=FALSE]; genes <- colnames(Xtr)
source("R/paca_dedup.R")
.seqacc19 <- { s <- Filter(function(co) co$sampInfo$dataset[1] == "PACA_AU_seq", dv)
  if (length(s)) { co <- s[[1]]; kp <- if (!is.null(co$samp_keeps)) co$samp_keeps else which(co$sampInfo$keep==1)
    paca_accession(rownames(co$sampInfo)[kp]) } else character(0) }
proj <- function(E,w){ g<-intersect(names(w),rownames(E)); as.numeric(drop(t(E[g,,drop=FALSE])%*%w[g])) }
cidx <- function(s,t,e) unname(survival::concordance(Surv(t,e)~s)$concordance)
transfer <- function(w){ ci<-nev<-c(); nm<-c()
  for(co in dv){ keep<-if(!is.null(co$samp_keeps))co$samp_keeps else which(co$sampInfo$keep==1)
    ids<-rownames(co$sampInfo)[keep]; dsn<-co$sampInfo$dataset[1]
    ded<-if(dsn=="PACA_AU_array") !(paca_accession(ids) %in% .seqacc19) else rep(TRUE,length(ids))
    E<-co$ex[,keep,drop=FALSE];tt<-co$sampInfo$time[keep];ee<-co$sampInfo$event[keep]
    ok<-is.finite(tt)&tt>0&is.finite(ee)&ded
    ci<-c(ci,cidx(proj(E[,ok,drop=FALSE],w),tt[ok],ee[ok]));nev<-c(nev,sum(ee[ok]));nm<-c(nm,dsn) }
  c(setNames(ci,nm), POOLED=sum(ci*nev)/sum(nev)) }

tg <- get_top_genes(W,270)$top_genes
Wtr <- matrix(0,length(genes),3,dimnames=list(genes,colnames(W)))
for(j in 1:3){ gj<-intersect(tg[[j]],genes); Wtr[gj,j]<-W[gj,j] }
w_D1 <- Wtr[,1]
Zt <- sapply(1:3,function(j) proj(ex,Wtr[,j])); colnames(Zt)<-c("D1","D2","D3")
cox_lp <- coxph(Surv(y,d)~Zt[,1]+Zt[,2]+Zt[,3]); b_lp <- coef(cox_lp)
w_LP <- setNames(as.numeric(Wtr%*%b_lp),genes)

dtr <- list(x=t(Xtr), y=y, censoring.status=d, featurenames=genes)
sp  <- superpc.train(dtr, type="survival")
superpc_pick <- function(seed){ set.seed(seed)
  cv <- superpc.cv(sp,dtr,n.threshold=NTHR,n.fold=5,min.features=5,max.features=length(genes))
  bi <- which(cv$scor==max(cv$scor,na.rm=TRUE),arr.ind=TRUE)[1,]
  list(ncomp=as.integer(bi[1]), thr=cv$thresholds[bi[2]], scor_by_ncomp=apply(cv$scor,1,max,na.rm=TRUE)) }
superpc_dir <- function(pk){ pr<-superpc.predict(sp,dtr,dtr,threshold=pk$thr,n.components=pk$ncomp,prediction.type="continuous")
  V<-as.matrix(pr$v.pred); scr<-as.numeric(V %*% coef(coxph(Surv(y,d) ~ V)))
  keepg<-abs(sp$feature.scores)>=pk$thr
  w<-setNames(rep(0,length(genes)),genes); w[keepg]<-apply(Xtr[,keepg,drop=FALSE],2,function(g)cor(g,scr))
  list(w=w, genes=genes[keepg], score=scr) }
set.seed(1); cvg <- cv.glmnet(Xtr,Surv(y,d),family="cox",alpha=1,nfolds=10)
w_SC <- setNames(as.numeric(coef(cvg,s="lambda.min")),genes)
pk1 <- superpc_pick(1); sd1 <- superpc_dir(pk1); w_SP <- sd1$w
al <- function(w) if (cor(proj(ex,w),proj(ex,w_D1))<0) -w else w
M <- lapply(list("DeSurv full LP"=w_LP,"DeSurv D1"=w_D1,"Sparse Cox"=w_SC,"Supervised PCA"=w_SP), al)

## ---- identity check against the reported C-index comparison ---------------
ref <- readRDS("results/desurv_vs_supervised_tuned.rds")
val <- t(sapply(M, transfer))
ident <- data.frame(predictor = rownames(val), pooled_C_recomputed = val[, "POOLED"],
                    pooled_C_reported = ref$prediction[rownames(val), "POOLED"])
ident$match <- round(ident$pooled_C_recomputed, 3) == ident$pooled_C_reported  # code/19 stores round(., 3)
print(ident, row.names = FALSE)
if (!all(ident$match)) stop("Rebuilt predictors do not reproduce the reported pooled C-indices; aborting.")
stopifnot(pk1$ncomp == ref$stability$sp_ncomp[1], length(sd1$genes) == ref$stability$sp_ngene[1],
          sum(w_SC != 0) == ref$stability$scox_ngene[1])

## ---- decomposition ---------------------------------------------------------
S  <- cov(t(H)); Hsd <- sqrt(diag(S))
decomp <- function(w){
  w <- w[genes]; a <- as.numeric(t(W[genes,]) %*% w)
  sc <- proj(ex, w); dd <- data.frame(sc, H1=H[1,], H2=H[2,], H3=H[3,])
  R2 <- function(f) summary(lm(f, dd))$r.squared
  r123 <- R2(sc~H1+H2+H3); r1 <- R2(sc~H1); r12 <- R2(sc~H1+H2); r13 <- R2(sc~H1+H3); r23 <- R2(sc~H2+H3)
  sens <- abs(a)*Hsd
  data.frame(n_genes = sum(w != 0),
             R2_program_space = r123,
             semipartial_R2_D2D3_given_D1 = r123 - r1,
             semipartial_R2_D1_given_D2D3 = r123 - r23,
             semipartial_R2_D2_given_D1D3 = r123 - r13,
             semipartial_R2_D3_given_D1D2 = r123 - r12,
             within_share_D1 = sens[1]/sum(sens), within_share_D2 = sens[2]/sum(sens), within_share_D3 = sens[3]/sum(sens),
             cor_H1 = cor(sc,H[1,]), cor_H2 = cor(sc,H[2,]), cor_H3 = cor(sc,H[3,])) }
dec <- do.call(rbind, lapply(names(M), function(nm) cbind(predictor = nm, decomp(M[[nm]]))))
rownames(dec) <- NULL

## sparse Cox genes by dominant program (H-mean scaled shares), for the record
C <- sweep(W,2,rowMeans(H),`*`); P <- C/pmax(rowSums(C),1e-300)
dom <- setNames(apply(P,1,which.max), rownames(W)); spec <- setNames(apply(P,1,max), rownames(W))
scg <- names(w_SC)[w_SC != 0]
scox_tab <- data.frame(gene = scg, beta_oriented = M[["Sparse Cox"]][scg], dominant = dom[scg], spec_max = round(spec[scg],3),
                       in_D1_exemplar = scg %in% tg$factor1, in_D3_exemplar = scg %in% tg$factor3, row.names = NULL)
sp_tab <- data.frame(gene = sd1$genes, w_oriented = M[["Supervised PCA"]][sd1$genes], dominant = dom[sd1$genes],
                     spec_max = round(spec[sd1$genes],3), row.names = NULL)

res <- list(meta = list(source = "code/19 seed-1 reference fits, rebuilt and identity-checked",
                        n_train = ncol(ex), n_genes = length(genes), NTHR = NTHR,
                        superpc = list(ncomp = pk1$ncomp, threshold = pk1$thr, n_genes = length(sd1$genes)),
                        sparse_cox_n_genes = sum(w_SC != 0),
                        desurv_lp_coef = b_lp, desurv_beta_fit = as.numeric(fit$beta),
                        H_cor = cor(t(H))),
            identity = ident, decomposition = dec, sparse_cox_genes = scox_tab, superpc_genes = sp_tab,
            weights = M)
saveRDS(res, "results/tuned_predictor_decomposition.rds")
cat("\n== DeSurv LP Cox coefficients on Z1..Z3 (re-estimated):", signif(b_lp,3), "\n")
cat("== superpc seed 1: ncomp", pk1$ncomp, " thr", round(pk1$thr,2), " genes", length(sd1$genes), "; sparse Cox genes", sum(w_SC!=0), "\n\n")
print(format(dec, digits = 3), row.names = FALSE)
cat("\nSparse Cox genes: sign x dominant program\n"); print(table(sign = sign(scox_tab$beta_oriented), dominant = scox_tab$dominant))
cat("\nsuperpc genes:\n"); print(sp_tab, row.names = FALSE)
