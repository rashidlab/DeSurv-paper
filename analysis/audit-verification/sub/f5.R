suppressMessages({ library(survival); library(superpc); library(DeSurv) })
dat <- readRDS("results/tar_data_filtered_tcgacptac.rds"); ex <- dat$ex; si <- dat$sampInfo; y <- si$time; d <- si$event
fit <- readRDS("results/tar_fit_desurv_tcgacptac.rds"); W <- fit$W
dv <- readRDS("results/data_val_filtered_tcgacptac.rds")
genes <- rownames(W); Xtr <- t(ex[genes,,drop=FALSE]); gsd <- apply(Xtr,2,sd); Xtr <- Xtr[,gsd>0,drop=FALSE]; genes <- colnames(Xtr)
cat("genes", length(genes), "\n")
source("R/paca_dedup.R"); source("R/get_top_genes.R")
.seqacc19 <- { s <- Filter(function(co) co$sampInfo$dataset[1] == "PACA_AU_seq", dv); co <- s[[1]]; kp <- if (!is.null(co$samp_keeps)) co$samp_keeps else which(co$sampInfo$keep==1); paca_accession(rownames(co$sampInfo)[kp]) }
cidx <- function(s,t,e) unname(survival::concordance(Surv(t,e)~s)$concordance)
proj <- function(E,w){ g<-intersect(names(w),rownames(E)); as.numeric(drop(t(E[g,,drop=FALSE])%*%w[g])) }
dtr <- list(x=t(Xtr), y=y, censoring.status=d, featurenames=genes)
sp <- superpc.train(dtr, type="survival")
set.seed(1); cv <- superpc.cv(sp,dtr,n.threshold=20,n.fold=5,min.features=5,max.features=length(genes))
lr<-cv$scor; lr[is.na(lr)]<- -Inf; bi<-which(lr==max(lr),arr.ind=TRUE)[1,]; ncomp<-as.integer(bi[1]); thr<-cv$thresholds[bi[2]]
keepg <- abs(sp$feature.scores)>=thr
cat("ncomp", ncomp, "thr", thr, "n genes", sum(keepg), "\n")
pr <- superpc.predict(sp,dtr,dtr,threshold=thr,n.components=ncomp,prediction.type="continuous")
V <- as.matrix(pr$v.pred); bcox <- coef(coxph(Surv(y,d)~V)); scr <- as.numeric(V %*% bcox)
w <- setNames(rep(0,length(genes)),genes); w[keepg] <- apply(Xtr[,keepg,drop=FALSE],2,function(g)cor(g,scr))
sur_tr <- proj(ex,w); cat("cor(native train, surrogate train)", cor(scr, sur_tr), "\n")
Zt1 <- proj(ex, { tg <- get_top_genes(W,270)$top_genes; v<-setNames(rep(0,length(genes)),genes); g1<-intersect(tg[[1]],genes); v[g1]<-W[g1,1]; v })
sgn <- if (cor(sur_tr, Zt1) < 0) -1 else 1; cat("orientation flip for surrogate:", sgn, " native-vs-D1 cor sign", sign(cor(scr,Zt1)), "\n")
res <- NULL
for (co in dv) { keep<-if(!is.null(co$samp_keeps))co$samp_keeps else which(co$sampInfo$keep==1)
  ids<-rownames(co$sampInfo)[keep]; dsn<-co$sampInfo$dataset[1]
  ded<-if(dsn=="PACA_AU_array") !(paca_accession(ids) %in% .seqacc19) else rep(TRUE,length(ids))
  E<-co$ex[,keep,drop=FALSE]; tt<-co$sampInfo$time[keep]; ee<-co$sampInfo$event[keep]; ok<-is.finite(tt)&tt>0&is.finite(ee)&ded
  E <- E[,ok,drop=FALSE]; miss <- setdiff(genes, rownames(E))
  dte <- list(x=E[genes,,drop=FALSE], y=tt[ok], censoring.status=ee[ok], featurenames=genes)
  pv <- superpc.predict(sp,dtr,dte,threshold=thr,n.components=ncomp,prediction.type="continuous")
  nat <- as.numeric(as.matrix(pv$v.pred) %*% bcox)
  s_sur <- sgn*proj(E,w)
  res <- rbind(res, data.frame(ds=dsn, n=sum(ok), ev=sum(ee[ok]), miss=length(miss), native=cidx(-nat,tt[ok],ee[ok]), surrogate=cidx(s_sur,tt[ok],ee[ok]))) }
res <- rbind(res, data.frame(ds="POOLED", n=sum(res$n), ev=sum(res$ev), miss=NA, native=sum(res$native*res$ev)/sum(res$ev), surrogate=sum(res$surrogate*res$ev)/sum(res$ev)))
print(res, digits=7)
c19 <- readRDS("results/desurv_vs_supervised_tuned.rds"); cat("cached names:", names(c19), "\n")
v <- c19$val; if (!is.null(v)) print(v["Supervised PCA",,drop=FALSE], digits=7)
