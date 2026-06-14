#!/usr/bin/env Rscript
# 16_desurv_vs_supervised.R
# ---------------------------------------------------------------------------
# DeSurv vs supervised-only methods (Nature Cancer revision).
#
# Question: what does DeSurv provide over survival-supervised-only methods? Truth =
# clinical outcome, so methods are compared on outcome-grounded properties, NOT on
# agreement with any one method's gene list.
#
# FINDINGS (this script):
#  - Prediction is PARITY among supervised methods (DeSurv full LP ~0.62, sparse Cox,
#    Cox-PLS, supervised PCA all ~0.60-0.64); unsupervised baselines fail (~0.49-0.55).
#    So prediction is NOT the differentiator.
#  - Reproducibility is NOT a DeSurv advantage: on cross-cohort gene transport DeSurv
#    (~0.11) ties the dense supervised methods Cox-PLS/SupPCA; only pure sparse outcome-
#    selection (sparse Cox, ~0.03) is clearly worse. (An earlier "within->cross reversal"
#    was a gene-wise-z-scoring artifact, removed by the unified within-sample-rank space.)
#  - The differentiator therefore is COMPARTMENT RESOLUTION (source-attributable programs,
#    validated elsewhere by single-cell + spatial), which dense supervised scores lack.
#
# Four parts:
#  (1) PREDICTION PARITY: freeze each method's leading prognostic axis on TCGA+CPTAC,
#      project to the 5 validation cohorts, compare frozen transfer C-index.
#  (2) REPRODUCIBILITY (within): refit on 80% subsamples (without replacement) of training,
#      mean pairwise top-100 gene Jaccard.
#  (3) ALPHA-SWEEP: vary DeSurv's supervision weight alpha; prognosis (CV) vs gene stability.
#  (4) REPRODUCIBILITY (cross-cohort): refit each method independently per cohort
#      (composition differs) and measure top-100 gene transport across cohorts.
#
# DeSurv loads from the R-4.5 library under R 4.6 (compiled deps compatible).
# ---------------------------------------------------------------------------

.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(survival); library(glmnet); library(DeSurv) })
set.seed(1)

# quick-mode contract: DESURV_QUICK=TRUE reduces iteration knobs for a fast smoke test
QUICK <- toupper(Sys.getenv("DESURV_QUICK")) %in% c("TRUE", "1", "YES")
MAXIT  <- if (QUICK) 60L else 200L   # DeSurv fit iterations
N_SUB2 <- if (QUICK) 4L  else 25L    # Part 2 within-cohort subsamples
N_SUB3 <- if (QUICK) 3L  else 15L    # Part 3 alpha-sweep subsamples
NTOP <- 270L   # gene-screen / top-gene count, matching DeSurv's top-270 (PGENES) convention
if (QUICK) cat("[DESURV_QUICK] reduced knobs: MAXIT=60, N_SUB2=4, N_SUB3=3, 3-point alpha grid\n")

dat <- readRDS("results/tar_data_filtered_tcgacptac.rds")
ex  <- dat$ex; si <- dat$sampInfo
fit <- readRDS("results/tar_fit_desurv_tcgacptac.rds"); W <- fit$W; Wg <- rownames(W)
sh  <- intersect(Wg, rownames(ex))
refD1 <- W[sh, 1]                                  # reference D1 gene loadings
HP  <- readRDS("results/tar_params_best_tcgacptac.rds")

# uniform scoring: within-sample rank over sh, then project on a gene-weight vector
rankM <- function(E){ g <- intersect(sh, rownames(E))
  xr <- apply(as.matrix(E[g, , drop=FALSE]), 2, rank, ties.method="average"); rownames(xr) <- g; xr }
score <- function(xr, w){ g <- intersect(names(w), rownames(xr))
  as.numeric(drop(t(xr[g, , drop=FALSE]) %*% w[g])) }
# FROZEN-model transfer C-index: direct rank concordance of the (training-aligned, protective)
# score against validation outcome. survival::concordance() does NOT fit a Cox model on
# validation -- it is a direct concordance estimator -- and for our protective-aligned scores it
# returns the correctly-oriented discrimination (>0.5 for a transferring protective score).
cidx  <- function(s, time, event){ unname(survival::concordance(Surv(time, event) ~ s)$concordance) }
topG  <- function(w, n=100) names(sort(abs(w), decreasing=TRUE))[seq_len(min(n,sum(w!=0)))]
jacc  <- function(a,b){ length(intersect(a,b))/length(union(a,b)) }
matchD1 <- function(Wfit){                          # pick + sign-align the factor closest to reference D1
  cc <- apply(Wfit[sh,,drop=FALSE], 2, function(v) suppressWarnings(cor(v, refD1, method="spearman", use="complete.obs")))
  j <- which.max(abs(cc)); v <- Wfit[sh, j]; if (cc[j] < 0) v <- -v; setNames(v, sh) }

## --- method fitters: each returns a named gene-weight vector over sh ----------
# Training matrix for the SUPERVISED-ONLY comparators (Sparse Cox, Cox-PLS, Supervised PCA,
# Unsup. PCA) = within-sample ranks over sh (samples x genes) -- the SAME representation used to
# score validation (rankM), so learned weights apply directly to rankM(E) with no train/validation
# transform mismatch. glmnet standardizes internally (coefs on input scale); PCA steps center below.
# (DeSurv is NOT fit via Xtrain: it uses its cached consensus W / desurv_fit on the rank-transformed
# expression and scores by the standard rank-then-project convention -- self-consistent with the
# manuscript, so the "unified rank space" statement scopes to the comparators, not DeSurv.)
Xtrain <- function(E, idx=NULL){ M <- t(rankM(E))
  if(!is.null(idx)) M <- M[idx,,drop=FALSE]; gsd <- apply(M,2,sd); M[, gsd>0, drop=FALSE] }

fit_desurv <- function(E, y, d, alpha){
  f <- desurv_fit(X=E, y=y, d=d, k=3, alpha=alpha, lambda=HP$lambda, nu=HP$nu,
                  lambdaW=0, lambdaH=0, maxit=MAXIT, imaxit=MAXIT, ninit=1)
  matchD1(f$W) }
fit_scox <- function(Xg, y, d, lambda=NULL){
  if (is.null(lambda)) lambda <- cv.glmnet(Xg, Surv(y,d), family="cox", alpha=1, nfolds=10)$lambda.min
  b <- as.numeric(coef(glmnet(Xg, Surv(y,d), family="cox", alpha=1, lambda=lambda)))
  setNames(b, colnames(Xg)) }
fit_pls <- function(Xg, y, d){ r <- residuals(coxph(Surv(y,d) ~ 1), type="martingale")
  w <- as.numeric(t(Xg) %*% r); setNames(w/sqrt(sum(w^2)), colnames(Xg)) }
fit_spca <- function(Xg, y, d){ z <- apply(Xg, 2, function(g) summary(coxph(Surv(y,d)~g))$coefficients[1,"z"])
  keep <- order(abs(z), decreasing=TRUE)[1:min(NTOP,ncol(Xg))]
  pcl <- prcomp(Xg[,keep], center=TRUE, scale.=FALSE)$rotation[,1]
  w <- setNames(rep(0, ncol(Xg)), colnames(Xg)); w[keep] <- pcl; w }
fit_pca <- function(Xg, y, d){ setNames(prcomp(Xg, center=TRUE, scale.=FALSE)$rotation[,1], colnames(Xg)) }

align <- function(w) if (suppressWarnings(cor(w[sh], refD1, use="complete.obs")) < 0) -w else w

## ===========================================================================
## (1) PREDICTION PARITY on the 5 validation cohorts
## ===========================================================================
cat("== Part 1: out-of-sample prediction ==\n")
y <- si$time; d <- si$event
Xg_full <- Xtrain(ex)
# Full DeSurv 3-factor linear predictor (the FAIR comparator), computed in-script:
# project training onto all 3 reference factors, fit a Cox LP, fold beta back into
# one gene-weight vector w_full = sum_j beta_j * W[,j] so it scores like the others.
Zt <- sapply(1:3, function(j) score(rankM(ex), setNames(W[sh, j], sh)))
b_lp <- coef(coxph(Surv(y, d) ~ Zt[,1] + Zt[,2] + Zt[,3]))
w_full <- setNames(as.numeric(W[sh, 1:3] %*% b_lp), sh)
W_meth <- list(
  "DeSurv (full LP)" = align(w_full),
  "DeSurv (D1)"      = align(setNames(W[sh, 1], sh)),   # consensus D1 (same source as full LP), not a fresh ninit=1 fit
  "Sparse Cox"       = align(fit_scox(Xg_full, y, d)),
  "Cox-PLS"          = align(fit_pls(Xg_full, y, d)),
  "Supervised PCA"   = align(fit_spca(Xg_full, y, d)),
  "Unsup. NMF (a=0)" = align(fit_desurv(ex, y, d, 0)),
  "Unsup. PCA (PC1)" = align(fit_pca(Xg_full, y, d)))

dv <- readRDS("results/data_val_filtered_tcgacptac.rds")
val_ci <- list(); pooled <- list()
for (nm in names(W_meth)) {
  w <- W_meth[[nm]]; cis <- c(); nev <- c(); ps <- data.frame()
  for (co in dv) { keep <- if(!is.null(co$samp_keeps)) co$samp_keeps else which(co$sampInfo$keep==1)
    E <- co$ex[, keep, drop=FALSE]; tt <- co$sampInfo$time[keep]; ee <- co$sampInfo$event[keep]
    ok <- is.finite(tt) & tt>0 & is.finite(ee); E<-E[,ok,drop=FALSE]; tt<-tt[ok]; ee<-ee[ok]
    s <- score(rankM(E), w); cis <- c(cis, cidx(s, tt, ee)); nev <- c(nev, sum(ee))
    ps <- rbind(ps, data.frame(z=as.numeric(scale(s)), time=tt, event=ee, ds=co$sampInfo$dataset[1])) }
  names(cis) <- sapply(dv, function(z) z$sampInfo$dataset[1])
  hr <- summary(coxph(Surv(time,event) ~ z + strata(ds), data=ps))  # association test only (HR), not used for C
  val_ci[[nm]] <- cis
  pooled[[nm]] <- c(pooled_C = sum(cis * nev) / sum(nev),   # event-weighted mean of frozen per-cohort C
                    HR = hr$conf.int["z","exp(coef)"], p = hr$coefficients["z","Pr(>|z|)"]) }
ci_tab <- do.call(rbind, val_ci); pool_tab <- do.call(rbind, pooled)
print(round(cbind(ci_tab, pool_tab), 3))

## ===========================================================================
## (2) REPRODUCIBILITY under data perturbation (80% subsampling, without replacement)
## ===========================================================================
cat("\n== Part 2: representation stability (top-100 gene Jaccard across subsamples) ==\n")
B <- N_SUB2; n <- nrow(Xg_full)
# shared sparse-Cox gene importance (path to >=N nonzero, |coef|) -- SAME rule in Part 2 and Part 4
scox_imp <- function(Xc, t, e, N=100){ gg <- glmnet(Xc, Surv(t,e), family="cox", alpha=1)
  nz <- apply(gg$beta!=0, 2, sum); j <- which(nz>=N)[1]; if (is.na(j)) j <- ncol(gg$beta)
  setNames(abs(as.numeric(gg$beta[, j])), colnames(Xc)) }
subs <- lapply(1:B, function(b) sample(n, round(0.8*n)))
sig_across <- function(fitter){ tops <- lapply(subs, function(idx){
    Xs <- Xtrain(ex, idx); w <- fitter(ex, Xs, idx); topG(w, 100) }); tops }
fitters <- list(
  "DeSurv (D1)"    = function(E,Xs,idx) fit_desurv(E[,idx,drop=FALSE], y[idx], d[idx], HP$alpha),
  "Sparse Cox"     = function(E,Xs,idx) scox_imp(Xs, y[idx], d[idx]),
  "Cox-PLS"        = function(E,Xs,idx) fit_pls(Xs, y[idx], d[idx]),
  "Supervised PCA" = function(E,Xs,idx) fit_spca(Xs, y[idx], d[idx]))
stab <- sapply(names(fitters), function(nm){ tops <- sig_across(fitters[[nm]])
  pj <- combn(length(tops), 2, function(ij) jacc(tops[[ij[1]]], tops[[ij[2]]])); mean(pj) })
cat(sprintf("  %-16s mean top-100 Jaccard = %.3f\n", names(stab), stab))

## ===========================================================================
## (3) ALPHA-SWEEP: prognosis (CV, from grid) and stability vs alpha
## ===========================================================================
cat("\n== Part 3: alpha-sweep (k=3) ==\n")
if (!file.exists("results/cv_grid/cv_grid_summary.csv"))
  stop("results/cv_grid/cv_grid_summary.csv not found -- run code/06_cv_grid.R first ",
       "(it is produced earlier in `make all`).")
grid <- read.csv("results/cv_grid/cv_grid_summary.csv")
g3 <- grid[grid$k==3 & is.na(grid$ntop), c("alpha","mean_cindex")]
# grid-aligned alphas so prognosis lookup and stability refit use the SAME value per row;
# 0.35 is the grid point nearest the BO-selected alpha (0.334)
bo_grid <- g3$alpha[which.min(abs(g3$alpha - HP$alpha))]
alphas <- if (QUICK) sort(unique(c(0, bo_grid, 0.7))) else sort(unique(c(0, 0.1, 0.2, bo_grid, 0.5, 0.7, 0.9)))
B3 <- N_SUB3; subs3 <- lapply(1:B3, function(b) sample(n, round(0.8*n)))
asweep <- data.frame()
for (a in alphas) {
  prog <- g3$mean_cindex[which.min(abs(g3$alpha - a))]
  tops <- lapply(subs3, function(idx) topG(fit_desurv(ex[,idx,drop=FALSE], y[idx], d[idx], a), 100))
  st <- mean(combn(length(tops),2, function(ij) jacc(tops[[ij[1]]], tops[[ij[2]]])))
  # D1 recovery strength at this alpha (full-data fit, loading cor w/ reference D1)
  rec <- suppressWarnings(cor(fit_desurv(ex, y, d, a)[sh], refD1, method="spearman"))
  asweep <- rbind(asweep, data.frame(alpha=a, cv_cindex=prog, stability=st, d1_recovery=abs(rec)))
  cat(sprintf("  alpha=%.2f  CV C=%.3f  stability=%.3f  D1-recovery=%.2f\n", a, prog, st, abs(rec))) }

## ===========================================================================
## (4) CROSS-COHORT REPRODUCIBILITY: refit each method independently on every
## cohort (composition DIFFERS across cohorts) and measure gene-program transport
## (top-100 Jaccard ACROSS cohorts). This is the proper test of whether a
## supervised-only method's gene SELECTION is composition-confounded: within-
## cohort resampling (Part 2) holds composition fixed; cross-cohort does not.
## ===========================================================================
cat("\n== Part 4: cross-cohort gene-program reproducibility ==\n")
cohorts <- list(Train = list(ex=ex, t=y, e=d))
for (co in dv) { keep <- if (!is.null(co$samp_keeps)) co$samp_keeps else which(co$sampInfo$keep==1)
  E <- co$ex[, keep, drop=FALSE]; tt <- co$sampInfo$time[keep]; ee <- co$sampInfo$event[keep]
  ok <- is.finite(tt) & tt>0 & is.finite(ee)
  cohorts[[co$sampInfo$dataset[1]]] <- list(ex=E[,ok,drop=FALSE], t=tt[ok], e=ee[ok]) }
topimp <- function(meth, E, t, e, N=100) {
  Xc <- Xtrain(E)
  w <- switch(meth,
    "DeSurv (D1)"    = abs(fit_desurv(E, t, e, HP$alpha)),
    "Sparse Cox"     = scox_imp(Xc, t, e, N),
    "Cox-PLS"        = abs(fit_pls(Xc, t, e)),
    "Supervised PCA" = abs(fit_spca(Xc, t, e)))
  names(sort(w, decreasing=TRUE))[seq_len(min(N, sum(w>0)))] }
xcohort <- data.frame()
for (mn in c("DeSurv (D1)","Sparse Cox","Cox-PLS","Supervised PCA")) {
  tops <- lapply(names(cohorts), function(cn) topimp(mn, cohorts[[cn]]$ex, cohorts[[cn]]$t, cohorts[[cn]]$e))
  pj <- combn(length(tops), 2, function(ij) jacc(tops[[ij[1]]], tops[[ij[2]]]))
  xcohort <- rbind(xcohort, data.frame(method=mn, cross_cohort_jaccard=mean(pj), within_cohort_jaccard=unname(stab[mn])))
  cat(sprintf("  %-16s cross-cohort Jaccard=%.3f  (within-cohort=%.3f)\n", mn, mean(pj), unname(stab[mn]))) }

## --- save ---
res <- list(meta=list(n=n, B_part2=B, B_part3=B3, bo_alpha=round(HP$alpha,2), bo_grid=bo_grid,
                      cohorts=sapply(cohorts, function(z) length(z$t))),
            val_cindex=ci_tab, pooled=pool_tab, stability=stab, alpha_sweep=asweep,
            cross_cohort=xcohort)
saveRDS(res, "results/desurv_vs_supervised_stats.rds")
cat("\nSaved -> results/desurv_vs_supervised_stats.rds\n")

## --- figure ---
pdf("figures/fig_desurv_vs_supervised.pdf", width=12, height=4.3)
par(mfrow=c(1,3), mar=c(5,4.5,3.6,1.2), mgp=c(2.6,0.7,0))
# A: per-cohort validation C-index (parity among supervised; unsupervised fails; heterogeneity)
mcols <- c("DeSurv (full LP)"="#67001F","DeSurv (D1)"="#B2182B","Sparse Cox"="#1B7837","Cox-PLS"="#5AAE61",
           "Supervised PCA"="#A6DBA0","Unsup. NMF (a=0)"="#999999","Unsup. PCA (PC1)"="#CCCCCC")
plot(NA, xlim=c(1,ncol(ci_tab)), ylim=c(0.45,0.72), xaxt="n", xlab="", ylab="Per-cohort validation C-index",
     main="A  Out-of-sample prediction (per cohort)")
axis(1, at=1:ncol(ci_tab), labels=gsub("_.*","",colnames(ci_tab)), cex.axis=0.8)
abline(h=0.5, lty=3, col="grey70")
for (m in rownames(ci_tab)) lines(1:ncol(ci_tab), ci_tab[m,], type="b", pch=19, cex=0.8, col=mcols[m], lwd=ifelse(grepl("Unsup",m),1,2))
legend("bottomright", names(mcols), col=mcols, lwd=2, bty="n", cex=0.62)
# B: within-cohort vs cross-cohort reproducibility (the reversal)
xc <- res$cross_cohort; bw <- t(as.matrix(xc[,c("within_cohort_jaccard","cross_cohort_jaccard")]))
bp <- barplot(bw, beside=TRUE, names.arg=rep("",ncol(bw)), col=c("#BDBDBD","#B2182B"), border="grey25",
        ylim=c(0,0.5), ylab="Gene-program Jaccard", main="B  Reproducibility: resample vs transport")
text(colMeans(bp), -0.02, xc$method, srt=20, adj=1, xpd=TRUE, cex=0.75)
legend("topright", c("within-cohort (resample)","cross-cohort (transport)"), fill=c("#BDBDBD","#B2182B"), bty="n", cex=0.75)
# C: alpha-sweep dual axis
a<-asweep
plot(a$alpha, a$cv_cindex, type="b", pch=19, col="#2166AC", ylim=range(a$cv_cindex)+c(-.01,.01),
     xlab=expression(alpha~"(supervision weight)"), ylab="CV C-index (prognosis)",
     main="C  Why the blend: prognosis vs stability")
par(new=TRUE); plot(a$alpha, a$stability, type="b", pch=17, col="#1B7837", axes=FALSE, xlab="", ylab="")
axis(4); mtext("within-cohort stability (Jaccard)", side=4, line=2.2, cex=0.7)
abline(v=res$meta$bo_grid, lty=2, col="grey40"); text(res$meta$bo_grid, par("usr")[3], "BO", pos=4, cex=0.8, col="grey30")
legend("right", c("prognosis (CV C)","stability"), col=c("#2166AC","#1B7837"), pch=c(19,17), bty="n", cex=0.75)
invisible(dev.off())
cat("Saved -> figures/fig_desurv_vs_supervised.pdf\n")
