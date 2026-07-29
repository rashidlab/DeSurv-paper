#!/usr/bin/env Rscript
# code/11_treated_cohort_analysis.R
# ---------------------------------------------------------------------------
# Per-cohort treated-disease translational analysis for DeSurv.
#
# REPRODUCIBILITY NOTE: the RASH-ACCEPT and Linehan cohorts are restricted-access
# bulk RNA-seq series held OUTSIDE this repository (under ~/research/rif-research/,
# gitignored). This script reads them, computes all per-cohort summary statistics,
# and writes a single derived object, results/treated_cohort_stats.rds, which IS
# tracked and from which the manuscript loads every treated-cohort numeral via
# `r load_result("treated_cohort_stats")$...` / `read_result(...)`. External users without the raw data
# cannot regenerate the .rds but can inspect/verify the derived values it contains.
#
# Scoring convention (matches main-paper validation, code/05): within-sample rank over the
# shared trained-W namespace, then project on the top-270-per-factor union Wtilde (PGENES);
# z-standardize each factor WITHIN cohort. D1 is reported MARGINAL (basal-classical axis,
# collinear with PurIST; adjusting over-adjusts the same axis); D2 marginal + PurIST+DeCAF-
# adjusted; D3 marginal. Cohorts are NOT pooled (settings differ; D1 heterogeneous).
# NOTE: a D2 random-effects meta is still COMPUTED and stored ($meta) for reference, but it is
# NOT used in the manuscript -- under the top-270 projection the D2 between-cohort heterogeneity
# is moderate (I^2 approx 49%), so the homogeneity that would justify pooling does not hold;
# the manuscript reports per-cohort results plus the full-vs-top-gene concordance only.
# ---------------------------------------------------------------------------
suppressMessages({library(survival)})

# Run from the repository root: paths below are relative (consistent with the rest of the pipeline,
# which is invoked from the repo root via make / run_pipeline.R).
fit  <- readRDS(file.path("results", "tar_fit_desurv_tcgacptac.rds"))
W <- fit$W; Wg <- rownames(W)
canon_path <- Sys.glob("~/research/rif-research/purist-twist/data/external/response_master_canonical.rds")
stopifnot(length(canon_path) == 1)
canon <- readRDS(canon_path)
# DeCAF/PurIST subtype-call objects (paths from the project's shared data-path config)
suppressWarnings(try(source("~/research/rif-research/shared/data_paths.R"), silent = TRUE))
if (!exists("COHORT")) stop("Restricted external data-path config (data_paths.R) not found; ",
  "this script requires the restricted-access RASH-ACCEPT/Linehan data environment. ",
  "Use the tracked results/treated_cohort_stats.rds instead.")

# Projection gene set: top-270-per-factor UNION (the validation convention; code/05 +
# tar_params_best$ntop=270). Rank is taken over the full Wg namespace, then the projection is
# restricted to PGENES, exactly as code/05 does (dv$ex rank-transformed, then ex[top_genes,] %*% W[top_genes,]).
source("R/get_top_genes.R")
PGENES <- unique(unlist(lapply(get_top_genes(fit$W, 270)$top_genes, as.character)))

clean <- function(ex){ rownames(ex) <- gsub("^X(?=[0-9_])", "", rownames(ex), perl = TRUE); ex }
rankX <- function(ex){ ex <- clean(ex); sh <- intersect(Wg, rownames(ex))
  xr <- apply(ex[sh, , drop = FALSE], 2, rank, ties.method = "average"); rownames(xr) <- sh; xr }
projZ <- function(xr, genes, j){ g <- intersect(genes, rownames(xr))
  if (!length(g)) return(setNames(rep(NA_real_, ncol(xr)), colnames(xr)))
  setNames(drop(t(xr[g, , drop = FALSE]) %*% W[g, j, drop = FALSE]), colnames(xr)) }

rash_study <- with(canon$RASH_ACCEPT$extras, setNames(as.character(study), sample))
mk <- function(co, rp, ov = NULL){
  xr <- rankX(t(canon[[co]]$X)); cl <- canon[[co]]$clin; sc <- readRDS(rp)$subtypeCall
  df <- data.frame(sample = colnames(xr), D1 = projZ(xr, PGENES, 1), D2 = projZ(xr, PGENES, 2),
                   D3 = projZ(xr, PGENES, 3),
                   os = cl$os_months[match(colnames(xr), cl$sample)],
                   death = cl$death[match(colnames(xr), cl$sample)],
                   purist = cl$pub_purist[match(colnames(xr), cl$sample)],
                   stringsAsFactors = FALSE)
  df <- merge(df, sc[, intersect(c("sampID", "DeCAF"), names(sc))], by.x = "sample", by.y = "sampID", all.x = TRUE)
  if ("DeCAF" %in% names(df)) df$DeCAF[df$DeCAF == "permCAF"] <- "proCAF"
  df$cohort <- if (is.null(ov)) co else ov[df$sample]
  list(df = df, xr = xr)
}
L  <- mk("Linehan", COHORT$Linehan)
RA <- mk("RASH_ACCEPT", COHORT$Yeh, rash_study)
D <- rbind(L$df, RA$df); D <- D[!is.na(D$os) & !is.na(D$death) & D$os > 0, ]
for (c0 in unique(D$cohort)) for (z in c("D1","D2","D3")) { i <- D$cohort == c0; D[i, z] <- as.numeric(scale(D[i, z])) }
COHS <- c("Linehan","Rash","Accept")

fmtcox <- function(d, sv, adj){
  d$z <- as.numeric(scale(d[[sv]]))
  if (adj){ d <- d[complete.cases(d[, c("purist","DeCAF")]), ]; ext <- ""
    for (a in c("purist","DeCAF")){ d[[a]] <- factor(d[[a]]); if (nlevels(droplevels(d[[a]])) >= 2) ext <- paste0(ext, "+", a) }
  } else ext <- ""
  if (nrow(d) < 10 || sum(d$death) < 5) return(NULL)
  s <- summary(coxph(as.formula(paste0("Surv(os,death)~z", ext)), d))
  data.frame(hr = s$conf.int["z",1], lo = s$conf.int["z",3], hi = s$conf.int["z",4],
             p = s$coef["z",5], n = nrow(d), ev = sum(d$death))
}

## ---- (1) per-cohort OS: D1 marginal; D2,D3 marginal + adjusted ----
os <- list()
for (c0 in COHS){ d <- D[D$cohort == c0, ]
  for (z in c("D1","D2","D3")) for (adj in c(FALSE, TRUE)){
    r <- fmtcox(d, z, adj); if (is.null(r)) next
    r$cohort <- c0; r$program <- z; r$model <- if (adj) "adjusted" else "marginal"
    os[[length(os)+1]] <- r } }
os <- do.call(rbind, os)

## ---- (2) basal-classical classifier triangulation (Classical-vs-Basal OS HR) ----
mof <- with(canon$RASH_ACCEPT$extras, setNames(as.character(moffitt), sample))
bail<- with(canon$RASH_ACCEPT$extras, setNames(as.character(bailey),  sample))
D$moffitt <- mof[D$sample]; D$bailey <- bail[D$sample]
classHR <- function(d, var){
  d <- d[!is.na(d[[var]]) & nzchar(d[[var]]), ]   # drop missing classifier calls BEFORE coding (grepl(.,NA)->FALSE would miscode them)
  if (!nrow(d)) return(NULL)
  if (var == "bailey") g <- factor(ifelse(grepl("squam", d$bailey, ignore.case=TRUE), "Basal","Classical"), levels=c("Basal","Classical"))
  else g <- factor(ifelse(grepl("class", d[[var]], ignore.case=TRUE), "Classical","Basal"), levels=c("Basal","Classical"))
  d$g <- g
  if (nlevels(droplevels(d$g)) < 2 || sum(d$death) < 5 || min(table(d$g)) < 3) return(NULL)
  s <- summary(coxph(Surv(os,death)~g, d))
  data.frame(hr=s$conf.int[1,1], lo=s$conf.int[1,3], hi=s$conf.int[1,4], p=s$coef[1,5], n_basal=sum(d$g=="Basal"))
}
axis <- list()
for (sp in list(c("Linehan","purist"),c("Rash","purist"),c("Accept","purist"),
                c("Rash","moffitt"),c("Accept","moffitt"),c("Rash","bailey"),c("Accept","bailey"))){
  r <- classHR(D[D$cohort==sp[1],], sp[2]); if (!is.null(r)){ r$cohort<-sp[1]; r$classifier<-c(purist="PurIST",moffitt="Moffitt",bailey="Bailey")[sp[2]]; axis[[length(axis)+1]]<-r } }
axis <- do.call(rbind, axis)
# D1 separates PurIST (AUC), per cohort
d1auc <- sapply(COHS, function(c0){ d<-D[D$cohort==c0 & !is.na(D$purist), ]; y<-d$purist=="Classical"
  if (length(unique(y))<2) return(NA); r<-rank(d$D1); (mean(r[y])-(sum(y)+1)/2)/sum(!y) })

## ---- (3) subtype prevalence per cohort ----
prev <- do.call(rbind, lapply(COHS, function(c0){ d<-D[D$cohort==c0, ]
  data.frame(cohort=c0, n=nrow(d), events=sum(d$death),
             n_basal=sum(d$purist=="Basal-like",na.rm=TRUE),
             pct_basal=round(100*mean(d$purist=="Basal-like",na.rm=TRUE),0)) }))

## ---- (4) PFS (RASH-ACCEPT) ----
exR <- canon$RASH_ACCEPT$extras; ZR <- cbind(D1=projZ(RA$xr,PGENES,1),D2=projZ(RA$xr,PGENES,2),D3=projZ(RA$xr,PGENES,3))
pfs <- list()
for (z in c("D1","D2","D3")){ dp<-data.frame(score=as.numeric(scale(ZR[,z])),
    pfs=exR$ttp_months[match(rownames(ZR),exR$sample)], ev=exR$ttp_event[match(rownames(ZR),exR$sample)])
  dp<-dp[is.finite(dp$pfs)&is.finite(dp$ev)&dp$pfs>0,]; if (nrow(dp)<10||sum(dp$ev)<5) next
  s<-summary(coxph(Surv(pfs,ev)~score,dp)); pfs[[z]]<-data.frame(program=z,hr=s$conf.int[1,1],lo=s$conf.int[1,3],hi=s$conf.int[1,4],p=s$coef[1,5],n=nrow(dp),ev=sum(dp$ev)) }
pfs <- do.call(rbind, pfs)

## ---- (5) D2 treatment-induced (Linehan paired pre/post) ----
exL <- canon$Linehan$extras; ZL <- projZ(L$xr,PGENES,2); clL<-canon$Linehan$clin
dl <- data.frame(sample=names(ZL), D2=as.numeric(scale(ZL)), pid=clL$patient_id[match(names(ZL),clL$sample)],
                 pp=exL$pre_post[match(names(ZL),exL$sample)])
pr <- merge(dl[dl$pp==1,], dl[dl$pp==2,], by="pid", suffixes=c("_pre","_post"))
paired <- list(n=nrow(pr), wilcox_p=wilcox.test(pr$D2_pre,pr$D2_post,paired=TRUE)$p.value,
               mean_delta=mean(pr$D2_post-pr$D2_pre),
               # paired pre/post D2 program scores for Fig S13B (derived z-scores, no raw
               # expression and no patient identifiers; rows are paired within patient)
               points=data.frame(D2_pre=pr$D2_pre, D2_post=pr$D2_post))

## ---- (6) D2 sub-component decomposition (cohort-stratified; exploratory) ----
GRP <- list(
  ECMcore   = c("COL12A1","COL8A1","COL4A4","COL11A1","COL10A1","COL5A1","COL5A2","COL6A3","PDGFRA","DDR2","INHBA","ADAMTS12","CDH11","FBN1","FGF7","LAMA2","LAMA4","WISP1","SVEP1","PAPPA","ITGBL1","POSTN"),
  Macrophage= c("CD163","MSR1","CYBB","F13A1","AOAH","FCGR2A","FCGR3A","CCL5","OSMR","MRC1","C1QA","C1QB"),
  Neural    = c("SLIT2","SEMA5A","PLXNC1","ZEB2","NFIB","SEMA3C","NRP1","NRP2","ROBO1"))
addg <- function(M){ df<-M$df; for(g in names(GRP)) df[[g]]<-projZ(M$xr,intersect(GRP[[g]],PGENES),2)[df$sample]; df }
DD <- rbind(addg(L), addg(RA)); DD <- DD[!is.na(DD$os)&!is.na(DD$death)&DD$os>0,]
for (c0 in unique(DD$cohort)) for (g in names(GRP)) { i<-DD$cohort==c0; DD[i,g]<-as.numeric(scale(DD[i,g])) }
da <- DD[complete.cases(DD[,c("purist","DeCAF")]),]; da$purist<-factor(da$purist); da$DeCAF<-factor(da$DeCAF)
decomp <- list()
for (g in names(GRP)){ s<-summary(coxph(as.formula(paste0("Surv(os,death)~",g,"+purist+DeCAF+strata(cohort)")),da))
  decomp[[length(decomp)+1]]<-data.frame(component=g,model="alone",hr=s$conf.int[g,1],lo=s$conf.int[g,3],hi=s$conf.int[g,4],p=s$coef[g,5]) }
sj<-summary(coxph(Surv(os,death)~ECMcore+Macrophage+Neural+purist+DeCAF+strata(cohort),da))
for (g in names(GRP)) decomp[[length(decomp)+1]]<-data.frame(component=g,model="joint",hr=sj$conf.int[g,1],lo=sj$conf.int[g,3],hi=sj$conf.int[g,4],p=sj$coef[g,5])
decomp <- do.call(rbind, decomp); decomp$n <- nrow(da); decomp$ev <- sum(da$death)

## ---- (7) D2 random-effects meta (secondary sensitivity) + leave-one-out ----
dl_meta <- function(b,se){ w<-1/se^2; mu0<-sum(w*b)/sum(w); Q<-sum(w*(b-mu0)^2); dfr<-length(b)-1
  tau2<-max(0,(Q-dfr)/(sum(w)-sum(w^2)/sum(w))); ws<-1/(se^2+tau2); mu<-sum(ws*b)/sum(ws); se_re<-sqrt(1/sum(ws))
  list(hr=exp(mu),lo=exp(mu-1.96*se_re),hi=exp(mu+1.96*se_re),I2=max(0,(Q-dfr)/Q)*100,p=2*pnorm(-abs(mu/se_re))) }
d2adj <- do.call(rbind, lapply(COHS, function(c0){ r<-fmtcox(D[D$cohort==c0,],"D2",TRUE); data.frame(cohort=c0,b=log(r$hr),se=(log(r$hi)-log(r$lo))/(2*1.96)) }))
meta_all <- dl_meta(d2adj$b, d2adj$se)
meta_loo <- lapply(COHS, function(c0){ k<-d2adj[d2adj$cohort!=c0,]; m<-dl_meta(k$b,k$se); data.frame(dropped=c0,hr=m$hr,lo=m$lo,hi=m$hi,p=m$p) })
meta_loo <- do.call(rbind, meta_loo)
meta <- list(all=data.frame(hr=meta_all$hr,lo=meta_all$lo,hi=meta_all$hi,I2=meta_all$I2,p=meta_all$p), loo=meta_loo)

## ---- (8) D1 vs GATA6 ISH orthogonal validation (COMPASS) ----
xrC <- rankX(t(canon$COMPASS_50$X)); D1c <- as.numeric(scale(projZ(xrC,PGENES,1)))
g6 <- suppressWarnings(as.numeric(as.character(canon$COMPASS_50$extras$gata6_ish[match(colnames(xrC),canon$COMPASS_50$extras$sample)])))
ok <- !is.na(g6)
gata6 <- list(rho = cor(D1c[ok], g6[ok], method="spearman"),
              p   = cor.test(D1c[ok], g6[ok], method="spearman")$p.value,
              n   = sum(ok),
              group_means = sapply(sort(unique(g6[ok])), function(lv) mean(D1c[ok][g6[ok]==lv])),
              # raw per-sample points for Fig 2E (scaled D1 score vs GATA6 ISH level);
              # aggregate-level derived values (no raw expression), so tracked in the .rds
              points = data.frame(D1 = D1c[ok], gata6 = g6[ok]))

## ---- (9) concordance: full-W vs top-270 projection (robustness) ----
.ps <- function(M, genes, j) projZ(M$xr, genes, j)
concordance <- sapply(1:3, function(j){
  a <- c(.ps(L, Wg, j), .ps(RA, Wg, j)); b <- c(.ps(L, PGENES, j), .ps(RA, PGENES, j))
  cor(a, b, use = "complete.obs") })
names(concordance) <- c("D1","D2","D3")

## ---- (10) O'Kane/COMPASS treated-metastatic D1 transportability (per ACTUAL arm) ----
# Frozen-projection transportability test in an independently treated metastatic cohort.
# Specimens are laser-capture microdissected (epithelial-enriched), so this tests the
# tumor-associated D1 direction, not the stromal programs. Treatment arms are kept
# SEPARATE (FFX, GA, GA/experimental) and never pooled (per-arm HRs + arm-stratified
# common effect with treatment-specific baseline hazards + factor x arm interaction).
okane <- local({
  ok <- canon$OKane; cl <- ok$clin
  xrk <- rankX(t(ok$X))
  Z1  <- projZ(xrk, PGENES, 1)                 # named by sample (colnames of xrk)
  keep <- cl$treatment %in% c("FFX", "GA", "GA/experimental") &
          is.finite(cl$os_months) & cl$os_months > 0 & !is.na(cl$death)
  # align the projected score to the clinical rows by SAMPLE NAME, then z-score
  # (do not rely on incidental row/column ordering matching between X and clin)
  dd <- data.frame(os = cl$os_months, ev = cl$death,
                   arm = factor(cl$treatment, levels = c("FFX", "GA", "GA/experimental")),
                   D1 = as.numeric(scale(Z1[cl$sample])))[keep, ]
  # D1 is already z-scaled over the cohort; use it directly so per-arm and
  # stratified HRs are all "per cohort SD" on a common scale (no per-arm rescaling).
  per_arm <- do.call(rbind, lapply(levels(dd$arm), function(a) {
    z <- dd[dd$arm == a, ]; s <- summary(coxph(Surv(os, ev) ~ D1, z))
    data.frame(arm = a, n = nrow(z), events = sum(z$ev),
               hr = exp(s$coef[1]), lo = s$conf.int[1, 3], hi = s$conf.int[1, 4])
  }))
  ss <- summary(coxph(Surv(os, ev) ~ D1 + strata(arm), dd))
  int_p <- anova(coxph(Surv(os, ev) ~ D1 + arm, dd),
                 coxph(Surv(os, ev) ~ D1 * arm, dd))[2, "Pr(>|Chi|)"]
  list(per_arm = per_arm,
       stratified = c(hr = exp(ss$coef[1]), lo = ss$conf.int[1, 3],
                      hi = ss$conf.int[1, 4], p = ss$coef[1, 5]),
       interaction_p = int_p, n = nrow(dd), events = sum(dd$ev),
       n_shared_genes = nrow(xrk))
})
cat(sprintf("OKane D1: arm-stratified HR=%.2f (%.2f-%.2f) P=%.3g | interaction P=%.2f | n=%d ev=%d\n",
    okane$stratified["hr"], okane$stratified["lo"], okane$stratified["hi"],
    okane$stratified["p"], okane$interaction_p, okane$n, okane$events))

## ---- assemble + save ----
treated_cohort_stats <- list(
  os = os, axis = axis, d1_auc = d1auc, prevalence = prev, pfs = pfs,
  paired = paired, decomp = decomp, meta = meta, gata6 = gata6, okane = okane,
  concordance = concordance, ntop = 270,
  cohorts = list(Linehan = "borderline-resectable/locally-advanced; FOLFIRINOX +/- CCR2 inhibitor (PF-04136309)",
                 Rash = "metastatic; gemcitabine + erlotinib", Accept = "metastatic; gemcitabine +/- afatinib",
                 OKane = "O'Kane/COMPASS metastatic (laser-capture microdissected); FOLFIRINOX or gemcitabine/nab-paclitaxel"),
  generated_note = "Derived from restricted-access external cohorts (O'Kane/COMPASS, RASH-ACCEPT, Linehan) via code/11_treated_cohort_analysis.R; only aggregate-level and paired derived program scores (no raw expression or patient identifiers) are stored; raw data not in repo.")
saveRDS(treated_cohort_stats, file.path("results", "treated_cohort_stats.rds"))
cat("Wrote results/treated_cohort_stats.rds\n")
print(os[, c("cohort","program","model","hr","lo","hi","p")])
cat("\nclassifier axis:\n"); print(axis[,c("cohort","classifier","hr","lo","hi","p")])
cat("\nD2 meta all:", sprintf("%.2f (%.2f-%.2f) I2=%.0f%%", meta$all$hr, meta$all$lo, meta$all$hi, meta$all$I2), "\n")
cat("GATA6 rho:", round(gata6$rho,2), "p:", signif(gata6$p,2), "n:", gata6$n, "\n")
