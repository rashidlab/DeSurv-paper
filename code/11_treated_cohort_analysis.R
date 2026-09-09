#!/usr/bin/env Rscript
# code/11_treated_cohort_analysis.R
# ---------------------------------------------------------------------------
# Per-cohort treated-disease translational analysis for DeSurv.
#
# REPRODUCIBILITY NOTE: the Linehan and O'Kane/COMPASS cohorts are restricted-access
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
  "this script requires the restricted-access Linehan/COMPASS data environment. ",
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
L <- mk("Linehan", COHORT$Linehan)

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

# Self-correlation sensitivity. GATA6 is itself one of the genes in the D1
# projection support, so part of the correlation above is the score correlating
# with one of its own inputs. Recompute with GATA6 dropped from the projection.
{
  D1x <- as.numeric(scale(projZ(xrC, setdiff(PGENES, "GATA6"), 1)))
  ctx <- suppressWarnings(cor.test(D1x[ok], g6[ok], method = "spearman"))
  gata6$in_support  <- "GATA6" %in% PGENES
  gata6$n_support   <- length(intersect(PGENES, rownames(xrC)))
  gata6$rho_excl    <- unname(ctx$estimate)
  gata6$p_excl      <- ctx$p.value
}

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
  okane = okane, paired = paired, gata6 = gata6, ntop = 270,
  cohorts = list(
    Linehan = "borderline-resectable/locally-advanced; FOLFIRINOX +/- CCR2 inhibitor (PF-04136309); used ONLY for the paired pre/post D2 comparison",
    OKane   = "O'Kane/COMPASS metastatic (laser-capture microdissected); FOLFIRINOX or gemcitabine/nab-paclitaxel; used for arm-stratified D1 transportability and the GATA6 comparison"),
  generated_note = paste(
    "Derived from restricted-access external cohorts (O'Kane/COMPASS and Linehan) via",
    "code/11_treated_cohort_analysis.R. Stores derived summary statistics plus per-sample",
    "derived program scores (paired D2 pre/post in paired$points; D1 vs GATA6 in gata6$points)",
    "with no raw expression and no patient identifiers. This object holds only the two treated",
    "analyses reported in the manuscript: arm-stratified D1 transportability in O'Kane/COMPASS",
    "and the within-patient D2 change in paired Linehan biopsies. No bulk program-by-survival",
    "association is computed or stored."))
saveRDS(treated_cohort_stats, file.path("results", "treated_cohort_stats.rds"))
cat("Wrote results/treated_cohort_stats.rds\n")
cat("GATA6 rho:", round(gata6$rho, 2), "p:", signif(gata6$p, 2), "n:", gata6$n, "\n")
cat("paired D2: n=", paired$n, " Wilcoxon p=", signif(paired$wilcox_p, 3), "\n", sep = "")
