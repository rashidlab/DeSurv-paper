#!/usr/bin/env Rscript
# 14_d1_variance_partition.R
# ---------------------------------------------------------------------------
# D1 variance partition (Nature Cancer revision) -- defuses the central tension:
# "is D1 just a re-derivation of the PurIST (Classical/Basal) and DeCAF
#  (restCAF/proCAF) classifiers, or does it carry novel, coupled signal?"
#
# Two complementary partitions over the five pooled external validation cohorts:
#   (A) VARIANCE: decompose the D1 score into the fraction explained by PurIST,
#       by DeCAF, by their interaction, and the residual (novel) fraction.
#       Because D1 is a *coupling* (Classical tumour + restCAF stroma), the
#       informative quantity is how much of D1 NEITHER classifier explains.
#   (B) PROGNOSIS: nested stratified Cox -- does D1 improve survival prediction
#       beyond PurIST + DeCAF (partial LR test, adjusted HR, C-index gain)?
# In-repo only (val_latent_desurv + data_val sampInfo). No new data.
# ---------------------------------------------------------------------------

suppressMessages({ library(survival) })
setwd("/home/naimrashid/Downloads/DeSurv-paper-clean")

val  <- readRDS("results/val_latent_desurv_tcgacptac.rds")
dval <- readRDS("results/data_val_filtered_tcgacptac.rds")

## --- assemble pooled per-sample table: D1/D2/D3 + survival + PurIST + DeCAF ---
rows <- list()
for (e in val) {
  Z <- as.matrix(e$latent); surv <- e$survival
  ds_match <- Filter(function(d) d$sampInfo$dataset[1] == e$dataset, dval)
  if (!length(ds_match)) next
  si <- ds_match[[1]]$sampInfo; si$sample_id <- rownames(si)
  ids <- rownames(surv)
  pur <- si$PurIST[match(ids, si$sample_id)]
  dec <- si$DeCAF [match(ids, si$sample_id)]
  dec[dec == "permCAF"] <- "proCAF"
  rows[[e$dataset]] <- data.frame(
    sample_id = ids, dataset = e$dataset,
    D1 = Z[,1], D2 = Z[,2], D3 = Z[,3],
    time = surv$time, event = surv$event,
    PurIST = pur, DeCAF = dec, stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows); rownames(df) <- NULL
df <- df[is.finite(df$time) & is.finite(df$D1) & !is.na(df$PurIST) & !is.na(df$DeCAF), ]
# z-score D1 within cohort (remove scale/level differences between datasets)
df$D1z <- ave(df$D1, df$dataset, FUN = function(x) as.numeric(scale(x)))
df$PurIST <- factor(df$PurIST); df$DeCAF <- factor(df$DeCAF)
cat(sprintf("pooled n = %d across %d cohorts; PurIST levels {%s}; DeCAF levels {%s}\n",
            nrow(df), length(unique(df$dataset)),
            paste(levels(df$PurIST),collapse=","), paste(levels(df$DeCAF),collapse=",")))

## ===========================================================================
## (A) Variance partition of the D1 score
##     Baseline removes between-cohort variance; classifiers explain the rest.
## ===========================================================================
r2 <- function(form) summary(lm(form, data = df))$r.squared
R2_base <- r2(D1z ~ dataset)                                   # cohort only
R2_pur  <- r2(D1z ~ dataset + PurIST)
R2_dec  <- r2(D1z ~ dataset + DeCAF)
R2_pd   <- r2(D1z ~ dataset + PurIST + DeCAF)
R2_full <- r2(D1z ~ dataset + PurIST * DeCAF)
# variance of D1 (beyond cohort) attributable to each component
denom <- 1 - R2_base                                           # within-cohort variance
sp_pur <- (R2_pd  - R2_dec)  / denom                           # semipartial: PurIST | DeCAF
sp_dec <- (R2_pd  - R2_pur)  / denom                           # semipartial: DeCAF | PurIST
sp_int <- (R2_full - R2_pd)  / denom                           # interaction
joint  <- (R2_pd  - R2_base) / denom                           # PurIST+DeCAF together (additive)
explained_total <- (R2_full - R2_base) / denom                 # all classifier terms
residual_novel  <- 1 - explained_total                         # D1 variance NEITHER explains

## ===========================================================================
## (B) Prognostic partition: does D1 add beyond PurIST + DeCAF?
## ===========================================================================
m0 <- coxph(Surv(time,event) ~ PurIST + DeCAF + strata(dataset), data = df)
m1 <- coxph(Surv(time,event) ~ D1z + PurIST + DeCAF + strata(dataset), data = df)
mm <- coxph(Surv(time,event) ~ D1z + strata(dataset), data = df)            # D1 marginal
lr <- anova(m0, m1)                                            # partial LR test for D1
lr_chisq <- lr$Chisq[2]; lr_p <- lr$`Pr(>|Chi|)`[2]
s1 <- summary(m1)
hr_adj <- s1$conf.int["D1z", c("exp(coef)","lower .95","upper .95")]
p_adj  <- s1$coefficients["D1z","Pr(>|z|)"]
sm <- summary(mm); hr_marg <- sm$conf.int["D1z", c("exp(coef)","lower .95","upper .95")]
c_m0 <- summary(m0)$concordance["C"]; c_m1 <- summary(m1)$concordance["C"]

res <- list(
  meta = list(n = nrow(df), n_events = sum(df$event), n_cohorts = length(unique(df$dataset))),
  variance = list(R2_base=R2_base, R2_purist=R2_pur, R2_decaf=R2_dec, R2_pd=R2_pd, R2_full=R2_full,
                  semipartial_purist=sp_pur, semipartial_decaf=sp_dec, interaction=sp_int,
                  joint_classifiers=joint, explained_total=explained_total,
                  residual_novel=residual_novel),
  prognosis = list(d1_adj_hr=unname(hr_adj[1]), d1_adj_lo=unname(hr_adj[2]), d1_adj_hi=unname(hr_adj[3]),
                   d1_adj_p=p_adj, d1_marg_hr=unname(hr_marg[1]), d1_marg_lo=unname(hr_marg[2]),
                   d1_marg_hi=unname(hr_marg[3]), lr_chisq=lr_chisq, lr_p=lr_p,
                   c_base=unname(c_m0), c_with_d1=unname(c_m1), c_gain=unname(c_m1 - c_m0)))
saveRDS(res, "results/d1_variance_partition_stats.rds")

fp <- function(p) ifelse(p<1e-4,"< 0.0001",sprintf("%.4f",p))
cat("\n================= D1 VARIANCE PARTITION =================\n")
cat(sprintf("pooled n=%d (%d events), %d cohorts\n\n", res$meta$n, res$meta$n_events, res$meta$n_cohorts))
cat("(A) Variance of D1 (within-cohort) explained by the classifiers:\n")
cat(sprintf("    PurIST alone (| DeCAF)  : %5.1f%%\n", 100*sp_pur))
cat(sprintf("    DeCAF  alone (| PurIST) : %5.1f%%\n", 100*sp_dec))
cat(sprintf("    PurIST x DeCAF interact : %5.1f%%\n", 100*sp_int))
cat(sprintf("    both classifiers (total): %5.1f%%\n", 100*explained_total))
cat(sprintf("    RESIDUAL / novel        : %5.1f%%  <- D1 variance neither classifier explains\n\n", 100*residual_novel))
cat("(B) Prognostic value of D1 beyond PurIST + DeCAF (stratified Cox):\n")
cat(sprintf("    D1 marginal HR per SD     : %.2f (%.2f-%.2f)\n", hr_marg[1],hr_marg[2],hr_marg[3]))
cat(sprintf("    D1 adjusted HR per SD     : %.2f (%.2f-%.2f), p %s\n", hr_adj[1],hr_adj[2],hr_adj[3], fp(p_adj)))
cat(sprintf("    partial LR chi-sq (1 df)  : %.1f, p %s\n", lr_chisq, fp(lr_p)))
cat(sprintf("    C-index: PurIST+DeCAF %.3f -> +D1 %.3f  (gain %+.3f)\n", c_m0, c_m1, c_m1-c_m0))
cat("\nSaved -> results/d1_variance_partition_stats.rds\n")
