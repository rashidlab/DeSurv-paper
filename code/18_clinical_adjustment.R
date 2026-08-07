#!/usr/bin/env Rscript
# 18_clinical_adjustment.R
# ---------------------------------------------------------------------------
# Reviewer request: does the DeSurv risk score (and the D1 axis) remain
# prognostic after adjustment for standard clinical prognostic variables
# (stage, grade, age, sex, nodal status, resection margin) and for tumor
# purity, per cohort?
#
# Scores: validation cohorts use the exact latent/LP from val_latent (full-W
# projection, same object the manuscript reports); training cohorts (TCGA,
# CPTAC) are scored by the identical projection Z = X^T W, LP = Z beta.
# Clinical covariates are joined from the original per-cohort objects
# (data/original/<cohort>.rds $sampInfo, positionally aligned to $ex columns).
# Covariate encodings were inspected per cohort; unparseable/absent covariates
# are simply dropped from that cohort's model. Puleo "developmental stage"
# (all "adult") and "ffpeblock age" (FFPE block age, not patient age) are NOT
# used.
# ---------------------------------------------------------------------------

suppressMessages({ library(survival) })
set.seed(1)

fit  <- readRDS("results/tar_fit_desurv_tcgacptac.rds")
W    <- fit$W; beta <- as.numeric(fit$beta)

## ---- 1. DeSurv scores per cohort ------------------------------------------
score_list <- list()
# validation cohorts: exact scores from val_latent
vv <- readRDS("results/val_latent_desurv_tcgacptac.rds")
for (e in vv) score_list[[e$dataset]] <- data.frame(
  id = rownames(e$latent), risk = as.numeric(e$risk_score),
  D1 = as.numeric(e$latent[, 1]), time = e$survival$time, event = e$survival$event,
  dataset = e$dataset, stringsAsFactors = FALSE)
# training cohorts: mirror extract_val_latent projection (full W)
tr <- readRDS("results/tar_data_filtered_tcgacptac.rds")
g  <- intersect(rownames(W), rownames(tr$ex))
Ztr  <- t(tr$ex[g, , drop = FALSE]) %*% W[g, , drop = FALSE]
lptr <- as.numeric(Ztr %*% beta)
si0  <- tr$sampInfo
for (d in unique(si0$dataset)) { idx <- which(si0$dataset == d)
  score_list[[d]] <- data.frame(id = colnames(tr$ex)[idx], risk = lptr[idx],
    D1 = Ztr[idx, 1], time = si0$time[idx], event = si0$event[idx],
    dataset = d, stringsAsFactors = FALSE) }

## ---- 2. original clinical (sampInfo aligned to ex columns) -----------------
get_orig <- function(cohort) {
  o <- readRDS(paste0("data/original/", cohort, ".rds"))
  si <- o$sampInfo; si$`.__exid` <- colnames(o$ex); si
}

## ---- 3. harmonizers -------------------------------------------------------
ord_stage <- function(x) { x <- toupper(trimws(as.character(x)))
  x <- sub("^STAGE\\s*", "", x)
  ifelse(grepl("^IV", x), 4L, ifelse(grepl("^III", x), 3L,
    ifelse(grepl("^II", x), 2L, ifelse(grepl("^I", x), 1L, NA_integer_)))) }
ord_grade <- function(x) { x <- tolower(trimws(as.character(x)))
  ifelse(grepl("^1|well", x), 1L, ifelse(grepl("^2|moder", x), 2L,
    ifelse(grepl("^3|poor", x), 3L, ifelse(grepl("^4|undiff", x), 3L, NA_integer_)))) }
bin_sex    <- function(x) { x <- tolower(trimws(as.character(x)))
  factor(ifelse(grepl("^m", x), "M", ifelse(grepl("^f", x), "F", NA_character_)), levels = c("F","M")) }
node_pos   <- function(x) { x <- toupper(trimws(as.character(x)))
  ifelse(grepl("N0", x), 0L, ifelse(grepl("N1|N2|N3", x), 1L, NA_integer_)) }
margin_pos <- function(x) { x <- tolower(trimws(as.character(x)))
  ifelse(grepl("r1|positive", x), 1L, ifelse(grepl("r0|negative", x), 0L, NA_integer_)) }
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# per-cohort covariate spec: name -> (column, transform)
SPEC <- list(
  TCGA_PAAD = list(stage=c("AJCC.pathologic.tumor.stage",  "ord_stage"),
                   grade=c("Grade","ord_grade"), age=c("Age.at.initial.pathologic.diagnosis","num"),
                   sex=c("Gender","bin_sex"), nodal=c("Pathology.N.stage","node_pos"),
                   purity=c("ABSOLUTE.Purity","num")),
  CPTAC     = list(stage=c("tumor_stage_pathological","ord_stage"), age=c("age","num"),
                   sex=c("sex","bin_sex"), nodal=c("pathologic_staging_regional_lymph_nodes_pn","node_pos")),
  PACA_AU_array = list(stage=c("AJCC.Pathology.Stage","ord_stage"), grade=c("Tumour.Grade","ord_grade"),
                   age=c("Age.at.Diagnosis.in.Years","num"), sex=c("Gender","bin_sex"),
                   purity=c("TumorPurity","num")),
  PACA_AU_seq   = list(stage=c("AJCC.Pathology.Stage","ord_stage"), grade=c("Tumour.Grade","ord_grade"),
                   age=c("Age.at.Diagnosis.in.Years","num"), sex=c("Gender","bin_sex"),
                   purity=c("TumorPurity","num")),
  Dijk      = list(grade=c("differentiation grade","ord_grade"), age=c("age (year)","num"),
                   sex=c("sex","bin_sex"), nodal=c("lymph node metastasis ratio","num")),
  Moffitt_GEO_array = list(sex=c("sex","bin_sex"), margin=c("Margin","margin_pos"),
                   purity=c("purity","num")),
  Puleo_array = list(sex=c("Characteristics.sex.","bin_sex"),
                   margin=c("Characteristics.resection.margin.","margin_pos")))

## ---- 4. per-cohort Cox models ---------------------------------------------
tf <- function(fun, v) get(fun)(v)
fit_hr <- function(df, extra) {   # adjusted HR per SD of DeSurv risk score
  df$rz <- as.numeric(scale(df$risk))
  keep <- c("time","event","rz", extra)
  d <- df[stats::complete.cases(df[, keep, drop=FALSE]), keep, drop=FALSE]
  # drop covariates that are constant after NA filtering
  extra2 <- extra[vapply(extra, function(k) length(unique(d[[k]])) > 1, logical(1))]
  if (nrow(d) < 15 || sum(d$event) < 8) return(NULL)
  f <- as.formula(paste("Surv(time,event) ~ rz", if(length(extra2)) paste("+", paste(extra2, collapse="+")) else ""))
  m <- tryCatch(coxph(f, data=d), error=function(e) NULL); if (is.null(m)) return(NULL)
  s <- summary(m)$coefficients["rz", , drop=FALSE]
  ci <- summary(m)$conf.int["rz", , drop=FALSE]
  list(n=nrow(d), events=sum(d$event), covars=extra2,
       hr=ci[,"exp(coef)"], lo=ci[,"lower .95"], hi=ci[,"upper .95"], p=s[,"Pr(>|z|)"])
}

rows <- list()
for (ds in names(score_list)) {
  sc <- score_list[[ds]]
  orig <- tryCatch(get_orig(orig_file <- ds), error=function(e) NULL)
  # build harmonized covariate frame joined by ex id
  clin <- data.frame(id = sc$id, stringsAsFactors=FALSE)
  spec <- SPEC[[ds]]
  matched <- 0
  if (!is.null(orig) && !is.null(spec)) {
    key <- orig$`.__exid`
    mi  <- match(sub("_seq$", "", sc$id), key)   # PACA_AU_seq score IDs carry a _seq suffix
    matched <- sum(!is.na(mi))
    for (nm in names(spec)) { col <- spec[[nm]][1]; fun <- spec[[nm]][2]
      raw <- if (col %in% names(orig)) orig[[col]][mi] else rep(NA, nrow(sc))
      clin[[nm]] <- tf(fun, raw) }
  }
  df <- cbind(sc, clin[,-1,drop=FALSE])
  avail <- names(spec)[vapply(names(spec), function(k) sum(!is.na(df[[k]])) >= 15, logical(1))]
  clin_covs <- setdiff(avail, "purity")
  m_unadj <- fit_hr(df, character(0))
  m_full  <- if (length(clin_covs)) fit_hr(df, clin_covs) else NULL
  m_pur   <- if ("purity" %in% avail) fit_hr(df, "purity") else NULL
  # Correlation of the scores themselves with tumor purity. The HR comparison
  # above shows the ASSOCIATION is not driven by purity; this shows how far the
  # SCORE tracks purity in the first place, which is the more direct question.
  m_purcor <- NULL
  if ("purity" %in% avail) {
    ok <- is.finite(df$purity)
    m_purcor <- list(
      n     = sum(ok),
      D1    = suppressWarnings(cor(df$D1[ok],   df$purity[ok], method = "spearman")),
      D1_p  = suppressWarnings(cor.test(df$D1[ok], df$purity[ok], method = "spearman")$p.value),
      risk  = suppressWarnings(cor(df$risk[ok], df$purity[ok], method = "spearman")))
  }
  rows[[ds]] <- list(dataset=ds, n=nrow(sc), matched_clinical=matched,
                     covariates_available=avail, unadjusted=m_unadj,
                     clinical_adjusted=m_full, purity_adjusted=m_pur,
                     purity_cor=m_purcor)
}

saveRDS(rows, "results/desurv_clinical_adjustment.rds")

## ---- 5. print summary -----------------------------------------------------
fmt <- function(m) if (is.null(m)) "  --" else sprintf("HR=%.2f (%.2f-%.2f) p=%.1e [n=%d,e=%d]", m$hr,m$lo,m$hi,m$p,m$n,m$events)
cat("\n=== DeSurv risk score: adjusted HR per SD ===\n")
for (r in rows) {
  cat(sprintf("\n%-18s matched-clinical=%d  covars={%s}\n", r$dataset, r$matched_clinical, paste(r$covariates_available, collapse=",")))
  cat("   unadjusted        :", fmt(r$unadjusted), "\n")
  cat("   clinical-adjusted :", fmt(r$clinical_adjusted), if(!is.null(r$clinical_adjusted)) paste0("  {",paste(r$clinical_adjusted$covars,collapse=","),"}") else "", "\n")
  cat("   purity-adjusted   :", fmt(r$purity_adjusted), "\n")
}
cat("\nSaved -> results/desurv_clinical_adjustment.rds\n")
