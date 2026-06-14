#!/usr/bin/env Rscript
# 12_single_cell_validation.R
# ---------------------------------------------------------------------------
# Tier-1 single-cell validation of the DeSurv programs (Nature Cancer revision).
#
# Question (the load-bearing one): are the DeSurv factors *real cell biology*,
# not just a survival-optimal direction in bulk? We project the trained,
# frozen W onto the Elyada (2019) PDAC scRNA atlas -- WITHOUT retraining and
# WITHOUT any survival objective in play -- and ask:
#   (1) State mapping. Do D1/D2/D3 localize to the expected cell states?
#         D1 (Classical/restCAF axis) -> Classical-malignant + restCAF-leaning CAFs
#         D2 (proCAF stroma)          -> proCAF-leaning CAFs (iCAF-like)
#         D3 (Basal-like)             -> Basal-like malignant cells
#   (2) Compartment concordance. Does D2 track the DeCAF proCAF signature and
#         D1's stromal arm track the restCAF signature, at single-cell resolution?
#   (3) Co-occurrence within patients. Do Classical-malignant cells and
#         restCAF stroma co-occur in the SAME patients -- the tumour-stroma
#         coupling that defines D1?
#
# Scoring convention is identical to the manuscript validation (code/05, code/11):
#   within-cell rank over the shared trained-W gene namespace, then project on
#   the top-270-per-factor union W-tilde (PGENES). No Seurat dependency: the
#   Elyada objects are read as serialized S4 and their dgCMatrix / data.frame
#   slots are pulled via attributes (R 4.6 cannot load the 4.5-built Seurat).
# ---------------------------------------------------------------------------

suppressMessages({ library(Matrix); library(survival) })
set.seed(1)

ROOT   <- "/home/naimrashid/Downloads/DeSurv-paper-clean"
SC_DIR <- "/home/naimrashid/Downloads/DeSurv-paper/data/derv"
setwd(ROOT)

## --- trained DeSurv fit + top-270 projection basis (manuscript convention) ---
fit <- readRDS(file.path("results", "tar_fit_desurv_tcgacptac.rds"))
W <- fit$W; Wg <- rownames(W)                      # factors = columns 1/2/3 = D1/D2/D3
source("R/get_top_genes.R")
PGENES <- unique(unlist(lapply(get_top_genes(fit$W, 270)$top_genes, as.character)))

clean <- function(ex){ rownames(ex) <- gsub("^X(?=[0-9_])", "", rownames(ex), perl = TRUE); ex }
# rank each cell (column) over the shared W namespace, exactly as code/11 rankX()
rankX <- function(ex){ ex <- clean(ex); sh <- intersect(Wg, rownames(ex))
  xr <- apply(as.matrix(ex[sh, , drop = FALSE]), 2, rank, ties.method = "average")
  rownames(xr) <- sh; xr }
projZ <- function(xr, genes, j){ g <- intersect(genes, rownames(xr))
  if (!length(g)) return(setNames(rep(NA_real_, ncol(xr)), colnames(xr)))
  setNames(drop(t(xr[g, , drop = FALSE]) %*% W[g, j, drop = FALSE]), colnames(xr)) }

## --- read a Seurat .Rds without Seurat: slots are attributes -----------------
read_seurat <- function(path){
  o  <- readRDS(path)
  A  <- attr(o, "assays")
  md <- attr(o, "meta.data")
  dat <- attr(A$RNA, "data")                       # log-normalized RNA, genes x cells
  list(data = dat, meta = md)
}

## --- DeCAF single-cell CAF reference signatures (from inst/single_cell.R) -----
proCAF_sig  <- c("IGFL2","NOX4","VSNL1","BICD1","NPR3","ETV1","ITGA11","CNIH3","COL11A1")
restCAF_sig <- c("CHRDL1","OGN","PI16","ANK2","ABCA8","TGFBR3","FBLN5","SCARA5","KIAA1217")
# mean log-norm expression of a signature, z-scored across cells
sig_score <- function(dat, sig){ g <- intersect(sig, rownames(dat))
  if (length(g) < 3) return(rep(NA_real_, ncol(dat)))
  s <- colMeans(as.matrix(dat[g, , drop = FALSE])); as.numeric(scale(s)) }

zc <- function(v) as.numeric(scale(v))             # z across cells
fmtp <- function(p) ifelse(p < 1e-4, "< 0.0001", sprintf("%.4f", p))

cat("== loading Elyada CAF + tumour objects (no Seurat) ==\n")
caf <- read_seurat(file.path(SC_DIR, "Elyada_caf_umap.Rds"))
tum <- read_seurat(file.path(SC_DIR, "Elyada_PDAC_umap.Rds"))

## ===========================================================================
## (1) CAF compartment: D2 = proCAF program
## ===========================================================================
cxr  <- rankX(caf$data)
cD1  <- zc(projZ(cxr, PGENES, 1)); cD2 <- zc(projZ(cxr, PGENES, 2)); cD3 <- zc(projZ(cxr, PGENES, 3))
cpro <- sig_score(caf$data, proCAF_sig)
crest<- sig_score(caf$data, restCAF_sig)
clab <- as.character(caf$meta$label)               # iCAF / myCAF / apCAF
csamp<- as.character(caf$meta$sample)
ccond<- as.character(caf$meta$condition)

# D2 vs DeCAF proCAF / restCAF signatures (Spearman, cell level)
caf_concord <- list(
  D2_vs_proCAF  = cor.test(cD2, cpro,  method = "spearman"),
  D2_vs_restCAF = cor.test(cD2, crest, method = "spearman"),
  D1_vs_restCAF = cor.test(cD1, crest, method = "spearman"),
  D1_vs_proCAF  = cor.test(cD1, cpro,  method = "spearman"))

# D2 across Elyada CAF subtypes (proCAF is inflammatory/iCAF-leaning)
caf_by_label <- aggregate(cbind(D1=cD1, D2=cD2, D3=cD3, proCAF=cpro, restCAF=crest),
                          by = list(label = clab), FUN = function(x) mean(x, na.rm=TRUE))
kw_D2  <- kruskal.test(cD2 ~ factor(clab))
# proCAF-leaning (iCAF) vs restCAF-leaning (myCAF) contrast on D2
icaf <- clab == "iCAF"; mycaf <- clab == "myCAF"
wil_D2_iVm <- if (sum(icaf) > 2 && sum(mycaf) > 2)
  wilcox.test(cD2[icaf], cD2[mycaf]) else NULL

## ===========================================================================
## (2) Tumour compartment: D1 = Classical, D3 = Basal-like
## ===========================================================================
txr  <- rankX(tum$data)
tD1  <- zc(projZ(txr, PGENES, 1)); tD2 <- zc(projZ(txr, PGENES, 2)); tD3 <- zc(projZ(txr, PGENES, 3))
tlab <- as.character(tum$meta$label_broad)         # Classical / Basal-like
tsamp<- as.character(tum$meta$sample)
tcond<- as.character(tum$meta$condition)
tmalig <- as.character(tum$meta$malig)             # Malignant / Normal

tum_by_label <- aggregate(cbind(D1=tD1, D2=tD2, D3=tD3),
                          by = list(label = tlab), FUN = function(x) mean(x, na.rm=TRUE))
cl <- tlab == "Classical"; bl <- tlab == "Basal-like"
wil_D1_CvB <- if (sum(cl) > 2 && sum(bl) > 2) wilcox.test(tD1[cl], tD1[bl]) else NULL
wil_D3_CvB <- if (sum(cl) > 2 && sum(bl) > 2) wilcox.test(tD3[cl], tD3[bl]) else NULL

## ===========================================================================
## (3) Co-occurrence WITHIN patients (the D1 coupling test)
##     PDAC samples only; join tumour + CAF compartments by `sample`.
##     x-axis: fraction of malignant cells that are Classical (tumour obj)
##     y-axis: mean restCAF signal in that patient's CAFs (CAF obj)
## ===========================================================================
pdac_t <- tcond == "PDAC" & tmalig == "Malignant"
pdac_c <- ccond == "PDAC"
# per-patient Classical fraction among malignant cells
clf <- tapply(tlab[pdac_t] == "Classical", tsamp[pdac_t], mean)
# per-patient restCAF signal (DeCAF restCAF sig) and D1 in CAFs
restp <- tapply(crest[pdac_c], csamp[pdac_c], mean, na.rm = TRUE)
D2p   <- tapply(cD2[pdac_c],   csamp[pdac_c], mean, na.rm = TRUE)
D1cp  <- tapply(cD1[pdac_c],   csamp[pdac_c], mean, na.rm = TRUE)
common <- intersect(names(clf), names(restp))
cooc <- data.frame(sample = common,
                   classical_frac = as.numeric(clf[common]),
                   restCAF_caf    = as.numeric(restp[common]),
                   D1_caf         = as.numeric(D1cp[common]),
                   D2_caf         = as.numeric(D2p[common]))
cooc_test_rest <- if (nrow(cooc) >= 4)
  cor.test(cooc$classical_frac, cooc$restCAF_caf, method = "spearman") else NULL
cooc_test_D1 <- if (nrow(cooc) >= 4)
  cor.test(cooc$classical_frac, cooc$D1_caf, method = "spearman") else NULL

## ===========================================================================
## Save
## ===========================================================================
res <- list(
  meta = list(n_caf = ncol(caf$data), n_tum = ncol(tum$data),
              n_patients = length(unique(c(csamp, tsamp))),
              n_pdac_patients = nrow(cooc), ntop = 270,
              W_genes = length(Wg), pgenes = length(PGENES),
              proCAF_genes_found = length(intersect(proCAF_sig, rownames(caf$data))),
              restCAF_genes_found= length(intersect(restCAF_sig, rownames(caf$data)))),
  caf_concord  = lapply(caf_concord, function(t) c(rho = unname(t$estimate), p = t$p.value)),
  caf_by_label = caf_by_label,
  caf_D2_kruskal = c(stat = unname(kw_D2$statistic), p = kw_D2$p.value),
  caf_D2_iCAF_vs_myCAF = if (!is.null(wil_D2_iVm)) c(W = unname(wil_D2_iVm$statistic), p = wil_D2_iVm$p.value,
                          iCAF_mean = mean(cD2[icaf]), myCAF_mean = mean(cD2[mycaf])) else NA,
  tum_by_label = tum_by_label,
  tum_D1_Class_vs_Basal = if (!is.null(wil_D1_CvB)) c(W = unname(wil_D1_CvB$statistic), p = wil_D1_CvB$p.value,
                          Classical_mean = mean(tD1[cl]), Basal_mean = mean(tD1[bl])) else NA,
  tum_D3_Class_vs_Basal = if (!is.null(wil_D3_CvB)) c(W = unname(wil_D3_CvB$statistic), p = wil_D3_CvB$p.value,
                          Classical_mean = mean(tD3[cl]), Basal_mean = mean(tD3[bl])) else NA,
  cooccurrence = cooc,
  cooc_test_restCAF = if (!is.null(cooc_test_rest)) c(rho = unname(cooc_test_rest$estimate), p = cooc_test_rest$p.value) else NA,
  cooc_test_D1      = if (!is.null(cooc_test_D1))   c(rho = unname(cooc_test_D1$estimate),   p = cooc_test_D1$p.value)   else NA
)
saveRDS(res, file.path("results", "sc_validation_stats.rds"))

## per-cell scores for figure (compact)
cells <- rbind(
  data.frame(compartment = "CAF", state = clab, sample = csamp,
             D1 = cD1, D2 = cD2, D3 = cD3, proCAF = cpro, restCAF = crest,
             stringsAsFactors = FALSE),
  data.frame(compartment = "Malignant", state = tlab, sample = tsamp,
             D1 = tD1, D2 = tD2, D3 = tD3, proCAF = NA, restCAF = NA,
             stringsAsFactors = FALSE))
saveRDS(cells, file.path("results", "sc_validation_cells.rds"))

## ===========================================================================
## Console report
## ===========================================================================
cat("\n================= SINGLE-CELL VALIDATION (Elyada 2019) =================\n")
cat(sprintf("CAFs: %d cells | Tumour: %d cells | patients: %d (PDAC overlap: %d)\n",
            res$meta$n_caf, res$meta$n_tum, res$meta$n_patients, res$meta$n_pdac_patients))
cat(sprintf("proCAF sig genes found: %d/%d | restCAF: %d/%d | PGENES in SC namespace ok\n\n",
            res$meta$proCAF_genes_found, length(proCAF_sig),
            res$meta$restCAF_genes_found, length(restCAF_sig)))

cat("(1) CAF concordance with DeCAF single-cell signatures (Spearman, cell-level):\n")
for (nm in names(res$caf_concord)) cat(sprintf("    %-16s rho = % .3f   p %s\n",
    nm, res$caf_concord[[nm]]["rho"], fmtp(res$caf_concord[[nm]]["p"])))
cat("\n    Mean DeSurv scores by Elyada CAF subtype:\n")
print(format(caf_by_label, digits = 2))
cat(sprintf("\n    D2 across CAF subtypes (Kruskal-Wallis): p %s\n", fmtp(res$caf_D2_kruskal["p"])))
if (is.numeric(res$caf_D2_iCAF_vs_myCAF))
  cat(sprintf("    D2 iCAF (%.2f) vs myCAF (%.2f): Wilcoxon p %s\n",
      res$caf_D2_iCAF_vs_myCAF["iCAF_mean"], res$caf_D2_iCAF_vs_myCAF["myCAF_mean"],
      fmtp(res$caf_D2_iCAF_vs_myCAF["p"])))

cat("\n(2) Tumour compartment, mean DeSurv scores by malignant subtype:\n")
print(format(tum_by_label, digits = 2))
if (is.numeric(res$tum_D1_Class_vs_Basal))
  cat(sprintf("    D1 Classical (%.2f) vs Basal (%.2f): Wilcoxon p %s\n",
      res$tum_D1_Class_vs_Basal["Classical_mean"], res$tum_D1_Class_vs_Basal["Basal_mean"],
      fmtp(res$tum_D1_Class_vs_Basal["p"])))
if (is.numeric(res$tum_D3_Class_vs_Basal))
  cat(sprintf("    D3 Classical (%.2f) vs Basal (%.2f): Wilcoxon p %s\n",
      res$tum_D3_Class_vs_Basal["Classical_mean"], res$tum_D3_Class_vs_Basal["Basal_mean"],
      fmtp(res$tum_D3_Class_vs_Basal["p"])))

cat("\n(3) Within-patient co-occurrence (PDAC samples):\n")
print(format(cooc, digits = 3))
if (is.numeric(res$cooc_test_restCAF))
  cat(sprintf("\n    Classical-malignant fraction vs restCAF stroma: rho = %.3f, p %s (n=%d)\n",
      res$cooc_test_restCAF["rho"], fmtp(res$cooc_test_restCAF["p"]), nrow(cooc)))
if (is.numeric(res$cooc_test_D1))
  cat(sprintf("    Classical-malignant fraction vs CAF D1 score:   rho = %.3f, p %s (n=%d)\n",
      res$cooc_test_D1["rho"], fmtp(res$cooc_test_D1["p"]), nrow(cooc)))
cat("\nSaved -> results/sc_validation_stats.rds\n")
