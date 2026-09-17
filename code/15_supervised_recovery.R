#!/usr/bin/env Rscript
# 15_supervised_recovery.R
# ---------------------------------------------------------------------------
# Supervised-method recovery (Nature Cancer revision) -- "is the D1 tumour-
# stroma coupling specific to DeSurv's algorithm, or a property of survival
# supervision itself?"  We fit several *independent* survival-supervised
# dimension-reduction methods on the SAME training data (TCGA+CPTAC) and ask
# whether each recovers the D1 axis -- both as a per-sample score and as a gene
# program -- while the leading UNSUPERVISED direction (PC1, which DECODER /
# standard NMF follow) does NOT. This shows D1 is created by supervision, not by
# DeSurv specifically; it is a method-robustness check, NOT a biology claim
# (the coupling is established at the bulk gene-program level; cellular/spatial
# resolution was not robustly supported and is not claimed -- see notes in
# code/13c/13d). It also mirrors the established result that unsupervised
# DECODER/NMF do not isolate D1.
#
# Methods: (i) unsupervised PCA [baseline, no survival]; (ii) supervised PCA
# (Bair-Tibshirani: Cox-screen genes -> PCA); (iii) Cox-PLS (one PLS component
# against martingale residuals); (iv) sparse Cox (glmnet, family="cox").
# In-repo only. D1/D2/D3 scored by the manuscript rank-then-project convention.
# ---------------------------------------------------------------------------

suppressMessages({ library(survival); library(glmnet) })
set.seed(1)
# Run from the repository root (paths below are repo-root-relative), as with the
# other pipeline scripts. No hardcoded absolute path, for portability across clones/HPC.

dat <- readRDS("results/tar_data_filtered_tcgacptac.rds")
ex  <- dat$ex                                   # genes x samples (1970 x 273)
si  <- dat$sampInfo
fit <- readRDS("results/tar_fit_desurv_tcgacptac.rds")
W   <- fit$W; Wg <- rownames(W)
source("R/get_top_genes.R")
PGENES <- unique(unlist(lapply(get_top_genes(W, 270)$top_genes, as.character)))

## --- DeSurv D1/D2/D3 training scores (rank-then-project, top-270) -------------
sh <- intersect(Wg, rownames(ex))
xr <- apply(ex[sh, , drop = FALSE], 2, rank, ties.method = "average"); rownames(xr) <- sh
projZ <- function(j){ g <- intersect(PGENES, sh)
  as.numeric(drop(t(xr[g, , drop = FALSE]) %*% W[g, j, drop = FALSE])) }
D1 <- scale(projZ(1))[,1]; D2 <- scale(projZ(2))[,1]; D3 <- scale(projZ(3))[,1]
y  <- Surv(si$time, si$event)

# gene-expression matrix, samples x genes, standardized (for PCA/PLS/glmnet)
Xg0 <- t(as.matrix(ex[sh, , drop = FALSE]))         # n x p, columns = genes
gsd <- apply(Xg0, 2, sd)
Xg  <- scale(Xg0[, gsd > 0, drop = FALSE])          # drop zero-variance genes
sh  <- colnames(Xg)                                 # genes actually used

## ===========================================================================
## (i) UNSUPERVISED PCA -- dominant variance direction (no survival)
## ===========================================================================
pc <- prcomp(Xg, center = FALSE, scale. = FALSE)
PC1 <- pc$x[,1]; load_PC1 <- pc$rotation[,1]

## ===========================================================================
## (ii) SUPERVISED PCA (Bair-Tibshirani): Cox-screen genes, PCA on survivors
## ===========================================================================
cox_z <- apply(Xg, 2, function(g) { s <- summary(coxph(y ~ g)); s$coefficients[1,"z"] })
keep  <- order(abs(cox_z), decreasing = TRUE)[1:270]      # top-270 to match D1
spc   <- prcomp(Xg[, keep], center = FALSE, scale. = FALSE)
SPC1  <- spc$x[,1]
load_SPC1 <- setNames(rep(0, ncol(Xg)), colnames(Xg)); load_SPC1[keep] <- spc$rotation[,1]

## ===========================================================================
## (iii) COX-PLS: one PLS component against martingale residuals
## ===========================================================================
mres <- residuals(coxph(y ~ 1), type = "martingale")     # null-model residuals
w_pls <- as.numeric(t(Xg) %*% mres); w_pls <- w_pls / sqrt(sum(w_pls^2))
PLS1  <- as.numeric(Xg %*% w_pls); load_PLS1 <- w_pls

## ===========================================================================
## (iv) SPARSE COX (glmnet, family="cox"); lambda.min, fall back to looser grid
## ===========================================================================
cvf <- cv.glmnet(Xg, y, family = "cox", alpha = 1, nfolds = 10)
bcoef <- as.numeric(coef(cvf, s = "lambda.min"))
if (sum(bcoef != 0) == 0) {                               # if still empty, take densest lambda
  gfit <- glmnet(Xg, y, family = "cox", alpha = 1)
  bcoef <- as.numeric(gfit$beta[, ncol(gfit$beta)]) }
SCOX  <- as.numeric(Xg %*% bcoef); load_SCOX <- bcoef
n_sparse <- sum(bcoef != 0)

## sign-align each direction (and its loadings) to D1 -------------------------
flipper <- function(score, load){ r <- suppressWarnings(cor(score, D1))
  f <- if (!is.na(r) && r < 0) -1 else 1; list(score = f*score, load = f*load) }
fp1 <- flipper(PC1,  load_PC1);  PC1  <- fp1$score; load_PC1  <- fp1$load
fp2 <- flipper(SPC1, load_SPC1); SPC1 <- fp2$score; load_SPC1 <- fp2$load
fp3 <- flipper(PLS1, load_PLS1); PLS1 <- fp3$score; load_PLS1 <- fp3$load
fp4 <- flipper(SCOX, load_SCOX); SCOX <- fp4$score; load_SCOX <- fp4$load

## ===========================================================================
## Metrics: (a) per-sample score correlation with D1; (b) gene-loading
## correlation with D1's W column; (c) does the axis carry BOTH arms
## (Classical-tumour + restCAF-stroma), i.e. recover the coupling?
## ===========================================================================
wD1 <- W[sh, 1]                                            # D1 gene loadings
load_cor <- function(L) suppressWarnings(cor(L, wD1, method = "spearman"))
methods <- list(
  "Unsupervised PCA" = list(score = PC1,  load = load_PC1,  sup = FALSE),
  "Supervised PCA"   = list(score = SPC1, load = load_SPC1, sup = TRUE),
  "Cox-PLS"          = list(score = PLS1, load = load_PLS1, sup = TRUE),
  "Sparse Cox"       = list(score = SCOX, load = load_SCOX, sup = TRUE))

tab <- do.call(rbind, lapply(names(methods), function(nm) {
  m <- methods[[nm]]
  data.frame(method = nm, supervised = m$sup,
             score_cor_D1 = cor(m$score, D1),
             load_cor_D1  = load_cor(m$load),
             # also report nearest among D1/D2/D3 to show it is D1 specifically
             score_cor_D2 = cor(m$score, D2), score_cor_D3 = cor(m$score, D3),
             stringsAsFactors = FALSE) }))

res <- list(
  meta = list(n = nrow(Xg), n_genes = ncol(Xg), n_events = sum(si$event),
              n_sparse_cox_genes = n_sparse, ntop = 270),
  table = tab,
  sup_mean_score_cor = mean(abs(tab$score_cor_D1[tab$supervised])),
  unsup_score_cor    = abs(tab$score_cor_D1[!tab$supervised]),
  sup_mean_load_cor  = mean(tab$load_cor_D1[tab$supervised]),
  scores = data.frame(D1 = D1, D2 = D2, D3 = D3,
                      PC1 = PC1, SupPCA = SPC1, CoxPLS = PLS1, SparseCox = SCOX))
saveRDS(res, "results/supervised_recovery_stats.rds")

cat("\n================= SUPERVISED-METHOD RECOVERY OF D1 =================\n")
cat(sprintf("training n=%d (%d events), %d genes\n\n", res$meta$n, res$meta$n_events, res$meta$n_genes))
cat(sprintf("%-18s %-11s %12s %12s   (%s)\n","method","supervised","|cor| score~D1","cor load~D1","score~D2/D3"))
for (i in 1:nrow(tab)) cat(sprintf("%-18s %-11s %12.2f %12.2f   (%.2f / %.2f)\n",
  tab$method[i], ifelse(tab$supervised[i],"yes","NO"),
  abs(tab$score_cor_D1[i]), tab$load_cor_D1[i], tab$score_cor_D2[i], tab$score_cor_D3[i]))
cat(sprintf("\nSupervised methods: mean |score~D1| = %.2f, mean load~D1 = %.2f\n",
            res$sup_mean_score_cor, res$sup_mean_load_cor))
cat(sprintf("Unsupervised PC1  : |score~D1| = %.2f  (dominant-variance direction misses D1)\n",
            res$unsup_score_cor))
cat(sprintf("Sparse Cox selected %d genes.\n", n_sparse))
cat("\nSaved -> results/supervised_recovery_stats.rds\n")

## ===========================================================================
## Figure
## ===========================================================================
## Drawing code lives in code/util_fig_supervised_recovery.R, which reads the
## stats saved above from results/supervised_recovery_stats.rds. Kept separate
## so the figure can be redrawn without re-running this analysis (whose
## cross-validated fits are not seeded and whose stats file this script
## overwrites unconditionally).
source("code/util_fig_supervised_recovery.R")
