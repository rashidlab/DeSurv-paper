#!/usr/bin/env Rscript
# 27_program_specificity.R
# ---------------------------------------------------------------------------
# Source specificity of gene-level supervised predictors vs DeSurv programs.
#
# Question: genes chosen because their bulk expression predicts survival
# (sparse Cox, Cox-screened supervised PCA, Cox-PLS) -- are they preferentially
# loaded on a single DeSurv program, or shared across D1-D3?  And how much does
# each gene-level predictor's score move with D2/D3 conditional on D1?
#
# Bulk expression of gene g under the fit is X_g ~ sum_j W_gj H_j.  The share of
# gene g's fitted expression attributable to program j is
#     p_gj = W_gj * mean(H_j) / sum_j' W_gj' * mean(H_j')
# (scaling by mean(H_j) removes the arbitrary column scale of W).  Specificity
# = max_j p_gj (fraction on the dominant program); we also report a normalized
# entropy 1 - H(p_g)/log(k).
#
# NOTE on circularity: DeSurv's published top-270 lists are selected by the
# exemplar rule in R/get_top_genes.R (column-max-normalized loading minus the
# max over other programs), i.e. BY specificity.  We therefore also report a
# non-exemplar DeSurv reference (top-270 by raw D1 loading) so the comparison
# is not only "a specificity-selected list is specific".
#
# Supervised gene sets reproduce code/15_supervised_recovery.R exactly
# (same data, same standardization, set.seed(1) before cv.glmnet).
# Output: results/program_specificity_stats.rds (not a cache_or_compute step).
# ---------------------------------------------------------------------------

suppressMessages({ library(survival); library(glmnet) })

dat <- readRDS("results/tar_data_filtered_tcgacptac.rds")
ex  <- dat$ex; si <- dat$sampInfo
fit <- readRDS("results/tar_fit_desurv_tcgacptac.rds")
W   <- fit$W; H <- fit$H; k <- ncol(W)
stopifnot(identical(rownames(W), rownames(ex)))
source("R/get_top_genes.R")
NTOP <- 270

## --- program-share matrix and specificity ------------------------------------
C   <- sweep(W, 2, rowMeans(H), `*`)                 # genes x k, contribution scale
P   <- C / pmax(rowSums(C), .Machine$double.eps)      # share of fitted expression
spec_max <- apply(P, 1, max)
ent      <- apply(P, 1, function(p) { p <- p[p > 0]; -sum(p * log(p)) })
spec_ent <- 1 - ent / log(k)
dom      <- apply(P, 1, which.max)
names(spec_max) <- names(spec_ent) <- names(dom) <- rownames(W)

## --- DeSurv gene sets ---------------------------------------------------------
tg <- get_top_genes(W, NTOP)$top_genes
sets <- list(
  "DeSurv D1 (exemplar top-270)" = as.character(tg$factor1),
  "DeSurv D2 (exemplar top-270)" = as.character(tg$factor2),
  "DeSurv D3 (exemplar top-270)" = as.character(tg$factor3),
  "DeSurv D1 (top-270 by raw loading)" =
    rownames(W)[order(W[, 1], decreasing = TRUE)[1:NTOP]])

## --- supervised gene sets (replicates code/15) --------------------------------
sh  <- rownames(W)
Xg0 <- t(as.matrix(ex[sh, , drop = FALSE]))
gsd <- apply(Xg0, 2, sd)
Xg  <- scale(Xg0[, gsd > 0, drop = FALSE]); sh <- colnames(Xg)
mu  <- attr(Xg, "scaled:center"); sdv <- attr(Xg, "scaled:scale")
y   <- Surv(si$time, si$event)

cox_z <- apply(Xg, 2, function(g) summary(coxph(y ~ g))$coefficients[1, "z"])
keep  <- order(abs(cox_z), decreasing = TRUE)[1:NTOP]
spc   <- prcomp(Xg[, keep], center = FALSE, scale. = FALSE)
load_SPC1 <- setNames(rep(0, ncol(Xg)), sh); load_SPC1[keep] <- spc$rotation[, 1]

mres  <- residuals(coxph(y ~ 1), type = "martingale")
w_pls <- as.numeric(t(Xg) %*% mres); w_pls <- w_pls / sqrt(sum(w_pls^2))
load_PLS1 <- setNames(w_pls, sh)

set.seed(1)
cvf   <- cv.glmnet(Xg, y, family = "cox", alpha = 1, nfolds = 10)
bcoef <- as.numeric(coef(cvf, s = "lambda.min"))
if (sum(bcoef != 0) == 0) {
  gfit <- glmnet(Xg, y, family = "cox", alpha = 1)
  bcoef <- as.numeric(gfit$beta[, ncol(gfit$beta)]) }
load_SCOX <- setNames(bcoef, sh)
n_sparse  <- sum(bcoef != 0)

sets[["Sparse Cox (nonzero coefficients)"]] <- sh[bcoef != 0]
sets[["Supervised PCA (Cox-screened top-270)"]] <- sh[keep]
sets[["Cox-PLS (top-270 by |weight|)"]] <- sh[order(abs(w_pls), decreasing = TRUE)[1:NTOP]]
sets[["All genes (background)"]] <- rownames(W)

## --- specificity summary per set ----------------------------------------------
bg <- spec_max[rownames(W)]
summ <- do.call(rbind, lapply(names(sets), function(nm) {
  g <- sets[[nm]]; s <- spec_max[g]; e <- spec_ent[g]
  data.frame(set = nm, n_genes = length(g),
             median_spec = median(s), mean_spec = mean(s),
             frac_spec_gt_0.5 = mean(s > 0.5), frac_spec_gt_0.6 = mean(s > 0.6),
             median_entropy_spec = median(e),
             dom_D1 = mean(dom[g] == 1), dom_D2 = mean(dom[g] == 2), dom_D3 = mean(dom[g] == 3),
             p_vs_background = if (nm == "All genes (background)") NA_real_ else
               wilcox.test(s, bg[setdiff(rownames(W), g)])$p.value,
             stringsAsFactors = FALSE) }))
rownames(summ) <- NULL

pair <- function(a, b) wilcox.test(spec_max[sets[[a]]], spec_max[sets[[b]]])$p.value
tests <- data.frame(
  comparison = c("D1 exemplar vs Sparse Cox", "D1 exemplar vs SupPCA screen",
                 "D1 exemplar vs Cox-PLS", "D1 raw-loading vs Sparse Cox",
                 "D1 raw-loading vs SupPCA screen", "D1 raw-loading vs Cox-PLS"),
  p = c(pair("DeSurv D1 (exemplar top-270)", "Sparse Cox (nonzero coefficients)"),
        pair("DeSurv D1 (exemplar top-270)", "Supervised PCA (Cox-screened top-270)"),
        pair("DeSurv D1 (exemplar top-270)", "Cox-PLS (top-270 by |weight|)"),
        pair("DeSurv D1 (top-270 by raw loading)", "Sparse Cox (nonzero coefficients)"),
        pair("DeSurv D1 (top-270 by raw loading)", "Supervised PCA (Cox-screened top-270)"),
        pair("DeSurv D1 (top-270 by raw loading)", "Cox-PLS (top-270 by |weight|)")))

## --- predictor decomposition into program contributions -------------------------
## score = b' X_std = sum_g b_g (X_g - mu_g)/sd_g  ~  sum_j a_j H_j + const,
## a_j = sum_g b_g W_gj / sd_g.  With S = Cov(H) (programs are correlated in the
## training fit; cor(H2,H3) ~ -0.98), the fraction of the score's fitted
## variance that moves with D2/D3 AT FIXED D1 is  a_{23}' Cov(H23 | H1) a_{23} / a'Sa.
## DeSurv scores use the pipeline convention (R/predict_validation_scores.R):
## Z_j = sum_{g in top_j} W_gj X_g on raw expression, risk = sum_j beta_j Z_j.
S <- cov(t(H))
bZ <- function(j) { b <- setNames(rep(0, length(sh)), sh)
  gl <- intersect(sets[[j]], sh); b[gl] <- W[gl, j] * sdv[gl]; b }
b_risk <- Reduce(`+`, lapply(seq_len(k), function(j) as.numeric(fit$beta)[j] * bZ(j)))
decomp <- function(bstd) {
  a <- as.numeric(t(W[sh, ]) %*% (bstd[sh] / sdv[sh]))
  tot <- drop(t(a) %*% S %*% a)
  S23g1 <- S[-1, -1] - S[-1, 1, drop = FALSE] %*% t(S[-1, 1, drop = FALSE]) / S[1, 1]
  S1g23 <- S[1, 1] - S[1, -1, drop = FALSE] %*% solve(S[-1, -1]) %*% S[-1, 1, drop = FALSE]
  sc <- as.numeric(Xg %*% bstd[sh]); d <- data.frame(sc, H1 = H[1, ], H2 = H[2, ], H3 = H[3, ])
  r1 <- summary(lm(sc ~ H1, d))$r.squared; r23 <- summary(lm(sc ~ H2 + H3, d))$r.squared
  r123 <- summary(lm(sc ~ H1 + H2 + H3, d))$r.squared
  data.frame(share_D1 = abs(a[1]) * sqrt(S[1,1]) / sum(abs(a) * sqrt(diag(S))),
             frac_var_D2D3_given_D1 = drop(t(a[-1]) %*% S23g1 %*% a[-1]) / tot,
             frac_var_D1_given_D2D3 = drop(a[1]^2 * S1g23) / tot,
             R2_program_space = r123, semipartial_R2_D2D3 = r123 - r1,
             semipartial_R2_D1 = r123 - r23, cor_H1 = cor(sc, H[1, ])) }
preds <- list("DeSurv Z1 (D1 score)" = bZ(1), "DeSurv risk score" = b_risk,
              "Sparse Cox" = load_SCOX, "Supervised PCA" = load_SPC1, "Cox-PLS" = load_PLS1)
dec <- do.call(rbind, lapply(names(preds), function(nm) cbind(predictor = nm, decomp(preds[[nm]]))))
rownames(dec) <- NULL

## --- sparse Cox genes: sign vs dominant program ----------------------------------
scox_tab <- data.frame(gene = names(load_SCOX)[load_SCOX != 0], beta = load_SCOX[load_SCOX != 0])
scox_tab$dominant <- dom[scox_tab$gene]; scox_tab$spec_max <- spec_max[scox_tab$gene]
scox_tab$in_D1_exemplar <- scox_tab$gene %in% sets[[1]]; scox_tab$in_D3_exemplar <- scox_tab$gene %in% sets[[3]]
rownames(scox_tab) <- NULL

PGENES <- unique(unlist(lapply(tg, as.character)))
## --- gene-level overlap: where do supervised genes sit in the program space? ----
ovl <- do.call(rbind, lapply(c("Sparse Cox (nonzero coefficients)",
                               "Supervised PCA (Cox-screened top-270)",
                               "Cox-PLS (top-270 by |weight|)"), function(nm) {
  g <- sets[[nm]]
  data.frame(set = nm,
             in_D1_exemplar = mean(g %in% sets[[1]]), in_D2_exemplar = mean(g %in% sets[[2]]),
             in_D3_exemplar = mean(g %in% sets[[3]]),
             in_no_exemplar_list = mean(!g %in% PGENES), stringsAsFactors = FALSE) }))

res <- list(meta = list(n = ncol(ex), n_genes = nrow(W), k = k, ntop = NTOP,
                        n_sparse_cox_genes = n_sparse, hyper = fit$hyper[c("alpha","lambda","nu")]),
            gene_spec = data.frame(gene = rownames(W), spec_max = spec_max,
                                   spec_entropy = spec_ent, dominant = dom, P),
            sets = sets, summary = summ, tests = tests, decomposition = dec, overlap = ovl,
            sparse_cox_genes = scox_tab, H_cor = cor(t(H)))
saveRDS(res, "results/program_specificity_stats.rds")

cat("\n=============== PROGRAM SPECIFICITY OF GENE SETS ===============\n")
cat(sprintf("n=%d, %d genes, k=%d, sparse Cox selected %d genes\n\n", ncol(ex), nrow(W), k, n_sparse))
print(format(summ[, c("set","n_genes","median_spec","frac_spec_gt_0.5","frac_spec_gt_0.6","dom_D1","dom_D2","dom_D3","p_vs_background")], digits = 3), row.names = FALSE)
cat("\nPairwise Wilcoxon on specificity:\n"); print(format(tests, digits = 3), row.names = FALSE)
cat("\n=============== PREDICTOR DECOMPOSITION INTO D1/D2/D3 ===============\n")
print(format(dec, digits = 3), row.names = FALSE)
cat("\nSparse Cox genes, sign x dominant program:\n"); print(table(sign = sign(scox_tab$beta), dominant = scox_tab$dominant))
cat("\nSupervised gene sets vs DeSurv exemplar lists:\n"); print(format(ovl, digits = 3), row.names = FALSE)
