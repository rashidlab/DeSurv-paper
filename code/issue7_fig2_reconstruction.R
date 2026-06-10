## issue7_fig2_reconstruction.R
## Reproduces every quantitative result in GitHub issue #7
## (Fig 2 reconstruction-Shapley x-axis, all-gene 2D correlation, N1 relabel,
##  cross-method projection R^2, overall reconstruction R^2).
##
## Run from the repo root:  Rscript code/issue7_fig2_reconstruction.R
## Inputs (PDAC TCGA+CPTAC, k=3):
##   results/tar_fit_desurv_tcgacptac.rds   (DeSurv fit: $W, $H)
##   results/fit_std_desurvk_tcgacptac.rds  (standard NMF fit: $W, $H)
##   results/tar_data_filtered_tcgacptac.rds ($ex = genes x samples expression)

f <- readRDS("results/tar_fit_desurv_tcgacptac.rds")   # DeSurv
s <- readRDS("results/fit_std_desurvk_tcgacptac.rds")   # NMF k=3
d <- readRDS("results/tar_data_filtered_tcgacptac.rds")

X     <- as.matrix(d$ex)
genes <- rownames(X)
Wd <- as.matrix(f$W); Hd <- as.matrix(f$H); rownames(Wd) <- genes
Wn <- as.matrix(s$W); Hn <- as.matrix(s$H); rownames(Wn) <- genes
Xnorm2 <- sum(X^2)

## ===========================================================================
## §1  Reconstruction Shapley value (Fig 2C x-axis)
## v(S) = ||X||^2 - ||X - sum_{k in S} w_k h_k^T||^2 ; phi_k = avg marginal gain.
## Efficiency axiom: sum_k phi_k = v(N), so normalized shares sum to 100%.
## ===========================================================================
recon_shapley <- function(W, H) {
  K  <- ncol(W)
  Rk <- lapply(seq_len(K), function(j) W[, j] %o% H[j, ])
  v  <- function(S) if (!length(S)) 0 else Xnorm2 - sum((X - Reduce(`+`, Rk[S]))^2)
  subs <- do.call(c, lapply(0:K, function(m) combn(K, m, simplify = FALSE)))
  phi <- numeric(K)
  for (k in seq_len(K)) for (S in Filter(function(z) !(k %in% z), subs)) {
    w <- factorial(length(S)) * factorial(K - length(S) - 1) / factorial(K)
    phi[k] <- phi[k] + w * (v(c(S, k)) - v(S))
  }
  list(phi = phi, share = phi / sum(phi), vN = v(seq_len(K)))
}
sh_n <- recon_shapley(Wn, Hn); sh_d <- recon_shapley(Wd, Hd)
cat("=== §1  Reconstruction Shapley (normalized share, sums to 100%) ===\n")
cat(sprintf("NMF    N1/N2/N3 : %s   (sum=%.3f, vN=%.4g)\n",
            paste0(sprintf("%.1f%%", 100 * sh_n$share), collapse = " / "),
            sum(sh_n$share), sh_n$vN))
cat(sprintf("DeSurv D1/D2/D3 : %s   (sum=%.3f, vN=%.4g)\n\n",
            paste0(sprintf("%.1f%%", 100 * sh_d$share), collapse = " / "), sum(sh_d$share), sh_d$vN))

## ===========================================================================
## §3  All-gene Spearman correlation (Fig 2D)
## ===========================================================================
stopifnot(identical(rownames(Wn), rownames(Wd)))
C <- cor(Wn, Wd, method = "spearman")
dimnames(C) <- list(paste0("N", 1:3), paste0("D", 1:3))
cat("=== §3  All-gene Spearman correlation (rows=NMF, cols=DeSurv; n=", nrow(Wn), " genes) ===\n", sep = "")
print(round(C, 2)); cat("\n")

## ===========================================================================
## §4  Marker loadings for the N1 "classical" -> "tumor" relabel
## (each factor scaled to max = 1).  N1<->D2 (stroma) correlation is in C above.
## ===========================================================================
sc <- function(W) sweep(W, 2, apply(W, 2, max), "/")
Wns <- sc(Wn); Wds <- sc(Wd)
sets <- list(
  `shared epithelial` = c("KRT8","KRT18","KRT19","CDH1","EPCAM","MUC1","KRT7","CEACAM5","CEACAM6","ELF3","CLDN4","CLDN3","KRT23"),
  `classical`         = c("GATA6","TFF1","TFF2","TFF3","REG4","AGR2","LYZ","CTSE","ANXA10","FAM3D","BTNL8","VSIG2","CLRN3"),
  `basal`             = c("KRT17","KRT5","KRT6A","KRT6C","KRT14","KRT15","S100A2","TP63","SPRR1B","SPRR3","LY6D","FAM83A","GPR87","DHRS9"))
cat("=== §4  Mean marker loading (scaled to max=1) ===\n")
cat(sprintf("%-18s %6s %6s %6s %6s %6s\n", "set", "N1", "D3", "D2", "N3", "D1"))
for (nm in names(sets)) {
  g <- intersect(sets[[nm]], genes)
  cat(sprintf("%-18s %6.2f %6.2f %6.2f %6.2f %6.2f\n", nm,
              mean(Wns[g, 1]), mean(Wds[g, 3]), mean(Wds[g, 2]), mean(Wns[g, 3]), mean(Wds[g, 1])))
}
cat(sprintf("N1<->D2 (stroma) Spearman correlation: %.2f\n\n", C["N1", "D2"]))

## ===========================================================================
## §5  Cross-method projection R^2 (both directions): how reconstructable is
## each factor's gene-loading vector from the *other* method's basis (OLS R^2).
## ===========================================================================
capt <- function(src, basis) vapply(seq_len(ncol(src)),
  function(j) summary(lm(src[, j] ~ basis))$r.squared, numeric(1))
r2_N <- capt(Wn, Wd); r2_D <- capt(Wd, Wn)
cat("=== §5  Projection R^2 (factor reconstructable from the OTHER basis) ===\n")
cat(sprintf("NMF  -> DeSurv basis : N1=%.2f  N2=%.2f  N3=%.2f\n", r2_N[1], r2_N[2], r2_N[3]))
cat(sprintf("DeSurv -> NMF basis  : D1=%.2f  D2=%.2f  D3=%.2f   (D1 ~ novel)\n\n", r2_D[1], r2_D[2], r2_D[3]))

## ===========================================================================
## §6  Overall reconstruction R^2 (whole matrix) for NMF vs DeSurv
## ===========================================================================
r2_overall <- function(W, H) 1 - sum((X - W %*% H)^2) / Xnorm2
cat("=== §6  Overall reconstruction R^2 ===\n")
cat(sprintf("NMF    : R^2 = %.3f  (residual fraction %.3f)\n", r2_overall(Wn, Hn), 1 - r2_overall(Wn, Hn)))
cat(sprintf("DeSurv : R^2 = %.3f  (residual fraction %.3f)\n", r2_overall(Wd, Hd), 1 - r2_overall(Wd, Hd)))
