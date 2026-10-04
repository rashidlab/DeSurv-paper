suppressMessages({library(DeSurv); library(survival)}); source("R/load_data_internal.R"); source("R/get_top_genes.R"); source("R/preprocess_helpers.R")
dat <- readRDS("results/tar_data_filtered_tcgacptac.rds"); genes <- rownames(dat$ex)
fit <- readRDS("results/tar_fit_desurv_tcgacptac.rds"); sup <- unique(unlist(get_top_genes(fit$W,270)$top_genes))
dv <- readRDS("results/data_val_filtered_tcgacptac.rds"); vl <- readRDS("results/val_latent_desurv_tcgacptac.rds")
cat("val_latent names:", names(vl), "\n")
C <- function(s,t,e) unname(concordance(Surv(t,e)~s, reverse=TRUE)$concordance)
for (i in seq_along(dv)) { co <- dv[[i]]; ds <- co$sampInfo$dataset[1]
  lp <- as.numeric(t(co$ex[sup,]) %*% fit$W[sup,] %*% fit$beta)
  v <- vl[[which(sapply(vl, function(e) e$dataset)==ds)]]
  agree <- max(abs(lp - as.numeric(v$risk_score)))
  raw <- load_data_internal(dataname=ds); av <- intersect(genes, rownames(raw$ex))
  alt <- suppressMessages(preprocess_validation_data(raw, genes=av, ngene=NULL, method_trans_train="rank", dataname=ds, zero_fill_missing=FALSE))
  stopifnot(identical(colnames(alt$ex), colnames(co$ex)))
  sa <- intersect(sup, rownames(alt$ex)); lpa <- as.numeric(t(alt$ex[sa,]) %*% fit$W[sa,] %*% fit$beta)
  t <- co$sampInfo$time; e <- co$sampInfo$event
  d <- lp - lpa
  cat(sprintf("%-18s |lp-cached| %.1e  C zero-fill %.7f  C drop-missing %.7f  diff %+.5f  sd(lp diff)/sd(lp) %.2e  cor %.6f\n", ds, agree, C(lp,t,e), C(lpa,t,e), C(lpa,t,e)-C(lp,t,e), sd(d)/sd(lp), cor(lp,lpa))) }
