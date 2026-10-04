D <- path.expand("~/Downloads/DeSurv-paper-recompute/analysis/audit-verification/sub")
suppressPackageStartupMessages(library(Rcpp)); sourceCpp(file.path(D,"cox101.cpp"), cacheDir=file.path(D,"rcppcache"))
.libPaths(c(Sys.getenv("DESURV_120_LIB"), .libPaths()))
s <- readRDS(file.path(D,"lp_sorted.rds")); y <- s$y; d <- s$d
tied <- y %in% y[duplicated(y)]
cat("tied values", sum(duplicated(y)), "distinct tied values", length(unique(y[tied])), " subjects at tied times", sum(tied), " events at tied times", sum(d[tied]), "\n")
a <- cox_stats_101(s$lp,d,y); b <- DeSurv:::cox_stats_cpp(s$lp,d,y)
e <- a$g_z-b$g_z
cat("at fitted lp: subjects with g differing", sum(abs(e)>1e-12), " max|dg|", max(abs(e)), " ||dg||/||g120||", sqrt(sum(e^2))/sqrt(sum(b$g_z^2)), " sum g101", sum(a$g_z), "\n")
cat("w rel diff", sqrt(sum((a$w-b$w)^2))/sqrt(sum(b$w^2)), " loglik equal:", all.equal(a$loglik,b$loglik), "\n")
