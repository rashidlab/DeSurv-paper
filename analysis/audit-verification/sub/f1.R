D <- path.expand("~/Downloads/DeSurv-paper-recompute/analysis/audit-verification/sub")
Sys.setenv(R_LIBS="")
suppressPackageStartupMessages({library(Rcpp); library(survival)})
sourceCpp(file.path(D,"cox101.cpp"), cacheDir=file.path(D,"rcppcache"))
.libPaths(c(Sys.getenv("DESURV_120_LIB"), .libPaths()))
library(DeSurv); cat("1.2.0 loaded:", as.character(packageVersion("DeSurv")), find.package("DeSurv"), "\n")
cs120 <- DeSurv:::cox_stats_cpp
# independent Breslow loglik
bres <- function(eta,d,t){ sum(sapply(which(d==1), function(i) eta[i]-log(sum(exp(eta[t>=t[i]]))))) }
fdg <- function(eta,d,t,h=1e-6) sapply(seq_along(eta), function(i){e1<-e2<-eta;e1[i]<-e1[i]+h;e2[i]<-e2[i]-h;(bres(e1,d,t)-bres(e2,d,t))/(2*h)})
fdw <- function(eta,d,t,h=1e-4) sapply(seq_along(eta), function(i){e1<-e2<-eta;e1[i]<-e1[i]+h;e2[i]<-e2[i]-h;-(bres(e1,d,t)-2*bres(eta,d,t)+bres(e2,d,t))/h^2})
# Example 1: audit toy at eta=0
t <- c(1,1,2,3,4,5); d <- c(1,1,1,0,1,0); eta <- rep(0,6)
a <- cox_stats_101(eta,d,t); b <- cs120(eta,d,t)
cat("\nToy eta=0\n"); print(rbind(g101=a$g_z, g120=b$g_z, fd=fdg(eta,d,t)), digits=6)
cat("sum g101", sum(a$g_z), " sum g120", sum(b$g_z), "\n")
cat("loglik 101/120/indep:", a$loglik, b$loglik, bres(eta,d,t), "\n")
# Example 2: random tied data, many groups
set.seed(42); n <- 40; t <- sort(sample(1:12, n, TRUE)); d <- rbinom(n,1,.7); eta <- rnorm(n)
a <- cox_stats_101(eta,d,t); b <- cs120(eta,d,t); g <- fdg(eta,d,t); w <- fdw(eta,d,t)
cat("\nRandom tied n=40, unique times", length(unique(t)), "\n")
cat("max|g-fd|: 101", max(abs(a$g_z-g)), " 120", max(abs(b$g_z-g)), "\n")
cat("max|w-fd|: 101", max(abs(a$w-w)), " 120", max(abs(b$w-w)), "\n")
cat("loglik diff 101", a$loglik-bres(eta,d,t), " 120", b$loglik-bres(eta,d,t), "\n")
cat("sum g (shift invariance, should be 0): 101", sum(a$g_z), " 120", sum(b$g_z), "\n")
cat("subjects with wrong g in 101:", sum(abs(a$g_z-g)>1e-6), "of", n, "\n")
# Example 3: untied data, both should agree
t3 <- sort(runif(n)); a3 <- cox_stats_101(eta,d,t3); b3 <- cs120(eta,d,t3)
cat("\nUntied: max|g101-g120|", max(abs(a3$g_z-b3$g_z)), " max|g101-fd|", max(abs(a3$g_z-fdg(eta,d,t3))), "\n")
# Example 4: within-tie permutation invariance of gradient (map back to subjects)
perm <- c(2,1,3:n); # swap first two (check they are tied)
cat("t[1]==t[2]:", t[1]==t[2], "\n")
a4 <- cox_stats_101(eta[perm],d[perm],t[perm]); b4 <- cs120(eta[perm],d[perm],t[perm])
cat("perm change in subject gradient: 101", max(abs(a4$g_z[order(perm)]-a$g_z)), " 120", max(abs(b4$g_z[order(perm)]-b$g_z)), "\n")
# check vs coxph score at beta=0 with one covariate
