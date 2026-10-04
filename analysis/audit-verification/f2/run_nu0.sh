cd /home/naimrashid/Downloads/DeSurv-paper-clean
HERE=~/Downloads/DeSurv-paper-recompute/analysis/audit-verification/f2
export R_LIBS=${DESURV_120_LIB:?set to a library with DeSurv 1.2.0}
Rscript -e 'suppressMessages(library(DeSurv)); stopifnot(packageVersion("DeSurv")=="1.2.0"); f<-readRDS("results/tar_fit_desurv_tcgacptac.rds"); fit<-desurv_fit(f$data$X, f$data$y, f$data$d, k=12, alpha=0.658658, lambda=0.74881, nu=0, lambdaW=0, lambdaH=0, seed=1, tol=1e-7, maxit=4000, ninit=30, verbose=FALSE); saveRDS(fit, "'"$HERE"'/fit_k12_nu0.rds")'
Rscript $HERE/f2_check.R $HERE/fit_k12_nu0.rds "inline 1.2.0 k12, eval-40 config with nu=0" > $HERE/out_k12_nu0.txt 2>&1
Rscript $HERE/f2_reopt.R $HERE/fit_k12_nu0.rds "inline k12 nu=0 (1.2.0)" $HERE/reopt_k12_nu0.rds 162 400 > $HERE/reopt_k12_nu0.txt 2>&1
Rscript $HERE/f2_val.R $HERE/reopt_k12_nu0.rds 162 $HERE/fit_k12_nu0.rds > $HERE/val_k12_nu0.txt 2>&1
echo DONE
