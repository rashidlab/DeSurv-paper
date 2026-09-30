# Rev 3: the CV folds are drawn ONCE from `fold_seed` and reused for every nu, so the
# nu comparison is on common folds and the result is reproducible (previously
# cv.glmnet drew new unseeded folds for each nu). The caller's RNG state is restored.
fit_cox_model = function(XtW,df,nfold=5,fold_seed=20261202L){
  # param grid
  nu_grid = seq(0,1,by=.05)
  
  # fixed fold assignment, drawn without disturbing the caller's RNG stream
  had_seed = exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed = get(".Random.seed", envir = globalenv(), inherits = FALSE)
  set.seed(fold_seed)
  foldid = sample(rep(seq_len(nfold), length.out = nrow(XtW)))
  if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else
    rm(".Random.seed", envir = globalenv())
  
  # loop over nu and fit glmnet
  cv = data.frame(matrix(ncol=3,nrow=0))
  cvfits = list()
  for(nu in nu_grid){
    cvfit = cv.glmnet(x=XtW,
                      y=Surv(df$time,df$event),
                      family="cox",
                      type.measure = "C",
                      foldid = foldid,
                      alpha = nu)
    temp = data.frame(lambda=cvfit$lambda,cvm=cvfit$cvm,cvsd=cvfit$cvsd)
    temp$nu = nu
    cv = rbind(cv,temp)
    cvfits[[as.character(nu)]] = cvfit
  }
  
  tops = cv %>% group_by(nu) %>%
    slice_max(order_by = cvm,n=1) %>%
    ungroup()
  best=tops[which.max(tops$cvm),]
  sd = best$cvm - best$cvsd
  selected = tops %>% filter(cvm > sd) %>% slice_max(order_by = lambda,n=1)
  mod = cvfits[[as.character(selected$nu)]]
  beta = coef(mod,s=selected$lambda)
  
  beta
}
