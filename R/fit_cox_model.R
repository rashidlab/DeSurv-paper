fit_cox_model = function(XtW,df,nfold=5){
  # param grid
  nu_grid = seq(0,1,by=.05)
  
  
  # loop over nu and fit glmnet
  cv = data.frame(matrix(ncol=3,nrow=0))
  cvfits = list()
  for(nu in nu_grid){
    cvfit = cv.glmnet(x=XtW,
                      y=Surv(df$time,df$event),
                      family="cox",
                      type.measure = "C",
                      nfolds = nfold,
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
