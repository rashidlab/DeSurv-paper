comp_genes = function(results1,results2,pairs=NULL,type='all',save_dir){
  library(ggvenn)
  
  if(type=="all") type="factor"
  if(type=="surv") type="survival"
  
  if(type=="factor"){
    tops1=results1$tops
    tops2=results2$tops
  }else if(type=="tumor"){
    cols1 = c(grep('basal',results1$labels$label,ignore.case = TRUE),
              grep('classical',results1$labels$label,ignore.case = TRUE))
    tops1 = data.frame(tumor=unlist(results1$tops[,cols1]))
    
    cols2 = c(grep('basal',results2$labels$label,ignore.case = TRUE),
              grep('classical',results2$labels$label,ignore.case = TRUE))
    tops2 = data.frame(tumor=unlist(results2$tops[,cols2]))
    
    pairs=data.frame(alpha.0_factor.number=1,alpha.0_factor.name="",
                     alpha.best_factor.number=1,alpha.0_factor.name="")
  }else if(type=="stroma"){
    cols1 = c(grep('caf',results1$labels$label,ignore.case = TRUE),
              grep('stroma',results1$labels$label,ignore.case = TRUE),
              grep('MS',results1$labels$label,ignore.case = FALSE))
    tops1 = data.frame(tumor=unlist(results1$tops[,cols1]))
    
    cols2 = c(grep('caf',results2$labels$label,ignore.case = TRUE),
              grep('stroma',results2$labels$label,ignore.case = TRUE),
              grep('MS',results2$labels$label,ignore.case = FALSE))
    tops2 = data.frame(tumor=unlist(results2$tops[,cols2]))
    
    pairs=data.frame(alpha.0_factor.number=1,alpha.0_factor.name="",
                     alpha.best_factor.number=1,alpha.0_factor.name="")
  }else if(type=="survival"){
    cols1 = which(results1$fit_cox$beta != 0)
    tops1 = data.frame(survival=unlist(results1$tops[,cols1]))
    
    cols2 = which(results2$fit_cox$beta != 0)
    tops2 = data.frame(survival=unlist(results2$tops[,cols2]))
    
    pairs=data.frame(alpha.0_factor.number=1,alpha.0_factor.name="",
                     alpha.best_factor.number=1,alpha.0_factor.name="")
  }
  
  if(ncol(tops1)>0 | ncol(tops2)>0){
    #loop over factors
    for(i in 1:max(ncol(tops1),ncol(tops2))){
      sets = vector("list",2)
      sets[[1]] = if(ncol(tops1)>0) tops1[,pairs$alpha.0_factor.number[i],drop=TRUE] else vector("character")
      sets[[2]] = if(ncol(tops2)>0) tops2[,pairs$alpha.best_factor.number[i],drop=TRUE] else vector("character")
      if(type=="factor"){
        names(sets) = c(paste0("NMF, factor",pairs$alpha.0_factor.number[i],"\n",pairs$alpha.0_factor.name[i]),
                        paste0("alpha=best, factor",pairs$alpha.best_factor.number[i],"\n",pairs$alpha.best_factor.name[i]))
        save_name=file.path(save_dir,
                            paste0("compare.genes_factor",
                                   pairs$alpha.0_factor.number[i],".png"))
      }else{
        names(sets) = c(paste0("NMF, ",type),
                        paste0("alpha=best, ",type))
        save_name=file.path(save_dir,paste0("compare.genes_",type,".png"))
      }
      
      png(save_name,width=800,height=800)
      print(ggvenn(sets, show_elements=TRUE,label_sep="\n",auto_scale = TRUE,
                   text_size = 4)+coord_cartesian(clip="off"))
      dev.off()
    }
  }else{
    warning("no factors of this type detected")
  }
  
  
}


comp_models = function(results_a0,results_best){
  
  if(results_a0$ntop != results_best$ntop){
    stop("Number of top genes for the two models must match")
  }
  ntop=basename(results_a0$model_save_dir)
  save_dir = dirname(dirname(results_a0$model_save_dir))
  mod1 = basename(dirname(results_a0$model_save_dir))
  mod2 = basename(dirname(results_best$model_save_dir))
  save_dir = file.path(save_dir,"comps",paste0(mod1,"_VS_",mod2),ntop)
  if(!dir.exists(save_dir)){
    dir.create(save_dir,recursive = TRUE)
  }
  
  pairs=comp_alphas(results_a0,results_best,save_dir)
  
  comp_genes(results_a0,results_best,pairs,save_dir=save_dir)
  comp_genes(results_a0,results_best,type="tumor",save_dir=save_dir)
  comp_genes(results_a0,results_best,type="stroma",save_dir=save_dir)
  comp_genes(results_a0,results_best,type="surv",save_dir=save_dir)
  
  return(pairs)
}

### compare factors between 2 alphas
comp_alphas = function(results1,results2,save_dir){
  library(pheatmap)
  library(dplyr)
  library(clue)
  
  ora1 = results1$labels
  ora2 = results2$labels
  
  
  
  tops1 = results1$tops
  tops2 = results2$tops
  
  colnames(tops1) = ora1$label
  colnames(tops2) = ora2$label
  
  W1 = results1$fit_cox$W
  W2 = results2$fit_cox$W
  
  colnames(W1) = ora1$label
  colnames(W2) = ora2$label
  
  
  
  result <- matrix(nrow = ncol(tops1), ncol = ncol(tops2),
                   dimnames = list(colnames(tops1), colnames(tops2)))
  
  for (i in seq_along(tops1)) {
    for (j in seq_along(tops2)) {
      result[i, j] <- length(intersect(tops1[[i]], tops2[[j]]))  # or store the actual intersect() result
    }
  }
  
  p1 = pheatmap(cor(W1,W2),cluster_rows = FALSE,cluster_cols = FALSE,
                display_numbers = TRUE,number_color = 'black',fontsize_number = 10,
                main="Correlation of gene weights between NMF (rows) and best alpha (cols)")
  p2 = pheatmap(result,cluster_rows = FALSE,cluster_cols = FALSE,
                display_numbers = TRUE,number_color = 'black',fontsize_number = 10,
                main="Number of top genes overlap between NMF (rows) and best alpha (cols)")
  
  png(file.path(save_dir,"gene.overlap.bw.alphas.png"), width = 700, height = 600)
  print(p2)
  dev.off()
  
  png(file.path(save_dir,"cor.W.bw.alphas.png"), width = 700, height = 600)
  print(p1)
  dev.off()
  
  
  assignment =as.numeric(solve_LSAP(result,maximum = TRUE))
  
  dat <- data.frame(
    alpha.0_factor.number = seq_along(assignment),
    alpha.0_factor.name = colnames(W1)[seq_along(assignment)],
    alpha.best_factor.number = assignment,
    alpha.best_factor.name = colnames(W2)[assignment]
  )
  
  # return a data frame with factors for results1 in one col and factors for results 2 in other col
  return(dat)
}
