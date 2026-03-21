# Writes data/derv/<datasets_joined>_formatted.rds and returns its path
load_data <- function(datasets) {
  stopifnot(length(datasets) >= 1, is.character(datasets))
  
  datasets_key <- paste0(datasets, collapse = ".")

  if (length(datasets) == 1) {
    data <- load_data_internal(dataname = datasets)
  } else if (length(datasets) > 1) {
    dat_list <- lapply(datasets, function(ds) load_data_internal(dataname = ds))
    
    # take the union of all feature IDs
    all_genes <- Reduce(union, lapply(dat_list, function(x) x$featInfo))
    
    # build aligned expression matrices with outer join
    ex_list <- lapply(dat_list, function(x) {
      m <- matrix(0, nrow = length(all_genes), ncol = ncol(x$ex),
                  dimnames = list(all_genes, colnames(x$ex)))
      m[x$featInfo, ] <- x$ex
      m
    })
    
    data <- list()
    data$ex <- do.call(cbind, ex_list)
    data$featInfo <- all_genes
    
    # row-bind sample info
    data$sampInfo <- do.call("rbind", lapply(dat_list, \(x) x$sampInfo))
    data$dataname <- datasets_key
  } else {
    stop("`datasets` must be a dataset name or a character vector of dataset names.")
  }
  
  data$samp_keeps <- which(data$sampInfo$keep == 1L)
  
  return(data)
}
