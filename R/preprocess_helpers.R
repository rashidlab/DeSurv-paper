# R/preprocess_helpers.R
# Preprocessing wrappers for training and validation data.

infer_dataset_name <- function(dataset, fallback = "validation") {
  dataname <- dataset$dataname
  if (!is.null(dataname) && nzchar(dataname)) {
    return(dataname)
  }
  dataset_ids <- unique(dataset$sampInfo$dataset)
  dataset_ids <- dataset_ids[!is.na(dataset_ids)]
  if (length(dataset_ids)) {
    dataset_ids[[1]]
  } else {
    fallback
  }
}

preprocess_training_data <- function(data, ngene, method_trans_train) {
  keep_idx <- data$samp_keeps
  X        <- data$ex
  y        <- data$sampInfo$time
  d        <- data$sampInfo$event
  dataset  <- data$sampInfo$dataset
  if (!is.null(keep_idx)) {
    if (is.logical(keep_idx)) keep_idx <- which(keep_idx)
    if (!length(keep_idx)) stop("No samples selected after applying `samp_keeps`.")
    X       <- X[, keep_idx, drop = FALSE]
    y       <- y[keep_idx]
    d       <- d[keep_idx]
    dataset <- dataset[keep_idx]
  }
  prep <- DeSurv::preprocess_data(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    samp_keeps         = NULL,
    ngene              = ngene,
    method_trans_train = method_trans_train,
    verbose            = FALSE
  )
  prep$dataname <- data$dataname
  prep
}

preprocess_validation_data <- function(dataset, genes = NULL, ngene = NULL,
                                       method_trans_train, dataname,
                                       transform_target = NULL,
                                       zero_fill_missing = FALSE) {
  if (is.null(genes) && is.null(ngene)) {
    stop("Provide either genes or ngene for validation preprocessing.")
  }
  keep_idx <- dataset$samp_keeps
  if (is.null(keep_idx) && !is.null(dataset$sampInfo$keep)) {
    keep_idx <- which(dataset$sampInfo$keep == 1)
  }
  X           <- dataset$ex
  y           <- dataset$sampInfo$time
  d           <- dataset$sampInfo$event
  dataset_ids <- dataset$sampInfo$dataset
  if (!is.null(keep_idx)) {
    if (is.logical(keep_idx)) keep_idx <- which(keep_idx)
    if (!length(keep_idx)) stop("No samples selected after applying `samp_keeps`.")
    X           <- X[, keep_idx, drop = FALSE]
    y           <- y[keep_idx]
    d           <- d[keep_idx]
    dataset_ids <- dataset_ids[keep_idx]
  }
  if (isTRUE(zero_fill_missing) && !is.null(genes)) {
    missing_genes <- setdiff(genes, rownames(X))
    if (length(missing_genes) > 0L) {
      message(sprintf(
        "Zero-filling %d of %d requested genes not found in %s: %s%s",
        length(missing_genes), length(genes), dataname,
        paste(head(missing_genes, 5L), collapse = ", "),
        if (length(missing_genes) > 5L) ", ..." else ""
      ))
      zero_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(X),
                         dimnames = list(missing_genes, colnames(X)))
      X <- rbind(X, zero_mat)
    }
  }
  args <- list(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset_ids,
    samp_keeps         = NULL,
    method_trans_train = method_trans_train,
    verbose            = FALSE
  )
  if (!is.null(genes)) args$genes <- genes else args$ngene <- ngene
  if (!is.null(transform_target)) args$transform_target <- transform_target
  prep <- do.call(DeSurv::preprocess_data, args)
  prep$dataname <- dataname
  meta <- dataset$sampInfo
  if (!is.null(meta) && nrow(meta)) {
    prep_ids <- colnames(prep$ex)
    if (is.null(prep_ids) && !is.null(rownames(prep$sampInfo))) prep_ids <- rownames(prep$sampInfo)
    if (is.null(prep_ids) && "ID" %in% names(prep$sampInfo))    prep_ids <- prep$sampInfo$ID
    if (!is.null(prep_ids) && length(prep_ids)) {
      if (!"ID" %in% names(meta) && !is.null(rownames(meta))) meta$ID <- rownames(meta)
      if ("ID" %in% names(meta)) {
        idx      <- match(prep_ids, meta$ID)
        meta_use <- meta[idx, , drop = FALSE]
        rownames(meta_use) <- prep_ids
        prep$sampInfo <- as.data.frame(prep$sampInfo, stringsAsFactors = FALSE)
        rownames(prep$sampInfo) <- prep_ids
        extra_cols <- setdiff(names(meta_use), names(prep$sampInfo))
        if (length(extra_cols)) prep$sampInfo[extra_cols] <- meta_use[extra_cols]
      }
    }
  }
  prep
}
