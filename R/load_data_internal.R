load_data_internal = function(dataname) {

  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  ### load data
  dat  = readRDS(paste0("data/original/", dataname, ".rds"))
  surv = readRDS(paste0("data/original/", dataname, ".survival_data.rds"))
  caf  = read.csv(paste0("data/original/", dataname, "_subtype.csv"))
  
  ### dataset specific edits
  if (dataname == "Dijk") {
    samp_keeps = 1:nrow(surv)
    dat$ex     = log2(dat$ex + 1)
    
  } else if (dataname == "CPTAC") {
    samp_keeps = which(dat$sampInfo$histology_diagnosis == "PDAC")
    
  } else if (dataname == "TCGA_PAAD") {
    samp_keeps       = which(!is.na(dat$sampInfo$Grade))
    rownames(dat$ex) = dat$featInfo$SYMBOL
    dat$ex           = log2(dat$ex + 1)
    
  } else if (dataname == "Grunwald") {
    samp_keeps = 30:nrow(dat$sampInfo)
    
  } else if (dataname == "Hayashi") {
    samp_keeps = which(dat$sampInfo$`Primary/_Met` %in% "Primary")
    
  } else if (dataname == "Linehan") {
    samp_keeps       = which(!is.na(dat$sampInfo$Treatment))
    rownames(dat$ex) = dat$featInfo$SYMBOL
    dat$ex           = log2(dat$ex + 1)
    
  } else if (dataname == "Moffitt_GEO_array") {
    samp_keeps = which(dat$sampInfo$specimen_type %in% c("Primary"))
    
  } else if (dataname == "Olive") {
    samp_keeps = which(dat$sampInfo$tumor %in% c("stroma") | dat$sampInfo$compartment %in% "Stroma")
    
  } else if (dataname == "PACA_AU_array") {
    samp_keeps       = which((dat$sampInfo$Sample.type %in% c("Primary tumour")) &
                             (dat$sampInfo$HistoSubtype %in% c("Pancreatic Ductal Adenocarcinoma")))
    # combine duplicated genes
    dat$ex$symbol    = dat$featInfo$SYMBOL
    dat$ex           = dplyr::as_tibble(dat$ex) |>
                       dplyr::group_by(.data$symbol) |>
                       dplyr::summarise(dplyr::across(dplyr::everything(), median), .groups = "drop")
    dat$ex           = as.data.frame(dat$ex[!is.na(dat$ex$symbol), ])
    rownames(dat$ex) = dat$ex$symbol
    dat$ex$symbol    = NULL
    
  } else if (dataname == "PACA_AU_seq") {
    samp_keeps       = which((dat$sampInfo$Sample.type %in% c("Primary tumour")) &
                             (dat$sampInfo$HistoSubtype %in% c("Pancreatic Ductal Adenocarcinoma")))
    
    # combine duplicated genes
    dat$ex$symbol    = dat$featInfo$SYMBOL
    dat$ex           = dplyr::as_tibble(dat$ex) |>
                       dplyr::group_by(.data$symbol) |>
                       dplyr::summarise(dplyr::across(dplyr::everything(), median), .groups = "drop")
    dat$ex           = as.data.frame(dat$ex[!is.na(dat$ex$symbol), ])
    rownames(dat$ex) = dat$ex$symbol
    dat$ex$symbol    = NULL
    dat$ex           = log2(dat$ex + 1)
    
    # ensure uniqueness vs array IDs
    colnames(dat$ex) = paste0(colnames(dat$ex), "_seq")
    caf$ID           = paste0(caf$ID, "_seq")
    surv$sampID      = paste0(surv$sampID,"_seq") 
    
  } else if (dataname == "Puleo_array") {
    samp_keeps = 1:nrow(dat$sampInfo)
    
  } else {
    stop("dataset not currently supported")
  }
  
  #### find observations with valid survival info
  surv_keeps = which(surv$time > 0 & !is.na(surv$time) & !is.na(surv$event))
  samp_keeps = intersect(samp_keeps,surv_keeps)
  
  #### Build sampInfo
  surv$ID                   = surv$sampID
  sampInfo                  = dplyr::left_join(surv, caf, by = "ID")
  sampInfo$event            = as.numeric(sampInfo$event)
  sampInfo$dataset          = dataname
  sampInfo$keep             = 0
  sampInfo$keep[samp_keeps] = 1
  
  ### Data checks
  # survival data valid
  if(!all(sampInfo$event[sampInfo$keep==1] %in% c(0,1))){
    stop("Event indicators must be 0 or 1")
  }
  if(!all(sampInfo$time[sampInfo$keep==1] > 0)){
    stop("Survival times must be positive values")
  }
  
  # Expression matrix must have sample IDs as column names
  if (is.null(colnames(dat$ex))) {
    stop("Expression matrix `dat$ex` must have column names (sample IDs).")
  }
  
  # sampInfo must have ID column, non-missing
  if (!"ID" %in% names(sampInfo)) {
    stop("`sampInfo` is missing 'ID' column after join.")
  }
  if (anyNA(sampInfo$ID)) {
    bad = unique(sampInfo$ID[is.na(sampInfo$ID)])
    stop("`sampInfo$ID` contains NA values; cannot align samples.")
  }
  
  # Duplicates in either source → error
  dup_ex  = unique(colnames(dat$ex)[duplicated(colnames(dat$ex))])
  if (length(dup_ex)) {
    stop(sprintf("Duplicate sample IDs in expression matrix columns: %s", paste(dup_ex, collapse = ", ")))
  }
  dup_si  = unique(sampInfo$ID[duplicated(sampInfo$ID)])
  if (length(dup_si)) {
    stop(sprintf("Duplicate sample IDs in sampInfo$ID: %s", paste(dup_si, collapse = ", ")))
  }
  
  # Sets must match exactly (difference → error)
  set_ex = sort(colnames(dat$ex))
  set_si = sort(as.character(sampInfo$ID))
  if (!identical(set_ex, set_si)) {
    only_in_ex = setdiff(set_ex, set_si)
    only_in_si = setdiff(set_si, set_ex)
    
    if (length(only_in_ex)) msg = sprintf("IDs only in expression: %s", paste(only_in_ex, collapse = ", "))
    if (length(only_in_si)) msg = sprintf("IDs only in sampInfo: %s", paste(only_in_si, collapse = ", "))
    stop(paste(c("Column names of expression matrix do not match sampInfo$ID.", msg), collapse = " "))
  }
  
  # Reorder expression matrix columns to match sampInfo$ID order
  if (!identical(colnames(dat$ex), as.character(sampInfo$ID))) {
    dat$ex = dat$ex[, as.character(sampInfo$ID), drop = FALSE]
  }
  
  # -------------- Package output --------------
  data = list(
    ex         = as.matrix(dat$ex),
    sampInfo   = sampInfo,
    featInfo   = rownames(dat$ex),
    dataname   = dataname
  )
  
  return(data)
}
