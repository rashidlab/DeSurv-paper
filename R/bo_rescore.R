# R/bo_rescore.R: fresh-partition re-evaluation of the rank-BO finalists (rev 3 amendment 1b,
# docs/plans/2026-10-03-rev3-amendment-1b-fresh-rescoring.md). Sourced by code/03d.
#
# rescore_finalists() repeats the per-finalist call of DeSurv::desurv_bo_finalists()
# (DeSurv 1.2.0, R/desurv_bo_finalists.R:176-193: .desurv_merge_args(data_args, cv_args,
# cv_overrides, point), cv_only = TRUE, desurv_cv()) with the cv_args stored in the finalists
# object, replacing only cv_seeds and seed_bank_base. The package's cv_overrides refuses those
# two, hence this copy. No cv_overrides here: parallel_grid / ncores_grid are the stored ones
# (the BO job's, which desurv_bo_finalists() also used). With the stored seeds this reproduces
# the finalists' mean_cindex exactly (quick-mode check, docs/plans/2026-10-03-amendment1b-impl-report.md).
# Defines functions only; reads and writes no files.

# fin: a desurv_bo_finalists object (one rank); tar_data: the training data the BO used.
# Returns list(configs, fold_cindex, diagnostics, cv_args): configs gain mean_cindex_out
# (mean held-out C over partition x fold) and train_cindex_out; fold_cindex is long
# (finalist_id, eval_id, k, rep, fold, cindex); diagnostics keeps the per-fold seed record.
rescore_finalists <- function(fin, tar_data, cv_seeds, seed_bank_base, verbose = TRUE) {
  cv_args <- fin$cv_args
  cv_args$cv_seeds <- as.integer(cv_seeds)
  cv_args$seed_bank_base <- as.integer(seed_bank_base)
  data_args <- list(X = tar_data$ex, y = tar_data$sampInfo$time, d = tar_data$sampInfo$event,
                    dataset = tar_data$sampInfo$dataset, samp_keeps = tar_data$samp_keeps)
  configs <- fin$configs
  configs$mean_cindex_out <- NA_real_
  configs$train_cindex_out <- NA_real_
  fold_rows <- diag_rows <- vector("list", nrow(configs))
  for (i in seq_len(nrow(configs))) {
    point <- as.list(configs[i, fin$param_names, drop = FALSE])
    args <- DeSurv:::.desurv_merge_args(data_args, cv_args, list(), point)
    args$cv_only <- TRUE
    if (verbose) {
      message(sprintf("  finalist %d/%d (eval %d, k = %d), cv_seeds %s, seed_bank_base %d",
                      i, nrow(configs), configs$eval_id[i], configs$k[i],
                      paste(cv_args$cv_seeds, collapse = ","), cv_args$seed_bank_base))
    }
    res <- do.call(DeSurv::desurv_cv, args)
    sf <- res$summary_fold
    if (!"rep" %in% names(sf)) sf$rep <- 1L
    fold_rows[[i]] <- data.frame(
      finalist_id = configs$finalist_id[i], eval_id = configs$eval_id[i],
      k = configs$k[i], rep = sf$rep, fold = sf$fold,
      cindex = sf$mean_cindex_fold, row.names = NULL
    )
    configs$mean_cindex_out[i] <- res$summary$mean_cindex[1L]
    dg <- res$diagnostics
    if (is.null(dg) || !all(c("cv_seed", "seed_first", "rep", "fold") %in% names(dg))) {
      stop("desurv_cv() diagnostics lack cv_seed / seed_first / rep / fold; ",
           "cannot verify which partitions and seed banks were used.")
    }
    diag_rows[[i]] <- data.frame(
      finalist_id = configs$finalist_id[i], rep = dg$rep, fold = dg$fold,
      cv_seed = dg$cv_seed, seed_first = dg$seed_first, seed_last = dg$seed_last,
      n_seed_fits = dg$n_seed_fits, train_cindex = as.numeric(dg$train_cindex),
      val_cindex = as.numeric(dg$val_cindex), status = dg$status, row.names = NULL
    )
    configs$train_cindex_out[i] <- mean(diag_rows[[i]]$train_cindex, na.rm = TRUE)
  }
  list(configs = configs, fold_cindex = do.call(rbind, fold_rows),
       diagnostics = do.call(rbind, diag_rows), cv_args = cv_args)
}
