# reconstruct.R — published-procedure reconstruction (amendment 1c, section 3).
#
# Comparison analysis only. Does NOT touch code/03d, code/03e, Step 5 wiring, or results/.
# Never sources code/00_helpers.R and never calls cache_or_compute() (avoids the inverted
# DESURV_RECOMPUTE default entirely — see CLAUDE.md).
#
# Three methods computed per input BO set:
#   1. "published (global-best, k replaced)" — the rule's selected k; non-k hyperparameters
#      from the OVERALL best evaluation (standardize_bo_params(bo$overall_best$params), k
#      replaced). This is exactly what code/03_bayesian_optimization.R:115-122 did.
#   2. "conditional-k" — same selected k; hyperparameters from the best evaluation WITHIN
#      that k (the rule's own per-k best row, including its alpha tie-break).
#   3. threshold / best_mean / best_se / k_best / reason, exactly as select_bo_k_by_cv_se()
#      returns them.
#
# Two input modes:
#   - May: cached results/desurv_bo_results_tcgacptac.rds (single desurv_cv_bo_refine object,
#     one search spanning k = 2..12, uneven evaluations per k).
#   - Rev 3 pooled: a directory of per-rank desurv_bo_results_k{K}_tcgacptac.rds objects
#     (K = 2..12), each a desurv_cv_bo object for ONE fixed k, pooled into one
#     history/diagnostics structure the published functions accept unmodified.
#
# R: tested under plain Rscript (DeSurv 1.0.1, site library) for the May object, and under
# R_LIBS=.../lib-rev3 (DeSurv 1.2.0) for the rev 3 rank objects (reading needs the matching
# DeSurv attached for S3 methods on class desurv_cv_bo / desurv_cv_bo_refine, though this
# script itself only touches $history / $diagnostics / $overall_best / $cv_args, which are
# plain data.frames / lists).

`%||%` <- function(a, b) if (is.null(a)) b else a

# published_rule.R lives next to this file. Resolve the directory robustly whether this
# script is run via `Rscript reconstruct.R`, sourced from another script, or sourced
# interactively; fall back to the conventional relative path from the repo root.
.reconstruct_find_dir <- function() {
  frame_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(frame_path)) return(dirname(normalizePath(frame_path)))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  "analysis/published-procedure-reconstruction"
}
source(file.path(.reconstruct_find_dir(), "published_rule.R"))
source(file.path(.reconstruct_find_dir(), "restore_diagnostics.R"))

# =============================================================================
# Adapter: pool per-rank desurv_cv_bo objects (one fixed k each) into a single
# structure that select_bo_k_by_cv_se() / collect_bo_diagnostics() /
# compute_bo_eval_se() accept WITHOUT any change to their bodies.
#
# Why an adapter is needed (read collect_bo_diagnostics() carefully before trusting this):
#   - A per-rank object (class desurv_cv_bo) carries $history (no "k"/"k_grid" column --
#     k is fixed for the whole object, stored only in $cv_args$k_grid) and a top-level
#     $diagnostics data.frame (which DOES carry a "k" column, fold-level). This hits the
#     FIRST branch of collect_bo_diagnostics() ("if (!is.null(bo_results$diagnostics))"),
#     unmodified, no adapter needed for that branch itself.
#   - The May object (class desurv_cv_bo_refine) carries $history WITH a "k_grid" column and
#     NO top-level $diagnostics; it instead has $runs, each with its own $diagnostics. This
#     hits the SECOND branch of collect_bo_diagnostics() ("runs <- bo_results$runs"),
#     unmodified.
#   - So collect_bo_diagnostics() already handles both shapes. The only real pooling problems:
#     (a) history needs a "k" column added per rank (copied from cv_args$k_grid) before
#         rbind, so select_bo_k_by_cv_se()'s k_col <- intersect(c("k","k_grid"), ...) finds it;
#     (b) eval_id is only unique WITHIN a rank (1..40 in every rank object); pooled naively,
#         eval_id=5 in rank k=2 would collide with eval_id=5 in rank k=3, corrupting the
#         left_join()s inside select_bo_k_by_cv_se()/collect_bo_diagnostics() (by = "eval_id").
#   - Fix for (b): offset eval_id by a per-rank constant (rank_index * 100000, ranks are
#     k = 2..12 so 11 ranks never need more than 40 ids each -- 100000 is a safe stride) and
#     apply the IDENTICAL offset to both $history$eval_id and $diagnostics$eval_id for that
#     rank. This changes nothing the SE computation depends on: compute_bo_eval_se() groups
#     exclusively by (eval_id, fold) then by eval_id, and the offset is applied uniformly to
#     every row carrying a given original eval_id, so fold membership and the eval_id lookup
#     key stay perfectly aligned; only the numeric label changes. No run_id is introduced or
#     relied upon (none of the pooled objects have one, and eval_id is now globally unique, so
#     collect_bo_diagnostics()'s run_id fallback path is simply never exercised).
#   - "Inherited vs new" evaluations: in the rank objects, a stage-level artifact of the
#     rank-stratified BO run is that the evaluations inherited from the step 4 main BO are
#     present in $history (with a mean_cindex) but have NO rows in $diagnostics (no fold-level
#     detail was preserved for them). We do not fabricate diagnostics for them; they simply
#     cannot contribute a c_se and collect_bo_diagnostics()/compute_bo_eval_se() leave c_se
#     NA for those eval_ids after the left_join, exactly as the published code would if an
#     evaluation's fold detail were missing. We report, per rank, how many eval_ids have no
#     corresponding diagnostics rows (== inherited) vs how many do (== new).
pool_rank_bo_results <- function(rank_objects) {
  # rank_objects: named list, names are k as character ("2".."12"), values are the
  # desurv_cv_bo objects read from desurv_bo_results_k{K}_tcgacptac.rds.
  stopifnot(length(rank_objects) >= 1)
  ks <- as.integer(names(rank_objects))
  stopifnot(!anyNA(ks), length(ks) == length(unique(ks)))

  offset_for <- function(k) (k - min(ks)) * 100000L

  inherited_report <- list()
  history_list <- list()
  diag_list <- list()

  for (nm in names(rank_objects)) {
    k <- as.integer(nm)
    obj <- rank_objects[[nm]]
    off <- offset_for(k)

    h <- obj$history
    stopifnot(is.data.frame(h), nrow(h) > 0)
    # k is fixed for the whole object; confirm cv_args agrees, then stamp it onto history.
    k_fixed <- obj$cv_args$k_grid
    if (!is.null(k_fixed) && !isTRUE(as.integer(k_fixed) == k)) {
      stop(sprintf(
        "pool_rank_bo_results: rank list name k=%d does not match cv_args$k_grid=%s",
        k, paste(k_fixed, collapse = ",")
      ))
    }
    h$k <- k
    h$eval_id_orig <- h$eval_id
    h$eval_id <- h$eval_id + off
    h$rank_k <- k

    d <- obj$diagnostics
    n_new <- 0L
    n_inherited <- 0L
    if (is.data.frame(d) && nrow(d) > 0) {
      d$eval_id_orig <- d$eval_id
      d$eval_id <- d$eval_id + off
      n_new <- length(unique(d$eval_id_orig))
      diag_list[[nm]] <- d
    }
    n_inherited <- length(setdiff(h$eval_id_orig, if (is.data.frame(d)) unique(d$eval_id_orig) else integer(0)))
    inherited_report[[nm]] <- data.frame(
      k = k,
      n_history_evals = nrow(h),
      n_new_evals_with_diagnostics = n_new,
      n_inherited_evals_without_diagnostics = n_inherited
    )

    history_list[[nm]] <- h
  }

  pooled_history <- do.call(rbind, history_list)
  pooled_diagnostics <- if (length(diag_list)) do.call(rbind, diag_list) else NULL
  inherited_tbl <- do.call(rbind, inherited_report)
  rownames(pooled_history) <- NULL
  if (!is.null(pooled_diagnostics)) rownames(pooled_diagnostics) <- NULL
  rownames(inherited_tbl) <- NULL

  # There is no single tracked "overall_best" for a pool of independent per-rank searches;
  # build it the same way select_bo_k_by_cv_se()/standardize_bo_params() would read it: the
  # pooled history row with the maximum mean_cindex among status=="ok" rows.
  ok <- !is.na(pooled_history$mean_cindex)
  if ("status" %in% names(pooled_history)) ok <- ok & pooled_history$status == "ok"
  ok_hist <- pooled_history[ok, , drop = FALSE]
  best_row <- ok_hist[which.max(ok_hist$mean_cindex), , drop = FALSE]
  param_cols <- intersect(
    c("alpha_grid", "lambda_grid", "nu_grid", "ntop", "lambdaW_grid", "lambdaH_grid", "k_grid", "k"),
    names(best_row)
  )
  overall_best_params <- as.list(best_row[1, param_cols])
  # overall_best$params in the real BO objects is a named numeric vector with a "k_grid"
  # entry (not "k"); match that shape so standardize_bo_params()'s sub("_grid$", "", ...)
  # behaves identically.
  if ("k" %in% names(overall_best_params) && !"k_grid" %in% names(overall_best_params)) {
    names(overall_best_params)[names(overall_best_params) == "k"] <- "k_grid"
  }
  overall_best_params <- unlist(overall_best_params)

  list(
    history = pooled_history,
    diagnostics = pooled_diagnostics,
    overall_best = list(params = overall_best_params, mean_cindex = best_row$mean_cindex[[1]]),
    inherited_report = inherited_tbl
  )
}

# =============================================================================
# Per-k table: best raw C per k, its eval_id, its unpaired SE (c_se), the threshold, and
# whether that k clears the threshold. Duplicates (but does not alter) the grouping /
# tie-break logic internal to select_bo_k_by_cv_se() because that function returns only the
# scalar winning-k summary, not the full per-k table the report needs. We verify below (in
# the May validation run) that this reproduces identical k_selected / k_best / best_mean /
# best_se / threshold, so the duplication is checked, not just asserted.
compute_per_k_table <- function(bo_results) {
  history_df <- bo_results$history
  ok <- !is.na(history_df$mean_cindex)
  if ("status" %in% names(history_df)) ok <- ok & history_df$status == "ok"
  history_df <- history_df[ok, , drop = FALSE]

  diagnostics <- collect_bo_diagnostics(bo_results, history_df)
  eval_se <- compute_bo_eval_se(diagnostics)
  if (is.data.frame(eval_se) && nrow(eval_se)) {
    join_cols <- intersect(c("run_id", "eval_id"), names(history_df))
    join_cols <- join_cols[join_cols %in% names(eval_se)]
    if (!length(join_cols)) join_cols <- "eval_id"
    history_df <- dplyr::left_join(history_df, eval_se, by = join_cols)
  }
  if (!"c_se" %in% names(history_df)) history_df$c_se <- NA_real_

  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  alpha_col <- intersect(c("alpha", "alpha_grid"), names(history_df))[1]
  groups <- split(seq_len(nrow(history_df)), history_df[[k_col]])
  best_idx <- vapply(groups, function(idx) {
    if (is.na(alpha_col)) {
      idx[order(-history_df$mean_cindex[idx])[1]]
    } else {
      idx[order(-history_df$mean_cindex[idx], history_df[[alpha_col]][idx])[1]]
    }
  }, integer(1))
  best_per_k <- history_df[best_idx, , drop = FALSE]
  best_per_k$k_numeric <- suppressWarnings(as.numeric(as.character(best_per_k[[k_col]])))
  best_per_k <- best_per_k[order(best_per_k$k_numeric), , drop = FALSE]

  k_sel <- select_bo_k_by_cv_se(bo_results)
  best_per_k$clears_threshold <- is.finite(best_per_k$mean_cindex) &
    is.finite(k_sel$lcb_threshold) &
    best_per_k$mean_cindex >= k_sel$lcb_threshold

  list(per_k = best_per_k, selection = k_sel, k_col = k_col, alpha_col = alpha_col)
}

# =============================================================================
# conditional-k parameters: non-k hyperparameters from the best row WITHIN the selected k
# (the per-k table's row for k_selected), rather than the global-best row.
extract_conditional_k_params <- function(per_k_table, k_selected) {
  row <- per_k_table[per_k_table$k_numeric == k_selected, , drop = FALSE]
  if (!nrow(row)) stop("extract_conditional_k_params: no per-k row for k_selected=", k_selected)
  row <- row[1, , drop = FALSE]
  param_cols <- intersect(
    c("alpha_grid", "alpha", "lambda_grid", "lambda", "nu_grid", "nu", "ntop",
      "lambdaW_grid", "lambdaH_grid"),
    names(row)
  )
  out <- as.list(row[1, param_cols])
  names(out) <- sub("_grid$", "", names(out))
  out <- out[!duplicated(names(out))]
  out$k <- as.integer(k_selected)
  out$eval_id <- if ("eval_id_orig" %in% names(row)) row$eval_id_orig[[1]] else row$eval_id[[1]]
  out
}

# =============================================================================
# Top-level driver: run both methods + the raw selection on one bo_results object (either the
# May desurv_cv_bo_refine object, or a pool_rank_bo_results() output).
run_reconstruction <- function(bo_results, label, require_finite_best_se = FALSE) {
  if (isTRUE(require_finite_best_se)) {
    # FIX 3: never let the rev 3 reconstruction silently fall back to "pick the best rank"
    # because the pooled overall-best row happened to be an inherited evaluation with no SE.
    # The May path does not set this flag and is unaffected.
    assert_overall_best_has_finite_se(bo_results)
  }
  k_sel <- select_bo_k_by_cv_se(bo_results)
  per_k_info <- compute_per_k_table(bo_results)
  stopifnot(identical(k_sel$k_selected, per_k_info$selection$k_selected))
  stopifnot(identical(k_sel$k_best, per_k_info$selection$k_best))
  stopifnot(isTRUE(all.equal(k_sel$best_mean, per_k_info$selection$best_mean)))
  stopifnot(isTRUE(all.equal(k_sel$best_se, per_k_info$selection$best_se)))
  stopifnot(isTRUE(all.equal(k_sel$lcb_threshold, per_k_info$selection$lcb_threshold)))

  published_params <- extract_published_params(bo_results, k_sel)
  conditional_params <- if (!is.null(k_sel$k_selected)) {
    extract_conditional_k_params(per_k_info$per_k, k_sel$k_selected)
  } else {
    NULL
  }

  list(
    label = label,
    selection = k_sel,
    per_k_table = per_k_info$per_k,
    published_global_best_params = published_params,
    conditional_k_params = conditional_params
  )
}

write_reconstruction <- function(result, out_dir, stem) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  saveRDS(result, file.path(out_dir, paste0(stem, ".rds")))
  summary_row <- data.frame(
    label = result$label,
    k_selected = result$selection$k_selected %||% NA_integer_,
    k_best = result$selection$k_best %||% NA_integer_,
    best_mean = result$selection$best_mean,
    best_se = result$selection$best_se,
    lcb_threshold = result$selection$lcb_threshold,
    reason = result$selection$reason
  )
  write.csv(summary_row, file.path(out_dir, paste0(stem, "_selection.csv")), row.names = FALSE)
  write.csv(result$per_k_table, file.path(out_dir, paste0(stem, "_per_k.csv")), row.names = FALSE)
  write.csv(
    as.data.frame(result$published_global_best_params),
    file.path(out_dir, paste0(stem, "_published_global_best_params.csv")),
    row.names = FALSE
  )
  if (!is.null(result$conditional_k_params)) {
    write.csv(
      as.data.frame(result$conditional_k_params),
      file.path(out_dir, paste0(stem, "_conditional_k_params.csv")),
      row.names = FALSE
    )
  }
  invisible(file.path(out_dir, stem))
}

# =============================================================================
# Mode 1: May (cached results/desurv_bo_results_tcgacptac.rds, read-only).
run_may_validation <- function(results_dir, out_dir) {
  bo <- readRDS(file.path(results_dir, "desurv_bo_results_tcgacptac.rds"))
  result <- run_reconstruction(bo, "May historical (cached desurv_bo_results_tcgacptac)")
  write_reconstruction(result, out_dir, "may_published_procedure_reconstruction")
  result
}

# =============================================================================
# Mode 2: rev 3 pooled (directory of desurv_bo_results_k{K}_tcgacptac.rds, K = 2..12).
#
# main_bo_path: step4_main_desurv_bo_results_tcgacptac.rds -- FIX 1, restores diagnostics for
#   evaluations each rank inherited from the step 4 main BO. Required: the rev 3 reconstruction
#   must not run with unrestored inherited evaluations in the pool (FIX 3 would stop it anyway
#   if the overall best happened to be one of them, but we always restore first so a *different*
#   best-row choice doesn't quietly dodge the check).
# noise_results_path: optional seed-bank noise job .rds (e.g. noise_results_3500640.rds) --
#   FIX 2, restores the two historically-relevant k=3 configs from their own re-run at the BO's
#   seed bank, if a rank's inherited row matches one exactly.
run_rev3_pooled <- function(input_dir, out_dir, ks = 2:12, main_bo_path, noise_results_path = NULL) {
  files <- file.path(input_dir, sprintf("desurv_bo_results_k%d_tcgacptac.rds", ks))
  present <- file.exists(files)
  if (!all(present)) {
    stop(
      "run_rev3_pooled: missing rank objects for k = ",
      paste(ks[!present], collapse = ", "),
      " in ", input_dir
    )
  }
  stopifnot(!missing(main_bo_path), file.exists(main_bo_path))
  main_bo <- readRDS(main_bo_path)
  noise_results <- if (!is.null(noise_results_path)) readRDS(noise_results_path) else NULL

  rank_objects <- setNames(lapply(files, readRDS), as.character(ks))
  restore_reports <- list()
  for (nm in names(rank_objects)) {
    k <- as.integer(nm)
    obj <- rank_objects[[nm]]
    fix1 <- restore_inherited_diagnostics_from_main(k, obj$history, obj$diagnostics, main_bo)
    obj$diagnostics <- fix1$diagnostics
    fix2 <- if (!is.null(noise_results)) {
      restore_noise_bank_diagnostics(k, obj$history, obj$diagnostics, noise_results)
    } else {
      list(diagnostics = obj$diagnostics, report = data.frame())
    }
    obj$diagnostics <- fix2$diagnostics
    rank_objects[[nm]] <- obj
    restore_reports[[nm]] <- list(k = k, fix1_report = fix1$report, fix2_report = fix2$report)
  }

  pooled <- pool_rank_bo_results(rank_objects)
  result <- run_reconstruction(
    pooled, "Rev 3 pooled published-procedure reconstruction (k2-k12)",
    require_finite_best_se = TRUE
  )
  result$inherited_report <- pooled$inherited_report
  result$restore_reports <- restore_reports

  overall_best_eval <- pooled$history$eval_id[which.max(ifelse(
    !is.na(pooled$history$mean_cindex) &
      (if ("status" %in% names(pooled$history)) pooled$history$status == "ok" else TRUE),
    pooled$history$mean_cindex, -Inf
  ))]
  result$may_comparable_se <- may_comparable_se_table(
    pooled, result$per_k_table, overall_best_eval, result$selection$best_mean
  )

  stem_path <- write_reconstruction(result, out_dir, "rev3_published_procedure_reconstruction")
  write.csv(pooled$inherited_report, paste0(stem_path, "_inherited_report.csv"), row.names = FALSE)
  write.csv(result$may_comparable_se$per_rep, paste0(stem_path, "_may_comparable_se.csv"), row.names = FALSE)
  result
}
