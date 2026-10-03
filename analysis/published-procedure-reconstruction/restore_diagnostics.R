# restore_diagnostics.R — adapter-only fixes (owner-approved follow-up, 2026-10-03).
#
# Comparison analysis only. The published functions in published_rule.R are NOT touched by
# anything in this file. No results/ writes, no code/00_helpers.R, no cache_or_compute().
#
# FIX 1: rank objects carry history rows "inherited" from the step 4 main BO
# (step4_main_desurv_bo_results_tcgacptac.rds, md5 ad25d5342a77572f71c5dd5e15182a27) with no
# fold-level diagnostics of their own. If the pooled overall-best evaluation happened to be one
# of these, select_bo_k_by_cv_se() would see a non-finite c_se and silently fall back to
# "pick the best rank" (reason = "missing_se"), which is an assembly artifact of the pooling,
# not a property of the rev 3 BO surface. This file restores the missing fold diagnostics for
# inherited rows by matching them to the exact main-BO evaluation that produced them (same
# recipe, partitions, seed bank) and transplanting its fold-level rows, under hard checks.
#
# FIX 2: the two historically-relevant k = 3 configurations (published config; May eval 136
# config) were independently re-run at seed bank 4,000,000 (the BO's own bank) in the seed-bank
# noise job (job 3500640, analysis/longleaf/bench/noise_bo_eval.R). If a rank's inherited row
# matches one of those configs exactly, transplant fold C from there instead, under the same
# hard check. If no matching inherited row exists (e.g. no k=3 rank object is available yet),
# this is a no-op that reports what the noise file supports.
#
# FIX 3: a guard that stops the rev 3 reconstruction outright if the pooled overall-best row
# (as select_bo_k_by_cv_se() itself would pick it) still lacks a finite c_se after FIX 1 / FIX
# 2. The May path is untouched by this guard.

this_dir <- "analysis/published-procedure-reconstruction"
if (!exists("select_bo_k_by_cv_se")) source(file.path(this_dir, "published_rule.R"))

# ---- FIX 1: restore inherited diagnostics from the step 4 main BO --------------------------

#' Identify rank-history rows with no fold diagnostics of their own.
inherited_rows <- function(rank_history, rank_diagnostics) {
  diag_ids <- if (is.data.frame(rank_diagnostics) && nrow(rank_diagnostics)) {
    unique(rank_diagnostics$eval_id)
  } else {
    integer(0)
  }
  rank_history[!(rank_history$eval_id %in% diag_ids), , drop = FALSE]
}

#' Exact-value match of one inherited row to the step 4 main BO history (k + alpha/lambda/nu/ntop).
#' HARD CHECK: exactly one match, or stop().
match_main_bo_eval <- function(row, main_history) {
  hit <- main_history$k_grid == row$k &
    main_history$alpha_grid == row$alpha_grid &
    main_history$lambda_grid == row$lambda_grid &
    main_history$nu_grid == row$nu_grid &
    main_history$ntop == row$ntop
  hit[is.na(hit)] <- FALSE
  idx <- which(hit)
  if (length(idx) != 1L) {
    stop(sprintf(
      "match_main_bo_eval: expected exactly one main-BO match for inherited eval_id=%s (k=%s), found %d",
      row$eval_id, row$k, length(idx)
    ))
  }
  main_history[idx, , drop = FALSE]
}

#' Recompute the fold-aggregated mean exactly as compute_bo_eval_se() would (group by fold,
#' mean within fold across every row regardless of rep/init_id, then mean across folds).
recompute_mean_from_fold_rows <- function(diag_rows) {
  fm <- tapply(diag_rows$val_cindex, diag_rows$fold, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    mean(x)
  })
  mean(fm, na.rm = TRUE)
}

#' Restore diagnostics for a rank's inherited evaluations from the step 4 main BO object.
#' `rank_history`/`rank_diagnostics` use the rank's OWN (un-offset) eval_id; offsetting for
#' pooling happens afterwards in pool_rank_bo_results(), unaffected by this step.
#' Returns the augmented diagnostics data.frame (original rows + transplanted rows), with two
#' new provenance columns (source, prov_id) added to every row (NA for pre-existing rows).
restore_inherited_diagnostics_from_main <- function(rank_k, rank_history, rank_diagnostics, main_bo) {
  main_history <- main_bo$history
  # Limitation (noted, not generalized; owner 2026-10-03): reads runs[[1]] only. The step 4 main BO object
  # (md5 ad25d534) has exactly one run holding all 150 evaluations' diagnostics.
  main_diag <- main_bo$runs[[1]]$diagnostics
  stopifnot(is.data.frame(main_history), is.data.frame(main_diag))

  base <- if (is.data.frame(rank_diagnostics) && nrow(rank_diagnostics)) rank_diagnostics else {
    rank_diagnostics[0, , drop = FALSE]
  }
  if (!"source" %in% names(base)) base$source <- NA_character_
  if (!"prov_id" %in% names(base)) base$prov_id <- NA_character_

  # (rank_history for a per-rank object has no "k" column of its own at this point in the
  #  pipeline -- it is added later by pool_rank_bo_results(); every inherited row belongs to
  #  this rank's fixed k, which is why rank_k is passed in explicitly and stamped onto each
  #  row inside the loop below rather than filtered here.)
  inh <- inherited_rows(rank_history, rank_diagnostics)

  report <- data.frame(
    eval_id = integer(0), main_eval_id = integer(0),
    history_mean = numeric(0), main_history_mean = numeric(0),
    recomputed_mean = numeric(0), n_fold_rows = integer(0), all_checks_pass = logical(0)
  )
  if (!nrow(inh)) {
    return(list(diagnostics = base, report = report))
  }

  transplants <- list()
  for (i in seq_len(nrow(inh))) {
    row <- inh[i, , drop = FALSE]
    row$k <- rank_k
    main_match <- match_main_bo_eval(row, main_history)

    # HARD CHECK 1: the matched main eval's mean_cindex equals the rank-history mean_cindex exactly.
    if (!isTRUE(main_match$mean_cindex[[1]] == row$mean_cindex[[1]])) {
      stop(sprintf(
        "restore_inherited_diagnostics_from_main: mean_cindex mismatch for inherited eval_id=%s (k=%s): rank history %.15f vs main-BO eval_id=%s %.15f",
        row$eval_id, rank_k, row$mean_cindex[[1]], main_match$eval_id[[1]], main_match$mean_cindex[[1]]
      ))
    }

    fold_rows <- main_diag[main_diag$eval_id == main_match$eval_id[[1]], , drop = FALSE]
    if (!nrow(fold_rows)) {
      stop(sprintf(
        "restore_inherited_diagnostics_from_main: no fold diagnostics found in main BO for matched eval_id=%s",
        main_match$eval_id[[1]]
      ))
    }

    # HARD CHECK 2: mean recomputed from the transplanted fold rows (same aggregation
    # compute_bo_eval_se() uses: mean within fold, then mean across folds) equals mean_cindex
    # exactly (tolerance 1e-12).
    recomputed <- recompute_mean_from_fold_rows(fold_rows)
    ok <- isTRUE(all.equal(recomputed, row$mean_cindex[[1]], tolerance = 1e-12))
    if (!ok) {
      stop(sprintf(
        "restore_inherited_diagnostics_from_main: recomputed fold mean %.15f does not match history mean_cindex %.15f for eval_id=%s (main eval_id=%s)",
        recomputed, row$mean_cindex[[1]], row$eval_id, main_match$eval_id[[1]]
      ))
    }

    fold_rows$eval_id <- row$eval_id  # relabel to the rank's own eval_id (pre-offset)
    fold_rows$source <- "step4_main"
    fold_rows$prov_id <- as.character(main_match$eval_id[[1]])
    transplants[[length(transplants) + 1L]] <- fold_rows

    report <- rbind(report, data.frame(
      eval_id = row$eval_id, main_eval_id = main_match$eval_id[[1]],
      history_mean = row$mean_cindex[[1]], main_history_mean = main_match$mean_cindex[[1]],
      recomputed_mean = recomputed, n_fold_rows = nrow(fold_rows), all_checks_pass = TRUE
    ))
  }

  all_cols <- union(names(base), names(transplants[[1]]))
  for (nm in setdiff(all_cols, names(base))) base[[nm]] <- NA
  for (j in seq_along(transplants)) {
    for (nm in setdiff(all_cols, names(transplants[[j]]))) transplants[[j]][[nm]] <- NA
    transplants[[j]] <- transplants[[j]][, all_cols, drop = FALSE]
  }
  base <- base[, all_cols, drop = FALSE]
  out <- rbind(base, do.call(rbind, transplants))
  rownames(out) <- NULL
  list(diagnostics = out, report = report)
}

# ---- FIX 2: restore from the seed-bank noise job (job 3500640) ----------------------------

#' Extract fold-level C for every config x bank in a noise_bo_eval.R output object, in
#' diagnostics-schema form (val_cindex instead of cindex; label/bank kept as provenance).
noise_results_to_diag_table <- function(noise_results) {
  rows <- lapply(noise_results$out, function(o) {
    if (!isTRUE(o$ok)) return(NULL)
    cr <- o$res$cv_results
    cr$val_cindex <- cr$cindex
    cr$noise_label <- as.character(o$cfg$label)
    cr$noise_bank <- o$cfg$bank
    cr$ntop <- o$cfg$ntop[[1]]  # cv_results itself carries no ntop column; it's a preprocessing arg
    cr$reported_mean_cindex <- o$res$summary$mean_cindex[[1]]
    cr
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  do.call(rbind, rows)
}

#' Restore diagnostics for inherited rows that exactly match a noise-job config at the given
#' bank (default 4,000,000, the BO's own bank). No-op (with a report) if no inherited row in
#' this rank matches any noise config -- does NOT fabricate anything.
restore_noise_bank_diagnostics <- function(rank_k, rank_history, rank_diagnostics, noise_results, bank = 4000000L) {
  noise_tbl <- noise_results_to_diag_table(noise_results)
  noise_tbl <- noise_tbl[noise_tbl$noise_bank == bank, , drop = FALSE]

  base <- if (is.data.frame(rank_diagnostics) && nrow(rank_diagnostics)) rank_diagnostics else {
    rank_diagnostics[0, , drop = FALSE]
  }
  if (!"source" %in% names(base)) base$source <- NA_character_
  if (!"prov_id" %in% names(base)) base$prov_id <- NA_character_

  inh <- inherited_rows(rank_history, base)
  inh$k <- rank_k

  report <- data.frame(
    eval_id = integer(0), noise_label = character(0), history_mean = numeric(0),
    recomputed_mean = numeric(0), matched = logical(0), transplanted = logical(0)
  )
  if (!nrow(inh) || !nrow(noise_tbl)) {
    return(list(diagnostics = base, report = report))
  }

  transplants <- list()
  for (i in seq_len(nrow(inh))) {
    row <- inh[i, , drop = FALSE]
    labels_here <- unique(noise_tbl$noise_label)
    matched_label <- NA_character_
    for (lbl in labels_here) {
      cfg_rows <- noise_tbl[noise_tbl$noise_label == lbl, , drop = FALSE]
      cfg <- cfg_rows[1, , drop = FALSE]
      if (isTRUE(cfg$k[[1]] == row$k[[1]]) &&
            isTRUE(cfg$alpha[[1]] == row$alpha_grid[[1]]) &&
            isTRUE(cfg$lambda[[1]] == row$lambda_grid[[1]]) &&
            isTRUE(cfg$nu[[1]] == row$nu_grid[[1]]) &&
            isTRUE(cfg$ntop[[1]] == row$ntop[[1]])) {
        matched_label <- lbl
        break
      }
    }
    if (is.na(matched_label)) next

    fold_rows <- noise_tbl[noise_tbl$noise_label == matched_label, , drop = FALSE]
    recomputed <- recompute_mean_from_fold_rows(fold_rows)
    ok <- isTRUE(all.equal(recomputed, row$mean_cindex[[1]], tolerance = 1e-12))
    if (!ok) {
      report <- rbind(report, data.frame(
        eval_id = row$eval_id, noise_label = matched_label, history_mean = row$mean_cindex[[1]],
        recomputed_mean = recomputed, matched = TRUE, transplanted = FALSE
      ))
      next  # do NOT manufacture an SE; leave NA.
    }
    fold_rows$eval_id <- row$eval_id
    fold_rows$source <- paste0("noise_3500640_", matched_label)
    fold_rows$prov_id <- matched_label
    transplants[[length(transplants) + 1L]] <- fold_rows
    report <- rbind(report, data.frame(
      eval_id = row$eval_id, noise_label = matched_label, history_mean = row$mean_cindex[[1]],
      recomputed_mean = recomputed, matched = TRUE, transplanted = TRUE
    ))
  }

  if (!length(transplants)) {
    return(list(diagnostics = base, report = report))
  }
  all_cols <- union(names(base), unique(unlist(lapply(transplants, names))))
  for (nm in setdiff(all_cols, names(base))) base[[nm]] <- NA
  base <- base[, all_cols, drop = FALSE]
  for (j in seq_along(transplants)) {
    for (nm in setdiff(all_cols, names(transplants[[j]]))) transplants[[j]][[nm]] <- NA
    transplants[[j]] <- transplants[[j]][, all_cols, drop = FALSE]
  }
  out <- rbind(base, do.call(rbind, transplants))
  rownames(out) <- NULL
  list(diagnostics = out, report = report)
}

# ---- FIX 3: guard against the missing_se fallback being reached silently ------------------

#' stop() if the pooled overall-best row (as select_bo_k_by_cv_se() itself would pick it) has
#' no finite c_se. Intended to be called on the rev 3 pool only, after FIX 1 / FIX 2 have been
#' applied; the May path never calls this.
assert_overall_best_has_finite_se <- function(bo_results) {
  k_sel <- select_bo_k_by_cv_se(bo_results)
  if (identical(k_sel$reason, "missing_se") || !is.finite(k_sel$best_se %||% NA_real_)) {
    stop(
      "assert_overall_best_has_finite_se: the rev 3 overall-best evaluation has no finite ",
      "c_se (reason='", k_sel$reason, "'). This is almost certainly an inherited evaluation ",
      "missing fold diagnostics (an assembly artifact), not a result. Restore its diagnostics ",
      "(FIX 1 / FIX 2) before calling the published rule, or report the gap -- do not let the ",
      "reconstruction silently fall back to 'pick the best rank'."
    )
  }
  invisible(k_sel)
}

# ---- Descriptive-only diagnostic: "May-comparable SE" (never a selection rule) -------------

#' For one evaluation's diagnostics rows, compute compute_bo_eval_se()'s SE formula separately
#' within each of the 3 rev-3 partitions (the "rep" column), instead of pooling all 15
#' partition x fold units into one SE as compute_bo_eval_se() does. Descriptive only: this is
#' NOT a selection rule and must be labelled as such everywhere it is reported.
compute_may_comparable_se <- function(diag_rows) {
  stopifnot(all(c("rep", "fold", "val_cindex") %in% names(diag_rows)))
  reps <- sort(unique(diag_rows$rep))
  rows <- lapply(reps, function(r) {
    sub <- diag_rows[diag_rows$rep == r, , drop = FALSE]
    fm <- tapply(sub$val_cindex, sub$fold, function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      mean(x)
    })
    fm <- fm[is.finite(fm)]
    se <- if (length(fm) <= 1L) NA_real_ else stats::sd(fm) / sqrt(length(fm))
    data.frame(rep = r, se_5fold = se, n_folds = length(fm), fold_mean_of_means = mean(fm))
  })
  out <- do.call(rbind, rows)
  out$label <- "May-comparable SE (descriptive only -- NOT a selection rule)"
  out
}

#' Build the full descriptive table: per-rep SE/threshold, plus which per-k best rows would
#' clear each of the 3 single-partition thresholds.
may_comparable_se_table <- function(bo_results, per_k_table, eval_id, best_mean) {
  diagnostics <- collect_bo_diagnostics(bo_results, bo_results$history)
  diag_rows <- diagnostics[diagnostics$eval_id == eval_id, , drop = FALSE]
  se_tbl <- compute_may_comparable_se(diag_rows)
  se_tbl$lcb_threshold_5fold <- best_mean - se_tbl$se_5fold
  clears <- lapply(se_tbl$lcb_threshold_5fold, function(th) {
    per_k_table$k_numeric[is.finite(per_k_table$mean_cindex) & per_k_table$mean_cindex >= th]
  })
  list(
    per_rep = se_tbl,
    mean_se_5fold = mean(se_tbl$se_5fold, na.rm = TRUE),
    range_se_5fold = range(se_tbl$se_5fold, na.rm = TRUE),
    clears_by_rep = clears,
    note = "Descriptive only. May SE came from one 5-fold partition; rev 3 pools 3 partitions x 5 folds = 15 units. This table isolates whether a k change comes from the mean curve or from SE tightening; it is never used to select k."
  )
}
