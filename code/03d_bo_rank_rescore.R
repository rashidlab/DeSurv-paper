#!/usr/bin/env Rscript
# code/03d_bo_rank_rescore.R — fresh-partition rescoring of the rank-BO finalists and the
# cross-rank selection on it (rev 3 amendment 1b; role per amendment 1c: robustness / sensitivity
# analysis of rank selection, never a Step 5 input)
#
# docs/plans/2026-10-03-rev3-amendment-1b-fresh-rescoring.md and
# docs/plans/2026-10-03-rev3-amendment-1c-roles-and-reconstruction.md (APPROVED). desurv_bo_finalists()
# re-evaluates on the BO's own partitions and seed banks, so its C equals the BO value and the
# cross-rank 1-SE rule compares selection-biased maxima. Here every finalist of code/03b is
# re-evaluated with the identical estimator and settings (its stored cv_args) on FRESH CV:
# partition seeds 20261111-20261113 and per-fold seed banks from 5,000,000. Nothing else changes.
#
# Two modes:
#   rank   DESURV_BO_K=<K>             one rank's finalists -> desurv_bo_rescore_k{K}_tcgacptac
#   select DESURV_BO_RESCORE_SELECT=1  paired 1-SE across ranks (DESURV_BO_RANKS, default 2:12,
#                                      parsed as in code/03c) on the fresh fold C
#
# Rank inputs:   tar_data_tcgacptac, desurv_bo_finalists_k{K}_tcgacptac,
#                desurv_bo_results_k{K}_tcgacptac (data fingerprint only)
# Rank output:   desurv_bo_rescore_k{K}_tcgacptac   list(configs, fold_cindex, diagnostics,
#                                                   cv_args, cv_args_bo, k): configs keep the
#                                                   BO-partition mean_cindex and gain mean_cindex_fresh
# Select inputs: tar_data_tcgacptac and, per K, desurv_bo_rescore_k{K}, desurv_bo_finalists_k{K},
#                desurv_bo_results_k{K}
# Select outputs (all _fresh names, so none can be mistaken for a pipeline selection; amendment 1c):
#   desurv_bo_finalists_fresh_tcgacptac pooled desurv_bo_finalists object on the FRESH fold C
#   tar_k_selection_fresh_tcgacptac     k_selection_from_finalists()   } attr selection_rule =
#   tar_params_best_fresh_tcgacptac     select_bo_params_finalists()   } "fresh_partition_paired_1se"
#   tar_data_filtered_fresh_tcgacptac   as code/03c (at B's selected ngene)
#   bo_rank_selection_side_by_side_tcgacptac (.rds + .csv, + _finalists.csv): per rank, B (fresh,
#                                       robustness / sensitivity), A (BO partitions), C (best GP mean), B minus A

message("=== Step 3d: fresh-partition rescoring of the rank-BO finalists (rev 3 amendment 1b) ===")
source("code/00_helpers.R")
library(DeSurv)
source("R/bo_helpers.R")
source("R/reported_settings.R")
source("R/bo_setup.R")
source("R/preprocess_helpers.R")
source("R/bo_rescore.R")

# ── Fresh CV (amendment 1b) ────────────────────────────────────────────────
# BO: partitions 20261101-3, banks from 4,000,000; sensitivity reruns 3,000,000+;
# 20261201-2 used elsewhere. Quick mode: one partition, as R/bo_setup.R's quick branch.
if (CONFIG$quick) {
  FRESH_CV_SEEDS <- 20261111L
} else {
  FRESH_CV_SEEDS <- c(20261111L, 20261112L, 20261113L)
}
FRESH_SEED_BANK_BASE <- 5000000L
SELECTION_RULE <- "fresh_partition_paired_1se"

SELECT_MODE <- identical(Sys.getenv("DESURV_BO_RESCORE_SELECT"), "1")
if (SELECT_MODE && nzchar(Sys.getenv("DESURV_BO_K"))) {
  stop("Set either DESURV_BO_K (rank mode) or DESURV_BO_RESCORE_SELECT=1 (select mode), not both.")
}

tar_data <- load_precomputed("tar_data_tcgacptac")

if (!SELECT_MODE) {
  # ════════════════════════════ rank mode ════════════════════════════════
  BO_K <- as.integer(Sys.getenv("DESURV_BO_K", NA))
  if (is.na(BO_K) || BO_K < 2L || BO_K > BO_K_UPPER) {
    stop("Set DESURV_BO_K to an integer in 2..", BO_K_UPPER, " (got '",
         Sys.getenv("DESURV_BO_K"), "'), or DESURV_BO_RESCORE_SELECT=1.")
  }
  message("  k = ", BO_K, "; fresh cv_seeds ", paste(FRESH_CV_SEEDS, collapse = ","),
          ", seed_bank_base ", FRESH_SEED_BANK_BASE)
  fin <- load_precomputed(sprintf("desurv_bo_finalists_k%d_tcgacptac", BO_K))
  bo  <- load_precomputed(sprintf("desurv_bo_results_k%d_tcgacptac", BO_K))
  stopifnot(inherits(fin, "desurv_bo_finalists"), all(fin$configs$k == BO_K),
            as.integer(fin$cv_args$k_grid) == BO_K,
            all(c("cv_seeds", "seed_bank_base") %in% names(fin$cv_args)))
  if (length(intersect(fin$cv_args$cv_seeds, FRESH_CV_SEEDS)) ||
      fin$cv_args$seed_bank_base == FRESH_SEED_BANK_BASE) {
    stop("The fresh partitions or seed bank coincide with the BO's.")
  }
  # The data must be the data the BO (and so the finalists) were run on.
  fp <- DeSurv:::.desurv_bo_data_fingerprint(tar_data$ex, tar_data$sampInfo$time,
                                             tar_data$sampInfo$event, tar_data$sampInfo$dataset,
                                             tar_data$samp_keeps)
  if (!isTRUE(all.equal(fp, bo$data_fingerprint, tolerance = 1e-12))) {
    stop("tar_data differs from the data the rank BO was run on.")
  }

  rescore <- cache_or_compute(sprintf("desurv_bo_rescore_k%d_tcgacptac", BO_K), {
    rs <- rescore_finalists(fin, tar_data, FRESH_CV_SEEDS, FRESH_SEED_BANK_BASE)
    configs <- fin$configs
    configs$mean_cindex_fresh <- rs$configs$mean_cindex_out
    configs$delta_fresh_vs_bo <- configs$mean_cindex_fresh - configs$mean_cindex
    configs$train_cindex_fresh <- rs$configs$train_cindex_out
    list(configs = configs, fold_cindex = rs$fold_cindex, diagnostics = rs$diagnostics,
         cv_args = rs$cv_args, cv_args_bo = fin$cv_args, param_names = fin$param_names,
         k = BO_K, fresh_cv_seeds = FRESH_CV_SEEDS, fresh_seed_bank_base = FRESH_SEED_BANK_BASE)
  })

  # ── HARD CHECKS ──
  # (a) The fresh partitions and seed banks were really used: every diagnostic row's cv_seed
  #     is a fresh seed (and every fresh seed occurs), every seed_first is in the 5,000,000
  #     bank, and equals that bank's first seed for its partition x fold.
  dg <- rescore$diagnostics
  nf <- as.integer(rescore$cv_args$nfolds)
  stride <- if (is.null(rescore$cv_args$seed_bank_stride)) 10000L else rescore$cv_args$seed_bank_stride  # desurv_cv() default
  expect_first <- vapply(seq_len(nrow(dg)), function(i) {
    DeSurv::desurv_cv_seed_bank(dg$rep[i], dg$fold[i], nf, dg$n_seed_fits[i],
                                base = FRESH_SEED_BANK_BASE, stride = stride)[1L]
  }, integer(1))
  chk_a <- nrow(dg) > 0L && all(dg$cv_seed %in% FRESH_CV_SEEDS) &&
    setequal(dg$cv_seed, FRESH_CV_SEEDS) && all(dg$seed_first >= FRESH_SEED_BANK_BASE) &&
    all(dg$seed_first == expect_first)
  if (!chk_a) stop("HARD CHECK (a) FAILED: diagnostics do not show the fresh partitions / seed banks.")
  message(sprintf("  check (a) PASS: %d diagnostic rows; cv_seed in {%s}; seed_first %d..%d (all = bank start); status not ok: %d",
                  nrow(dg), paste(sort(unique(dg$cv_seed)), collapse = ","),
                  min(dg$seed_first), max(dg$seed_first), sum(dg$status != "ok")))
  # (b) The arguments used differ from the stored ones only in cv_seeds and seed_bank_base.
  drop_seed <- function(a) a[setdiff(names(a), c("cv_seeds", "seed_bank_base"))]
  chk_b <- identical(drop_seed(rescore$cv_args), drop_seed(fin$cv_args)) &&
    identical(rescore$cv_args$cv_seeds, FRESH_CV_SEEDS) &&
    identical(rescore$cv_args$seed_bank_base, FRESH_SEED_BANK_BASE)
  if (!chk_b) stop("HARD CHECK (b) FAILED: used cv_args differ from the stored ones beyond the seeds.")
  message("  check (b) PASS: cv_args identical to the stored ones apart from cv_seeds and seed_bank_base")
  # (c) length(FRESH_CV_SEEDS) x nfolds fold rows per finalist (15 in full mode), one per unit.
  n_units <- length(FRESH_CV_SEEDS) * nf
  per_fin <- table(factor(rescore$fold_cindex$finalist_id, levels = fin$configs$finalist_id))
  units_ok <- all(tapply(paste(rescore$fold_cindex$rep, rescore$fold_cindex$fold),
                         rescore$fold_cindex$finalist_id, anyDuplicated) == 0L)
  if (!all(per_fin == n_units) || !units_ok || (!CONFIG$quick && n_units != 15L)) {
    stop("HARD CHECK (c) FAILED: fold rows per finalist ", paste(per_fin, collapse = "/"),
         ", expected ", n_units, ".")
  }
  message(sprintf("  check (c) PASS: %s fold rows per finalist (expected %d)",
                  paste(per_fin, collapse = "/"), n_units))

  print(rescore$configs[, c("finalist_id", "eval_id", fin$param_names, "gp_mean", "mean_cindex",
                            "mean_cindex_fresh", "delta_fresh_vs_bo", "train_cindex_bo",
                            "train_cindex_fresh")])
  message("=== Step 3d (rank mode, k = ", BO_K, ") complete ===")
} else {
  # ════════════════════════════ select mode ══════════════════════════════
  parse_ranks <- function(s) {                      # as code/03c
    parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
    out <- unlist(lapply(parts, function(p) {
      ab <- as.integer(strsplit(p, ":", fixed = TRUE)[[1]])
      if (length(ab) == 2L) seq(ab[1], ab[2]) else ab
    }))
    if (!length(out) || anyNA(out) || anyDuplicated(out)) stop("Bad DESURV_BO_RANKS: '", s, "'.")
    sort(out)
  }
  RANKS <- parse_ranks(Sys.getenv("DESURV_BO_RANKS", "2:12"))
  message("  select mode; ranks: ", paste(RANKS, collapse = ", "))

  fins <- lapply(RANKS, function(k) load_precomputed(sprintf("desurv_bo_finalists_k%d_tcgacptac", k)))
  bos  <- lapply(RANKS, function(k) load_precomputed(sprintf("desurv_bo_results_k%d_tcgacptac", k)))
  rss  <- lapply(RANKS, function(k) load_precomputed(sprintf("desurv_bo_rescore_k%d_tcgacptac", k)))
  names(fins) <- names(bos) <- names(rss) <- RANKS

  # ── Consistency: as code/03c, plus every rescore is of this rank's finalists on the fresh CV ──
  drop_exec <- c("k_grid", "parallel_grid", "ncores_grid")
  strip <- function(a) { a <- a[setdiff(names(a), drop_exec)]; a[sort(names(a))] }
  for (i in seq_along(RANKS)) {
    f <- fins[[i]]; r <- rss[[i]]; k <- RANKS[i]
    stopifnot(inherits(f, "desurv_bo_finalists"), all(f$configs$k == k),
              all(f$fold_cindex$k == k), as.integer(f$cv_args$k_grid) == k,
              identical(f$param_names, fins[[1]]$param_names),
              identical(strip(f$cv_args), strip(fins[[1]]$cv_args)),
              r$k == k, identical(r$configs$eval_id, f$configs$eval_id),
              identical(r$configs$finalist_id, f$configs$finalist_id),
              identical(r$cv_args$cv_seeds, FRESH_CV_SEEDS),
              identical(r$cv_args$seed_bank_base, FRESH_SEED_BANK_BASE),
              all(r$diagnostics$cv_seed %in% FRESH_CV_SEEDS),
              all(r$diagnostics$seed_first >= FRESH_SEED_BANK_BASE),
              all(r$fold_cindex$k == k))
  }

  # ── Pool (as code/03c): finalist ids unique across ranks; desurv_select_k_paired_1se() keys
  #    finalists and their paired units by finalist_id. The same ids for A and B. ──
  pool <- function(cfg_list, fc_list) {
    configs <- do.call(rbind, cfg_list)
    configs$rank_finalist_id <- configs$finalist_id
    configs$finalist_id <- seq_len(nrow(configs))
    id_map <- paste(configs$k, configs$rank_finalist_id)
    fold_cindex <- do.call(rbind, fc_list)
    fold_cindex$rank_finalist_id <- fold_cindex$finalist_id
    fold_cindex$finalist_id <- match(paste(fold_cindex$k, fold_cindex$rank_finalist_id), id_map)
    stopifnot(!anyNA(fold_cindex$finalist_id))
    rownames(configs) <- rownames(fold_cindex) <- NULL
    configs$k_grid <- configs$k
    list(configs = configs, fold_cindex = fold_cindex)
  }
  B <- pool(lapply(rss, `[[`, "configs"), lapply(rss, `[[`, "fold_cindex"))
  A <- pool(lapply(fins, `[[`, "configs"), lapply(fins, `[[`, "fold_cindex"))
  stopifnot(identical(A$configs$eval_id, B$configs$eval_id), identical(A$configs$k, B$configs$k))
  # B's pooled object: the fresh C is the finalists' C (mean_cindex); the BO-partition value
  # is kept as mean_cindex_bo_partition.
  B$configs$mean_cindex_bo_partition <- B$configs$mean_cindex
  B$configs$mean_cindex <- B$configs$mean_cindex_fresh
  param_names <- c("k_grid", fins[[1]]$param_names)
  sel_A <- DeSurv::desurv_select_k_paired_1se(A$fold_cindex, A$configs)

  desurv_bo_finalists <- cache_or_compute("desurv_bo_finalists_fresh_tcgacptac", {
    cv_args <- rss[[1]]$cv_args
    cv_args$k_grid <- NULL
    structure(
      list(
        configs = B$configs,
        fold_cindex = B$fold_cindex,
        selection = DeSurv::desurv_select_k_paired_1se(B$fold_cindex, B$configs),
        param_names = param_names,
        cv_args = cv_args,                   # the FRESH CV; k_grid per rank in configs
        n_per_k = fins[[1]]$n_per_k,
        ranks = RANKS,
        selection_rule = SELECTION_RULE,
        source = "rank-stratified BO finalists rescored on fresh CV, rev 3 amendment 1b (code/03d)"
      ),
      class = "desurv_bo_finalists"
    )
  })
  print(desurv_bo_finalists)

  # ── Selection objects, as code/03c, tagged with the rule ──
  tar_k_selection <- cache_or_compute("tar_k_selection_fresh_tcgacptac", {
    x <- k_selection_from_finalists(desurv_bo_finalists)
    attr(x, "selection_rule") <- SELECTION_RULE
    x
  })
  tar_params_best <- cache_or_compute("tar_params_best_fresh_tcgacptac", {
    x <- select_bo_params_finalists(desurv_bo_finalists)
    attr(x, "selection_rule") <- SELECTION_RULE
    x
  })
  stopifnot(identical(tar_params_best$k, tar_k_selection$k_selected),
            identical(attr(tar_k_selection, "selection_rule"), SELECTION_RULE),
            identical(attr(tar_params_best, "selection_rule"), SELECTION_RULE))

  tar_data_filtered <- cache_or_compute("tar_data_filtered_fresh_tcgacptac", {
    preprocess_training_data(data = tar_data, ngene = bo_ngene(tar_params_best),
                             method_trans_train = "rank")
  })

  # ── Side-by-side: B (fresh, robustness / sensitivity), A (BO partitions), C (best GP mean, descriptive) ──
  sel_B <- desurv_bo_finalists$selection
  hp <- fins[[1]]$param_names
  best_cols <- function(sel, cfg, tag) {
    pk <- sel$per_k
    row <- cfg[match(pk$finalist_id, cfg$finalist_id), , drop = FALSE]
    out <- data.frame(k = pk$k, rank_finalist_id = row$rank_finalist_id, eval_id = pk$eval_id,
                      cv = pk$cv_mean, mean_diff = pk$mean_diff, se_diff = pk$se_diff,
                      n_pairs = pk$n_pairs, within_1se = pk$within_1se,
                      selected = pk$k == sel$k_selected, row[, hp, drop = FALSE],
                      train_cindex_bo = row$train_cindex_bo, row.names = NULL)
    names(out)[-1] <- paste0(tag, "_", names(out)[-1])
    out
  }
  side <- cache_or_compute("bo_rank_selection_side_by_side_tcgacptac", {
    sb <- merge(best_cols(sel_B, B$configs, "B"), best_cols(sel_A, A$configs, "A"), by = "k")
    C_gp <- vapply(as.character(sb$k), function(kk) {
      h <- bos[[kk]]$history; ok <- h$status == "ok" & is.finite(h$mean_cindex)
      max(h$gp_mean[ok])
    }, numeric(1))
    sb$C_best_gp_mean <- unname(C_gp)
    sb$C_top <- sb$k == sb$k[which.max(sb$C_best_gp_mean)]
    # B minus A: per rank's best finalist under each evaluation, and for B's best finalist
    # itself (same configuration, fresh minus BO-partition C).
    sb$B_minus_A_best <- sb$B_cv - sb$A_cv
    sb$B_minus_A_same_finalist <- sb$B_cv -
      B$configs$mean_cindex_bo_partition[match(paste(sb$k, sb$B_rank_finalist_id),
                                               paste(B$configs$k, B$configs$rank_finalist_id))]
    sb$selected_by <- vapply(seq_len(nrow(sb)), function(i) {
      s <- c(if (sb$B_selected[i]) "B", if (sb$A_selected[i]) "A")
      if (length(s)) paste(s, collapse = "+") else ""
    }, character(1))
    fin_tab <- data.frame(B$configs[, c("finalist_id", "k", "rank_finalist_id", "eval_id", "source",
                                        hp, "gp_mean", "train_cindex_bo", "train_cindex_fresh")],
                          A_cv = B$configs$mean_cindex_bo_partition,
                          B_cv = B$configs$mean_cindex_fresh,
                          B_minus_A = B$configs$mean_cindex_fresh - B$configs$mean_cindex_bo_partition,
                          row.names = NULL)
    attr(sb, "finalists") <- fin_tab
    attr(sb, "selection") <- data.frame(
      rule = c("B fresh_partition_paired_1se (robustness / sensitivity)", "A bo_partition_paired_1se",
               "C best GP mean (descriptive, no SE)"),
      k_selected = c(sel_B$k_selected, sel_A$k_selected, sb$k[sb$C_top]),
      k_best = c(sel_B$k_best, sel_A$k_best, NA), runner_up_k = c(sel_B$runner_up_k, sel_A$runner_up_k, NA))
    attr(sb, "note") <- paste(
      "B: fresh CV (cv_seeds", paste(FRESH_CV_SEEDS, collapse = ","), "; seed banks from",
      FRESH_SEED_BANK_BASE, "), paired 1-SE over all ranks' finalists; robustness / sensitivity",
      "analysis of rank selection (amendment 1c); never a Step 5 input.",
      "A: the same finalists on the BO's own partitions (their stored fold C), same rule; comparison.",
      "C: max GP posterior mean in the rank BO history; descriptive only.",
      "B_minus_A_best mixes selection optimism with fresh-CV sampling variation (amendment 1b).",
      "train_cindex_bo: training C at the BO evaluation; NA for the k = 3 seed rows.")
    sb
  })
  utils::write.csv(side, file.path(RESULTS_DIR, "bo_rank_selection_side_by_side_tcgacptac.csv"),
                   row.names = FALSE)
  utils::write.csv(attr(side, "finalists"),
                   file.path(RESULTS_DIR, "bo_rank_selection_side_by_side_finalists_tcgacptac.csv"),
                   row.names = FALSE)
  print(side)
  print(attr(side, "selection"))
  if (sel_A$k_selected != sel_B$k_selected) {
    message("  NOTE: A (BO partitions) selects k = ", sel_A$k_selected, ", B (fresh) selects k = ",
            sel_B$k_selected, "; both reported (robustness / sensitivity, amendment 1c).")
  }
  message("  Fresh-partition selection: k=", tar_params_best$k, " (runner-up k=",
          tar_k_selection$runner_up_k, "), alpha=", round(tar_params_best$alpha, 3))
  message("=== Step 3d (select mode) complete ===")
}
