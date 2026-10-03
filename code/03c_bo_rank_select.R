#!/usr/bin/env Rscript
# code/03c_bo_rank_select.R — cross-rank selection after the rank-stratified BO (rev 3 amendment 1)
#
# Pools the per-rank finalists of code/03b_bo_rank.R (each re-scored on the identical
# repeated CV), applies the paired 1-SE rule across ranks, and writes the same selection
# objects code/03's main stage writes after its finalists.
#
# Env:     DESURV_BO_RANKS  ranks to pool, "2:12" (default) or a comma list such as "2,3,4"
# Inputs:  desurv_bo_finalists_k{K}_tcgacptac, desurv_bo_results_k{K}_tcgacptac (each K),
#          tar_data_tcgacptac
# Outputs: desurv_bo_finalists_tcgacptac   pooled desurv_bo_finalists object (k_grid tuned
#                                          across ranks, finalist ids unique)
#          tar_k_selection_tcgacptac       k_selection_from_finalists()
#          tar_params_best_tcgacptac       select_bo_params_finalists()
#          tar_data_filtered_tcgacptac     preprocess_training_data(ngene = BO_FIXED$ngene)
#          bo_rank_summary_tcgacptac       list(per_k, finalists): evaluations, best CV C,
#                                          finalists' rescored CV C and training C

message("=== Step 3c: cross-rank selection (rev 3 amendment 1) ===")
source("code/00_helpers.R")
library(DeSurv)
source("R/bo_helpers.R")
source("R/reported_settings.R")
source("R/bo_setup.R")
source("R/preprocess_helpers.R")

parse_ranks <- function(s) {
  parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
  out <- unlist(lapply(parts, function(p) {
    ab <- as.integer(strsplit(p, ":", fixed = TRUE)[[1]])
    if (length(ab) == 2L) seq(ab[1], ab[2]) else ab
  }))
  if (!length(out) || anyNA(out) || anyDuplicated(out)) stop("Bad DESURV_BO_RANKS: '", s, "'.")
  sort(out)
}
RANKS <- parse_ranks(Sys.getenv("DESURV_BO_RANKS", "2:12"))
message("  ranks: ", paste(RANKS, collapse = ", "))

tar_data <- load_precomputed("tar_data_tcgacptac")
fins <- lapply(RANKS, function(k) load_precomputed(sprintf("desurv_bo_finalists_k%d_tcgacptac", k)))
bos  <- lapply(RANKS, function(k) load_precomputed(sprintf("desurv_bo_results_k%d_tcgacptac", k)))
names(fins) <- names(bos) <- RANKS

# ── Consistency: one rank per object, identical scoring arguments across ranks ──
drop_exec <- c("k_grid", "parallel_grid", "ncores_grid")
strip <- function(a) { a <- a[setdiff(names(a), drop_exec)]; a[sort(names(a))] }
for (i in seq_along(RANKS)) {
  f <- fins[[i]]; k <- RANKS[i]
  stopifnot(inherits(f, "desurv_bo_finalists"), all(f$configs$k == k),
            all(f$fold_cindex$k == k), as.integer(f$cv_args$k_grid) == k,
            identical(f$param_names, fins[[1]]$param_names),
            identical(strip(f$cv_args), strip(fins[[1]]$cv_args)))
}

# ── Pool: renumber finalist ids uniquely (rank order, then the rank's own order) ──
configs <- do.call(rbind, lapply(fins, function(f) f$configs))
configs$rank_finalist_id <- configs$finalist_id
configs$finalist_id <- seq_len(nrow(configs))
id_map <- paste(configs$k, configs$rank_finalist_id)
fold_cindex <- do.call(rbind, lapply(fins, function(f) f$fold_cindex))
fold_cindex$rank_finalist_id <- fold_cindex$finalist_id
fold_cindex$finalist_id <- match(paste(fold_cindex$k, fold_cindex$rank_finalist_id), id_map)
stopifnot(!anyNA(fold_cindex$finalist_id))
rownames(configs) <- rownames(fold_cindex) <- NULL
# eval_id stays the rank-local BO evaluation id (unique only within k); the 1-SE rule
# uses it only as a within-k tie-break.

# k is a tuned coordinate of the pooled object, as in the main BO's finalists, so that
# select_bo_params_finalists() returns k with the other hyperparameters.
configs$k_grid <- configs$k
param_names <- c("k_grid", fins[[1]]$param_names)

desurv_bo_finalists <- cache_or_compute("desurv_bo_finalists_tcgacptac", {
  cv_args <- fins[[1]]$cv_args
  cv_args$k_grid <- NULL
  structure(
    list(
      configs = configs,
      fold_cindex = fold_cindex,
      selection = DeSurv::desurv_select_k_paired_1se(fold_cindex, configs),
      param_names = param_names,
      cv_args = cv_args,                     # common to all ranks; k_grid per rank in configs
      n_per_k = fins[[1]]$n_per_k,
      ranks = RANKS,
      source = "rank-stratified BO, rev 3 amendment 1 (code/03b + code/03c)"
    ),
    class = "desurv_bo_finalists"
  )
})
print(desurv_bo_finalists)

# ── Selection objects, exactly as code/03's main stage after its finalists ──
# Tagged with the rule (amendment 1b): code/04 refuses these objects when
# DESURV_REQUIRE_SELECTION_RULE asks for the fresh-partition rule (code/03d).
tar_k_selection <- cache_or_compute("tar_k_selection_tcgacptac", {
  x <- k_selection_from_finalists(desurv_bo_finalists)
  attr(x, "selection_rule") <- "bo_partition_paired_1se"
  x
})

tar_params_best <- cache_or_compute("tar_params_best_tcgacptac", {
  x <- select_bo_params_finalists(desurv_bo_finalists)
  attr(x, "selection_rule") <- "bo_partition_paired_1se"
  x
})
stopifnot(identical(tar_params_best$k, tar_k_selection$k_selected))

tar_data_filtered <- cache_or_compute("tar_data_filtered_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = bo_ngene(tar_params_best),
                           method_trans_train = "rank")
})

# ── Summary table ──────────────────────────────────────────────────────────
bo_rank_summary <- cache_or_compute("bo_rank_summary_tcgacptac", {
  cv_mean <- vapply(configs$finalist_id, function(id) {
    x <- fold_cindex$cindex[fold_cindex$finalist_id == id]; mean(x[is.finite(x)])
  }, numeric(1))
  per_k <- do.call(rbind, lapply(seq_along(RANKS), function(i) {
    b <- bos[[i]]; rb <- attr(b, "rank_bo"); h <- b$history
    ok <- h$status == "ok" & is.finite(h$mean_cindex)
    cf <- configs[configs$k == RANKS[i], , drop = FALSE]
    cvk <- cv_mean[configs$k == RANKS[i]]
    data.frame(
      k = RANKS[i],
      n_evals_ok = sum(ok),
      n_inherited = rb$n_inherited, n_seeds = rb$n_seeds,
      n_new = sum(ok) - rb$n_existing,
      best_bo_cv = max(h$mean_cindex[ok]),
      best_gp_mean = max(h$gp_mean[ok]),
      finalist_cv_best = max(cvk),
      finalist_cv = paste(sprintf("%.4f", cvk), collapse = "/"),
      finalist_train = paste(sprintf("%.3f", cf$train_cindex_bo), collapse = "/"),
      row.names = NULL
    )
  }))
  sel <- desurv_bo_finalists$selection$per_k
  per_k <- merge(per_k, sel[, c("k", "mean_diff", "se_diff", "within_1se")], by = "k", all.x = TRUE)
  per_k$selected <- per_k$k == tar_k_selection$k_selected
  per_k$runner_up <- per_k$k %in% tar_k_selection$runner_up_k
  finalists <- data.frame(configs[, c("finalist_id", "k", "rank_finalist_id", "eval_id", "source",
                                      "main_eval_id", fins[[1]]$param_names, "gp_mean",
                                      "bo_mean_cindex", "mean_cindex", "train_cindex_bo")],
                          cv_mean_units = cv_mean, row.names = NULL)
  list(per_k = per_k, finalists = finalists,
       note = paste("train_cindex_bo: mean training C of the reported fit over partition x fold",
                    "at the BO evaluation of that configuration (same partitions and seed bank);",
                    "NA for the k = 3 seed rows, evaluated outside the BOs. Diagnostic only."))
})
print(bo_rank_summary$per_k)

message("  Rank-stratified BO: k=", tar_params_best$k, " (runner-up k=", tar_k_selection$runner_up_k,
        "), alpha=", round(tar_params_best$alpha, 3))
message("=== Step 3c complete ===")
