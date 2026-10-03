# test_rev3_mechanics.R — MECHANICS TEST ONLY for the rev 3 pooling path, incl. the
# owner-approved follow-up fixes (FIX 1/2/3 + descriptive "May-comparable SE").
#
# Per the brief: do NOT obtain or run on real rank objects other than the k2 one already
# staged. No real 11-rank cross-rank selection is performed or reported.
this_dir <- "analysis/published-procedure-reconstruction"
source(file.path(this_dir, "reconstruct.R"))

in_dir <- "/tmp/claude-1000/-home-naimrashid-Downloads-DeSurv-paper-clean/4d6f63f8-6e3b-4abd-8ce0-f3388c26d037/scratchpad/recon_inputs"
k2_path <- file.path(in_dir, "desurv_bo_results_k2_tcgacptac.rds")
main_path <- file.path(in_dir, "step4_main_desurv_bo_results_tcgacptac.rds")
noise_path <- file.path(in_dir, "noise_results_3500640.rds")
stopifnot(file.exists(k2_path), file.exists(main_path), file.exists(noise_path))

k2 <- readRDS(k2_path)
main_bo <- readRDS(main_path)
noise_results <- readRDS(noise_path)

cat("md5 of step4 main BO object (must match PROVENANCE ad25d5342a77572f71c5dd5e15182a27):\n")
cat(tools::md5sum(main_path), "\n")

# ============================================================================
cat("\n=== FIX 1: restore inherited diagnostics from the step 4 main BO (real k2 object) ===\n")
before_diag_ids <- unique(k2$diagnostics$eval_id)
inh_before <- inherited_rows(k2$history, k2$diagnostics)
cat("inherited eval_ids before restoration:", paste(inh_before$eval_id, collapse = ","), "\n")
stopifnot(setequal(inh_before$eval_id, c(1, 2, 3, 4)))

fix1 <- restore_inherited_diagnostics_from_main(2L, k2$history, k2$diagnostics, main_bo)
cat("\nFIX 1 report (matched main eval, means, hard checks):\n")
print(fix1$report)
stopifnot(nrow(fix1$report) == 4L, all(fix1$report$all_checks_pass))
stopifnot(isTRUE(all.equal(fix1$report$history_mean, fix1$report$main_history_mean, tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(fix1$report$history_mean, fix1$report$recomputed_mean, tolerance = 1e-12)))

after_diag_ids <- unique(fix1$diagnostics$eval_id)
cat("\neval_ids with diagnostics: before =", length(before_diag_ids),
    " after FIX 1 =", length(after_diag_ids), "(expect +4)\n")
stopifnot(length(after_diag_ids) == length(before_diag_ids) + 4L)
stopifnot(all(c(1, 2, 3, 4) %in% after_diag_ids))

k2_restored <- k2
k2_restored$diagnostics <- fix1$diagnostics

cat("\nper-k table for k2 now that inherited evals 1-4 have finite c_se:\n")
pooled_k2_only <- pool_rank_bo_results(list("2" = k2_restored))
per_k_fix1 <- compute_per_k_table(pooled_k2_only)
print(per_k_fix1$per_k[, c("eval_id_orig", "mean_cindex", "c_se")])
stopifnot(all(is.finite(per_k_fix1$per_k$c_se[per_k_fix1$per_k$eval_id_orig %in% c(1, 2, 3, 4)])))

# ============================================================================
cat("\n=== FIX 2: seed-bank noise job (3500640) for the k=3 seed configs ===\n")
noise_tbl <- noise_results_to_diag_table(noise_results)
cat("configs x banks present in the noise file:\n")
print(unique(noise_tbl[, c("noise_label", "noise_bank")]))
b4_published <- noise_tbl[noise_tbl$noise_label == "B_published_k3" & noise_tbl$noise_bank == 4000000, ]
b4_old <- noise_tbl[noise_tbl$noise_label == "C_oldBObest_k3" & noise_tbl$noise_bank == 4000000, ]
cat("\nB_published_k3 @ bank 4e6: n_fold_rows=", nrow(b4_published),
    " recomputed_mean=", recompute_mean_from_fold_rows(b4_published),
    " reported_mean_cindex=", b4_published$reported_mean_cindex[1], "\n")
cat("C_oldBObest_k3 @ bank 4e6: n_fold_rows=", nrow(b4_old),
    " recomputed_mean=", recompute_mean_from_fold_rows(b4_old),
    " reported_mean_cindex=", b4_old$reported_mean_cindex[1], "\n")
stopifnot(isTRUE(all.equal(recompute_mean_from_fold_rows(b4_published), b4_published$reported_mean_cindex[1], tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(recompute_mean_from_fold_rows(b4_old), b4_old$reported_mean_cindex[1], tolerance = 1e-12)))
cat("\nSelf-consistency hard check (compute_bo_eval_se-style fold aggregation reproduces the\n")
cat("noise job's own res$summary$mean_cindex) PASSES for both configs. No real k=3 rank object\n")
cat("is available to apply the transplant end-to-end (no inherited row to match against); the\n")
cat("function is exercised below on k2 (no k=3 rows present there) to confirm it is a correct\n")
cat("no-op rather than fabricating anything.\n")

fix2_noop <- restore_noise_bank_diagnostics(2L, k2$history, k2$diagnostics, noise_results)
cat("\nFIX 2 applied to k2 (no k=3 inherited rows exist in k2): report nrow =",
    nrow(fix2_noop$report), "(expect 0, since none of k2's inherited rows are k=3)\n")
stopifnot(nrow(fix2_noop$report) == 0L)
stopifnot(identical(dim(fix2_noop$diagnostics), dim(k2$diagnostics) + c(0L, ncol(fix2_noop$diagnostics) - ncol(k2$diagnostics))))

# Mechanics check with a synthetic inherited row at k=3 matching B_published_k3 exactly, to
# show the matching + hard-check + transplant logic itself works (SYNTHETIC: this is not a
# real rank history row).
cat("\n*** SYNTHETIC CHECK (not real rev 3 data): a fabricated k=3 'inherited' row with\n")
cat("*** alpha/lambda/nu set to B_published_k3's exact values and mean_cindex set to that\n")
cat("*** config's own recomputed mean, to confirm restore_noise_bank_diagnostics() transplants\n")
cat("*** and passes its hard check when a true match exists.\n")
synthetic_h <- k2$history[1, , drop = FALSE]
synthetic_h$eval_id <- 999L
synthetic_h$alpha_grid <- b4_published$alpha[1]
synthetic_h$lambda_grid <- b4_published$lambda[1]
synthetic_h$nu_grid <- b4_published$nu[1]
synthetic_h$ntop <- b4_published$ntop[1]
synthetic_h$mean_cindex <- b4_published$reported_mean_cindex[1]
synthetic_history <- rbind(k2$history, synthetic_h)
fix2_synth <- restore_noise_bank_diagnostics(3L, synthetic_history, k2$diagnostics, noise_results)
print(fix2_synth$report)
stopifnot(nrow(fix2_synth$report) == 1L, isTRUE(fix2_synth$report$transplanted[1]))

# ============================================================================
cat("\n=== FIX 3: guard against the missing_se fallback ===\n")
cat("Synthetic case: a pool whose best row is an UNRESTORED inherited eval (no diagnostics).\n")
k2_unrestored_best <- k2
# Force the global best row to be an inherited (no-diagnostics) eval by inflating its mean.
k2_unrestored_best$history$mean_cindex[k2_unrestored_best$history$eval_id == 1] <- 0.99
pooled_bad <- pool_rank_bo_results(list("2" = k2_unrestored_best))
caught <- tryCatch({
  assert_overall_best_has_finite_se(pooled_bad)
  "NO_ERROR"
}, error = function(e) conditionMessage(e))
cat("assert_overall_best_has_finite_se() on the rigged pool:\n", caught, "\n")
stopifnot(!identical(caught, "NO_ERROR"), grepl("missing_se", caught))

caught2 <- tryCatch({
  run_reconstruction(pooled_bad, "rigged", require_finite_best_se = TRUE)
  "NO_ERROR"
}, error = function(e) conditionMessage(e))
cat("run_reconstruction(..., require_finite_best_se=TRUE) on the rigged pool:\n", caught2, "\n")
stopifnot(!identical(caught2, "NO_ERROR"))
cat("Control: the SAME rigged pool with require_finite_best_se=FALSE (May-path default) does\n")
cat("NOT stop (documents the exact failure mode FIX 3 closes):\n")
ctrl <- run_reconstruction(pooled_bad, "rigged-control", require_finite_best_se = FALSE)
cat("reason:", ctrl$selection$reason, " k_selected:", ctrl$selection$k_selected, "\n")
stopifnot(identical(ctrl$selection$reason, "missing_se"))

# ============================================================================
cat("\n=== Descriptive-only: May-comparable SE (mechanics on k2's own best row) ===\n")
pooled_restored_k2 <- pool_rank_bo_results(list("2" = k2_restored))
per_k_restored <- compute_per_k_table(pooled_restored_k2)
best_eval_id <- per_k_restored$per_k$eval_id[which.max(per_k_restored$per_k$mean_cindex)]
diag_all <- collect_bo_diagnostics(pooled_restored_k2, pooled_restored_k2$history)
diag_best <- diag_all[diag_all$eval_id == best_eval_id, ]
may_se <- compute_may_comparable_se(diag_best)
cat("best eval_id (k2 pool):", best_eval_id, "\n")
print(may_se)
cat("mean SE_5fold across the 3 partitions:", mean(may_se$se_5fold, na.rm = TRUE), "\n")
cat("range SE_5fold:", paste(range(may_se$se_5fold, na.rm = TRUE), collapse = " - "), "\n")
stopifnot(nrow(may_se) == 3L, all(is.finite(may_se$se_5fold)))
cat("(Label carried in the table itself:)", unique(may_se$label), "\n")

cat("\n=== Directory-of-11 acceptance check (argument validation only, not run on real data) ===\n")
res <- tryCatch(
  run_rev3_pooled(in_dir, "/tmp/should_not_be_used", ks = 2:12, main_bo_path = main_path, noise_results_path = noise_path),
  error = function(e) e
)
cat("run_rev3_pooled(in_dir, ks=2:12) on the real (incomplete) input dir correctly errors:\n")
cat(conditionMessage(res), "\n")
stopifnot(inherits(res, "error"), grepl("missing rank objects", conditionMessage(res)))

cat("\nAll mechanics checks passed (FIX 1, FIX 2 self-consistency + synthetic transplant, FIX 3, descriptive SE).\n")
