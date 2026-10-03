#!/usr/bin/env Rscript
# code/03e_published_config.R — Step 5 parameters: the prespecified PUBLISHED model configuration
# (rev 3 amendment 1c, docs/plans/2026-10-03-rev3-amendment-1c-roles-and-reconstruction.md)
#
# Writes tar_params_best_tcgacptac for the Step 5 fit (code/04) and nothing else. Not a selection:
# the values are constants from the audited May analysis, not optimized here or conditionally at
# k = 3. Provenance: the May BO's overall best evaluation (k = 7, C 0.6545342172; cached
# results/desurv_bo_results_tcgacptac.rds, md5 9d6dab3b...) with k replaced by the k = 3 that the
# published 1-SE rule selected (public:code/03:111-122). Identical in value to the published
# results/tar_params_best_tcgacptac.rds (md5 8106019c...).
#
# Run only through slurm/run_rev3_step5.sh, which gives it an empty results directory, exports
# DESURV_REQUIRE_MODEL_ROLE=published_configuration_k3 and runs code/04 immediately after.
# Refuses to overwrite: an existing tar_params_best here would be a stale A/B/other object.

source("code/00_helpers.R")

MODEL_ROLE <- "published_configuration_k3"
path <- file.path(RESULTS_DIR, "tar_params_best_tcgacptac.rds")
if (file.exists(path)) {
  stop(path, " already exists; Step 5 needs a fresh results directory (no stale selection objects).",
       call. = FALSE)
}
if (file.exists(file.path(RESULTS_DIR, "tar_k_selection_tcgacptac.rds"))) {
  stop("tar_k_selection_tcgacptac.rds is present in ", RESULTS_DIR,
       "; Step 5 uses no selection object.", call. = FALSE)
}

tar_params_best <- list(
  k      = 3L,
  alpha  = 0.333935556177224480,
  lambda = 0.349448737941329357,
  nu     = 0.055651216804573779,
  ntop   = 270
)
attr(tar_params_best, "model_role") <- MODEL_ROLE
attr(tar_params_best, "source") <- "published_May_global_BO_optimum_with_k_replaced_by_3"
attr(tar_params_best, "note") <- paste(
  "Hyperparameters of the May BO's overall best evaluation (k = 7, CV C 0.6545) with k set to 3,",
  "as the published pipeline did. Not optimized at k = 3; never cross-validated at k = 3 in May.",
  "Rev 3 amendment 1c: prespecified published model configuration re-fit under the rev 3 estimator.")

saveRDS(tar_params_best, path)
chk <- readRDS(path)
stopifnot(identical(attr(chk, "model_role"), MODEL_ROLE), identical(chk$k, 3L),
          isTRUE(all.equal(unlist(chk[c("alpha", "lambda", "nu", "ntop")]),
                           c(alpha = 0.33393555617722448, lambda = 0.34944873794132936,
                             nu = 0.055651216804573779, ntop = 270), tolerance = 0)))
message("  wrote ", path, ": k=3, alpha=", signif(chk$alpha, 7), ", lambda=", signif(chk$lambda, 7),
        ", nu=", signif(chk$nu, 7), ", ntop=", chk$ntop, "; model_role=", MODEL_ROLE)
