# ---------------------------------------------------------------------------
# PACA-AU de-duplication for COMBINED (pooled) validation analyses.
#
# PACA-AU array and RNA-seq are two platforms profiling a partially overlapping
# set of the SAME ICGC patients (shared SA... accession; RNA-seq sample IDs
# carry a "_seq" suffix). For any patient profiled on both platforms we keep ONE
# record in combined/pooled analyses, preferentially the RNA-seq profile and
# excluding the corresponding array profile -- consistent with the preprocessing
# rule used in the previously published DeCAF analysis (Peng et al.).
#
# Platform-specific (per-dataset) analyses are UNAFFECTED and retain both
# records: the correction concerns combined patient-level inference only.
#
# Single source of truth: sourced by every pooled-analysis script and by the
# manuscript (via paper/load_precomputed.R). Do not re-implement the rule
# elsewhere.
# ---------------------------------------------------------------------------

# Stable ICGC patient accession (strip the RNA-seq "_seq" suffix).
paca_accession <- function(id) sub("_seq$", "", as.character(id))

# Logical keep-mask for COMBINED analyses: drop PACA_AU_array rows whose
# accession also appears in PACA_AU_seq; retain everything else.
paca_combined_keep <- function(dataset, id) {
  dataset <- as.character(dataset); id <- as.character(id)
  seq_acc <- paca_accession(id[dataset == "PACA_AU_seq"])
  drop <- dataset == "PACA_AU_array" & (paca_accession(id) %in% seq_acc)
  !drop
}

# Filter a patient-level data.frame to the combined (de-duplicated) set and
# assert one row per unique patient accession. `id_col` must be platform-tagged
# (array IDs bare, seq IDs "..._seq").
dedup_combined <- function(df, dataset_col = "dataset", id_col = "id") {
  keep <- paca_combined_keep(df[[dataset_col]], df[[id_col]])
  out  <- df[keep, , drop = FALSE]
  pid  <- paca_accession(out[[id_col]])
  if (anyDuplicated(pid)) {
    stop(sprintf("dedup_combined: %d duplicate patient accession(s) remain after de-duplication.",
                 sum(duplicated(pid))))
  }
  out
}

# Expected combined-cohort counts (build-time safeguard against pipeline drift).
DESURV_EXPECTED_COMBINED_PATIENTS <- 570L
DESURV_EXPECTED_COMBINED_EVENTS   <- 388L
