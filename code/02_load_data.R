#!/usr/bin/env Rscript
# code/02_load_data.R — Load and validate all datasets
#
# Inputs:  data/original/*.rds, *.survival_data.rds, *_subtype.csv
# Outputs: results/precomputed/tar_data_tcgacptac.rds
#
# Training: TCGA_PAAD + CPTAC (pooled)
# Validation: Dijk, Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq, Puleo_array

message("=== Step 2: Load Data ===")
source("code/00_helpers.R")
source("R/load_data.R")
source("R/load_data_internal.R")

# ── Load pooled training data ─────────────────────────────────────────────
tar_data_tcgacptac <- cache_or_compute("tar_data_tcgacptac", {
  load_data(c("TCGA_PAAD", "CPTAC"))
})

message(sprintf("  Training data: %d genes x %d samples (%d events)",
                nrow(tar_data_tcgacptac$ex),
                ncol(tar_data_tcgacptac$ex),
                sum(tar_data_tcgacptac$sampInfo$event[tar_data_tcgacptac$samp_keeps])))

message("=== Step 2 complete ===")
