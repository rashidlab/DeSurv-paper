# ---------------------------------------------------------------------------
# Fixed, biologically prespecified reference-program panel displayed in
# main-text Fig. 2A/B (identical rows and order in BOTH panels), chosen
# independently of the observed rank-biserial values to avoid a second,
# outcome-dependent selection step. Ordered by biological family:
#   classical -> basal-like -> activated/proCAF/ECM -> restCAF/normal
#   -> immune -> exocrine/endocrine controls.
#
# Single source of truth, sourced by code/09a_figures.R (production pipeline),
# code/util_fig2_standalone.R (standalone preview), and
# code/util_rb_enrichment_table.R (SI source table), so `make all` reproduces
# the exact manuscript figure rather than the observed-effect fallback.
# ---------------------------------------------------------------------------
fig2_display_sigs <- c(
  # classical (tumor)
  "Moffitt_Classical", "PurIST_Classical", "Collisson_Classical",
  # basal-like (tumor)
  "Moffitt_Basal-like", "PurIST_BasalLike", "Bailey_Squamous",
  # activated / proCAF / ECM / desmoplastic (stroma)
  "Moffitt_Activated", "SCISSORS_proCAF", "deCAF_proCAF", "Maurer_ECM",
  # restCAF / normal stroma
  "Moffitt_Normal", "SCISSORS_restCAF", "deCAF_restCAF",
  # immune
  "DECODER_Immune", "Moffitt_Immune",
  # exocrine / endocrine controls
  "DECODER_Exocrine", "DECODER_Endocrine")
