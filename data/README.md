# Data

This directory should contain the input datasets for the DeSurv analysis.

## Required files

The following datasets are needed in `data/original/`:

### Training cohorts
- `TCGA_PAAD.rds`, `TCGA_PAAD.survival_data.rds`, `TCGA_PAAD_subtype.csv`
- `CPTAC.rds`, `CPTAC.survival_data.rds`, `CPTAC_subtype.csv`

### Validation cohorts
- `Dijk.rds`, `Dijk.survival_data.rds`, `Dijk_subtype.csv`
- `Moffitt_GEO_array.rds`, `Moffitt_GEO_array.survival_data.rds`, `Moffitt_GEO_array_subtype.csv`
- `PACA_AU_array.rds`, `PACA_AU_array.survival_data.rds`, `PACA_AU_array_subtype.csv`
- `PACA_AU_seq.rds`, `PACA_AU_seq.survival_data.rds`, `PACA_AU_seq_subtype.csv`
- `Puleo_array.rds`, `Puleo_array.survival_data.rds`, `Puleo_array_subtype.csv`

### Reference data
- `cmbSubtypes.RData` (combined molecular subtypes for sensitivity analysis)

## Data sources

- **TCGA_PAAD:** The Cancer Genome Atlas (https://www.cancer.gov/tcga)
- **CPTAC:** Clinical Proteomic Tumor Analysis Consortium
- **Moffitt:** GEO accession GSE71729
- **Puleo:** ArrayExpress E-MTAB-6134
- **Dijk:** ArrayExpress E-MTAB-6830
- **PACA-AU:** ICGC data portal, EGA study EGAS00001000154

## Download

Data files are included in this repository.
