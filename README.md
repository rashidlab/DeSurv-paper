# DeSurv Paper — Reproducible Analysis

Reproducible analysis for: **"Survival-guided matrix factorization identifies reproducible prognostic programs in pancreatic cancer"**

By Amber M. Young, Alisa Yurovsky, Xianlu Laura Peng, Didong Li, Jen Jen Yeh, and Naim U. Rashid.

## Quick Start

### Compile the manuscript from pre-computed results (~2 min)

```bash
git clone https://github.com/rashidlab/DeSurv-paper.git
cd DeSurv-paper
make paper                       # Compile paper + supplement
```

### Smoke test the full pipeline (~10 min, any laptop)

```bash
make quick
```

### Full re-computation (HPC recommended)

```bash
make all NCORES=8
# Or on Slurm:
sbatch slurm/run_full_pipeline.sh
```

## Repository Structure

```
DeSurv-paper/
├── code/                          # Analysis pipeline (numbered scripts)
│   ├── 00_helpers.R               #   Shared config, caching utilities
│   ├── 01_install.R               #   Install DeSurv package
│   ├── 02_load_data.R             #   Load TCGA+CPTAC training data
│   ├── 03_bayesian_optimization.R #   Hyperparameter search via BO
│   ├── 04_fit_models.R            #   Multi-start fitting + consensus init
│   ├── 05_external_validation.R   #   Project to 5 validation cohorts
│   ├── 06_sensitivity_analysis.R  #   K-sensitivity grid analysis
│   ├── 07_simulations.R           #   Simulation studies (3 scenarios)
│   ├── 08_figures.R               #   Figure generation
│   └── 09_render_paper.R          #   Compile manuscript
│
├── R/                             # Helper functions
│   ├── simulation_functions/      #   Simulation data generation
│   └── *.R                        #   Data loading, BO, validation, figures
│
├── paper/                         # Manuscript source
│   ├── paper.Rmd                  #   Main document
│   ├── si_appendix.Rmd            #   SI Appendix
│   ├── load_precomputed.R         #   Loads results from results/precomputed/
│   └── *.Rmd                      #   Child sections
│
├── results/
│   ├── precomputed/               # Pre-computed analysis objects (Zenodo)
│   └── cv_grid/                   # K-sensitivity analysis tables
│
├── figures/                       # Static PDF figures
├── data/                          # Input datasets (see data/README.md)
├── slurm/                         # HPC job scripts
│
├── Makefile                       # make quick / make all / make paper
└── run_pipeline.R                 # Alternative: Rscript run_pipeline.R --quick
```

## Reproducing the Analysis

There are three ways to use this repository, depending on your goal:

### 1. Compile the manuscript (no recomputation)

The `results/precomputed/` directory contains all intermediate analysis objects. The paper `.Rmd` files load these directly via `readRDS()`, so you can compile the manuscript without running any analysis:

```bash
make paper
```

### 2. Quick mode (verify pipeline runs end-to-end)

Quick mode runs every pipeline step with reduced iterations (~10 minutes on a laptop):

```bash
make quick
```

This uses `DESURV_QUICK=TRUE` to set: 4 BO iterations (vs 100), 3 initialization seeds (vs 100), 2 simulation replicates (vs 100). Results won't match production, but every step completes successfully.

### 3. Full reproduction (requires HPC or patience)

```bash
make all NCORES=8
```

Or step-by-step:
```bash
Rscript run_pipeline.R --full --ncores 8
```

Or individual steps:
```bash
Rscript run_pipeline.R --full --step 3 --only  # Just BO
```

## DeSurv Package

This analysis requires the [DeSurv](https://github.com/rashidlab/DeSurv) R package:

```r
devtools::install_github("rashidlab/DeSurv")
```

## Data Availability

Input data files and pre-computed results are included in this repository. See `data/README.md` for details on each dataset and its source (TCGA, CPTAC, GEO, ArrayExpress, ICGC).

## Development History

Development history prior to submission is archived at the original repositories:
- Paper: https://github.com/ayoung31/DeSurv-paper
- Package: https://github.com/ayoung31/DeSurv

## License

MIT License. See [LICENSE](LICENSE).
