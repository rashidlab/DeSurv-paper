# DeSurv Paper — Reproducible Analysis

Reproducible analysis for: **"Survival-guided matrix factorization identifies reproducible prognostic programs in pancreatic cancer"**

By Amber M. Young, Alisa Yurovsky, Xianlu Laura Peng, Didong Li, Jen Jen Yeh, and Naim U. Rashid.

## Prerequisites

- **R >= 4.3** with a C++ compiler (for RcppArmadillo)
- **Required R packages:** Install the DeSurv package first, then remaining dependencies:

```r
# Install DeSurv
devtools::install_github("rashidlab/DeSurv")

# Other key packages (install if missing)
install.packages(c("survival", "NMF", "ggplot2", "cowplot", "dplyr",
                    "survminer", "rmarkdown", "optparse"))

# Bioconductor packages
BiocManager::install(c("preprocessCore", "ComplexHeatmap"))
```

Or run: `Rscript code/01_install.R`

## Quick Start

### Compile the manuscript from pre-computed results (~2 min)

```bash
git clone https://github.com/rashidlab/DeSurv-paper.git
cd DeSurv-paper
make paper
```

This compiles `paper/paper.pdf` and `paper/si_appendix.pdf` from pre-computed analysis objects — no pipeline execution needed.

### Smoke test the full pipeline (~10-15 min, any laptop)

```bash
make quick
```

Runs every analysis step end-to-end with reduced iterations (4 BO iterations vs 100, 3 initialization seeds vs 100, 2 simulation replicates vs 100). Proves the code works on your system. Results will differ from production due to reduced iterations.

### Full production re-computation

**Locally (single machine, multi-core):**
```bash
make all NCORES=8
```

**On a Slurm HPC cluster:**
```bash
sbatch slurm/run_full_pipeline.sh
```

**Step by step (to run or re-run individual steps):**
```bash
Rscript run_pipeline.R --full --ncores 8              # All steps
Rscript run_pipeline.R --full --step 3 --only          # Just step 3 (BO)
Rscript run_pipeline.R --full --step 7 --only --ncores 30  # Just simulations
```

Full re-computation replaces the pre-computed results with freshly computed ones. BO and simulations have stochastic elements, so exact numerical values may differ slightly but conclusions are unchanged.

## Pipeline Steps

| Step | Script | What it does | Runtime (full) | Runtime (quick) |
|------|--------|-------------|----------------|-----------------|
| 1 | `code/01_install.R` | Install DeSurv and check dependencies | 1 min | 1 min |
| 2 | `code/02_load_data.R` | Load TCGA+CPTAC training data | 1 min | 1 min |
| 3 | `code/03_bayesian_optimization.R` | BO hyperparameter search (DeSurv, alpha=0, elbow-k) | 4-8 hours | 3 min |
| 4 | `code/04_fit_models.R` | Multi-start fitting with consensus initialization | 2-4 hours | 2 min |
| 5 | `code/05_external_validation.R` | Project to 5 independent PDAC cohorts | 10 min | 2 min |
| 6 | `code/07_simulations.R` | 3 scenarios x 100 replicates x 4 methods | 6-12 hours | 10 min |
| 8 | `code/08_figures.R` | Generate manuscript figures | 5 min | uses cached |
| 9 | `code/09_render_paper.R` | Compile paper + SI appendix | 2 min | 2 min |

Steps 3, 4, and 7 benefit from multiple cores. DeSurv handles within-step parallelism internally via `parallel::mclapply` — set `NCORES` to control this.

## How It Works

Each script uses a `cache_or_compute()` pattern: if a result file exists in `results/precomputed/`, it loads it; if not (or if `DESURV_RECOMPUTE=TRUE`), it computes and saves the result. This means:

- **`make paper`** — loads all 30 pre-computed results, renders the manuscript
- **`make quick`** — forces recomputation with reduced parameters, overwrites cached results
- **`make all`** — forces full production recomputation

Environment variables that control behavior:
- `DESURV_QUICK=TRUE` — reduced iterations for smoke testing
- `DESURV_RECOMPUTE=TRUE` — recompute even if cached results exist
- `DESURV_NCORES=N` — number of cores for parallel steps (default: 1)

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
│   ├── 07_simulations.R           #   Simulation studies (3 scenarios)
│   ├── 08_figures.R               #   Figure generation
│   ├── 09_render_paper.R          #   Compile manuscript
│   └── sim_helpers.R              #   Simulation infrastructure functions
│
├── R/                             # Helper functions
│   ├── simulation_functions/      #   Simulation data generation
│   └── *.R                        #   Data loading, BO, validation, figures
│
├── paper/                         # Manuscript source
│   ├── paper.Rmd                  #   Main document
│   ├── si_appendix.Rmd            #   SI Appendix
│   ├── load_precomputed.R         #   Loads results from results/precomputed/
│   └── *.Rmd                      #   Child sections (intro, methods, results, discussion)
│
├── results/
│   ├── precomputed/               # Pre-computed analysis objects (30 RDS files)
│   └── cv_grid/                   # K-sensitivity analysis tables
│
├── data/original/                 # Input datasets (7 PDAC cohorts)
├── figures/                       # Static PDF figures
├── slurm/                         # HPC job scripts
│
├── Makefile                       # make quick / make all / make paper
└── run_pipeline.R                 # Rscript run_pipeline.R --quick / --full
```

## DeSurv Package

This analysis requires the [DeSurv](https://github.com/rashidlab/DeSurv) R package:

```r
devtools::install_github("rashidlab/DeSurv")
```

DeSurv provides survival-driven nonnegative matrix factorization with Bayesian optimization for hyperparameter selection. See the package repository for documentation and vignettes.

## Data Availability

All input data files and pre-computed results are included in this repository. See `data/README.md` for provenance of each dataset:

- **TCGA-PAAD** — The Cancer Genome Atlas
- **CPTAC** — Clinical Proteomic Tumor Analysis Consortium
- **Moffitt** — GEO GSE71729
- **Puleo** — ArrayExpress E-MTAB-6134
- **Dijk** — ArrayExpress E-MTAB-6830
- **PACA-AU** — ICGC EGA EGAS00001000154

## Development History

Development history prior to submission is archived at the original repositories:
- Paper: https://github.com/ayoung31/DeSurv-paper
- Package: https://github.com/ayoung31/DeSurv

## License

MIT License. See [LICENSE](LICENSE).
