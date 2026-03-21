# DeSurv-paper Makefile
#
# Usage:
#   make from-precomputed   # Figures + paper from pre-computed results (~2 min)
#   make quick              # Smoke test full pipeline end-to-end (~10 min)
#   make all NCORES=8       # Full pipeline re-computation (hours, HPC recommended)
#   make paper              # Just compile the manuscript

NCORES ?= 1

.PHONY: all quick from-precomputed paper clean install

# ── Default: from pre-computed results ────────────────────────────────────
from-precomputed: paper

# ── Quick mode (smoke test) ───────────────────────────────────────────────
quick:
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/01_install.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/02_load_data.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/03_bayesian_optimization.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/04_fit_models.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/05_external_validation.R
	@echo "=== Quick mode complete. Pipeline runs end-to-end. ==="

# ── Full pipeline ─────────────────────────────────────────────────────────
all:
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/01_install.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/02_load_data.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/03_bayesian_optimization.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/04_fit_models.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/05_external_validation.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/07_simulations.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/08_figures.R
	Rscript code/09_render_paper.R
	@echo "=== Full pipeline complete ==="

# ── Paper only ────────────────────────────────────────────────────────────
paper:
	Rscript code/09_render_paper.R

# ── Install DeSurv ────────────────────────────────────────────────────────
install:
	Rscript code/01_install.R

clean:
	rm -rf results/precomputed/*.rds
	@echo "Cleaned pre-computed results. Static figures and cv_grid results preserved."
