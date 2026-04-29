# DeSurv-paper Makefile
#
# Usage:
#   make from-precomputed   # Figures + paper from pre-computed results (~2 min)
#   make quick              # Smoke test full pipeline end-to-end (~10 min)
#   make all NCORES=8       # Full pipeline re-computation (hours, HPC recommended)
#   make paper              # Just compile the manuscript
#
# ntop-specific runs (steps 3-5, 8-9 only — assumes step 2 data exists):
#   make ntop NTOP=150 NCORES=8                              # Fixed ntop
#   make ntop-bo NTOP_LOWER=50 NTOP_UPPER=250 NCORES=8      # BO-tuned ntop

NCORES ?= 1
NTOP ?=
NTOP_LOWER ?=
NTOP_UPPER ?=

.PHONY: all quick from-precomputed paper clean install ntop ntop-bo

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
	DESURV_NTOP_LOWER=50 DESURV_NTOP_UPPER=300 Rscript code/09_render_paper.R

# ── Install DeSurv ────────────────────────────────────────────────────────
install:
	Rscript code/01_install.R

# ── ntop-specific pipeline (steps 3-5, 8-9) ─────────────────────────────
ntop:
	@test -n "$(NTOP)" || (echo "Error: NTOP not set. Usage: make ntop NTOP=150" && exit 1)
	DESURV_NTOP=$(NTOP) DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/03_bayesian_optimization.R
	DESURV_NTOP=$(NTOP) DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/04_fit_models.R
	DESURV_NTOP=$(NTOP) DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/05_external_validation.R
	DESURV_NTOP=$(NTOP) DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/08_figures.R
	DESURV_NTOP=$(NTOP) Rscript code/09_render_paper.R
	@echo "=== ntop=$(NTOP) pipeline complete ==="

ntop-bo:
	@test -n "$(NTOP_LOWER)" || (echo "Error: NTOP_LOWER not set. Usage: make ntop-bo NTOP_LOWER=50 NTOP_UPPER=250" && exit 1)
	@test -n "$(NTOP_UPPER)" || (echo "Error: NTOP_UPPER not set." && exit 1)
	DESURV_NTOP_LOWER=$(NTOP_LOWER) DESURV_NTOP_UPPER=$(NTOP_UPPER) DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/03_bayesian_optimization.R
	DESURV_NTOP_LOWER=$(NTOP_LOWER) DESURV_NTOP_UPPER=$(NTOP_UPPER) DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/04_fit_models.R
	DESURV_NTOP_LOWER=$(NTOP_LOWER) DESURV_NTOP_UPPER=$(NTOP_UPPER) DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/05_external_validation.R
	DESURV_NTOP_LOWER=$(NTOP_LOWER) DESURV_NTOP_UPPER=$(NTOP_UPPER) DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/08_figures.R
	DESURV_NTOP_LOWER=$(NTOP_LOWER) DESURV_NTOP_UPPER=$(NTOP_UPPER) Rscript code/09_render_paper.R
	@echo "=== ntop BO [$(NTOP_LOWER),$(NTOP_UPPER)] pipeline complete ==="

clean:
	rm -rf results/precomputed/*.rds
	@echo "Cleaned pre-computed results. Static figures and cv_grid results preserved."
