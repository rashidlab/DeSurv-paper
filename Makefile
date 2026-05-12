# DeSurv-paper Makefile
#
# Pipeline steps (run in order):
#   01  install packages
#   02  load + preprocess data
#   03  Bayesian optimisation (ntop tuned in [50, 300])
#   04  fit final models
#   05  external validation
#   06  CV grid search (HPC only, ~24-48h — separate make target)
#   07  simulations (HPC: SLURM array via slurm/run_simulations.sh)
#   07b merge simulation array results (HPC only)
#   08  cutpoint analysis
#   09a main figures
#   09b SI figures
#   09c simulation figures (needs step 07/07b output)
#   10  render paper
#
# Usage:
#   make all NCORES=8         # Full pipeline, simulations run locally (~6-12h)
#   make main NCORES=8        # Steps 01-05, 08, 09a, 09b (no simulations or render)
#   make from-precomputed     # Figures + paper from cached results (~2 min)
#   make quick                # Smoke test full pipeline end-to-end (~10 min)
#   make paper                # Just compile the manuscript
#   make cv-grid NCORES=32    # CV grid search only (also runs as part of make all)
#
# On HPC, use slurm/submit_all.sh to run main pipeline, simulation array,
# and CV grid search in parallel, then merge results and render.
#
# Default: DESURV_RECOMPUTE=TRUE (recompute everything).
# Use DESURV_RECOMPUTE=FALSE to load from cache.

NCORES ?= 1
# Use a portable Rscript invocation by default; override with `make RSCRIPT=...`
RSCRIPT ?= Rscript

.PHONY: all main quick from-precomputed paper clean install cv-grid

# ── Default: from pre-computed results ────────────────────────────────────
from-precomputed: paper

# ── Quick mode (smoke test) ───────────────────────────────────────────────
quick:
	DESURV_QUICK=TRUE DESURV_NCORES=1 $(RSCRIPT) code/01_install.R
	DESURV_QUICK=TRUE DESURV_NCORES=1 $(RSCRIPT) code/02_load_data.R
	DESURV_QUICK=TRUE DESURV_NCORES=1 $(RSCRIPT) code/03_bayesian_optimization.R
	DESURV_QUICK=TRUE DESURV_NCORES=1 $(RSCRIPT) code/04_fit_models.R
	DESURV_QUICK=TRUE DESURV_NCORES=1 $(RSCRIPT) code/05_external_validation.R
	@echo "=== Quick mode complete. Pipeline runs end-to-end. ==="

# ── Main pipeline (no simulations — for HPC parallel workflow) ────────────
main:
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/01_install.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/02_load_data.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/03_bayesian_optimization.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/04_fit_models.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/05_external_validation.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/08_cutpoint_analysis.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/09a_figures.R
	@echo "=== Main pipeline complete (figures 3-4 ready) ==="

# ── Full pipeline (simulations run locally) ───────────────────────────────
all:
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/01_install.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/02_load_data.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/03_bayesian_optimization.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/04_fit_models.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/05_external_validation.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/06_cv_grid.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/07_simulations.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/08_cutpoint_analysis.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/09a_figures.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/09b_si_figures.R
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/09c_sim_figures.R
	DESURV_RECOMPUTE=FALSE $(RSCRIPT) code/10_render_paper.R
	@echo "=== Full pipeline complete ==="

# ── Paper only ────────────────────────────────────────────────────────────
paper:
	DESURV_RECOMPUTE=FALSE $(RSCRIPT) code/10_render_paper.R

# ── Install DeSurv ────────────────────────────────────────────────────────
install:
	$(RSCRIPT) code/01_install.R

# ── CV grid search (k × alpha × ntop) — HPC only ─────────────────────────
cv-grid:
	DESURV_NCORES=$(NCORES) $(RSCRIPT) code/06_cv_grid.R

clean:
	rm -rf results/*.rds
	@echo "Cleaned pre-computed results. Static figures and cv_grid results preserved."
