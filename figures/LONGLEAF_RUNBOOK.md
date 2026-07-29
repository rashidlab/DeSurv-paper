# Longleaf figure-regeneration run-book (Nature Cancer resubmission)

Regenerating the main figures needs the populated `results/precomputed/` store, which exists
**only on Longleaf** (Amber's pipeline runs), plus the DeSurv package and Amber's ggplot2. In the
`DeSurv-paper-clean` local checkout `results/precomputed/` is empty, so figures cannot be rebuilt
there. This run-book is the exact sequence to do it on the cluster.

Work through it top to bottom. The per-figure content changes (what each panel should show) live in
`figures/FIGURE_REGENERATION.md`; this file is the **environment + build + verify** wrapper.

---

## 0. Manuscript-figure ↔ file ↔ producer map (numbering is permuted — do not trust the filename)

| Manuscript | File (`FIGURE_DIR`) | Built by | Changes needed |
|---|---|---|---|
| Fig 1 (schema) | `model_schematic_final.pdf` | manual (PowerPoint) | label fixes only — see FIGURE_REGENERATION.md |
| Fig 2 (biology) | `fig3_tcgacptac.pdf` | `code/09a_figures.R` | **add panel E** (GATA6) |
| Fig 3 (validation) | `fig4_tcgacptac.pdf` | `code/09a_figures.R` | **3C: NMF KM → supervised heatmap** (code already in 09a) |
| Fig 4 (simulations) | `fig2_tcgacptac.pdf` | `code/09c_sim_figures.R` | **4D scenario gradient**, 4B relabel |

So a full main-figure rebuild = run **both** `code/09a_figures.R` (Figs 2, 3) **and**
`code/09c_sim_figures.R` (Fig 4). Fig 4 additionally needs the step-07 simulation outputs
(`sim_figs_by_scenario`) already in the store.

---

## 1. Get onto the cluster and into the right checkout

```bash
ssh <onyen>@longleaf.unc.edu
cd <path-to-Amber's DeSurv-paper checkout>        # the one whose results/precomputed/ is populated
git fetch origin
git checkout naim/nature-resubmission
git pull origin naim/nature-resubmission          # pulls splot_cutpoint recovery (71c319f) + FIGURE_REGENERATION fixes
```

Confirm the store is present (this is what the clean checkout lacks):

```bash
ls results/precomputed/tar_fit_desurv_tcgacptac.rds results/precomputed/data_val_filtered_tcgacptac.rds
# both must exist; if empty, you are in the clean checkout, not the pipeline checkout
```

## 2. Environment

```bash
module load r/4.4.0
Rscript -e 'cat("R", as.character(getRversion()), "| ggplot2", as.character(packageVersion("ggplot2")),
               "| DeSurv", as.character(packageVersion("DeSurv")), "\n")'
```

Ensure the figure dependencies are installed (safe to re-run; only installs what's missing):

```bash
Rscript code/01_install.R      # includes DiceKriging, superpc, survminer, cowplot, glmnet, etc.
```

`splot_cutpoint` is now in the repo (`R/figure_plot_helpers.R`, sourced by `code/09a`), so the
previously-fatal "could not find function splot_cutpoint" no longer occurs.

## 3. ggplot2-version guard (only if the check above shows ggplot2 ≥ 4.0)

`code/09a_figures.R:233` uses `geom_errorbarh(...)`, which errors under ggplot2 ≥ 4.0. If step 2
reports ggplot2 4.x, replace that call with a horizontal `geom_errorbar`:

```r
# code/09a_figures.R ~L233, forest panel — from:
#   geom_errorbarh(aes(xmin = lower, xmax = upper), ...)
# to:
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", ...)
```

(If ggplot2 is < 4.0 — the version Amber built with — skip this; `geom_errorbarh` is fine.)

## 4. Build the figures

```bash
Rscript code/09a_figures.R        # → fig3_tcgacptac.pdf (Fig 2), fig4_tcgacptac.pdf (Fig 3)
Rscript code/09c_sim_figures.R    # → fig2_tcgacptac.pdf (Fig 4)
```

Watch for the "Saved fig3_tcgacptac.pdf" / "Saved fig4_tcgacptac.pdf" messages. Fig 3 panel C should
now be the tuned-supervised-vs-DeSurv heatmap (the `km_block_4` block was already swapped in 09a);
sanity-check its values against `results/desurv_vs_supervised_tuned.rds$axis_decomposition`
(Penalized/Sparse Cox: 0.71 / 0.06 / 0.10; Supervised PCA: 0.46 / 0.14 / **0.31**).

## 5. Content edits that are NOT automatic (apply per FIGURE_REGENERATION.md)

These need code/caption edits beyond re-running the builders:

- **Fig 2E (GATA6 scatter).** `results/treated_cohort_stats.rds$gata6` stores only `rho/p/n`, not the
  raw points. Either save `D1c[ok]` / `g6[ok]` from `code/11_treated_cohort_analysis.R:165–167` into
  that cache and add panel E in `code/09a`, or build the panel directly from the COMPASS data.
  Then append the panel-E sentence to the `fig-pdac` caption in `paper/04_results_REVISED.Rmd`.
- **Fig 3C caption.** In the `fig-val` chunk (`paper/04_results_REVISED.Rmd`), change the "(B,C)
  Kaplan–Meier … DeSurv (B) and standard NMF (C)" wording to the panel-B-KM / panel-C-heatmap wording
  given in FIGURE_REGENERATION.md.
- **Fig 4D + 4B.** In `code/09c_sim_figures.R`, replace the k×α tuning-surface panel with the
  scenario gradient (prognostic → mixed → null; data already loaded as `sim_figs_by_scenario`, ids
  `R0_easy` / `R_mixed` / `R00_null`), move the tuning surface to the SI, and relabel 4B away from
  "Precision" (title "Recovery of true prognostic genes"). Update the `fig-sim` caption.
- **Fig 1 (manual).** In PowerPoint: `Z = XᵀW → Z = WᵀX`; confirm `α ∈ [0,1)`; optional lowercase
  panel letters. Re-export `model_schematic_final.pdf`.

## 6. Re-render and verify

```bash
make paper
```

Then confirm the manuscript is clean:

```bash
grep -c '??' paper/paper.pdf.log 2>/dev/null   # or check render stdout: 0 unresolved "??" refs
```

Manually confirm the three new panel references now resolve in `paper/paper.pdf`: **Fig 2E** (GATA6),
**Fig 3C** (supervised heatmap), **Fig 4D** (scenario gradient), and that there are 0 "??" and 0
undefined citations in the render log.

## 7. Commit the regenerated assets

```bash
git add figures/*.pdf paper/04_results_REVISED.Rmd code/09a_figures.R code/09c_sim_figures.R \
        code/11_treated_cohort_analysis.R results/treated_cohort_stats.rds
git commit -m "Regenerate main figures (Fig 2E GATA6, Fig 3C supervised heatmap, Fig 4D scenario gradient)"
git push origin naim/nature-resubmission
```

---

### Quick reference: minimal path if you only want Fig 3 (the 3C swap)

```bash
git pull origin naim/nature-resubmission
module load r/4.4.0
# apply §3 ggplot2 guard if ggplot2 >= 4.0
Rscript code/09a_figures.R       # rebuilds fig4_tcgacptac.pdf with the new 3C heatmap
# edit the fig-val caption (§5), then: make paper
```
