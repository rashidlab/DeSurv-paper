# Longleaf figure-regeneration run-book (Nature Cancer resubmission)

> **STALE (Aug 2026).** This file describes a figure redesign that has since been
> completed differently, so treat it as history, not instructions. Two things in
> particular no longer hold. (1) **Figure numbering shifted**: the treated-cohort
> figure is now Fig. 4 and the simulations figure is now Fig. 5, because the
> treated section was merged into the external-validation section. References to
> "Fig 4 (simulations)" below mean what is now Fig. 5. (2) **Fig. 3 panel C no
> longer exists**: the dichotomized-risk-group KM panel was removed (see
> docs/LESSONS.md L7) and the supervised-score heatmap is now panel b, so the
> "3C" instructions below are superseded. **The mapping tables below are also wrong**
> (`FIGURE_REGENERATION.md` names `code/09a_figures.R` for the simulations figure,
> which is built by `code/09c_sim_figures.R`, and neither table lists the treated
> figure). The verified current mapping is:
>
> | Figure | Asset | Script |
> |---|---|---|
> | 1 (schema) | `figures/model_schematic_final.pdf` | manual (PowerPoint) |
> | 2 (factor structure) | `figures/fig3_tcgacptac.pdf` | `code/09a_figures.R` |
> | 3 (external validation) | `figures/fig4_tcgacptac.pdf` | `code/09a_figures.R` |
> | 4 (treated cohorts) | `figures/fig_treated_context.pdf` | `code/util_fig_treated_context.R` |
> | 5 (simulations) | `figures/fig2_tcgacptac.pdf` | `code/09c_sim_figures.R` |
>
> Note the asset filenames are permuted against the figure numbers.

The cached model store lives in the **flat `results/` directory** (`RESULTS_DIR <- "results"`,
`code/00_helpers.R:28` — the `results/precomputed/` path in the stale CLAUDE.md doc is not used).
The `DeSurv-paper-clean` checkout on this machine **already contains that store** (59 `.rds`), so
**Fig 3 can be, and has been, rebuilt locally** via `code/util_fig3_standalone_preview.R` (no DeSurv
/ DiceKriging / ggplot2-downgrade needed — see that script). Longleaf is still the cleanest place to
regenerate the **full** set in one pass, because `code/09a_figures.R` also recomputes cutpoint KMs
through `DeSurv::`, Fig 4 needs the step-07 simulation store, and Amber's env has the exact ggplot2
Amber built with. This run-book is the **environment + build + verify** wrapper; the per-panel
content changes live in `figures/FIGURE_REGENERATION.md`.

If you only need Fig 3, skip the cluster: run `Rscript code/util_fig3_standalone_preview.R` locally
(output `figures/standalone_preview/fig3_full.pdf`) and drop it in as `fig4_tcgacptac.pdf`.

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
cd <path-to-Amber's DeSurv-paper checkout>        # the one whose flat results/ store is populated
git fetch origin
git checkout naim/nature-resubmission
git pull origin naim/nature-resubmission          # pulls splot_cutpoint recovery (71c319f) + FIGURE_REGENERATION fixes
```

Confirm the store is present (flat `results/`, not `results/precomputed/`):

```bash
ls results/tar_fit_desurv_tcgacptac.rds results/data_val_filtered_tcgacptac.rds results/sim_figs_by_scenario.rds
# the first two drive Figs 2-3; sim_figs_by_scenario drives Fig 4
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
make paper 2>&1 | tee /tmp/desurv_render.log     # make paper -> code/10_render_paper.R (3-pass)
```

`code/10_render_paper.R` deletes the LaTeX aux/`.log` files on completion, so there is no persistent
`paper/paper.pdf.log` to grep — check the render **stdout** instead:

```bash
grep -nE '\?\?|undefined (reference|citation)|LaTeX Warning: (Reference|Citation)' /tmp/desurv_render.log
# expect no matches: 0 unresolved "??" cross-refs, 0 undefined citations
```

Then open `paper/paper.pdf` and confirm the three new panel references resolve: **Fig 2E** (GATA6),
**Fig 3C** (supervised heatmap), **Fig 4D** (scenario gradient).

## 7. Commit the regenerated assets

```bash
git add figures/*.pdf paper/04_results_REVISED.Rmd code/09a_figures.R code/09c_sim_figures.R \
        code/11_treated_cohort_analysis.R results/treated_cohort_stats.rds
git commit -m "Regenerate main figures (Fig 2E GATA6, Fig 3C supervised heatmap, Fig 4D scenario gradient)"
git push origin naim/nature-resubmission
```

---

### Quick reference: Fig 3 only (the 3C swap) — no cluster needed

```bash
Rscript code/util_fig3_standalone_preview.R      # -> figures/standalone_preview/fig3_full.pdf
cp figures/standalone_preview/fig3_full.pdf figures/fig4_tcgacptac.pdf   # drop-in as manuscript Fig 3
# edit the fig-val caption (§5), then: make paper
```

This reproduces `code/09a`'s Fig-3 output (A forest + B KM + C heatmap) from the cached `results/`
store using default-library ggplot2/survminer/cowplot — the fastest path when only Fig 3 changed.
