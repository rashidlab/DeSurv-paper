# Main-figure regeneration instructions (Nature Cancer resubmission)

> **STALE (Aug 2026).** This file describes a figure redesign that has since been
> completed differently, so treat it as history, not instructions. Two things in
> particular no longer hold. (1) **Figure numbering shifted**: the treated-cohort
> figure is now Fig. 4 and the simulations figure is now Fig. 5, because the
> treated section was merged into the external-validation section. References to
> "Fig 4 (simulations)" below mean what is now Fig. 5. (2) **Fig. 3 panel C no
> longer exists**: the dichotomized-risk-group KM panel was removed (see
> docs/LESSONS.md L7) and the supervised-score heatmap is now panel b, so the
> "3C" instructions below are superseded. Current mapping of figure number to
> source file is in the table below, which is still correct.

The Results text was restructured (branch `naim/nature-resubmission`, commits `d2911ba`, `dc44361`)
into: **learn → discover → validate → distinguish → explain**. The prose now references figure
panels that the current figure PDFs do not yet contain. This file specifies every figure change so
the display items match the text.

The main figures are pre-generated PDFs included via `include_graphics`; the manuscript render
(`make paper`) does **not** rebuild them. They are produced by `code/09a_figures.R`:

| Manuscript | File | Producer |
|---|---|---|
| Fig 1 (schema) | `figures/model_schematic_final.pdf` | manual (PowerPoint) |
| Fig 2 (pdac)   | `figures/fig3_tcgacptac.pdf` | `code/09a_figures.R` |
| Fig 3 (val)    | `figures/fig4_tcgacptac.pdf` | `code/09a_figures.R` |
| Fig 4 (sim)    | `figures/fig2_tcgacptac.pdf` | `code/09a_figures.R` |

## Environment blockers (must fix before `code/09a_figures.R` will run)

Running `code/09a_figures.R` in the `DeSurv-paper-clean` checkout currently fails, independent of any
figure change:

1. **`splot_cutpoint` is undefined.** It is *called* at `code/09a_figures.R:325,331` but defined
   nowhere in the repo, in any gitignored file, or in the installed `DeSurv` package. The figures were
   last built in an environment where this helper existed — it needs to be restored/located and
   sourced (or the two calls replaced with the current cutpoint-plot helper).
2. **`DiceKriging`** was not installed (it is in `code/01_install.R` but was missing); `install.packages("DiceKriging")` fixes it.
3. **ggplot2 4.0.3 deprecations.** `geom_errorbarh()` now errors; replace with
   `geom_errorbar(..., orientation = "y")`. Check other 4.0 breakages while building.

---

## Fig 1 — DeSurv schema (MANUAL, PowerPoint)

`figures/model_schematic_final.pdf` (source `.pptx`). Two label fixes:
- Panel A score equation: **`Z = XᵀW` → `Z = WᵀX`** (the caption already says `Z = W^T X`; the image is wrong).
- Confirm the supervision-parameter label reads **α ∈ [0,1)**.
- (Optional, Nature style) lowercase panel letters.

No prose change needed.

---

## Fig 2 — biological factor structure: ADD panel E (GATA6)

`figures/fig3_tcgacptac.pdf`. Keep A–D; **add panel E**: scatter of the **D1 score vs GATA6 RNA-ISH
in COMPASS**, Spearman **ρ = 0.68, P < 0.001, n = 33**.

- **Data:** the raw per-sample points are computed in `code/11_treated_cohort_analysis.R:165–167`
  (`D1c` = scaled D1 score, `g6` = GATA6 ISH). The cache `results/treated_cohort_stats.rds$gata6`
  currently stores only `rho/p/n/group_means`, **not** the raw points — so either (a) save `D1c[ok]`
  and `g6[ok]` into that cache from `code/11`, or (b) rebuild the panel directly in `code/11`/`code/09a`
  from the COMPASS data.
- **Caption** (`fig-pdac` chunk, `paper/04_results_REVISED.Rmd`): append
  `(E) DeSurv D1 score versus GATA6 RNA in situ hybridization in COMPASS (Spearman correlation; n = 33).`
- **Prose already references it:** Section 2 ¶1 cites `Fig. \ref{fig:pdac}E`.
- **Provenance:** COMPASS is restricted-access; the Data Availability statement already lists it.

---

## Fig 3 — external validation: REDESIGN (biggest change)

`figures/fig4_tcgacptac.pdf`, built in `code/09a_figures.R` (~L454–518). Current layout: A = two forests
(DeSurv D1–D3 / NMF N1–N3), B = DeSurv KM, C = NMF KM. **New layout:**

- **3A** — keep the forest panel (DeSurv D1–D3 + matched-rank standard NMF). This is the continuous
  frozen-linear-predictor transfer result and is the visual center.
- **3B** — keep the DeSurv pooled KM.
- **3C** — **REPLACE the NMF KM** with a heatmap of the tuned supervised comparators' correspondence
  with the three DeSurv programs (this supports Section 4C). Values (|Spearman r|, from
  `results/desurv_vs_supervised_tuned.rds$axis_decomposition`):

  |               | D1   | D2 (proCAF) | D3   |
  |---------------|------|------|------|
  | Penalized Cox | 0.71 | **0.06** | 0.10 |
  | Supervised PCA| 0.46 | **0.14** | 0.31 |

  Message: supervised risk scores concentrate on the D1 direction and are near-orthogonal to the
  separately-resolved D2 stromal program. **Do NOT add a DeSurv-LP row** — DeSurv's trained score has
  β(D2)=0, so it would also read ≈0 on D2 and blur the point; the DeSurv contrast is that D1/D2/D3 exist
  as the *columns*.

  **Drop-in code** (already written into `code/09a_figures.R`, replacing the `km_block_4` block — verify
  once the script runs):
  ```r
  .axd <- readRDS("results/desurv_vs_supervised_tuned.rds")$axis_decomposition
  .hm_df <- data.frame(
    method  = factor(rep(c("Supervised PCA", "Penalized Cox"), each = 3),
                     levels = c("Supervised PCA", "Penalized Cox")),
    program = factor(rep(c("D1", "D2", "D3"), 2), levels = c("D1", "D2", "D3")),
    r = c(abs(as.numeric(.axd["Supervised PCA", c("D1","D2","D3")])),
          abs(as.numeric(.axd["Sparse Cox",     c("D1","D2","D3")]))))
  fig_supcorr_hm <- ggplot(.hm_df, aes(program, method, fill = r)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", r)), size = km_text_size / ggplot2::.pt) +
    scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1),
                        name = expression("|" * italic(r) * "|")) +
    labs(x = NULL, y = NULL, title = "Supervised score vs DeSurv program") + theme_pnas +
    theme(plot.title = element_text(size = 9), legend.position = "right")
  km_block_4 <- plot_grid(stack_surv(fig_median_survival_desurv, "DeSurv"),
                          ggdraw(km_legend_grob), fig_supcorr_hm,
                          nrow = 3, labels = c("B", "", "C"), rel_heights = c(5, 0.6, 4.2))
  ```
- **Caption** (`fig-val` chunk): change the "(B,C) Kaplan–Meier … DeSurv (B) and standard NMF (C)"
  wording to: `(B) Pooled DeSurv Kaplan–Meier curve. (C) Correspondence (|Spearman r|) of tuned
  supervised risk scores (penalized Cox, supervised PCA) with the three DeSurv programs (D1–D3).`
- **Prose already references it:** Section 4C cites `Fig. \ref{fig:val}C` as the supervised comparison.

---

## Fig 4 — simulations only

`figures/fig2_tcgacptac.pdf`. Make it **entirely simulation-based**:

- **4A** — test-set C-index, DeSurv vs NMF (keep). Numbers: median 0.837 vs 0.724, Δ 0.113, P<0.001.
- **4B** — recovery of true prognostic genes (keep the metric). **Relabel** the panel/axis away from
  "Precision": title `Recovery of true prognostic genes`, y-axis `Proportion of selected genes from the
  true prognostic program`. (Metric = fraction of selected genes that are true-program genes; only
  factor 1 is survival-driving, β = c(2,0,0).)
- **4C** — selected-rank distribution (keep).
- **4D** — **REPLACE** the current PDAC (k×α) tuning-surface panel (which is *not* a simulation) with the
  **scenario gradient**: performance across **prognostic → mixed → null**, visually delivering "benefit
  grows as variance and prognosis diverge." Data in `results/sim_figs_by_scenario.rds` (scenarios
  `R0_easy`, `R_mixed`, `R00_null`). **Move the tuning surface to the SI** (model-selection sensitivity).
- **Caption** (`fig-sim` chunk): update panel D from the GP tuning surface to the scenario-gradient
  description; note the tuning surface now lives in the SI.

---

## After regenerating: re-render and check

```
Rscript code/09a_figures.R      # rebuilds fig2/3/4_tcgacptac.pdf
make paper                      # re-includes them; verify 0 "??" and 0 undefined citations
```
Confirm the panel references now resolve: Fig 2E (GATA6), Fig 3C (supervised heatmap), Fig 4D (scenario
gradient). All numeric values above are verified against the current caches.
