# Stage 26 — Heaps refit, figure regeneration, Table 5 lookup

Python computes, R draws. The Python scripts in `scripts/` apply the
non-viral exclusion, refit Heaps' law, run the PCA and parse the Scoary
output, writing one CSV per panel into `02_results/`. The R scripts in `R/`
read those CSVs and produce the figures, using the theme and palette taken
from your Supplementary Figure S1 map script.

Nothing here needs SLURM. The heaviest step is 100 permutations over a
178 × 290 matrix, which takes seconds on the login node.

## Which figures need redrawing

**Three.** Everything else in the paper is unaffected.

| Figure | Status | Why |
|---|---|---|
| **Figure 1** | **redraw** | panels A, B and D change with the 178-orthogroup matrix (cloud 50 → 22); C and E are redrawn too so the whole figure comes from one matrix |
| **Figure 4A** | **redraw** | three African associations withdrawn |
| **Figure 5B** | **redraw** | extended with the ASR, masked and partition rows |
| Figure 2 | unchanged | Proksee circular alignment, independent of the orthogroup matrix |
| Figure 3 A/B/C | unchanged | recombination and fastBAPS results are unrevised |
| Figure 4B | unchanged | heatmap of named orthogroups; none of the withdrawn three appears in it |
| Figure 5A, 5C | unchanged | counts and the cutoff sweep are unrevised |
| Supp. Fig. S1 | unchanged | the 336-genome map predates the filtering |
| Supp. Fig. S2 | unchanged | the exclusion counts (10 / 32 / 4) are unrevised |
| Supp. Fig. S3 | unchanged | PastML reconstruction is unrevised |
| Supp. Fig. S4 | unchanged | root-to-tip regression is unrevised |
| Supp. Fig. S5 | unchanged | the 60-gene selection landscape is unrevised |

The new root-inference result is carried by Table 1, so it needs no figure.
If you would like one — a small panel showing the three rootings with their
AU and ELW values, plus the ten balanced replicates — say so and I will add
`R/figureS6.R`.

**One thing to check in Supplementary Figure S1 while you are there.** Its
caption says "Gold represents SPPV, Red represents GTPV", but the map script
assigns `SPPV = "#871746"` (maroon) and `GTPV = "#F7B978"` (light orange).
The caption and the figure disagree about which species is which colour.
Either the caption or the palette needs correcting before submission.

## Run it

```bash
export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
cd "$PANG/revision_2026-11/26_Figure_Regeneration"

conda activate iqtree_env        # or any env with pandas and numpy
python3 -c "import scipy" || conda install -y -c conda-forge scipy

# --- Python: compute -------------------------------------------------------
python3 scripts/01_nonviral_screen.py    # MUST print 28
python3 scripts/02_refit_heaps.py        # item 1 — k and gamma
python3 scripts/03_figure1_data.py       # panel data for Figure 1
python3 scripts/04_scoary_table.py       # item 5 — the Table 5 numbers

# --- R: draw ---------------------------------------------------------------
Rscript R/figure1.R                      # item 2
Rscript R/figure4a.R                     # item 3
Rscript R/figure5b.R                     # item 4
```

R needs `ggplot2`, `dplyr`, `patchwork` and `ggrepel`. If any are missing:

```r
install.packages(c("ggplot2", "dplyr", "patchwork", "ggrepel"))
```

Dependencies: steps 2 and 3 need step 1; `figure1.R` needs steps 2 and 3;
`figure4a.R` needs step 4; `figure5b.R` needs nothing but
`01_inputs/forest_rows.csv`, which ships with this stage.

## The theme

`R/00_theme.R` holds everything shared, so adjusting it changes all three
figures at once:

- `theme_cpx()` — `theme_minimal`, serif throughout, `#2C3E50` text, bold
  legend titles, a thin dark panel border, no minor grid. Taken from your map
  script, with axes kept since these panels need them.
- `SPECIES_COLS` — your exact three hex values: `#002F57`, `#871746`,
  `#F7B978`.
- `COMPARTMENT_COLS` — a sequential ramp anchored on the LSDV navy, so
  Figure 1D reads as one family rather than four unrelated hues.
- `CATEGORY_COLS`, `FOREST_COLS` — the same palette extended for the volcano
  and the forest plot.
- `save_fig()` — writes PDF via `cairo_pdf` and PNG at 600 dpi, as your map
  script does, falling back to the standard device if cairo is unavailable.
- `legend_inside()` — picks the right legend API for your ggplot2 version;
  3.5 deprecated `legend.position = c(x, y)`.

## What each step prints that you should check

| Step | Check |
|---|---|
| 1 | non-viral count is **28**, pan-genome after exclusion is **178** |
| 2 | γ should come out near **0.0351** — tell me if it does not |
| 2 | the caption string printed at the end, ready to paste into Figure 1B |
| 3 | partition is **132 / 10 / 14 / 22**; a mismatch is flagged explicitly |
| 3 | PC1 and PC2 percentages, against the Figure 1 caption |
| 4 | significant pairs should be **91**, across **16** traits and **15** orthogroups |
| 4 | the OG0000151 and OG0000004 values, which fill the two `% CONFIRM` gaps |
| figure4a.R | three rows withdrawn, 5 convergent points remain |
| figure5b.R | 13 rows, 5 of them without confidence intervals |

## Step 1 is the one that can go wrong

Everything downstream depends on knowing which 28 orthogroups are excluded.
Step 1 writes that list once and the rest read it. It refuses to write a list
that is not 28 long, because the manuscript states 206 − 28 = 178.

If you already have the 28 IDs from the earlier screen, that is the fastest
and safest route:

```bash
printf 'OG0000xxx\nOG0000yyy\n' > 01_inputs/nonviral_orthogroups.txt
python3 scripts/01_nonviral_screen.py          # verifies and reports
```

Failing that it re-derives the list with DIAMOND. **If it prints a count
other than 28, send it to me rather than overriding** — a different number
means the screen has drifted and the manuscript's 178 would be wrong.

## Step 4 needs to see both Scoary runs

It globs `*.results.csv` under `$PANG/11_Scoary` and
`$PANG/revision_2026-11/22_Scoary_Vaccine`. `OG0000004` appears only in the
280-genome sensitivity run, so if that lives elsewhere:

```bash
SCOARY_DIRS="$PANG/11_Scoary:/path/to/the/280/run" \
  python3 scripts/04_scoary_table.py
```

## Figure 5B and the missing intervals

Five rows are drawn as open diamonds with no error bars, because stage 25
saved aggregate counts rather than the per-gene table a bootstrap needs. The
caption says so. If you want intervals there too, say so and I will patch
stage 25 step 4 to dump the per-gene counts — about ten lines, then one
re-run of that step.

Row values live in `01_inputs/forest_rows.csv`. Edit that file and re-run to
redraw; nothing is recomputed. Leave `lo` and `hi` empty for a row with no
interval.

## If you prefer your own plotting code

Every panel's numbers are written as CSV in `02_results/`, so you can ignore
the R scripts entirely and feed those into whatever you already have:

| File | Panel |
|---|---|
| `fig1_panelA_matrix_long.csv` | heatmap, long format, present cells only |
| `fig1_panelA_species.csv` | column order and species per genome |
| `accumulation_curves.csv` | pan and core curves, mean and sd per N |
| `fig1_panelC_spectrum.csv` | gene frequency spectrum |
| `fig1_panelD_partition.csv` | compartment counts and percentages |
| `fig1_panelE_pca.csv` | PC1/PC2 per genome with variance explained |
| `fig4_panelA_points.csv` | every plotted volcano point |
| `03_summary/heaps_fit.csv` | k, γ, R² for both fitting methods |
| `03_summary/table5_rows.csv` | the Table 5 values |

## Files

| Path | Content |
|---|---|
| `scripts/01_nonviral_screen.py` | fixes the 28-orthogroup exclusion list |
| `scripts/02_refit_heaps.py` | accumulation curves and the Heaps fit |
| `scripts/03_figure1_data.py` | panel data for Figure 1 |
| `scripts/04_scoary_table.py` | Scoary parse, Table 5 rows |
| `R/00_theme.R` | shared theme, palette, saving helper |
| `R/figure1.R` | Figure 1, all five panels |
| `R/figure4a.R` | Figure 4A volcano |
| `R/figure5b.R` | Figure 5B forest plot |
| `01_inputs/forest_rows.csv` | editable Figure 5B rows |
| `04_figures/` | Figure_1, Figure_4A, Figure_5B as PDF and PNG |

## Tested

All four Python scripts and all three R scripts were run end to end against a
synthetic 290-genome, 206-orthogroup dataset before being sent, and the three
figures were inspected. The numbers they produced there are meaningless —
only the mechanics were tested. Six bugs were found and fixed that way:
overlapping species labels on the Figure 1 strip, the Heaps equation sitting
on the core-genome line, an invalid plotmath expression that aborted
`figure4a.R`, the ∞ and → glyphs failing to render in the serif PDF font, a
stray `NA` row in the forest plot, and the forest rows coming out in the
wrong order because a discrete scale was ordering its levels unexpectedly
(now replaced with explicit numeric positions).
