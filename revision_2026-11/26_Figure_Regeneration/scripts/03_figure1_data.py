#!/usr/bin/env python3
"""Stage 26, step 3 — write the data behind Figure 1, without plotting.

The figures are drawn in R (see R/figure1.R) so that they match the theme of
your existing map script. This step does the computation only: it applies the
non-viral exclusion, partitions the pan-genome, builds the frequency spectrum
and runs the accessory-genome PCA, writing one CSV per panel.

Outputs
  02_results/fig1_panelA_matrix_long.csv   genome x orthogroup, long format
  02_results/fig1_panelA_species.csv       column order and species per genome
  02_results/fig1_panelC_spectrum.csv      gene frequency spectrum
  02_results/fig1_panelD_partition.csv     core/soft-core/shell/cloud counts
  02_results/fig1_panelE_pca.csv           PC1/PC2 per genome
  03_summary/03_figure1_data.log
"""
import os
import sys

def _locate_common():
    env = os.environ.get("COMMON_DIR")
    if env and os.path.exists(os.path.join(env, "_common.py")):
        return env
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(8):
        if os.path.exists(os.path.join(d, "_common.py")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    sys.exit(
        "could not find _common.py in any parent of\n  "
        + os.path.abspath(__file__)
        + "\nIt belongs in the revision_2026-11 directory alongside the "
          "numbered stages.\nIf it lives somewhere else, point at it:\n"
          "  export COMMON_DIR=/path/to/the/directory/holding/_common.py")


sys.path.insert(0, _locate_common())
from _common import stage_dir, load_pa, groups                    # noqa: E402

import numpy as np                                                # noqa: E402
import pandas as pd                                               # noqa: E402

D = stage_dir("26_Figure_Regeneration")
NONVIRAL = os.path.join(D, "01_inputs", "nonviral_orthogroups.txt")
LOG = open(os.path.join(D, "03_summary", "03_figure1_data.log"), "w")

# cloud orthogroups present in fewer genomes than this are dropped from the
# panel A heatmap only, to keep it readable. The submitted figure did the same.
PANEL_A_MIN_GENOMES = int(os.environ.get("PANEL_A_MIN_GENOMES", "2"))


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


if not os.path.exists(NONVIRAL):
    sys.exit(f"{NONVIRAL} not found; run 01_nonviral_screen.py first")

drop = {l.strip() for l in open(NONVIRAL) if l.strip()}
pa = load_pa()
say(f"matrix before exclusion: {pa.shape[0]} orthogroups x {pa.shape[1]} genomes")
pa = pa.loc[[i for i in pa.index if i not in drop]]
n_og, n_gen = pa.shape
say(f"matrix after excluding {len(drop)} non-viral: {n_og} x {n_gen}")

G = groups()
sp_of = {g: s for s, ids in G.items() for g in ids}
species = pd.Series({g: sp_of.get(g, "other") for g in pa.columns})
say("genomes per species: " +
    ", ".join(f"{s}={int((species == s).sum())}"
              for s in ["LSDV", "SPPV", "GTPV", "other"]))

# ---- partition --------------------------------------------------------------
freq = pa.sum(axis=1) / n_gen
part = pd.Series(index=pa.index, dtype=object)
part[freq >= 0.99] = "core"
part[(freq >= 0.95) & (freq < 0.99)] = "soft-core"
part[(freq >= 0.15) & (freq < 0.95)] = "shell"
part[freq < 0.15] = "cloud"
counts = part.value_counts().reindex(["core", "soft-core", "shell", "cloud"])
say("\npartition:")
for k, v in counts.items():
    say(f"  {k:<10} {int(v):>4}  ({100*v/n_og:.1f}%)")
say(f"  {'total':<10} {n_og:>4}")
say("manuscript states: core 132, soft-core 10, shell 14, cloud 22 = 178")
if list(counts.values) != [132, 10, 14, 22]:
    say("  *** MISMATCH with the manuscript. Report before drawing. ***")

pd.DataFrame({"compartment": counts.index, "orthogroups": counts.values,
              "percent": (100 * counts.values / n_og).round(1)}).to_csv(
    os.path.join(D, "02_results", "fig1_panelD_partition.csv"), index=False)

# ---- panel A ----------------------------------------------------------------
SP_ORDER = ["LSDV", "SPPV", "GTPV"]
col_order = sorted(
    pa.columns,
    key=lambda g: (SP_ORDER.index(species[g]) if species[g] in SP_ORDER
                   else len(SP_ORDER), g))
keep = pa.index[pa.sum(axis=1) >= PANEL_A_MIN_GENOMES]
sub = pa.loc[keep, col_order]
row_order = list(sub.sum(axis=1).sort_values(ascending=False).index)
sub = sub.loc[row_order]
say(f"\npanel A: {len(sub)} of {n_og} orthogroups "
    f"(present in >= {PANEL_A_MIN_GENOMES} genomes), "
    f"{n_og - len(sub)} ultra-rare excluded from the heatmap")

long = sub.reset_index().melt(id_vars="Gene", var_name="genome",
                              value_name="present")
long["og_rank"] = long.Gene.map({g: i for i, g in enumerate(row_order)})
long["genome_rank"] = long.genome.map({g: i for i, g in enumerate(col_order)})
long = long[long.present == 1]          # only present cells; absent is background
long.to_csv(os.path.join(D, "02_results", "fig1_panelA_matrix_long.csv"),
            index=False)
say(f"panel A long format: {len(long)} present cells written")

pd.DataFrame({"genome": col_order,
              "genome_rank": range(len(col_order)),
              "species": [species[g] for g in col_order]}).to_csv(
    os.path.join(D, "02_results", "fig1_panelA_species.csv"), index=False)

# ---- panel C ----------------------------------------------------------------
n_present = pa.sum(axis=1)
edges = np.arange(0, n_gen + 11, 10)
hist, _ = np.histogram(n_present, bins=edges)
pd.DataFrame({"bin_low": edges[:-1], "bin_high": edges[1:],
              "bin_mid": (edges[:-1] + edges[1:]) / 2,
              "orthogroups": hist}).to_csv(
    os.path.join(D, "02_results", "fig1_panelC_spectrum.csv"), index=False)

# ---- panel E ----------------------------------------------------------------
acc = pa.loc[part.isin(["soft-core", "shell", "cloud"])]
X = acc.values.T.astype(float)
X = X - X.mean(axis=0)
U, S, Vt = np.linalg.svd(X, full_matrices=False)
pc = U[:, :2] * S[:2]
var = (S ** 2) / (S ** 2).sum() * 100
pd.DataFrame({"genome": acc.columns, "PC1": pc[:, 0], "PC2": pc[:, 1],
              "species": [species[g] for g in acc.columns],
              "PC1_var": round(float(var[0]), 1),
              "PC2_var": round(float(var[1]), 1)}).to_csv(
    os.path.join(D, "02_results", "fig1_panelE_pca.csv"), index=False)
say(f"\npanel E: {acc.shape[0]} accessory orthogroups, "
    f"PC1={var[0]:.1f}% PC2={var[1]:.1f}%")

say("\nFigure 1 caption numbers:")
say(f"  panel A: {len(sub)} orthogroups shown, "
    f"{n_og - len(sub)} ultra-rare cloud excluded")
say(f"  panel D: core {int(counts['core'])}, soft-core "
    f"{int(counts['soft-core'])}, shell {int(counts['shell'])}, "
    f"cloud {int(counts['cloud'])}")
say(f"  panel E: PC1 = {var[0]:.1f}%, PC2 = {var[1]:.1f}%")
say("\nnow draw it:  Rscript R/figure1.R")
LOG.close()
