#!/usr/bin/env python3
"""
Identify which significance filter reproduces the 92 gene-trait pairs across
16 traits and 18 orthogroups reported in the manuscript.

Reads only; writes nothing.

Run from: /cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
"""
import os, glob
import numpy as np
import pandas as pd

os.chdir("/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis")
RESULTS_DIR = "11_Scoary/scoary_results_tree"

files = sorted(glob.glob(os.path.join(RESULTS_DIR, "*.results.csv"))) \
        or sorted(glob.glob(os.path.join(RESULTS_DIR, "*.csv")))
print(f"{len(files)} per-trait result files\n")

frames = []
for f in files:
    d = pd.read_csv(f)
    d["trait"] = os.path.basename(f).split(".")[0]
    frames.append(d)
D = pd.concat(frames, ignore_index=True)
print(f"Total rows written by Scoary across all traits: {len(D)}")

for c in ["Naive_p", "Bonferroni_p", "Benjamini_H_p", "Empirical_p",
          "Best_pairwise_comp_p", "Worst_pairwise_comp_p",
          "Max_supporting_pairs", "Max_opposing_pairs", "Odds_ratio"]:
    if c in D.columns:
        D[c] = pd.to_numeric(D[c], errors="coerce")

TARGET = (92, 16, 18)

def report(name, mask):
    sub = D[mask]
    trio = (len(sub), sub.trait.nunique(), sub.Gene.nunique())
    flag = "   <<< MATCHES MANUSCRIPT" if trio == TARGET else ""
    print(f"{name:58s} pairs={trio[0]:4d}  traits={trio[1]:3d}  orthogroups={trio[2]:3d}{flag}")

print("\n--- single filters ---")
report("Empirical_p < 0.05", D.Empirical_p < 0.05)
report("Bonferroni_p < 0.05", D.Bonferroni_p < 0.05)
report("Benjamini_H_p < 0.05", D.Benjamini_H_p < 0.05)
report("Naive_p < 0.05", D.Naive_p < 0.05)
report("Worst_pairwise_comp_p < 0.05", D.Worst_pairwise_comp_p < 0.05)
report("Best_pairwise_comp_p < 0.05", D.Best_pairwise_comp_p < 0.05)

print("\n--- combinations ---")
report("Empirical < 0.05 AND Bonferroni < 0.05", (D.Empirical_p < 0.05) & (D.Bonferroni_p < 0.05))
report("Empirical < 0.05 AND Benjamini_H < 0.05", (D.Empirical_p < 0.05) & (D.Benjamini_H_p < 0.05))
report("Empirical < 0.05 AND Worst_pairwise < 0.05", (D.Empirical_p < 0.05) & (D.Worst_pairwise_comp_p < 0.05))
report("Empirical < 0.05 AND Best_pairwise < 0.05", (D.Empirical_p < 0.05) & (D.Best_pairwise_comp_p < 0.05))
report("Bonferroni < 0.05 AND Worst_pairwise < 0.05", (D.Bonferroni_p < 0.05) & (D.Worst_pairwise_comp_p < 0.05))
report("Empirical <= 0.05 AND Bonferroni <= 0.05", (D.Empirical_p <= 0.05) & (D.Bonferroni_p <= 0.05))
report("Empirical < 0.05 AND Max_supporting_pairs >= 1", (D.Empirical_p < 0.05) & (D.Max_supporting_pairs >= 1))
report("Empirical < 0.05 AND Max_supporting >= Max_opposing", (D.Empirical_p < 0.05) & (D.Max_supporting_pairs >= D.Max_opposing_pairs))

print("\n--- study-wide Bonferroni on the naive p ---")
report("Naive_p < 1.27e-5", D.Naive_p < 1.27e-5)
report("Naive_p < 1.28e-5", D.Naive_p < 1.28e-5)
report("Naive_p < 1.27e-5 AND Empirical < 0.05", (D.Naive_p < 1.27e-5) & (D.Empirical_p < 0.05))

print("\n--- Tier 1 criterion as described in the paper ---")
t1 = (D.Empirical_p < 0.05) & (D.Max_supporting_pairs >= 3) & \
     ((D.Max_supporting_pairs - D.Max_opposing_pairs) > 0)
report("Empirical < 0.05, pairs >= 3, net positive", t1)
print("\nTier 1 rows under that definition:")
cols = [c for c in ["trait", "Gene", "Odds_ratio", "Empirical_p", "Bonferroni_p",
                    "Max_supporting_pairs", "Max_opposing_pairs"] if c in D.columns]
print(D[t1][cols].sort_values(["trait", "Empirical_p"]).to_string(index=False))
