#!/usr/bin/env python3
"""Root-to-tip regression against collection date, stratified by continent.

A negative or absent slope is a result, not a failure: it indicates that
sampling date and geography are confounded, which precludes a molecular clock.
The per-stratum output makes that confounding visible rather than averaging it
away.
"""
import argparse
import numpy as np, pandas as pd
from Bio import Phylo
from scipy import stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--metadata", required=True,
                    help="TSV with columns id, Continent, Year")
    ap.add_argument("--species-filter", default=None,
                    help="restrict to one species, e.g. LSDV")
    ap.add_argument("--out-distances", required=True)
    ap.add_argument("--out-regression", required=True)
    a = ap.parse_args()

    tree = Phylo.read(a.tree, "newick")
    meta = pd.read_csv(a.metadata, sep="\t")
    meta.columns = [c.strip() for c in meta.columns]
    if a.species_filter and "Species" in meta.columns:
        meta = meta[meta.Species == a.species_filter]

    dist = {t.name: tree.distance(tree.root, t) for t in tree.get_terminals()}
    df = meta.assign(root_to_tip=meta.id.map(dist)).dropna(subset=["root_to_tip"])
    df["Year"] = pd.to_numeric(df.Year, errors="coerce")
    df = df.dropna(subset=["Year"])
    df.to_csv(a.out_distances, sep="\t", index=False)

    def fit(sub, name):
        if len(sub) < 5:
            return {"subset": name, "n": len(sub), "slope": None,
                    "r_squared": None, "p_value": None}
        r = stats.linregress(sub.Year, sub.root_to_tip)
        return {"subset": name, "n": len(sub), "slope": r.slope,
                "r_squared": r.rvalue ** 2, "p_value": r.pvalue}

    rows = [fit(df, "all")]
    rows.append(fit(df[df.Year >= 2000], "post_2000"))
    for c, sub in df.groupby("Continent"):
        rows.append(fit(sub, f"continent_{c}"))

    out = pd.DataFrame(rows)
    out.to_csv(a.out_regression, sep="\t", index=False)
    print(out.to_string(index=False))

    if len(df.Continent.unique()) > 1:
        groups = [g.root_to_tip.values for _, g in df.groupby("Continent")]
        F, p = stats.f_oneway(*groups)
        print(f"\nANOVA, root-to-tip by continent: F = {F:.2f}, p = {p:.4g}")

if __name__ == "__main__":
    main()
