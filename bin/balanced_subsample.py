#!/usr/bin/env python3
"""Balanced sub-sampling replicates for the phylogeographic root test.

Public archives over-represent recent Eurasian outbreaks. Equalising
continental representation to the least-sampled continent tests whether an
inferred ancestral root survives that imbalance or is an artefact of it.
"""
import argparse
import pandas as pd
from Bio import Phylo

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--column", default="Continent")
    ap.add_argument("--replicates", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20260817)
    a = ap.parse_args()

    meta = pd.read_csv(a.metadata, sep="\t")
    meta.columns = [c.strip() for c in meta.columns]
    tips = {t.name for t in Phylo.read(a.tree, "newick").get_terminals()}
    meta = meta[meta.id.isin(tips)]

    n_min = meta[a.column].value_counts().min()
    print(f"{a.column} counts:\n{meta[a.column].value_counts().to_string()}")
    print(f"normalising each stratum to n = {n_min}")

    for rep in range(a.replicates):
        keep = (meta.groupby(a.column, group_keys=False)
                    .apply(lambda g: g.sample(n_min, random_state=a.seed + rep)))
        tree = Phylo.read(a.tree, "newick")
        for t in list(tree.get_terminals()):
            if t.name not in set(keep.id):
                tree.prune(t)
        Phylo.write(tree, f"rep{rep}_tree.nwk", "newick")
        keep.to_csv(f"rep{rep}_metadata.tsv", sep="\t", index=False)
        print(f"  rep{rep}: {len(keep)} tips")

if __name__ == "__main__":
    main()
