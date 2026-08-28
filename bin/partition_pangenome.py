#!/usr/bin/env python3
"""Partition an OrthoFinder gene-count matrix into core / soft-core / shell / cloud.

Also writes the binary presence/absence matrix used by Scoary and the list of
strictly single-copy orthogroups used for the supermatrix and codon analyses.
"""
import argparse, sys
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-counts", required=True)
    ap.add_argument("--core", type=float, default=0.99)
    ap.add_argument("--soft-core", type=float, default=0.95)
    ap.add_argument("--shell", type=float, default=0.15)
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-single-copy", required=True)
    ap.add_argument("--out-partitions", required=True)
    a = ap.parse_args()

    df = pd.read_csv(a.gene_counts, sep="\t", index_col=0)
    df = df.drop(columns=[c for c in ["Total"] if c in df.columns])
    n = df.shape[1]
    if n == 0:
        sys.exit("ERROR: no genome columns found in the gene-count matrix.")

    binary = (df > 0).astype(int)
    binary.to_csv(a.out_matrix)

    freq = binary.sum(axis=1) / n

    def label(f):
        if f >= a.core:      return "core"
        if f >= a.soft_core: return "soft_core"
        if f >= a.shell:     return "shell"
        return "cloud"

    part = pd.DataFrame({
        "orthogroup": freq.index,
        "n_genomes": binary.sum(axis=1).values,
        "frequency": freq.round(4).values,
        "partition": [label(f) for f in freq],
        "max_copies": df.max(axis=1).values,
        "single_copy_everywhere": ((df == 1).all(axis=1)).values,
    })
    part.to_csv(a.out_partitions, sep="\t", index=False)

    sc = part.loc[part.single_copy_everywhere, "orthogroup"]
    sc.to_csv(a.out_single_copy, index=False, header=False)

    counts = part.partition.value_counts()
    print(f"genomes: {n}   orthogroups: {len(part)}")
    for k in ["core", "soft_core", "shell", "cloud"]:
        c = int(counts.get(k, 0))
        print(f"  {k:<10} {c:5d}  ({100*c/len(part):.1f}%)")
    print(f"  single-copy in every genome: {len(sc)}")

if __name__ == "__main__":
    main()
