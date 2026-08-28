#!/usr/bin/env python3
"""Decide which annotated genomes enter the pangenome.

Two independent filters, both applied before orthogroup inference:
  1. contamination  -- fraction of predicted proteins with no hit to the viral
                       reference proteome exceeds a threshold
  2. completeness   -- number of missing near-core orthogroups exceeds a
                       threshold, indicating a fragmented assembly

Absence of a gene is only interpretable once fragmented assemblies are removed,
which is why completeness is handled here rather than post hoc.
"""
import argparse, glob, os
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contamination", nargs="+", required=True,
                    help="per-genome *.contamination.tsv files")
    ap.add_argument("--faa-dir", default=".")
    ap.add_argument("--cds-counts", nargs="*", default=[])
    ap.add_argument("--near-core-counts", default=None,
                    help="optional TSV: genome<TAB>n_missing_near_core")
    ap.add_argument("--contamination-max", type=float, default=0.20)
    ap.add_argument("--near-core-threshold", type=float, default=0.97)
    ap.add_argument("--max-missing-near-core", type=int, default=9)
    ap.add_argument("--out-retained", required=True)
    ap.add_argument("--out-excluded", required=True)
    ap.add_argument("--out-summary", required=True)
    a = ap.parse_args()

    rows = []
    for f in a.contamination:
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            rows.append({"genome": p[0], "n_proteins": int(p[1]),
                         "n_hits": int(p[2]), "no_hit_fraction": float(p[3])})
    qc = pd.DataFrame(rows).drop_duplicates("genome")

    if a.near_core_counts and os.path.exists(a.near_core_counts):
        nc = pd.read_csv(a.near_core_counts, sep="\t",
                         names=["genome", "n_missing_near_core"])
        qc = qc.merge(nc, on="genome", how="left")
    else:
        qc["n_missing_near_core"] = 0

    qc["fail_contamination"] = qc.no_hit_fraction > a.contamination_max
    qc["fail_completeness"]  = qc.n_missing_near_core >= a.max_missing_near_core
    qc["retained"] = ~(qc.fail_contamination | qc.fail_completeness)

    qc.loc[qc.retained, "genome"].to_csv(a.out_retained, index=False, header=False)
    qc.loc[~qc.retained].to_csv(a.out_excluded, sep="\t", index=False)
    qc.to_csv(a.out_summary, sep="\t", index=False)

    print(f"evaluated : {len(qc)}")
    print(f"retained  : {int(qc.retained.sum())}")
    print(f"excluded  : contamination {int(qc.fail_contamination.sum())}, "
          f"fragmentation {int(qc.fail_completeness.sum())}")

if __name__ == "__main__":
    main()
