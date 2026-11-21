#!/usr/bin/env python3
"""Concatenate per-orthogroup alignments into a supermatrix with a partition file.

Sequence headers are reduced to the genome identifier so that every alignment
contributes the same taxon set. Any orthogroup that is not present in every
genome is reported and skipped, since a partitioned analysis requires a
complete matrix.
"""
import argparse, os, re
from collections import defaultdict
import pandas as pd

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name: seqs[name] = "".join(buf)
    return seqs

def genome_id(header):
    # Prokka locus tags look like GENOME_00123; strip the trailing gene number
    return re.sub(r"_\d+$", "", header)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alignments", nargs="+", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-partitions", required=True)
    ap.add_argument("--out-stats", required=True)
    a = ap.parse_args()

    parsed, taxa = {}, None
    for path in sorted(a.alignments):
        og = os.path.basename(path).split(".")[0]
        s = {genome_id(h): v for h, v in read_fasta(path).items()}
        parsed[og] = s
        taxa = set(s) if taxa is None else taxa & set(s)

    taxa = sorted(taxa)
    if not taxa:
        raise SystemExit("ERROR: no taxon is present in every alignment.")

    blocks, start, stats = [], 1, []
    cat = defaultdict(list)
    for og in sorted(parsed):
        s = parsed[og]
        if not all(t in s for t in taxa):
            stats.append({"orthogroup": og, "included": False, "length": 0})
            continue
        L = len(next(iter(s.values())))
        if any(len(s[t]) != L for t in taxa):
            stats.append({"orthogroup": og, "included": False, "length": 0})
            continue
        for t in taxa:
            cat[t].append(s[t])
        blocks.append((og, start, start + L - 1))
        stats.append({"orthogroup": og, "included": True, "length": L})
        start += L

    with open(a.out_fasta, "w") as fh:
        for t in taxa:
            fh.write(f">{t}\n{''.join(cat[t])}\n")

    with open(a.out_partitions, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for og, s, e in blocks:
            fh.write(f"    charset {og} = {s}-{e};\n")
        fh.write("end;\n")

    st = pd.DataFrame(stats)
    st.to_csv(a.out_stats, sep="\t", index=False)
    print(f"taxa: {len(taxa)}   partitions: {len(blocks)}   "
          f"supermatrix length: {start - 1}")
    dropped = st.loc[~st.included, "orthogroup"].tolist()
    if dropped:
        print(f"skipped (incomplete): {len(dropped)} -> {dropped[:5]}")

if __name__ == "__main__":
    main()
