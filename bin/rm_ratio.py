#!/usr/bin/env python3
"""Recombination-to-mutation ratio (r/m) from Gubbins output.

r/m is the number of SNPs introduced by recombination divided by the number
introduced by point mutation, summed over the phylogeny. It is reported here
relative to the alignable core alignment, not to the full reference genome, so
the denominator of every percentage is stated explicitly.
"""
import argparse, re
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="*.recombination_predictions.gff")
    ap.add_argument("--embl", default=None, help="*.branch_base_reconstruction.embl")
    ap.add_argument("--alignment", required=True, help="core alignment FASTA")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    aln_len = 0
    with open(a.alignment) as fh:
        seq = []
        for line in fh:
            if line.startswith(">"):
                if seq: break
                continue
            seq.append(line.strip())
        aln_len = len("".join(seq))

    blocks, snps_recomb = [], 0
    for line in open(a.gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        start, end = int(f[3]), int(f[4])
        blocks.append((start, end))
        m = re.search(r'snp_count[= ]"?(\d+)', f[8])
        if m:
            snps_recomb += int(m.group(1))

    merged, covered = [], 0
    for s, e in sorted(blocks):
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    covered = sum(e - s + 1 for s, e in merged)

    snps_mutation = 0
    if a.embl:
        try:
            txt = open(a.embl).read()
            snps_mutation = len(re.findall(r'FT\s+variation', txt))
        except OSError:
            pass

    rm = (snps_recomb / snps_mutation) if snps_mutation else float("nan")

    pd.DataFrame([{
        "n_blocks": len(blocks),
        "merged_blocks": len(merged),
        "bp_covered": covered,
        "core_alignment_bp": aln_len,
        "fraction_core_covered": round(covered / aln_len, 4) if aln_len else None,
        "snps_from_recombination": snps_recomb,
        "snps_from_mutation": snps_mutation,
        "r_over_m": round(rm, 4) if rm == rm else None,
    }]).to_csv(a.out, sep="\t", index=False)

    print(f"blocks {len(blocks)}  covered {covered} bp "
          f"({100*covered/aln_len:.1f}% of the {aln_len} bp core alignment)")
    print(f"r/m = {rm:.4f}" if rm == rm else "r/m not computed (no mutation count)")

if __name__ == "__main__":
    main()
