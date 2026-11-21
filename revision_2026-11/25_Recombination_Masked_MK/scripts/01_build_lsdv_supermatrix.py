#!/usr/bin/env python3
"""Stage 25, step 1 — LSDV-only supermatrix for recombination detection.

Why this stage exists. Restricting the MK ingroup to the 139 non-recombinant
LSDV genomes removes the adaptive signal entirely (alpha = -0.10, CMH
p = 0.91), while the 101 genomes descending from the Node_147 recombination
event retain it (alpha = 0.52, p = 0.028). Removing genomes is the wrong
test, though: it reclassifies sites from polymorphism to divergence and it
throws away three quarters of the standing variation. The right test masks
the recombinant *regions* and keeps every genome.

The manuscript's Gubbins run was performed on the Parsnp whole-genome
alignment, whose coordinates do not map onto the per-gene codon alignments
without going back through the Parsnp LCB blocks. Running Gubbins directly on
the 60-gene supermatrix avoids that entirely: the predicted intervals come
back in supermatrix coordinates, and the partition file maps them onto genes
and codons exactly.

Output: 01_inputs/lsdv_super_nt.fasta, 01_inputs/partitions_nt.txt
"""
import os
import sys
import glob

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import CODON_DIR, CORE_TREE, groups, stage_dir      # noqa: E402
from Bio import SeqIO, Phylo                                     # noqa: E402

D = stage_dir("25_Recombination_Masked_MK")
INP = os.path.join(D, "01_inputs")

G = groups()
tips = [t.name for t in Phylo.read(CORE_TREE, "newick").get_terminals()]
lsdv = [t for t in tips if t in G["LSDV"]]
files = sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta"))
if not files:
    sys.exit(f"no codon alignments under {CODON_DIR}")

seqs = {t: [] for t in lsdv}
parts, pos = [], 1
for f in files:
    og = os.path.basename(f).replace("_codon.fasta", "")
    rec = {x.id: str(x.seq).upper() for x in SeqIO.parse(f, "fasta")}
    L = len(next(iter(rec.values())))
    for t in lsdv:
        seqs[t].append(rec.get(t, "-" * L))
    parts.append(f"{og}\t{pos}\t{pos + L - 1}\t{L}")
    pos += L

with open(f"{INP}/lsdv_super_nt.fasta", "w") as fh:
    for t in lsdv:
        fh.write(f">{t}\n{''.join(seqs[t])}\n")
with open(f"{INP}/partitions_nt.txt", "w") as fh:
    fh.write("orthogroup\tstart\tend\tlength\n")
    fh.write("\n".join(parts) + "\n")
print(f"{len(lsdv)} LSDV genomes x {pos - 1} nt from {len(files)} orthogroups")
print(f"written to {INP}/lsdv_super_nt.fasta")
