#!/usr/bin/env python3
"""Stage 24, step 1 — nucleotide supermatrix for outgroup-free rooting.

Rooting under a non-reversible substitution model needs nucleotide data: the
non-reversible DNA models carry directional information that a reversible
model, and a protein model, do not. The 60 single-copy core orthogroups are
concatenated in a fixed order, with a partition file, and taxa missing from a
gene are filled with gaps.

Output: 01_inputs/super_nt.fasta, 01_inputs/partitions_nt.txt
"""
import os
import sys
import glob

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import CODON_DIR, CORE_TREE, stage_dir          # noqa: E402
from Bio import SeqIO, Phylo                                 # noqa: E402

D = stage_dir("24_Root_Inference")
INP = os.path.join(D, "01_inputs")

tips = [t.name for t in Phylo.read(CORE_TREE, "newick").get_terminals()]
files = sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta"))
if not files:
    sys.exit(f"no codon alignments under {CODON_DIR}")

seqs = {t: [] for t in tips}
parts, pos = [], 1
for f in files:
    og = os.path.basename(f).replace("_codon.fasta", "")
    rec = {x.id: str(x.seq).upper() for x in SeqIO.parse(f, "fasta")}
    L = len(next(iter(rec.values())))
    for t in tips:
        seqs[t].append(rec.get(t, "-" * L))
    parts.append(f"DNA, {og} = {pos}-{pos + L - 1}")
    pos += L

with open(f"{INP}/super_nt.fasta", "w") as fh:
    for t in tips:
        fh.write(f">{t}\n{''.join(seqs[t])}\n")
with open(f"{INP}/partitions_nt.txt", "w") as fh:
    fh.write("\n".join(parts) + "\n")
print(f"{len(tips)} taxa x {pos - 1} nt from {len(files)} orthogroups")
