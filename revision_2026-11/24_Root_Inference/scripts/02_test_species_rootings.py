#!/usr/bin/env python3
"""Stage 24, step 2 — outgroup-free test of the three species-level rootings.

The manuscript polarizes the McDonald-Kreitman test by treating the MRCA of
the SPPV+GTPV clade as the ancestral state for LSDV, which is valid only if
LSDV is the earliest-diverging lineage. An outgroup-rooted phylogeny cannot
test its own rooting assumption, and the published ancient-genome analysis
reports the root as unresolved and prefers SPPV as basal. This step tests the
three rootings directly, using a non-reversible model, which infers the root
from the data rather than from an assumed outgroup.

The same unrooted topology is rooted three ways -- LSDV basal, SPPV basal,
GTPV basal -- and the three rooted trees are compared by the approximately
unbiased test with 10,000 RELL replicates.

Outputs
  01_inputs/{LSDV,SPPV,GTPV}_basal.nwk, three_rootings.nwk
  02_results/au_test_three_rootings.iqtree
  03_summary/species_rootings.csv
"""
import os
import re
import sys
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import CORE_TREE, groups, stage_dir, iqtree_bin   # noqa: E402
import pandas as pd                                            # noqa: E402
from Bio import Phylo                                          # noqa: E402
from Bio.Phylo.Newick import Clade                             # noqa: E402

MODEL = os.environ.get("ROOT_MODEL", "UNREST+G4")   # fully general non-reversible
NREP = int(os.environ.get("ROOT_ZB", "10000"))

D = stage_dir("24_Root_Inference")
INP, RES, SUM = (os.path.join(D, x) for x in ("01_inputs", "02_results", "03_summary"))
ALN = f"{INP}/super_nt.fasta"
if not os.path.exists(ALN):
    sys.exit("run 01_build_nt_supermatrix.py first")

G = groups()


def root_off(tree, target, path):
    """Root `tree` on the branch subtending `target`, keeping the root bifurcating."""
    tree.root_with_outgroup(tree.common_ancestor(sorted(target)))
    kids = list(tree.root.clades)
    inside = [c for c in kids if {x.name for x in c.get_terminals()} <= target]
    outside = [c for c in kids if c not in inside]
    if len(kids) > 2:                 # Biopython can leave a trifurcation
        tree.root.clades = [Clade(branch_length=0.0, clades=inside)] + outside
    Phylo.write(tree, path, "newick")
    return [{k: int(v) for k, v in
             pd.Series([s for x in c.get_terminals()
                        for s in [next((sp for sp, m in G.items()
                                        if x.name in m), "?")]]).value_counts().items()}
            for c in tree.root.clades]


SPECS = [("LSDV_basal", G["SPPV"] | G["GTPV"]),
         ("SPPV_basal", G["SPPV"]),
         ("GTPV_basal", G["GTPV"])]
order = []
with open(f"{INP}/three_rootings.nwk", "w") as zf:
    for name, target in SPECS:
        p = f"{INP}/{name}.nwk"
        comp = root_off(Phylo.read(CORE_TREE, "newick"), set(target), p)
        print(f"{name:<12} root sides: {comp}")
        zf.write(open(p).read().strip() + "\n")
        order.append(name)

pre = f"{RES}/au_test_three_rootings"
subprocess.run([iqtree_bin(), "-s", ALN, "-m", MODEL, "-z",
                f"{INP}/three_rootings.nwk", "-n", "0", "-zb", str(NREP),
                "-au", "--prefix", pre, "-T", "2", "-redo"],
               capture_output=True, text=True)

txt = open(pre + ".iqtree").read()
blk = txt.split("USER TREES")[1]
rows = []
for ln in blk.splitlines():
    m = re.match(r"\s*(\d+)\s+(-[\d.]+)\s+([\d.]+)\s+([\d.e-]+)\s*([+-])\s+"
                 r"([\d.e-]+)\s*[+-]\s+([\d.e-]+)\s*[+-]\s+([\d.e-]+)\s*[+-]\s+"
                 r"([\d.e-]+)\s*([+-])", ln)
    if m:
        rows.append(dict(rooting=order[int(m.group(1)) - 1],
                         logL=float(m.group(2)), deltaL=float(m.group(3)),
                         bp_RELL=float(m.group(4)), c_ELW=float(m.group(8)),
                         p_AU=float(m.group(9)),
                         in_95_set=(m.group(10) == "+")))
df = pd.DataFrame(rows)
df.to_csv(f"{SUM}/species_rootings.csv", index=False)
print(f"\nmodel {MODEL}, {NREP} RELL replicates\n")
print(df.to_string(index=False))
