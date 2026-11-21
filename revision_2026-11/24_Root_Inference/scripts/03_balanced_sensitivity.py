#!/usr/bin/env python3
"""Stage 24, step 3 — balanced-sampling sensitivity for the root test.

The dataset holds 240 LSDV, 34 SPPV and 16 GTPV genomes. Rooting under a
non-reversible model can be pulled toward the densely sampled group, so the
three species-level rootings are re-tested on replicates that take an equal
number of genomes from each species (16, the size of the smallest group).

Outputs
  02_results/balanced_roottest.csv
  03_summary/balanced_summary.csv
"""
import os
import re
import random
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import CORE_TREE, groups, stage_dir, iqtree_bin   # noqa: E402
import pandas as pd                                            # noqa: E402
from Bio import SeqIO, Phylo                                   # noqa: E402
from Bio.Phylo.Newick import Clade                             # noqa: E402

MODEL = os.environ.get("ROOT_MODEL", "UNREST+G4")
NREP = int(os.environ.get("ROOT_ZB", "10000"))
NSUB = int(os.environ.get("ROOT_REPS", "10"))

D = stage_dir("24_Root_Inference")
INP, RES, SUM = (os.path.join(D, x) for x in ("01_inputs", "02_results", "03_summary"))
W = os.path.join(RES, "balanced"); os.makedirs(W, exist_ok=True)

G = groups()
N = min(len(v) for v in G.values())
rec = {r.id: str(r.seq) for r in SeqIO.parse(f"{INP}/super_nt.fasta", "fasta")}


def root_off(tree, target, path):
    tree.root_with_outgroup(tree.common_ancestor(sorted(target)))
    kids = list(tree.root.clades)
    inside = [c for c in kids if {x.name for x in c.get_terminals()} <= target]
    outside = [c for c in kids if c not in inside]
    if len(kids) > 2:
        tree.root.clades = [Clade(branch_length=0.0, clades=inside)] + outside
    Phylo.write(tree, path, "newick")


rows = []
for rep in range(NSUB):
    rnd = random.Random(1000 + rep)
    keep = sum((rnd.sample(sorted(G[s]), N) for s in ["LSDV", "SPPV", "GTPV"]), [])
    ks = set(keep)
    aln = f"{W}/r{rep}.fa"
    with open(aln, "w") as fh:
        for k in keep:
            fh.write(f">{k}\n{rec[k]}\n")

    t = Phylo.read(CORE_TREE, "newick")
    for tip in [x for x in t.get_terminals() if x.name not in ks]:
        t.prune(tip)
    pruned = f"{W}/r{rep}_pruned.nwk"
    Phylo.write(t, pruned, "newick")

    order = ["LSDV_basal", "SPPV_basal", "GTPV_basal"]
    targets = {"LSDV_basal": (G["SPPV"] | G["GTPV"]) & ks,
               "SPPV_basal": G["SPPV"] & ks,
               "GTPV_basal": G["GTPV"] & ks}
    zf = f"{W}/r{rep}_three.nwk"
    with open(zf, "w") as fh:
        for name in order:
            p = f"{W}/r{rep}_{name}.nwk"
            root_off(Phylo.read(pruned, "newick"), targets[name], p)
            fh.write(open(p).read().strip() + "\n")

    pre = f"{W}/r{rep}_z"
    subprocess.run([iqtree_bin(), "-s", aln, "-m", MODEL, "-z", zf, "-n", "0",
                    "-zb", str(NREP), "-au", "--prefix", pre, "-T", "2", "-redo"],
                   capture_output=True, text=True)
    blk = open(pre + ".iqtree").read().split("USER TREES")[1]
    vals = {}
    for ln in blk.splitlines():
        m = re.match(r"\s*(\d+)\s+(-[\d.]+)\s+([\d.]+)\s+([\d.e-]+)\s*[+-]\s+"
                     r"([\d.e-]+)\s*[+-]\s+([\d.e-]+)\s*[+-]\s+([\d.e-]+)\s*[+-]\s+"
                     r"([\d.e-]+)\s*([+-])", ln)
        if m:
            vals[order[int(m.group(1)) - 1]] = dict(
                logL=float(m.group(2)), deltaL=float(m.group(3)),
                bp_RELL=float(m.group(4)), c_ELW=float(m.group(7)),
                p_AU=float(m.group(8)), in_95_set=(m.group(9) == "+"))
    if not vals:
        print(f"replicate {rep}: could not parse the test table"); continue
    best = min(vals, key=lambda k: vals[k]["deltaL"])
    print(f"rep {rep}: best={best:<11} " + "  ".join(
        f"{k.split('_')[0]} ELW={v['c_ELW']:.3f} AU={v['p_AU']:.3f}"
        for k, v in vals.items()))
    for k, v in vals.items():
        rows.append(dict(replicate=rep, rooting=k, best=(k == best), **v))

df = pd.DataFrame(rows)
df.to_csv(f"{RES}/balanced_roottest.csv", index=False)
s = df.groupby("rooting").agg(times_best=("best", "sum"),
                              mean_ELW=("c_ELW", "mean"),
                              mean_AU=("p_AU", "mean"),
                              in_95_set=("in_95_set", "sum")).round(3)
s.to_csv(f"{SUM}/balanced_summary.csv")
print(f"\n=== {N} genomes per species, {NSUB} replicates ===")
print(s.to_string())
