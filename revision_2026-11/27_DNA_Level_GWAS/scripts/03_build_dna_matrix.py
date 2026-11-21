#!/usr/bin/env python3
"""Stage 27, step 3 — build the DNA-level presence/absence matrix.

Same calling rule as stage 20: an orthogroup is present in a genome when its
best tblastn hit covers at least MIN_COV of the query at MIN_PID identity,
summing non-overlapping HSPs within a contig so that a gene split across a
contig break is not scored absent.

Writes the matrix in the same layout as the Scoary input it replaces, so
step 4 can run Scoary with the arguments used for the original analysis and
nothing else changes between the two runs.

Orthogroups with no query sequence are carried over from the annotation
matrix unchanged, and flagged in the comparison so they are not mistaken
for DNA-level calls.

Outputs
  02_results/presence_dna.csv              orthogroup x genome, 0/1
  02_results/annotation_vs_dna.csv         per-orthogroup comparison
  03_summary/gene_presence_absence_dna.csv Scoary input
  03_summary/core_tree_pruned.newick       tree pruned to matching taxa
  03_summary/traits_matched.csv            traits subset to matching taxa
  03_summary/03_build_dna_matrix.log
"""
import os
import sys


def _locate_common():
    env = os.environ.get("COMMON_DIR")
    if env and os.path.exists(os.path.join(env, "_common.py")):
        return env
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(8):
        if os.path.exists(os.path.join(d, "_common.py")):
            return d
        p = os.path.dirname(d)
        if p == d:
            break
        d = p
    sys.exit("could not find _common.py above " + os.path.abspath(__file__)
             + "\nSet COMMON_DIR to the directory that holds it.")


sys.path.insert(0, _locate_common())
from _common import (PA_CSV, TRAITS, CORE_TREE, META_COLS,        # noqa: E402
                     stage_dir, load_pa, is_sra)

import pandas as pd                                               # noqa: E402
from Bio import Phylo                                             # noqa: E402

MIN_COV = float(os.environ.get("MIN_COV", "0.80"))
MIN_PID = float(os.environ.get("MIN_PID", "80"))

D = stage_dir("27_DNA_Level_GWAS")
HITS = os.path.join(D, "02_results", "tblastn_hits.tsv")
LOG = open(os.path.join(D, "03_summary", "03_build_dna_matrix.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


if not os.path.exists(HITS):
    sys.exit(f"{HITS} not found; run 02_run_tblastn.sbatch first")

say(f"thresholds: query coverage >= {MIN_COV:.0%}, identity >= {MIN_PID:.0f}%")

ann_pa = load_pa()
genomes = list(ann_pa.columns)
say(f"annotation matrix: {ann_pa.shape[0]} orthogroups x {len(genomes)} genomes")

cols = ["qseqid", "sseqid", "pident", "length", "qlen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"]
bl = pd.read_csv(HITS, sep="\t", names=cols)
bl["genome"] = bl.sseqid.str.split("|").str[0]
bl = bl[bl.pident >= MIN_PID]
say(f"{len(bl)} HSPs pass the identity filter, "
    f"{bl.qseqid.nunique()} orthogroups, {bl.genome.nunique()} genomes")


def covered(sub):
    """Non-overlapping query positions covered by the HSPs of one contig."""
    tot = 0
    cs = ce = None
    for s, e in sorted(zip(sub.qstart, sub.qend)):
        if cs is None:
            cs, ce = s, e
        elif s <= ce + 1:
            ce = max(ce, e)
        else:
            tot += ce - cs + 1
            cs, ce = s, e
    if cs is not None:
        tot += ce - cs + 1
    return tot


rec = []
for (q, g, _c), sub in bl.groupby([bl.qseqid, bl.genome, bl.sseqid]):
    rec.append((q, g, covered(sub) / sub.qlen.iloc[0]))
cov = (pd.DataFrame(rec, columns=["Orthogroup", "genome", "qcov"])
       .sort_values("qcov").groupby(["Orthogroup", "genome"]).last())
say(f"{len(cov)} orthogroup-genome pairs with a hit")

present = (cov.qcov >= MIN_COV).unstack(fill_value=False)
dna = pd.DataFrame(0, index=ann_pa.index, columns=genomes, dtype=int)
shared_g = [g for g in genomes if g in present.columns]
shared_o = [o for o in ann_pa.index if o in present.index]
dna.loc[shared_o, shared_g] = present.loc[shared_o, shared_g].astype(int).values

no_query = [o for o in ann_pa.index if o not in present.index]
if no_query:
    say(f"\n{len(no_query)} orthogroups had no tblastn hit at all; "
        f"carried over from the annotation matrix:")
    say("  " + ", ".join(no_query[:20]) + (" ..." if len(no_query) > 20 else ""))
    dna.loc[no_query] = ann_pa.loc[no_query]

dna.to_csv(os.path.join(D, "02_results", "presence_dna.csv"))

# ---- comparison -------------------------------------------------------------
sra = [g for g in genomes if is_sra(g)]
gb = [g for g in genomes if not is_sra(g)]
rows = []
for og in ann_pa.index:
    a, d = ann_pa.loc[og], dna.loc[og]
    rows.append(dict(
        Orthogroup=og, ann_n=int(a.sum()), dna_n=int(d.sum()),
        gained=int(((d == 1) & (a == 0)).sum()),
        lost=int(((d == 0) & (a == 1)).sum()),
        dna_SRA=int(d[sra].sum()), dna_GB=int(d[gb].sum()),
        carried_over=og in no_query))
cmp = pd.DataFrame(rows)
cmp.to_csv(os.path.join(D, "02_results", "annotation_vs_dna.csv"), index=False)

tested = cmp[~cmp.carried_over]
unchanged = int(((tested.gained == 0) & (tested.lost == 0)).sum())
say("\n=== how much the matrix changed ===")
say(f"orthogroups tested at the DNA level : {len(tested)}")
say(f"  unchanged                         : {unchanged}")
say(f"  gained genomes                    : {int((tested.gained > 0).sum())}")
say(f"  lost genomes                      : {int((tested.lost > 0).sum())}")
say(f"total annotation-based presences    : {int(ann_pa.values.sum())}")
say(f"total DNA-level presences           : {int(dna.values.sum())}")

big = tested[tested.gained >= 50].sort_values("gained", ascending=False)
if len(big):
    say(f"\n{len(big)} orthogroups gained 50 or more genomes "
        f"(annotation dropout on a large scale):")
    say(f"  {'Orthogroup':<13}{'ann':>6}{'DNA':>6}{'gained':>8}")
    for r in big.head(25).itertuples():
        say(f"  {r.Orthogroup:<13}{r.ann_n:>6}{r.dna_n:>6}{r.gained:>8}")

# compartment shift
def classify(m):
    f = m.sum(axis=1) / m.shape[1]
    return pd.cut(f, [-.01, .15, .95, .99, 1.01],
                  labels=["cloud", "shell", "soft-core", "core"])


say("\n=== pan-genome compartments ===")
say(f"  {'compartment':<12}{'annotation':>12}{'DNA':>8}")
ca, cd = classify(ann_pa).value_counts(), classify(dna).value_counts()
for k in ["core", "soft-core", "shell", "cloud"]:
    say(f"  {k:<12}{int(ca.get(k, 0)):>12}{int(cd.get(k, 0)):>8}")

# ---- Scoary inputs ----------------------------------------------------------
tree = Phylo.read(CORE_TREE, "newick")
tips = {t.name for t in tree.get_terminals() if t.name}
traits = pd.read_csv(TRAITS, index_col=0)
keep = sorted(set(dna.columns) & tips & set(traits.index))
say(f"\n=== Scoary inputs ===")
say(f"matrix genomes {len(dna.columns)}, tree tips {len(tips)}, "
    f"trait rows {len(traits)}  ->  {len(keep)} in common")
if len(keep) < 50:
    sys.exit("too few genomes shared between matrix, tree and traits; "
             "check that the names match")

for t in tree.get_terminals():
    if t.name and t.name not in keep:
        tree.prune(t)
Phylo.write(tree, os.path.join(D, "03_summary", "core_tree_pruned.newick"),
            "newick")
traits.loc[keep].to_csv(os.path.join(D, "03_summary", "traits_matched.csv"))

# mirror the original Scoary input layout exactly
orig = pd.read_csv(PA_CSV, low_memory=False)
meta = [c for c in META_COLS if c in orig.columns and c != "Gene"]
out = pd.concat([
    pd.DataFrame({"Gene": dna.index}, index=dna.index),
    orig.set_index("Gene")[meta].reindex(dna.index),
    dna[keep],
], axis=1)
out.to_csv(os.path.join(D, "03_summary", "gene_presence_absence_dna.csv"),
           index=False)
say(f"written: 03_summary/gene_presence_absence_dna.csv "
    f"({dna.shape[0]} orthogroups x {len(keep)} genomes)")
say("written: 03_summary/core_tree_pruned.newick")
say("written: 03_summary/traits_matched.csv")
say("\nnext:  sbatch scripts/04_run_scoary.sbatch")
LOG.close()
