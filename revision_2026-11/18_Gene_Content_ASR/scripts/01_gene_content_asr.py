#!/usr/bin/env python3
"""Stage 18 — ancestral gene-content reconstruction.

Reviewer 2, comment 1: the frequency distribution of orthogroups describes
presence/absence across sampled genomes but says nothing about the timing or
direction of gene gain and loss, so "transient acquisition" cannot be read
off it. This stage supplies the evolutionary evidence the reviewer asks for.

The 206-orthogroup presence/absence matrix is treated as a binary alignment
(one site per orthogroup) on the fixed core-genome ML tree, and IQ-TREE
reconstructs marginal ancestral gene content at every internal node under a
two-state substitution model. Each branch is then classified as a gain
(0 -> 1) or a loss (1 -> 0) and, separately, as terminal or internal.

The whole reconstruction is repeated on the 237 GenBank genomes alone,
because stage 19 shows the low-frequency compartment is inflated by the
de novo SRA assemblies.

Outputs
  02_results/gene_content_events_all290.csv
  02_results/gene_content_events_genbank237.csv
  03_summary/gene_content_summary.csv
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import (CORE_TREE, ANNOT_ALL, load_pa, stage_dir,   # noqa: E402
                     run_iqtree, iqtree_bin, is_sra)
import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402
from Bio import Phylo                                            # noqa: E402

MODEL = os.environ.get("PA_MODEL", "GTR2+FO")
D = stage_dir("18_Gene_Content_ASR")
LOG = open(os.path.join(D, "03_summary", "01_gene_content_asr.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


say(f"IQ-TREE: {iqtree_bin()} | binary model: {MODEL}")
PA = load_pa()
ann = pd.read_csv(ANNOT_ALL).set_index("Orthogroup")["product"]
say(f"presence/absence matrix: {PA.shape[0]} orthogroups x {PA.shape[1]} genomes")


def classify(freq, n):
    return pd.cut(freq / n, [-.01, .15, .95, .99, 1.01],
                  labels=["cloud", "shell", "soft_core", "core"])


def reconstruct(tag, cols):
    sub = PA[cols]
    sub = sub.loc[sub.sum(axis=1) > 0]            # drop orthogroups absent here
    n = len(cols)
    cls = classify(sub.sum(axis=1), n)
    work = os.path.join(D, "02_results", tag)
    os.makedirs(work, exist_ok=True)

    t = Phylo.read(CORE_TREE, "newick")
    for tip in [x for x in t.get_terminals() if x.name not in set(cols)]:
        t.prune(tip)
    tre = f"{work}/tree.nwk"
    Phylo.write(t, tre, "newick")

    tips = [x.name for x in t.get_terminals()]
    phy = f"{work}/pa.phy"
    with open(phy, "w") as fh:
        fh.write(f"{len(tips)} {sub.shape[0]}\n")
        for x in tips:
            fh.write(f"{x}  {''.join(sub[x].astype(str))}\n")

    run_iqtree(["-s", phy, "-te", tre, "-m", MODEL, "-asr", "-T", "1",
                "--prefix", f"{work}/pa", "--quiet", "-redo"],
               log=os.path.join(D, "logs", f"iqtree_{tag}.log"))
    if not os.path.exists(f"{work}/pa.state"):
        raise SystemExit(f"IQ-TREE produced no .state for {tag}; see logs/")

    st = pd.read_csv(f"{work}/pa.state", sep="\t", comment="#")
    st["call"] = (st["p_1"] >= 0.5).astype(int)   # resolves '-' ties by posterior
    anc = st.pivot(index="Node", columns="Site", values="call")

    t2 = Phylo.read(f"{work}/pa.treefile", "newick")
    for cl in t2.get_nonterminals():
        if cl.name is None and cl.confidence is not None:
            cl.name = str(cl.confidence)

    def state(node, site):
        return (int(sub.iloc[site - 1][node.name]) if node.is_terminal()
                else int(anc.loc[node.name, site]))

    rows = []
    for site in range(1, sub.shape[0] + 1):
        og = sub.index[site - 1]
        gt = gi = lt = li = 0
        for parent in t2.get_nonterminals():
            ps = state(parent, site)
            for child in parent.clades:
                cs = state(child, site)
                if cs == ps:
                    continue
                term = child.is_terminal()
                if ps == 0 and cs == 1:
                    gt += term; gi += (not term)
                else:
                    lt += term; li += (not term)
        rows.append(dict(Orthogroup=og, n_present=int(sub.loc[og].sum()),
                         pangenome_class=str(cls[og]),
                         gain_terminal=gt, gain_internal=gi,
                         loss_terminal=lt, loss_internal=li,
                         product=str(ann.get(og, ""))[:80]))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(D, "02_results",
                            f"gene_content_events_{tag}.csv"), index=False)

    agg = res.groupby("pangenome_class", observed=True)[
        ["gain_terminal", "gain_internal", "loss_terminal",
         "loss_internal"]].sum()
    agg["n_orthogroups"] = res.groupby("pangenome_class", observed=True).size()
    agg["gains"] = agg.gain_terminal + agg.gain_internal
    agg["losses"] = agg.loss_terminal + agg.loss_internal
    agg["pct_gains_terminal"] = (
        100 * agg.gain_terminal / agg.gains.replace(0, np.nan)).round(1)
    agg.insert(0, "dataset", tag)
    say(f"\n=== {tag}: {n} genomes, {sub.shape[0]} orthogroups present ===")
    say(agg[["dataset", "n_orthogroups", "gains", "losses", "gain_terminal",
             "gain_internal", "pct_gains_terminal"]].to_string())
    return agg.reset_index()


cols_all = list(PA.columns)
cols_gb = [c for c in PA.columns if not is_sra(c)]
out = pd.concat([reconstruct("all290", cols_all),
                 reconstruct("genbank237", cols_gb)], ignore_index=True)
out.to_csv(os.path.join(D, "03_summary", "gene_content_summary.csv"),
           index=False)
say("\nwritten: 03_summary/gene_content_summary.csv")
LOG.close()
