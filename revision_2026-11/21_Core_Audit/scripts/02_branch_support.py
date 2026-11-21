#!/usr/bin/env python3
"""Stage 21, step 2 — branch support for Figure 3.

Reviewer 2, comment 4: support values must be shown on the Figure 3
phylogeny. Before they can be, one ambiguity has to be resolved: the Results
text describes Figure 3A as the IQ-TREE supermatrix tree (60 genes,
UFBoot 1000) while the figure caption describes the RAxML tree built inside
Gubbins on the recombination-filtered alignment. Those are two different
trees and only the first carries support values.

This script reports, for both trees, how many internal splits they share, so
the right one can be chosen deliberately:

  * if Figure 3A is the IQ-TREE tree, the .contree already has UFBoot values
    and nothing needs rerunning - only the caption has to be corrected
  * if Figure 3A is the Gubbins tree, support has to be generated on the
    recombination-filtered alignment; use 03_ufboot_gubbins.sbatch

It also writes the IQ-TREE support values in an iTOL-ready file.

Outputs
  02_results/tree_split_comparison.txt
  02_results/itol_branch_support.txt
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import CORE_CONTREE, ROOT, stage_dir                 # noqa: E402
from Bio import Phylo                                             # noqa: E402

GUBBINS = f"{ROOT}/06_Recombination/02_Gubbins_Results/lsdv_gubbins.final_tree_renamed.nwk"
D = stage_dir("21_Core_Audit")


def splits(tree):
    tips = {t.name for t in tree.get_terminals()}
    out = {}
    for cl in tree.get_nonterminals():
        s = frozenset(t.name for t in cl.get_terminals())
        if 1 < len(s) < len(tips):
            key = min(s, frozenset(tips - s), key=lambda x: (len(x), sorted(x)[0]))
            out[key] = cl
    return out, tips


lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


iq = Phylo.read(CORE_CONTREE, "newick")
iq_splits, iq_tips = splits(iq)
say(f"IQ-TREE consensus tree : {len(iq_tips)} tips, {len(iq_splits)} internal splits")
say(f"  {CORE_CONTREE}")

sup = [c.confidence for c in iq.get_nonterminals() if c.confidence is not None]
if sup:
    say(f"  UFBoot values present on {len(sup)} nodes "
        f"(median {sorted(sup)[len(sup)//2]:.0f}, "
        f"{sum(1 for s in sup if s >= 95)} nodes >= 95)")
else:
    say("  WARNING: no support values found on the consensus tree")

if os.path.exists(GUBBINS):
    gb = Phylo.read(GUBBINS, "newick")
    gb_splits, gb_tips = splits(gb)
    say(f"\nGubbins/RAxML tree     : {len(gb_tips)} tips, "
        f"{len(gb_splits)} internal splits")
    say(f"  {GUBBINS}")
    g_sup = [c.confidence for c in gb.get_nonterminals()
             if c.confidence is not None]
    say(f"  support values on this tree: "
        f"{len(g_sup) if g_sup else 0} nodes"
        f"{'' if g_sup else '  <- this is why Figure 3 has none'}")
    shared = set(iq_splits) & set(gb_splits)
    say(f"\nshared tips   : {len(iq_tips & gb_tips)}")
    say(f"shared splits : {len(shared)} "
        f"({100*len(shared)/max(len(gb_splits),1):.1f}% of the Gubbins splits)")

    # the fraction alone is misleading: what matters is whether the two trees
    # disagree about well-supported nodes or only about unresolved ones
    conf = {k: (c.confidence if c.confidence is not None else None)
            for k, c in iq_splits.items()}
    hi = [k for k, v in conf.items() if v is not None and v >= 95]
    lo = [k for k, v in conf.items() if v is not None and v < 95]
    say(f"\n  splits with UFBoot >= 95 : {len(hi):>4}  shared with Gubbins "
        f"{len(set(hi) & shared):>4} "
        f"({100*len(set(hi) & shared)/max(len(hi),1):.0f}%)")
    say(f"  splits with UFBoot <  95 : {len(lo):>4}  shared with Gubbins "
        f"{len(set(lo) & shared):>4} "
        f"({100*len(set(lo) & shared)/max(len(lo),1):.0f}%)")
    say(f"  splits with no value     : {len(conf) - len(hi) - len(lo):>4}")
    say("\nIf disagreement is concentrated in the low-support splits, the two")
    say("trees agree wherever the data resolve anything and the cleanest fix is")
    say("to show the IQ-TREE tree, which already carries UFBoot, and correct the")
    say("caption. Report honestly that most internal LSDV nodes are weakly")
    say("supported - that is consistent with the near-absent within-LSDV")
    say("diversity reported in stages 09 and 13, and none of the manuscript's")
    say("claims depend on the internal LSDV topology.")
    say("Only if well-supported splits genuinely conflict is a rerun needed:")
    say("then use 03_ufboot_gubbins.sbatch.")
else:
    say(f"\nGubbins tree not found at {GUBBINS}")

with open(os.path.join(D, "02_results", "tree_split_comparison.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")

# ---- iTOL dataset: UFBoot values as branch labels --------------------------
out = os.path.join(D, "02_results", "itol_branch_support.txt")
with open(out, "w") as fh:
    fh.write("DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\tUFBoot\n"
             "COLOR\t#444444\nMARGIN\t2\nSHOW_INTERNAL\t1\nDATA\n")
    n = 0
    for cl in iq.get_nonterminals():
        if cl.confidence is None:
            continue
        tips = sorted(t.name for t in cl.get_terminals())
        if len(tips) < 2:
            continue
        fh.write(f"{tips[0]}|{tips[-1]}\t{cl.confidence:.0f}\t-1\t#444444"
                 f"\tnormal\t1\t0\n")
        n += 1
say(f"\nwritten: 02_results/itol_branch_support.txt ({n} labelled nodes)")
with open(os.path.join(D, "02_results", "tree_split_comparison.txt"), "a") as fh:
    fh.write(lines[-1] + "\n")
