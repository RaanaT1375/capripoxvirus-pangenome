#!/usr/bin/env python3
"""Stage 17, step 1 — marginal ML ancestral-state reconstruction per core gene.

Reviewer 2, comment 3: consensus codons from extant SPPV/GTPV do not
necessarily represent the ancestral state. This step reconstructs ancestral
states explicitly under a reversible nucleotide model, so substitutions,
reversions and parallel mutations inside the outgroup lineages are modelled
rather than averaged over.

For each of the 60 single-copy core orthogroups the per-gene tree is rooted
on the SPPV+GTPV clade and IQ-TREE reconstructs marginal ancestral states at
every internal node (-asr). Branch lengths are re-optimised on the fixed
topology (-te), exactly as HyPhy does in stages 10/13/15/16.

Output: 02_results/asr_runs/<OG>.state and <OG>.treefile
"""
import os
import sys
import glob

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import (CODON_DIR, TREE_DIR, groups, stage_dir,  # noqa: E402
                     run_iqtree, iqtree_bin)
import pandas as pd                                            # noqa: E402
from Bio import Phylo                                          # noqa: E402

MODEL = os.environ.get("ASR_MODEL", "GTR+F+G4")

D = stage_dir("17_ASR_MK")
OUT = os.path.join(D, "02_results", "asr_runs")
os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(D, "logs", "01_run_asr.log")

G = groups()
OUTG = G["SPPV"] | G["GTPV"]
print(f"IQ-TREE: {iqtree_bin()}\nmodel  : {MODEL}\noutput : {OUT}\n", flush=True)

rows = []
for aln in sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta")):
    og = os.path.basename(aln).replace("_codon.fasta", "")
    tre = f"{TREE_DIR}/{og}_pruned.nwk"
    if not os.path.exists(tre):
        rows.append((og, "missing_tree", 0, 0)); continue

    t = Phylo.read(tre, "newick")
    tips = {x.name for x in t.get_terminals()}
    og_tips = sorted(tips & OUTG)
    if len(og_tips) < 2:
        rows.append((og, "too_few_outgroup", len(tips), len(og_tips))); continue

    # root on the outgroup clade so that the LSDV stem is an edge at the root
    t.root_with_outgroup(t.common_ancestor(og_tips))
    rooted = f"{OUT}/{og}_rooted.nwk"
    Phylo.write(t, rooted, "newick")

    pre = f"{OUT}/{og}"
    run_iqtree(["-s", aln, "-te", rooted, "-m", MODEL, "-asr",
                "-T", "1", "--prefix", pre, "--quiet", "-redo"], log=LOG)
    ok = os.path.exists(pre + ".state")
    rows.append((og, "ok" if ok else "FAILED", len(tips), len(og_tips)))
    print(f"{og}\t{rows[-1][1]}\t{len(tips)} tips ({len(og_tips)} outgroup)",
          flush=True)

st = pd.DataFrame(rows, columns=["Orthogroup", "status", "n_tips",
                                 "n_outgroup_tips"])
st.to_csv(os.path.join(D, "03_summary", "asr_run_status.csv"), index=False)
n_ok = int((st.status == "ok").sum())
print(f"\ndone: {n_ok}/{len(st)} genes reconstructed")
if n_ok != len(st):
    print("WARNING: some genes failed - see", LOG)
    sys.exit(1)
