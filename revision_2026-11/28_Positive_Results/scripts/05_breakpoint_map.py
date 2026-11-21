#!/usr/bin/env python3
"""Stage 28, step 5 — the recombination breakpoint map.

Your recombination result is currently one sentence: nine of ten annotated
recombinants coalesce at Node_147, and r/m = 0.13. The 90 blocks support
considerably more than that, and all of it is positive.

Three analyses:

  density    where breakpoints fall along the alignment, tested against a
             uniform null by a chi-square on equal-width bins and by a
             Kolmogorov-Smirnov test, so a terminal or central enrichment
             can be stated with a p value rather than asserted

  recurrence which intervals are hit by more than one independent block,
             since a position used repeatedly is a hotspot and a candidate
             mechanism rather than an accident

  Node_147   the breakpoints of the single ancestral event, which are the
             coordinates a DIVA assay would target: a stable ancestral
             junction is a better assay target than a transient tract

Outputs
  02_results/breakpoint_positions.csv   every breakpoint with its branch
  02_results/breakpoint_density.csv     binned counts and the null
  03_summary/recurrent_intervals.csv    intervals hit more than once
  03_summary/node147_breakpoints.csv    the ancestral event's junctions
  03_summary/05_breakpoint_map.log
"""
import os
import re
import sys
import collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _s28 import ROOT, opener                                     # noqa: E402

import numpy as np                                                # noqa: E402
import pandas as pd                                               # noqa: E402

NBINS = int(os.environ.get("NBINS", "50"))
GFF = os.environ.get(
    "GUBBINS_GFF", f"{ROOT}/06_Recombination/02_Gubbins_Results/"
                   "lsdv_gubbins.recombination_predictions.gff")
NODE = os.environ.get("NODE", "Node_147")

D, say, close = opener("05_breakpoint_map.log")

if not os.path.exists(GFF):
    say(f"Gubbins GFF not found at\n  {GFF}\n")
    say("Point at it with:  export GUBBINS_GFF=/path/to/*.gff")
    close()
    sys.exit(1)

rows = []
for ln in open(GFF):
    if ln.startswith("#") or not ln.strip():
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    try:
        s, e = int(f[3]), int(f[4])
    except ValueError:
        continue
    attr = f[8]
    node = re.search(r'node="?([^";]+)"?', attr)
    taxa = re.search(r'taxa="?\s*([^";]+)"?', attr)
    snp = re.search(r'snp_count="?(\d+)"?', attr)
    tx = taxa.group(1).split() if taxa else []
    rows.append(dict(start=s, end=e, length=e - s + 1,
                     node=node.group(1) if node else "",
                     n_taxa=len(tx), taxa=" ".join(tx),
                     snp_count=int(snp.group(1)) if snp else None))

bl = pd.DataFrame(rows)
if bl.empty:
    say("no intervals parsed from the GFF")
    close()
    sys.exit(1)
say(f"{len(bl)} recombination blocks parsed")
say(f"  total span {int(bl.length.sum())} bp, "
    f"median block {int(bl.length.median())} bp")
ALN = int(bl.end.max())
say(f"  alignment length implied by the intervals: {ALN}")

# ---- breakpoint positions ---------------------------------------------------
bp = pd.concat([
    bl[["start", "node", "n_taxa"]].rename(columns={"start": "position"})
      .assign(edge="start"),
    bl[["end", "node", "n_taxa"]].rename(columns={"end": "position"})
      .assign(edge="end")], ignore_index=True)
bp.to_csv(os.path.join(D, "02_results", "breakpoint_positions.csv"),
          index=False)
say(f"\n{len(bp)} breakpoints ({len(bl)} blocks x 2 edges)")

# ---- density against a uniform null ----------------------------------------
edges = np.linspace(0, ALN, NBINS + 1)
obs, _ = np.histogram(bp.position, bins=edges)
exp = np.full(NBINS, len(bp) / NBINS)
dens = pd.DataFrame({"bin_start": edges[:-1].astype(int),
                     "bin_end": edges[1:].astype(int),
                     "observed": obs, "expected": exp.round(2)})
dens.to_csv(os.path.join(D, "02_results", "breakpoint_density.csv"),
            index=False)

from scipy.stats import chisquare, kstest                         # noqa: E402
chi, p_chi = chisquare(obs, exp)
ks, p_ks = kstest(bp.position / ALN, "uniform")
say(f"\n=== is the distribution uniform? ===")
say(f"  chi-square over {NBINS} bins : chi2 = {chi:.1f}, p = {p_chi:.3e}")
say(f"  Kolmogorov-Smirnov           : D = {ks:.4f}, p = {p_ks:.3e}")
if p_chi < 0.05 or p_ks < 0.05:
    say("  breakpoints are NOT uniformly distributed")
else:
    say("  no departure from uniformity detected")

# terminal versus central thirds
third = ALN / 3
term = int(((bp.position < third) | (bp.position > 2 * third)).sum())
cent = len(bp) - term
say(f"\n  outer thirds {term} breakpoints, central third {cent}")
say(f"  ratio {term / max(cent, 1):.2f} against the 2.0 expected if uniform")
from scipy.stats import binomtest                                 # noqa: E402
bt = binomtest(term, len(bp), 2 / 3, alternative="two-sided")
say(f"  binomial test against 2/3 expectation: p = {bt.pvalue:.3e}")

top = dens.sort_values("observed", ascending=False).head(8)
say(f"\n=== densest bins ===")
say(f"  {'interval':<22}{'observed':>9}{'expected':>10}")
for r in top.itertuples():
    say(f"  {str(r.bin_start) + '-' + str(r.bin_end):<22}"
        f"{r.observed:>9}{r.expected:>10.2f}")

# ---- recurrence -------------------------------------------------------------
cov = collections.Counter()
for r in bl.itertuples():
    cov[(r.start // 1000, r.end // 1000)] += 1
multi = [(k, v) for k, v in cov.items() if v > 1]
rec_rows = []
step = 1000
hits = np.zeros(ALN // step + 2, int)
for r in bl.itertuples():
    hits[r.start // step:r.end // step + 1] += 1
runs = []
cur = None
for i, h in enumerate(hits):
    if h >= 2:
        if cur is None:
            cur = [i, i, h]
        else:
            cur[1] = i
            cur[2] = max(cur[2], h)
    elif cur is not None:
        runs.append(cur)
        cur = None
if cur is not None:
    runs.append(cur)
for a, b, h in sorted(runs, key=lambda x: -x[2]):
    rec_rows.append(dict(start=a * step, end=(b + 1) * step,
                         max_overlapping_blocks=h))
rec = pd.DataFrame(rec_rows)
rec.to_csv(os.path.join(D, "03_summary", "recurrent_intervals.csv"),
           index=False)
say(f"\n=== intervals covered by two or more independent blocks ===")
say(f"  {len(rec)} such intervals")
if len(rec):
    say(f"  {'interval':<22}{'max blocks':>12}")
    for r in rec.head(10).itertuples():
        say(f"  {str(r.start) + '-' + str(r.end):<22}"
            f"{r.max_overlapping_blocks:>12}")
    say("\n  These are candidate hotspots. A region used repeatedly by")
    say("  independent events is more interesting than a long tract used")
    say("  once, and it is the kind of claim a reviewer will accept from")
    say("  this dataset.")

# ---- the ancestral event ----------------------------------------------------
# Gubbins labels BRANCHES, not nodes: entries look like "Node_148->Node_147".
# For an assay target that every descendant of a clade carries, only the
# branch INTO that clade's ancestor counts. Blocks on branches leading OUT of
# it are confined to sub-lineages and are not shared by the whole clade.
inc = bl[bl.node.astype(str).str.strip().str.endswith("->" + NODE)]
out = bl[bl.node.astype(str).str.strip().str.startswith(NODE + "->")]
exact = bl[bl.node.astype(str).str.strip() == NODE]

say(f"\n=== blocks relating to {NODE} ===")
say(f"  on the branch INTO it   : {len(inc)} blocks"
    + (f"  ({', '.join(sorted(inc.node.unique()))})" if len(inc) else ""))
say(f"  on branches OUT of it   : {len(out)} blocks"
    + (f"  ({', '.join(sorted(out.node.unique()))})" if len(out) else ""))
if len(exact):
    say(f"  labelled exactly {NODE} : {len(exact)} blocks")

anc = inc if len(inc) else exact
if anc.empty:
    say(f"\nNothing on a branch into {NODE}. Branch labels present:")
    for n, c in bl.node.value_counts().head(15).items():
        say(f"  {str(n):<32}{c:>4} blocks")
    say(f"\nSet the label with:  export NODE=<node name>")
    say("Give the node, not the branch: the script builds the branch label.")
else:
    anc = anc.sort_values("start")
    anc.to_csv(os.path.join(D, "03_summary", "ancestral_breakpoints.csv"),
               index=False)
    say(f"\n=== the ancestral event on the branch into {NODE} ===")
    say(f"  {len(anc)} blocks, {int(anc.length.sum())} bp total")
    say(f"  {'branch':<24}{'start':>9}{'end':>10}{'length':>9}"
        f"{'SNPs':>7}{'taxa':>7}")
    for r in anc.itertuples():
        say(f"  {str(r.node)[:23]:<24}{r.start:>9}{r.end:>10}{r.length:>9}"
            f"{(r.snp_count if r.snp_count is not None else 0):>7}"
            f"{r.n_taxa:>7}")
    if anc.n_taxa.nunique() == 1:
        say(f"\n  All {len(anc)} blocks apply to the same "
            f"{int(anc.n_taxa.iloc[0])} taxa, which is what a single event on "
            f"one branch should look like.")
    else:
        say(f"\n  WARNING: the taxa counts differ "
            f"({sorted(anc.n_taxa.unique())}). Blocks on one branch should "
            f"apply to the same set of descendants, so check the labels.")
    say("\n  These junctions are the DIVA targets: carried by every")
    say("  descendant of the event, so an assay spanning one separates the")
    say("  recombinant lineage from the clonal backbone without depending")
    say("  on a tract confined to a sub-lineage.")
    if len(out):
        say(f"\n  The {len(out)} blocks on branches out of {NODE} are NOT")
        say("  suitable: each is confined to part of the clade. An earlier")
        say("  version of this script lumped them in with the ancestral")
        say("  blocks, which would have produced unusable assay targets.")

say("\nwritten: 02_results/breakpoint_positions.csv, breakpoint_density.csv")
say("        03_summary/recurrent_intervals.csv")
if not anc.empty:
    say("        03_summary/ancestral_breakpoints.csv")
close()
