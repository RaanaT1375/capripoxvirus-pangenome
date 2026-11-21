#!/usr/bin/env python3
"""Stage 28, step 6 — turn 4,290 diagnostic sites into a usable panel.

Step 4 found plenty of perfectly diagnostic columns, but its "recommended
panel" was not usable as it stood, for three reasons this step fixes.

  1. Every site it picked sat in the first few hundred columns. With
     thousands of sites tied at zero exceptions, the sort returned whichever
     came first. A panel should be spread along the genome so that one
     assembly gap or one recombinant tract cannot take out every marker.

  2. None was checked against the recombination intervals. A marker inside a
     tract that moves between lineages is not a stable assay target, and the
     Node_147 markers are precisely the ones at risk, since the event
     defining that lineage is what moves.

  3. BAPS cluster 2 is SPPV and cluster 3 is GTPV, so those groupings
     duplicate the species rows exactly (1,220 and 729 sites in both). They
     are dropped to avoid reporting the same markers twice.

Selection rule: keep perfectly diagnostic sites only, drop any site inside
or within FLANK bp of a recombination interval, then pick N_PER sites per
grouping spaced as evenly along the alignment as the survivors allow.

Outputs
  03_summary/panel_final.csv            the recommended panel
  02_results/sites_recombination_ok.csv all perfect sites outside tracts
  03_summary/06_panel_filter.log
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _s28 import ROOT, opener                                     # noqa: E402

import numpy as np                                                # noqa: E402
import pandas as pd                                               # noqa: E402

N_PER = int(os.environ.get("N_PER", "5"))
FLANK = int(os.environ.get("FLANK", "200"))
GFF = os.environ.get(
    "GUBBINS_GFF", f"{ROOT}/06_Recombination/02_Gubbins_Results/"
                   "lsdv_gubbins.recombination_predictions.gff")
# BAPS clusters that duplicate a species grouping
DROP = {"baps:BAPS_cluster_2": "species:SPPV",
        "baps:BAPS_cluster_3": "species:GTPV"}

D, say, close = opener("06_panel_filter.log")
SITES = os.path.join(D, "02_results", "diagnostic_sites.csv")

if not os.path.exists(SITES):
    say(f"{SITES} not found; run 04_diagnostic_markers.py first")
    close()
    sys.exit(1)

s = pd.read_csv(SITES)
say(f"{len(s)} candidate sites across {s.grouping.nunique()} groupings")

# ---- drop duplicated groupings ---------------------------------------------
for dup, same in DROP.items():
    if dup in set(s.grouping) and same in set(s.grouping):
        a = set(s.column[s.grouping == dup])
        b = set(s.column[s.grouping == same])
        ov = len(a & b) / max(len(a | b), 1)
        say(f"  {dup} and {same}: {ov:.0%} of sites identical "
            f"-> dropping {dup}")
        s = s[s.grouping != dup]

perf = s[s.perfect].copy()
say(f"\n{len(perf)} perfectly diagnostic sites retained")

# ---- recombination intervals ------------------------------------------------
iv = []
if os.path.exists(GFF):
    for ln in open(GFF):
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split("\t")
        if len(f) < 5:
            continue
        try:
            iv.append((int(f[3]), int(f[4])))
        except ValueError:
            continue
    say(f"{len(iv)} recombination intervals read from the Gubbins GFF")
else:
    say(f"WARNING: {GFF} not found. Sites CANNOT be screened against")
    say("recombinant tracts, so the panel below is provisional.")

if iv:
    aln = max(e for _, e in iv)
    masked = np.zeros(aln + FLANK + 2, bool)
    for a, b in iv:
        masked[max(0, a - FLANK):min(len(masked), b + FLANK + 1)] = True
    say(f"masked span including {FLANK} bp flanks: "
        f"{int(masked.sum())} of {aln} columns "
        f"({100 * masked.sum() / aln:.1f}%)")
    inside = perf.column.apply(lambda c: bool(masked[c])
                               if c < len(masked) else False)
    say(f"\nsites inside or near a recombinant tract: {int(inside.sum())} "
        f"of {len(perf)} ({100 * inside.mean():.1f}%)")
    for g, sub in perf.groupby("grouping"):
        n_in = int(inside[sub.index].sum())
        say(f"  {g:<32}{n_in:>6} of {len(sub):<6} removed")
    perf = perf[~inside]
    say(f"\n{len(perf)} sites survive the recombination screen")

perf.to_csv(os.path.join(D, "02_results", "sites_recombination_ok.csv"),
            index=False)

# ---- spread along the alignment --------------------------------------------
def spread(sub, n):
    """Pick n sites as evenly spaced along the alignment as possible."""
    sub = sub.sort_values("column")
    if len(sub) <= n:
        return sub
    lo, hi = sub.column.min(), sub.column.max()
    targets = np.linspace(lo, hi, n)
    picked, used = [], set()
    for t in targets:
        cand = sub[~sub.index.isin(used)]
        if cand.empty:
            break
        j = (cand.column - t).abs().idxmin()
        picked.append(j)
        used.add(j)
    return sub.loc[picked].sort_values("column")


rows = []
for g, sub in perf.groupby("grouping"):
    rows.append(spread(sub, N_PER))
if not rows:
    say("\nno sites survive. Relax the screen:")
    say("  FLANK=0 python3 scripts/06_panel_filter.py")
    close()
    sys.exit(1)
panel = pd.concat(rows).sort_values(["grouping", "column"])
panel.to_csv(os.path.join(D, "03_summary", "panel_final.csv"), index=False)

say(f"\n=== final panel: up to {N_PER} sites per grouping, "
    f"outside recombinant tracts, spread along the alignment ===")
say(f"  {'grouping':<32}{'column':>9}{'group':>7}{'other':>7}")
for g, sub in panel.groupby("grouping"):
    for r in sub.itertuples():
        say(f"  {r.grouping:<32}{r.column:>9}{r.group_base:>7}"
            f"{r.other_base:>7}")
    say(f"  {'':<32}{'spacing:':>9} "
        f"{', '.join(str(int(x)) for x in np.diff(sub.column.values))}")

say(f"\n{len(panel)} sites in the panel across "
    f"{panel.grouping.nunique()} groupings")

missing = sorted(set(s.grouping) - set(panel.grouping))
if missing:
    say(f"\ngroupings with no usable site: {', '.join(missing)}")
    say("Report these as negative: no stable diagnostic marker exists for")
    say("them in the core alignment, which for the vaccine grouping is")
    say("expected, since vaccine strains are not a single lineage.")

say("\nwritten: 03_summary/panel_final.csv")
say("        02_results/sites_recombination_ok.csv")
close()
