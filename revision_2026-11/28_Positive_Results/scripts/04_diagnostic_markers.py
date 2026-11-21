#!/usr/bin/env python3
"""Stage 28, step 4 — a minimal diagnostic SNP panel from the core alignment.

Entirely positive, translational, and the data are already on disk. For each
grouping of interest the script finds alignment columns at which one group
is fixed for one base and every other genome is fixed for a different base,
which are the columns a diagnostic assay would target.

Four groupings:
  species        LSDV / SPPV / GTPV
  vaccine        vaccine-derived versus field, within LSDV
  recombinant    Node_147 descendants versus the clonal backbone, within LSDV
  baps           each fastBAPS cluster against the rest

Columns are reported in core-alignment coordinates and, where the Parsnp
XMFA is available, translated to NC_003027.1 genome coordinates so they can
be used to design primers. A perfectly diagnostic site tolerating no
exceptions is reported first; sites tolerating one or two exceptions are
reported separately, since a marker that fails in one of 240 genomes is
often still useful and a strict rule can return nothing.

Outputs
  02_results/diagnostic_sites.csv      every site, with its grouping
  03_summary/diagnostic_panel.csv      the recommended minimal panel
  03_summary/04_diagnostic_markers.log
"""
import os
import sys
import collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _s28 import ROOT, opener                                     # noqa: E402
from _common import TRAITS, groups, vaccine_ids                   # noqa: E402

import numpy as np                                                # noqa: E402
import pandas as pd                                               # noqa: E402
from Bio import SeqIO                                             # noqa: E402

MAX_EXC = int(os.environ.get("MAX_EXCEPTIONS", "2"))
CORE = os.environ.get(
    "CORE_ALN", f"{ROOT}/06_Recombination/01_Parsnp_Alignment/"
                "parsnp_core_alignment.fasta")

D, say, close = opener("04_diagnostic_markers.log")

if not os.path.exists(CORE):
    say(f"core alignment not found at\n  {CORE}\n")
    say("Point at it with:  export CORE_ALN=/path/to/alignment.fasta")
    close()
    sys.exit(1)

recs = {r.id.split("|")[0].split(".fna")[0]: str(r.seq).upper()
        for r in SeqIO.parse(CORE, "fasta")}
say(f"core alignment: {len(recs)} sequences")
L = len(next(iter(recs.values())))
say(f"alignment length: {L} columns")

G = groups()
VACC = vaccine_ids()
names = [n for n in recs if n in set().union(*G.values())]
if len(names) < 50:
    # names in the alignment may carry suffixes; try a looser match
    allg = set().union(*G.values())
    remap = {}
    for n in recs:
        for g in allg:
            if n.startswith(g) or g.startswith(n):
                remap[n] = g
                break
    recs = {remap[n]: s for n, s in recs.items() if n in remap}
    names = list(recs)
    say(f"matched {len(names)} sequence names to the metadata by prefix")
say(f"usable genomes: {len(names)}")

M = np.array([list(recs[n]) for n in names], dtype="U1")
idx = {n: i for i, n in enumerate(names)}

# ---- groupings --------------------------------------------------------------
lsdv = [n for n in names if n in G["LSDV"]]
tests = []
for sp in ("LSDV", "SPPV", "GTPV"):
    ing = [n for n in names if n in G[sp]]
    if ing:
        tests.append((f"species:{sp}", ing, [n for n in names if n not in G[sp]]))

vac = [n for n in lsdv if n in VACC]
fld = [n for n in lsdv if n not in VACC]
if vac and fld:
    tests.append(("vaccine:LSDV_vaccine", vac, fld))

if os.path.exists(TRAITS):
    tr = pd.read_csv(TRAITS, index_col=0)
    for col, label in [("Ancestral_recombinant", "recombinant:Node_147"),
                       ("No_recombination", "recombinant:clonal_backbone")]:
        if col in tr.columns:
            pos = [n for n in lsdv if n in tr.index and tr.loc[n, col] == 1]
            neg = [n for n in lsdv if n in tr.index and tr.loc[n, col] == 0]
            if len(pos) >= 3 and len(neg) >= 3:
                tests.append((label, pos, neg))
    for c in [c for c in tr.columns if c.startswith("BAPS_cluster")]:
        pos = [n for n in names if n in tr.index and tr.loc[n, c] == 1]
        neg = [n for n in names if n in tr.index and tr.loc[n, c] == 0]
        if len(pos) >= 3 and len(neg) >= 3:
            tests.append((f"baps:{c}", pos, neg))
else:
    say(f"note: {TRAITS} not found; only species and vaccine groupings tested")

say(f"\n{len(tests)} groupings to scan")

# ---- scan -------------------------------------------------------------------
ACGT = set("ACGT")
rows = []
for label, pos, neg in tests:
    pi = [idx[n] for n in pos]
    ni = [idx[n] for n in neg]
    found = 0
    for c in range(L):
        colp = M[pi, c]
        coln = M[ni, c]
        pc = collections.Counter(x for x in colp if x in ACGT)
        nc = collections.Counter(x for x in coln if x in ACGT)
        if not pc or not nc:
            continue
        pb, pn = pc.most_common(1)[0]
        exc_in = len(colp) - pn
        shared = nc.get(pb, 0)
        if shared > MAX_EXC or exc_in > MAX_EXC:
            continue
        nb = nc.most_common(1)[0][0]
        if nb == pb:
            continue
        rows.append(dict(grouping=label, column=c + 1, group_base=pb,
                         other_base=nb, n_group=len(pos), n_other=len(neg),
                         exceptions_in_group=exc_in,
                         exceptions_outside=shared,
                         perfect=(exc_in == 0 and shared == 0)))
        found += 1
    say(f"  {label:<34}{found:>6} candidate sites "
        f"({sum(1 for r in rows if r['grouping'] == label and r['perfect'])}"
        f" perfect)")

if not rows:
    say("\nno diagnostic sites found. Try relaxing the tolerance:")
    say("  MAX_EXCEPTIONS=5 python3 scripts/04_diagnostic_markers.py")
    close()
    sys.exit(1)

sites = pd.DataFrame(rows)

# ---- translate to reference coordinates where possible ----------------------
XMFA = f"{ROOT}/06_Recombination/01_Parsnp_Alignment/parsnp.xmfa"
if os.path.exists(XMFA):
    say(f"\nreference coordinates available from {os.path.basename(XMFA)};")
    say("the stage 25 step 4 translation can be reused to convert columns to")
    say("NC_003027.1 positions. Not performed here to keep this step")
    say("independent of that stage; the column numbers are unambiguous.")
else:
    say(f"\n{os.path.basename(XMFA)} not present, so sites are reported in")
    say("core-alignment coordinates only.")

sites.to_csv(os.path.join(D, "02_results", "diagnostic_sites.csv"),
             index=False)

# ---- minimal panel ----------------------------------------------------------
panel = []
for label, sub in sites.groupby("grouping"):
    best = sub.sort_values(["perfect", "exceptions_in_group",
                            "exceptions_outside"],
                           ascending=[False, True, True]).head(3)
    panel.append(best)
panel = pd.concat(panel)
panel.to_csv(os.path.join(D, "03_summary", "diagnostic_panel.csv"),
             index=False)

say(f"\n=== recommended panel: up to 3 sites per grouping ===")
say(f"  {'grouping':<34}{'column':>8}{'group':>7}{'other':>7}  exceptions")
for r in panel.itertuples():
    say(f"  {r.grouping:<34}{r.column:>8}{r.group_base:>7}{r.other_base:>7}  "
        f"{r.exceptions_in_group} in / {r.exceptions_outside} out"
        f"{'  PERFECT' if r.perfect else ''}")

n_perf = int(sites.perfect.sum())
say(f"\n{len(sites)} candidate sites in total, {n_perf} perfectly diagnostic")
say("\nwritten: 02_results/diagnostic_sites.csv")
say("        03_summary/diagnostic_panel.csv")
say("\nBefore reporting these, check that each site falls outside the")
say("recombinant tracts from stage 06: a marker inside a region that moves")
say("between lineages is not a stable assay target. Step 5 produces the")
say("interval list needed for that check.")
close()
