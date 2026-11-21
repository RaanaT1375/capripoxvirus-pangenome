#!/usr/bin/env python3
"""Stage 28, step 1 — is selective constraint structured by functional class?

Your inter-species omega estimates are valid: the codon models failed only
WITHIN LSDV, where there is no synonymous variation to fit them to. Across
the genus the 60 single-copy core genes gave a median omega of 0.169 with
every gene below one, and that number is currently reported only as a
preamble to a negative result. It is a positive finding in its own right if
the variation among genes is structured rather than arbitrary.

This classifies each core gene by functional class from its product
annotation and tests whether omega differs across classes
(Kruskal-Wallis, then pairwise Mann-Whitney with Holm correction).

The expected pattern is a gradient: transcription and replication machinery
under the strongest constraint, host-interaction genes the weakest. If that
holds, it is one sentence in the abstract and a figure.

Override the automatic classification by writing
01_inputs/functional_categories.csv with columns orthogroup,category.

Outputs
  02_results/gene_constraint.csv       per-gene omega, product, category
  03_summary/constraint_by_class.csv   per-class medians and the tests
  04_figures/constraint_by_class.pdf   boxplot, if matplotlib is present
  03_summary/01_constraint_structure.log
"""
import os
import sys
import glob

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _s28 import (ROOT, classify, find_selection_table, opener)   # noqa: E402

import pandas as pd                                               # noqa: E402

D, say, close = opener("01_constraint_structure.log")
CAT_OVERRIDE = os.path.join(D, "01_inputs", "functional_categories.csv")

# ---- omega per gene ---------------------------------------------------------
src = os.environ.get("SELECTION_TABLE")
og_col = os.environ.get("OG_COL")
om_col = os.environ.get("OMEGA_COL")
if not src:
    src, og_col, om_col = find_selection_table()
    if not src:
        say("Could not find a per-gene selection table with an omega column.")
        say("Candidates inspected:")
        for f in (om_col or [])[:25]:
            say("  " + os.path.relpath(f, ROOT))
        say("\nPoint at the right file directly:")
        say("  SELECTION_TABLE=/path/table.csv OG_COL=Orthogroup "
            "OMEGA_COL=omega_MG94 \\")
        say("    python3 scripts/01_constraint_structure.py")
        close()
        sys.exit(1)
say(f"selection table: {os.path.relpath(src, ROOT) if src.startswith(ROOT) else src}")
say(f"  orthogroup column = {og_col}, omega column = {om_col}")

t = pd.read_csv(src)
df = t[[og_col, om_col]].rename(columns={og_col: "orthogroup",
                                         om_col: "omega"})
df = df.dropna().drop_duplicates(subset="orthogroup")
df["omega"] = pd.to_numeric(df.omega, errors="coerce")
df = df.dropna()
say(f"{len(df)} genes with an omega estimate")
say(f"  median {df.omega.median():.4f}, range "
    f"{df.omega.min():.4f}-{df.omega.max():.4f}")
if (df.omega >= 1).any():
    say(f"  note: {int((df.omega >= 1).sum())} gene(s) at or above 1")

# ---- product annotation -----------------------------------------------------
prod = {}
for pat in [f"{ROOT}/00_Metadata/*.csv", f"{ROOT}/99_Supplementary/**/*.csv"]:
    for f in sorted(glob.glob(pat, recursive=True)):
        try:
            a = pd.read_csv(f, nrows=3)
        except Exception:                                         # noqa: BLE001
            continue
        cols = {c.lower(): c for c in a.columns}
        oc = next((cols[c] for c in cols if c in ("orthogroup", "gene")), None)
        pc = next((cols[c] for c in cols
                   if c in ("product", "annotation", "description")), None)
        if oc and pc:
            full = pd.read_csv(f)
            for r in full[[oc, pc]].dropna().itertuples(index=False):
                prod.setdefault(str(r[0]), str(r[1]))
if prod:
    say(f"{len(prod)} product annotations collected")
else:
    say("no product annotations found; every gene will be Unassigned unless "
        "you supply 01_inputs/functional_categories.csv")
df["product"] = df.orthogroup.map(prod).fillna("")

if os.path.exists(CAT_OVERRIDE):
    ov = pd.read_csv(CAT_OVERRIDE)
    m = dict(zip(ov.iloc[:, 0].astype(str), ov.iloc[:, 1].astype(str)))
    auto = df["product"].map(classify)
    df["category"] = [m.get(o, a) for o, a in zip(df.orthogroup, auto)]
    say(f"using the override for {sum(o in m for o in df.orthogroup)} genes")
else:
    df["category"] = df["product"].map(classify)

df = df.sort_values("omega")
df.to_csv(os.path.join(D, "02_results", "gene_constraint.csv"), index=False)

say("\n=== genes per class ===")
for c, n in df.category.value_counts().items():
    sub = df[df.category == c]
    say(f"  {c:<22}{n:>4}   median omega {sub.omega.median():.4f}")

un = int((df.category == "Unassigned").sum())
if un:
    say(f"\n{un} genes could not be classified from their annotation. "
        f"They are kept in the table but excluded from the test.")
    for r in df[df.category == "Unassigned"].itertuples():
        say(f"  {r.orthogroup:<13}{str(r.product)[:60]}")
    say("Classify them by hand in 01_inputs/functional_categories.csv if you "
        "want them included.")

# ---- tests ------------------------------------------------------------------
test = df[df.category != "Unassigned"]
groups = [g.omega.values for _, g in test.groupby("category") if len(g) >= 3]
names = [c for c, g in test.groupby("category") if len(g) >= 3]
rows = []
if len(groups) < 2:
    say("\nfewer than two classes with at least 3 genes; no test performed")
else:
    from scipy.stats import kruskal, mannwhitneyu
    H, p = kruskal(*groups)
    say(f"\n=== Kruskal-Wallis across {len(groups)} classes "
        f"({', '.join(names)}) ===")
    say(f"  H = {H:.3f}, p = {p:.3e}")
    rows.append(dict(comparison="Kruskal-Wallis (all classes)",
                     statistic=H, p_raw=p, p_holm=p, n1=len(test), n2=None))
    if p < 0.05:
        say("  constraint differs significantly among functional classes")
    else:
        say("  no significant difference among classes")

    pw = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            U, pp = mannwhitneyu(groups[i], groups[j], alternative="two-sided")
            pw.append((names[i], names[j], U, pp, len(groups[i]),
                       len(groups[j])))
    pw.sort(key=lambda x: x[3])
    m = len(pw)
    say("\n=== pairwise Mann-Whitney, Holm-corrected ===")
    prev = 0.0
    for k, (a, b, U, pp, n1, n2) in enumerate(pw):
        ph = min(1.0, max(prev, pp * (m - k)))
        prev = ph
        flag = " *" if ph < 0.05 else ""
        say(f"  {a:<22} vs {b:<22} p = {pp:.4f}  Holm {ph:.4f}{flag}")
        pw_row = dict(comparison=f"{a} vs {b}", statistic=U, p_raw=pp,
                      p_holm=ph, n1=n1, n2=n2)
        rows.append(pw_row)

summ = (test.groupby("category").omega
        .agg(["count", "median", "mean", "min", "max"]).reset_index())
summ.to_csv(os.path.join(D, "03_summary", "constraint_by_class.csv"),
            index=False)
pd.DataFrame(rows).to_csv(
    os.path.join(D, "03_summary", "constraint_tests.csv"), index=False)
say("\nwritten: 03_summary/constraint_by_class.csv, constraint_tests.csv")
say("        02_results/gene_constraint.csv")
say("\nThe R script R/constraint.R draws the figure from "
    "02_results/gene_constraint.csv.")
close()
