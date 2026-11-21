#!/usr/bin/env python3
"""
Build the balanced sub-sampling summary table (Supplementary File S7).

For the full analysis and for each of the ten balanced sub-sampling replicates,
extract the MPPA ancestral state(s) assigned to the root node and the marginal
posterior probability of each continent at that node.

Reads only; writes two files into 07_Phylogeography/summary/.

Run from: /cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
"""
import os, glob, re, sys
import pandas as pd

os.chdir("/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis")
BASE = "07_Phylogeography"
OUT = os.path.join(BASE, "summary")
os.makedirs(OUT, exist_ok=True)

MARG = "marginal_probabilities.character_Continent.model_F81.tab"
COMB = "combined_ancestral_states.tab"

runs = [("full", os.path.join(BASE, "pastml_output_v3"))]
for d in sorted(glob.glob(os.path.join(BASE, "sensitivity_analysis", "rep*_pastml")),
                key=lambda x: int(re.search(r"rep(\d+)", x).group(1))):
    runs.append((os.path.basename(d).replace("_pastml", ""), d))
print(f"Runs found: {[r[0] for r in runs]}\n")


def load(path):
    for sep in ["\t", ","]:
        try:
            d = pd.read_csv(path, sep=sep)
            if d.shape[1] > 1:
                return d
        except Exception:
            pass
    return None


def find_root(df):
    """Return the row of df corresponding to the tree root."""
    idcol = df.columns[0]
    ids = df[idcol].astype(str)
    for pat in [r"^root$", r"^ROOT$", r"root"]:
        hit = df[ids.str.match(pat, case=False, na=False)]
        if len(hit):
            return hit.iloc[0], idcol
    # PastML often names the root n1 or n0
    for name in ["n1", "n0", "ROOT"]:
        hit = df[ids == name]
        if len(hit):
            return hit.iloc[0], idcol
    return None, idcol


rows = []
first = True
for label, d in runs:
    mpath, cpath = os.path.join(d, MARG), os.path.join(d, COMB)
    rec = {"replicate": label}

    marg = load(mpath) if os.path.exists(mpath) else None
    if marg is not None:
        if first:
            print(f"Columns in {MARG}: {list(marg.columns)}")
            print(f"First 5 node ids: {marg[marg.columns[0]].head().tolist()}\n")
            first = False
        r, idcol = find_root(marg)
        if r is not None:
            rec["root_node_id"] = r[idcol]
            for c in marg.columns[1:]:
                v = pd.to_numeric(pd.Series([r[c]]), errors="coerce").iloc[0]
                if v == v:
                    rec[f"P_{c}"] = round(float(v), 4)
        else:
            rec["root_node_id"] = "NOT FOUND"

    comb = load(cpath) if os.path.exists(cpath) else None
    if comb is not None:
        idcol = comb.columns[0]
        ids = comb[idcol].astype(str)
        sel = comb[ids.str.contains("root", case=False, na=False)]
        if not len(sel) and "root_node_id" in rec:
            sel = comb[ids == str(rec["root_node_id"])]
        if len(sel):
            statecol = comb.columns[-1] if len(comb.columns) > 1 else idcol
            states = sorted(set(str(x) for x in sel[statecol].dropna()))
            rec["MPPA_root_state"] = " | ".join(states)
            rec["n_states_at_root"] = len(states)

    probs = {k[2:]: v for k, v in rec.items() if k.startswith("P_")}
    if probs:
        top = max(probs, key=probs.get)
        rec["highest_probability_continent"] = top
        rec["highest_probability"] = probs[top]
        rec["resolved"] = "ambiguous" if rec.get("n_states_at_root", 1) > 1 else top

    rows.append(rec)

res = pd.DataFrame(rows)
front = [c for c in ["replicate", "root_node_id", "MPPA_root_state", "n_states_at_root",
                     "highest_probability_continent", "highest_probability", "resolved"]
         if c in res.columns]
res = res[front + [c for c in res.columns if c not in front]]
res.to_csv(f"{OUT}/S7_subsampling_root_states.csv", index=False)

print(res.to_string(index=False))
print(f"\nWrote {OUT}/S7_subsampling_root_states.csv")

reps = res[res.replicate != "full"]
if len(reps) and "MPPA_root_state" in reps.columns:
    print("\n--- summary across the ten replicates ---")
    print(reps["MPPA_root_state"].value_counts().to_string())
    if "n_states_at_root" in reps.columns:
        amb = int((reps.n_states_at_root > 1).sum())
        print(f"\nambiguous (more than one state retained): {amb}/{len(reps)}")
    afr = int(reps["MPPA_root_state"].fillna("").str.contains("Africa").sum())
    print(f"replicates in which Africa was retained at the root: {afr}/{len(reps)}")
    print("\nCHECK against the manuscript: 5/10 ambiguous, Africa resolved in 0/10.")
