#!/usr/bin/env python3
"""Summarise BUSTED output and record whether codon models are applicable.

The diagnostics matter as much as the p-values. In a clonal population with
almost no synonymous variation the likelihood surface degenerates, and a
significant result is then a numerical artefact rather than evidence of
selection. This script records, per gene, the number of distinct protein
haplotypes, the proportion of sites with dS approximately zero, and whether the
likelihood-ratio statistic collapsed to zero.

Portable reimplementation of bin/legacy/parse_busted.py.
"""
import argparse, glob, json, os
import numpy as np, pandas as pd

def uniq_proteins(path):
    seqs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return len(set(seqs))

def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-dir", required=True)
    ap.add_argument("--alignments", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-diagnostics", required=True)
    a = ap.parse_args()

    rows = []
    for jf in sorted(glob.glob(os.path.join(a.json_dir, "*_BUSTED.json"))):
        og = os.path.basename(jf).replace("_BUSTED.json", "")
        try:
            d = json.load(open(jf))
        except Exception:
            rows.append({"orthogroup": og, "status": "unparsable"}); continue

        tr = d.get("test results", {})
        lrt, p = tr.get("LRT"), tr.get("p-value")

        fits = d.get("fits", {})
        unc  = fits.get("Unconstrained model", {})
        rd   = unc.get("Rate Distributions", {})
        srv  = rd.get("Synonymous site-to-site rates", {})
        n_srv = len(srv) if isinstance(srv, dict) else np.nan

        omega_mean, omega_max = np.nan, np.nan
        w = rd.get("Test") or rd.get("Shared")
        if isinstance(w, dict):
            vals = [(v.get("omega"), v.get("proportion")) for v in w.values()
                    if isinstance(v, dict)]
            vals = [(o, pr) for o, pr in vals if o is not None and pr is not None]
            if vals:
                omega_mean = float(sum(o * pr for o, pr in vals))
                omega_max  = float(max(o for o, _ in vals))

        mg = fits.get("MG94xREV with separate rates for branch sets", {})
        omega_mg94 = np.nan
        mgrd = mg.get("Rate Distributions", {})
        if isinstance(mgrd, dict):
            for v in mgrd.values():
                if isinstance(v, dict) and "omega" in v:
                    omega_mg94 = float(v["omega"]); break

        aln = os.path.join(a.alignments, f"{og}_codon.fasta")
        n_uniq = uniq_proteins(aln) if os.path.exists(aln) else np.nan

        rows.append({"orthogroup": og, "status": "ok", "LRT": lrt, "p_value": p,
                     "omega_mean": omega_mean, "omega_max": omega_max,
                     "omega_MG94": omega_mg94, "n_SRV_classes": n_srv,
                     "n_unique_protein_seqs": n_uniq,
                     "LRT_is_zero": (lrt == 0) if lrt is not None else None})

    df = pd.DataFrame(rows)
    ok = df.status == "ok"
    if ok.any() and df.loc[ok, "p_value"].notna().any():
        sub = df.loc[ok & df.p_value.notna()]
        df.loc[sub.index, "q_value_BH"] = bh(sub.p_value.values)
        df["significant"] = df.q_value_BH < a.fdr
    df.to_csv(a.out_summary, index=False)

    n = int(ok.sum())
    diag = pd.DataFrame([{
        "genes_evaluated": n,
        "genes_with_LRT_zero": int(df.LRT_is_zero.fillna(False).sum()),
        "median_unique_protein_seqs": float(df.n_unique_protein_seqs.median(skipna=True)),
        "median_omega_MG94": float(df.omega_MG94.median(skipna=True)),
        "genes_significant_FDR": int(df.get("significant", pd.Series(dtype=bool)).sum()),
        "verdict": ("codon models degenerate: interpret with a population-genetic "
                    "framework instead"
                    if n and df.LRT_is_zero.fillna(False).sum() > n / 3
                    else "codon models appear identifiable"),
    }])
    diag.to_csv(a.out_diagnostics, sep="\t", index=False)
    print(diag.to_string(index=False))

if __name__ == "__main__":
    main()
