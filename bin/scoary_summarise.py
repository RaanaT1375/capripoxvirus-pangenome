#!/usr/bin/env python3
"""Per-trait summary of a Scoary run, including traits with no significant hit.

An association is retained only if it satisfies BOTH the per-trait Bonferroni
correction and the permutation empirical p-value. Either filter alone is
substantially more permissive. Convergent (Tier 1) associations additionally
require at least the given number of independent supporting pairs with net
positive support, which is what separates adaptation from lineage expansion.

Portable reimplementation of bin/legacy/build_S6_pangwas_table.py.
"""
import argparse, glob, os
import numpy as np, pandas as pd

def pick(cols, *cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c.lower() in low: return low[c.lower()]
    for c in cands:
        for lc, orig in low.items():
            if c.lower() in lc: return orig
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--traits", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--bonferroni", type=float, default=0.05)
    ap.add_argument("--empirical", type=float, default=0.05)
    ap.add_argument("--min-pairs", type=int, default=3)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-pairs", required=True)
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.results_dir, "*.results.csv"))) \
            or sorted(glob.glob(os.path.join(a.results_dir, "*.csv")))
    traits = pd.read_csv(a.traits, index_col=0)

    mat = pd.read_csv(a.matrix, index_col=0, low_memory=False)
    bincols = [c for c in mat.columns
               if set(pd.unique(mat[c].astype(str))) <= {"0", "1"}]
    m = mat[bincols].astype(int) if bincols else mat.notna().astype(int)
    rs, ng = m.sum(axis=1), m.shape[1]
    n_variable = int(((rs > 0) & (rs < ng)).sum())

    rows, sig_all = [], []
    for tname in traits.columns:
        col = traits[tname].dropna()
        rec = {"trait": tname,
               "n_genomes_positive": int((col == 1).sum()),
               "n_genomes_negative": int((col == 0).sum()),
               "n_orthogroups_tested": n_variable,
               "statistical_method": "Fisher's exact with phylogenetic pairwise comparisons",
               "multiple_testing_correction":
                   f"per-trait Bonferroni < {a.bonferroni} AND empirical p < {a.empirical}",
               "n_significant": 0, "n_passing_bonferroni": 0, "n_tier1_convergent": 0,
               "best_orthogroup": None, "best_odds_ratio": None,
               "best_empirical_p": None, "best_bonferroni_p": None,
               "best_supporting_pairs": None, "best_opposing_pairs": None,
               "outcome": "no significant association"}

        match = [f for f in files
                 if os.path.basename(f).split(".")[0].lower() == str(tname).lower()] \
                or [f for f in files if str(tname).lower() in os.path.basename(f).lower()]

        if match:
            d = pd.read_csv(match[0])
            g   = pick(d.columns, "Gene", "Orthogroup")
            emp = pick(d.columns, "Empirical_p")
            bon = pick(d.columns, "Bonferroni_p")
            orr = pick(d.columns, "Odds_ratio")
            sup = pick(d.columns, "Max_supporting_pairs")
            opp = pick(d.columns, "Max_opposing_pairs")

            mask = pd.Series(True, index=d.index)
            if emp: mask &= pd.to_numeric(d[emp], errors="coerce") < a.empirical
            if bon: mask &= pd.to_numeric(d[bon], errors="coerce") < a.bonferroni
            sig = d[mask]

            rec["n_significant"] = int(len(sig))
            if bon:
                rec["n_passing_bonferroni"] = int(
                    (pd.to_numeric(d[bon], errors="coerce") < a.bonferroni).sum())
            if sup and opp:
                s = pd.to_numeric(sig[sup], errors="coerce")
                o = pd.to_numeric(sig[opp], errors="coerce")
                rec["n_tier1_convergent"] = int(((s >= a.min_pairs) & ((s - o) > 0)).sum())

            if len(sig):
                key = emp or bon
                best = sig.loc[pd.to_numeric(sig[key], errors="coerce").idxmin()]
                rec.update(best_orthogroup=best.get(g), best_odds_ratio=best.get(orr),
                           best_empirical_p=best.get(emp), best_bonferroni_p=best.get(bon),
                           best_supporting_pairs=best.get(sup),
                           best_opposing_pairs=best.get(opp),
                           outcome=f"{len(sig)} significant gene-trait pair(s)")
                keep = [c for c in [g, orr, emp, bon, sup, opp] if c]
                tmp = sig[keep].copy(); tmp.insert(0, "trait", tname)
                sig_all.append(tmp)
        else:
            rec["outcome"] = "no Scoary output located for this trait"

        rows.append(rec)
        print(f"{tname:<26} pos={rec['n_genomes_positive']:>4} "
              f"sig={rec['n_significant']:>3} tier1={rec['n_tier1_convergent']:>3}")

    pd.DataFrame(rows).to_csv(a.out_summary, index=False)
    if sig_all:
        s = pd.concat(sig_all, ignore_index=True)
        s.to_csv(a.out_pairs, index=False)
        gc = [c for c in s.columns if c.lower() in ("gene", "orthogroup")]
        print(f"\npairs={len(s)}  traits={s.trait.nunique()}"
              + (f"  orthogroups={s[gc[0]].nunique()}" if gc else ""))
    else:
        pd.DataFrame().to_csv(a.out_pairs, index=False)

if __name__ == "__main__":
    main()
