#!/usr/bin/env python3
"""MK v2: frequency cutoffs, CMH stratified test, per-gene DoS sign test."""
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import permutations
from scipy.stats import fisher_exact, chi2, binomtest, wilcoxon
import pandas as pd, numpy as np, glob, os, collections

ALN = "10_Selection_Pressure/01_codon_alignments"
VALID = set("ACGT"); _tr = {}

def tr(c):
    if c not in _tr: _tr[c] = str(Seq(c).translate())
    return _tr[c]

def path_counts(a, b):
    pos = [i for i in range(3) if a[i] != b[i]]
    if not pos: return 0.0, 0.0
    np_, s_t, n_t = 0, 0.0, 0.0
    for order in set(permutations(pos)):
        cur, s, n, ok = a, 0, 0, True
        for p in order:
            nxt = cur[:p] + b[p] + cur[p+1:]
            if tr(nxt) == "*": ok = False; break
            if tr(nxt) == tr(cur): s += 1
            else: n += 1
            cur = nxt
        if ok: np_ += 1; s_t += s; n_t += n
    if np_ == 0:
        d = len(pos)
        return (d, 0.0) if tr(a) == tr(b) else (0.0, d)
    return s_t/np_, n_t/np_

def clean(c): return len(c) == 3 and set(c) <= VALID and tr(c) != "*"

def scan_gene(path, ingroup, outgroup, min_cov=0.90, out_cons=0.90):
    """Return per-variant polymorphism records + fixed-difference totals."""
    recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}
    ing = [recs[i] for i in ingroup if i in recs]
    out = [recs[i] for i in outgroup if i in recs]
    if not ing or not out: return None
    L = len(next(iter(recs.values()))) // 3
    poly, Dn, Ds = [], 0.0, 0.0
    for k in range(L):
        sl = slice(3*k, 3*k+3)
        ic = [c for c in (s[sl] for s in ing) if clean(c)]
        oc = [c for c in (s[sl] for s in out) if clean(c)]
        if len(ic) < min_cov*len(ing) or len(oc) < min_cov*len(out): continue
        icc, occ = collections.Counter(ic), collections.Counter(oc)
        anc, an = occ.most_common(1)[0]
        if an/len(oc) < out_cons: continue
        major = icc.most_common(1)[0][0]
        if len(icc) > 1:
            for v, cnt in icc.items():
                if v == major: continue
                s, n = path_counts(major, v)
                poly.append((s, n, cnt/len(ic)))
        elif major != anc:
            s, n = path_counts(anc, major); Ds += s; Dn += n
    return poly, Dn, Ds, len(ing), len(out), L

def counts_at(poly, Dn, Ds, cutoff):
    Ps = sum(s for s, n, f in poly if f >= cutoff)
    Pn = sum(n for s, n, f in poly if f >= cutoff)
    return Pn, Ps, Dn, Ds

def cmh(tables):
    """Cochran-Mantel-Haenszel. tables: list of [[Dn,Ds],[Pn,Ps]] (rounded)."""
    A = E = V = 0.0; num = den = 0.0
    for (a, b), (c, d) in tables:
        n = a+b+c+d
        if n == 0 or (a+b) == 0 or (c+d) == 0 or (a+c) == 0 or (b+d) == 0: continue
        A += a; E += (a+b)*(a+c)/n
        V += (a+b)*(c+d)*(a+c)*(b+d)/(n*n*(n-1)) if n > 1 else 0
        num += a*d/n; den += b*c/n
    if V <= 0 or den == 0: return np.nan, np.nan, 0
    stat = (abs(A-E)-0.5)**2 / V
    return chi2.sf(stat, 1), num/den, len(tables)

if __name__ == "__main__":
    sp = pd.read_csv("00_Metadata/species_assignment.csv")
    grp = {s: sp.loc[sp.species == s, "Name"].tolist() for s in ["LSDV","SPPV","GTPV"]}
    ann = pd.read_csv("10_Selection_Pressure/04_summary_tables/orthogroup_annotations.csv")
    CUTS = [0.0, 0.01, 0.05, 0.10, 0.15]

    for ogsp in ["SPPV", "GTPV"]:
        scans = {}
        for f in sorted(glob.glob(f"{ALN}/*_codon.fasta")):
            og = os.path.basename(f).replace("_codon.fasta", "")
            r = scan_gene(f, grp["LSDV"], grp[ogsp])
            if r: scans[og] = r

        print(f"\n{'='*80}\n=== LSDV vs {ogsp}  ({len(scans)} ژن) ===")
        print(f"\n--- α و NI به تفکیک آستانهٔ فراوانی ---")
        print(f"{'cutoff':>7} {'Pn':>8} {'Ps':>8} {'Dn':>8} {'Ds':>8} "
              f"{'NI':>7} {'alpha':>7} {'CMH_p':>10} {'MH_OR':>7}")
        for cut in CUTS:
            tot = np.zeros(4); tabs = []
            for og, (poly, Dn, Ds, *_ ) in scans.items():
                Pn, Ps, Dn_, Ds_ = counts_at(poly, Dn, Ds, cut)
                tot += [Pn, Ps, Dn_, Ds_]
                tabs.append([[round(Dn_), round(Ds_)], [round(Pn), round(Ps)]])
            Pn, Ps, Dn, Ds = tot
            NI = (Pn/Ps)/(Dn/Ds) if Ps and Ds and Dn else np.nan
            p_cmh, mhor, _ = cmh(tabs)
            print(f"{cut:>7.2f} {Pn:>8.1f} {Ps:>8.1f} {Dn:>8.1f} {Ds:>8.1f} "
                  f"{NI:>7.3f} {1-NI:>7.3f} {p_cmh:>10.3e} {mhor:>7.3f}")

        # جدول کامل به‌ازای ژن در آستانهٔ ۵٪ (بدون حذف ژن‌های بدون تفاوت ثابت)
        allrows = []
        for og, (poly, Dn, Ds, nin, nout, L) in scans.items():
            Pn, Ps, Dn_, Ds_ = counts_at(poly, Dn, Ds, 0.05)
            allrows.append(dict(Orthogroup=og, Pn=Pn, Ps=Ps,
                                Dn=Dn_, Ds=Ds_, n_codons=L))
        pd.DataFrame(allrows).to_csv(
            f"14_MK_Test/02_results/mk_v2_vs_{ogsp}_allgenes.csv", index=False)
        print(f"  جدول کامل به‌ازای ژن: {len(allrows)} ژن")

        # آزمون علامت روی DoS در آستانهٔ ۵٪
        rows = []
        for og, (poly, Dn, Ds, nin, nout, L) in scans.items():
            Pn, Ps, Dn_, Ds_ = counts_at(poly, Dn, Ds, 0.05)
            if (Dn_+Ds_) == 0 or (Pn+Ps) == 0: continue
            dos = Dn_/(Dn_+Ds_) - Pn/(Pn+Ps)
            p = fisher_exact([[round(Dn_), round(Ds_)], [round(Pn), round(Ps)]])[1]
            rows.append(dict(Orthogroup=og, Pn=Pn, Ps=Ps, Dn=Dn_, Ds=Ds_,
                             DoS=dos, fisher_p=p, n_codons=L))
        d = pd.DataFrame(rows).merge(ann[["Orthogroup","product"]],
                                     on="Orthogroup", how="left")
        d = d.sort_values("fisher_p").reset_index(drop=True)
        _bh = d.fisher_p.values * len(d) / np.arange(1, len(d)+1)
        d["q_BH"] = np.clip(np.minimum.accumulate(_bh[::-1])[::-1], 0, 1)
        d.to_csv(f"14_MK_Test/02_results/mk_v2_vs_{ogsp}.csv", index=False)

        pos = int((d.DoS > 0).sum()); tot_n = int((d.DoS != 0).sum())
        bt = binomtest(pos, tot_n, 0.5)
        wp = wilcoxon(d.DoS[d.DoS != 0])[1]
        print(f"\n--- آزمون علامت DoS (آستانهٔ ۵٪، {tot_n} ژن) ---")
        print(f"  DoS>0 در {pos}/{tot_n} ژن | binomial p={bt.pvalue:.4f} "
              f"| Wilcoxon p={wp:.4f} | میانه DoS={d.DoS.median():+.4f}")
        print(f"\n--- ۸ ژن با بالاترین DoS (شمارش کافی) ---")
        big = d[(d.Dn+d.Ds >= 5) & (d.Pn+d.Ps >= 5)].nlargest(8, "DoS")
        print(big[["Orthogroup","Pn","Ps","Dn","Ds","DoS","fisher_p","q_BH","product"]]
              .to_string(index=False, max_colwidth=40))
