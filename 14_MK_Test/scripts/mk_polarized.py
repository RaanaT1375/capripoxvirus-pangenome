#!/usr/bin/env python3
"""MK قطبی‌شده: حالت اجدادی = کدونی که SPPV و GTPV بر آن توافق دارند.
تفاوت‌های ثابت به شاخهٔ LSDV نسبت داده می‌شوند. درون‌گروه: LSDV بدون واکسن."""
import sys, glob, os, collections
sys.path.insert(0, "14_MK_Test/scripts")
from mk_test_v2 import tr, path_counts, clean, cmh
from Bio import SeqIO
from scipy.stats import fisher_exact, binomtest, wilcoxon
import pandas as pd, numpy as np

ALN, CUT, CONS, COV = "10_Selection_Pressure/01_codon_alignments", 0.05, 0.90, 0.90

def consensus(seqs, sl):
    c = [x for x in (s[sl] for s in seqs) if clean(x)]
    if not c: return None
    cod, n = collections.Counter(c).most_common(1)[0]
    return cod if n/len(c) >= CONS else None

def polarized_gene(path, ing_ids, sppv, gtpv):
    r = {x.id: str(x.seq).upper() for x in SeqIO.parse(path, "fasta")}
    ing = [r[i] for i in ing_ids if i in r]
    s1  = [r[i] for i in sppv    if i in r]
    s2  = [r[i] for i in gtpv    if i in r]
    if not ing or not s1 or not s2: return None
    L = len(next(iter(r.values()))) // 3
    Pn = Ps = Dn = Ds = 0.0; n_anc = 0
    for k in range(L):
        sl = slice(3*k, 3*k+3)
        a1, a2 = consensus(s1, sl), consensus(s2, sl)
        if a1 is None or a2 is None or a1 != a2: continue   # حالت اجدادی نامطمئن
        n_anc += 1
        ic = [x for x in (s[sl] for s in ing) if clean(x)]
        if len(ic) < COV*len(ing): continue
        icc = collections.Counter(ic); major = icc.most_common(1)[0][0]
        if len(icc) > 1:
            for v, cnt in icc.items():
                if v == major or cnt/len(ic) < CUT: continue
                s_, n_ = path_counts(major, v); Ps += s_; Pn += n_
        elif major != a1:
            s_, n_ = path_counts(a1, major); Ds += s_; Dn += n_
    return dict(Pn=Pn, Ps=Ps, Dn=Dn, Ds=Ds, n_anc_sites=n_anc, n_codons=L)

sp = pd.read_csv("00_Metadata/species_assignment.csv")
grp = {s: sp.loc[sp.species == s, "Name"].tolist() for s in ["LSDV","SPPV","GTPV"]}
vacc = set(pd.read_csv("11_Scoary/01_inputs/traits.csv", index_col=0)
             .query("Vaccine == 1").index)
ing_ids = [g for g in grp["LSDV"] if g not in vacc]
ann = pd.read_csv("10_Selection_Pressure/04_summary_tables/orthogroup_annotations.csv")
print(f"درون‌گروه: {len(ing_ids)} ژنوم LSDV (بدون {len(vacc)} واکسن)\n")

rows = []
for f in sorted(glob.glob(f"{ALN}/*_codon.fasta")):
    d = polarized_gene(f, ing_ids, grp["SPPV"], grp["GTPV"])
    if d: d["Orthogroup"] = os.path.basename(f).replace("_codon.fasta",""); rows.append(d)

df = pd.DataFrame(rows)
Pn, Ps, Dn, Ds = df[["Pn","Ps","Dn","Ds"]].sum()
NI = (Pn/Ps)/(Dn/Ds)
tabs = [[[round(r.Dn), round(r.Ds)], [round(r.Pn), round(r.Ps)]] for r in df.itertuples()]
p_cmh, mhor, _ = cmh(tabs)
print(f"=== MK قطبی‌شده (شاخهٔ LSDV، cutoff {CUT}) ===")
print(f"  سایت‌های با حالت اجدادی مطمئن: {df.n_anc_sites.sum():.0f} "
      f"از {df.n_codons.sum():.0f} کدون ({100*df.n_anc_sites.sum()/df.n_codons.sum():.1f}%)")
print(f"  Pn={Pn:.1f}  Ps={Ps:.1f}  Dn={Dn:.1f}  Ds={Ds:.1f}")
print(f"  NI={NI:.3f} | alpha={1-NI:.3f} | CMH p={p_cmh:.3e} | MH_OR={mhor:.3f}")

df["DoS"] = df.Dn/(df.Dn+df.Ds) - df.Pn/(df.Pn+df.Ps)
df["fisher_p"] = [fisher_exact([[round(r.Dn),round(r.Ds)],[round(r.Pn),round(r.Ps)]])[1]
                  for r in df.itertuples()]
d = df.dropna(subset=["DoS"]).sort_values("fisher_p").reset_index(drop=True)
_bh = d.fisher_p.values*len(d)/np.arange(1, len(d)+1)
d["q_BH"] = np.clip(np.minimum.accumulate(_bh[::-1])[::-1], 0, 1)
d = d.merge(ann[["Orthogroup","product"]], on="Orthogroup", how="left")
d.to_csv("14_MK_Test/02_results/mk_polarized_LSDV_branch.csv", index=False)

nz = d[d.DoS.notna() & (d.DoS != 0)]
print(f"\n  DoS>0 در {int((nz.DoS>0).sum())}/{len(nz)} ژن | "
      f"binomial p={binomtest(int((nz.DoS>0).sum()), len(nz), 0.5).pvalue:.4f} | "
      f"Wilcoxon p={wilcoxon(nz.DoS)[1]:.4f} | میانه={nz.DoS.median():+.4f}")
print(f"\n--- ۱۰ ژن برتر بر اساس DoS (شمارش ≥۵ در هر ردیف) ---")
big = d[(d.Dn+d.Ds >= 5) & (d.Pn+d.Ps >= 5)].nlargest(10, "DoS")
print(big[["Orthogroup","Pn","Ps","Dn","Ds","DoS","fisher_p","q_BH","product"]]
      .to_string(index=False, max_colwidth=40))
