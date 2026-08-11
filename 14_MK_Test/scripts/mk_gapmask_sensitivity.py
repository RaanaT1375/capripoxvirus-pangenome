#!/usr/bin/env python3
"""مرحلهٔ ۱۴ — پنجمین آزمون استحکام MK: ماسک ستون‌های گپ‌دار.
نتیجه‌ای که تولید می‌کند: alpha از 0.604 (بدون ماسک) به 0.577 (گپ>۵٪)
و 0.570 (گپ>۱٪) می‌رسد؛ Dn فقط از 74.8 به 69.8 افت می‌کند."""
import sys, glob, collections
sys.path.insert(0, "14_MK_Test/scripts")
from mk_test_v2 import path_counts, clean, cmh
from Bio import SeqIO
import pandas as pd, numpy as np

ALN = "10_Selection_Pressure/01_codon_alignments"
sp = pd.read_csv("00_Metadata/species_assignment.csv")
grp = {s: sp.loc[sp.species == s, "Name"].tolist() for s in ["LSDV", "SPPV", "GTPV"]}
vacc = set(pd.read_csv("11_Scoary/01_inputs/traits.csv", index_col=0).query("Vaccine==1").index)
ing_ids = [g for g in grp["LSDV"] if g not in vacc]

def cons(seqs, sl, thr=0.90):
    c = [x for x in (s[sl] for s in seqs) if clean(x)]
    if not c: return None
    cod, n = collections.Counter(c).most_common(1)[0]
    return cod if n / len(c) >= thr else None

for MAXGAP in [1.00, 0.05, 0.01]:
    tot = np.zeros(4); tabs = []; nmask = ntot = 0
    for f in sorted(glob.glob(f"{ALN}/*_codon.fasta")):
        r = {x.id: str(x.seq).upper() for x in SeqIO.parse(f, "fasta")}
        ing = [r[i] for i in ing_ids if i in r]
        s1 = [r[i] for i in grp["SPPV"] if i in r]
        s2 = [r[i] for i in grp["GTPV"] if i in r]
        L = len(next(iter(r.values()))) // 3
        Pn = Ps = Dn = Ds = 0.0
        for k in range(L):
            ntot += 1
            sl = slice(3*k, 3*k+3)
            if sum(1 for s in r.values() if "-" in s[sl]) / len(r) > MAXGAP:
                nmask += 1; continue
            a1, a2 = cons(s1, sl), cons(s2, sl)
            if a1 is None or a2 is None or a1 != a2: continue
            ic = [x for x in (s[sl] for s in ing) if clean(x)]
            if len(ic) < 0.9 * len(ing): continue
            icc = collections.Counter(ic); major = icc.most_common(1)[0][0]
            if len(icc) > 1:
                for v, cnt in icc.items():
                    if v == major or cnt / len(ic) < 0.05: continue
                    s_, n_ = path_counts(major, v); Ps += s_; Pn += n_
            elif major != a1:
                s_, n_ = path_counts(a1, major); Ds += s_; Dn += n_
        tot += [Pn, Ps, Dn, Ds]
        tabs.append([[round(Dn), round(Ds)], [round(Pn), round(Ps)]])
    Pn, Ps, Dn, Ds = tot
    NI = (Pn/Ps) / (Dn/Ds)
    p, orr, _ = cmh(tabs)
    print(f"MAXGAP={MAXGAP:.2f} | ماسک {nmask}/{ntot} | Pn={Pn:.1f} Ps={Ps:.1f} "
          f"Dn={Dn:.1f} Ds={Ds:.1f} | alpha={1-NI:.3f} | CMH p={p:.2e} | OR={orr:.3f}")
