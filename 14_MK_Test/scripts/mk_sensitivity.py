#!/usr/bin/env python3
"""حساسیت MK: حذف سویه‌های واکسن و حذف ژنوم‌های SRA."""
import sys, glob, os
sys.path.insert(0, "14_MK_Test/scripts")
from mk_test_v2 import scan_gene, counts_at, cmh
import pandas as pd, numpy as np

ALN = "10_Selection_Pressure/01_codon_alignments"
sp  = pd.read_csv("00_Metadata/species_assignment.csv")
lsdv = sp.loc[sp.species == "LSDV", "Name"].tolist()
grp = {s: sp.loc[sp.species == s, "Name"].tolist() for s in ["SPPV", "GTPV"]}

tr = pd.read_csv("11_Scoary/01_inputs/traits.csv", index_col=0)
print("ستون‌های traits.csv:", tr.columns.tolist())
vcol = [c for c in tr.columns if "vacc" in c.lower()]
vacc = set(tr.index[tr[vcol[0]] == 1]) if vcol else set()
print(f"ستون واکسن: {vcol} → {len(vacc)} سویه\n")

sets = {
    "همه LSDV":        lsdv,
    "بدون واکسن":      [g for g in lsdv if g not in vacc],
    "بدون SRA":        [g for g in lsdv if not str(g).startswith("SRR")],
    "بدون هر دو":      [g for g in lsdv if g not in vacc and not str(g).startswith("SRR")],
}

print(f"{'زیرمجموعه':<14}{'n':>5}{'outgrp':>8}{'cut':>6}"
      f"{'Pn':>8}{'Ps':>8}{'alpha':>8}{'CMH_p':>11}{'MH_OR':>7}")
for label, ing in sets.items():
    for ogsp in ["SPPV", "GTPV"]:
        scans = {}
        for f in sorted(glob.glob(f"{ALN}/*_codon.fasta")):
            r = scan_gene(f, ing, grp[ogsp])
            if r: scans[os.path.basename(f)] = r
        for cut in [0.05, 0.10]:
            tot = np.zeros(4); tabs = []
            for poly, Dn, Ds, *_ in scans.values():
                Pn, Ps, Dn_, Ds_ = counts_at(poly, Dn, Ds, cut)
                tot += [Pn, Ps, Dn_, Ds_]
                tabs.append([[round(Dn_), round(Ds_)], [round(Pn), round(Ps)]])
            Pn, Ps, Dn, Ds = tot
            NI = (Pn/Ps)/(Dn/Ds) if Ps and Ds and Dn else np.nan
            p, orr, _ = cmh(tabs)
            print(f"{label:<14}{len(ing):>5}{ogsp:>8}{cut:>6.2f}"
                  f"{Pn:>8.1f}{Ps:>8.1f}{1-NI:>8.3f}{p:>11.2e}{orr:>7.3f}")
