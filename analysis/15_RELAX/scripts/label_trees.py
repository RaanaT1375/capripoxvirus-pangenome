#!/usr/bin/env python3
"""مرحلهٔ ۱۵ — برچسب‌گذاری {Test}=کلاد LSDV و {Reference}=SPPV+GTPV.
ورودی : 10_Selection_Pressure/02_pruned_trees/OG*_pruned.nwk
خروجی : 15_RELAX/01_labeled_trees/OG*_labeled.nwk  (۶۰ فایل)
بازبینی: در OG0000120 باید {Test}=472 و {Reference}=99 باشد.
اجرا از ریشهٔ پروژه با orthofinder_env."""
from Bio import Phylo
import pandas as pd, glob, os

sp = pd.read_csv("00_Metadata/species_assignment.csv")
smap = dict(zip(sp.Name, sp.species))
lab = lambda s: "Test" if s == "LSDV" else "Reference"

n_ok = 0
for f in sorted(glob.glob("10_Selection_Pressure/02_pruned_trees/*_pruned.nwk")):
    og = os.path.basename(f).replace("_pruned.nwk", "")
    tree = Phylo.read(f, "newick")

    # پاس ۱: ترکیب گونه‌ای گره‌های داخلی را با نام‌های اصلی حساب کن
    internal = {}
    for cl in tree.get_nonterminals():
        s = {smap.get(t.name) for t in cl.get_terminals()}
        internal[id(cl)] = ("Test" if s == {"LSDV"}
                            else ("Reference" if not (s & {"LSDV"}) else None))

    # پاس ۲: برچسب‌گذاری. confidence پاک می‌شود تا bootstrap نوشته نشود
    for cl in tree.get_terminals():
        cl.name = f"{cl.name}{{{lab(smap.get(cl.name))}}}"
    for cl in tree.get_nonterminals():
        cl.confidence = None
        cl.name = f"{{{internal[id(cl)]}}}" if internal[id(cl)] else None

    Phylo.write(tree, f"15_RELAX/01_labeled_trees/{og}_labeled.nwk", "newick")
    n_ok += 1

print(f"✓ {n_ok} درخت برچسب‌گذاری شد")
s = open("15_RELAX/01_labeled_trees/OG0000120_labeled.nwk").read()
print(f"بازبینی OG0000120: Test={s.count('{Test}')} Reference={s.count('{Reference}')}")
