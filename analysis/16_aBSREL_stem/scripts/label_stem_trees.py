#!/usr/bin/env python3
"""مرحلهٔ ۱۶ — برچسب {Stem} روی یک یال مشخص.

پیش‌فرض (بدون آرگومان) یالِ منتهی به کلاد SPPV+GTPV را برچسب می‌زند. در
درخت بدون ریشه این همان یالی است که LSDV را از برون‌گروه جدا می‌کند —
یعنی دقیقاً یالی که MK قطبی‌شده واگرایی را به آن نسبت می‌دهد.

با آرگومان SPPV یا GTPV، کنترل اختصاصیت ساخته می‌شود (یالِ منتهی به آن
کلاد) و خروجی در 04_specificity ذخیره می‌شود.

    python3 16_aBSREL_stem/scripts/label_stem_trees.py            # اصلی
    python3 16_aBSREL_stem/scripts/label_stem_trees.py SPPV       # کنترل

هر فایل خروجی باید دقیقاً یک {Stem} داشته باشد.
اجرا از ریشهٔ پروژه با orthofinder_env."""
from Bio import Phylo
import pandas as pd, glob, os, sys

target = sys.argv[1] if len(sys.argv) > 1 else "OUTGROUP"
want = {"SPPV", "GTPV"} if target == "OUTGROUP" else {target}
outdir = ("16_aBSREL_stem/01_labeled_trees" if target == "OUTGROUP"
          else "16_aBSREL_stem/04_specificity")
suffix = "_stem" if target == "OUTGROUP" else f"_{target}stem"
os.makedirs(outdir, exist_ok=True)

sp = pd.read_csv("00_Metadata/species_assignment.csv")
smap = dict(zip(sp.Name, sp.species))

ok, skipped = 0, []
for f in sorted(glob.glob("10_Selection_Pressure/02_pruned_trees/*_pruned.nwk")):
    og = os.path.basename(f).replace("_pruned.nwk", "")
    tree = Phylo.read(f, "newick")
    tgt = next((cl for cl in tree.get_nonterminals()
                if {smap.get(t.name) for t in cl.get_terminals()} == want), None)
    if tgt is None:
        skipped.append(og); continue
    for cl in tree.get_nonterminals():
        cl.confidence = None
        cl.name = "{Stem}" if cl is tgt else None
    Phylo.write(tree, f"{outdir}/{og}{suffix}.nwk", "newick")
    ok += 1

print(f"✓ {ok} درخت برچسب خورد (هدف: {target})")
if skipped:
    print(f"⚠ کلاد پیدا نشد در: {skipped}")
