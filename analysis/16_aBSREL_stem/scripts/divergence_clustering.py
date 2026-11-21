#!/usr/bin/env python3
"""مرحلهٔ ۱۶ — آیا جایگزینی‌های ثابت روی یال LSDV در طول ژن پراکنده‌اند
یا در یک بلوک جمع شده‌اند؟

این همان تحلیلی است که ادعای انتخاب روی OG0000120 را باطل کرد. منطقش:
انتخاب متنوع‌کننده جایگزینی‌ها را پراکنده می‌کند، ورود یک قطعهٔ واگرا
(نوترکیبی، یا ناحیهٔ ایندل‌دار غیرقابل‌الاینمنت) آن‌ها را متمرکز می‌کند.

اعداد گزارش‌شده:
  فقط ۹ ژن پنج جایگزینی ثابت یا بیشتر دارند
  تنها OG0000120 پس از FDR خوشه‌ای است: KS D=0.802، q=3.0e-09،
  هر ۱۴ جایگزینی در کدون ۱۷۴ تا ۲۰۲ از ۲۱۷
  میانه KS_D = 0.396 و میانه span = 0.735 → خوشه‌ای بودن عمومی نیست

دو حالت اجرا:
    python3 .../divergence_clustering.py               # همهٔ ۶۰ ژن
    python3 .../divergence_clustering.py --gaps OG0000120   # بازرسی یک ژن
حالت دوم هیستوگرام موقعیت‌ها و درصد گپ را در پنجره‌های ۲۰ کدونی می‌دهد.
اجرا از ریشهٔ پروژه با orthofinder_env."""
import sys, collections, glob, os
sys.path.insert(0, "14_MK_Test/scripts")
from mk_test_v2 import tr, clean
from Bio import SeqIO
from scipy.stats import kstest
import numpy as np, pandas as pd

ALN = "10_Selection_Pressure/01_codon_alignments"
sp = pd.read_csv("00_Metadata/species_assignment.csv")
grp = {s: sp.loc[sp.species == s, "Name"].tolist() for s in ["LSDV", "SPPV", "GTPV"]}
vacc = set(pd.read_csv("11_Scoary/01_inputs/traits.csv",
                       index_col=0).query("Vaccine==1").index)
ing_ids = [g for g in grp["LSDV"] if g not in vacc]

def cons(seqs, sl, thr=0.90):
    c = [x for x in (s[sl] for s in seqs) if clean(x)]
    if not c:
        return None
    cod, n = collections.Counter(c).most_common(1)[0]
    return cod if n / len(c) >= thr else None

def fixed_positions(f):
    """موقعیت کدون‌هایی که LSDV تک‌شکل است و با حالت اجدادی مشترک
    SPPV/GTPV فرق دارد، تفکیک‌شده به مترادف و غیرمترادف."""
    r = {x.id: str(x.seq).upper() for x in SeqIO.parse(f, "fasta")}
    ing = [r[i] for i in ing_ids if i in r]
    s1 = [r[i] for i in grp["SPPV"] if i in r]
    s2 = [r[i] for i in grp["GTPV"] if i in r]
    L = len(next(iter(r.values()))) // 3
    syn, non = [], []
    for k in range(L):
        sl = slice(3 * k, 3 * k + 3)
        a1, a2 = cons(s1, sl), cons(s2, sl)
        if a1 is None or a2 is None or a1 != a2:
            continue
        ic = [x for x in (s[sl] for s in ing) if clean(x)]
        if len(ic) < 0.9 * len(ing):
            continue
        icc = collections.Counter(ic)
        major = icc.most_common(1)[0][0]
        if len(icc) > 1 or major == a1:
            continue
        (non if tr(major) != tr(a1) else syn).append(k)
    return r, L, syn, non

if "--gaps" in sys.argv:
    og = sys.argv[sys.argv.index("--gaps") + 1]
    r, L, syn, non = fixed_positions(f"{ALN}/{og}_codon.fasta")
    pos = sorted(syn + non)
    print(f"{og}: {len(non)} غیرمترادف + {len(syn)} مترادف")
    if pos:
        print(f"محدوده: کدون {pos[0]}–{pos[-1]} از {L}")
        print(f"هیستوگرام ۱۰ پنجره‌ای: "
              f"{list(np.histogram(pos, bins=10, range=(0, L))[0])}")
    print("درصد گپ در پنجره‌های ۲۰ کدونی:")
    for st in range(0, L, 20):
        en = min(st + 20, L)
        g = sum(s[3 * st:3 * en].count("-") for s in r.values())
        print(f"  کدون {st:3d}-{en:3d}: {100 * g / (len(r) * (en - st) * 3):5.1f}%")
    sys.exit()

rows = []
for f in sorted(glob.glob(f"{ALN}/*_codon.fasta")):
    og = os.path.basename(f).replace("_codon.fasta", "")
    _, L, syn, non = fixed_positions(f)
    pos = sorted(syn + non)
    if len(pos) < 5:
        continue
    ks = kstest(np.array(pos) / L, "uniform")
    rows.append(dict(Orthogroup=og, n_fixed=len(pos), n_nonsyn=len(non),
                     n_codons=L, KS_D=ks.statistic, KS_p=ks.pvalue,
                     span_frac=(pos[-1] - pos[0] + 1) / L))

d = pd.DataFrame(rows).sort_values("KS_p").reset_index(drop=True)
n = len(d)
d["q_BH"] = np.clip(np.minimum.accumulate(
    (d.KS_p.values * n / np.arange(1, n + 1))[::-1])[::-1], 0, 1)
ann = pd.read_csv("10_Selection_Pressure/04_summary_tables/orthogroup_annotations.csv")
d = d.merge(ann[["Orthogroup", "product"]], on="Orthogroup", how="left")
d.to_csv("16_aBSREL_stem/03_summary/divergence_clustering.csv", index=False)

print(f"ژن با ۵ جایگزینی ثابت یا بیشتر: {n} از ۶۰")
print(f"خوشه‌ای معنادار (q<0.05): {(d.q_BH < 0.05).sum()}")
print(f"میانه KS_D {d.KS_D.median():.3f} | میانه span {d.span_frac.median():.3f}")
print(d[["Orthogroup", "n_fixed", "n_codons", "KS_D", "KS_p", "q_BH",
         "span_frac", "product"]].to_string(index=False, max_colwidth=34))
