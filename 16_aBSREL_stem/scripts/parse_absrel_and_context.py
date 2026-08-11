#!/usr/bin/env python3
"""مرحلهٔ ۱۶ — پارس aBSREL + بستر طول شاخه زیر مدل نوکلئوتیدی GTR.

اعداد گزارش‌شده که این اسکریپت تولید می‌کند:
  ۴۹ ژن قابل آزمون، ۲ معنادار (هر دو رد شده)
  ۱۰ ژن با tested=0 چون روی این یال Dn=Ds=0 دارند
  میانه omega پایه روی یال = 0.146
  سهم یال LSDV از کل طول درخت: میانه 0.083، و OG0000120 با 0.437 رتبهٔ ۱

نکته: کلید Global MG94xREV در JSONهای این نسخهٔ HyPhy وجود ندارد؛ بستر
طول شاخه از Nucleotide GTR گرفته می‌شود که مدل‌مستقل‌تر هم هست.
اجرا از ریشهٔ پروژه با orthofinder_env."""
import json, glob, os
import pandas as pd, numpy as np

def branch_info(path):
    d = json.load(open(path))
    ba = d["branch attributes"]["0"]
    tgt = next((n for n, a in ba.items()
                if a.get("Corrected P-value") is not None), None)
    if tgt is None:
        return None
    gtr = {n: a["Nucleotide GTR"] for n, a in ba.items()
           if a.get("Nucleotide GTR") is not None}
    if tgt not in gtr:
        return None
    v = np.array(list(gtr.values()))
    a = ba[tgt]
    rd = a.get("Rate Distributions") or []
    pos = [(w, wt) for w, wt in rd if w > 1]
    return dict(branch=tgt, gtr=gtr[tgt], frac=gtr[tgt] / v.sum(),
                rank=int((v > gtr[tgt]).sum()) + 1, n_branches=len(v),
                omega_base=a.get("Baseline MG94xREV omega ratio"),
                omega_pos=max([w for w, _ in pos], default=np.nan),
                weight_pos=sum(wt for _, wt in pos) if pos else 0.0,
                LRT=a.get("LRT"), p_corr=a.get("Corrected P-value"),
                n_codons=d["input"]["number of sites"])

rows, untested = [], []
for p in sorted(glob.glob("16_aBSREL_stem/02_results/*_ABSREL.json")):
    og = os.path.basename(p).replace("_ABSREL.json", "")
    try:
        i = branch_info(p)
    except Exception as e:
        print(f"✗ {og}: {str(e)[:40]}"); continue
    if i is None:
        untested.append(og); continue
    i["Orthogroup"] = og
    rows.append(i)

d = pd.DataFrame(rows)
ann = pd.read_csv("10_Selection_Pressure/04_summary_tables/orthogroup_annotations.csv")
d = d.merge(ann[["Orthogroup", "product"]], on="Orthogroup", how="left")
d.to_csv("16_aBSREL_stem/03_summary/absrel_stem_final.csv", index=False)
d.sort_values("frac", ascending=False).to_csv(
    "16_aBSREL_stem/03_summary/stem_branch_gtr.csv", index=False)

print(f"قابل آزمون: {len(d)} | آزمون‌نشده (Dn=Ds=0): {len(untested)}")
print(f"  {untested}")
print(f"معنادار (p_corr<0.05): {(d.p_corr < 0.05).sum()}")
print(f"میانه omega پایه روی یال: {d.omega_base.median():.3f}")
print(f"میانه سهم یال از کل درخت: {d.frac.median():.3f} | "
      f"چارک سوم {d.frac.quantile(.75):.3f}")
print(f"Spearman(p_corr، طول ژن): "
      f"{d[['p_corr','n_codons']].corr(method='spearman').iloc[0,1]:.3f}")
print("\nژن‌های معنادار:")
print(d[d.p_corr < 0.05][["Orthogroup", "omega_base", "omega_pos",
                          "weight_pos", "p_corr", "product"]]
      .to_string(index=False, max_colwidth=38))
print("\n۵ ژن با بیشترین سهم یال LSDV:")
print(d.nlargest(5, "frac")[["Orthogroup", "gtr", "frac", "rank",
                             "omega_base", "product"]]
      .to_string(index=False, max_colwidth=34))
