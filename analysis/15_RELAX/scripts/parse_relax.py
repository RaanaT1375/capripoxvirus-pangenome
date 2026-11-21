#!/usr/bin/env python3
"""مرحلهٔ ۱۵ — پارس RELAX + کنترل کیفیت برازش.
اعداد گزارش‌شده که این اسکریپت تولید می‌کند:
  ۳۶ ژن پارس‌شده از ۵۹ اجرا (۲۴ شکست، ۴۱٪)
  ۱۶ ژن p<0.05 با تقسیم ۱۱ تشدید / ۵ سست‌شدن؛ میانه K=1.242
  omega_max>100 در ۱۱ ژن؛ Fisher p=1.00؛ Mann-Whitney p=0.043
  با فیلتر سخت (0.01<K<10 و omega_max<10) چهار ژن باقی می‌ماند
اجرا از ریشهٔ پروژه با orthofinder_env."""
import json, glob, os
from scipy.stats import fisher_exact, mannwhitneyu
import pandas as pd, numpy as np

rows = []
for p in sorted(glob.glob("15_RELAX/02_results/*_RELAX.json")):
    og = os.path.basename(p).replace("_RELAX.json", "")
    try:
        d = json.load(open(p))
    except Exception:
        continue
    tr = d.get("test results", {})
    fits = d["fits"]["RELAX alternative"]["Rate Distributions"]
    t, r = fits.get("Test", {}), fits.get("Reference", {})
    tv = tuple(round(t[k]["omega"], 6) for k in sorted(t))
    rv = tuple(round(r[k]["omega"], 6) for k in sorted(r))
    rows.append(dict(Orthogroup=og,
                     K=tr.get("relaxation or intensification parameter"),
                     p=tr.get("p-value"), LRT=tr.get("LRT"),
                     identical_dists=(tv == rv),
                     omega_max=max(max(tv), max(rv))))

d = pd.DataFrame(rows)
d["degen"] = d.omega_max > 100
d["sig"] = d.p < 0.05
d["direction"] = np.where(d.K > 1, "intensification", "relaxation")
d.to_csv("15_RELAX/03_summary/relax_summary.csv", index=False)

print(f"پارس‌شده: {len(d)} | p<0.05: {d.sig.sum()} "
      f"({d[d.sig].direction.value_counts().to_dict()})")
print(f"میانه K: {d.K.median():.3f} | omega_max>100 در {d.degen.sum()} ژن")
print(f"توزیع test و reference یکسان در {d.identical_dists.sum()} ژن")

ct = pd.crosstab(d.degen, d.sig)
if ct.shape == (2, 2):
    print(f"معناداری در برابر دژنراسیون: Fisher p={fisher_exact(ct.values)[1]:.4f}")
    print(f"K: دژنره {d[d.degen].K.median():.3f} در برابر بقیه "
          f"{d[~d.degen].K.median():.3f} | "
          f"Mann-Whitney p={mannwhitneyu(d[d.degen].K, d[~d.degen].K)[1]:.4f}")

ok = d[d.sig & d.K.between(0.01, 10) & (d.omega_max < 10)].sort_values("p")
ann = pd.read_csv("10_Selection_Pressure/04_summary_tables/orthogroup_annotations.csv")
ok = ok.merge(ann[["Orthogroup", "product"]], on="Orthogroup", how="left")
print(f"\nقابل تفسیر (0.01<K<10، omega_max<10): {len(ok)} ژن")
print(ok[["Orthogroup", "K", "p", "omega_max", "product"]]
      .to_string(index=False, max_colwidth=40))
