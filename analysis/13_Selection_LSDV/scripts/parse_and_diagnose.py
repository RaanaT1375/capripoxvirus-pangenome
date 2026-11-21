#!/usr/bin/env python3
"""مرحلهٔ ۱۳ — پارس BUSTED (LSDV-only) + تشخیص dS.
اعدادی که تولید می‌کند و در گزارش استفاده شده‌اند:
  - ۱۱ ژن q<0.05 که همه از سایت‌های dS≈0 می‌آیند
  - prop_dS_near0 بین ۰.۸۸ و ۰.۹۶ برای هر ۱۱ ژن
  - LRT=0 در ۳۰ از ۵۶ ژن قابل آزمون
اجرا از ریشهٔ پروژه با orthofinder_env."""
import json, glob, os
import pandas as pd, numpy as np

rows = []
for path in sorted(glob.glob("13_Selection_LSDV/03_busted_results/*_BUSTED.json")):
    og = os.path.basename(path).replace("_BUSTED.json", "")
    d = json.load(open(path))
    tr = d["test results"]
    ud = d["fits"]["Unconstrained model"]["Rate Distributions"]
    rd = ud["Test"]
    om = np.array([rd[k]["omega"] for k in sorted(rd)])
    pr = np.array([rd[k]["proportion"] for k in sorted(rd)])

    mg = d["fits"]["MG94xREV with separate rates for branch sets"]["Rate Distributions"]
    vv = mg.get("non-synonymous/synonymous rate ratio for *test*", [])
    if vv:
        v = np.array([x[0] for x in vv], float); w = np.array([x[1] for x in vv], float)
        omg = float((v * w).sum() / w.sum())
    else:
        omg = np.nan

    srv = ud.get("Synonymous site-to-site rates")
    if srv:
        r = np.array([srv[k]["rate"] for k in sorted(srv)])
        p = np.array([srv[k]["proportion"] for k in sorted(srv)])
        near0, dsmin, dsmax = float(p[r < 0.05].sum()), float(r.min()), float(r.max())
    else:
        near0 = dsmin = dsmax = np.nan

    inp = d.get("input", {})
    rows.append(dict(Orthogroup=og, N_sequences=inp.get("number of sequences"),
                     N_codons=inp.get("number of sites"),
                     LRT=tr.get("LRT"), p_value=tr.get("p-value"),
                     omega_max=float(om.max()), omega_MG94=omg,
                     prop_sites_omega_gt1=float(pr[om > 1].sum()),
                     prop_dS_near0=near0, dS_min=dsmin, dS_max=dsmax))

df = pd.DataFrame(rows)
div = pd.read_csv("13_Selection_LSDV/diversity_LSDV_only.csv")
df = df.merge(div[["Orthogroup", "uniq_aa", "poly_nt_sites", "poly_aa_sites"]], on="Orthogroup")

testable = df[df.poly_nt_sites > 0].sort_values("p_value").reset_index(drop=True)
excluded = df[df.poly_nt_sites == 0]
n = len(testable)
bh = testable.p_value.values * n / np.arange(1, n + 1)
testable["q_value_BH"] = np.clip(np.minimum.accumulate(bh[::-1])[::-1], 0, 1)

ann = pd.read_csv("10_Selection_Pressure/04_summary_tables/orthogroup_annotations.csv")
testable = testable.merge(ann[["Orthogroup", "product"]], on="Orthogroup", how="left")
testable.to_csv("13_Selection_LSDV/04_summary_tables/busted_lsdv_final.csv", index=False)

print(f"قابل آزمون: {n} | کنار گذاشته (صفر سایت چندشکلی): {len(excluded)} "
      f"({', '.join(excluded.Orthogroup)})")
print(f"معنادار q<0.05: {(testable.q_value_BH < 0.05).sum()}")
print(f"LRT = 0 در {(testable.LRT == 0).sum()} از {n} ژن")
print(f"میانه omega_MG94: {testable.omega_MG94.median():.4f}")
print(f"میانه uniq_aa: {df.uniq_aa.median():.0f} از ۲۴۰ | ژن با uniq_aa<5: {(df.uniq_aa<5).sum()}")
sig = testable[testable.q_value_BH < 0.05]
print("\nprop_dS_near0 در ژن‌های معنادار:")
print(sig[["Orthogroup","prop_dS_near0","dS_max","poly_aa_sites","omega_MG94","q_value_BH"]]
      .to_string(index=False))
