#!/usr/bin/env python3
"""
parse_busted.py — پس‌پردازش کامل نتایج BUSTED مرحلهٔ ۱۰

از پوشهٔ 03_busted_results/ فایل‌های JSON را می‌خواند و جدول نهایی مقاله
را می‌سازد. تمام مراحلی که به‌صورت تعاملی اجرا شدند در این فایل تجمیع
شده‌اند تا کل مرحله با یک دستور بازتولید شود.

اجرا (از داخل پوشهٔ 10_Selection_Pressure):
    conda activate orthofinder_env
    python3 scripts/parse_busted.py

خروجی (در 04_summary_tables/):
    busted_summary.csv                    نتایج خام BUSTED + تصحیح BH
    orthogroup_annotations.csv            انوتیشن Prokka هر orthogroup
    busted_annotated.csv                  نتایج + انوتیشن
    busted_annotated_with_divergence.csv  + طول درخت و واگرایی
    busted_final_table.csv                جدول نهایی مقاله
"""
import json, glob, os, sys
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

ROOT     = Path("..")                                  # Pangenome_Analysis/
RESULTS  = ROOT / "03_OrthoFinder_Results/OrthoFinder_Results_v3/Results_Aug04"
JSON_DIR = Path("03_busted_results")
OUT      = Path("04_summary_tables"); OUT.mkdir(exist_ok=True)

# ══════════════════════════════════════════════════════════════════
# ۱) پارس فایل‌های JSON
# ══════════════════════════════════════════════════════════════════
rows, failed = [], []

for path in sorted(glob.glob(str(JSON_DIR / "*_BUSTED.json"))):
    og = os.path.basename(path).replace("_BUSTED.json", "")
    try:
        d = json.load(open(path))
    except Exception as e:
        failed.append((og, f"JSON نامعتبر: {e}"))
        continue

    tr = d["test results"]

    # توزیع ω در مدل غیرمقید (شاخهٔ Test)
    rd = d["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    cls    = [(rd[k]["omega"], rd[k]["proportion"]) for k in sorted(rd)]
    omegas = np.array([c[0] for c in cls]); props = np.array([c[1] for c in cls])

    # ω کلی از مدل MG94 — ساختار [[value, weight], ...]
    mg  = d["fits"]["MG94xREV with separate rates for branch sets"]["Rate Distributions"]
    key = "non-synonymous/synonymous rate ratio for *test*"
    vv  = mg.get(key, [])
    if vv:
        v = np.array([x[0] for x in vv], float); w = np.array([x[1] for x in vv], float)
        omega_mg94 = float((v * w).sum() / w.sum())
    else:
        omega_mg94 = np.nan

    # نرخ‌های مترادف سایت‌به‌سایت (SRV)
    srv   = d["fits"]["Unconstrained model"]["Rate Distributions"].get(
                "Synonymous site-to-site rates", {})
    rates = np.array([srv[k]["rate"] for k in sorted(srv)], float)
    sprop = np.array([srv[k]["proportion"] for k in sorted(srv)], float)
    prop_dS0 = float(sprop[rates < 0.01].sum()) if len(rates) else np.nan

    # مجموع طول شاخه‌ها در هر مدل
    ba = d.get("branch attributes", {}).get("0", {})
    def total_len(model_key):
        return sum(v for br in ba.values()
                   if isinstance(v := br.get(model_key), (int, float)))

    inp = d.get("input", {})
    rows.append({
        "Orthogroup": og,
        "N_sequences": inp.get("number of sequences"),
        "N_codons": inp.get("number of sites"),
        "LRT": tr.get("LRT"),
        "p_value": tr.get("p-value"),
        "omega_mean": float((omegas * props).sum()),
        "omega_max": float(omegas.max()),
        "prop_sites_omega_gt1": float(props[omegas > 1].sum()),
        "omega_MG94": omega_mg94,
        "prop_dS_near_zero": prop_dS0,
        "n_SRV_classes": len(rates),
        "treelen_MG94": total_len("MG94xREV with separate rates for branch sets"),
        "treelen_GTR": total_len("Nucleotide GTR"),
        "n_branches": len(ba),
    })

df = pd.DataFrame(rows)
print(f"فایل‌های پارس‌شده: {len(df)}")
for og, msg in failed:
    print(f"  مشکل‌دار: {og} — {msg}")
if df.empty:
    sys.exit("هیچ نتیجه‌ای پارس نشد")

# ══════════════════════════════════════════════════════════════════
# ۲) تصحیح Benjamini-Hochberg
# ══════════════════════════════════════════════════════════════════
df = df.sort_values("p_value").reset_index(drop=True)
n  = len(df)
bh = df["p_value"].values * n / np.arange(1, n + 1)
df["q_value_BH"] = np.minimum(np.minimum.accumulate(bh[::-1])[::-1], 1.0)
df["Selection"] = np.where(
    df["q_value_BH"] < 0.05, "Episodic diversifying (positive)",
    np.where(df["omega_mean"] < 1, "Purifying", "Neutral / inconclusive"))
df["divergence_per_codon"] = df["treelen_MG94"] / df["N_codons"]

df.drop(columns=["omega_MG94", "prop_dS_near_zero", "n_SRV_classes",
                 "treelen_MG94", "treelen_GTR", "n_branches",
                 "divergence_per_codon"]).to_csv(OUT / "busted_summary.csv", index=False)

# ══════════════════════════════════════════════════════════════════
# ۳) انوتیشن Prokka (رأی اکثریت میان ۲۹۰ ژنوم)
# ══════════════════════════════════════════════════════════════════
tsv_paths = {}
with open(ROOT / "06_Recombination/nucleotide_file_map.tsv") as fh:
    next(fh)
    for line in fh:
        gid, fna = line.rstrip("\n").split("\t")
        tsv_paths[gid] = fna[:-4] + ".tsv"

og_table = pd.read_csv(RESULTS / "Orthogroups/Orthogroups.tsv",
                       sep="\t").set_index("Orthogroup")
genome_cols = og_table.columns.tolist()
sco_list = [l.strip() for l in
            open(ROOT / "04_Pangenome_Statistics/strict_single_copy_orthogroups_v3.txt")]

ann_cache = {}
def get_ann(gid):
    if gid not in ann_cache:
        p, d = tsv_paths.get(gid), {}
        if p and Path(p).exists():
            t = pd.read_csv(p, sep="\t", dtype=str).fillna("")
            for _, r in t.iterrows():
                d[r["locus_tag"]] = (r.get("gene", ""), r.get("product", ""), r.get("COG", ""))
        ann_cache[gid] = d
    return ann_cache[gid]

ann_rows = []
for og in sco_list:
    if og not in og_table.index:
        continue
    r = og_table.loc[og]
    prods, genes, cogs, found = Counter(), Counter(), Counter(), 0
    for g in genome_cols:
        if pd.isna(r[g]):
            continue
        a = get_ann(g).get(str(r[g]).strip())
        if a is None:
            continue
        found += 1
        gene, product, cog = a
        if product: prods[product] += 1
        if gene:    genes[gene] += 1
        if cog:     cogs[cog] += 1
    top = prods.most_common(1)[0] if prods else ("", 0)
    ann_rows.append({
        "Orthogroup": og, "n_annotated": found, "product": top[0],
        "product_agreement_pct": round(100 * top[1] / found, 1) if found else 0,
        "n_distinct_products": len(prods),
        "gene": genes.most_common(1)[0][0] if genes else "",
        "COG": cogs.most_common(1)[0][0] if cogs else "",
    })

ann_df = pd.DataFrame(ann_rows)
ann_df.to_csv(OUT / "orthogroup_annotations.csv", index=False)
print(f"انوتیشن: {len(ann_df)} orthogroup | "
      f"hypothetical: {ann_df['product'].str.contains('hypothetical', case=False).sum()}")

# ══════════════════════════════════════════════════════════════════
# ۴) ادغام و جدول نهایی
# ══════════════════════════════════════════════════════════════════
m = df.merge(ann_df, on="Orthogroup", how="left").sort_values("q_value_BH")
m.to_csv(OUT / "busted_annotated.csv", index=False)
m.to_csv(OUT / "busted_annotated_with_divergence.csv", index=False)
m.to_csv(OUT / "busted_final_table.csv", index=False)

# ══════════════════════════════════════════════════════════════════
# ۵) گزارش
# ══════════════════════════════════════════════════════════════════
pd.set_option("display.max_colwidth", 45)
print("\n=== ژن‌های معنادار (q < 0.05)، مرتب بر اساس omega_MG94 ===")
print(m[m["q_value_BH"] < 0.05].sort_values("omega_MG94", ascending=False)[
    ["Orthogroup","N_codons","omega_max","prop_sites_omega_gt1",
     "omega_MG94","q_value_BH","product"]].to_string(index=False))

print("\n=== ۱۵ ژن برتر بر اساس omega_MG94 ===")
top = m.nlargest(15, "omega_MG94").copy()
top["sig"] = top["q_value_BH"] < 0.05
print(top[["Orthogroup","omega_MG94","q_value_BH","sig","product"]].to_string(index=False))

print("\n=== توزیع و همبستگی‌ها ===")
print(m["Selection"].value_counts().to_string())
s  = m.dropna(subset=["omega_MG94"])
lg = np.log10(s["omega_max"].clip(lower=1e-9))
for c in ["omega_MG94", "divergence_per_codon", "prop_dS_near_zero"]:
    print(f"  Spearman(log10 omega_max , {c}) = {lg.corr(s[c], method='spearman'):.3f}")

print(f"\nذخیره شد در {OUT}/")
