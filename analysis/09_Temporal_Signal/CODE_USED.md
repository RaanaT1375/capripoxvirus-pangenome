# کدهای استفاده‌شده — مرحلهٔ ۰۹: سیگنال زمانی

## ۱) استخراج تاریخ از نام Isolate (روش اولیه، ناکافی)

```python
import pandas as pd, re

xl = pd.read_excel("05_Phylogeny/iTOL_Inputs_and_Scripts/Supplementary_File1.xlsx",
                   sheet_name=None, header=2)

def extract_year_from_name(text):
    if pd.isna(text): return None
    m = re.findall(r'(19[5-9]\d|20[0-2]\d)', str(text))
    return float(m[0]) if m else None

# نتیجه: فقط 156 از 290 (53.8%) -- ناکافی، به سراغ NCBI رفتیم
```

## ۲) پرس‌وجوی مستقیم NCBI Entrez برای تاریخ جمع‌آوری واقعی

```python
from Bio import Entrez, SeqIO
from pathlib import Path
import time, json

Entrez.email = "xpkk@zhaw.ch"   # الزام NCBI: ایمیل واقعی

genome_ids = [f.stem for f in Path("02_OrthoFinder_Input/OrthoFinder_Input_v2").glob("*.faa")]
genbank_ids = [g for g in genome_ids if not g.startswith("SRR")]

results = {}
failed = []

for i, acc in enumerate(genbank_ids):
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        collection_date = None
        for feature in record.features:
            if feature.type == "source":
                cd = feature.qualifiers.get("collection_date", [None])[0]
                if cd:
                    collection_date = cd
                break

        results[acc] = collection_date
        time.sleep(0.34)   # محدودیت نرخ NCBI بدون API key: حداکثر 3 در ثانیه
    except Exception as e:
        failed.append((acc, str(e)))
        time.sleep(0.34)

with open("09_Temporal_Signal/genbank_collection_dates.json", "w") as f:
    json.dump(results, f, indent=2)
with open("09_Temporal_Signal/genbank_failed.json", "w") as f:
    json.dump(failed, f, indent=2)

# نتیجه: 201 از 237 موفق، 0 خطا
```

## ۳) پارس فرمت‌های ناهمگون تاریخ به سال اعشاری

```python
MONTHS = {"Jan":1,"Feb":2,"Mar":3,"Apr":4,"May":5,"Jun":6,
          "Jul":7,"Aug":8,"Sep":9,"Oct":10,"Nov":11,"Dec":12}

def parse_to_decimal_year(raw):
    if not raw or raw == "None":
        return None
    raw = str(raw).strip()

    # بازهٔ دو سال: "2021/2023" -> میانگین
    if "/" in raw:
        parts = raw.split("/")
        years = [int(p) for p in parts if re.fullmatch(r"\d{4}", p)]
        if len(years) == 2:
            return sum(years) / 2
        return None

    # تاریخ کامل: "2017-06-12"
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", raw)
    if m:
        y, mo, d = int(m[1]), int(m[2]), int(m[3])
        return y + (mo - 1)/12 + d/365

    # روز-ماه-سال: "21-Sep-2022"  <-- رایج‌ترین فرمت در این دیتاست
    m = re.fullmatch(r"(\d{1,2})-([A-Za-z]{3})-(\d{4})", raw)
    if m and m[2] in MONTHS:
        d, mo, y = int(m[1]), MONTHS[m[2]], int(m[3])
        return y + (mo - 1)/12 + d/365

    # ماه-سال: "Sep-2021"
    m = re.fullmatch(r"([A-Za-z]{3})-(\d{4})", raw)
    if m and m[1] in MONTHS:
        return int(m[2]) + (MONTHS[m[1]] - 1)/12

    # فقط سال: "2022"
    if re.fullmatch(r"\d{4}", raw):
        return float(raw)

    return None
```

## ۴) ترکیب همهٔ منابع و تطبیق با ۲۹۰ ژنوم نهایی

```python
import json
import pandas as pd
from pathlib import Path

with open("09_Temporal_Signal/genbank_collection_dates.json") as f:
    genbank_raw = json.load(f)

genbank_years = {acc: parse_to_decimal_year(raw)
                 for acc, raw in genbank_raw.items()
                 if parse_to_decimal_year(raw)}

xl = pd.read_excel("05_Phylogeny/iTOL_Inputs_and_Scripts/Supplementary_File1.xlsx",
                   sheet_name=None, header=2)

# SRA: استخراج از نام (نه efetch مستقیم -- نیاز به BioSample/elink دارد)
df_r = xl['Raw_Data']
df_r = df_r[df_r['Run'].notna()]
sra_years = {}
for _, row in df_r.iterrows():
    run = str(row['Run']).strip()
    for col in ['Isolate', 'Experiment_title', 'Description', 'SampleName']:
        y = extract_year_from_name(row.get(col))
        if y:
            sra_years[run] = y
            break

# GenBank fallback: فقط برای accession هایی که efetch تاریخ نداد
df_a = xl['Assembeled']
df_a = df_a[df_a['Accession Number'].notna()]
genbank_fallback = {}
for _, row in df_a.iterrows():
    acc = str(row['Accession Number']).strip()
    if acc not in genbank_years:
        y = extract_year_from_name(row.get('Isolate'))
        if y:
            genbank_fallback[acc] = y

# اولویت: efetch واقعی > نام Isolate
all_years = {**genbank_fallback, **sra_years, **genbank_years}

# تطبیق با ۲۹۰ ژنوم نهایی (شامل حل نمونه‌های _Merged)
genome_ids = [f.stem for f in Path("02_OrthoFinder_Input/OrthoFinder_Input_v2").glob("*.faa")]
matched = {}
for gid in genome_ids:
    if gid in all_years:
        matched[gid] = all_years[gid]
    elif "_Merged" in gid:
        for part in gid.replace("_Merged", "").split("_"):
            if part in all_years:
                matched[gid] = all_years[part]
                break

with open("09_Temporal_Signal/all_tip_years.json", "w") as f:
    json.dump(matched, f, indent=2)

# نتیجه: 234 از 290 (80.7%)
```

## ۵) Root-to-tip Regression

```python
import json
from Bio import Phylo
from scipy import stats
import pandas as pd

with open("09_Temporal_Signal/all_tip_years.json") as f:
    tip_years = json.load(f)

# نکتهٔ مهم: درخت اصلاح‌شدهٔ نوترکیبی (Gubbins)، نه درخت خام Parsnp
# چون این تحلیل به فرض «تجمع یکنواخت جهش با زمان» تکیه دارد که
# نوترکیبی آن را نقض می‌کند
tree = Phylo.read("07_Phylogeography/lsdv_rerooted_with_SPPV.nwk", "newick")
sppv = next((t for t in tree.get_terminals() if t.name == "NC_004002.1"), None)
if sppv:
    tree.prune(sppv)   # حذف outgroup از خودِ تحلیل

root = tree.root
root_to_tip = {tip.name: tree.distance(root, tip) for tip in tree.get_terminals()}

# فقط تاکسون‌هایی که هم سال دارند هم در درخت هستند
paired = [(tip_years[t], root_to_tip[t]) for t in root_to_tip if t in tip_years]

years = [p[0] for p in paired]
dists = [p[1] for p in paired]

slope, intercept, r_value, p_value, std_err = stats.linregress(years, dists)

print(f"n={len(paired)}  R²={r_value**2:.4f}  p={p_value:.4e}  slope={slope:.4e}")

out_df = pd.DataFrame(paired, columns=["year", "root_to_tip_distance"])
out_df.to_csv("09_Temporal_Signal/root_to_tip_data.csv", index=False)

with open("09_Temporal_Signal/regression_stats.json", "w") as f:
    json.dump({"slope": slope, "intercept": intercept, "r_squared": r_value**2,
              "p_value": p_value, "std_err": std_err, "n": len(paired)}, f, indent=2)
```

## ۶) رسم نمودار

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json

df = pd.read_csv("root_to_tip_data.csv")
with open("regression_stats.json") as f:
    stats = json.load(f)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(df["year"], df["root_to_tip_distance"], alpha=0.5, s=20, color="#2E6E8E")

x_line = np.linspace(df["year"].min(), df["year"].max(), 100)
y_line = stats["slope"] * x_line + stats["intercept"]
ax.plot(x_line, y_line, color="#C1440E", linewidth=2,
        label=f"R²={stats['r_squared']:.4f}, p={stats['p_value']:.3f}")

ax.set_xlabel("Collection Year")
ax.set_ylabel("Root-to-tip Distance")
ax.set_title("Temporal Signal Analysis (n=234)")
ax.legend()
plt.tight_layout()
plt.savefig("root_to_tip_regression.png", dpi=150)
```
