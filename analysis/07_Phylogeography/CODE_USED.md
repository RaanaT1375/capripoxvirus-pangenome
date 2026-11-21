# کدهای استفاده‌شده — مرحلهٔ ۰۷: فیلوژئوگرافی

## ۱) ساخت جدول متادیتای جغرافیایی

```python
# استخراج کشور و قاره از فایل‌های iTOL مرحلهٔ نوترکیبی
country_map = {}
with open("06_Recombination/03_iTOL_Recombination/itol_country_labels_FINAL.txt") as fh:
    in_data = False
    for line in fh:
        line = line.rstrip("\n")
        if line.strip() == "DATA":
            in_data = True
            continue
        if not in_data or not line.strip():
            continue
        parts = line.split("\t")
        country_map[parts[0]] = parts[1]

# قاره: باید از کد رنگ به برچسب ترجمه شود
legend_colors, legend_labels = [], []
continent_data = {}
with open("06_Recombination/03_iTOL_Recombination/itol_continent_FINAL.txt") as fh:
    in_data = False
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("LEGEND_COLORS"):
            legend_colors = line.split("\t")[1:]
        if line.startswith("LEGEND_LABELS"):
            legend_labels = line.split("\t")[1:]
        if line.strip() == "DATA":
            in_data = True
            continue
        if not in_data or not line.strip():
            continue
        parts = line.split("\t")
        continent_data[parts[0]] = parts[1]

color_to_label = dict(zip(legend_colors, legend_labels))
all_tips = set(country_map) | set(continent_data)

with open("07_Phylogeography/geo_metadata.csv", "w") as fh:
    fh.write("tip,Country,Continent\n")
    for tip in sorted(all_tips):
        country = country_map.get(tip, '')
        continent = color_to_label.get(continent_data.get(tip, ''), '')
        fh.write(f"{tip},{country},{continent}\n")

print(f"تعداد ژنوم: {len(all_tips)}")
```

## ۲) بررسی توزیع قاره‌ای

```python
import csv
from collections import Counter

c = Counter()
with open('geo_metadata.csv') as fh:
    for row in csv.DictReader(fh):
        c[row['Continent']] += 1

for k, v in c.most_common():
    print(f'  {k}: {v} ({v/sum(c.values())*100:.1f}%)')
```

## ۳) اعتبارسنجی انتخاب outgroup

```python
from Bio import Phylo

tree = Phylo.read("../06_Recombination/02_Gubbins_Results/lsdv_gubbins.final_tree_renamed.nwk",
                  "newick")

# فاصلهٔ هر outgroup کاندید تا چند سویهٔ LSDV
sample = [x.name for x in tree.get_terminals()
          if x.name not in ("NC_004002.1", "NC_004003.1")][:5]

for s in sample:
    print(f"  SPPV <-> {s}: {tree.distance('NC_004002.1', s):.5f}")
    print(f"  GTPV <-> {s}: {tree.distance('NC_004003.1', s):.5f}")

# مقایسه با درخت پروتئینی (بدون تأثیر Gubbins)
t_prot = Phylo.read("../05_Phylogeny/04_tree/core_genome_ML_tree_v3.contree", "newick")
print(f"درخت پروتئینی — SPPV↔GTPV: {t_prot.distance('NC_004002.1','NC_004003.1'):.5f}")
```

## ۴) ریشه‌دهی با outgroup

```python
from Bio import Phylo

tree = Phylo.read("../06_Recombination/02_Gubbins_Results/lsdv_gubbins.final_tree_renamed.nwk",
                  "newick")

# ریشه‌دهی با SPPV
og = next(t for t in tree.get_terminals() if t.name == "NC_004002.1")
tree.root_with_outgroup(og)
Phylo.write(tree, "lsdv_rerooted_with_SPPV.nwk", "newick")

# حذف outgroup از تحلیل (چون خودش بخشی از سؤال منشأ LSDV نیست)
tree2 = Phylo.read("lsdv_rerooted_with_SPPV.nwk", "newick")
sppv = next(t for t in tree2.get_terminals() if t.name == "NC_004002.1")
tree2.prune(sppv)
Phylo.write(tree2, "lsdv_ingroup_only_rooted.nwk", "newick")

print(f"تعداد برگ نهایی: {len(tree2.get_terminals())}")
```

## ۵) اجرای PastML

```bash
conda activate pastml_env

pastml \
    --tree lsdv_ingroup_only_rooted.nwk \
    --data geo_metadata.csv \
    --data_sep , \
    --id_index 0 \
    --columns Country Continent \
    --prediction_method MPPA \
    --work_dir pastml_output_v3 \
    --html_compressed pastml_output_v3/tree_geography_v3.html
```

پارامترها:
- `--columns` : ستون‌های متادیتا برای بازسازی
- `--prediction_method MPPA` : اجازهٔ حالت‌های چندگانه در گره‌های مبهم
  (در مقابل MAP که همیشه یک حالت انتخاب می‌کند و رویداد کاذب می‌سازد)
- `--html_compressed` : تجسم تعاملی

## ۶) استخراج حالت ریشه و شمارش رویدادهای انتقال

**نکتهٔ مهم:** فایل `combined_ancestral_states.tab` فرمت «طولانی» دارد —
وقتی گره‌ای مبهم است، PastML چند سطر جداگانه برای همان گره می‌نویسد.
پارس ساده (یک سطر = یک گره) نتیجهٔ غلط می‌دهد.

```python
from collections import Counter, defaultdict
from Bio import Phylo

# پارس صحیح: تجمیع همهٔ سطرهای هر گره در یک مجموعه
node_continents = defaultdict(set)
with open("pastml_output_v3/combined_ancestral_states.tab") as fh:
    fh.readline()   # skip header
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 1 and parts[1].strip():
            node_continents[parts[0]].add(parts[1].strip())

print(f"حالت اجدادی ریشه: {sorted(node_continents.get('root', set()))}")

# شمارش رویدادهای انتقال، فقط روی یال‌هایی که هر دو سر تک‌حالته‌اند
tree = Phylo.read("pastml_output_v3/named.tree_lsdv_ingroup_only_rooted.nwk",
                  "newick")
transitions = Counter()
resolved, ambiguous = 0, 0

def walk(clade, parent_id=None):
    global resolved, ambiguous
    st = node_continents.get(clade.name, set())
    if parent_id is not None:
        ps = node_continents.get(parent_id, set())
        if len(st) == 1 and len(ps) == 1:
            p, c = next(iter(ps)), next(iter(st))
            if p != c:
                transitions[(p, c)] += 1
            resolved += 1
        else:
            ambiguous += 1
    for ch in clade.clades:
        walk(ch, clade.name)

walk(tree.root, None)

print(f"یال قابل‌تفسیر: {resolved}  |  مبهم: {ambiguous}")
for (p, c), n in transitions.most_common():
    print(f"  {p} -> {c}: {n}")
```

## ۷) تحلیل حساسیت — نمونه‌برداری متوازن

هدف: بررسی این‌که آیا حالت ریشه به عدم‌توازن نمونه‌برداری حساس است.

```python
import csv, random, subprocess
from collections import Counter, defaultdict
from pathlib import Path
from Bio import Phylo

# بارگذاری متادیتا
geo = {}
with open("geo_metadata.csv") as fh:
    for row in csv.DictReader(fh):
        if row["Continent"]:
            geo[row["tip"]] = row["Continent"]

by_continent = defaultdict(list)
for tip, cont in geo.items():
    by_continent[cont].append(tip)

# اندازهٔ نمونه = کمترین تعداد در بین قاره‌ها (متوازن‌سازی)
n_per_continent = min(len(v) for v in by_continent.values())
print(f"نمونه به‌ازای هر قاره: {n_per_continent}")

outdir = Path("sensitivity_analysis")
outdir.mkdir(exist_ok=True)
root_states = []

for rep in range(10):
    random.seed(rep)   # seed متفاوت برای هر تکرار
    sampled = []
    for cont, tips in by_continent.items():
        sampled.extend(random.sample(tips, n_per_continent))

    # هرس درخت به همین زیرمجموعه
    tree = Phylo.read("lsdv_ingroup_only_rooted.nwk", "newick")
    for t in [x for x in tree.get_terminals() if x.name not in set(sampled)]:
        tree.prune(t)

    tree_path = outdir / f"rep{rep}_tree.nwk"
    Phylo.write(tree, tree_path, "newick")

    meta_path = outdir / f"rep{rep}_meta.csv"
    with open(meta_path, "w") as fh:
        fh.write("tip,Continent\n")
        for tip in sampled:
            fh.write(f"{tip},{geo[tip]}\n")

    work_dir = outdir / f"rep{rep}_pastml"
    subprocess.run([
        "pastml", "--tree", str(tree_path), "--data", str(meta_path),
        "--data_sep", ",", "--id_index", "0", "--columns", "Continent",
        "--prediction_method", "MPPA", "--work_dir", str(work_dir)
    ], capture_output=True, text=True)

    # خواندن حالت ریشه
    states_file = work_dir / "combined_ancestral_states.tab"
    root_set = set()
    if states_file.exists():
        with open(states_file) as fh:
            fh.readline()
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if parts[0] == "root" and len(parts) > 1 and parts[1].strip():
                    root_set.add(parts[1].strip())
    root_states.append(tuple(sorted(root_set)))
    print(f"  تکرار {rep}: ریشه = {sorted(root_set)}")

print("\n=== خلاصه ===")
for state, count in Counter(root_states).most_common():
    print(f"  {list(state)}: {count}/10")
```

## ۸) تطبیق نوترکیبی با جغرافیا

```python
import csv
from collections import defaultdict, Counter

geo = {}
with open("geo_metadata.csv") as fh:
    for row in csv.DictReader(fh):
        geo[row["tip"]] = row

recomb_status = {}
with open("../06_Recombination/03_iTOL_Recombination/itol_recombination_status.txt") as fh:
    in_data = False
    for line in fh:
        line = line.rstrip("\n")
        if line.strip() == "DATA":
            in_data = True
            continue
        if not in_data or not line.strip():
            continue
        tip, color = line.split("\t")
        recomb_status[tip] = color

color_labels = {
    "#8B5E3C": "Ancestral recombinant clade",
    "#A63446": "Independent recent recombination",
    "#D8D4C8": "No major recombination detected",
}

continent_by_group = defaultdict(Counter)
for tip, color in recomb_status.items():
    label = color_labels.get(color, color)
    continent = geo.get(tip, {}).get("Continent", "Unknown")
    continent_by_group[label][continent] += 1

for group, counts in continent_by_group.items():
    print(f"=== {group} ===")
    for cont, n in counts.most_common():
        print(f"   {cont}: {n}")
```
