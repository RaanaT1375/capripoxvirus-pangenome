# کدهای استفاده‌شده — مرحلهٔ ۰۶: تحلیل نوترکیبی

## ۱) جمع‌آوری ژنوم‌های نوکلئوتیدی

```python
import subprocess, shutil
from pathlib import Path

# لیست ۲۹۰ ژنوم نهایی از پوشهٔ ورودی OrthoFinder
genome_ids = [f.stem for f in
              Path("02_OrthoFinder_Input/OrthoFinder_Input_v2").glob("*.faa")]

outdir = Path("06_Recombination/00_nucleotide_genomes")
outdir.mkdir(parents=True, exist_ok=True)
found, missing = [], []

# جست‌وجوی فایل .fna متناظر در درخت پروژه
for gid in genome_ids:
    result = subprocess.run(
        ["find", "/cfs/earth/scratch/xpkk/Raana/Cow",
         "/cfs/earth/scratch/xpkk/Raana/Sheep",
         "/cfs/earth/scratch/xpkk/Raana/Goat",
         "-type", "f", "-iname", f"{gid}.fna"],
        capture_output=True, text=True
    )
    hits = [l for l in result.stdout.strip().split("\n") if l]
    if hits:
        found.append((gid, hits[0]))
    else:
        missing.append(gid)

# کپی + ساخت فایل نگاشت
with open("06_Recombination/nucleotide_file_map.tsv", "w") as fh:
    fh.write("genome_id\tsource_path\n")
    for gid, path in found:
        shutil.copy(path, outdir / f"{gid}.fasta")
        fh.write(f"{gid}\t{path}\n")

print(f"پیدا‌شده: {len(found)}  |  پیدا‌نشده: {len(missing)}")
```

نمونه‌های `_Merged` اصلاح‌شده مسیر متفاوتی دارند:

```python
CORRECTED = Path(".../Merged_Amplicons_CORRECTED")
mapping = {
    "SRR26747378_SRR26747377_Merged": "SRR26747378_SRR26747377_Merged_CORRECTED",
    "SRR26747380_SRR26747379_Merged": "SRR26747380_SRR26747379_Merged_CORRECTED",
    "SRR26747382_SRR26747381_Merged": "SRR26747382_SRR26747381_Merged_CORRECTED",
    "SRR26747384_SRR26747383_Merged": "SRR26747384_SRR26747383_Merged_CORRECTED",
    "SRR26747386_SRR26747385_Merged": "SRR26747386_SRR26747385_Merged_CORRECTED",
    "SRR26747387_SRR26747376_Merged": "SRR26747387_SRR26747376_Merged_CORRECTED",
}

with open("06_Recombination/nucleotide_file_map.tsv", "a") as fh:
    for final_name, folder_name in mapping.items():
        fna = CORRECTED / folder_name / "Prokka_Output" / f"{folder_name}.fna"
        if fna.exists():
            shutil.copy(fna, outdir / f"{final_name}.fasta")
            fh.write(f"{final_name}\t{fna}\n")
```

## ۲) الاینمنت کامل ژنومی (Parsnp)

```bash
conda activate parsnp_env

REF="/cfs/earth/scratch/xpkk/Raana/LSDV_Ref_Genome/Cow_LSDV_Ref.fasta"

# مهم: پوشهٔ خروجی نباید از قبل وجود داشته باشد
rm -rf 01_Parsnp_Alignment

parsnp -r "$REF" \
       -d 00_nucleotide_genomes \
       -o 01_Parsnp_Alignment \
       -p 30 \
       -c
```

پارامترها:
- `-r` : ژنوم رفرنس
- `-d` : پوشهٔ ژنوم‌های ورودی
- `-o` : پوشهٔ خروجی (نباید از قبل موجود باشد)
- `-p` : تعداد thread
- `-c` : اجبار مشارکت همهٔ ژنوم‌ها در الاینمنت

### تبدیل به FASTA
```bash
harvesttools -i 01_Parsnp_Alignment/parsnp.ggr \
             -M 01_Parsnp_Alignment/parsnp_core_alignment.fasta

# بررسی
grep -c '^>' 01_Parsnp_Alignment/parsnp_core_alignment.fasta

# تأیید یکسان بودن طول همهٔ توالی‌ها (پیش‌نیاز Gubbins)
awk '/^>/{if(seq)print length(seq); seq=""} !/^>/{seq=seq $0} END{print length(seq)}' \
    01_Parsnp_Alignment/parsnp_core_alignment.fasta | sort -u | wc -l
# خروجی باید 1 باشد
```

## ۳) اجرای Gubbins (SLURM)

```bash
#!/bin/bash
#SBATCH --job-name=Gubbins_v3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=gubbins_v3_run.log

module load USS/2022
module load slurm
source ~/.bashrc
conda activate gubbins_fix   # نه gubbins_env (باگ pkg_resources)

cd .../06_Recombination/02_Gubbins_Results

run_gubbins.py \
    --prefix lsdv_gubbins \
    --threads 30 \
    --tree-builder raxml \
    --iterations 5 \
    ../01_Parsnp_Alignment/parsnp_core_alignment.fasta
```

## ۴) تحلیل آماری نتایج نوترکیبی

```python
import pandas as pd

df = pd.read_csv("lsdv_gubbins.per_branch_statistics.csv", sep="\t")
df['is_internal'] = df['Node'].astype(str).str.startswith('Node')

# --- r/m سراسری درخت ---
total_in  = df['Number of SNPs Inside Recombinations'].sum()
total_out = df['Number of SNPs Outside Recombinations'].sum()
print(f"SNP داخل: {total_in}  |  خارج: {total_out}  |  r/m = {total_in/total_out:.2f}")

# --- اتحاد بازه‌های نوترکیبی (بدون شمارش مضاعف) ---
intervals = []
with open("lsdv_gubbins.recombination_predictions.gff") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        intervals.append((int(parts[3]), int(parts[4])))

intervals.sort()
merged = []
for s, e in intervals:
    if merged and s <= merged[-1][1] + 1:
        merged[-1] = (merged[-1][0], max(merged[-1][1], e))
    else:
        merged.append([s, e])

union = sum(e - s + 1 for s, e in merged)
print(f"نواحی: {len(intervals)}  |  اتحاد: {union} bp ({union/110321*100:.1f}%)")

# --- شاخه‌های انتهایی با بیشترین نوترکیبی مختص‌به‌خود ---
# نکته: ستون Cumulative تجمعی از ریشه است — نباید استفاده شود
terminal = df[~df['is_internal']].sort_values(
    'Bases in Recombinations Excluding Gaps', ascending=False)
print(terminal[['Node', 'Number of Recombination Blocks',
                'Bases in Recombinations Excluding Gaps', 'r/m']].head(10))

# --- شاخه‌های داخلی (اجدادی) ---
internal = df[df['is_internal']].sort_values(
    'Bases in Recombinations Excluding Gaps', ascending=False)
print(internal[['Node', 'Number of Recombination Blocks',
                'Bases in Recombinations Excluding Gaps']].head(5))
```

## ۵) اعتبارسنجی متقابل با برچسب‌های پیشین

```python
from Bio import Phylo

tree = Phylo.read("lsdv_gubbins.node_labelled.final_tree.tre", "newick")

recombinant_labeled = ["MH646674.1", "OL752713.2", "MZ577076.1", "MZ577075.1",
                       "MZ577074.1", "MZ577073.1", "MW732649.1", "MW355944.1",
                       "MT992618.1", "MT134042.1"]
top_internal = ["Node_143", "Node_146", "Node_144", "Node_147", "Node_142"]

for target in top_internal:
    clade = next((c for c in tree.find_clades(name=target)), None)
    if clade is None:
        continue
    descendants = [t.name.replace(".fasta", "") for t in clade.get_terminals()]
    overlap = [s for s in recombinant_labeled if s in descendants]
    print(f"{target}: {len(descendants)} نواده — {len(overlap)} سویه: {overlap}")
```

## ۶) پاک‌سازی درخت نهایی

```python
from Bio import Phylo
import re

tree = Phylo.read("lsdv_gubbins.final_tree.tre", "newick")

# حذف برگ رفرنس
ref_tip = next((t for t in tree.get_terminals() if "ref" in t.name.lower()), None)
if ref_tip:
    tree.prune(ref_tip)

# حذف پسوند .fasta از نام برگ‌ها
for t in tree.get_terminals():
    t.name = re.sub(r"\.fasta$", "", t.name)

Phylo.write(tree, "lsdv_gubbins.final_tree_renamed.nwk", "newick")
print(f"تعداد برگ نهایی: {len(tree.get_terminals())}")
```

## ۷) ساخت لایهٔ iTOL وضعیت نوترکیبی

```python
from Bio import Phylo
from collections import Counter

# استخراج نوادگان Node_147 (کلاد نوترکیبی اجدادی)
tree_labelled = Phylo.read("02_Gubbins_Results/lsdv_gubbins.node_labelled.final_tree.tre",
                           "newick")
node147 = next(tree_labelled.find_clades(name="Node_147"))
clade_members = set(t.name.replace(".fasta", "") for t in node147.get_terminals())
clade_members.discard("Cow_LSDV_Ref.fasta.ref")

tree_final = Phylo.read("02_Gubbins_Results/lsdv_gubbins.final_tree_renamed.nwk",
                        "newick")
all_genomes = [t.name for t in tree_final.get_terminals()]

INDEPENDENT = {"MT134042.1"}
colors = {
    "Ancestral recombinant clade":      "#8B5E3C",
    "Independent recent recombination": "#A63446",
    "No major recombination detected":  "#D8D4C8",
}

rows = {}
for g in all_genomes:
    if g in INDEPENDENT:
        cat = "Independent recent recombination"
    elif g in clade_members:
        cat = "Ancestral recombinant clade"
    else:
        cat = "No major recombination detected"
    rows[g] = colors[cat]

header = [
    "DATASET_COLORSTRIP", "SEPARATOR TAB",
    "DATASET_LABEL\tRecombination Status",
    "COLOR\t#000000", "STRIP_WIDTH\t30", "MARGIN\t2", "SHOW_INTERNAL\t0",
    "LEGEND_TITLE\tRecombination (Gubbins)",
    "LEGEND_SHAPES\t1\t1\t1",
    f"LEGEND_COLORS\t{colors['Ancestral recombinant clade']}\t"
    f"{colors['Independent recent recombination']}\t"
    f"{colors['No major recombination detected']}",
    "LEGEND_LABELS\tAncestral recombinant clade\t"
    "Independent recent recombination\tNo major recombination detected",
    "DATA",
]

with open("03_iTOL_Recombination/itol_recombination_status.txt", "w") as fh:
    for line in header:
        fh.write(line + "\n")
    for g, color in rows.items():
        fh.write(f"{g}\t{color}\n")

label_map = {v: k for k, v in colors.items()}
print(Counter(label_map[c] for c in rows.values()))
```

## ۸) ساخت بقیهٔ فایل‌های iTOL

```bash
python3 ../05_Phylogeny/iTOL_Inputs_and_Scripts/prepare_itol_datasets.py \
    --tree 02_Gubbins_Results/lsdv_gubbins.final_tree_renamed.nwk \
    --supplementary ../05_Phylogeny/iTOL_Inputs_and_Scripts/Supplementary_File1.xlsx \
    --itol-dir ../05_Phylogeny/iTOL_Inputs_and_Scripts \
    --outdir 03_iTOL_Recombination
```
