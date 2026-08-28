# کدهای استفاده‌شده — مرحلهٔ ۰۱: کنترل کیفیت و فیلترینگ

## فیلتر ۱ — غربالگری آلودگی

### ۱.۱ ساخت دیتابیس رفرنس (DIAMOND v2.1.10)
```bash
diamond makedb \
    --in AF325528.1.faa \
    -d /tmp/lsdv_ref_db \
    --quiet
```
- `--in` : پروتئوم رفرنس LSDV (ایزولهٔ NI-2490)
- `-d`   : مسیر خروجی دیتابیس باینری

### ۱.۲ جست‌وجوی هر ژنوم در برابر رفرنس
```bash
for f in OrthoFinder_Input/*.faa; do
    diamond blastp \
        -q "$f" \
        -d /tmp/lsdv_ref_db \
        --outfmt 6 qseqid sseqid pident length qlen slen evalue \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --quiet \
        > /tmp/$(basename "$f" .faa)_vs_ref.tsv
done
```
- `--outfmt 6`        : فرمت جدولی BLAST
- `--max-target-seqs 1`: فقط بهترین hit برای هر پروتئین
- `--evalue 1e-5`     : آستانهٔ معناداری

### ۱.۳ محاسبهٔ آمار QC
```python
import pandas as pd
import os
from pathlib import Path

records = []
for faa_file in Path("OrthoFinder_Input").glob("*.faa"):
    genome_id = faa_file.stem
    n_proteins = sum(1 for line in open(faa_file) if line.startswith(">"))

    tsv_file = f"/tmp/{genome_id}_vs_ref.tsv"
    if os.path.exists(tsv_file) and os.path.getsize(tsv_file) > 0:
        hits = pd.read_csv(tsv_file, sep='\t', header=None,
                           names=['qseqid','sseqid','pident','length',
                                  'qlen','slen','evalue'])
        n_hit = len(hits)
        avg_pident = hits['pident'].mean()
    else:
        n_hit, avg_pident = 0, 0

    pct_no_hit = (1 - n_hit / n_proteins) * 100 if n_proteins > 0 else 100

    records.append({
        'genome_id': genome_id,
        'n_proteins': n_proteins,
        'n_hit': n_hit,
        'pct_no_hit': round(pct_no_hit, 2),
        'avg_pident': round(avg_pident, 2)
    })

df = pd.DataFrame(records)
print(df.describe(percentiles=[.5, .75, .9, .95, .99]))
df.to_csv("full_contamination_qc.tsv", sep='\t', index=False)
```

### ۱.۴ حذف ژنوم‌های آلوده
```python
import shutil, os

threshold = 20.0   # بر پایهٔ شکاف بین صدک ۹۰ (8.86%) و صدک ۹۵ (38.73%)
contaminated = df[df['pct_no_hit'] > threshold]['genome_id'].tolist()

with open("excluded_samples_list.txt", "w") as f:
    for g in contaminated:
        f.write(g + "\n")

os.makedirs("Excluded_Contaminated_Genomes", exist_ok=True)
for g in contaminated:
    src = f"OrthoFinder_Input/{g}.faa"
    if os.path.exists(src):
        shutil.move(src, "Excluded_Contaminated_Genomes/")

print(f"حذف‌شده: {len(contaminated)}  |  باقی‌مانده: {323 - len(contaminated)}")
```

### ۱.۵ اعتبارسنجی آستانه با گونه‌های خواهر
```bash
for ref_acc in NC_004002.1 NC_004003.1; do
    diamond blastp \
        -q "${ref_acc}.faa" \
        -d /tmp/lsdv_ref_db \
        --outfmt 6 qseqid sseqid pident \
        --max-target-seqs 1 --evalue 1e-5 --quiet \
        | awk -v acc="$ref_acc" \
              '{sum+=$3; n++} END{printf "%s: %d hits, avg_pident=%.1f%%\n", acc, n, sum/n}'
done
```

---

## فیلتر ۲ — غربالگری ژنوم ناقص

### ۲.۱ اجرای OrthoFinder روی ۲۹۷ ژنوم
```bash
orthofinder \
    -f OrthoFinder_Input \
    -t 30 -a 30 \
    -S diamond \
    -og \
    -o OrthoFinder_Results
```

### ۲.۲ شمارش ژن‌های نزدیک-هسته غایب
```python
import pandas as pd

df = pd.read_csv('Orthogroups/Orthogroups.GeneCount.tsv',
                 sep='\t').set_index('Orthogroup')
if 'Total' in df.columns:
    df = df.drop(columns=['Total'])

presence = (df > 0)
n_genomes = presence.shape[1]

# تعریف نزدیک-هسته: حضور در ≥97% ژنوم‌ها
near_core = presence[presence.sum(axis=1) >= (0.97 * n_genomes)]
print(f"تعداد orthogroup نزدیک-هسته: {len(near_core)}")

n_missing_per_genome = (~near_core).sum(axis=0).sort_values(ascending=False)
print(n_missing_per_genome.describe(percentiles=[.5, .9, .95, .97, .99]))
print(n_missing_per_genome.head(15))
```

### ۲.۳ تأیید ماهیت با موقعیت ژنومی
```python
import re

og_table = pd.read_csv('Orthogroups/Orthogroups.tsv',
                       sep='\t').set_index('Orthogroup')
ref_col = 'AF325528.1'

# نگاشت هر orthogroup به موقعیت ژنومی (از locus tag رفرنس)
og_to_position = {}
for og in near_core.index:
    cell = og_table.loc[og, ref_col]
    if pd.isna(cell):
        continue
    first_id = str(cell).split(',')[0].strip()
    m = re.search(r'_(\d+)$', first_id)
    if m:
        og_to_position[og] = int(m.group(1))

target_genomes = ['SRR26747376_SRR26747377_Merged',
                  'SRR26747386_SRR26747387_Merged',
                  'PP034673.1', 'OR602866.1',
                  'SRR19090746_SRR19090747_Merged']

for g in target_genomes:
    missing_ogs = near_core.index[~near_core[g]].tolist()
    positions = sorted([og_to_position[og] for og in missing_ogs
                        if og in og_to_position])
    print(f"--- {g}: {len(positions)} ژن غایب ---")
    print(f"  موقعیت‌ها: {positions}")
```

### ۲.۴ تأیید با تعداد کل پروتئین
```bash
for f in PP034673.1.faa OR602866.1.faa SRR19090746_SRR19090747_Merged.faa \
         SRR26747376_SRR26747377_Merged.faa SRR26747386_SRR26747387_Merged.faa; do
    n=$(grep -c '^>' "$f")
    echo "${f}: ${n} پروتئین (رفرنس سالم ≈ 155)"
done
```

### ۲.۵ حذف ژنوم‌های ناقص
```python
threshold = 9   # صدک ۹۷ به بالا
incomplete = n_missing_per_genome[n_missing_per_genome >= threshold]
incomplete.to_csv('incomplete_genomes_list.tsv', sep='\t',
                  header=['n_missing_near_core'])
print(f"نامزد حذف: {len(incomplete)}")
```

```bash
mkdir -p Excluded_Incomplete_Genomes

INCOMPLETE="SRR26747376_SRR26747377_Merged SRR26747386_SRR26747387_Merged \
PP034673.1 SRR19090748 OR602866.1 SRR19090746_SRR19090747_Merged \
OR886082.1 OR393168.1 PQ472735.1 SRR34490580"

for sample in $INCOMPLETE; do
    cp "OrthoFinder_Input/${sample}.faa" "Excluded_Incomplete_Genomes/"
    rm "OrthoFinder_Input/${sample}.faa"
done
```

---

## فیلتر ۳ — اصلاح خطای ادغام نادرست

### ۳.۱ حذف کایمراهای باقی‌مانده
```bash
WRONG_MERGES="SRR26747378_SRR26747379_Merged SRR26747380_SRR26747381_Merged \
SRR26747382_SRR26747383_Merged SRR26747384_SRR26747385_Merged"

for d in $WRONG_MERGES; do
    rm -rf "Merged_Amplicons/${d}"
    rm -f  "Merged_Amplicons/FAA_Results/${d}.faa"
done
```

### ۳.۲ ادغام صحیح، اسمبلی و انوتیشن مجدد
```bash
#!/bin/bash
THREADS=30
MEMORY=120
BASE_DIR=".../PCR_data"
OUT_DIR="$BASE_DIR/Merged_Amplicons_CORRECTED"
FAA_OUT="$OUT_DIR/FAA_Results"
mkdir -p "$FAA_OUT"

# نگاشت صحیح: نام_ایزوله  run_P1  run_P2
declare -A PAIRS=(
    ["ALB1000"]="SRR26747387 SRR26747376"
    ["GRC715"]="SRR26747386 SRR26747385"
    ["MKD5000"]="SRR26747384 SRR26747383"
    ["MKD5011"]="SRR26747382 SRR26747381"
    ["SRB4592"]="SRR26747380 SRR26747379"
    ["SRB5778"]="SRR26747378 SRR26747377"
)

for isolate in "${!PAIRS[@]}"; do
    read -r run1 run2 <<< "${PAIRS[$isolate]}"
    sample_name="${run1}_${run2}_Merged_CORRECTED"
    sample_dir="$OUT_DIR/$sample_name"
    mkdir -p "$sample_dir/merged_reads"

    # ادغام reads پاک‌شده از میزبان
    cat "$BASE_DIR/${run1}/host_depleted/${run1}_host_removed_R1.fastq.gz" \
        "$BASE_DIR/${run2}/host_depleted/${run2}_host_removed_R1.fastq.gz" \
        > "$sample_dir/merged_reads/${sample_name}_R1.fastq.gz"
    cat "$BASE_DIR/${run1}/host_depleted/${run1}_host_removed_R2.fastq.gz" \
        "$BASE_DIR/${run2}/host_depleted/${run2}_host_removed_R2.fastq.gz" \
        > "$sample_dir/merged_reads/${sample_name}_R2.fastq.gz"

    # اسمبلی
    out_assembly="$sample_dir/Spades_Assembly"
    mkdir -p "$out_assembly"
    spades.py \
        --pe1-1 "$sample_dir/merged_reads/${sample_name}_R1.fastq.gz" \
        --pe1-2 "$sample_dir/merged_reads/${sample_name}_R2.fastq.gz" \
        --only-assembler \
        --threads $THREADS --memory $MEMORY \
        -o "$out_assembly" > "$out_assembly/spades.log" 2>&1

    # فیلتر کنتیگ‌های زیر 500bp
    awk 'BEGIN{RS=">";FS="\n"} NR>1{
        seq=""; for(i=2;i<=NF;i++) seq=seq""$i;
        if(length(seq)>=500) print ">"$1"\n"seq
    }' "$out_assembly/contigs.fasta" \
        > "$out_assembly/${sample_name}_final_contigs.fasta"

    # انوتیشن
    prokka_out="$sample_dir/Prokka_Output"
    rm -rf "$prokka_out"
    prokka --outdir "$prokka_out" --prefix "$sample_name" \
        "$out_assembly/${sample_name}_final_contigs.fasta" \
        --cpus $THREADS --force --kingdom Viruses > /dev/null 2>&1

    cp "$prokka_out/${sample_name}.faa" "$FAA_OUT/"
done
```

### ۳.۳ QC نمونه‌های اصلاح‌شده
```bash
REF_FAA="AF325528.1.faa"
diamond makedb --in "$REF_FAA" -d /tmp/lsdv_ref_db_v2 --quiet

for f in Merged_Amplicons_CORRECTED/FAA_Results/*.faa; do
    base=$(basename "$f" .faa)
    diamond blastp -q "$f" -d /tmp/lsdv_ref_db_v2 \
        --outfmt 6 qseqid sseqid pident length qlen slen evalue \
        --max-target-seqs 1 --evalue 1e-5 --quiet > /tmp/${base}_vs_ref.tsv

    n_total=$(grep -c '^>' "$f")
    n_hit=$(wc -l < /tmp/${base}_vs_ref.tsv)
    avg_pident=$(awk '{sum+=$3; n++} END{printf "%.1f", sum/n}' /tmp/${base}_vs_ref.tsv)
    echo "${base}: کل=${n_total}  hit=${n_hit}  pident=${avg_pident}%"
done
```

### ۳.۴ افزودن به دیتاست نهایی
```bash
CORRECTED_DIR=".../Merged_Amplicons_CORRECTED/FAA_Results"

for f in "$CORRECTED_DIR"/*.faa; do
    base=$(basename "$f" "_CORRECTED.faa")
    cp "$f" "OrthoFinder_Input_v2/${base}.faa"
done

echo "تعداد نهایی: $(ls OrthoFinder_Input_v2/*.faa | wc -l)"
```

---

## افزودن نمونهٔ جاافتاده (SRR26792425)

```bash
FAA=".../Short_Reads_Illumina/SRR26792425/Prokka_Output/SRR26792425.faa"
REF="OrthoFinder_Input_v2/AF325528.1.faa"

diamond makedb --in "$REF" -d /tmp/lsdv_qc_ref --quiet
diamond blastp -q "$FAA" -d /tmp/lsdv_qc_ref \
    --outfmt 6 qseqid sseqid pident length qlen \
    --max-target-seqs 1 --evalue 1e-5 --quiet > /tmp/srr26792425_vs_ref.tsv

n_total=$(grep -c '^>' "$FAA")
n_hit=$(wc -l < /tmp/srr26792425_vs_ref.tsv)
avg_pident=$(awk '{sum+=$3;n++}END{printf "%.1f",sum/n}' /tmp/srr26792425_vs_ref.tsv)
echo "کل=${n_total}  hit=${n_hit}  pident=${avg_pident}%"

# QC سالم بود → افزودن
cp "$FAA" "OrthoFinder_Input_v2/SRR26792425.faa"
```
