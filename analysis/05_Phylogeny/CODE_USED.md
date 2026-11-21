# کدهای استفاده‌شده — مرحلهٔ ۰۵: فیلوژنی

## ۱) استخراج توالی ژن‌های تک‌کپی

```python
import pandas as pd
from Bio import SeqIO
from pathlib import Path

results_dir = Path("03_OrthoFinder_Results/OrthoFinder_Results_v3/Results_Aug04")
outdir = Path("05_Phylogeny/01_extracted_single_copy")
outdir.mkdir(parents=True, exist_ok=True)

# لیست ۶۰ orthogroup تک‌کپی
sco_list = [l.strip() for l in
            open("04_Pangenome_Statistics/strict_single_copy_orthogroups_v3.txt")]

og_table = pd.read_csv(results_dir / "Orthogroups/Orthogroups.tsv",
                       sep="\t").set_index("Orthogroup")
genome_cols = og_table.columns.tolist()

written = 0
for og in sco_list:
    if og not in og_table.index:
        continue
    row = og_table.loc[og]

    # نگاشت شناسهٔ ژن → شناسهٔ ژنوم
    gene_to_genome = {}
    for genome in genome_cols:
        cell = row[genome]
        if pd.isna(cell):
            continue
        gene_to_genome[str(cell).strip()] = genome

    seq_file = results_dir / "Orthogroup_Sequences" / f"{og}.fa"
    if not seq_file.exists():
        print(f"WARNING: {seq_file} not found")
        continue

    # تغییر هدر هر توالی به شناسهٔ ژنوم
    records_out = []
    for record in SeqIO.parse(seq_file, "fasta"):
        genome = gene_to_genome.get(record.id)
        if genome is None:
            continue
        record.id = genome
        record.description = ""
        records_out.append(record)

    # فقط اگر همهٔ ژنوم‌ها حاضر باشند بنویس
    if len(records_out) == len(genome_cols):
        SeqIO.write(records_out, outdir / f"{og}.fasta", "fasta")
        written += 1
    else:
        print(f"WARNING: {og} has {len(records_out)} seqs, "
              f"expected {len(genome_cols)}")

print(f"[DONE] {written} فایل تک‌ژنی نوشته شد")
```

## ۲) الاینمنت با MAFFT

```bash
for f in 05_Phylogeny/01_extracted_single_copy/*.fasta; do
    base=$(basename "$f" .fasta)
    mafft --auto --thread 4 --quiet "$f" > "05_Phylogeny/02_aligned/${base}.aln"
done

echo "تعداد الاینمنت‌ها: $(ls 05_Phylogeny/02_aligned/*.aln | wc -l)"
```

## ۳) ساخت سوپرماتریکس + فایل پارتیشن

```python
from Bio import SeqIO
from pathlib import Path

aln_dir = Path("05_Phylogeny/02_aligned")
aln_files = sorted(aln_dir.glob("*.aln"))

# جمع‌آوری همهٔ شناسه‌های ژنوم
all_genomes = set()
per_file = {}
for f in aln_files:
    recs = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}
    per_file[f] = recs
    all_genomes.update(recs.keys())
all_genomes = sorted(all_genomes)

# الحاق افقی + ثبت مرز پارتیشن‌ها
supermatrix = {g: [] for g in all_genomes}
partitions = []
cursor = 1
for f in aln_files:
    recs = per_file[f]
    aln_len = len(next(iter(recs.values())))
    for g in all_genomes:
        # اگر ژنومی در این ژن نبود، با gap پر کن (در عمل رخ نمی‌دهد)
        supermatrix[g].append(recs.get(g, "-" * aln_len))
    partitions.append(f"AUTO, {f.stem} = {cursor}-{cursor+aln_len-1}")
    cursor += aln_len

out_fasta = Path("05_Phylogeny/03_supermatrix/supermatrix.fasta")
with open(out_fasta, "w") as fh:
    for g in all_genomes:
        fh.write(f">{g}\n{''.join(supermatrix[g])}\n")

with open("05_Phylogeny/03_supermatrix/partitions.txt", "w") as fh:
    fh.write("\n".join(partitions) + "\n")

print(f"[DONE] {len(all_genomes)} تاکسون × {cursor-1} موقعیت")
```

## ۴) اجرای IQ-TREE (SLURM)

```bash
#!/bin/bash
#SBATCH --job-name=IQTREE_LSDV_v3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=iqtree_v3_run.log

module load USS/2022
module load slurm

cd /cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis/05_Phylogeny

/cfs/earth/scratch/xpkk/tools/bin/iqtree \
    -s 03_supermatrix/supermatrix.fasta \
    -p 03_supermatrix/partitions.txt \
    -m MFP \
    -bb 1000 \
    -nt AUTO \
    -ntmax 30 \
    -pre 04_tree/core_genome_ML_tree_v3 \
    -redo
```

پارامترها:
- `-s`  : فایل سوپرماتریکس
- `-p`  : فایل پارتیشن (مدل مستقل برای هر ژن)
- `-m MFP` : ModelFinder Plus — انتخاب خودکار بهترین مدل
- `-bb 1000` : ۱۰۰۰ تکرار ultrafast bootstrap
- `-nt AUTO -ntmax 30` : تشخیص خودکار تعداد thread بهینه

## ۵) تأیید درخت نهایی

```bash
sacct -j <JOB_ID> --format=JobID,JobName,State,ExitCode,Elapsed

grep "Input data:" 04_tree/core_genome_ML_tree_v3.iqtree

python3 -c "
from Bio import Phylo
t = Phylo.read('04_tree/core_genome_ML_tree_v3.contree', 'newick')
print(f'تعداد برگ: {len(t.get_terminals())}')
"
```

## ۶) ساخت فایل‌های iTOL

```bash
python3 iTOL_Inputs_and_Scripts/prepare_itol_datasets.py \
    --tips-from-fasta 03_supermatrix/supermatrix.fasta \
    --supplementary iTOL_Inputs_and_Scripts/Supplementary_File1.xlsx \
    --itol-dir iTOL_Inputs_and_Scripts \
    --outdir 06_iTOL_Final

# یا پس از اتمام درخت (نتیجهٔ یکسان):
python3 iTOL_Inputs_and_Scripts/prepare_itol_datasets.py \
    --tree 04_tree/core_genome_ML_tree_v3.contree \
    --supplementary iTOL_Inputs_and_Scripts/Supplementary_File1.xlsx \
    --itol-dir iTOL_Inputs_and_Scripts \
    --outdir 06_iTOL_Final
```

اسکریپت `prepare_itol_datasets.py` کارهای زیر را انجام می‌دهد:
1. خواندن نام تاکسون‌ها از درخت یا سوپرماتریکس
2. فیلتر کردن فایل‌های خام iTOL به همین تاکسون‌ها
3. حل خودکار نمونه‌های `_Merged` (که در متادیتا دو رکورد جدا دارند)
4. ساخت دیتاست گونهٔ ویروسی و برچسب خوانا از فایل Supplementary
5. اعمال پالت رنگی هماهنگ روی همهٔ دیتاست‌ها
