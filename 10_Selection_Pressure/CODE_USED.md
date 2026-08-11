# مرحلهٔ ۱۰ — کدهای استفاده‌شده (فشار انتخابی)

## نصب محیط‌ها
```bash
conda create -n hyphy_env   -c conda-forge -c bioconda hyphy  -y   # HyPhy 2.5.101
conda create -n pal2nal_env -c bioconda -c conda-forge pal2nal -y  # pal2nal v14
```

## اجرای اول BUSTED (Job 644026)
```bash
#!/bin/bash
#SBATCH --job-name=BUSTED_LSDV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=busted_run.log

module load USS/2022
module load slurm
source ~/.bashrc
conda activate hyphy_env

BASE="/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis/10_Selection_Pressure"

for codon_file in "$BASE"/01_codon_alignments/*_codon.fasta; do
    og=$(basename "$codon_file" _codon.fasta)
    tree_file="$BASE/02_pruned_trees/${og}_pruned.nwk"
    out_json="$BASE/03_busted_results/${og}_BUSTED.json"

    if [[ -f "$out_json" ]]; then
        echo "SKIP (موجود): $og"
        continue
    fi

    echo "=== $og ==="
    hyphy busted \
        --alignment "$codon_file" \
        --tree "$tree_file" \
        --output "$out_json" \
        CPU=4 \
        > "$BASE/03_busted_results/${og}_BUSTED.log" 2>&1

    if [[ -f "$out_json" ]]; then
        echo "  ✅ موفق"
    else
        echo "  ❌ شکست -- بررسی ${og}_BUSTED.log"
    fi
done

echo "=================================================="
echo "🏁 BUSTED finished for all genes"
echo "=================================================="
```

## اجرای مجدد ۴ ژن ناموفق (Job 647755)
تفاوت کلیدی: افزودن `ENV=TOLERATE_NUMERICAL_ERRORS=1`
```bash
#!/bin/bash
#SBATCH --job-name=BUSTED_retry
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=06:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=10_Selection_Pressure/busted_retry.log

module load USS/2022
module load slurm
source ~/.bashrc
conda activate hyphy_env

BASE="/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis/10_Selection_Pressure"

for og in OG0000070 OG0000075 OG0000091 OG0000121; do
    echo "=== $og ==="
    # فایل خالی قبلی را حذف کن تا با نتیجهٔ جدید اشتباه نشود
    rm -f "$BASE/03_busted_results/${og}_BUSTED.json"

    hyphy busted \
        --alignment "$BASE/01_codon_alignments/${og}_codon.fasta" \
        --tree "$BASE/02_pruned_trees/${og}_pruned.nwk" \
        --output "$BASE/03_busted_results/${og}_BUSTED.json" \
        CPU=4 \
        ENV=TOLERATE_NUMERICAL_ERRORS=1 \
        > "$BASE/03_busted_results/${og}_BUSTED_retry.log" 2>&1

    if [[ -s "$BASE/03_busted_results/${og}_BUSTED.json" ]]; then
        echo "  ✅ موفق ($(stat -c%s "$BASE/03_busted_results/${og}_BUSTED.json") بایت)"
    else
        echo "  ❌ باز هم شکست — بررسی ${og}_BUSTED_retry.log"
        tail -20 "$BASE/03_busted_results/${og}_BUSTED_retry.log"
    fi
done

echo "=================================================="
echo "🏁 retry finished"
echo "=================================================="
```

## پس‌پردازش
اسکریپت تجمیع‌شدهٔ Python: `scripts/parse_busted.py`
