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
