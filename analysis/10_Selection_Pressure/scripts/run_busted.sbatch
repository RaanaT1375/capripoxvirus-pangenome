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
