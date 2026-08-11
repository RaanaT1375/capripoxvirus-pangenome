#!/bin/bash
#SBATCH --job-name=MEME_LSDV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=24:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=12_MEME/logs/meme_run.log

module load USS/2022
module load slurm
source ~/.bashrc
conda activate hyphy_env

BASE="/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis"
SEL="$BASE/10_Selection_Pressure"
OUT="$BASE/12_MEME/01_results"

for og in OG0000120 OG0000123 OG0000110 OG0000116 OG0000103; do
    echo "=========== $og ==========="
    out_json="$OUT/${og}_MEME.json"
    [[ -s "$out_json" ]] && { echo "  SKIP (موجود)"; continue; }

    hyphy meme \
        --alignment "$SEL/01_codon_alignments/${og}_codon.fasta" \
        --tree      "$SEL/02_pruned_trees/${og}_pruned.nwk" \
        --output    "$out_json" \
        CPU=6 \
        ENV=TOLERATE_NUMERICAL_ERRORS=1 \
        > "$OUT/${og}_MEME.log" 2>&1

    if [[ -s "$out_json" ]]; then
        echo "  ✅ موفق ($(stat -c%s "$out_json") بایت)"
    else
        echo "  ❌ شکست"
        tail -20 "$OUT/${og}_MEME.log"
    fi
done

echo "=================================================="
echo "🏁 MEME finished"
echo "=================================================="
