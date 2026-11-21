#!/bin/bash
#SBATCH --job-name=BUSTED_LSDV240
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=13_Selection_LSDV/logs/busted_lsdv_run.log

module load USS/2022
module load slurm
source ~/.bashrc
conda activate hyphy_env

BASE="/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis/13_Selection_LSDV"
n_ok=0; n_fail=0

for aln in "$BASE"/01_codon_alignments/*_codon.fasta; do
    og=$(basename "$aln" _codon.fasta)
    out="$BASE/03_busted_results/${og}_BUSTED.json"
    [[ -s "$out" ]] && { echo "$og: SKIP"; continue; }

    hyphy busted \
        --alignment "$aln" \
        --tree      "$BASE/02_pruned_trees/${og}_pruned.nwk" \
        --output    "$out" \
        CPU=4 \
        ENV=TOLERATE_NUMERICAL_ERRORS=1 \
        > "$BASE/03_busted_results/${og}_BUSTED.log" 2>&1

    if [[ -s "$out" ]]; then
        echo "$og: ✅ ($(stat -c%s "$out") بایت)"; ((n_ok++))
    else
        echo "$og: ❌"; ((n_fail++)); tail -12 "$BASE/03_busted_results/${og}_BUSTED.log"
    fi
done

echo "=================================================="
echo "موفق: $n_ok | ناموفق: $n_fail"
echo "🏁 BUSTED LSDV-only finished"
echo "=================================================="
