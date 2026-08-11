#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --job-name=absrel_stem
#SBATCH --time=12:00:00
#SBATCH --array=1-60%12
#SBATCH --output=16_aBSREL_stem/logs/absrel_%A_%a.out
#SBATCH --error=16_aBSREL_stem/logs/absrel_%A_%a.err
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch

set -u
cd "$SLURM_SUBMIT_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hyphy_env

OG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 16_aBSREL_stem/01_labeled_trees/gene_list.txt)
OUT=16_aBSREL_stem/02_results/${OG}_ABSREL.json

if [[ -s "$OUT" ]]; then
    echo "$OG: از قبل موجود، پرش"; exit 0
fi

echo "شروع $OG در $(date)"
hyphy absrel \
    --alignment 10_Selection_Pressure/01_codon_alignments/${OG}_codon.fasta \
    --tree      16_aBSREL_stem/01_labeled_trees/${OG}_stem.nwk \
    --branches  Stem \
    --output    "$OUT"
rc=$?
if [[ -s "$OUT" ]]; then echo "$OG: ✅ پایان $(date)"; else echo "$OG: ❌ rc=$rc"; fi
exit $rc
