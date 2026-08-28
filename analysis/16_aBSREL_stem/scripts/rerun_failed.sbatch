#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --job-name=absrel_rerun
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#SBATCH --output=16_aBSREL_stem/logs/rerun_%A_%a.out
#SBATCH --error=16_aBSREL_stem/logs/rerun_%A_%a.err

set -u
cd "$SLURM_SUBMIT_DIR"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hyphy_env

GENES=(OG0000070 OG0000075 OG0000091 OG0000121)
OG=${GENES[$((SLURM_ARRAY_TASK_ID-1))]}
rm -f 16_aBSREL_stem/02_results/${OG}_ABSREL.json

export TOLERATE_NUMERICAL_ERRORS=1
hyphy absrel \
    --alignment 10_Selection_Pressure/01_codon_alignments/${OG}_codon.fasta \
    --tree      16_aBSREL_stem/01_labeled_trees/${OG}_stem.nwk \
    --branches  Stem \
    --output    16_aBSREL_stem/02_results/${OG}_ABSREL.json
echo "$OG: rc=$?"
