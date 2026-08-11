#!/bin/bash
#SBATCH --job-name=Scoary_tree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=8
#SBATCH --mem=30GB
#SBATCH --time=06:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=11_Scoary/logs/scoary_tree_run.log

module load USS/2022
module load slurm
source ~/.bashrc
conda activate scoary_env

BASE="/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis"
cd "$BASE/11_Scoary"

rm -rf scoary_results_tree

scoary \
    --genes      01_inputs/gene_presence_absence.csv \
    --traits     01_inputs/traits.csv \
    --newicktree "$BASE/05_Phylogeny/04_tree/core_genome_ML_tree_v3.contree" \
    --outdir     scoary_results_tree \
    --start_col  15 \
    --no-time \
    -e 1000 \
    --threads 8

echo "کد خروج: $?"
echo "تعداد فایل نتیجه: $(ls scoary_results_tree 2>/dev/null | wc -l)"
echo "=================================================="
echo "🏁 Scoary with tree finished"
echo "=================================================="
