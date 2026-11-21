#!/bin/bash
#SBATCH --job-name=rev27_scoary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --output=logs/scoary_%j.out
#SBATCH --error=logs/scoary_%j.err
#
# Submit from the stage directory, so that logs/ already exists:
#   cd $PANG/revision_2026-11/27_DNA_Level_GWAS
#   sbatch scripts/04_run_scoary.sbatch
#
# The same Scoary invocation as the original analysis -- tree-constrained
# pairwise comparisons, 1000 permutations -- run on the DNA-level matrix
# instead of the annotation-based one. Nothing else differs between the two
# runs, so any change in the result comes from the matrix.

set -euo pipefail
export PANG="${PANG:-/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis}"
STAGE="$PANG/revision_2026-11/27_DNA_Level_GWAS"
OUT="$STAGE/02_results/scoary_results_dna"
mkdir -p "$OUT"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scoary_env       # the pinned Python 3.6 environment

# --start_col is the first genome column, i.e. one past the metadata block.
# Step 3 writes the same metadata columns as the original input, so the
# value used for the original run applies unchanged here.
START_COL="${START_COL:-15}"

scoary \
  --genes   "$STAGE/03_summary/gene_presence_absence_dna.csv" \
  --traits  "$STAGE/03_summary/traits_matched.csv" \
  --newicktree "$STAGE/03_summary/core_tree_pruned.newick" \
  --start_col "$START_COL" \
  -e 1000 \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --outdir "$OUT"

echo
echo "result files:"
ls -1 "$OUT" | head -30
