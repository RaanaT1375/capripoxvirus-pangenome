#!/bin/bash
#SBATCH --job-name=rev25_gubbins
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=08:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --output=logs/gubbins_%j.out
#SBATCH --error=logs/gubbins_%j.err
#
# Submit from the stage directory, so that logs/ exists:
#   cd $PANG/revision_2026-11/25_Recombination_Masked_MK
#   sbatch scripts/02_run_gubbins.sbatch
#
# Gubbins is run on the 60-gene LSDV supermatrix rather than on the Parsnp
# whole-genome alignment, so that the predicted recombinant intervals come
# back in coordinates that map directly onto the codon alignments used for
# the MK test.

set -euo pipefail
export PANG="${PANG:-/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis}"
STAGE="$PANG/revision_2026-11/25_Recombination_Masked_MK"
cd "$STAGE/02_results"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gubbins_env

run_gubbins.py \
  --prefix lsdv_core_gubbins \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --filter-percentage 50 \
  --min-snps 3 \
  "$STAGE/01_inputs/lsdv_super_nt.fasta"

echo "exit=$?"
ls -la lsdv_core_gubbins.* | head -20
echo
echo "recombinant intervals predicted: $(grep -vc '^##' lsdv_core_gubbins.recombination_predictions.gff || true)"
