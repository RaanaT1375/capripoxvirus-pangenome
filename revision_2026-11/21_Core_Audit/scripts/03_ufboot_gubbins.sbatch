#!/bin/bash
#SBATCH --job-name=rev21_ufboot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --output=logs/ufboot_%j.out
#SBATCH --error=logs/ufboot_%j.err
#
# Only run this if 02_branch_support.py shows that Figure 3A really is the
# Gubbins tree and that it differs from the IQ-TREE topology. It puts UFBoot
# support on the recombination-filtered alignment, which is what the figure
# caption currently claims is plotted.
#
# Submit from the stage directory so that logs/ exists:
#   cd $PANG/revision_2026-11/21_Core_Audit && sbatch scripts/03_ufboot_gubbins.sbatch

set -euo pipefail
export PANG="${PANG:-/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis}"
STAGE="$PANG/revision_2026-11/21_Core_Audit"
ALN="$PANG/06_Recombination/02_Gubbins_Results/lsdv_gubbins.filtered_polymorphic_sites.fasta"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate iqtree_env

iqtree -s "$ALN" \
       -m GTR+ASC+G4 \
       -bb 1000 -alrt 1000 \
       -nt "${SLURM_CPUS_PER_TASK:-8}" \
       -pre "$STAGE/02_results/gubbins_filtered_ufboot" \
       -redo

echo "consensus tree with support: $STAGE/02_results/gubbins_filtered_ufboot.contree"
