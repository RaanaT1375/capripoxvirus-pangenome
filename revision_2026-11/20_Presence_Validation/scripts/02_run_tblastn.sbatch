#!/bin/bash
#SBATCH --job-name=rev20_tblastn
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=04:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --output=logs/tblastn_%j.out
#SBATCH --error=logs/tblastn_%j.err
# Submit from the stage directory, so that logs/ already exists:
#   cd $PANG/revision_2026-11/20_Presence_Validation
#   sbatch scripts/02_run_tblastn.sbatch

set -euo pipefail
export PANG="${PANG:-/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis}"
STAGE="$PANG/revision_2026-11/20_Presence_Validation"

tblastn \
  -query   "$STAGE/01_inputs/queries.faa" \
  -db      "$STAGE/01_inputs/genomes_db" \
  -out     "$STAGE/02_results/tblastn_hits.tsv" \
  -outfmt  '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' \
  -evalue  1e-5 \
  -max_target_seqs 10000 \
  -num_threads "${SLURM_CPUS_PER_TASK:-8}"

echo "hits: $(wc -l < "$STAGE/02_results/tblastn_hits.tsv")"
