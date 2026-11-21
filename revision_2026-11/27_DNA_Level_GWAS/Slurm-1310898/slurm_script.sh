#!/bin/bash
#SBATCH --job-name=rev27_tblastn
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=08:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --output=logs/tblastn_%j.out
#SBATCH --error=logs/tblastn_%j.err
#
# Submit from the stage directory, so that logs/ already exists:
#   cd $PANG/revision_2026-11/27_DNA_Level_GWAS
#   sbatch scripts/02_run_tblastn.sbatch
#
# Six times the queries of stage 20 against the same database. Stage 20 took
# well under a minute, so this should too; the eight-hour wall time is only
# there so a slow queue cannot truncate the search.

set -euo pipefail
export PANG="${PANG:-/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis}"
STAGE="$PANG/revision_2026-11/27_DNA_Level_GWAS"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prokka_env

# the database is either this stage's own or stage 20's, recorded by step 1
DB="$STAGE/01_inputs/genomes_db"
if [ -f "$STAGE/01_inputs/genomes_db.reused" ]; then
  DB="$(cat "$STAGE/01_inputs/genomes_db.reused")"
fi
echo "database: $DB"

tblastn \
  -query   "$STAGE/01_inputs/queries.faa" \
  -db      "$DB" \
  -out     "$STAGE/02_results/tblastn_hits.tsv" \
  -outfmt  '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' \
  -evalue  1e-5 \
  -max_target_seqs 10000 \
  -num_threads "${SLURM_CPUS_PER_TASK:-8}"

echo "hits: $(wc -l < "$STAGE/02_results/tblastn_hits.tsv")"
echo "queries with at least one hit: $(cut -f1 "$STAGE/02_results/tblastn_hits.tsv" | sort -u | wc -l)"
