#!/bin/bash
#SBATCH --job-name=revision_2026_11
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --time=08:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --output=logs/revision_%j.out
#SBATCH --error=logs/revision_%j.err
#
# Runs stages 17, 18, 19 and 21. Stage 20 needs BLAST+ and is submitted
# separately (see README).
#
#   cd $PANG/revision_2026-11 && mkdir -p logs && sbatch run_all.sbatch

set -euo pipefail
export PANG="${PANG:-/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis}"
REV="$PANG/revision_2026-11"
cd "$REV"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate iqtree_env 2>/dev/null || true
export IQTREE="${IQTREE:-$(command -v iqtree2 || command -v iqtree)}"
echo "PANG=$PANG"
echo "IQTREE=$IQTREE"
python3 -c "import Bio,pandas,numpy,scipy;print('python deps ok')"

echo; echo "################ 17 - ASR + McDonald-Kreitman ################"
python3 17_ASR_MK/scripts/01_run_asr.py
python3 17_ASR_MK/scripts/02_mk_asr.py

echo; echo "################ 19 - provenance audit #######################"
python3 19_Provenance_Audit/scripts/01_provenance_audit.py

echo; echo "################ 18 - ancestral gene content #################"
python3 18_Gene_Content_ASR/scripts/01_gene_content_asr.py

echo; echo "################ 21 - core audit + branch support ############"
python3 21_Core_Audit/scripts/01_excluded_core_table.py
python3 21_Core_Audit/scripts/02_branch_support.py

echo; echo "done. summaries:"
find "$REV" -path '*/03_summary/*' -name '*.csv' -newermt '-1 day' | sort
