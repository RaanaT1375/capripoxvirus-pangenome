#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=30
#SBATCH --mem=60GB
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --job-name=relax_lsdv
#SBATCH --time=24:00:00
#SBATCH --output=15_RELAX/logs/relax_%j.out
#SBATCH --error=15_RELAX/logs/relax_%j.err

set -u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hyphy_env

cd "$SLURM_SUBMIT_DIR"
ALN=10_Selection_Pressure/01_codon_alignments
TRE=15_RELAX/01_labeled_trees
OUT=15_RELAX/02_results
LOG=15_RELAX/logs/relax_run.log
: > "$LOG"

ok=0; fail=0
for f in "$ALN"/*_codon.fasta; do
    og=$(basename "$f" _codon.fasta)
    [[ -s "$OUT/${og}_RELAX.json" ]] && { echo "$og: پرش (موجود)" >> "$LOG"; continue; }
    hyphy relax --alignment "$f" --tree "$TRE/${og}_labeled.nwk" \
        --test Test --reference Reference --models Minimal \
        --output "$OUT/${og}_RELAX.json" > "15_RELAX/logs/${og}.log" 2>&1
    if [[ -s "$OUT/${og}_RELAX.json" ]]; then
        echo "$og: ✅ ($(stat -c%s "$OUT/${og}_RELAX.json") بایت)" >> "$LOG"; ok=$((ok+1))
    else
        echo "$og: ❌" >> "$LOG"; fail=$((fail+1))
    fi
done
echo "==================================================" >> "$LOG"
echo "موفق: $ok | ناموفق: $fail" >> "$LOG"
echo "🏁 RELAX finished" >> "$LOG"
