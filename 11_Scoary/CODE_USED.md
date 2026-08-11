# مرحلهٔ ۱۱ — کدهای استفاده‌شده (Scoary)

## نصب
```bash
conda create -n scoary_env -c conda-forge -c bioconda scoary -y
```

## اجرای Scoary (Job 645294)
نکته: پرچم permutation در این نسخه `-e 1000` است، نه `--permutations`.
```bash
#!/bin/bash
#SBATCH --job-name=Scoary_LSDV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=8
#SBATCH --mem=30GB
#SBATCH --time=06:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=scoary_run.log

module load USS/2022
module load slurm
source ~/.bashrc
conda activate scoary_env

cd /cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis/11_Scoary

rm -rf scoary_results

scoary \
    --genes gene_presence_absence.csv \
    --traits traits.csv \
    --outdir scoary_results \
    --start_col 15 \
    --no-time \
    -e 1000 \
    --threads 8

echo "کد خروج: $?"
echo "=================================================="
echo "🏁 Scoary finished"
echo "=================================================="
```
