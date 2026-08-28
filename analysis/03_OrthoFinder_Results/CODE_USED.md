# کدهای استفاده‌شده — مرحلهٔ ۰۳: OrthoFinder

## اسکریپت SLURM

```bash
#!/bin/bash
#SBATCH --job-name=OrthoFinder_LSDV_v3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=xpkk@zhaw.ch
#SBATCH --cpus-per-task=30
#SBATCH --mem=120GB
#SBATCH --time=12:00:00
#SBATCH --partition=earth-3
#SBATCH --constraint=rhel8
#SBATCH --output=OrthoFinder_v3_run.log

module load USS/2022
module load slurm
source ~/.bashrc
conda activate orthofinder_env

cd /cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis

orthofinder \
    -f 02_OrthoFinder_Input/OrthoFinder_Input_v2 \
    -t 30 \
    -a 30 \
    -S diamond \
    -og \
    -o 03_OrthoFinder_Results/OrthoFinder_Results_v3
```

ارسال:
```bash
sbatch 03_OrthoFinder_Results/run_orthofinder_v3.sbatch
squeue -u xpkk
```

## تأیید خروجی

```bash
RESULTS_DIR="OrthoFinder_Results_v3/Results_Aug04"

# تعداد ستون‌ها (باید 292 باشد: Orthogroup + 290 ژنوم + Total)
head -1 "$RESULTS_DIR/Orthogroups/Orthogroups.GeneCount.tsv" | tr '\t' '\n' | wc -l

# تعداد orthogroup
tail -n +2 "$RESULTS_DIR/Orthogroups/Orthogroups.GeneCount.tsv" | wc -l

# درصد تخصیص ژن‌ها
total=0
for f in ../02_OrthoFinder_Input/OrthoFinder_Input_v2/*.faa; do
    n=$(grep -c '^>' "$f")
    total=$((total + n))
done
unassigned=$(($(wc -l < "$RESULTS_DIR/Orthogroups/Orthogroups_UnassignedGenes.tsv") - 1))

python3 -c "print(f'کل={$total}  تخصیص‌نیافته={$unassigned}  درصد={(1-$unassigned/$total)*100:.2f}%')"
```

## تأیید وضعیت Job

```bash
sacct -j <JOB_ID> --format=JobID,JobName,State,ExitCode,Elapsed
grep -iE "error|fail|traceback" OrthoFinder_v3_run.log
```
