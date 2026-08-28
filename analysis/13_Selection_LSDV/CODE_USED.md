# مرحلهٔ ۱۳ — BUSTED روی ۲۴۰ ژنوم LSDV

محیط: `hyphy_env` برای اجرا، `orthofinder_env` برای پارس | اجرا از ریشهٔ پروژه

## ۱) ساخت زیرمجموعهٔ LSDV
الاینمنت‌های ۲۹۰ ژنومی مرحلهٔ ۱۰ به ۲۴۰ ژنوم LSDV هرس شدند
(فهرست از `00_Metadata/lsdv_only_genomes.txt`)، و درخت‌ها متناظراً.

- خروجی: `01_codon_alignments/OG*_codon.fasta` (۲۴۰ توالی)
- خروجی: `02_pruned_trees/OG*_pruned.nwk`
- خروجی: `diversity_LSDV_only.csv` — شمارش تنوع به‌ازای ژن:
  `uniq_aa`، `poly_nt_sites`، `poly_aa_sites`

## ۲) اجرای BUSTED
```bash
sbatch 13_Selection_LSDV/scripts/run_busted_lsdv.sbatch
```
۶۰ از ۶۰ ژن موفق. (برخلاف مرحلهٔ ۱۰ که ۴ ژن نیاز به
`ENV=TOLERATE_NUMERICAL_ERRORS=1` داشتند.)

## ۳) پارس و تشخیص
```bash
conda activate orthofinder_env
python3 13_Selection_LSDV/scripts/parse_and_diagnose.py
```

این اسکریپت علاوه بر جدول نتایج، چهار تشخیصی را می‌سازد که ۱۱ ژن
معنادار را رد می‌کند:

- `prop_dS_near0` — سهم سایت‌هایی که نرخ مترادفشان زیر ۰.۰۵ است،
  از `fits → Unconstrained model → Rate Distributions →
  Synonymous site-to-site rates`
- `dS_min` و `dS_max` — دامنهٔ نرخ مترادف
- `LRT == 0` — شمارش ژن‌های غیرقابل آزمون
- پیوست `uniq_aa` از `diversity_LSDV_only.csv`

خروجی: `04_summary_tables/busted_lsdv_final.csv`

## اسکریپت‌ها
| فایل | کار |
|---|---|
| `scripts/run_busted_lsdv.sbatch` | اجرای BUSTED روی ۶۰ ژن |
| `scripts/parse_and_diagnose.py` | پارس + چهار تشخیصی dS |
