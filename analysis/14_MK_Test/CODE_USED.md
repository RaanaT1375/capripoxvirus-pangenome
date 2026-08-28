# مرحلهٔ ۱۴ — آزمون McDonald–Kreitman

محیط: `orthofinder_env` | اجرا از ریشهٔ پروژه

## ورودی
- `10_Selection_Pressure/01_codon_alignments/OG*_codon.fasta` (۶۰ ژن، ۲۹۰ توالی)
- `00_Metadata/species_assignment.csv` — ستون‌های `Name` و `species`
- `11_Scoary/01_inputs/traits.csv` — ستون `Vaccine` برای حذف ۲۰ سویهٔ واکسن

## ترتیب اجرا

```bash
conda activate orthofinder_env

# ۱) نسخهٔ پایه — دو برون‌گروه، جدول ۲×۲، فیشر
python3 14_MK_Test/scripts/mk_test.py

# ۲) نسخهٔ دوم — آستانهٔ فراوانی، CMH طبقه‌بندی‌شده، آزمون علامت DoS
python3 14_MK_Test/scripts/mk_test_v2.py

# ۳) حساسیت — حذف واکسن و حذف ژنوم‌های SRA
python3 14_MK_Test/scripts/mk_sensitivity.py

# ۴) نتیجهٔ اصلی — قطبی‌شده به شاخهٔ LSDV
python3 14_MK_Test/scripts/mk_polarized.py

# ۵) آزمون استحکام پنجم — ماسک ستون‌های گپ‌دار
python3 14_MK_Test/scripts/mk_gapmask_sensitivity.py
```

`mk_test_v2.py` توابع مشترک را صادر می‌کند (`tr`، `path_counts`، `clean`،
`cmh`، `scan_gene`، `counts_at`) و بقیهٔ اسکریپت‌ها با
`sys.path.insert(0, "14_MK_Test/scripts")` از آن import می‌کنند.

## تصمیم‌های روش‌شناختی پیاده‌شده در کد

| تصمیم | جای پیاده‌سازی |
|---|---|
| دو برون‌گروه جداگانه (SPPV n=34، GTPV n=16) | حلقهٔ `for ogsp in [...]` |
| حالت اجدادی = کدون اکثریت برون‌گروه، اجماع ≥۹۰٪ | `cons()` / `out_cons` |
| سایت چندشکل فقط P شمرده شود (محافظه‌کارانه) | شاخهٔ `if len(icc) > 1` |
| مسیریابی Nei–Gojobori برای کدون‌های چنداختلافی | `path_counts()` |
| میانگین‌گیری روی مسیرهای بدون کدون استاپ | فیلتر `tr(nxt) == "*"` |
| حذف واریانت‌های زیر ۵٪ | `counts_at(..., cutoff)` |
| حذف ۲۰ سویهٔ واکسن | `ing_ids` |
| CMH طبقه‌بندی‌شده به‌ازای ژن (نه جمع ساده) | `cmh()` |
| قطبی‌سازی: اجدادی = توافق SPPV و GTPV | `polarized_gene()` |

## قطبی‌سازی
معتبر است چون کلاد SPPV+GTPV در `05_Phylogeny/04_tree/core_genome_ML_tree_v3.contree`
با bootstrap=100 تأیید شد. سایتی که دو برون‌گروه بر آن توافق ندارند کنار
گذاشته می‌شود؛ ۱۴٬۴۱۲ از ۱۶٬۲۱۵ کدون (۸۸.۹٪) باقی می‌ماند.

## خروجی‌ها
| فایل | محتوا |
|---|---|
| `02_results/mk_vs_SPPV.csv`، `mk_vs_GTPV.csv` | نسخهٔ پایه |
| `02_results/mk_v2_vs_SPPV.csv`، `mk_v2_vs_GTPV.csv` | با آستانهٔ ۵٪ و DoS |
| `02_results/mk_polarized_LSDV_branch.csv` | **نتیجهٔ اصلی** |
| `01_inputs/alignment_species_composition.csv` | ترکیب گونه‌ای هر الاینمنت |

## اشکالات شناخته‌شده
- `Series.clip` در این نسخهٔ pandas پارامتر `upper=` می‌گیرد نه `max=`
  (نسخهٔ اول `mk_test_v2.py` به همین دلیل وسط اجرا متوقف شد).
