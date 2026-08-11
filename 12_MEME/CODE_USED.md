# مرحلهٔ ۱۲ — MEME

محیط: `hyphy_env` | اجرا از ریشهٔ پروژه

## ورودی
- `10_Selection_Pressure/01_codon_alignments/OG*_codon.fasta` (۲۹۰ توالی)
- `10_Selection_Pressure/02_pruned_trees/OG*_pruned.nwk`

## اجرا
پنج ژن با بالاترین omega در BUSTED مرحلهٔ ۱۰ انتخاب شدند:
OG0000120، OG0000123، OG0000110، OG0000116، OG0000103

```bash
for OG in OG0000120 OG0000123 OG0000110 OG0000116 OG0000103; do
    hyphy meme \
        --alignment 10_Selection_Pressure/01_codon_alignments/${OG}_codon.fasta \
        --tree      10_Selection_Pressure/02_pruned_trees/${OG}_pruned.nwk \
        --output    12_MEME/01_results/${OG}_MEME.json
done
```

## اسکریپت‌ها
| فایل | کار |
|---|---|
| `scripts/` | اجرای MEME و پارس نتایج |

## نکتهٔ تفسیری
هنگام پارس، ستون `alpha` (نرخ مترادف) هر سایت معنادار باید بررسی شود.
`alpha` نزدیک صفر یعنی سیگنال از فروریختن مخرج می‌آید نه از انتخاب —
همان تشخیصی که در مرحلهٔ ۱۳ کمّی شد. هر دو سایت معنادار این مرحله
`alpha ≈ 0` داشتند.
