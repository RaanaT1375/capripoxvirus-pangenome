# مرحلهٔ ۱۶ — aBSREL روی یال LSDV

محیط: `orthofinder_env` برای برچسب‌گذاری و پارس، `hyphy_env` برای اجرا
اجرا از ریشهٔ پروژه.

## منطق طراحی

برخلاف مرحلهٔ ۱۵ که کل کلاد LSDV را آزمون می‌کند، اینجا فقط **یک یال**
برچسب می‌خورد: یالی که LSDV را از کلاد SPPV+GTPV جدا می‌کند. این همان
یالی است که MK قطبی‌شده واگرایی را به آن نسبت می‌دهد، پس مرحلهٔ ۱۶ همان
فرضیه را با درست‌نمایی می‌آزماید به‌جای شمارش.

در درخت بدون ریشه، این یال با پیدا کردن گرهی که دقیقاً برگ‌های
SPPV و GTPV را دارد مشخص می‌شود.

## ۱) برچسب‌گذاری

    conda activate orthofinder_env
    python3 16_aBSREL_stem/scripts/label_stem_trees.py

خروجی در `01_labeled_trees/OG*_stem.nwk`. هر فایل باید دقیقاً یک
`{Stem}` داشته باشد:

    grep -o '{Stem}' 16_aBSREL_stem/01_labeled_trees/OG0000120_stem.nwk | wc -l

کنترل اختصاصیت با همان اسکریپت ساخته می‌شود:

    python3 16_aBSREL_stem/scripts/label_stem_trees.py SPPV
    python3 16_aBSREL_stem/scripts/label_stem_trees.py GTPV

که یالِ منتهی به کلاد SPPV یا GTPV را برچسب می‌زند و در `04_specificity`
ذخیره می‌کند.

## ۲) اجرا

aBSREL برای هر ژن ۳ تا ۱۳ دقیقه می‌گیرد (مدل پایه روی کل درخت ۲۹۰ برگی
برازش می‌شود حتی وقتی فقط یک شاخه آزمون می‌شود)، پس به‌صورت job array
اجرا شد نه متوالی:

    sbatch 16_aBSREL_stem/scripts/run_absrel.sbatch

آرایهٔ `1-60%12` با ۸ CPU و ۱۶GB به‌ازای هر تسک. Job 647986، حدود ۳۲ دقیقه.

فرمان اصلی:

    hyphy absrel --alignment <codon.fasta> --tree <stem.nwk> \
        --branches Stem --output <OG>_ABSREL.json

### چهار ژن مقاوم

OG0000070، OG0000075، OG0000091، OG0000121 در هر تحلیل HyPhy شکست
می‌خورند و JSON خالی می‌دهند. راه‌حل، **به‌صورت آرگومان خط فرمان**:

    hyphy absrel ENV=TOLERATE_NUMERICAL_ERRORS=1 --alignment ... --tree ...

متغیر محیطی export شده (`export TOLERATE_NUMERICAL_ERRORS=1`) کار نمی‌کند.
این نکته یک بار وقت تلف کرد.

## ۳) پارس

    conda activate orthofinder_env
    python3 16_aBSREL_stem/scripts/parse_absrel_and_context.py

علاوه بر p و omega، بستر مقایسه هم می‌سازد: طول همان یال زیر مدل خالص
نوکلئوتیدی GTR، سهمش از کل طول درخت، و رتبه‌اش بین ۳۰۰ شاخه. این ستون‌ها
مدل-مستقل‌اند و برای تشخیص آرتیفکت لازم بودند.

نکته: کلید `Global MG94xREV` در JSONهای این نسخهٔ HyPhy وجود ندارد؛
`Nucleotide GTR` استفاده می‌شود.

## ۴) آزمون خوشه‌ای بودن — اسکریپت تعیین‌کننده

    python3 16_aBSREL_stem/scripts/divergence_clustering.py

این تحلیلی است که ادعای انتخاب روی OG0000120 را باطل کرد. منطقش ساده
است: انتخاب متنوع‌کننده جایگزینی‌ها را در طول ژن پراکنده می‌کند، ولی
ورود یک قطعهٔ واگرا یا یک ناحیهٔ غیرقابل‌الاینمنت آن‌ها را در یک بلوک
جمع می‌کند. آزمون کولموگروف-اسمیرنوف در برابر توزیع یکنواخت این دو را
از هم جدا می‌کند.

بازرسی یک ژن خاص همراه با توزیع گپ:

    python3 16_aBSREL_stem/scripts/divergence_clustering.py --gaps OG0000120

که هیستوگرام موقعیت جایگزینی‌ها و درصد گپ در پنجره‌های ۲۰ کدونی را
می‌دهد — همان دو خروجی که نشان داد بلوک OG0000120 دقیقاً روی
گپ‌دارترین ناحیهٔ ژن افتاده است.

## اسکریپت‌ها

| فایل | کار |
|---|---|
| scripts/label_stem_trees.py | برچسب {Stem}؛ با آرگومان SPPV/GTPV کنترل می‌سازد |
| scripts/run_absrel.sbatch | job array روی ۶۰ ژن |
| scripts/rerun_failed.sbatch | اجرای مجدد ۴ ژن مقاوم |
| scripts/parse_absrel_and_context.py | پارس + بستر طول شاخه زیر GTR |
| scripts/divergence_clustering.py | آزمون KS + بازرسی گپ |

## خروجی‌ها

| فایل | محتوا |
|---|---|
| 03_summary/absrel_stem_final.csv | نتیجهٔ اصلی به‌ازای ژن |
| 03_summary/stem_branch_gtr.csv | همان، مرتب بر اساس سهم یال |
| 03_summary/divergence_clustering.csv | آزمون خوشه‌ای بودن |
| 04_specificity/ | کنترل یال SPPV و GTPV برای OG0000120 |

## ۴) آزمون خوشه‌ای بودن — اسکریپت تعیین‌کننده

    python3 16_aBSREL_stem/scripts/divergence_clustering.py

این تحلیلی است که ادعای انتخاب روی OG0000120 را باطل کرد. منطقش ساده
است: انتخاب متنوع‌کننده جایگزینی‌ها را در طول ژن پراکنده می‌کند، ولی
ورود یک قطعهٔ واگرا یا یک ناحیهٔ غیرقابل‌الاینمنت آن‌ها را در یک بلوک
جمع می‌کند. آزمون کولموگروف-اسمیرنوف در برابر توزیع یکنواخت این دو را
از هم جدا می‌کند.

بازرسی یک ژن خاص همراه با توزیع گپ:

    python3 16_aBSREL_stem/scripts/divergence_clustering.py --gaps OG0000120

که هیستوگرام موقعیت جایگزینی‌ها و درصد گپ در پنجره‌های ۲۰ کدونی را
می‌دهد — همان دو خروجی که نشان داد بلوک OG0000120 دقیقاً روی
گپ‌دارترین ناحیهٔ ژن افتاده است.

## اسکریپت‌ها

| فایل | کار |
|---|---|
| scripts/label_stem_trees.py | برچسب {Stem}؛ با آرگومان SPPV/GTPV کنترل می‌سازد |
| scripts/run_absrel.sbatch | job array روی ۶۰ ژن |
| scripts/rerun_failed.sbatch | اجرای مجدد ۴ ژن مقاوم |
| scripts/parse_absrel_and_context.py | پارس + بستر طول شاخه زیر GTR |
| scripts/divergence_clustering.py | آزمون KS + بازرسی گپ |

## خروجی‌ها

| فایل | محتوا |
|---|---|
| 03_summary/absrel_stem_final.csv | نتیجهٔ اصلی به‌ازای ژن |
| 03_summary/stem_branch_gtr.csv | همان، مرتب بر اساس سهم یال |
| 03_summary/divergence_clustering.csv | آزمون خوشه‌ای بودن |
| 04_specificity/ | کنترل یال SPPV و GTPV برای OG0000120 |
