# مرحلهٔ ۱۵ — RELAX

محیط: `orthofinder_env` برای برچسب‌گذاری و پارس، `hyphy_env` برای اجرا
اجرا از ریشهٔ پروژه.

## ۱) برچسب‌گذاری درخت‌ها

    conda activate orthofinder_env
    python3 15_RELAX/scripts/label_trees.py

هر برگ و هر گره داخلی برچسب می‌گیرد: `{Test}` برای کلاد LSDV و
`{Reference}` برای کلاد SPPV+GTPV. گره‌های مختلط بدون برچسب می‌مانند.

دو نکتهٔ پیاده‌سازی:

- ترکیب گونه‌ای گره‌های داخلی باید در یک پاس اول و با نام‌های اصلی برگ‌ها
  حساب شود. اگر اول برگ‌ها را تغییر نام دهید، جست‌وجو در
  species_assignment در پاس دوم شکست می‌خورد.
- `cl.confidence = None` لازم است، وگرنه مقدار bootstrap جای برچسب
  نوشته می‌شود.

بازبینی: در OG0000120_labeled.nwk باید {Test}=472 و {Reference}=99 باشد.

## ۲) اجرا

    sbatch 15_RELAX/scripts/run_relax.sbatch

فرمان اصلی داخل اسکریپت:

    hyphy relax --alignment <codon.fasta> --tree <labeled.nwk> \
        --test Test --reference Reference --models Minimal \
        --output <OG>_RELAX.json

Job 648048، زمان کل ۱:۲۷:۱۷.

اشکال SLURM که یک بار کار را کشت: مسیر `--output` در سربرگ نباید پوشهٔ
ناموجود داشته باشد. SLURM فایل خروجی را پیش از اجرای اسکریپت باز می‌کند،
پس `mkdir` داخل اسکریپت دیر است. Job 647985 به همین علت در یک ثانیه
FAILED شد.

## ۳) پارس

    conda activate orthofinder_env
    python3 15_RELAX/scripts/parse_relax.py

علاوه بر K و p، دو کنترل کیفیت هم می‌سازد:

- `identical_dists` — آیا توزیع نرخ Test و Reference عیناً یکی درآمده؟
  اگر بله، بهینه‌ساز دو مجموعه شاخه را از هم تفکیک نکرده. اینجا صفر بود.
- `degen` یعنی `omega_max > 100` — برازش‌های منفجرشده.

و آزمون می‌کند که آیا معناداری با دژنراسیون همبسته است (فیشر) و آیا
دژنراسیون K را بالا می‌برد (من-ویتنی).

## نرخ شکست: ۲۴ از ۵۹

خطای تکرارشونده هنگام برازش مدل K != 1 :

    relax._renormalize_with_weights  ->  Max(USie_UMg.mean, 0.001)
    Attempting to operate on an undefined value

این فروریختن شبکهٔ اولیهٔ توزیع نرخ است، نه کرش تصادفی. ژن‌هایی که همگرا
شدند دقیقاً آن‌هایی‌اند که توزیعشان نریخت، پس ۳۶ ژن پارس‌شده نمونهٔ
تصادفی نیستند و نتیجه به ژنوم تعمیم‌پذیر نیست.

با ENV=TOLERATE_NUMERICAL_ERRORS=1 می‌شد بعضی را نجات داد، ولی نتیجه‌ای
که از نادیده‌گرفتن خطای عددی به دست بیاید قابل استناد نیست — عمداً انجام
نشد.

## اسکریپت‌ها

| فایل | کار |
|---|---|
| scripts/label_trees.py | برچسب {Test} و {Reference} |
| scripts/run_relax.sbatch | اجرای RELAX روی ۶۰ ژن |
| scripts/parse_relax.py | پارس + کنترل کیفیت برازش |

خروجی: `03_summary/relax_summary.csv`
