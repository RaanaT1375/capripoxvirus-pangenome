# کدهای استفاده‌شده — مرحلهٔ ۰۲: آماده‌سازی ورودی

## ساخت پوشهٔ ورودی نهایی

```bash
cd .../02_OrthoFinder_Input

mkdir -p OrthoFinder_Input_v2

# کپی ژنوم‌های سالم از نسخهٔ قبلی (۲۸۷ ژنوم)
cp OrthoFinder_Input_Final/*.faa OrthoFinder_Input_v2/

# حذف ۴ کایمرای باقی‌مانده
rm -f OrthoFinder_Input_v2/SRR26747378_SRR26747379_Merged.faa \
      OrthoFinder_Input_v2/SRR26747380_SRR26747381_Merged.faa \
      OrthoFinder_Input_v2/SRR26747382_SRR26747383_Merged.faa \
      OrthoFinder_Input_v2/SRR26747384_SRR26747385_Merged.faa

# افزودن ۶ نمونهٔ اصلاح‌شده (با حذف پسوند _CORRECTED)
CORRECTED_DIR=".../Merged_Amplicons_CORRECTED/FAA_Results"
for f in "$CORRECTED_DIR"/*.faa; do
    base=$(basename "$f" "_CORRECTED.faa")
    cp "$f" "OrthoFinder_Input_v2/${base}.faa"
done

# افزودن نمونهٔ جاافتاده
cp .../SRR26792425/Prokka_Output/SRR26792425.faa OrthoFinder_Input_v2/

echo "تعداد نهایی: $(ls OrthoFinder_Input_v2/*.faa | wc -l)"
```

## بررسی سلامت دیتاست

```bash
# تعداد کل ژن‌های ورودی
total=0
for f in OrthoFinder_Input_v2/*.faa; do
    n=$(grep -c '^>' "$f")
    total=$((total + n))
done
echo "کل ژن‌ها: $total"

# توزیع تعداد پروتئین
for f in OrthoFinder_Input_v2/*.faa; do
    echo "$(grep -c '^>' "$f")"
done | sort -n | awk '{
    a[NR]=$1
} END {
    print "کمینه:", a[1]
    print "میانه:", a[int(NR/2)]
    print "بیشینه:", a[NR]
}'
```
