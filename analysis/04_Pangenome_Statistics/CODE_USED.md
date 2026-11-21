# کدهای استفاده‌شده — مرحلهٔ ۰۴: آمار پان‌ژنوم

## ۱) طبقه‌بندی core / soft-core / shell / cloud

```python
import pandas as pd

RESULTS = "03_OrthoFinder_Results/OrthoFinder_Results_v3/Results_Aug04"

df = pd.read_csv(f"{RESULTS}/Orthogroups/Orthogroups.GeneCount.tsv",
                 sep='\t').set_index('Orthogroup')
if 'Total' in df.columns:
    df = df.drop(columns=['Total'])

presence = (df > 0)
n_genomes = presence.shape[1]
presence_fraction = presence.sum(axis=1) / n_genomes

def classify(f):
    if f >= 0.99:  return "Core (>=99%)"
    elif f >= 0.95: return "Soft-core (95-99%)"
    elif f >= 0.15: return "Shell (15-95%)"
    else:           return "Cloud (<15%)"

categories = presence_fraction.apply(classify)
print(f"=== ترکیب پان‌ژنوم ({n_genomes} ژنوم، {len(df)} orthogroup) ===")
print(categories.value_counts())

summary = pd.DataFrame({
    'n_genomes_present': presence.sum(axis=1),
    'fraction': presence_fraction,
    'category': categories
})
summary.to_csv('04_Pangenome_Statistics/pangenome_categories_v3.tsv', sep='\t')
```

## ۲) استخراج ژن‌های تک‌کپی-جهانی

```python
# تک‌کپی-جهانی: دقیقاً یک نسخه در هر ژنوم
strict_single_copy = df[(df == 1).all(axis=1)]
print(f"orthogroup تک‌کپی-جهانی: {len(strict_single_copy)}")

strict_single_copy.index.to_series().to_csv(
    '04_Pangenome_Statistics/strict_single_copy_orthogroups_v3.txt',
    index=False, header=False)
```

## ۳) منحنی Rarefaction و برازش قانون Heaps

```python
import numpy as np
from scipy.optimize import curve_fit

genomes = presence.columns.tolist()
n_genomes = len(genomes)

rng = np.random.default_rng(42)   # seed ثابت برای تکرارپذیری
records = []
for N in range(1, n_genomes + 1):
    for _ in range(100):          # 100 تکرار تصادفی به‌ازای هر N
        sampled = rng.choice(genomes, size=N, replace=False)
        sub = presence[sampled]
        records.append({
            "N": N,
            "pan_size":  int(sub.any(axis=1).sum()),   # اتحاد: حضور در ≥1
            "core_size": int(sub.all(axis=1).sum())    # اشتراک: حضور در همه
        })

curve_df = pd.DataFrame(records)
rare = curve_df.groupby("N").agg(
    pan_mean=("pan_size", "mean"),
    core_mean=("core_size", "mean")
).reset_index()
rare.to_csv('04_Pangenome_Statistics/rarefaction_summary_v3.tsv',
            sep='\t', index=False)

# برازش قانون Heaps: pan(N) = k * N^gamma
def heaps(N, k, gamma):
    return k * np.power(N, gamma)

popt, _ = curve_fit(heaps, rare["N"], rare["pan_mean"],
                    p0=[rare["pan_mean"].iloc[0], 0.5], maxfev=10000)
k_fit, gamma_fit = popt
openness = "OPEN" if gamma_fit > 0.05 else "CLOSED"

print(f"k = {k_fit:.3f}")
print(f"gamma = {gamma_fit:.4f}  -> {openness}")
print(f"pan-genome at N={n_genomes}: {rare['pan_mean'].iloc[-1]:.1f}")
print(f"strict core at N={n_genomes}: {rare['core_mean'].iloc[-1]:.1f}")
```

## ۴) بررسی outlier نهایی

```python
# تعریف نزدیک-هسته: حضور در ≥97% ژنوم‌ها
near_core = presence[presence.sum(axis=1) >= (0.97 * n_genomes)]
n_missing = (~near_core).sum(axis=0).sort_values(ascending=False)

print(f"orthogroup نزدیک-هسته: {len(near_core)}")
print(n_missing.describe(percentiles=[.5, .9, .95, .99]))
print("\n۱۰ ژنوم با بیشترین غیبت:")
print(n_missing.head(10))
```
