# کدهای استفاده‌شده — مرحلهٔ ۰۸: خوشه‌بندی جمعیتی

## ۱) نصب fastBAPS (از bioconda، باینری از‌پیش‌کامپایل‌شده)

```bash
conda create -n fastbaps_env \
    -c conda-forge -c bioconda -c defaults \
    r-fastbaps r-ape -y
conda activate fastbaps_env
R -e 'library(fastbaps); packageVersion("fastbaps")'
# نسخه: 1.0.8
```

## ۲) آماده‌سازی الاینمنت ورودی

```python
from Bio import SeqIO
import re

# مهم: از الاینمنت اصلاح‌شدهٔ Gubbins (نه Parsnp خام)
# نواحی نوترکیبی در این فایل حذف شده‌اند
records = []
for rec in SeqIO.parse(
    "../06_Recombination/02_Gubbins_Results/lsdv_gubbins.filtered_polymorphic_sites.fasta",
    "fasta"):
    name = re.sub(r"\.fasta$", "", rec.id)
    if "ref" in name.lower():
        continue
    rec.id = name
    rec.description = ""
    records.append(rec)

SeqIO.write(records, "filtered_polymorphic_sites_clean.fasta", "fasta")
print(f"تعداد توالی: {len(records)}")
# نتیجه: 290
```

## ۳) اجرای اصلی fastBAPS

```r
library(fastbaps)
library(ape)

sparse.data <- import_fasta_sparse_nt("filtered_polymorphic_sites_clean.fasta")
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")

cat("تعداد ژنوم:", ncol(sparse.data$snp.matrix), "\n")
cat("تعداد موقعیت SNP:", nrow(sparse.data$snp.matrix), "\n")

baps.hc <- fast_baps(sparse.data)
best.partition <- best_baps_partition(sparse.data, baps.hc)

df <- data.frame(
    genome_id = colnames(sparse.data$snp.matrix),
    cluster = best.partition
)
write.csv(df, "fastbaps_clusters.csv", row.names = FALSE)
print(table(best.partition))
```

## ۴) تحلیل حساسیت با روش‌های مختلف Prior

```r
library(fastbaps)

sparse.data <- import_fasta_sparse_nt("filtered_polymorphic_sites_clean.fasta")

for (prior_type in c("symmetric", "baps", "optimise.symmetric", "optimise.baps")) {
    tryCatch({
        sd.opt <- optimise_prior(sparse.data, type = prior_type)
        baps.hc <- fast_baps(sd.opt)
        best.part <- best_baps_partition(sd.opt, baps.hc)
        n_clusters <- length(unique(best.part))
        cat(sprintf("%s → %d خوشه | اندازه‌ها: %s\n",
            prior_type, n_clusters,
            paste(as.vector(table(best.part)), collapse=",")))
    }, error = function(e) {
        cat(sprintf("%s → خطا: %s\n", prior_type, e$message))
    })
}
```

## ۵) تطبیق با متادیتا و نوترکیبی

```python
import csv
from collections import defaultdict, Counter

clusters = {}
with open("fastbaps_clusters.csv") as fh:
    for row in csv.DictReader(fh):
        clusters[row["genome_id"]] = int(row["cluster"])

geo = {}
with open("../07_Phylogeography/geo_metadata.csv") as fh:
    for row in csv.DictReader(fh):
        geo[row["tip"]] = {"country": row["Country"], "continent": row["Continent"]}

recomb = {}
color_labels = {
    "#8B5E3C": "Ancestral_recombinant",
    "#A63446": "Independent_recombinant",
    "#D8D4C8": "No_major_recombination",
}
with open("../06_Recombination/03_iTOL_Recombination/itol_recombination_status.txt") as fh:
    in_data = False
    for line in fh:
        line = line.rstrip("\n")
        if line.strip() == "DATA":
            in_data = True; continue
        if not in_data or not line.strip(): continue
        tip, color = line.split("\t")
        recomb[tip] = color_labels.get(color, "Unknown")

for c in sorted(set(clusters.values())):
    members = [g for g, cl in clusters.items() if cl == c]
    continents = Counter(geo.get(g, {}).get("continent", "?") for g in members)
    recomb_status = Counter(recomb.get(g, "?") for g in members)
    countries = Counter(geo.get(g, {}).get("country", "?") for g in members)
    top3 = [f"{k}({v})" for k, v in countries.most_common(3)]
    print(f"خوشهٔ {c} ({len(members)} ژنوم)")
    print(f"  قاره: {dict(continents)}")
    print(f"  نوترکیبی: {dict(recomb_status)}")
    print(f"  کشور (۳ تا): {', '.join(top3)}")

# ذخیرهٔ جدول کامل
with open("fastbaps_annotated.csv", "w") as fh:
    fh.write("genome_id,cluster,continent,country,recombination_status\n")
    for g, c in clusters.items():
        continent = geo.get(g, {}).get("continent", "")
        country = geo.get(g, {}).get("country", "")
        r = recomb.get(g, "")
        fh.write(f"{g},{c},{continent},{country},{r}\n")
```
