#!/usr/bin/env Rscript
library(fastbaps)
library(ape)

sparse.data <- import_fasta_sparse_nt("filtered_polymorphic_sites_clean.fasta")
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")

cat("=== ساختار داده ===\n")
cat("تعداد ژنوم:", ncol(sparse.data$snp.matrix), "\n")
cat("تعداد موقعیت SNP:", nrow(sparse.data$snp.matrix), "\n")

baps.hc <- fast_baps(sparse.data)
best.partition <- best_baps_partition(sparse.data, baps.hc)

df <- data.frame(
    genome_id = colnames(sparse.data$snp.matrix),
    cluster = best.partition
)
write.csv(df, "fastbaps_clusters.csv", row.names = FALSE)

cat("\n=== تعداد خوشه‌های شناسایی‌شده ===\n")
print(table(best.partition))
