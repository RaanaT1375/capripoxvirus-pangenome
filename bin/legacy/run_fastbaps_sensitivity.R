#!/usr/bin/env Rscript
library(fastbaps)
library(ape)

sparse.data <- import_fasta_sparse_nt("filtered_polymorphic_sites_clean.fasta")

cat("=== مقایسه روش‌های مختلف Prior ===\n")
for (prior_type in c("symmetric", "baps", "optimise.symmetric", "optimise.baps")) {
    tryCatch({
        sd.opt <- optimise_prior(sparse.data, type = prior_type)
        baps.hc <- fast_baps(sd.opt)
        best.part <- best_baps_partition(sd.opt, baps.hc)
        n_clusters <- length(unique(best.part))
        cat(sprintf("  %s  →  %d خوشه  |  توزیع: %s\n",
            prior_type, n_clusters,
            paste(names(table(best.part)), collapse=",")))
        cat(sprintf("           اندازه‌ها: %s\n",
            paste(as.vector(table(best.part)), collapse=",")))
    }, error = function(e) {
        cat(sprintf("  %s  →  خطا: %s\n", prior_type, e$message))
    })
}

cat("\n=== ساخت درخت برای شرط‌گذاری (اختیاری) ===\n")
tree <- read.tree("../06_Recombination/02_Gubbins_Results/lsdv_gubbins.final_tree_renamed.nwk")
sparse.data.opt <- optimise_prior(sparse.data, type = "optimise.symmetric")
baps.hc.tree <- fast_baps(sparse.data.opt)
best.partition.tree <- best_baps_partition(sparse.data.opt, baps.hc.tree)
cat("نتیجه بدون شرط‌گذاری (مرجع):", length(unique(best.partition.tree)), "خوشه\n")
print(table(best.partition.tree))
