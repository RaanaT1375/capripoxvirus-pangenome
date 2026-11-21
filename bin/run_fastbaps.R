#!/usr/bin/env Rscript
# Bayesian hierarchical clustering of recombination-free SNPs.
# The number of clusters is not fixed a priori; it follows from the
# partitioning of the SNP matrix itself. Prior sensitivity is reported so that
# the stability of the partition can be judged rather than assumed.
#
# Portable version of bin/legacy/run_fastbaps.R
suppressPackageStartupMessages({ library(fastbaps); library(ape) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: run_fastbaps.R <snp_alignment.fasta> <clusters.csv> <priors.tsv>")

aln <- args[1]; out_clusters <- args[2]; out_priors <- args[3]

sparse <- import_fasta_sparse_nt(aln)
cat("genomes:", ncol(sparse$snp.matrix), " SNP sites:", nrow(sparse$snp.matrix), "\n")

sparse <- optimise_prior(sparse, type = "optimise.symmetric")
hc   <- fast_baps(sparse)
best <- best_baps_partition(sparse, hc)

write.csv(data.frame(genome_id = colnames(sparse$snp.matrix), cluster = best),
          out_clusters, row.names = FALSE)

# prior sensitivity: the partition should not depend on the prior chosen
res <- do.call(rbind, lapply(c("symmetric", "baps", "optimise.baps", "optimise.symmetric"),
  function(p) {
    s <- tryCatch(optimise_prior(import_fasta_sparse_nt(aln), type = p), error = function(e) NULL)
    if (is.null(s)) return(data.frame(prior = p, n_clusters = NA))
    k <- length(unique(best_baps_partition(s, fast_baps(s))))
    data.frame(prior = p, n_clusters = k)
  }))
write.table(res, out_priors, sep = "\t", row.names = FALSE, quote = FALSE)
print(res)
