# Reference files

These are not committed to the repository because of their size. Fetch them
before the first run:

```bash
# LSDV Neethling vaccine backbone, used for synteny and contamination screening
efetch -db nuccore -id AF325528.1 -format fasta > AF325528.1.fasta

# LSDV RefSeq proteome, used to guide Prokka annotation
datasets download genome accession GCF_000839805.1 --include protein
unzip -p ncbi_dataset.zip '*/protein.faa' > GCF_000839805.1_protein.faa
```

The bacterial minimap2 index is built separately from the NCBI RefSeq
representative bacterial genome collection:

```bash
minimap2 -t 32 -d Bacteria_Representative_DB.mmi Bacteria_Representative_DB.fna
```

Host Bowtie2 indices are built from the Ensembl top-level assemblies for
*Bos taurus*, *Ovis aries* and *Capra hircus*, and must be named
`Bos_taurus`, `Ovis_aries` and `Capra_hircus` inside the directory passed to
`--host_bowtie2_index`.
