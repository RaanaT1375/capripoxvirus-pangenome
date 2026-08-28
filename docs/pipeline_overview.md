# Pipeline overview

```mermaid
flowchart TB
    subgraph P["01-04  reads to annotated assemblies"]
        A1[fastp] --> A2[Bowtie2 host depletion]
        A2 --> A3["SPAdes --careful"]
        A3 --> A4[Prokka, reference-guided]
        A4 --> A5{"CDS > 200?"}
        A5 -- yes --> A6[minimap2 bacterial depletion]
        A6 --> A7["SPAdes --only-assembler"]
        A7 --> A8[Prokka]
        A5 -- no --> B0
        A8 --> B0
    end

    subgraph Q["05-06  QC and pangenome"]
        B0[annotated proteomes] --> B1[DIAMOND contamination screen]
        B1 --> B2[completeness filter]
        B2 --> B3[OrthoFinder]
        B3 --> B4[core / soft-core / shell / cloud]
        B3 --> B5["Heaps' law fit"]
    end

    subgraph R["07-08  phylogenomics"]
        B4 --> C1[MAFFT per orthogroup]
        C1 --> C2[concatenated supermatrix]
        C2 --> C3[IQ-TREE, partitioned]
    end

    subgraph S["09-13  population processes"]
        B2 --> D1[Parsnp core alignment]
        D1 --> D2[Gubbins]
        D2 --> D3[fastBAPS]
        D2 --> D4[PastML, full and sub-sampled]
        D2 --> D5[root-to-tip regression]
        C3 --> D6[Scoary, tree-corrected]
        B3 --> D6
    end

    subgraph T["14-16  selection"]
        C1 --> E1[pal2nal codon alignments]
        E1 --> E2[BUSTED / RELAX / aBSREL]
        E2 --> E3["dN/dS applicability diagnostics"]
        E1 --> E4[polarised McDonald-Kreitman]
        E4 --> E5[frequency cutoff sweep]
    end

    E3 -.->|"dS collapse: codon models uninterpretable"| E4
```

## Why the diagnostic arrow matters

The dashed edge is the analytical hinge of the study. The codon-model scans are
not run because they are expected to work; they are run to establish, with
recorded evidence, that they cannot work on a clonal population with almost no
synonymous variation. `SELECTION_DIAGNOSTICS` writes that evidence out per gene:
the number of distinct protein haplotypes, the proportion of sites with dS
approximately zero, and whether the likelihood-ratio statistic degenerated.
Only then does the population-genetic test become the primary result rather than
an alternative one.

## Stage inputs and outputs

| Stage | Key input | Key output | Tool |
|---|---|---|---|
| 01 | raw FASTQ | trimmed reads | fastp |
| 02a | trimmed reads | host-depleted reads | Bowtie2 |
| 02b | flagged reads | bacterially depleted reads | minimap2, samtools |
| 03 | depleted reads | filtered contigs | SPAdes |
| 04 | contigs | proteome, CDS, GFF | Prokka |
| 05 | proteomes | retained genome list | DIAMOND |
| 06 | retained proteomes | orthogroups, presence/absence | OrthoFinder |
| 07 | single-copy orthogroups | supermatrix, partitions | MAFFT |
| 08 | supermatrix | ML tree | IQ-TREE |
| 09 | genomes | core alignment, recombination blocks, r/m | Parsnp, Gubbins |
| 10 | recombination-free SNPs | genetic clusters | fastBAPS |
| 11 | filtered tree, metadata | ancestral states, sub-sampling control | PastML |
| 12 | filtered tree, dates | root-to-tip regression | Biopython |
| 13 | presence/absence, traits, tree | gene-trait associations | Scoary |
| 14 | codon alignments, tree | episodic and branch selection | HyPhy |
| 15 | codon alignments | polarised alpha, threshold sweep | in-house |

## Resource notes

Bacterial depletion maps against a ~100 GB minimap2 index and carries the
`process_high_memory` label. OrthoFinder and Gubbins are the two other stages
that dominate wall time on a full dataset. The HyPhy stages are numerous and
individually small, so they carry `errorStrategy = 'ignore'`: a per-gene
optimisation failure is a recorded diagnostic, not a pipeline failure, and the
diagnostics stage counts those failures explicitly.
