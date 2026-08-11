# Pan-genome and selection analysis of 290 Capripoxvirus genomes

Reproducible pipeline, intermediate results and documentation for a
comparative genomic analysis of the genus *Capripoxvirus*: **240 lumpy skin
disease virus (LSDV), 34 sheeppox virus (SPPV) and 16 goatpox virus (GTPV)**
genomes assembled from GenBank and SRA.

> **Read [`ANALYSIS_STATUS.md`](ANALYSIS_STATUS.md) before citing any number
> from this repository.** Several analyses were superseded during the project
> and each stage carries an explicit verdict.

## Headline results

**1. Present-day LSDV is evolutionarily static.** Across 240 genomes the
median number of distinct protein sequences per core gene is **6**; 27 of 60
genes have fewer than five. Purifying selection dominates (median
ω = 0.105) and there is no temporal signal in root-to-tip distance.

**2. Adaptive divergence happened deep in the past.** A McDonald–Kreitman
test polarized to the LSDV branch gives **α = 0.604** (CMH p = 0.009,
MH odds ratio 1.98), i.e. roughly 60% of amino-acid substitutions fixed on
that branch were driven by positive selection — during LSDV's divergence
from the SPPV/GTPV ancestor, not during its recent clonal expansion. The
estimate survives five independent robustness checks (two outgroups,
vaccine-strain removal, SRA-genome removal, polarization, gap masking) and
is a lower bound.

**3. Codon substitution models break down on this dataset.** This is
documented quantitatively rather than assumed — see below.

## Why codon models fail here — and how we caught it

Three independent failures, each quantified:

| Method | Failure mode | Evidence |
|---|---|---|
| BUSTED (240 LSDV) | dS collapses to zero | `prop_dS_near_zero` is 0.88–0.96 for **all 11** genes reaching q<0.05; LRT = 0 for 30 of 56 testable genes |
| RELAX | rate distribution collapses | 24 of 59 genes (41%) fail to converge; the survivors are not a random subset |
| aBSREL | mistakes an indel block for selection | see below |

The aBSREL case is the instructive one. `OG0000120` (Ig domain OX-2-like)
looked compelling from several angles: under the plain nucleotide GTR model
its LSDV stem branch is the longest of 300 branches in the tree, holding
43.7% of total tree length against a median of 0.082 across genes — and it
inverts the genome-wide pattern, where the LSDV stem is the *shortest* of
the three species stems. aBSREL returned a positively selected rate class at
ω = 697.6, LRT = 103.5, corrected p = 0, while the SPPV and GTPV stems gave
a single rate class and LRT = 0.

It is nonetheless an artefact. All 14 fixed differences fall in codons
174–202 of 217 (KS D = 0.802, p < 1e-9). The weight of the "positively
selected" class (0.141) matches that block's share of the gene
(29/217 = 0.134) almost exactly — the rate class *is* the block. And the
block sits on the gappiest part of the alignment (0% gaps in codons 20–160,
36% in 200–217), where LSDV carries a seven-residue indel and a divergent
hydrophobic C-terminal tail. dN/dS cannot be estimated across an indel-rich
region.

The clustering test that exposed this
(`16_aBSREL_stem/scripts/divergence_clustering.py`) also shows the problem is
gene-specific: of the nine genes with at least five fixed differences, only
this one is significantly clustered after FDR correction. The MK result is
therefore unaffected.

## Repository layout

Each stage directory contains `RESULTS_SUMMARY.txt` (key numbers, in Persian),
`CODE_USED.md` (commands and design decisions) and, from stage 09 onward,
`STATUS.md` (an English verdict: FINAL, SUPERSEDED, NEGATIVE or INCONCLUSIVE).

| Stage | Content |
|---|---|
| `01`–`04` | QC and filtering (323 → 290), OrthoFinder input, 206 orthogroups, pan-genome statistics |
| `05_Phylogeny` | ML tree from 60 single-copy genes (290 × 16,215, 1000 bootstraps) |
| `06_Recombination` | Parsnp + Gubbins (r/m = 0.13, 29.2% of the genome) |
| `07_Phylogeography` | PastML ancestral state reconstruction |
| `08_Population_Structure` | fastBAPS (two of five clusters are species, not LSDV lineages) |
| `09_Temporal_Signal` | root-to-tip regression — no temporal signal |
| `10_Selection_Pressure` | BUSTED on 290 genomes; **genus-level, not within-LSDV** |
| `11_Scoary` | gene–trait association; most hits are phylogenetic confounding |
| `12_MEME` | site-level selection — negative |
| `13_Selection_LSDV` | BUSTED on 240 LSDV genomes — negative + diagnostic |
| `14_MK_Test` | **McDonald–Kreitman — the primary result** |
| `15_RELAX` | selection intensity — inconclusive (41% non-convergence) |
| `16_aBSREL_stem` | episodic selection on the LSDV stem — negative |

`metadata/` holds species assignments, genome lists and length statistics.

## What is and is not included

Included: all analysis scripts, summary tables, phylogenetic trees, codon
alignments (stages 10 and 13), labelled trees, figures and documentation.

Not included, for size reasons: raw genome assemblies and protein FASTA
files, raw OrthoFinder and Gubbins output, raw HyPhy JSON results, and
cluster logs. All of these are regenerable from the scripts and the accession
lists in `metadata/`. Genome accessions are available from NCBI GenBank and
SRA.

## Reproducing an analysis

Scripts are run from the project root. Two conda environments are used —
`hyphy_env` for HyPhy analyses, `orthofinder_env` for everything else.

```bash
python3 13_Selection_LSDV/scripts/parse_and_diagnose.py
python3 14_MK_Test/scripts/mk_polarized.py
python3 16_aBSREL_stem/scripts/divergence_clustering.py
```

Two practical notes that cost time during this project:

- Four genes (`OG0000070`, `075`, `091`, `121`) fail in every HyPhy analysis.
  The fix is `hyphy <analysis> ENV=TOLERATE_NUMERICAL_ERRORS=1 ...` **as a
  command-line argument** — the exported environment variable does not work.
- In SLURM scripts, `--output` must not point into a directory that does not
  yet exist; SLURM opens the file before the script runs.

## Software versions

HyPhy 2.5.101, OrthoFinder 3.1.5, IQ-TREE 3.0.1, MAFFT 7.526, pal2nal 14.1,
Parsnp 2.1.5, Gubbins 3.4.3, PastML 1.9.50, fastBAPS 1.0.8, Scoary 1.6.16,
Prokka 1.15.6, Python 3.14 (Biopython 1.87, pandas 3.0.3, NumPy 2.5.1,
SciPy 1.18.0).

## License

Source code (`scripts/`) is released under the MIT License (`LICENSE`).
Data, result tables, alignments, trees and documentation are released under
CC BY 4.0 (`LICENSE-DATA`). Underlying sequences remain subject to the terms
of NCBI GenBank and SRA.

## Citation

See `CITATION.cff`, or use the GitHub "Cite this repository" button.
