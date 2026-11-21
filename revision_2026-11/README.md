# Revision analyses — November 2026

Analyses added in response to the Scientific Reports review of
*Pan-genomic Evaluation of Host Adaptation and Recombination-Mediated
Emergence in Capripoxviruses* (submission `ecb9a1ba`, decision
**Request Revision**, handling editor Xiuguang Mao).

These stages continue the numbering of the original pipeline (01–16) and
follow the same layout, so they drop into the existing repository without
disturbing anything already published.

```
revision_2026-11/
├── _common.py                 shared paths, IQ-TREE discovery, loaders
├── run_all.sbatch             stages 17, 18, 19, 21 in order
├── 17_ASR_MK/                 Reviewer 2 comment 3  — ancestral-state polarization
├── 18_Gene_Content_ASR/       Reviewer 2 comment 1  — gene gain/loss reconstruction
├── 19_Provenance_Audit/       Reviewer 2 comment 1  + editor (overstated conclusions)
├── 20_Presence_Validation/    protects the pan-GWAS result (DNA-level presence)
├── 21_Core_Audit/             Reviewer 2 comments 2 and 4
├── 22_Scoary_Clean/           pan-GWAS repeated on the decontaminated matrix
└── 23_Trait_Audit/            integrity of the trait matrix itself
```

Every stage uses the same five directories: `scripts/`, `01_inputs/`,
`02_results/`, `03_summary/`, `logs/`. Only `03_summary/` is meant to be read
by a human or cited in the manuscript; `02_results/` holds the full machine
output.

## Running

```bash
export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
cd "$PANG/revision_2026-11"
mkdir -p logs
sbatch run_all.sbatch
```

Stage 20 needs BLAST+ and runs in three steps:

```bash
cd "$PANG/revision_2026-11/20_Presence_Validation"
conda activate prokka_env                       # or any env with BLAST+
python3 scripts/01_build_db_and_queries.py
sbatch scripts/02_run_tblastn.sbatch
python3 scripts/03_summarise_validation.py      # after the job finishes
```

Requirements: IQ-TREE 2.x or 3.x on `PATH` (or `export IQTREE=/path/to/iqtree2`),
Python 3 with Biopython, pandas, numpy and scipy, and BLAST+ for stage 20.
Nothing in stages 01–16 is modified or re-run.

## What each stage answers

### 17 — McDonald–Kreitman with ML ancestral states

> *"consensus sequences from extant outgroup lineages do not necessarily
> represent the true ancestral codon states … substitutions, reversions, or
> parallel mutations occurring independently within the SPPV and GTPV
> lineages may result in incorrect polarization"* — Reviewer 2

Each per-gene tree is rooted on the SPPV+GTPV clade and IQ-TREE reconstructs
marginal ancestral states at every internal node. The ≥90 % consensus codon
is replaced by the reconstruction at `ANC_OUT`, the MRCA of the outgroup
clade; everything downstream — Nei–Gojobori pathway averaging, the 5 % minor
variant cutoff, the 90 % coverage filter, the per-gene CMH test — is the
stage-14 code unchanged, so the new α is directly comparable with the
published value. A second, stricter estimator counts divergence between
`ANC_OUT` and `ANC_LSDV`, localising substitutions to the LSDV stem edge.

Reference values from a validation run (α = 0.604 from stage 14 reproduced
exactly before any change):

| estimator | α | CMH *p* | confident codons |
|---|---|---|---|
| consensus ≥90 % (published) | 0.604 | 8.96e-03 | 88.9 % |
| ASR branch, PP ≥ 0.90 | **0.629** | 9.01e-03 | **98.8 %** |
| ASR branch, PP ≥ 0.50 | 0.658 | 6.66e-05 | 100 % |
| ASR branch, PP ≥ 0.95 | 0.682 | 5.01e-02 | 97.9 % |
| ASR branch, OG0000120 removed | 0.505 | 2.25e-03 | — |
| ASR stem (`ANC_OUT`→`ANC_LSDV`) | 0.457 | 3.36e-02 | — |

Gene bootstrap 95 % CI [0.25, 0.77], overlapping the published [0.33, 0.75].
The conclusion is unchanged under every variant and the fraction of sites
with a confident ancestral state rises from 88.9 % to 98.8 %.

### 18 — Ancestral gene-content reconstruction

> *"The presence of low-frequency orthogroups does not, in itself,
> demonstrate that these genes have been acquired transiently … the authors
> should either provide additional evolutionary evidence, such as ancestral
> gene-content reconstruction, or revise the statement"* — Reviewer 2

The 206-orthogroup presence/absence matrix is reconstructed as a binary
character set on the fixed core tree, and every branch is scored as a gain or
a loss, terminal or internal.

| class | gains | losses | % gains terminal |
|---|---|---|---|
| cloud | 228 | 11 | 89.9 |
| shell | 91 | 115 | 56.0 |
| soft-core | 4 | 39 | 25.0 |
| core | 0 | 23 | — |

The accessory compartment is dominated by gains, and those gains sit almost
entirely on terminal branches — the manuscript's statement now rests on
evidence rather than on the frequency distribution. The reconstruction is
repeated on the 237 GenBank genomes alone, where the pattern survives
(49 gains vs 11 losses, 75.5 % terminal).

### 19 — Sequence provenance as a confounder

Not requested by a reviewer, but it is the strongest available answer to the
editor's instruction that *"any overstated conclusions are rewritten and the
limitations of the work fully explained"*, and it is what makes the
Heaps' law claim defensible.

The dataset has two halves: 237 single-contig GenBank assemblies and 53
genomes assembled de novo here from SRA reads. The audit shows the effect is
confined to the low-frequency compartment, which is where an assembly or
annotation artefact would live and nowhere else:

- 30 of the 50 cloud orthogroups occur in **zero** GenBank genomes
- per-genome cloud load is 7.9× higher in SRA-derived genomes (3.58 vs 0.45,
  Mann–Whitney *p* = 2.4e-04), while shell load is identical (8.98 vs 9.19)
- γ falls from 0.0698 (all 290) to **0.0325** on GenBank genomes alone —
  below the manuscript's own γ > 0.05 threshold for an open pan-genome
- the three robust Africa-associated pan-GWAS hits (OG0000159, OG0000156,
  OG0000161, all hypothetical proteins) are present in zero GenBank genomes,
  and the Africa trait is itself enriched among SRA-derived genomes
  (28.3 % vs 15.6 %, *p* = 0.046)
- the Asia-associated hits run the opposite way (enriched in GenBank,
  absent from SRA) and are therefore not explained by this mechanism

Consequence for the manuscript: the *terminal, lineage-specific* character of
the accessory compartment stands, the *open pan-genome* classification does
not survive provenance control and should be withdrawn or heavily qualified,
and the Africa pan-GWAS signal should be withdrawn pending stage 20.

### 20 — DNA-level validation of candidate loci

OrthoFinder calls a gene absent when Prokka did not predict it, which
conflates true absence with annotation dropout in a fragmented draft. Each
candidate orthogroup is re-tested with `tblastn` directly against all 290
assemblies, independently of annotation. Candidates recovered in GenBank
genomes by DNA search survive and can be re-tested with the corrected
matrix; candidates still absent from every GenBank genome at the DNA level
are artefacts and must be withdrawn.

### 21 — Core-orthogroup audit and branch support

> *"the use of only 60 single-copy core genes may introduce a bias toward
> more conserved genomic components"* — Reviewer 2, comment 2

The exclusion is decomposed exactly:

| | count |
|---|---|
| core, present in ≥ 99 % of genomes | 132 |
| &nbsp;&nbsp;– present in all 290 but not single-copy everywhere | 51 |
| &nbsp;&nbsp;– not present in all 290 (each missing from 1–2 genomes) | 21 |
| &nbsp;&nbsp;= single-copy and universal, analysed | 60 |

The reviewer's concern is empirically correct and the Discussion must say so:
seven host-interaction or immunomodulatory orthogroups are among the 72
excluded, including kelch-like, two ankyrin-repeat proteins, IL-1
receptor-like, a putative host-range protein, a soluble interferon-γ receptor
and CD47-like. `03_summary/S_excluded_core_orthogroups.csv` is the
supplementary table that lets a reader verify this directly.

Comment 4 (support values on Figure 3) first needs an ambiguity resolved: the
Results text describes Figure 3A as the IQ-TREE supermatrix tree while the
caption describes the RAxML tree built inside Gubbins. Only the first carries
UFBoot values, which is why the figure has none. `02_branch_support.py`
compares the two topologies; if they agree, showing the IQ-TREE tree and
correcting the caption is sufficient and no rerun is needed.

### 23 — Trait-matrix audit

`traits.csv` drives every pan-GWAS association and had never been checked
against the rest of the dataset. It is structurally sound — 290 genomes, no
missing values, strictly binary, and the genome set agrees exactly with the
species assignment, the presence/absence matrix and the phylogeny tip labels.
Three trait groups are not mutually exclusive, and in each case the overlap is
a deliberate encoding choice that the Methods does not state:

- `Host_Wildlife` is a superset of `Host_Giraffe` (one genome, OR644284.1)
- `Recombinant` (10 designated strains) overlaps `Ancestral_recombinant` and
  `Independent_recombinant`; the clean partition is
  `Ancestral_recombinant` (100) + `Independent_recombinant` (1) +
  `No_recombination` (189) = 290
- the same 10 recombinant strains are neither `Vaccine` nor `Wild`

The last one has a real consequence. Scoary treats a 0 as a negative, not as
missing, so those 10 strains are scored as evidence *against* both the vaccine
and the wild-type trait. Both surviving Vaccine tier-1 hits therefore need a
sensitivity run with those genomes excluded. The audit also shows
`No_recombination` contains all 34 SPPV and all 16 GTPV genomes, so any
recombination association is tested against a negative set that is a quarter
other species.

## Reviewer comment → stage map

| comment | stage | status |
|---|---|---|
| R2 #1 — transient acquisition not demonstrated | 18, 19 | new analysis; claim supported, γ claim withdrawn |
| R2 #2 — 60 of 132 core genes, possible bias | 21 | acknowledged, supplementary table |
| R2 #3 — consensus is not the ancestral state | 17 | new analysis; α unchanged |
| R2 #4 — branch support on Figure 3 | 21 | caption/tree ambiguity resolved first |
| R1 — refocus on H3/H4, shorten | — | manuscript restructuring, no new analysis |
| editor — overstated conclusions, limitations | 19 | provenance audit |
