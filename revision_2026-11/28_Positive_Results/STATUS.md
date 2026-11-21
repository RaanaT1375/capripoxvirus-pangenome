# Stage 28 — five positive results from data already on disk

None of these needs new sequencing, new jobs of any size, or anything the
cluster does not already hold. Four run in seconds; one needs stage 17 to
have been run first.

## Run it

```bash
export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
cd "$PANG/revision_2026-11/28_Positive_Results"
conda activate iqtree_env

python3 scripts/01_constraint_structure.py     # constraint by functional class
python3 scripts/03_dna_pangenome.py            # the sequence-based pan-genome
python3 scripts/04_diagnostic_markers.py       # diagnostic SNP panel
python3 scripts/05_breakpoint_map.py           # breakpoint map and DIVA targets

conda activate r_figs
Rscript R/constraint.R                         # figure for step 1
```

Step 2 needs the ancestral states from stage 17, which is not currently on
the cluster. Unpack `revision_2026-11.zip` (sent alongside this), then:

```bash
conda activate iqtree_env
cd "$PANG/revision_2026-11/17_ASR_MK"
python3 scripts/01_run_asr.py                  # IQ-TREE ASR, 60 genes
cd "$PANG/revision_2026-11/28_Positive_Results"
python3 scripts/02_ancestral_substitutions.py  # name the fixed substitutions
```

## What each step gives you

**Step 1 — constraint is structured by function.** Your inter-species omega
estimates are valid; the codon models failed only *within* LSDV. At the
moment the median of 0.169 appears only as a preamble to a negative result.
This classifies the 60 core genes by functional class and tests whether
omega differs across them (Kruskal-Wallis, then pairwise Mann-Whitney with
Holm correction). If the expected gradient holds — replication and
transcription most constrained, host-interaction least — that is a sentence
in the abstract and a figure.

The classification comes from product annotation by keyword. Check it: the
script prints every gene it could not classify, and you can override the
whole thing with `01_inputs/functional_categories.csv` (columns
`orthogroup,category`). A reviewer will scrutinise the class assignment more
than the test, so it is worth ten minutes of your own judgement.

**Step 2 — name the substitutions.** This is the one I would prioritise.
Your positive result is currently a single number; stage 17 already
reconstructed an ancestral codon at 16,018 of 16,215 sites, so the specific
changes can be listed: gene, codon, ancestral residue, derived residue. With
alpha near 0.6 roughly three in five of the non-synonymous entries are
expected to be adaptive, though the test cannot say which. It also
classifies each change as radical or conservative by physicochemical class,
because that is the first thing a reviewer will ask.

The integer counts here should come close to the MK test's Dn = 74.8 and
Ds = 83.2 without matching exactly, since Nei-Gojobori splits multi-step
codons fractionally and this counts events whole. If they differ by more
than a few, tell me.

**Step 3 — the sequence-based pan-genome.** Stage 27's matrix has not yet
been used to describe the pan-genome itself. It is more closed than the
gene-call version (strict core 139 against 132, cloud 35 against 50), and
as far as I can tell it is the first sequence-based pan-genome for this
genus. This is the positive form of the pan-GWAS negative: the same
correction that destroys the associations improves the pan-genome
description.

**Step 4 — diagnostic SNP panel.** Finds alignment columns at which one
group is fixed for one base and everything else for another: species,
vaccine versus field, recombinant versus clonal, and each fastBAPS cluster.
Translational, entirely positive, and the data are already there. Sites
tolerating one or two exceptions are reported separately from perfect ones,
since a strict rule can return nothing. Cross-check the chosen sites against
step 5's intervals before reporting them — a marker inside a tract that
moves between lineages is not a stable assay target.

**Step 5 — breakpoint map.** Your 90 blocks currently support one sentence.
This tests whether breakpoints are uniformly distributed (chi-square and
Kolmogorov-Smirnov, plus a binomial test of outer thirds against the central
third), finds intervals hit by two or more independent blocks, and extracts
the Node_147 junctions. Those junctions are the DIVA targets worth naming: a
stable ancestral junction is a better assay target than a transient tract.

## Honest expectations

I cannot promise any of these returns a positive result — the pan-GWAS went
the other way tonight. Steps 2, 3 and 5 are near-certain to produce
reportable output because they are descriptive. Steps 1 and 4 are genuine
tests and could come back null: the classes might not differ, and there
might be no perfectly diagnostic site outside recombinant regions. If that
happens we report it. Each script prints what it found rather than only what
it hoped to find.

## Files

| Path | Content |
|---|---|
| `scripts/_s28.py` | shared helpers, the functional-class keyword rules |
| `scripts/01_constraint_structure.py` | omega by functional class, with tests |
| `scripts/02_ancestral_substitutions.py` | the fixed substitutions, named |
| `scripts/03_dna_pangenome.py` | partition and Heaps on both matrices |
| `scripts/04_diagnostic_markers.py` | diagnostic sites and a minimal panel |
| `scripts/05_breakpoint_map.py` | breakpoint density, hotspots, DIVA targets |
| `R/constraint.R` | the figure for step 1, in the stage 26 theme |
| `01_inputs/functional_categories.csv` | optional override for step 1 |

## Tested

All five scripts and the R figure were run end to end against synthetic
inputs before being sent: a selection table with a planted constraint
gradient, a core alignment with planted species-diagnostic columns, a
Gubbins GFF with 90 blocks of which nine are attributed to Node_147, and
stage 27's matrix. Step 1 recovered the planted gradient, step 4 recovered
the planted diagnostic columns including one with an exception, step 5
recovered the nine ancestral blocks, and step 2 exits with instructions when
stage 17 output is absent. The numbers from those runs are meaningless; only
the mechanics were tested.

## step 7 — substitution clustering screen (added 2026-11-20)

Reason: step 2's populated `substitutions` column showed that in some genes
every fixed difference falls inside a few codons. OG0000120's 16 changes all
sit in codons 177-205, which is the 29-codon unalignable C-terminal block the
manuscript already attributes its aBSREL signal to -- recovered independently.
OG0000067 shows the same signature (4 changes in codons 60-65, Ds = 0).

`07_substitution_clustering.py` applies, to every gene with >= 3
non-synonymous fixed differences, a permutation test on the largest number of
substitutions in any 30-codon window, a KS uniformity test, and a local gap
frequency check against the gene's own mean.

Scope: the screen governs which loci may be NAMED as candidates. It does not
touch the MK test -- removing clustered sites post hoc would bias alpha.
Report both totals and the screen.

### corrections after the first run (2026-11-20)

The first version had three defects, all found by reading its own output:

1. The gap rule had no absolute floor, so a ratio on a near-zero baseline
   flagged noise. OG0000122 was flagged at "4.4x gappier" = 0.23% against
   0.05%, under one gapped cell in three hundred; OG0000057 at "34x" =
   2.41% against 0.07%, about five cells in two hundred. Both were false
   positives. Fixed with GAP_FLOOR (default 5% local gap frequency).
2. The gap rule ran on genes below MIN_N. OG0000057 has n = 2, was never
   position-tested, and was flagged anyway. Now restricted to tested genes.
3. A fixed 30-codon window is wrong for short genes and missed a real
   cluster. OG0000067 is 89 codons, so a 30-codon window is a third of the
   sequence: p_window = 0.102 although its four substitutions span six
   codons. Replaced by an exact, scale-free span probability
   (p = 3.5e-04, Holm 2.4e-03), Holm-corrected over the tested genes. The
   window statistic is retained in the table for comparison only.

Flagged after the correction: OG0000120 and OG0000067, not OG0000122 or
OG0000057. Removal is therefore 16 of 60 non-synonymous fixed differences
(27%), not 21 (35%), leaving 44 across 18 genes -- confirm against the
re-run.

### the gap column is circular -- retracted

An earlier note here claimed the gap evidence contradicted the manuscript's
"29-codon unalignable block". That was wrong and is withdrawn.

Gap frequency at OG0000120's substitution columns is 0.0% against 14.5% for
the gene, but step 2 admits a site only when every codon position
reconstructs at PP >= 0.90 and the ingroup is strictly monomorphic. A column
gapped in a third of the taxa cannot satisfy that, so the filter SELECTS
gap-free columns and a low gap figure at them is guaranteed. It is not
evidence about the surrounding block.

The manuscript measures gaps the right way: a sliding window over the
29-codon block records 36.1% gaps, attributed to a species-specific
7-residue indel and a divergent hydrophobic tail. That stands. A misaligned
divergent tail produces both indels and apparent substitutions, and the
substitutions land on whichever columns align without gaps. The manuscript's
wording needs no change.

### what step 7 does contribute

1. It reproduces the manuscript's own KS test on OG0000120 independently
   (D = 0.804 here against D = 0.802 published) without being told the
   answer, which is a check on that published result.
2. It finds the same signature at OG0000067, which the manuscript attributes
   to synonymous-rate collapse but does not report as positionally
   clustered: 4 non-synonymous fixed differences inside 6 codons of 89
   (p = 3.5e-04, Holm 2.4e-03). A second, independent line of evidence for
   a claim the manuscript already makes.
3. It clears the other six tested genes: their substitutions are
   positionally unremarkable, so the 40 differences across 18 genes may be
   named.
4. OG0000111 is borderline -- raw p = 0.035, Holm 0.209, so it does not
   survive correction. Report it, because a reader with the CSV will ask.
