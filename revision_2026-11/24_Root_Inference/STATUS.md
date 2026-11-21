# Stage 24 — Outgroup-free inference of the capripoxvirus root

## Why

The McDonald–Kreitman analysis polarizes divergence by treating the MRCA of
the SPPV+GTPV clade as the ancestral state for LSDV. That is valid only if
LSDV is the earliest-diverging lineage of the genus. The manuscript supports
that rooting with a bootstrap value of 100, but a bootstrap measures
confidence in a split, not in the placement of the root, and an
outgroup-rooted phylogeny cannot test its own rooting assumption.

The point is live. The published ancient-genome analysis of sheeppox
([Mühlemann et al., 2026](https://www.biorxiv.org/content/10.64898/2026.02.23.707469v2.full))
reports that its core-genome tree only weakly supports LSDV as basal, prefers
SPPV as the earliest-diverging lineage, obtains materially different genus
ages under the two rootings (median 7,286 BP against about 5,000 BP), and
concludes that ancient LSDV or GTPV genomes are needed to settle it.

This stage tests the root from the data instead, under a non-reversible
substitution model, which carries directional information and therefore does
not need an outgroup.

## What was run

The 60 single-copy core orthogroups were concatenated into a nucleotide
supermatrix of 290 taxa and 48,645 sites. The same unrooted ML topology was
then rooted three ways — LSDV basal, SPPV basal, GTPV basal — and the three
rooted trees compared under `UNREST+G4`, the fully general 12-parameter
non-reversible model, by the approximately unbiased test with 10,000 RELL
replicates. The test was repeated on ten replicates that sample 16 genomes
from each species, because the full dataset is heavily unbalanced
(240 / 34 / 16) and non-reversible rooting can be pulled toward the densely
sampled group.

A separate unconstrained search over all 365 candidate root branches was also
run, for the limitation noted below.

## Results

Full dataset, 290 genomes:

| Rooting | logL | ΔlogL | RELL | ELW | AU *p* |
|---|---|---|---|---|---|
| **LSDV basal** | −92012.69 | 0 | **0.937** | **0.888** | **0.945** |
| SPPV basal | −92017.82 | 5.13 | 0.030 | 0.056 | 0.071 |
| GTPV basal | −92017.81 | 5.12 | 0.033 | 0.056 | 0.072 |

Balanced subsamples, 16 genomes per species, 10 replicates:

| Rooting | Times best | Mean ELW | Mean AU *p* |
|---|---|---|---|
| **LSDV basal** | **10 / 10** | **0.914** | **0.962** |
| SPPV basal | 0 / 10 | 0.043 | 0.050 |
| GTPV basal | 0 / 10 | 0.043 | 0.049 |

LSDV-basal rooting is preferred on the full dataset and in every balanced
replicate, carrying about 89–91% of the expected likelihood weight. The two
alternatives are excluded by the RELL bootstrap on the full data and sit at
or below the 5% AU boundary on balanced data. Balancing the sampling
strengthens rather than weakens the result, so the preference is not an
artefact of LSDV being over-represented.

## Limitation, which must be reported

The unconstrained search over all 365 candidate root branches finds its
optimum **inside** a clade rather than at any species split: the best branch
reaches logL −91992.24, about 20 log-likelihood units better than the best
species-level rooting at −92012.69. The non-reversible model therefore has
rooting signal that is not sharply localized at the deep divergences, and
part of what it is fitting is within-clade compositional structure. The
honest reading is that this analysis is **consistent with and supportive of**
LSDV-basal rooting among the three biologically meaningful alternatives, not
that it resolves the root of the genus. The ancient-genome paper's conclusion
that ancient LSDV or GTPV material is needed for that still stands.

## Consequence for the manuscript

The MK polarization rests on an assumption that now has independent support
from a method that does not share the assumption. The Discussion should state
the discrepancy with the ancient-genome analysis, present this test, and note
that under the disfavoured alternatives the branch to which adaptive
divergence is attributed would have to be reinterpreted, even though the test
itself — within-LSDV polymorphism against LSDV-to-outgroup divergence —
remains valid either way.

## Files

| Path | Content |
|---|---|
| `scripts/01_build_nt_supermatrix.py` | concatenates the 60 codon alignments |
| `scripts/02_test_species_rootings.py` | builds the three rootings, runs the AU test |
| `scripts/03_balanced_sensitivity.py` | the balanced replicate analysis |
| `01_inputs/{LSDV,SPPV,GTPV}_basal.nwk` | the three rooted topologies |
| `02_results/au_test_three_rootings.iqtree` | full IQ-TREE report |
| `02_results/all_branch_roottest.csv` | unconstrained test over all 365 branches |
| `02_results/balanced_roottest.csv` | per-replicate results |
| `03_summary/species_rootings.csv` | the main table |
| `03_summary/balanced_summary.csv` | the balanced summary |

## Reproducing

```bash
export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
cd "$PANG/revision_2026-11"
conda activate iqtree_env
python3 24_Root_Inference/scripts/01_build_nt_supermatrix.py
python3 24_Root_Inference/scripts/02_test_species_rootings.py
python3 24_Root_Inference/scripts/03_balanced_sensitivity.py
```

Note on models: `NONREV` is rejected as a model name by some IQ-TREE builds.
`UNREST` (equivalently `12.12`) is the fully general non-reversible model and
is used here.
