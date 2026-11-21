# Scripts

The repository carries two script sets, and the distinction matters.

## `bin/` — pipeline scripts

Portable command-line tools called by the Nextflow modules. Each takes explicit
arguments and writes to explicit paths, so a run is reproducible on any machine
without a fixed directory layout.

| Script | Stage | Purpose |
|---|---|---|
| `pangenome_qc.py` | 05 | contamination and completeness verdicts per genome |
| `partition_pangenome.py` | 06 | core / soft-core / shell / cloud partitioning, single-copy list |
| `heaps_law.py` | 06 | accumulation curve and Heaps' law fit |
| `concat_alignments.py` | 07 | supermatrix and NEXUS partition file |
| `rm_ratio.py` | 09 | recombination-to-mutation ratio from Gubbins output |
| `run_fastbaps.R` | 10 | Bayesian clustering with prior sensitivity |
| `balanced_subsample.py` | 11 | continentally balanced replicates for the root test |
| `root_to_tip.py` | 12 | root-to-tip regression, stratified by continent |
| `scoary_summarise.py` | 13 | per-trait pan-GWAS summary, including null results |
| `busted_summary.py` | 14 | BUSTED table plus dN/dS applicability diagnostics |

## `bin/legacy/` — the published analysis

The scripts that produced the results reported in the manuscript, preserved
verbatim. They expect the original project directory layout and are run from the
project root rather than through Nextflow. They are kept because they are the
provenance record: if a number in the paper is questioned, this is the code that
generated it.

Several `bin/` scripts are portable reimplementations of a `bin/legacy/`
counterpart — `scoary_summarise.py` of `build_S6_pangwas_table.py`, and
`busted_summary.py` of `parse_busted.py`. They implement the same criteria, but
they are not byte-identical, so before relying on a pipeline run for a published
figure, check the values against the archived results.

Several stages of the original analysis were run interactively rather than from
a saved script; the `bin/` implementations of those steps
(`pangenome_qc.py`, `partition_pangenome.py`, `heaps_law.py`,
`concat_alignments.py`, `rm_ratio.py`, `balanced_subsample.py`,
`root_to_tip.py`) are new code written for this pipeline. Validate them against
the published values before treating their output as a reproduction.

## Validation targets

A full run on the archived samplesheet should return:

| Quantity | Expected |
|---|---|
| Genomes entering the pangenome | 290 |
| Orthogroups | 206 |
| Strict core | 132 |
| Single-copy core | 60 |
| Heaps' gamma | 0.0679 |
| Recombination blocks | 90 |
| r/m | 0.13 |
| fastBAPS clusters | 5 |
| Significant gene-trait pairs | 92 across 16 traits, 18 orthogroups |
| Polarised alpha | 0.604 (95% CI 0.33-0.75) |
