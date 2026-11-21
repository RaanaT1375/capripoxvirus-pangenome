# Stage 27 — the pan-GWAS redone on sequence instead of gene calls

Stage 20 tested 34 orthogroups at the DNA level and showed the
annotation-based presence/absence matrix does not describe gene content.
The three strongest associations turned out to be present in 269–289 of 290
genomes once the sequence itself was searched:

| Orthogroup | Trait | Annotated in | Present at DNA level |
|---|---|---|---|
| OG0000157 | Asia | 26 | **289** |
| OG0000154 | Europe | 114 | **289** |
| OG0000152 | Vaccine | 91 | **269** |

This stage rebuilds the whole matrix that way and re-runs the association
analysis on it. Either real associations survive, or none does and the
negative result becomes clean.

## Run it

```bash
export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
cd "$PANG/revision_2026-11/27_DNA_Level_GWAS"

conda activate prokka_env           # BLAST+
python3 scripts/01_build_queries.py

sbatch scripts/02_run_tblastn.sbatch        # submit from this directory

# once it finishes
conda activate iqtree_env           # pandas, Biopython
python3 scripts/03_build_dna_matrix.py

sbatch scripts/04_run_scoary.sbatch         # submit from this directory

# once that finishes
conda activate iqtree_env
python3 scripts/05_compare_results.py
```

Step 1 reuses the BLAST database stage 20 already built, so it takes
seconds. Step 2 is six times the queries of stage 20, which ran in under a
minute, so expect a few minutes. Step 4 is the long one: Scoary with 1,000
permutations over 19 traits.

## What each step prints that you should check

| Step | Check |
|---|---|
| 1 | 206 queries written; if any orthogroup has no sequence file it says so |
| 3 | how many orthogroups gained or lost genomes, and the compartment shift |
| 3 | genomes in common between matrix, tree and traits — should be 290 |
| 5 | the convergent count, against 7 from the annotation-based matrix |
| 5 | the verdict on each of the old seven |

## The decision this stage settles

Step 5 prints one line per association from the old analysis, with one of
three verdicts.

**still convergent** — supported by sequence rather than by gene calls.
Report it as a biological result; it is stronger than before, because it has
survived a test the original did not face.

**significant, not convergent** — the presence/absence pattern is real but
confined to too few independent branches to separate convergence from shared
ancestry.

**no longer significant** — an artefact of annotation.

If nothing survives, that is not a failure. It makes the negative result
clean: a pan-GWAS on viral genomes annotated by heterogeneous pipelines
recovers annotation provenance rather than biology, demonstrated directly by
rebuilding the matrix from sequence. That is a methodological finding with
reach beyond this genus, and it belongs in the same part of the paper as the
codon-model result.

## Two things that could go wrong

**`--start_col`.** Step 3 writes the same 14 metadata columns as the
original Scoary input, so the genome columns begin at column 15 and the
original `--start_col 15` applies unchanged. If Scoary complains about trait
or gene names, check the header of
`03_summary/gene_presence_absence_dna.csv` and set `START_COL` accordingly:

```bash
START_COL=15 sbatch scripts/04_run_scoary.sbatch
```

**Taxon mismatch.** Scoary refuses to run when the tree tips, the gene
matrix and the trait table do not name the same genomes. Step 3 prunes all
three to their intersection and writes matched copies, which is why step 4
points at `03_summary/` rather than at the original inputs. The count it
prints should be 290; if it is much lower, the names have drifted and the
intersection is wrong.

## Calling rule

Identical to stage 20, so the two are comparable: an orthogroup is present
in a genome when its best tblastn hit covers at least 80% of the query
protein at 80% identity or better, summing non-overlapping HSPs within a
contig so a gene split across a contig break is not scored absent. Override
with `MIN_COV` and `MIN_PID` if you want to see how sensitive the result is
to the threshold — worth doing once the main run is in, since a reviewer
will ask.

## Files

| Path | Content |
|---|---|
| `scripts/01_build_queries.py` | one representative protein per orthogroup |
| `scripts/02_run_tblastn.sbatch` | tblastn, all orthogroups vs all assemblies |
| `scripts/03_build_dna_matrix.py` | DNA-level matrix, comparison, Scoary inputs |
| `scripts/04_run_scoary.sbatch` | Scoary on the DNA-level matrix |
| `scripts/05_compare_results.py` | new result, and the fate of the old seven |
| `02_results/presence_dna.csv` | orthogroup x genome, 0/1 |
| `02_results/annotation_vs_dna.csv` | per-orthogroup gained/lost |
| `03_summary/scoary_dna_convergent.csv` | the new convergent set |
| `03_summary/old_vs_new.csv` | verdict on each of the old seven |

## Tested

Steps 3 and 5 were run end to end against synthetic tblastn and Scoary
output before being sent, so the parsing, the matrix construction, the tree
pruning and the comparison logic are known to execute. Steps 1, 2 and 4 are
thin wrappers around BLAST+ and Scoary and were not run here, since they
need the real database and the pinned Python 3.6 environment.
