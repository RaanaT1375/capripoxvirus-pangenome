# Stage 25 — McDonald–Kreitman with recombinant regions masked

## The problem this stage exists to settle

The distance distribution of LSDV core haplotypes is bimodal, and the two
modes correspond exactly to recombination ancestry: the 139 genomes with no
inferred recombination (BAPS cluster 5) sit within 50 substitutions of the
consensus, while the 101 genomes descending from the Node_147 event (BAPS
clusters 1 and 4) sit 100–250 substitutions away.

Splitting the MK ingroup along that line removes the adaptive signal:

| Ingroup | n | Pn | Ps | Dn | Ds | α | CMH *p* |
|---|---|---|---|---|---|---|---|
| all non-vaccine LSDV | 227 | 146.0 | 379.0 | 69.8 | 67.2 | **0.629** | 9.0 × 10⁻³ |
| non-recombinant only | 133 | 29.7 | 39.3 | 108.8 | 159.2 | **−0.104** | 0.914 |
| recombinant-descended only | 94 | 91.0 | 206.0 | 78.8 | 85.2 | **0.522** | 2.8 × 10⁻² |

The mechanism is visible in the counts. Dropping the recombinant-descended
genomes cuts Ps by a factor of ten (379 → 39) but Pn only by a factor of five
(146 → 30): the polymorphism contributed by those lineages is
disproportionately synonymous, which is exactly what makes α positive. The
McDonald–Kreitman test assumes polymorphism and divergence arise from the
same mutational process. Recombination imports old variants from a divergent
donor lineage instead, and the manuscript itself reports that this affected
29.2% of the genome with 101 genomes descending from a single event.

Dropping genomes is nevertheless the wrong test. It discards three quarters
of the standing variation, and it reclassifies sites from polymorphism to
divergence — which is why Dn and Ds *rise* in the non-recombinant row. The
correct test masks the recombinant **regions** and keeps every genome.

## Why Gubbins is re-run here

The manuscript's Gubbins run was performed on the Parsnp whole-genome
alignment. Its GFF coordinates are Parsnp alignment coordinates, not genome
coordinates, so mapping them onto the per-gene codon alignments would require
going back through the Parsnp LCB blocks. Running Gubbins directly on the
60-gene LSDV supermatrix removes that step entirely: the predicted intervals
come back in supermatrix coordinates and the partition table maps them onto
genes and codons exactly.

A codon is masked when any of its three positions falls inside a predicted
recombinant interval **on any branch**. That is deliberately conservative: it
removes every site recombination is inferred to have touched anywhere in the
tree, not only in the genomes carrying it.

## How to run it

```bash
export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
cd "$PANG/revision_2026-11/25_Recombination_Masked_MK"

conda activate iqtree_env          # or any env with Biopython + pandas
python3 scripts/01_build_lsdv_supermatrix.py

sbatch scripts/02_run_gubbins.sbatch     # submit from this directory

# once the job finishes
python3 scripts/03_masked_mk.py
```

Step 3 reads the ancestral states produced by stage 17, so stage 17 must have
been run first. The three unmasked rows in the table above are recomputed by
step 3 as controls, so they should reproduce exactly; if they do not, the
inputs have changed.

## Step 3 is retracted — do not report its masked row

Step 3 ran Gubbins on the 60-gene concatenation, and that design is wrong.
With no physical adjacency between neighbouring genes, rate heterogeneity
*between* genes creates blocks of high and low SNP density that Gubbins reads
as recombination. The result was a mask that operates at whole-gene
granularity: 24,118 of 48,645 sites (49.6%), 8,054 of 16,215 codons (49.7%),
with many orthogroups masked at exactly 100% and 17 untouched entirely. It
removed 80% of the non-synonymous fixed differences (Dn 69.8 → 14.0) and
dropped Dn/Ds from 1.04 to 0.44.

That is circular: the mask deletes the fast-evolving genes wholesale, which
is precisely where adaptive divergence would sit. The masked value it
produced (α = 0.304, CMH *p* = 0.080) is therefore not a
recombination-corrected estimate and must not appear in the manuscript.

The three unmasked rows stand, because they partition genomes rather than
sites.

Step 4 replaces it, using the stage-06 Gubbins run on the Parsnp whole-genome
core alignment, where adjacency is real. Its coordinates are columns of
`parsnp_core_alignment.fasta` and are translated in three hops — core-alignment
column → `NC_003027.1` genome position via the reference row of each Parsnp
LCB → codon, by locating each reference CDS in the genome. The script
validates the translation before using it: the reference row rebuilt from the
XMFA must match the core alignment's reference row exactly, and it computes
nothing if they differ. It also reports a per-genome mask, using each
interval's `taxa` field, which is what Gubbins actually infers and which
discards far less data than a union mask.

```bash
python3 scripts/04_parsnp_masked_mk.py
```

The control row in that output must reproduce α = 0.629.

## How to read the result

If α stays positive and significant with recombinant codons masked while all
227 genomes are retained, the adaptive signal is not an artefact of
introgressed polymorphism, and the central claim of the manuscript stands —
strengthened, because the most obvious alternative explanation has been
tested and rejected.

If α collapses towards zero as it does when the recombinant genomes are
dropped, then what the original analysis measured was the synonymous-rich
polymorphism imported by recombination, not adaptive fixation along the LSDV
branch. In that case the manuscript's headline result has to be rewritten,
and the honest version is still publishable: it becomes a demonstration that
a polymorphism-to-divergence test applied to a recombining clonal virus
recovers the recombination history rather than adaptation, which is a
methodological result with reach beyond this genus.

Either way the number must be obtained before the Discussion is written.

## Files

| Path | Content |
|---|---|
| `scripts/01_build_lsdv_supermatrix.py` | 240 LSDV genomes × 48,645 nt, plus the partition table |
| `scripts/02_run_gubbins.sbatch` | Gubbins on that supermatrix |
| `scripts/03_masked_mk.py` | GFF → codon mask, then the four MK analyses |
| `01_inputs/lsdv_super_nt.fasta` | the alignment (regenerate; not committed) |
| `01_inputs/partitions_nt.txt` | orthogroup → supermatrix interval |
| `02_results/masked_codons.csv` | masked codons per orthogroup |
| `03_summary/mk_recombination_masked.csv` | the four-row result table |

## Note

The code path was validated end to end against a synthetic GFF before the
real Gubbins run existed; the three unmasked control rows reproduced the
reference values exactly. No synthetic file is retained in this stage — if a
`recombination_predictions.gff` is present, it came from step 2.
