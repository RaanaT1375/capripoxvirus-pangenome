# Caveats and interpretive limits

Every constraint on how the results in this repository may be read, collected
in one place. Each entry names the stage it applies to.

## Dataset-level

**The dataset is genus-level, not species-level.** 290 genomes = 240 LSDV +
34 SPPV + 16 GTPV. This was discovered mid-project. Any analysis run on all
290 genomes measures between-species divergence unless explicitly restricted.
Affects stages 06, 07, 08, 10, 11.

**Two heterogeneous halves.** 237 GenBank single-contig assemblies (median
150,497 bp) and 53 SRA multi-contig assemblies (up to 9.7 Mb, including
non-viral sequence). Removing the SRA half changes α by 0.002, so it is not
a source of bias for stage 14, but no population-level STR catalogue was
attempted for this reason.

**Vaccine strains.** 20 culture-passaged vaccine strains carry private
non-synonymous variants. They are excluded from the primary MK analysis.
Including them lowers α and introduces a spurious dependence on the
frequency cutoff.

## Stage-specific

**06 Recombination.** Gubbins GFF coordinates are Parsnp-alignment
coordinates, not reference-genome coordinates. Mapping recombinant regions to
genes requires converting via `parsnp.maf`/`xmfa`. r/m = 0.13 was computed at
genus level.

**07 Phylogeography.** Run on all three species. The root state is
unresolved regardless of sampling strategy. Continental transitions should
not be read as LSDV epidemiology.

**08 Population structure.** fastBAPS clusters 2 and 3 correspond to SPPV and
GTPV — they are species, not LSDV lineages.

**09 Temporal signal.** Use `*_LSDV_only*` files. The apparent negative slope
comes from sampling confounding: all 20 pre-2000 genomes are African, while
Asia (120) and Europe (44) are entirely post-2000. Restricted to 2000+, the
slope is zero (p = 0.887). No molecular clock; no BEAST/TMRCA was attempted.

**10 BUSTED (290 genomes).** Valid, but measures between-species divergence.
Must not be described as within-LSDV selection.

**11 Scoary.** All BAPS and host traits have `n_robust = 0` and
`max_pairs = 1` — zero independent evolutionary events. These associations
are phylogenetic confounding. `Host_Goat` is GTPV and `Host_Sheep` is SPPV,
so the "host-associated genes" are species markers. Only 11 of 92 pairs have
at least three supporting pairs; four have more opposing than supporting.

**12 MEME.** The only two sites reaching q<0.05 have synonymous rate near
zero and 1–2 supporting branches — the same dS-underflow pathology as
stage 13. Not independent confirmation.

**13 BUSTED (240 LSDV).** Negative and diagnostic. Median 6 distinct protein
sequences per gene across 240 genomes; effective sample size is ~6, not 240.
LRT = 0 for 30 of 56 testable genes, which HyPhy reports as p = 0.5 —
"untestable", not "rejected". The 4 genes with ω > 1 are all
non-significant and the significant genes all have low ω: the inverse of a
real diversifying-selection pattern.

**14 McDonald-Kreitman.** Genome-wide result; no single gene passes FDR and
no gene-specific claim may be made from it. Only 9 of 60 genes have ≥5 fixed
differences on the LSDV branch; 17 have none. α is a floor: slightly
deleterious variants and the conservative treatment of 239 sites that are
both polymorphic and divergent both depress it. Partly reuses the divergence
counts of stage 10, so the two are not fully independent.

**15 RELAX.** Inconclusive, and reported as a limitation rather than a
result. 24 of 59 genes (41%) fail with a rate-distribution collapse, so the
36 survivors are not a random subset. Of those, 16 reach p<0.05 with no
coherent direction (11 intensification, 5 relaxation). The test clade has
almost no divergence, which pushes the ω distribution to its bounds; K > 1
here is not biological.

**16 aBSREL.** Negative. Two genes reach significance and both are rejected:
`OG0000067` has dS = 0 across the entire gene, and `OG0000120` (OX-2) is an
alignment artefact — see below. 11 genes could not be tested because
Dn = Ds = 0 on the labelled branch.

## The OG0000120 (OX-2) case

This gene must **not** be cited as evidence of selection anywhere. It looked
compelling: longest stem branch of 300 in the GTR tree, ω = 697.6,
LRT = 103.5, corrected p = 0, with SPPV and GTPV stems giving LRT = 0.

It is an artefact. All 14 fixed differences fall in codons 174–202 of 217
(KS D = 0.802, p < 1e-9). The weight of the "positively selected" rate class
(0.141) matches that block's share of the gene (29/217 = 0.134): the rate
class is the block. The block sits on the gappiest region of the alignment
(0% gaps in codons 20–160, 29% in 180–200, 36% in 200–217), where LSDV
carries a 7-residue indel and a divergent hydrophobic tail. dN/dS is
undefined in an indel-containing region. Codons 0–20 are 93.5% gapped,
indicating inconsistent ORF start annotation.

Clustering is not a general problem: of the 9 genes with ≥5 fixed
differences, only OG0000120 is significantly clustered after FDR. The MK
result is therefore unaffected.

## Gene-level candidate — discussion only

`OG0000118` (LSDV135, putative IFN-α/β binding protein) is elevated in MK
against GTPV (p = 0.0085) and in RELAX (K = 4.67, p = 5e-06) with a clean
fit, and its divergence is dispersed rather than clustered (span = 0.795,
KS q = 0.135). It does not pass FDR on its own and is offered as a
hypothesis, not a finding.
