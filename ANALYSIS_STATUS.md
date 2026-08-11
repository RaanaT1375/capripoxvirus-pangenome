# Analysis Status — Capripoxvirus Pan-genome Project

Dataset: 290 genomes — **240 LSDV, 34 SPPV, 16 GTPV** (genus level, not species level).
This distinction was discovered mid-project and invalidates the original
interpretation of several stages. Read this file before citing any number.

Topology (core_genome_ML_tree_v3.contree): SPPV+GTPV form a clade, bootstrap 100.
LSDV is the outgroup to that clade.

| Stage | Status | Notes |
|---|---|---|
| 01 QC & filtering | FINAL | 323 -> 290 genomes |
| 02-04 Pan-genome (OrthoFinder) | FINAL, caveat | 206 orthogroups; some families split into allelic orthogroups (see 03 STATUS) |
| 05 Phylogeny (IQ-TREE) | FINAL | 60 single-copy genes, 290 x 16,215 |
| 06 Recombination (Gubbins) | FINAL | r/m=0.13; GFF coords are Parsnp-alignment coords, not genome coords |
| 07 Phylogeography (PastML) | FINAL, caveat | run on all 3 species; root state unresolved |
| 08 Population structure (fastBAPS) | FINAL, caveat | clusters 2 and 3 are SPPV and GTPV, i.e. species, not LSDV lineages |
| 09 Temporal signal | FINAL (revised) | use `*_LSDV_only*` files only |
| 10 BUSTED, 290 genomes | VALID but genus-level | do NOT interpret as within-LSDV selection |
| 11 Scoary | FINAL (tree version) | use `scoary_results_tree/` |
| 12 MEME | NEGATIVE | signal is a dS~0 artifact |
| 13 BUSTED, 240 LSDV | NEGATIVE + DIAGNOSTIC | basis of the codon-model limitation section |
| 14 McDonald-Kreitman | FINAL — primary result | alpha = 0.604 on the LSDV branch, CMH p = 0.009; five robustness checks passed |
| 15 RELAX | INCONCLUSIVE | 41% non-convergence; no consistent direction |
| 16 aBSREL on LSDV stem | NEGATIVE | 2/49 significant, both rejected as artefacts |

## Headline results
1. Within present-day LSDV there is almost no diversity (median 6 distinct
   proteins per gene across 240 genomes) and strong purifying selection
   (median omega = 0.105). No temporal signal.
2. Codon models (BUSTED, MEME) are not applicable within LSDV: ~90% of sites
   have dS ~ 0, and LRT = 0 for 30 of 56 testable genes.
3. McDonald-Kreitman on the LSDV branch gives alpha = 0.604 (CMH p = 0.009),
   i.e. adaptive amino-acid divergence happened during LSDV's divergence from
   the SPPV/GTPV ancestor, not during its recent clonal expansion. Robust to
   two independent outgroups, removal of vaccine strains, removal of
   SRA-derived genomes, polarization, and masking of gap-rich columns
   (alpha 0.570-0.604 across all of them).
4. The MK signal is genome-wide, not gene-specific. No single gene passes FDR.
5. Branch-level likelihood tests add nothing. RELAX fails to converge on 41%
   of genes and shows no consistent direction. aBSREL finds two significant
   genes on the LSDV stem and both are artefacts: one has dS = 0 across the
   whole gene, and OG0000120 (OX-2) has all 14 of its fixed differences
   packed into a 29-codon indel-rich block at the C-terminus, which is
   exactly the rate class the model called "positively selected".
6. Best gene-level candidate for discussion only: OG0000118 (LSDV135,
   putative IFN-alpha/beta binding protein) - elevated in MK vs GTPV
   (p = 0.0085) and RELAX (K = 4.67, p = 5e-06) with a clean fit and
   dispersed rather than clustered divergence. Does not pass FDR alone.
