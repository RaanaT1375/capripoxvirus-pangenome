# Analysis Status — Capripoxvirus Pan-genome Project

Dataset: 290 genomes — **240 LSDV, 34 SPPV, 16 GTPV** (genus level, not
species level). This distinction was discovered mid-project and invalidates
the original interpretation of several stages. Read this file before citing
any number.

Topology (`core_genome_ML_tree_v3.contree`): SPPV+GTPV form a clade,
bootstrap 100. LSDV is the outgroup to that clade.

| Stage | Status | Notes |
|---|---|---|
| 01 QC & filtering | FINAL | 323 -> 290 genomes |
| 02-04 Pan-genome (OrthoFinder) | FINAL, caveat | 206 orthogroups; some families split into allelic orthogroups |
| 05 Phylogeny (IQ-TREE) | FINAL | 60 single-copy genes, 290 x 16,215 |
| 06 Recombination (Gubbins) | FINAL, caveat | r/m=0.13; GFF coords are Parsnp-alignment coords |
| 07 Phylogeography (PastML) | FINAL, caveat | run on all 3 species; root state unresolved |
| 08 Population structure (fastBAPS) | FINAL, caveat | clusters 2 and 3 are SPPV and GTPV |
| 09 Temporal signal | FINAL (revised) | use `*_LSDV_only*` files only |
| 10 BUSTED, 290 genomes | VALID but genus-level | not within-LSDV selection |
| 11 Scoary | FINAL (tree version) | use `scoary_results_tree/` |
| 12 MEME | NEGATIVE | signal is a dS~0 artefact |
| 13 BUSTED, 240 LSDV | NEGATIVE + DIAGNOSTIC | basis of the codon-model limitation section |
| 14 McDonald-Kreitman | **FINAL — primary result** | alpha = 0.604, 95% CI 0.33-0.75 |
| 15 RELAX | INCONCLUSIVE | 24/59 genes failed; no coherent direction |
| 16 aBSREL, LSDV stem | NEGATIVE | 2 significant, both rejected as artefacts |

## Headline results
1. Within present-day LSDV there is almost no diversity (median 6 distinct
   proteins per gene across 240 genomes) and strong purifying selection
   (median omega = 0.105). No temporal signal.
2. Codon models (BUSTED, MEME) are not applicable within LSDV: ~90% of sites
   have dS ~ 0, and LRT = 0 for 30 of 56 testable genes.
3. McDonald-Kreitman on the LSDV branch gives alpha = 0.604
   (bootstrap 95% CI 0.33-0.75; CMH p = 0.009), i.e. adaptive amino-acid
   divergence happened during LSDV's divergence from the SPPV/GTPV ancestor,
   not during its recent clonal expansion.
4. The MK signal is genome-wide, not gene-specific. No single gene passes
   FDR. Seven independent robustness tests pass, including an APOBEC3
   context check that found no editing signature.
5. OG0000120 (OX-2) is an alignment artefact, not an adaptive signal: all 14
   fixed differences fall in codons 174-202, which is also the most gap-rich
   region of the gene (29-36% gaps). It must not be cited as evidence of
   selection anywhere in the manuscript.

## Central narrative
Adaptation in the deep past (alpha ~ 0.6 on the LSDV stem) versus stasis in
the present (6 haplotypes, omega = 0.105, no molecular clock).

## Status
All analyses complete (stages 01-16). Per-stage STATUS.md files are
authoritative for numbers. Remaining: final report, then public repository.
