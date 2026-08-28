# 10 Selection Pressure (BUSTED, 290 genomes) — VALID, GENUS-LEVEL

The computation is sound; the original interpretation was not. This alignment
set contains all three species, so significance here reflects divergence
BETWEEN LSDV, SPPV and GTPV, not adaptation within LSDV.

Do NOT write "positive selection in LSDV" from these results. For the
within-LSDV question see stage 13 (negative) and stage 14 (MK).

Results: 60 genes; 10 with q<0.05; 48 purifying; median omega_MG94 = 0.169;
no gene with omega > 1. Run with SRV (synonymous rate variation).
Note: 4 genes initially failed on numerical instability in the GTR fit and
succeeded with ENV=TOLERATE_NUMERICAL_ERRORS=1 (job 647755) — a convergence
issue, not a data problem.

This directory also holds the inputs reused by stages 13 and 14:
`01_codon_alignments/` (60 pal2nal codon alignments, 290 taxa) and
`02_pruned_trees/`.
