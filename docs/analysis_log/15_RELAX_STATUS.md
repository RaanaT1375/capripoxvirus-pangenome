# 15 RELAX — INCONCLUSIVE (report as a limitation, not a result)

Design: test set = LSDV clade (472 labelled branches), reference set =
SPPV+GTPV clade (99 branches), --models Minimal. Job 648048, 1h27m.

OUTCOME: 35 of 59 genes converged, 24 failed (41%). The failures are a
recurring HyPhy error inside relax._renormalize_with_weights
("Max(USie_UMg.mean,0.001)" on an undefined value) raised while fitting the
K != 1 model. This is the rate distribution collapsing, not a random crash,
so the 36 parsed genes are a biased subset and NOTHING here generalises to
the genome.

Of the 36 parsed: 16 have p<0.05, split 11 intensification / 5 relaxation —
no consistent direction. Median K = 1.242. omega_max > 100 in 11 genes.
Degeneracy does not explain significance (Fisher p=1.00), but it does
inflate K (median 2.08 vs 1.12, Mann-Whitney p=0.043).

Applying a strict filter (0.01 < K < 10 and omega_max < 10, i.e. excluding
boundary estimates of K such as 0.0 and 5.5e-08) leaves four genes, all
intensification: OG0000118 (K=4.67, p=5e-06), OG0000060 (K=2.57),
OG0000088 (K=1.78), OG0000117 (K=1.50).

CAVEAT THAT MATTERS MOST: the test branch set has almost no divergence
(median 6 distinct proteins per gene across 240 genomes). When a branch set
carries no substitutions, its omega distribution is pushed to the
boundaries, which RELAX reads as intensification. K > 1 here should not be
interpreted as biology.

CONCLUSION: no evidence for a consistent change in selection intensity on
the LSDV lineage. Report the 41% non-convergence rate alongside any number
taken from this stage.
