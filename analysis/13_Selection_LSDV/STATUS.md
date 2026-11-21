# 13 BUSTED on 240 LSDV genomes — NEGATIVE + DIAGNOSTIC

All 60 jobs completed; 4 genes had zero polymorphic nucleotide sites and were
excluded, leaving 56 testable. 11 genes reached q<0.05.

These 11 are NOT reportable as positive selection. Diagnostics:
  * prop_dS_near_zero is 0.88-0.96 for ALL 11 significant genes (dS_max 8-27),
    so the omega>1 rate class is populated by sites where dS -> 0.
  * median number of distinct protein sequences per gene is 6 out of 240;
    27 of the 60 genes have fewer than 5. Effective sample size is ~6, not 240.
  * LRT = 0 for 30 of 56 genes (HyPhy reports p = 0.5 in that case, which is
    "untestable", not "rejected").
  * The 4 genes with omega_MG94 > 1 are all non-significant, and the
    significant genes all have low omega — the inverse of a real
    diversifying-selection pattern.

VALUE OF THIS STAGE: it quantifies why codon models cannot be applied within
LSDV. Keep these outputs; they are the evidence for the limitations section.
The within-LSDV selection question is answered instead in stage 14.
