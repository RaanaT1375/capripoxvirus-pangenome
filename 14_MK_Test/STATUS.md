# 14 McDonald-Kreitman — FINAL, PRIMARY RESULT OF THE PROJECT

Why MK: within-LSDV dS is ~0 (see 13/STATUS.md), so any dN/dS is undefined.
MK contrasts LSDV polymorphism with divergence to an outgroup, where
divergence is ~3% and the denominator is well estimated.

Design (scripts/mk_test.py, mk_test_v2.py, mk_polarized.py):
  * two independent outgroups run separately (SPPV n=34, GTPV n=16)
  * ancestral state = outgroup majority codon, >=90% consensus required
  * sites polymorphic in LSDV counted as polymorphism only (conservative)
  * Nei-Gojobori pathway averaging for multi-nucleotide codon differences
  * variants below 5% frequency excluded (slightly-deleterious control)
  * 20 vaccine strains excluded (culture-passaged; private nonsyn variants)
  * Cochran-Mantel-Haenszel stratified by gene, not pooled counts

PRIMARY (polarized to the LSDV branch; ancestral = codon on which SPPV and
GTPV agree; valid because SPPV+GTPV are a clade, bootstrap 100):
  14,412 / 16,215 codons had a confident ancestral state (88.9%)
  Pn=106.5  Ps=299.5  Dn=74.8  Ds=83.2
  NI = 0.396,  alpha = 0.604,  CMH p = 0.0090,  MH odds ratio 1.98

FIVE INDEPENDENT ROBUSTNESS CHECKS, ALL PASSED:
  1. two outgroups agree: alpha 0.491 (SPPV) and 0.540 (GTPV), unpolarized
  2. removing the 53 SRA-derived genomes: alpha changes by 0.002
  3. excluding vaccine strains: alpha rises and becomes cutoff-independent
     (identical at 5% and 10%), showing the earlier cutoff-dependence was
     entirely a vaccine-clade effect
  4. polarization removes ~78% of fixed differences (they sit on the
     SPPV/GTPV branches) and alpha still rises to 0.604
  5. masking gap-rich alignment columns: alpha 0.604 -> 0.577 (>5% gaps)
     -> 0.570 (>1% gaps); Dn drops by only 5

  DoS > 0 in 33/48 genes vs SPPV (Wilcoxon p=0.0014); polarized 20/32
  (Wilcoxon p=0.013, binomial p=0.22)

REPORTING RULES:
  * This is a GENOME-WIDE result. No single gene passes FDR.
  * alpha is a LOWER bound: segregating slightly deleterious variants and
    the conservative treatment of 239 polymorphic-and-divergent sites both
    push it down.
  * Only 9 of 60 genes have >=5 fixed differences on the LSDV branch, and
    10 genes have none at all. Divergence on this branch is sparse.
  * Partially reuses the divergence counts of stage 10; not fully independent.
