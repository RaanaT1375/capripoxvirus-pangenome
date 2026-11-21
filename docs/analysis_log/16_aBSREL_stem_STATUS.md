# 16 aBSREL on the LSDV stem branch — NEGATIVE (one artefact, well documented)

Design: label only the edge separating LSDV from the SPPV+GTPV clade
({Stem}), then `hyphy absrel --branches Stem`. This is the same edge that
polarized MK attributes divergence to, tested by likelihood instead of
counting. Run as a 60-task array (job 647986, ~32 min).

Note: the four genes that fail everywhere (OG0000070, 075, 091, 121) needed
`hyphy absrel ENV=TOLERATE_NUMERICAL_ERRORS=1 ...` — as a command-line
argument. `export TOLERATE_NUMERICAL_ERRORS=1` does NOT work.

RESULTS: 49 genes testable, 2 significant. Both are rejected:
  * OG0000067 — baseline omega and the positive class both pinned at 1e10
    with weight 1.0, i.e. dS = 0 across the whole gene. Degenerate fit.
  * OG0000120 (Ig domain OX-2-like) — see below.

Eleven genes could not be tested at all: they have zero fixed differences on
the LSDV branch (independently confirmed from the MK counts). Median
baseline omega on the stem is 0.147 — purifying selection dominates even
across the species split.

THE OG0000120 STORY (worth keeping; it is the methodological centrepiece):
  Model-free evidence looked strong. Under the plain nucleotide GTR model
  the LSDV stem is 0.105217 subs/site, the single longest branch of 300,
  holding 43.7% of the whole tree; the median across 49 genes is 0.082 and
  the runner-up is 0.176. Genome-wide the LSDV stem is the SHORTEST of the
  three (0.00926 vs GTPV 0.01333 vs SPPV 0.02162), so this gene inverts the
  genomic pattern.
  aBSREL: baseline omega 0.590, a positively selected class at omega=697.6
  with weight 0.141, LRT=103.5, corrected p=0.
  Specificity control: the SPPV and GTPV stems give a single rate class and
  LRT=0 (p=1.0) — though those branches are 5-6x shorter, so that control
  has less power.

  AND YET IT IS AN ARTEFACT. All 14 fixed differences sit in codons
  174-202 of 217 (KS D=0.802, p<1e-9; window histogram
  [0,0,0,0,0,0,0,0,10,4]). The weight of the "positively selected" class
  (0.141) matches the block's share of the gene (29/217 = 0.134) almost
  exactly: the rate class IS the block. That block coincides with the
  gappiest part of the alignment (0% gaps in codons 20-160, 29% in 180-200,
  36% in 200-217), where LSDV carries a 7-residue indel and a divergent
  hydrophobic C-terminal tail. dN/dS cannot be estimated across an
  indel-rich region.
  (Codons 0-20 are 93.5% gaps — the ORF start is annotated inconsistently
  across sequences. Worth mentioning in Methods.)

  This clustering is NOT a general problem: of the 9 genes with >=5 fixed
  differences, only OG0000120 is significantly clustered after FDR. The MK
  result is therefore unaffected.

BEST REMAINING GENE-LEVEL CANDIDATE (for Discussion, not as a result):
  OG0000118, LSDV135 putative IFN-alpha/beta binding protein. Elevated in
  MK vs GTPV (rank 3, p=0.0085), in RELAX (K=4.67, p=5e-06, clean fit), and
  its divergence is dispersed rather than clustered (span 0.795, KS q=0.135)
  — the opposite of OG0000120. It does not pass FDR in MK on its own.
