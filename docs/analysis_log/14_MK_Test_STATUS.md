# 14 McDonald-Kreitman — FINAL, primary selection result

Rationale: within-LSDV dS is ~0, so any dN/dS ratio is undefined. MK avoids
this by contrasting LSDV polymorphism with divergence to an outgroup, where
divergence is ~3% and the denominator is well estimated.

## Design decisions (implemented in `scripts/`)
  * two independent outgroups run separately (SPPV n=34, GTPV n=16)
  * ancestral state = outgroup majority codon, >= 90% consensus required
  * sites polymorphic in LSDV are counted as polymorphism only (conservative)
  * Nei-Gojobori pathway averaging for codons differing at >1 position
  * variants below 5% frequency excluded (slightly-deleterious control)
  * vaccine strains (n=20) excluded from the primary analysis
  * Cochran-Mantel-Haenszel stratified by gene, not pooled counts

## PRIMARY RESULT — polarised to the LSDV branch
Ancestral = codon on which SPPV and GTPV agree; valid because SPPV+GTPV form
a clade with bootstrap 100. Ingroup n=227 (240 LSDV minus 20 vaccine).

  14,412 of 16,215 codons had a confident ancestral state (88.9%)
  Pn=106.5  Ps=299.5  Dn=74.8  Ds=83.2      NI=0.396
  alpha = 0.604   bootstrap 95% CI [0.33, 0.75]
  alpha_MH = 0.496  RBG 95% CI [0.18, 0.69]
  CMH p = 8.96e-03   OR_MH = 1.98 (95% CI 1.22-3.24, 56 non-empty strata)

Polarisation removes ~78% of fixed differences (Dn 343.1 -> 74.8): most
divergence in the unpolarised version sits on the SPPV/GTPV branches.
17 of 60 genes have zero fixed differences on the LSDV branch. This is why
the CI is wide; the lower bound 0.33 is the conservative reading.

## CONFIDENCE INTERVALS
Gene-level bootstrap, B=10,000, seed=20260817 (`scripts/mk_bootstrap_ci.py`).
Full table: `03_summary/alpha_confidence_intervals.csv`.
Per-replicate draws retained as `03_summary/*_alpha_bootstrap.npy`.

| analysis | alpha | bootstrap 95% CI | alpha_MH | RBG 95% CI |
|---|---|---|---|---|
| polarised (primary)    | 0.604 | 0.33 - 0.75  | 0.496 | 0.18 - 0.69 |
| no vaccine vs GTPV     | 0.540 | 0.35 - 0.66  | 0.485 | 0.30 - 0.62 |
| no vaccine vs SPPV     | 0.491 | 0.32 - 0.61  | 0.480 | 0.33 - 0.60 |
| gap-mask >5%           | 0.577 | 0.29 - 0.72  | 0.477 | 0.14 - 0.68 |
| gap-mask >1%           | 0.570 | 0.28 - 0.72  | 0.453 | 0.10 - 0.67 |
| no vaccine+no SRA GTPV | 0.442 | 0.14 - 0.63  | 0.444 | 0.25 - 0.59 |
| all LSDV vs GTPV       | 0.435 | 0.12 - 0.62  | 0.429 | 0.23 - 0.58 |
| no SRA vs GTPV         | 0.434 | 0.12 - 0.62  | 0.432 | 0.24 - 0.58 |
| no vaccine+no SRA SPPV | 0.399 | 0.12 - 0.58  | 0.421 | 0.26 - 0.55 |
| all LSDV vs SPPV       | 0.393 | 0.12 - 0.57  | 0.409 | 0.24 - 0.54 |
| no SRA vs SPPV         | 0.391 | 0.12 - 0.57  | 0.409 | 0.24 - 0.54 |
| no cutoff vs GTPV      | 0.228 | -0.13 - 0.45 | 0.110 | -0.16 - 0.32 |
| no cutoff vs SPPV      | 0.179 | -0.15 - 0.40 | 0.118 | -0.10 - 0.29 |

Bootstrap CIs vary in the third decimal between full-pipeline runs because
each file draws from a shared generator in glob order; the seed fixes the
stream, not the per-file allocation. Quote to two decimals.

## ROBUSTNESS — all seven tests pass
  1. two independent outgroups: alpha 0.491 (SPPV) and 0.540 (GTPV)
  2. removing the 53 SRA-derived genomes: alpha shifts by 0.002 / 0.001
  3. removing vaccine strains: alpha rises AND removes cutoff dependence.
     With vaccine: 0.393 (5%) -> 0.495 (10%) vs SPPV.
     Without:     0.491 (5%) -> 0.491 (10%). Identical.
  4. polarisation: removes ~78% of fixed differences, alpha still rises
  5. gap masking: 306/16,215 codons (>5%) and 614/16,215 (>1%) removed;
     alpha 0.604 -> 0.577 -> 0.570. Displacement << CI width.
     MAXGAP=1.00 control reproduces 0.604 exactly.
  6. APOBEC3 context: no enrichment in polymorphism or divergence.
     TC->TT 25/79 polymorphic (expected 0.319, 0.99x, p=0.56), 0/16 fixed;
     GA->AA 25/81 polymorphic (expected 0.405, 0.76x), 3/17 fixed.
     C>T + G>A = 14.7% of fixed differences; A>G + T>C = 50.9%.
     No masking applied. See `03_apobec_diagnostic/`.
  7. two independent implementations (`mk_test_v3_novaccine.py` and
     `mk_sensitivity.py`) give identical alpha (0.491 / 0.540). The duplicate
     output files are retained deliberately as a cross-validation record.

  DoS sign test: positive in 33/48 genes vs SPPV (Wilcoxon p=0.0014,
  binomial p=0.0133); polarised 20/32 (Wilcoxon p=0.0128, binomial p=0.2153).

## GENE COUNTS — 60 / 41 / 32 / 56
  60 = all single-copy core orthogroups; basis of the genome-wide alpha
  41 = genes with a defined DoS (17 have Dn+Ds=0, 2 have Pn+Ps=0)
  32 = genes with DoS != 0; the sign test runs on these only
  56 = non-empty strata in the Mantel-Haenszel calculation
  All four must appear with their definitions in the manuscript.

## REPORTING RULES
  * this is a GENOME-WIDE result; no single gene passes FDR
  * alpha is a floor, not a ceiling
  * only 9 of 60 genes have >= 5 fixed differences on the LSDV branch
  * partly reuses the same divergence counts as stage 10 -> not independent
  * phrase it as "fixed amino-acid differences", never "amino-acid
    substitutions" -- alpha concerns divergence, not polymorphism

## FILE MAP (`02_results/`)
  mk_polarized_LSDV_branch_allgenes.csv  <- PRIMARY (60 genes)
  mk_gapmask{100,05,01}_allgenes.csv     <- gap-mask series (100 = control)
  mk_v3_novacc_vs_{SPPV,GTPV}_allgenes.csv
  mk_subset_{all,novacc,nosra,noboth}_vs_{SPPV,GTPV}_allgenes.csv
  mk_v2_vs_{SPPV,GTPV}_allgenes.csv      <- with vaccine, 5% cutoff
  mk_v1_nocutoff_vs_{SPPV,GTPV}.csv      <- no frequency cutoff
Superseded 41- and 55-gene tables: `99_Backups/mk_partial_gene_tables/`.
