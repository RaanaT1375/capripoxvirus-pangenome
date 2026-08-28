# 09 Temporal Signal — FINAL (revised)

USE: `root_to_tip_LSDV_only.csv`, `regression_stats_LSDV_only.json`
SUPERSEDED: `root_to_tip_data.csv`, `regression_stats.json` (all 3 species, n=234)

The original run mixed LSDV with SPPV/GTPV. Root-to-tip distance for the LSDV
reference is 0.000003 vs 0.034 (SPPV) and 0.032 (GTPV) — a ~10^4 difference, so
between-species variance swamped any within-LSDV signal.

Revised result (LSDV only, n=205, 1954-2025): a significant NEGATIVE slope
(R2=0.111, p=1e-06) driven entirely by sampling — all 20 pre-2000 genomes are
African, while Asia (120) and Europe (44) are all post-2000. Restricting to
2000+ gives slope ~ 0 (n=185, R2=0.0001, p=0.887). ANOVA continent x distance:
F=6.29, p=0.0022; no within-continent regression is significant.

CONCLUSION: no temporal signal. BEAST / TMRCA not attempted.
