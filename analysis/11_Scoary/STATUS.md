# 11 Scoary — FINAL is the tree-aware run

USE: `scoary_results_tree/`   (run with -n and the core-genome ML tree)
SUPERSEDED: `_SUPERSEDED_scoary_results_notree/`

Tree-aware run: 92 significant gene-trait pairs, 16 traits, 18 orthogroups.
Bonferroni_p was identical in both runs (it comes from Fisher's test and is
tree-independent), which confirms the reruns are consistent;
Max_Pairwise_comparisons changed, which confirms the tree was actually used.

CRITICAL CAVEAT: every BAPS and host trait has n_robust = 0 and
max_pairs = 1 — zero independent evolutionary events. These associations are
entirely phylogenetic confounding. Specifically, Host_Goat is GTPV and
Host_Sheep is SPPV, so those "host-associated genes" are species markers.
Only 11 of 92 pairs have >= 3 supporting pairs. Four pairs have more opposing
than supporting pairs and should be rejected outright.
