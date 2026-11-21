# Stage 26 — Figure 4A (pan-GWAS volcano), drawn in R.
#
# Reads the parse written by scripts/04_scoary_table.py, so the plot and
# Table 5 cannot disagree.
#
#   export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
#   Rscript R/figure4a.R
#
# Points meeting the convergence criterion are ringed and labelled. Panel B
# of Figure 4 is a presence/absence heatmap of named orthogroups and is
# unchanged, so it is not redrawn here.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ggrepel)
})
.sd <- local({ a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd() })
source(file.path(.sd, "00_theme.R"))

P <- stage_paths()

# Nothing is withdrawn: the clean per-trait parse reproduces the
# submitted 92 pairs / 16 traits / 18 orthogroups. Add IDs here only
# if an association is dropped in a later revision.
WITHDRAWN <- character(0)
MINPAIRS  <- 3

# Trait names as Scoary actually writes them, so that nothing but a genuinely
# unrecognised trait falls into "Other". Matching is by prefix, so
# Host_Cattle, BAPS_cluster_3 and the rest land in the right family.
categorise <- function(trait) {
  t <- as.character(trait)
  ifelse(t %in% c("Africa", "Asia", "Europe"),               "Geography",
  ifelse(grepl("^Host_", t),                                 "Host",
  ifelse(grepl("^BAPS_cluster", t),                          "Population structure",
  ifelse(t %in% c("Ancestral_recombinant", "Independent_recombinant",
                  "No_recombination", "Recombinant"),        "Recombination",
  ifelse(t %in% c("Vaccine", "Wild", "Wild_type"),           "Clinical",
  ifelse(t %in% c("LSDV", "SPPV", "GTPV"),                   "Species",
                                                             "Other"))))))
}

df <- read.csv(need(file.path(P$summary, "scoary_significant.csv")))
if ("run" %in% names(df)) df <- df[df$run == "primary_290", ]
message(sprintf("%d significant pairs before withdrawal", nrow(df)))

pick <- function(d, ...) {
  for (n in c(...)) if (n %in% names(d)) return(n)
  stop("none of these columns found: ", paste(c(...), collapse = ", "))
}
c_bonf <- pick(df, "Bonferroni_p", "Bonferroni_p_value")
c_or   <- pick(df, "Odds_ratio", "Odds_Ratio")
c_sup  <- pick(df, "Max_supporting_pairs", "Best_supporting_pairs")
c_opp  <- pick(df, "Max_opposing_pairs", "Worst_opposing_pairs")

gone <- df[df$Gene %in% WITHDRAWN, ]
if (nrow(gone)) {
  message(sprintf("withdrawing %d rows for %s", nrow(gone),
                  paste(sort(unique(gone$Gene)), collapse = ", ")))
}
df <- df[!df$Gene %in% WITHDRAWN, ]
message(sprintf("%d pairs plotted, %d traits, %d orthogroups",
                nrow(df), dplyr::n_distinct(df$Trait),
                dplyr::n_distinct(df$Gene)))

orv <- df[[c_or]]
finite_max <- suppressWarnings(max(orv[is.finite(orv)], na.rm = TRUE))
if (!is.finite(finite_max)) finite_max <- 100
cap <- finite_max * 2
df$OR_infinite <- is.infinite(orv) | is.na(orv)
orv[df$OR_infinite] <- cap
orv <- pmax(orv, 1e-4)

df <- df |>
  mutate(
    log2OR    = log2(orv),
    neglog10p = -log10(pmax(.data[[c_bonf]], 1e-300)),
    category  = categorise(Trait),
    converged = .data[[c_sup]] >= MINPAIRS &
                (.data[[c_sup]] - .data[[c_opp]]) > 0,
    label     = ifelse(converged,
                       ifelse(OR_infinite,
                              paste0(Gene, "  (OR → ∞)"), Gene), NA)
  )
df$category <- factor(df$category,
                      levels = names(CATEGORY_COLS)[
                        names(CATEGORY_COLS) %in% df$category])
write.csv(df, file.path(P$results, "fig4_panelA_points.csv"), row.names = FALSE)
message(sprintf("%d points meet the convergence criterion", sum(df$converged)))

lab <- df[df$converged, ] |> group_by(Gene) |> slice_max(neglog10p, n = 1) |>
  ungroup()

p <- ggplot(df, aes(log2OR, neglog10p)) +
  geom_vline(xintercept = 0, colour = "grey75", linetype = "dotted",
             linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), colour = "grey75",
             linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(fill = category), shape = 21, size = 2.6, stroke = 0,
             alpha = 0.78) +
  geom_point(data = df[df$converged, ], shape = 21, size = 3.8,
             fill = NA, colour = INK, stroke = 0.6) +
  geom_text_repel(data = lab, aes(label = label), family = "serif",
                  size = 2.7, colour = INK, min.segment.length = 0,
                  segment.colour = "grey60", segment.size = 0.3,
                  box.padding = 0.45, max.overlaps = 20, seed = 1) +
  annotate("text", x = -Inf, y = -log10(0.05), hjust = -0.05, vjust = -0.6,
           label = "Bonferroni~italic(p)==0.05", parse = TRUE,
           family = "serif", size = 2.6, colour = "grey45") +
  scale_fill_manual(values = CATEGORY_COLS, name = "Phenotype",
                    drop = TRUE) +
  labs(x = expression("Effect size ("*log[2]*" odds ratio)"),
       y = expression(-log[10]*" Bonferroni-adjusted "*italic(p))) +
  theme_cpx()

if (any(df$OR_infinite)) {
  p <- p + labs(caption = paste(
    "Points marked OR \u2192 \u221E are plotted at the right-hand limit: the",
    "orthogroup is present in every genome carrying the phenotype.",
    "Points at the left-hand limit have an odds ratio of zero."))
  p <- p + theme(plot.caption = element_text(size = 7.5, colour = "grey40",
                                             hjust = 0, family = "serif"))
}

save_fig(p, file.path(P$figures, "Figure_4A"), width = 7.4, height = 5.6)
