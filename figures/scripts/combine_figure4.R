#!/usr/bin/env Rscript
# =====================================================================
# Combined Figure 4 (Panels A & B) - Unified Serif & Q1 Palette
# Panel A: Pan-GWAS Volcano Plot (Bonferroni-adjusted p, External Legend)
# Panel B: Presence/Absence Heatmap of Key Orthogroups (No clipped text)
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(grid)
})

FONT_FAM <- "serif"
out_dir  <- "revision_2026-11/final_figures"

# --- Input paths ---
scoary_file <- "revision_2026-11/26_Figure_Regeneration/03_summary/scoary_significant.csv"
traits_file <- Sys.glob("11_Scoary/01_inputs/*trait*.csv")[1]
gpa_file    <- Sys.glob("11_Scoary/01_inputs/*gene_presence_absence*.csv")[1]
gc_fallback <- "03_OrthoFinder_Results/OrthoFinder_Results_v3/Results_Aug04/Orthogroups/Orthogroups.GeneCount.tsv"

# =========================================================
# PANEL A: Pan-GWAS Volcano Plot
# =========================================================
df <- read.csv(scoary_file, stringsAsFactors = FALSE)
if ("run" %in% names(df)) df <- df[df$run == "primary_290", ]

pick_col <- function(d, ...) {
  for (n in c(...)) if (n %in% names(d)) return(n)
  stop("Missing expected column: ", paste(c(...), collapse = ", "))
}

c_bonf <- pick_col(df, "Bonferroni_p", "Bonferroni_p_value")
c_or   <- pick_col(df, "Odds_ratio", "Odds_Ratio")
c_sup  <- pick_col(df, "Max_supporting_pairs", "Best_supporting_pairs")
c_opp  <- pick_col(df, "Max_opposing_pairs", "Worst_opposing_pairs")

categorise <- function(trait) {
  t <- as.character(trait)
  ifelse(t %in% c("Africa", "Asia", "Europe"),                 "Geography",
  ifelse(grepl("^Host_", t) | t %in% c("Cattle", "Sheep", "Goat", "Wildlife", "Giraffe"), "Host",
  ifelse(grepl("^BAPS_cluster", t),                            "Population structure",
  ifelse(t %in% c("Ancestral_recombinant", "Independent_recombinant",
                  "No_recombination", "Recombinant"),          "Recombination",
  ifelse(t %in% c("Vaccine", "Wild", "Wild_type"),             "Clinical",
  ifelse(t %in% c("LSDV", "SPPV", "GTPV"),                     "Species",
                                                               "Other"))))))
}

orv <- df[[c_or]]
finite_max <- suppressWarnings(max(orv[is.finite(orv)], na.rm = TRUE))
if (!is.finite(finite_max)) finite_max <- 100
cap <- finite_max * 2

df <- df %>%
  mutate(
    OR_infinite = is.infinite(orv) | is.na(orv),
    or_clean    = ifelse(OR_infinite, cap, pmax(orv, 1e-4)),
    log2OR      = log2(or_clean),
    neglog10p   = -log10(pmax(.data[[c_bonf]], 1e-300)),
    category    = categorise(Trait),
    converged   = .data[[c_sup]] >= 3 & (.data[[c_sup]] - .data[[c_opp]]) > 0,
    label       = ifelse(converged,
                         ifelse(OR_infinite, paste0(Gene, " (OR \u2192 \u221E)"), Gene),
                         NA_character_)
  )

# Harmonized palette with Figures 1 & 3
PAL_PHENO <- c(
  "Geography"            = "#002F57",
  "Clinical"             = "#6B4E71",
  "Host"                 = "#CB9A9F",
  "Recombination"        = "#871746",
  "Population structure" = "#90A0BC",
  "Species"              = "#585D26",
  "Other"                = "#B8B8B8"
)

df$category <- factor(df$category, levels = intersect(names(PAL_PHENO), unique(df$category)))

# One clean label per converged orthogroup (at its most significant point)
lab_df <- df %>%
  filter(converged) %>%
  group_by(Gene) %>%
  slice_max(neglog10p, n = 1, with_ties = FALSE) %>%
  ungroup()

x_min <- min(df$log2OR, na.rm = TRUE)
x_max <- max(df$log2OR, na.rm = TRUE)

p_A <- ggplot(df, aes(x = log2OR, y = neglog10p)) +
  geom_vline(xintercept = 0, colour = "grey65", linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), colour = "#871746", linetype = "dashed", linewidth = 0.55) +
  geom_point(aes(fill = category), shape = 21, size = 3.0, colour = "#2B2B2B", stroke = 0.35, alpha = 0.85) +
  geom_point(data = filter(df, converged), shape = 21, size = 4.6, fill = NA, colour = "#111111", stroke = 0.8) +
  geom_text_repel(
    data = lab_df, aes(label = label),
    family = FONT_FAM, size = 3.2, fontface = "bold", colour = "#111111",
    box.padding = 0.6, point.padding = 0.35, min.segment.length = 0,
    segment.colour = "grey45", segment.size = 0.4, max.overlaps = 30, seed = 42
  ) +
  annotate("text", x = x_min + 0.3, y = -log10(0.05) + 1.8, hjust = 0,
           label = "Bonferroni~italic(p)==0.05", parse = TRUE,
           family = FONT_FAM, size = 3.1, colour = "#871746") +
  scale_fill_manual(values = PAL_PHENO, name = NULL, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  theme_bw(base_size = 12, base_family = FONT_FAM) +
  labs(
    x = expression("Effect size (" * log[2] * " odds ratio)"),
    y = expression(-log[10] * " Bonferroni-adjusted " * italic(p)),
    caption = "Ringed points meet the phylogenetic convergence criterion (>= 3 supporting pairs). Points at the right/left limits represent OR -> Inf and OR = 0, respectively."
  ) +
  theme(
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(color = "black", size = 10.5),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 10, face = "bold"),
    legend.key.size  = unit(0.4, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.85),
    plot.caption     = element_text(size = 8.8, hjust = 0, color = "#333333", family = FONT_FAM),
    plot.margin      = margin(8, 12, 6, 12)
  )

# =========================================================
# PANEL B: Presence/Absence Heatmap of 10 Key Orthogroups
# =========================================================
target_ogs <- c(
  "OG0000157", "OG0000155", "OG0000154",
  "OG0000159", "OG0000161", "OG0000156",
  "OG0000152", "OG0000004",
  "OG0000151", "OG0000153"
)

if (!is.na(gpa_file) && file.exists(gpa_file)) {
  gene_pa <- read.csv(gpa_file, stringsAsFactors = FALSE, check.names = FALSE)
  sub_pa  <- gene_pa[match(target_ogs, gene_pa$Gene), , drop = FALSE]
  pa_mat  <- sub_pa[, 15:ncol(sub_pa), drop = FALSE]
  pa_mat[is.na(pa_mat) | pa_mat == ""] <- 0
  pa_mat[pa_mat != 0] <- 1
  pa_matrix <- apply(pa_mat, 2, as.numeric)
  rownames(pa_matrix) <- target_ogs
} else {
  gc_raw <- read.delim(gc_fallback, check.names = FALSE, stringsAsFactors = FALSE)
  if ("Total" %in% colnames(gc_raw)) gc_raw$Total <- NULL
  sub_gc <- gc_raw[match(target_ogs, gc_raw$Orthogroup), -1, drop = FALSE]
  pa_matrix <- ifelse(as.matrix(sub_gc) > 0, 1, 0)
  rownames(pa_matrix) <- target_ogs
}

traits <- read.csv(traits_file, stringsAsFactors = FALSE, check.names = FALSE)
colnames(traits)[1] <- "Genome"

meta <- traits %>%
  mutate(
    Genome = trimws(Genome),
    Geography = case_when(
      Asia == 1   ~ "Asia",
      Europe == 1 ~ "Europe",
      Africa == 1 ~ "Africa",
      TRUE        ~ "Unknown"
    ),
    Condition = case_when(
      Vaccine == 1     ~ "Vaccine",
      Recombinant == 1 ~ "Recombinant",
      Wild == 1        ~ "Wild",
      TRUE             ~ "Unknown"
    )
  ) %>%
  select(Genome, Geography, Condition)

rownames(meta) <- meta$Genome
colnames(pa_matrix) <- trimws(colnames(pa_matrix))

common_genomes <- intersect(colnames(pa_matrix), rownames(meta))
pa_matrix <- pa_matrix[, common_genomes, drop = FALSE]
meta      <- meta[common_genomes, , drop = FALSE]

heat_colors <- structure(c("#FFFFFF", "#002F57"), names = c("0", "1"))
geo_colors  <- c("Asia" = "#871746", "Europe" = "#CB9A9F", "Africa" = "#495C7C", "Unknown" = "#E5E5E5")
cond_colors <- c("Recombinant" = "#1F3A5F", "Vaccine" = "#6B4E71", "Wild" = "#EAE6DF", "Unknown" = "#E5E5E5")

legend_params <- list(
  title_gp  = gpar(fontface = "bold", fontfamily = FONT_FAM, fontsize = 10.5),
  labels_gp = gpar(fontfamily = FONT_FAM, fontsize = 9.5)
)

ha <- HeatmapAnnotation(
  Geography = meta$Geography,
  Condition = meta$Condition,
  col = list(Geography = geo_colors, Condition = cond_colors),
  border = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontfamily = FONT_FAM, fontsize = 10.5, fontface = "bold"),
  annotation_legend_param = list(Geography = legend_params, Condition = legend_params),
  simple_anno_size = unit(0.38, "cm")
)

row_groups <- factor(
  c(rep("Geographic Adaptation", 6), rep("Allelic Variants & Vaccine", 4)),
  levels = c("Geographic Adaptation", "Allelic Variants & Vaccine")
)

ht <- Heatmap(
  pa_matrix,
  name = "Presence",
  col = heat_colors,
  top_annotation = ha,
  column_split = meta$Geography,
  row_split = row_groups,
  column_title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = FONT_FAM),
  row_title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = FONT_FAM),
  row_names_gp = gpar(fontsize = 9.5, fontfamily = FONT_FAM),
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  border = TRUE,
  rect_gp = gpar(col = NA),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = c(legend_params, list(at = c(0, 1), labels = c("Absent", "Present")))
)

# Capture heatmap cleanly with generous right padding so no legend text is clipped
pdf(NULL, width = 10.5, height = 5.6, family = FONT_FAM)
ht_grob <- grid.grabExpr(draw(ht, merge_legend = TRUE, padding = unit(c(4, 6, 4, 18), "mm")))
dev.off()
p_B <- wrap_elements(full = ht_grob)

# =========================================================
# COMBINE PANELS A & B INTO UNIFIED FIGURE 4
# =========================================================
fig4 <- p_A / p_B +
  plot_layout(heights = c(1.15, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", family = FONT_FAM))

ggsave(file.path(out_dir, "Figure_4_Final.pdf"), fig4, width = 10.5, height = 11.8, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_4_Final.png"), fig4, width = 10.5, height = 11.8, dpi = 300)
message("Saved unified Figure_4_Final.pdf and Figure_4_Final.png in: ", out_dir)
