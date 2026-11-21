#!/usr/bin/env Rscript
# =====================================================================
# Combined Figure 1 (Panels A-E) - Exact Original Style (Serif)
# =====================================================================
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ape)
  library(dplyr)
  library(tidyr)
  library(circlize)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

FONT_FAM  <- "serif"
out_dir   <- "revision_2026-11/final_figures"
tree_file <- "06_Recombination/02_Gubbins_Results/lsdv_gubbins.final_tree_renamed.nwk"
gc_file   <- "03_OrthoFinder_Results/OrthoFinder_Results_v3/Results_Aug04/Orthogroups/Orthogroups.GeneCount.tsv"
cat_file  <- "04_Pangenome_Statistics/pangenome_categories_v3.tsv"
rare_file <- "04_Pangenome_Statistics/rarefaction_summary_v3.tsv"
meta_file <- "00_Metadata/species_assignment.csv"

# =========================================================
# PANEL A: ComplexHeatmap + Dendrogram (captured as grob)
# =========================================================
tree <- read.tree(tree_file)
tree_ultra <- compute.brlen(tree)
tree_order <- tree_ultra$tip.label

counts <- read.delim(gc_file, check.names = FALSE)
rownames(counts) <- counts$Orthogroup
counts$Orthogroup <- NULL
if ("Total" %in% colnames(counts)) counts$Total <- NULL

mat <- t(as.matrix(counts) > 0) + 0
common_genomes <- intersect(rownames(mat), tree_order)
mat <- mat[common_genomes, ]
mat <- mat[tree_order[tree_order %in% common_genomes], ]

categories <- read.delim(cat_file) %>% arrange(desc(fraction))
common_genes <- intersect(colnames(mat), categories$Orthogroup)
mat <- mat[, common_genes]

meta_df <- read.csv(meta_file, stringsAsFactors = FALSE)
colnames(meta_df)[1:2] <- c("Genome", "Species")
meta_df$Genome  <- trimws(meta_df$Genome)
meta_df$Species <- trimws(meta_df$Species)

clean_rownames <- trimws(rownames(mat))
lineage_vector <- meta_df$Species[match(clean_rownames, meta_df$Genome)]
lineage_vector[is.na(lineage_vector)] <- "Unknown"

heat_colors    <- structure(c("#FFFFFF", "#002F57"), names = c("0", "1"))
species_colors <- c("LSDV" = "#90A0BC", "SPPV" = "#CB9A9F", "GTPV" = "#585D26", "Unknown" = "#871746")
cat_colors     <- c("Core (>=99%)" = "#90A0BC", "Soft-core (95-99%)" = "#585D26",
                    "Shell (15-95%)" = "#871746", "Cloud (<15%)" = "#CB9A9F")

left_ann <- rowAnnotation(
  Lineage = lineage_vector,
  col = list(Lineage = species_colors),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.38, "cm"),
  border = TRUE,
  annotation_legend_param = list(Lineage = list(
    title_gp  = gpar(fontface = "bold", fontfamily = FONT_FAM, fontsize = 11),
    labels_gp = gpar(fontfamily = FONT_FAM, fontsize = 10)
  ))
)

top_ann <- HeatmapAnnotation(
  Compartment = categories$category[match(colnames(mat), categories$Orthogroup)],
  col = list(Compartment = cat_colors),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.28, "cm"),
  border = TRUE,
  annotation_legend_param = list(Compartment = list(
    title     = "Pan-genome",
    title_gp  = gpar(fontface = "bold", fontfamily = FONT_FAM, fontsize = 11),
    labels_gp = gpar(fontfamily = FONT_FAM, fontsize = 10)
  ))
)

ht <- Heatmap(
  mat,
  name = "Presence",
  col = heat_colors,
  cluster_rows = as.dendrogram(as.hclust(tree_ultra)),
  row_dend_width = unit(3.0, "cm"),
  row_dend_gp = gpar(lwd = 0.7),
  left_annotation = left_ann,
  top_annotation = top_ann,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  border = TRUE,
  use_raster = TRUE,
  raster_quality = 5,
  show_heatmap_legend = FALSE
)

# Capture ComplexHeatmap cleanly without stray Rplots.pdf
pdf(NULL, width = 11, height = 5.5, family = FONT_FAM)
ht_grob <- grid.grabExpr(draw(ht, merge_legend = TRUE, padding = unit(c(2, 4, 2, 4), "mm")))
dev.off()
p_A <- wrap_elements(full = ht_grob)

# =========================================================
# PANEL B: Accumulation Curve
# =========================================================
rare <- read.delim(rare_file)
has_sd <- all(c("pan_sd", "core_sd") %in% colnames(rare))
rare_long <- rare %>%
  pivot_longer(cols = c(pan_mean, core_mean), names_to = "Metric", values_to = "Size") %>%
  mutate(
    SD = if (has_sd) ifelse(Metric == "pan_mean", rare$pan_sd[match(N, rare$N)], rare$core_sd[match(N, rare$N)]) else 0,
    Metric = ifelse(Metric == "pan_mean", "Pan-genome", "Core Genome")
  )

curve_colors <- c("Pan-genome" = "#871746", "Core Genome" = "#90A0BC")

p_B <- ggplot(rare_long, aes(x = N, y = Size, color = Metric, fill = Metric)) +
  { if (has_sd) geom_ribbon(aes(ymin = Size - SD, ymax = Size + SD), alpha = 0.2, color = NA) } +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = curve_colors) +
  scale_fill_manual(values = curve_colors) +
  theme_classic(base_size = 12, base_family = FONT_FAM) +
  labs(x = "Number of Genomes", y = "Number of Orthogroups") +
  annotate("text", x = 150, y = 188,
           label = "italic(n) == 140.71 %*% italic(N)^0.0679",
           parse = TRUE, family = FONT_FAM, size = 4.2, color = "black") +
  theme(
    axis.title         = element_text(face = "bold", size = 12),
    axis.text          = element_text(color = "black", size = 10.5),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 10.5, face = "bold"),
    legend.background  = element_rect(fill = "transparent", color = NA),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    plot.margin        = margin(8, 10, 4, 10)
  )

# =========================================================
# PANEL C: Frequency Spectrum
# =========================================================
categories$category <- factor(categories$category,
                              levels = c("Cloud (<15%)", "Shell (15-95%)", "Soft-core (95-99%)", "Core (>=99%)"))

p_C <- ggplot(categories, aes(x = n_genomes_present, fill = category)) +
  geom_histogram(binwidth = 3, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = cat_colors) +
  theme_bw(base_size = 12, base_family = FONT_FAM) +
  labs(x = "Number of Genomes", y = "Number of Gene Clusters") +
  scale_x_continuous(breaks = seq(0, 300, by = 50)) +
  theme(
    axis.title        = element_text(face = "bold", size = 12),
    axis.text         = element_text(color = "black", size = 10.5),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text       = element_text(size = 8.5, face = "bold"),
    legend.key.size   = unit(0.35, "cm"),
    panel.grid.major  = element_line(color = "grey90"),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", linewidth = 0.9),
    plot.margin       = margin(8, 10, 4, 10)
  )

# =========================================================
# PANEL D: Donut Chart
# =========================================================
pg_data <- data.frame(
  Category = c("Core (>=99%)", "Soft-core (95-99%)", "Shell (15-95%)", "Cloud (<15%)"),
  Count    = c(132, 10, 14, 50)
) %>%
  mutate(
    Fraction   = Count / sum(Count),
    Percentage = round(Fraction * 100, 1),
    Category   = factor(Category, levels = c("Cloud (<15%)", "Shell (15-95%)", "Soft-core (95-99%)", "Core (>=99%)"))
  )

p_D <- ggplot(pg_data, aes(x = 2, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "white", linewidth = 1.1) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(values = cat_colors) +
  xlim(0.5, 2.5) +
  theme_void(base_family = FONT_FAM) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            color = "white", family = FONT_FAM, fontface = "bold", size = 3.2) +
  annotate("text", x = 0.5, y = 0, label = "Total\n206",
           family = FONT_FAM, fontface = "bold", size = 5.5, color = "black") +
  theme(
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 10.5, face = "bold"),
    legend.key.size = unit(0.45, "cm"),
    plot.margin     = margin(4, 10, 4, 4)
  )

# =========================================================
# PANEL E: PCA + Marginal Densities
# =========================================================
mat_t <- t(as.matrix(counts) > 0) + 0
gene_vars <- apply(mat_t, 2, var)
mat_filtered <- mat_t[, gene_vars > 0]

pca_res <- prcomp(mat_filtered, scale. = TRUE)
pca_df  <- as.data.frame(pca_res$x)
pca_var <- summary(pca_res)$importance[2, ]
pc1_label <- paste0("PC1 (", round(pca_var[1] * 100, 1), "% Variance)")
pc2_label <- paste0("PC2 (", round(pca_var[2] * 100, 1), "% Variance)")

pca_df$Genome_clean <- trimws(rownames(pca_df))
pca_df <- left_join(pca_df, meta_df, by = c("Genome_clean" = "Genome"))
pca_df$Species[is.na(pca_df$Species)] <- "Unknown"

p_main <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Species, fill = Species)) +
  geom_point(size = 2.6, alpha = 0.85, shape = 16) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  theme_bw(base_size = 12, base_family = FONT_FAM) +
  labs(x = pc1_label, y = pc2_label) +
  theme(
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(color = "black", size = 10.5),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.9)
  )

p_top <- ggplot(pca_df, aes(x = PC1, fill = Species, color = Species)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  theme_void() + theme(legend.position = "none")

p_right <- ggplot(pca_df, aes(x = PC2, fill = Species, color = Species)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  coord_flip() + theme_void() + theme(legend.position = "none")

layout_E <- "
  A#
  BC
  DD
"

p_E_inner <- p_top + p_main + p_right + guide_area() +
  plot_layout(design = layout_E, widths = c(1, 0.22), heights = c(0.22, 1, 0.14), guides = "collect")

p_E <- wrap_elements(full = p_E_inner)

# =========================================================
# COMBINE ALL 5 PANELS INTO UNIFIED FIGURE 1
# =========================================================
design_main <- "
  AA
  BC
  DE
"

fig1 <- p_A + p_B + p_C + p_D + p_E +
  plot_layout(
    design  = design_main,
    heights = c(1.25, 0.92, 1.02),
    widths  = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", family = FONT_FAM))

ggsave(file.path(out_dir, "Figure_1_Final.pdf"), fig1, width = 11, height = 14, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_1_Final.png"), fig1, width = 11, height = 14, dpi = 300)
message("Saved unified Figure_1_Final.pdf and Figure_1_Final.png in: ", out_dir)
