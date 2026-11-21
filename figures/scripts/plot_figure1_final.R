#!/usr/bin/env Rscript
# =====================================================================
# Figure 1 (All 5 Panels + Combined) - Cluster Ready (r_figs)
# Exact Original Palette | Sans-Serif | Strictly External Legends
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(grid)
})

FONT_FAM <- "sans"

# Exact original palette
heat_colors    <- c("0" = "#FFFFFF", "1" = "#002F57")
species_colors <- c("LSDV" = "#90A0BC", "SPPV" = "#CB9A9F", "GTPV" = "#585D26", "Unknown" = "#871746")
cat_colors     <- c(
  "Core (>=99%)"       = "#90A0BC",
  "Soft-core (95-99%)" = "#585D26",
  "Shell (15-95%)"     = "#871746",
  "Cloud (<15%)"       = "#CB9A9F"
)

# Helper: extract tip order from Newick without ape
get_newick_tips <- function(tree_path) {
  txt <- paste(readLines(tree_path, warn = FALSE), collapse = "")
  txt <- gsub(":[0-9eE\\.\\-\\+]+", "", txt)
  txt <- gsub("\\)[^,\\)\\(;]*", ")", txt)
  txt <- gsub("[\\(\\);\\s]", "", txt)
  tips <- unlist(strsplit(txt, ","))
  tips[nzchar(tips)]
}

# --- 1. Load Cluster Data ---
tree_file <- "06_Recombination/02_Gubbins_Results/lsdv_gubbins.final_tree_renamed.nwk"
gc_file   <- "03_OrthoFinder_Results/OrthoFinder_Results_v3/Results_Aug04/Orthogroups/Orthogroups.GeneCount.tsv"
cat_file  <- "04_Pangenome_Statistics/pangenome_categories_v3.tsv"
meta_file <- "00_Metadata/species_assignment.csv"
out_dir   <- "revision_2026-11/final_figures"

counts <- read.delim(gc_file, check.names = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- counts$Orthogroup
counts$Orthogroup <- NULL
if ("Total" %in% colnames(counts)) counts$Total <- NULL

mat_bin <- ifelse(as.matrix(counts) > 0, 1L, 0L) # OGs x Genomes

categories <- read.delim(cat_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  arrange(desc(fraction)) %>%
  mutate(category = factor(category, levels = c(
    "Cloud (<15%)", "Shell (15-95%)", "Soft-core (95-99%)", "Core (>=99%)"
  )))

meta_df <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  rename(Genome = Name, Species = species) %>%
  mutate(
    Genome  = trimws(Genome),
    Species = ifelse(is.na(Species) | !(Species %in% names(species_colors)), "Unknown", trimws(Species))
  )

tree_order <- get_newick_tips(tree_file)
tree_order <- tree_order[tree_order %in% colnames(mat_bin)]
missing_g  <- setdiff(colnames(mat_bin), tree_order)
if (length(missing_g) > 0) tree_order <- c(tree_order, missing_g)

mat_bin <- mat_bin[, tree_order, drop = FALSE]

# =====================================================================
# PANEL A: Phylogenomic-Pangenome Heatmap
# =====================================================================
common_genes <- intersect(categories$Orthogroup, rownames(mat_bin))
sub_mat <- t(mat_bin[common_genes, , drop = FALSE]) # Genomes x OGs

hm_long <- as.data.frame(as.table(sub_mat), stringsAsFactors = FALSE)
colnames(hm_long) <- c("Genome", "Orthogroup", "Present")
hm_long <- hm_long %>%
  mutate(
    Genome     = factor(Genome, levels = tree_order),
    Orthogroup = factor(Orthogroup, levels = common_genes),
    Present    = factor(Present, levels = c("0", "1"))
  )

sp_bar_df <- tibble(Genome = factor(tree_order, levels = tree_order)) %>%
  left_join(meta_df, by = "Genome") %>%
  mutate(Species = factor(ifelse(is.na(Species), "Unknown", Species),
                          levels = c("LSDV", "SPPV", "GTPV", "Unknown")))

p_sp_bar <- ggplot(sp_bar_df, aes(x = 1, y = Genome, fill = Species)) +
  geom_tile() +
  scale_fill_manual(values = species_colors, name = "Lineage") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = "Genomes (n = 290, phylogeny-ordered)", x = NULL) +
  theme_void(base_family = FONT_FAM) +
  theme(
    axis.title.y    = element_text(angle = 90, size = 10, face = "bold", color = "black", margin = margin(r = 4)),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none",
    plot.margin     = margin(0, 2, 0, 2)
  )

og_bar_df <- categories %>%
  filter(Orthogroup %in% common_genes) %>%
  mutate(Orthogroup = factor(Orthogroup, levels = common_genes))

p_og_bar <- ggplot(og_bar_df, aes(x = Orthogroup, y = 1, fill = category)) +
  geom_tile() +
  scale_fill_manual(values = cat_colors, name = "Pan-genome") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void(base_family = FONT_FAM) +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 9.5, face = "bold", color = "black"),
    legend.text      = element_text(size = 9, color = "black"),
    legend.key.size  = unit(0.35, "cm"),
    legend.margin    = margin(b = 3),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin      = margin(0, 0, 2, 0)
  )

p_hm <- ggplot(hm_long, aes(x = Orthogroup, y = Genome, fill = Present)) +
  geom_raster() +
  scale_fill_manual(values = heat_colors, guide = "none") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = paste0("Orthologous Gene Clusters (n = ", length(common_genes), ")"), y = NULL) +
  theme_bw(base_size = 10, base_family = FONT_FAM) +
  theme(
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.7),
    axis.title.x = element_text(face = "bold", margin = margin(t = 4)),
    plot.margin  = margin(0, 2, 2, 0)
  )

p_A <- (plot_spacer() + p_og_bar + plot_layout(widths = c(0.032, 1))) /
       (p_sp_bar + p_hm + plot_layout(widths = c(0.032, 1))) +
       plot_layout(heights = c(0.06, 1))

ggsave(file.path(out_dir, "Figure1_PanelA_Heatmap.pdf"), p_A, width = 10, height = 6)
ggsave(file.path(out_dir, "Figure1_PanelA_Heatmap.png"), p_A, width = 10, height = 6, dpi = 300)

# =====================================================================
# PANEL B: Pan-genome Accumulation Curve (External Legend at Bottom)
# =====================================================================
set.seed(42)
n_gen  <- ncol(mat_bin)
n_perm <- 100
pan_mat  <- matrix(0, nrow = n_perm, ncol = n_gen)
core_mat <- matrix(0, nrow = n_perm, ncol = n_gen)

for (p in seq_len(n_perm)) {
  ord   <- sample.int(n_gen)
  m_ord <- mat_bin[, ord, drop = FALSE]
  cs    <- t(apply(m_ord, 1, cumsum))
  pan_mat[p, ]  <- colSums(cs > 0)
  core_mat[p, ] <- colSums(cs == rep(seq_len(n_gen), each = nrow(m_ord)))
}

rare_long <- rbind(
  data.frame(N = seq_len(n_gen), Size = colMeans(pan_mat),  SD = apply(pan_mat, 2, sd),  Metric = "Pan-genome"),
  data.frame(N = seq_len(n_gen), Size = colMeans(core_mat), SD = apply(core_mat, 2, sd), Metric = "Core Genome")
)

curve_colors <- c("Pan-genome" = "#871746", "Core Genome" = "#90A0BC")

p_B <- ggplot(rare_long, aes(x = N, y = Size, color = Metric, fill = Metric)) +
  geom_ribbon(aes(ymin = Size - SD, ymax = Size + SD), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.3) +
  scale_color_manual(values = curve_colors) +
  scale_fill_manual(values = curve_colors) +
  theme_bw(base_size = 11, base_family = FONT_FAM) +
  labs(x = "Number of Genomes", y = "Number of Orthogroups") +
  annotate("text", x = 150, y = 186,
           label = "italic(n) == 140.71 %*% italic(N)^0.0679",
           parse = TRUE, family = FONT_FAM, size = 4.0, color = "black") +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "black"),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.8)
  )

ggsave(file.path(out_dir, "Figure1_PanelB_Accumulation.pdf"), p_B, width = 6.5, height = 4.8)
ggsave(file.path(out_dir, "Figure1_PanelB_Accumulation.png"), p_B, width = 6.5, height = 4.8, dpi = 300)

# =====================================================================
# PANEL C: Pan-genome Frequency Spectrum (External Legend at Bottom)
# =====================================================================
p_C <- ggplot(categories, aes(x = n_genomes_present, fill = category)) +
  geom_histogram(binwidth = 3, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = cat_colors) +
  theme_bw(base_size = 11, base_family = FONT_FAM) +
  labs(x = "Number of Genomes", y = "Number of Gene Clusters") +
  scale_x_continuous(breaks = seq(0, 300, by = 50)) +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "black"),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 8.5, face = "bold"),
    legend.key.size  = unit(0.35, "cm"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.8)
  )

ggsave(file.path(out_dir, "Figure1_PanelC_Frequency.pdf"), p_C, width = 6.5, height = 4.8)
ggsave(file.path(out_dir, "Figure1_PanelC_Frequency.png"), p_C, width = 6.5, height = 4.8, dpi = 300)

# =====================================================================
# PANEL D: Pan-genome Composition Donut (External Legend at Right)
# =====================================================================
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
            color = "white", family = FONT_FAM, fontface = "bold", size = 3.0) +
  annotate("text", x = 0.5, y = 0, label = "Total\n206",
           family = FONT_FAM, fontface = "bold", size = 5.0, color = "black") +
  theme(
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 9.5, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.margin     = margin(4, 4, 4, 4)
  )

ggsave(file.path(out_dir, "Figure1_PanelD_Donut.pdf"), p_D, width = 6, height = 4.5)
ggsave(file.path(out_dir, "Figure1_PanelD_Donut.png"), p_D, width = 6, height = 4.5, dpi = 300)

# =====================================================================
# PANEL E: PCA with Marginal Densities (External Legend at Bottom)
# =====================================================================
mat_t <- t(mat_bin)
gene_vars <- apply(mat_t, 2, var)
mat_filtered <- mat_t[, gene_vars > 0, drop = FALSE]

pca_res <- prcomp(mat_filtered, scale. = TRUE)
pca_df  <- as.data.frame(pca_res$x)
pca_var <- summary(pca_res)$importance[2, ]
pc1_label <- paste0("PC1 (", round(pca_var[1] * 100, 1), "% Variance)")
pc2_label <- paste0("PC2 (", round(pca_var[2] * 100, 1), "% Variance)")

pca_df$Genome_clean <- trimws(rownames(pca_df))
pca_df <- left_join(pca_df, meta_df, by = c("Genome_clean" = "Genome"))
pca_df$Species[is.na(pca_df$Species)] <- "Unknown"

p_main <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Species, fill = Species)) +
  geom_point(size = 2.4, alpha = 0.85, shape = 16) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  theme_bw(base_size = 11, base_family = FONT_FAM) +
  labs(x = pc1_label, y = pc2_label) +
  theme(
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "black"),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 9.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.8)
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
  coord_flip() +
  theme_void() + theme(legend.position = "none")

layout_E <- "
  A#
  BC
  DD
"

p_E <- p_top + p_main + p_right + guide_area() +
  plot_layout(design = layout_E, widths = c(1, 0.22), heights = c(0.22, 1, 0.14), guides = "collect")

ggsave(file.path(out_dir, "Figure1_PanelE_PCA.pdf"), p_E, width = 6.5, height = 5.8)
ggsave(file.path(out_dir, "Figure1_PanelE_PCA.png"), p_E, width = 6.5, height = 5.8, dpi = 300)

# =====================================================================
# COMBINED FIGURE 1 (All 5 Panels)
# =====================================================================
bottom_grid <- (p_B | p_C) / (p_D | wrap_elements(full = p_E)) +
  plot_layout(heights = c(1, 1.15))

fig1_combined <- wrap_elements(full = p_A) / bottom_grid +
  plot_layout(heights = c(1.15, 1.65)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold", family = FONT_FAM))

ggsave(file.path(out_dir, "Figure_1_Combined.pdf"), fig1_combined, width = 10, height = 12.5)
ggsave(file.path(out_dir, "Figure_1_Combined.png"), fig1_combined, width = 10, height = 12.5, dpi = 300)
message("All 5 panels and Figure_1_Combined saved in: ", out_dir)
