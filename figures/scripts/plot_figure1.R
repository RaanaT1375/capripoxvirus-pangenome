#!/usr/bin/env Rscript
# =====================================================================
# Figure 1: Pan-genome architecture and population structure
# Soft Muted Q1 Palette | Sans-Serif | Strictly External Legends
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(ggrepel)
  library(grid)
})

source("revision_2026-11/final_figures/scripts/theme_common.R")

get_newick_tips <- function(tree_path) {
  txt <- paste(readLines(tree_path, warn = FALSE), collapse = "")
  txt <- gsub(":[0-9eE\\.\\-\\+]+", "", txt)
  txt <- gsub("\\)[^,\\)\\(;]*", ")", txt)
  txt <- gsub("[\\(\\);\\s]", "", txt)
  tips <- unlist(strsplit(txt, ","))
  tips[nzchar(tips)]
}

# --- 1. Load Data ---
gc_file   <- "03_OrthoFinder_Results/OrthoFinder_Results_v3/Results_Aug04/Orthogroups/Orthogroups.GeneCount.tsv"
cat_file  <- "04_Pangenome_Statistics/pangenome_categories_v3.tsv"
meta_file <- "00_Metadata/species_assignment.csv"
tree_file <- "05_Phylogeny/04_tree/core_genome_ML_tree_v3.treefile"
out_dir   <- "revision_2026-11/final_figures"

gc_raw <- read.delim(gc_file, check.names = FALSE, stringsAsFactors = FALSE)
if ("Total" %in% colnames(gc_raw)) gc_raw$Total <- NULL

og_ids  <- gc_raw$Orthogroup
mat_num <- as.matrix(gc_raw[, -1, drop = FALSE])
rownames(mat_num) <- og_ids
mat_bin <- ifelse(mat_num > 0, 1L, 0L)

cats <- read.delim(cat_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(category = factor(category, levels = c(
    "Core (>=99%)", "Soft-core (95-99%)", "Shell (15-95%)", "Cloud (<15%)"
  )))

meta <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  rename(Genome = Name, Species = species) %>%
  mutate(Species = ifelse(is.na(Species) | !(Species %in% names(PAL_SPECIES)), "Unknown", Species))

tip_order <- get_newick_tips(tree_file)
tip_order <- tip_order[tip_order %in% colnames(mat_bin)]
missing_tips <- setdiff(colnames(mat_bin), tip_order)
if (length(missing_tips) > 0) tip_order <- c(tip_order, missing_tips)
mat_bin <- mat_bin[, tip_order, drop = FALSE]

# --- 2. Panel A: Phylogenomic Presence-Absence Heatmap ---
og_keep <- cats %>%
  filter(n_genomes_present >= 2) %>%
  arrange(desc(n_genomes_present), Orthogroup)

sub_mat <- t(mat_bin[og_keep$Orthogroup, , drop = FALSE])
hm_long <- as.data.frame(as.table(sub_mat), stringsAsFactors = FALSE)
colnames(hm_long) <- c("Genome", "Orthogroup", "Present")
hm_long <- hm_long %>%
  mutate(
    Genome     = factor(Genome, levels = tip_order),
    Orthogroup = factor(Orthogroup, levels = og_keep$Orthogroup),
    Present    = factor(ifelse(Present == 1, "Present", "Absent"), levels = c("Present", "Absent"))
  )

sp_bar_df <- tibble(Genome = factor(tip_order, levels = tip_order)) %>%
  left_join(meta, by = "Genome") %>%
  mutate(Species = factor(ifelse(is.na(Species), "Unknown", Species),
                          levels = c("LSDV", "SPPV", "GTPV", "Unknown")))

p_sp_bar <- ggplot(sp_bar_df, aes(x = 1, y = Genome, fill = Species)) +
  geom_tile() +
  scale_fill_manual(values = PAL_SPECIES, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = "Genomes (n = 290, phylogeny-ordered)", x = NULL) +
  theme_void(base_family = "sans") +
  theme(
    axis.title.y = element_text(angle = 90, size = 9, face = "bold", color = "#222222", margin = margin(r = 5)),
    panel.border = element_rect(color = "#444444", fill = NA, linewidth = 0.5),
    plot.margin  = margin(0, 2, 0, 2)
  )

og_bar_df <- og_keep %>%
  mutate(Orthogroup = factor(Orthogroup, levels = og_keep$Orthogroup))

p_og_bar <- ggplot(og_bar_df, aes(x = Orthogroup, y = 1, fill = category)) +
  geom_tile() +
  scale_fill_manual(values = PAL_PANGENOME, name = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void(base_family = "sans") +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 8.5, face = "bold", color = "#333333"),
    legend.key.size  = unit(0.35, "cm"),
    legend.margin    = margin(b = 4),
    panel.border     = element_rect(color = "#444444", fill = NA, linewidth = 0.5),
    plot.margin      = margin(0, 0, 2, 0)
  )

p_hm <- ggplot(hm_long, aes(x = Orthogroup, y = Genome, fill = Present)) +
  geom_raster() +
  scale_fill_manual(values = c("Present" = "#3D5A80", "Absent" = "#FAF9F6"), name = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = paste0("Orthologous gene clusters (n = ", nrow(og_keep), "; singletons excluded)"), y = NULL) +
  theme_q1_sans(base_size = 9) +
  theme(
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    legend.position = "none",
    axis.title.x    = element_text(margin = margin(t = 4)),
    plot.margin     = margin(0, 2, 2, 0)
  )

panel_A <- (plot_spacer() + p_og_bar + plot_layout(widths = c(0.032, 1))) /
           (p_sp_bar + p_hm + plot_layout(widths = c(0.032, 1))) +
           plot_layout(heights = c(0.06, 1))

# --- 3. Panel B: Accumulation Curves (Legend strictly at bottom) ---
set.seed(42)
n_gen  <- ncol(mat_bin)
n_perm <- 100
pan_mat  <- matrix(0, nrow = n_perm, ncol = n_gen)
core_mat <- matrix(0, nrow = n_perm, ncol = n_gen)

for (p in seq_len(n_perm)) {
  ord <- sample.int(n_gen)
  m_ord <- mat_bin[, ord, drop = FALSE]
  cs <- t(apply(m_ord, 1, cumsum))
  pan_mat[p, ]  <- colSums(cs > 0)
  core_mat[p, ] <- colSums(cs == rep(seq_len(n_gen), each = nrow(m_ord)))
}

acc_df <- tibble(
  N         = seq_len(n_gen),
  pan_mean  = colMeans(pan_mat),
  pan_sd    = apply(pan_mat, 2, sd),
  core_mean = colMeans(core_mat),
  core_sd   = apply(core_mat, 2, sd)
)

panel_B <- ggplot(acc_df, aes(x = N)) +
  geom_ribbon(aes(ymin = pan_mean - pan_sd, ymax = pan_mean + pan_sd),
              fill = "#E07A5F", alpha = 0.20) +
  geom_ribbon(aes(ymin = core_mean - core_sd, ymax = core_mean + core_sd),
              fill = "#3D5A80", alpha = 0.20) +
  geom_line(aes(y = pan_mean, color = "Pan-genome (Total)"), linewidth = 1.1) +
  geom_line(aes(y = core_mean, color = "Strict core"), linewidth = 1.1) +
  annotate("text", x = 155, y = 183,
           label = expression(italic(n) == 140.71 %*% italic(N)^{0.0679}),
           size = 3.2, family = "sans", color = "#222222") +
  scale_color_manual(
    values = c("Pan-genome (Total)" = "#E07A5F", "Strict core" = "#3D5A80"),
    name = NULL
  ) +
  labs(x = "Number of genomes", y = "Number of orthogroups") +
  theme_q1_sans(base_size = 9) +
  theme(legend.position = "bottom")

# --- 4. Panel C: Gene Frequency Spectrum (No in-plot legend) ---
panel_C <- ggplot(cats, aes(x = n_genomes_present, fill = category)) +
  geom_histogram(binwidth = 3, boundary = 0, color = "#FAF9F6", linewidth = 0.2) +
  scale_fill_manual(values = PAL_PANGENOME, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Number of genomes present", y = "Number of orthogroups") +
  theme_q1_sans(base_size = 9)

# --- 5. Panel D: Donut Chart (Clean external repelled labels for small slices) ---
donut_df <- cats %>%
  count(category) %>%
  mutate(
    pct      = n / sum(n) * 100,
    ymax     = cumsum(pct),
    ymin     = lag(ymax, default = 0),
    ymid     = (ymin + ymax) / 2,
    is_small = pct < 10,
    short_cat = gsub(" \\(.*\\)", "", as.character(category)),
    lbl_in   = paste0(short_cat, "\n", n, " (", sprintf("%.1f%%", pct), ")"),
    lbl_out  = paste0(short_cat, ": ", n, " (", sprintf("%.1f%%", pct), "%)")
  )

panel_D <- ggplot(donut_df, aes(ymax = ymax, ymin = ymin, xmax = 3.8, xmin = 2.2, fill = category)) +
  geom_rect(color = "white", linewidth = 1.0) +
  geom_text(
    data = filter(donut_df, !is_small),
    aes(x = 3.0, y = ymid, label = lbl_in),
    color = "white", size = 2.8, fontface = "bold", family = "sans", lineheight = 0.95
  ) +
  geom_label_repel(
    data = filter(donut_df, is_small),
    aes(x = 3.8, y = ymid, label = lbl_out),
    nudge_x = 1.0, direction = "y", segment.size = 0.45, segment.color = "#444444",
    size = 2.7, fontface = "bold", family = "sans", color = "#222222",
    fill = "white", label.size = 0, min.segment.length = 0, box.padding = 0.25
  ) +
  annotate("text", x = 0, y = 0, label = paste0("Total\n", sum(donut_df$n)),
           size = 4.3, fontface = "bold", family = "sans", color = "#333333", lineheight = 0.95) +
  coord_polar(theta = "y", clip = "off") +
  xlim(c(0, 5.1)) +
  scale_fill_manual(values = PAL_PANGENOME, guide = "none") +
  theme_void(base_family = "sans") +
  theme(plot.margin = margin(0, 8, 0, 8))

# --- 6. Panel E: PCA of Accessory Genome (Legend strictly at bottom) ---
acc_mat  <- t(mat_bin)
var_cols <- apply(acc_mat, 2, var) > 0
pca_res  <- prcomp(acc_mat[, var_cols, drop = FALSE], center = TRUE, scale. = TRUE)
var_exp  <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

pca_df <- tibble(
  Genome = rownames(pca_res$x),
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2]
) %>%
  left_join(meta, by = "Genome") %>%
  mutate(Species = factor(ifelse(is.na(Species), "Unknown", Species),
                          levels = c("LSDV", "SPPV", "GTPV", "Unknown")))

panel_E <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(size = 2.2, alpha = 0.85, shape = 21, color = "white", stroke = 0.3) +
  scale_fill_manual(values = PAL_SPECIES, name = NULL, drop = TRUE) +
  labs(
    x = paste0("PC1 (", var_exp[1], "% variance)"),
    y = paste0("PC2 (", var_exp[2], "% variance)")
  ) +
  theme_q1_sans(base_size = 9) +
  theme(legend.position = "bottom")

# --- 7. Combine All Panels ---
bottom_grid <- (panel_B | panel_C) / (panel_D | panel_E) +
  plot_layout(heights = c(1, 1.05))

fig1 <- wrap_elements(full = panel_A) / bottom_grid +
  plot_layout(heights = c(1.20, 1.60)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 13, face = "bold", family = "sans", color = "#111111"))

ggsave(file.path(out_dir, "Figure_1.pdf"), fig1, width = 8.6, height = 10.8)
ggsave(file.path(out_dir, "Figure_1.png"), fig1, width = 8.6, height = 10.8, dpi = 300)
message("Saved Q1 Sans-Serif Figure_1.pdf and Figure_1.png successfully.")
