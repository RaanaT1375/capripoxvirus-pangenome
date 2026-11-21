#!/usr/bin/env Rscript
# =====================================================================
# Combined Figure 3 (Panels A, B, C) - Final 4-Ring Version
# Panel A: Rooted Fan Phylogeny + Numeric Bootstrap + 4 Slim Rings
#          (Viral Species, Recombination, Condition, Region)
# Panel B & C: Preserved exactly as approved
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ape)
  library(ggtree)
  library(ggnewscale)
  library(patchwork)
})

FONT_FAM <- "serif"
out_dir  <- "revision_2026-11/final_figures"

tree_file <- "05_Phylogeny/04_tree/core_genome_ML_tree_v3.treefile"
gff_file  <- Sys.glob("06_Recombination/02_Gubbins_Results/*.recombination_predictions.gff")[1]
baps_file <- "08_Population_Structure/fastbaps_annotated.csv"
sp_file   <- "00_Metadata/species_assignment.csv"
cond_file <- "00_Metadata/itol_condition_FINAL.txt"

# Harmonized Palettes
PAL_SPECIES <- c("LSDV" = "#90A0BC", "SPPV" = "#CB9A9F", "GTPV" = "#585D26", "Unknown" = "#871746")

PAL_RECOMB_RING <- c(
  "Ancestral recombinant"   = "#871746",
  "Independent recombinant" = "#D48C70",
  "No major recombination"  = "#DCE3EC"
)

PAL_RECOMB_BAR <- c(
  "Ancestral recombinant"   = "#871746",
  "Independent recombinant" = "#CB9A9F",
  "No major recombination"  = "#90A0BC"
)

PAL_CONDITION <- c(
  "Recombinant" = "#1F3A5F",
  "Vaccine"     = "#6B4E71",
  "Wild"        = "#EAE6DF",
  "Unknown"     = "#E5E5E5"
)

PAL_REGION <- c(
  "Africa"  = "#C29B48",
  "Asia"    = "#3A6B6B",
  "Europe"  = "#6D7B8D",
  "Unknown" = "#E5E5E5"
)

# --- Parse iTOL Condition file ---
cond_lines <- readLines(cond_file, warn = FALSE)
data_idx   <- which(trimws(cond_lines) == "DATA")[1]
cond_raw   <- cond_lines[(data_idx + 1):length(cond_lines)]
cond_raw   <- cond_raw[nzchar(trimws(cond_raw))]
cond_split <- strsplit(cond_raw, "\t")
cond_df <- data.frame(
  genome_id = trimws(sapply(cond_split, `[`, 1)),
  hex       = toupper(trimws(sapply(cond_split, `[`, 2))),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Condition = case_when(
      hex == "#3F482C" ~ "Recombinant",
      hex == "#B69065" ~ "Vaccine",
      hex == "#ECD7BC" ~ "Wild",
      TRUE             ~ "Unknown"
    )
  ) %>%
  select(genome_id, Condition)

# =========================================================
# PANEL A: Rooted Fan Phylogeny + Numeric Bootstrap + 4 Slim Rings
# =========================================================
tree <- read.tree(tree_file)

sp_df <- read.csv(sp_file, stringsAsFactors = FALSE)
colnames(sp_df)[1:2] <- c("genome_id", "Species")
sp_df$genome_id <- trimws(sp_df$genome_id)
sp_df$Species   <- trimws(sp_df$Species)

baps_df <- read.csv(baps_file, stringsAsFactors = FALSE) %>%
  mutate(
    genome_id = trimws(genome_id),
    Recombination = case_when(
      grepl("Ancestral", recombination_status, ignore.case = TRUE)   ~ "Ancestral recombinant",
      grepl("Independent", recombination_status, ignore.case = TRUE) ~ "Independent recombinant",
      TRUE ~ "No major recombination"
    ),
    Region = ifelse(is.na(continent) | continent == "", "Unknown", continent)
  )

# Root cleanly on SPPV + GTPV outgroup clade & ladderize
outgroup_tips <- sp_df$genome_id[sp_df$Species %in% c("SPPV", "GTPV")]
outgroup_tips <- intersect(outgroup_tips, tree$tip.label)
if (length(outgroup_tips) >= 2) {
  tree <- root(tree, outgroup = outgroup_tips, resolve.root = TRUE)
}
tree <- ladderize(tree, right = FALSE)

meta_all <- data.frame(genome_id = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(sp_df, by = "genome_id") %>%
  left_join(baps_df, by = "genome_id") %>%
  left_join(cond_df, by = "genome_id") %>%
  mutate(
    Species       = ifelse(is.na(Species), "Unknown", Species),
    Recombination = ifelse(is.na(Recombination), "No major recombination", Recombination),
    Condition     = ifelse(is.na(Condition), "Unknown", Condition),
    Region        = ifelse(is.na(Region), "Unknown", Region)
  )
rownames(meta_all) <- meta_all$genome_id

# Count descendant tips per internal node to label major branches cleanly
n_tips  <- length(tree$tip.label)
n_nodes <- tree$Nnode
desc_count <- integer(n_tips + n_nodes)
desc_count[1:n_tips] <- 1L
tr_post <- reorder(tree, "postorder")
for (i in seq_len(nrow(tr_post$edge))) {
  p <- tr_post$edge[i, 1]
  c <- tr_post$edge[i, 2]
  desc_count[p] <- desc_count[p] + desc_count[c]
}

bs_num <- suppressWarnings(as.numeric(tree$node.label))
node_df <- data.frame(
  node  = (n_tips + 1):(n_tips + n_nodes),
  bs    = bs_num,
  ndesc = desc_count[(n_tips + 1):(n_tips + n_nodes)]
) %>%
  mutate(
    bs_label = ifelse(!is.na(bs) & bs >= 70 & ndesc >= 6, as.character(round(bs)), NA_character_)
  )

p_tree_base <- ggtree(tree, layout = "fan", open.angle = 15, branch.length = "none",
                      size = 0.24, color = "#8C8C8C") %<+% meta_all %<+% node_df

p_tree_bs <- p_tree_base +
  geom_nodepoint(
    aes(subset = !isTip & !is.na(bs) & bs >= 95),
    color = "#002F57", fill = "#002F57", size = 0.9, shape = 21
  ) +
  geom_text2(
    aes(subset = !isTip & !is.na(bs_label), label = bs_label),
    size = 2.3, family = FONT_FAM, fontface = "bold", color = "#002F57",
    hjust = 1.15, vjust = -0.35
  )

# 4 Slim, evenly spaced rings (Inner -> Outer: Species, Recombination, Condition, Region)
ring1 <- meta_all[, "Species", drop = FALSE]
p_A <- gheatmap(p_tree_bs, ring1, offset = 0.8, width = 0.042, colnames = FALSE, color = NA) +
  scale_fill_manual(values = PAL_SPECIES, name = "Viral Species") +
  ggnewscale::new_scale_fill()

ring2 <- meta_all[, "Recombination", drop = FALSE]
p_A <- gheatmap(p_A, ring2, offset = 3.0, width = 0.042, colnames = FALSE, color = NA) +
  scale_fill_manual(values = PAL_RECOMB_RING, name = "Recombination") +
  ggnewscale::new_scale_fill()

ring3 <- meta_all[, "Condition", drop = FALSE]
p_A <- gheatmap(p_A, ring3, offset = 5.2, width = 0.042, colnames = FALSE, color = NA) +
  scale_fill_manual(values = PAL_CONDITION, name = "Condition") +
  ggnewscale::new_scale_fill()

ring4 <- meta_all[, "Region", drop = FALSE]
p_A <- gheatmap(p_A, ring4, offset = 7.4, width = 0.042, colnames = FALSE, color = NA) +
  scale_fill_manual(values = PAL_REGION, name = "Region") +
  labs(caption = "Rings (inner to outer): Viral Species, Recombination, Condition, Region. Numbers show ultrafast bootstrap (>= 70%); solid navy circles mark >= 95%.") +
  theme(
    text             = element_text(family = FONT_FAM),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 9.5, family = FONT_FAM),
    legend.text      = element_text(size = 8.5, family = FONT_FAM),
    legend.key.size  = unit(0.35, "cm"),
    legend.spacing.y = unit(0.20, "cm"),
    plot.caption     = element_text(size = 8.8, hjust = 0.5, family = FONT_FAM, color = "#222222"),
    plot.margin      = margin(-22, 4, -12, -20)
  )

# =========================================================
# PANEL B: Recombination Hotspots Map (Preserved)
# =========================================================
gff_data <- read.table(gff_file, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(gff_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

blocks <- gff_data %>%
  filter(start > 0 & end > 0) %>%
  mutate(track_y = -0.15 - (row_number() %% 4) * 0.05)

expanded_bp <- unlist(mapply(seq, blocks$start, blocks$end))
dens <- density(expanded_bp, bw = 2000, n = 1000, from = 0, to = 152000)
density_df <- data.frame(x = dens$x, y = dens$y)
density_df$y_norm <- density_df$y / max(density_df$y)

p_B <- ggplot() +
  geom_segment(data = density_df,
               aes(x = x, xend = x, y = 0, yend = y_norm, color = y_norm),
               linewidth = 1) +
  scale_color_gradientn(colors = c("#E5E5E5", "#90A0BC", "#871746"), guide = "none") +
  geom_line(data = density_df, aes(x = x, y = y_norm), color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -0.05, color = "grey30", linewidth = 1.1) +
  geom_rect(aes(xmin = 0, xmax = 2400, ymin = -0.07, ymax = -0.03), fill = "#585D26") +
  geom_rect(aes(xmin = 148600, xmax = 151000, ymin = -0.07, ymax = -0.03), fill = "#585D26") +
  geom_segment(data = blocks,
               aes(x = start, xend = end, y = track_y, yend = track_y),
               color = "#002F57", linewidth = 1.8, lineend = "round") +
  scale_x_continuous(
    limits = c(0, 155000),
    breaks = seq(0, 150000, by = 25000),
    labels = function(x) paste0(x / 1000, " kb"),
    expand = c(0.02, 0.02)
  ) +
  scale_y_continuous(
    limits = c(-0.35, 1.05),
    breaks = c(0, 0.5, 1.0),
    labels = c("0.0", "0.5", "1.0")
  ) +
  theme_minimal(base_family = FONT_FAM) +
  labs(x = "LSDV Genome Coordinates (kb)", y = "Relative Recombination Density") +
  theme(
    axis.title.x       = element_text(face = "bold", size = 11.5, margin = margin(t = 10)),
    axis.title.y       = element_text(face = "bold", size = 11.5, margin = margin(r = 10)),
    axis.text          = element_text(color = "black", size = 10),
    panel.grid.major.x = element_line(color = "grey85", linetype = "dotted"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(8, 15, 8, 15)
  )

# =========================================================
# PANEL C: fastBAPS Clusters by Continent (Preserved)
# =========================================================
baps_plot_df <- baps_df %>%
  filter(Region %in% c("Africa", "Asia", "Europe")) %>%
  mutate(
    Cluster = factor(paste0("Cluster ", cluster), levels = paste0("Cluster ", 1:5)),
    Recombination = factor(Recombination, levels = names(PAL_RECOMB_BAR))
  ) %>%
  count(Region, Cluster, Recombination)

p_C <- ggplot(baps_plot_df, aes(x = Cluster, y = n, fill = Recombination)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.72) +
  facet_wrap(~ Region, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = PAL_RECOMB_BAR, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_bw(base_size = 11.5, base_family = FONT_FAM) +
  labs(x = NULL, y = "Number of Genomes") +
  theme(
    axis.title.y       = element_text(face = "bold", size = 11.5),
    axis.text.x        = element_text(color = "black", size = 10, angle = 30, hjust = 1),
    axis.text.y        = element_text(color = "black", size = 10),
    strip.background   = element_rect(fill = "#222222", color = "black"),
    strip.text         = element_text(color = "white", face = "bold", size = 11),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size = 10, face = "bold"),
    legend.key.size    = unit(0.4, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.8),
    plot.margin        = margin(8, 15, 6, 15)
  )

# =========================================================
# COMBINE ALL 3 PANELS INTO UNIFIED FIGURE 3
# =========================================================
fig3 <- p_A / p_B / p_C +
  plot_layout(heights = c(2.45, 0.80, 0.90)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", family = FONT_FAM))

ggsave(file.path(out_dir, "Figure_3_Final.pdf"), fig3, width = 10.5, height = 15.6, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_3_Final.png"), fig3, width = 10.5, height = 15.6, dpi = 300)
message("Saved final 4-ring Figure_3_Final.pdf and Figure_3_Final.png in: ", out_dir)
