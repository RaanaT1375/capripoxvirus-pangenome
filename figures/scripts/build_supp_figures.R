#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Figures S2, S4, and S5 - Unified Serif & Q1 Palette
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(patchwork)
})

FONT_FAM <- "serif"
out_dir  <- "revision_2026-11/final_supplementary/figures"

base_theme <- function() {
  theme_bw(base_size = 12, base_family = FONT_FAM) +
    theme(
      axis.title       = element_text(face = "bold", size = 11.5),
      axis.text        = element_text(color = "black", size = 10),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", linewidth = 0.8),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 10, face = "bold")
    )
}

# =====================================================================
# SUPPLEMENTARY FIGURE S2: Pan-genomic QC Thresholds (Panels A, B, C)
# =====================================================================
qc_contam <- read.delim("01_QC_and_Filtering/full_contamination_qc.tsv", stringsAsFactors = FALSE)
qc_contam$Genome <- sub("\\.faa$", "", qc_contam$file)

qc_incomp <- read.delim("01_QC_and_Filtering/incomplete_genomes_list.tsv", stringsAsFactors = FALSE)
colnames(qc_incomp)[1] <- "Genome"

qc_df <- qc_contam %>%
  left_join(qc_incomp, by = "Genome") %>%
  mutate(
    n_missing_near_core = ifelse(is.na(n_missing_near_core), 0, n_missing_near_core),
    Status = case_when(
      pct_no_hit > 20          ~ "Excluded: Contamination (>20%)",
      n_missing_near_core >= 9 ~ "Excluded: Fragmentation (>=9 missing)",
      TRUE                     ~ "Retained (n = 290)"
    )
  )

PAL_QC <- c(
  "Retained (n = 290)"                    = "#002F57",
  "Excluded: Fragmentation (>=9 missing)" = "#C29B48",
  "Excluded: Contamination (>20%)"        = "#871746"
)

p_s2a <- ggplot(qc_df, aes(x = pct_no_hit)) +
  geom_histogram(bins = 35, fill = "#90A0BC", color = "black", linewidth = 0.35) +
  geom_vline(xintercept = 20, color = "#871746", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 22, y = Inf, vjust = 1.6, hjust = 0, label = "20% cutoff",
           family = FONT_FAM, fontface = "bold", size = 3.5, color = "#871746") +
  labs(x = "Unmapped proteome (%)", y = "Number of genomes") +
  base_theme()

p_s2b <- ggplot(qc_df, aes(x = n_missing_near_core)) +
  geom_histogram(bins = 30, fill = "#90A0BC", color = "black", linewidth = 0.35) +
  geom_vline(xintercept = 9, color = "#871746", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 10.5, y = Inf, vjust = 1.6, hjust = 0, label = "9 missing genes",
           family = FONT_FAM, fontface = "bold", size = 3.5, color = "#871746") +
  labs(x = "Missing near-core genes (n)", y = "Number of genomes") +
  base_theme()

p_s2c <- ggplot(qc_df, aes(x = pct_no_hit, y = n_missing_near_core, fill = Status)) +
  geom_vline(xintercept = 20, color = "#871746", linetype = "dashed", linewidth = 0.65) +
  geom_hline(yintercept = 9,  color = "#871746", linetype = "dashed", linewidth = 0.65) +
  geom_point(shape = 21, size = 2.8, color = "#222222", stroke = 0.35, alpha = 0.85) +
  scale_fill_manual(values = PAL_QC) +
  labs(x = "Unmapped proteome (%)", y = "Missing near-core genes (n)") +
  base_theme()

fig_s2 <- (p_s2a | p_s2b) / p_s2c +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold", family = FONT_FAM))

ggsave(file.path(out_dir, "Supplementary_Figure_S2.pdf"), fig_s2, width = 10.5, height = 8.2, device = cairo_pdf)
ggsave(file.path(out_dir, "Supplementary_Figure_S2.png"), fig_s2, width = 10.5, height = 8.2, dpi = 300)
message("[OK] Saved Supplementary_Figure_S2")

# =====================================================================
# SUPPLEMENTARY FIGURE S4: Root-to-Tip Regression (205 LSDV Isolates)
# =====================================================================
rtt <- read.csv("09_Temporal_Signal/root_to_tip_LSDV_only.csv", stringsAsFactors = FALSE)
if (!"Continent" %in% names(rtt) && file.exists("07_Phylogeography/geo_metadata.csv")) {
  geo <- read.csv("07_Phylogeography/geo_metadata.csv", stringsAsFactors = FALSE)
  colnames(geo)[1] <- colnames(rtt)[1]
  rtt <- left_join(rtt, geo[, 1:2], by = colnames(rtt)[1])
  colnames(rtt)[ncol(rtt)] <- "Continent"
}

col_year <- intersect(c("year", "Year", "collection_year"), names(rtt))[1]
col_rtt  <- intersect(c("root_to_tip", "distance", "rtt"), names(rtt))[1]
col_cont <- intersect(c("Continent", "continent", "Region"), names(rtt))[1]

rtt <- rtt %>%
  rename(Year = all_of(col_year), RTT = all_of(col_rtt), Continent = all_of(col_cont)) %>%
  filter(!is.na(Year), !is.na(RTT))

PAL_CONT <- c("Africa" = "#495C7C", "Asia" = "#871746", "Europe" = "#CB9A9F", "Unknown" = "#B8B8B8")

p_s4 <- ggplot(rtt, aes(x = Year, y = RTT)) +
  geom_smooth(method = "lm", se = TRUE, color = "#002F57", fill = "#90A0BC", alpha = 0.25, linewidth = 0.9) +
  geom_point(aes(fill = Continent), shape = 21, size = 2.8, color = "#222222", stroke = 0.35, alpha = 0.88) +
  scale_fill_manual(values = PAL_CONT) +
  annotate("text", x = min(rtt$Year) + 1, y = max(rtt$RTT) * 0.98, hjust = 0, vjust = 1,
           label = "All LSDV (n = 205): slope = -4.87e-5, R^2 = 0.111, p < 0.001\nPost-2000 (n = 185): slope = -5.40e-6, R^2 < 0.001, p = 0.887",
           family = FONT_FAM, size = 3.4, color = "#111111") +
  labs(x = "Collection year", y = "Root-to-tip genetic distance (substitutions/site)") +
  base_theme()

ggsave(file.path(out_dir, "Supplementary_Figure_S4.pdf"), p_s4, width = 8.5, height = 5.6, device = cairo_pdf)
ggsave(file.path(out_dir, "Supplementary_Figure_S4.png"), p_s4, width = 8.5, height = 5.6, dpi = 300)
message("[OK] Saved Supplementary_Figure_S4")

# =====================================================================
# SUPPLEMENTARY FIGURE S5: Core Genome Selection Landscape (60 Genes)
# =====================================================================
busted <- read.csv("10_Selection_Pressure/04_summary_tables/busted_final_table.csv", stringsAsFactors = FALSE) %>%
  mutate(
    neglog10q = -log10(pmax(q_value_BH, 1e-10)),
    Trajectory = case_when(
      q_value_BH < 0.05   ~ "Significant episodic signal (q < 0.05)",
      omega_MG94 >= 0.50  ~ "Relaxed constraint (omega >= 0.50)",
      TRUE                ~ "Purifying selection"
    ),
    Label = ifelse(Orthogroup %in% c("OG0000120", "OG0000067", "OG0000068", "OG0000058",
                                     "OG0000100", "OG0000116", "OG0000123", "OG0000110"),
                   Orthogroup, NA_character_)
  )

PAL_SEL <- c(
  "Purifying selection"                    = "#002F57",
  "Significant episodic signal (q < 0.05)" = "#871746",
  "Relaxed constraint (omega >= 0.50)"     = "#585D26"
)

p_s5 <- ggplot(busted, aes(x = omega_MG94, y = neglog10q, fill = Trajectory)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#871746", linewidth = 0.65) +
  geom_point(shape = 21, size = 3.2, color = "#222222", stroke = 0.4, alpha = 0.88) +
  geom_text_repel(aes(label = Label), family = FONT_FAM, size = 3.2, fontface = "bold",
                  color = "#111111", box.padding = 0.45, min.segment.length = 0,
                  segment.color = "grey50", max.overlaps = 20, seed = 42) +
  scale_fill_manual(values = PAL_SEL) +
  annotate("text", x = max(busted$omega_MG94) * 0.98, y = -log10(0.05) + 0.18, hjust = 1,
           label = "BH q = 0.05", family = FONT_FAM, size = 3.2, color = "#871746") +
  labs(
    x = expression("Global " * omega[MG94] * " (dN/dS)"),
    y = expression("BUSTED significance (-" * log[10] * " BH " * italic(q) * "-value)")
  ) +
  base_theme()

ggsave(file.path(out_dir, "Supplementary_Figure_S5.pdf"), p_s5, width = 8.6, height = 5.8, device = cairo_pdf)
ggsave(file.path(out_dir, "Supplementary_Figure_S5.png"), p_s5, width = 8.6, height = 5.8, dpi = 300)
message("[OK] Saved Supplementary_Figure_S5")
