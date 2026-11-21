#!/usr/bin/env Rscript
# =====================================================================
# Combined Figure 5 (Panels A, B, C) - Unified Serif & Q1 Palette
# Panel A: McDonald-Kreitman Test (from alpha_confidence_intervals.csv)
# Panel B: Forest Plot of alpha Sensitivity (from forest_rows.csv, no clipping)
# Panel C: MAF Cutoff Sensitivity Curve (External Legend)
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
})

FONT_FAM <- "serif"
out_dir  <- "revision_2026-11/final_figures"

ci_file     <- "14_MK_Test/03_summary/alpha_confidence_intervals.csv"
forest_file <- "revision_2026-11/26_Figure_Regeneration/01_inputs/forest_rows.csv"

# =========================================================
# PANEL A: McDonald-Kreitman Test Bar Chart
# =========================================================
ci_df <- read.csv(ci_file, stringsAsFactors = FALSE)
# Primary polarized LSDV branch row
prim <- ci_df[ci_df[[1]] == "mk_polarized_LSDV_branch_allgenes", , drop = FALSE]

mk_df <- data.frame(
  Class = factor(
    c("Intra-species\nPolymorphism (P)", "Intra-species\nPolymorphism (P)",
      "Inter-species\nFixed Divergence (D)", "Inter-species\nFixed Divergence (D)"),
    levels = c("Intra-species\nPolymorphism (P)", "Inter-species\nFixed Divergence (D)")
  ),
  Type  = factor(c("Synonymous", "Non-synonymous", "Synonymous", "Non-synonymous"),
                 levels = c("Synonymous", "Non-synonymous")),
  Count = c(prim$Ps[1], prim$Pn[1], prim$Ds[1], prim$Dn[1])
)

PAL_MK <- c("Synonymous" = "#90A0BC", "Non-synonymous" = "#871746")

p_A <- ggplot(mk_df, aes(x = Class, y = Count, fill = Type)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.58,
           color = "black", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", Count)),
            position = position_dodge(width = 0.68),
            vjust = -0.5, family = FONT_FAM, fontface = "bold", size = 3.5, color = "black") +
  annotate("text", x = 1.5, y = 228,
           label = "alpha %~~% 60*'%'~Adaptive", parse = TRUE,
           family = FONT_FAM, fontface = "bold", size = 4.0, color = "#871746") +
  annotate("segment", x = 1.32, xend = 1.68, y = 208, yend = 208,
           colour = "#333333", linewidth = 0.7,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  scale_fill_manual(values = PAL_MK, name = NULL) +
  scale_y_continuous(limits = c(0, 345), breaks = seq(0, 300, by = 50), expand = c(0, 0)) +
  theme_classic(base_size = 12, base_family = FONT_FAM) +
  labs(x = NULL, y = "Number of Sites (Counts)") +
  theme(
    axis.title.y     = element_text(face = "bold", size = 11.5, margin = margin(r = 8)),
    axis.text.x      = element_text(color = "black", size = 10.5, face = "bold"),
    axis.text.y      = element_text(color = "black", size = 10),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 10.5, face = "bold"),
    legend.key.size  = unit(0.4, "cm"),
    axis.line        = element_line(color = "black", linewidth = 0.7),
    plot.margin      = margin(10, 12, 6, 10)
  )

# =========================================================
# PANEL B: Forest Plot of alpha Sensitivity (No Right Clipping)
# =========================================================
rows <- read.csv(forest_file, stringsAsFactors = FALSE)
rows <- rows[!is.na(rows$alpha), ]

PAL_FOREST <- c(
  "Polarization"  = "#002F57",
  "Recombination" = "#871746",
  "Outgroup"      = "#3A6B6B",
  "Alignment"     = "#7A848D",
  "Dataset"       = "#C29B48"
)

group_order <- names(PAL_FOREST)[names(PAL_FOREST) %in% rows$group]
rows$group  <- factor(rows$group, levels = group_order)
rows        <- rows[order(rows$group, seq_len(nrow(rows))), ]
rows$has_ci <- !(is.na(rows$lo) | is.na(rows$hi))
rows$ypos   <- rev(seq_len(nrow(rows)))

rows$annot <- ifelse(
  is.na(rows$cmh_p),
  sprintf("%.3f", rows$alpha),
  sprintf("%.3f   p = %s", rows$alpha, formatC(rows$cmh_p, format = "g", digits = 2))
)

xlo    <- min(c(rows$lo, rows$alpha), na.rm = TRUE)
xhi    <- max(c(rows$hi, rows$alpha), na.rm = TRUE)
span   <- xhi - xlo
text_x <- xhi + span * 0.14

p_B <- ggplot(rows, aes(y = ypos, colour = group)) +
  geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(data = filter(rows, has_ci),
                aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.24, linewidth = 0.75) +
  geom_point(data = filter(rows, has_ci), aes(x = alpha),
             size = 2.8, shape = 21, fill = "white", stroke = 1.15) +
  geom_point(data = filter(rows, !has_ci), aes(x = alpha),
             size = 3.0, shape = 23, fill = "white", stroke = 1.15) +
  geom_text(aes(x = text_x, label = annot), hjust = 0,
            family = FONT_FAM, size = 3.0, colour = "#111111") +
  scale_colour_manual(values = PAL_FOREST, name = "Analysis", drop = TRUE) +
  scale_y_continuous(breaks = rows$ypos, labels = rows$label,
                     expand = expansion(add = 0.6)) +
  coord_cartesian(xlim = c(xlo - span * 0.05, xhi + span * 0.04), clip = "off") +
  theme_bw(base_size = 11.5, base_family = FONT_FAM) +
  labs(
    x = expression(alpha * ", adaptive fraction of fixed amino-acid differences"),
    y = NULL,
    caption = "Open diamonds mark analyses for which no bootstrap confidence interval was computed."
  ) +
  theme(
    axis.title.x       = element_text(face = "bold", size = 11, margin = margin(t = 8)),
    axis.text.x        = element_text(color = "black", size = 10),
    axis.text.y        = element_text(color = "black", size = 9.5, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.border       = element_rect(color = "black", linewidth = 0.8),
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.title       = element_text(face = "bold", size = 9.5),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.35, "cm"),
    plot.caption       = element_text(size = 8.5, colour = "#333333", hjust = 0,
                                      family = FONT_FAM, margin = margin(t = 6)),
    plot.margin        = margin(10, 135, 8, 8)
  ) +
  guides(colour = guide_legend(nrow = 1))

# =========================================================
# PANEL C: MAF Cutoff Sensitivity Curve (External Legend)
# =========================================================
maf_df <- data.frame(
  cutoff = rep(c(0, 1, 2, 5, 10, 15), 2),
  Group  = rep(c("Vaccine strains excluded", "Vaccine strains included"), each = 6),
  alpha  = c(0.272, 0.448, 0.471, 0.604, 0.608, 0.556,
             0.272, 0.448, 0.471, 0.495, 0.608, 0.556),
  lo     = c(-0.31, 0.00, 0.04, 0.33, 0.33, 0.25,
             -0.12, 0.00, 0.04, 0.07, 0.33, 0.25),
  hi     = c(0.55, 0.66, 0.68, 0.74, 0.75, 0.72,
             0.55, 0.66, 0.68, 0.70, 0.75, 0.72)
)

PAL_MAF <- c(
  "Vaccine strains excluded" = "#871746",
  "Vaccine strains included" = "#90A0BC"
)

p_C <- ggplot(maf_df, aes(x = cutoff, y = alpha, color = Group, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.6) +
  geom_vline(xintercept = 5, linetype = "dotted", color = "grey55", linewidth = 0.6) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(shape = Group), size = 2.8, color = "black", stroke = 0.6) +
  scale_color_manual(values = PAL_MAF, name = NULL) +
  scale_fill_manual(values = PAL_MAF, name = NULL) +
  scale_shape_manual(values = c("Vaccine strains excluded" = 21, "Vaccine strains included" = 22), name = NULL) +
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 15),
                     labels = paste0(c(0, 1, 2, 5, 10, 15), "%")) +
  scale_y_continuous(limits = c(-0.35, 0.82), breaks = seq(-0.2, 0.8, by = 0.2)) +
  theme_classic(base_size = 12, base_family = FONT_FAM) +
  labs(
    x = "Minor-allele frequency cutoff",
    y = expression(alpha * " (adaptive fixed differences)")
  ) +
  theme(
    axis.title       = element_text(face = "bold", size = 11.5),
    axis.text        = element_text(color = "black", size = 10),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 9.5, face = "bold"),
    legend.key.size  = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    axis.line        = element_line(color = "black", linewidth = 0.7),
    plot.margin      = margin(10, 12, 6, 10)
  )

# =========================================================
# COMBINE ALL 3 PANELS INTO UNIFIED FIGURE 5
# =========================================================
design_5 <- "
  AB
  CB
"

fig5 <- p_A + p_B + p_C +
  plot_layout(design = design_5, widths = c(0.95, 1.35), heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", family = FONT_FAM))

ggsave(file.path(out_dir, "Figure_5_Final.pdf"), fig5, width = 14.2, height = 8.8, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_5_Final.png"), fig5, width = 14.2, height = 8.8, dpi = 300)
message("Saved unified Figure_5_Final.pdf and Figure_5_Final.png in: ", out_dir)
