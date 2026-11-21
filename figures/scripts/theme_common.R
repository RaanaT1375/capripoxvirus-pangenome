# =====================================================================
# Unified Q1 Publication Theme & Soft Muted Palette (Sans-Serif)
# No in-plot legends; clean full borders; harmonized across all figures
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
})

# Soft Muted Q1 Palette (from user's original suite)
PAL_PANGENOME <- c(
  "Core (>=99%)"       = "#3D5A80",  # Steel Blue
  "Soft-core (95-99%)" = "#81B29A",  # Sage Green
  "Shell (15-95%)"     = "#F4A261",  # Sandy Orange
  "Cloud (<15%)"       = "#E07A5F"   # Dusty Rose
)

PAL_SPECIES <- c(
  "LSDV"    = "#3D5A80",  # Steel Blue
  "SPPV"    = "#E07A5F",  # Dusty Rose
  "GTPV"    = "#81B29A",  # Sage Green
  "Unknown" = "#9E9E9E"   # Neutral Grey
)

PAL_RECOMB <- c(
  "Ancestral recombinant"   = "#E07A5F",
  "Independent recombinant" = "#F4A261",
  "No major recombination"  = "#3D5A80"
)

PAL_REGION <- c(
  "Africa" = "#3D5A80",
  "Asia"   = "#E07A5F",
  "Europe" = "#81B29A"
)

# Clean Q1 theme (sans-serif, full border, external bottom legend)
theme_q1_sans <- function(base_size = 9.5, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title        = element_text(face = "bold", size = base_size + 1.5, hjust = 0, color = "#111111"),
      axis.title        = element_text(face = "bold", size = base_size, color = "#222222"),
      axis.text         = element_text(size = base_size - 1, color = "#222222"),
      panel.border      = element_rect(color = "#444444", fill = NA, linewidth = 0.65),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "#E9ECEF", linewidth = 0.4),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_size - 1, face = "bold", color = "#333333"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.key.size   = unit(0.35, "cm"),
      legend.margin     = margin(t = 2, b = 0),
      strip.background  = element_rect(fill = "#3D5A80", color = "#444444", linewidth = 0.5),
      strip.text        = element_text(color = "white", face = "bold", size = base_size - 0.5),
      plot.margin       = margin(6, 8, 6, 8)
    )
}
