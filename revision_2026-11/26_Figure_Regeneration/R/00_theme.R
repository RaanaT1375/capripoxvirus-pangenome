# Stage 26 — shared theme and palette for the regenerated figures.
#
# Taken from the Supplementary Figure S1 map script so that the regenerated
# figures sit alongside the ones you already have: theme_minimal base, serif
# throughout, #2C3E50 text, bold legend title, a thin dark panel border, and
# the same three species colours.
#
# source() this at the top of each figure script.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---- palette ---------------------------------------------------------------
# Exactly the hex values from your map script.
SPECIES_COLS <- c(
  "LSDV" = "#002F57",   # navy
  "SPPV" = "#871746",   # maroon
  "GTPV" = "#F7B978"    # terracotta / light orange
)

# Pan-genome compartments: a sequential ramp anchored on the LSDV navy, so
# panel D reads as one family rather than as four unrelated hues.
COMPARTMENT_COLS <- c(
  "core"      = "#002F57",
  "soft-core" = "#2E6189",
  "shell"     = "#7FA6C4",
  "cloud"     = "#D6E4EF"
)

# Trait categories in the volcano plot.
CATEGORY_COLS <- c(
  "Geography"            = "#002F57",
  "Clinical"             = "#871746",
  "Host"                 = "#F7B978",
  "Species"              = "#117A65",
  "Population structure" = "#99A3A4",
  "Recombination"        = "#5D6D7E",
  "Other"                = "#BDC3C7"
)

# Analysis families in the forest plot.
FOREST_COLS <- c(
  "Polarization"  = "#002F57",
  "Recombination" = "#871746",
  "Outgroup"      = "#117A65",
  "Alignment"     = "#99A3A4",
  "Dataset"       = "#F7B978"
)

INK <- "#2C3E50"

# ---- theme ------------------------------------------------------------------
# For ordinary panels with axes.
theme_cpx <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "serif") +
    theme(
      text             = element_text(family = "serif", colour = INK),
      plot.title       = element_text(size = base_size + 3, face = "bold",
                                      hjust = 0, colour = INK,
                                      margin = margin(b = 6)),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey30",
                                      hjust = 0, margin = margin(b = 8)),
      axis.title       = element_text(size = base_size - 1, colour = INK),
      axis.text        = element_text(size = base_size - 2, colour = INK),
      panel.background = element_rect(fill = "white", colour = "black",
                                      linewidth = 0.5),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position     = "right",
      legend.background   = element_rect(fill = "white", colour = NA),
      legend.key          = element_rect(fill = "white", colour = NA),
      legend.title        = element_text(face = "bold", size = base_size - 1),
      legend.text         = element_text(size = base_size - 2),
      legend.margin       = margin(6, 6, 6, 6),
      plot.margin         = margin(8, 10, 8, 8),
      plot.tag            = element_text(size = base_size + 4, face = "bold",
                                         colour = INK)
    )
}

# For heatmaps and other panels with no axes to speak of.
theme_cpx_bare <- function(base_size = 11) {
  theme_cpx(base_size) +
    theme(
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid.major = element_blank()
    )
}

# For the donut, which needs no panel furniture at all.
theme_cpx_void <- function(base_size = 11) {
  theme_void(base_size = base_size, base_family = "serif") +
    theme(
      text          = element_text(family = "serif", colour = INK),
      legend.title  = element_text(face = "bold", size = base_size - 1),
      legend.text   = element_text(size = base_size - 2),
      legend.key    = element_rect(fill = "white", colour = NA),
      plot.tag      = element_text(size = base_size + 4, face = "bold",
                                   colour = INK),
      plot.margin   = margin(8, 10, 8, 8)
    )
}

# ---- saving -----------------------------------------------------------------
# cairo_pdf where available, matching your map script; falls back to the
# standard pdf device on a build without cairo.
save_fig <- function(plot, stem, width, height, dpi = 600) {
  dir.create(dirname(stem), showWarnings = FALSE, recursive = TRUE)
  dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else NULL
  if (!is.null(dev)) {
    ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
           device = dev)
  } else {
    message("cairo not available; using the standard pdf device")
    ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  }
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi)
  message("written: ", stem, ".pdf and .png")
}

# ---- legend placement across ggplot2 versions -------------------------------
# ggplot2 3.5 deprecated legend.position = c(x, y) in favour of
# legend.position.inside. This returns whichever the installed version wants,
# so the scripts run unchanged on either.
legend_inside <- function(x, y, just = NULL) {
  if (is.null(just)) just <- c(x, y)
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside",
          legend.position.inside = c(x, y),
          legend.justification = just)
  } else {
    theme(legend.position = c(x, y), legend.justification = just)
  }
}

# ---- paths ------------------------------------------------------------------
# Directory of the running script, so source("00_theme.R") works from anywhere.
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}

# Resolve the stage directory from $PANG, as the Python scripts do.
stage_paths <- function() {
  pang <- Sys.getenv("PANG",
                     "/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis")
  d <- file.path(pang, "revision_2026-11", "26_Figure_Regeneration")
  if (!dir.exists(d)) {
    stop("stage directory not found: ", d,
         "\nSet PANG, e.g.  export PANG=/path/to/Pangenome_Analysis")
  }
  list(stage = d,
       inputs = file.path(d, "01_inputs"),
       results = file.path(d, "02_results"),
       summary = file.path(d, "03_summary"),
       figures = file.path(d, "04_figures"))
}

need <- function(path) {
  if (!file.exists(path)) {
    stop("missing input: ", path,
         "\nRun the matching Python step in scripts/ first.")
  }
  path
}
