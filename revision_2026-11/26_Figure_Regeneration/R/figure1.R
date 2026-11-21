# Stage 26 — Figure 1, drawn in R.
#
# Reads the CSVs written by scripts/03_figure1_data.py and
# scripts/02_refit_heaps.py. Run those first.
#
#   export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
#   Rscript R/figure1.R
#
# Panels A, B and D change because the 28 non-viral orthogroups are gone
# (cloud 50 -> 22). C and E are redrawn too so the whole figure comes from
# one matrix.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})
.sd <- local({ a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd() })
source(file.path(.sd, "00_theme.R"))

P <- stage_paths()

panelA_long <- read.csv(need(file.path(P$results, "fig1_panelA_matrix_long.csv")))
panelA_sp   <- read.csv(need(file.path(P$results, "fig1_panelA_species.csv")))
curves      <- read.csv(need(file.path(P$results, "accumulation_curves.csv")))
spectrum    <- read.csv(need(file.path(P$results, "fig1_panelC_spectrum.csv")))
partition   <- read.csv(need(file.path(P$results, "fig1_panelD_partition.csv")))
pca         <- read.csv(need(file.path(P$results, "fig1_panelE_pca.csv")))
heaps       <- read.csv(need(file.path(P$summary, "heaps_fit.csv")))

ll <- heaps[heaps$fit == "log-log OLS", ][1, ]
k  <- ll$k; gam <- ll$gamma
n_og <- sum(partition$orthogroups)
message(sprintf("Heaps: k = %.3f, gamma = %.4f", k, gam))
message(sprintf("pan-genome: %d orthogroups", n_og))

sp_levels <- c("LSDV", "SPPV", "GTPV")
panelA_sp$species <- factor(panelA_sp$species, levels = sp_levels)

# ---- A: species strip + presence/absence heatmap ---------------------------
# Inline labels only where the segment is wide enough to hold one; the
# narrow species rely on the legend instead, so nothing collides.
sp_label <- panelA_sp |>
  group_by(species) |>
  summarise(mid = mean(genome_rank), n = dplyr::n(), .groups = "drop") |>
  mutate(frac = n / nrow(panelA_sp)) |>
  filter(frac >= 0.08)

strip <- ggplot(panelA_sp, aes(x = genome_rank, y = 1, fill = species)) +
  geom_tile() +
  geom_text(data = sp_label, inherit.aes = FALSE,
            aes(x = mid, y = 1, label = sprintf("%s (%d)", species, n)),
            family = "serif", size = 2.8, colour = "white",
            fontface = "bold") +
  scale_fill_manual(values = SPECIES_COLS, name = "Species",
                    labels = sprintf(
                      "%s (%d)", levels(panelA_sp$species),
                      as.integer(table(panelA_sp$species)[
                        levels(panelA_sp$species)]))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_cpx_void() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 2, 0),
        legend.key.size = unit(8, "pt"),
        panel.background = element_blank(),
        plot.margin = margin(0, 10, 0, 8))

heat <- ggplot(panelA_long, aes(x = genome_rank, y = og_rank)) +
  geom_raster(fill = "#002F57") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = "Genomes", y = "Orthogroups") +
  theme_cpx_bare() +
  theme(panel.background = element_rect(fill = "#F2F4F6", colour = "black",
                                        linewidth = 0.5))

pA <- strip / heat + plot_layout(heights = c(0.045, 1))

# ---- B: accumulation curves + Heaps ----------------------------------------
curves$heaps <- k * curves$N ^ gam
heaps_lab <- sprintf("Heaps: italic(n)==%.1f*italic(N)^%.4f", k, gam)

pB <- ggplot(curves, aes(x = N)) +
  geom_ribbon(aes(ymin = pan_mean - pan_sd, ymax = pan_mean + pan_sd),
              fill = "#002F57", alpha = 0.16) +
  geom_line(aes(y = pan_mean, colour = "Pan-genome"), linewidth = 0.8) +
  geom_line(aes(y = heaps, colour = "Heaps' law fit"),
            linewidth = 0.55, linetype = "dashed") +
  geom_ribbon(aes(ymin = core_mean - core_sd, ymax = core_mean + core_sd),
              fill = "#99A3A4", alpha = 0.18) +
  geom_line(aes(y = core_mean, colour = "Strict core"), linewidth = 0.8) +
  scale_colour_manual(
    values = c("Pan-genome" = "#002F57", "Heaps' law fit" = "#871746",
               "Strict core" = "#7F8C8D"),
    breaks = c("Pan-genome", "Heaps' law fit", "Strict core"),
    name = NULL) +
  labs(x = expression("Genomes sampled ("*italic(N)*")"), y = "Orthogroups",
       subtitle = parse(text = heaps_lab)[[1]]) +
  theme_cpx() +
  legend_inside(0.97, 0.52, just = c(1, 0.5)) +
  theme(legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = NA),
        legend.key.height = unit(10, "pt"))

# ---- C: gene frequency spectrum --------------------------------------------
pC <- ggplot(spectrum, aes(x = bin_mid, y = orthogroups)) +
  geom_col(width = 9, fill = "#002F57", colour = "white", linewidth = 0.25) +
  labs(x = "Genomes in which an orthogroup is present",
       y = "Number of orthogroups") +
  theme_cpx()

# ---- D: compartment donut ---------------------------------------------------
partition$compartment <- factor(partition$compartment,
                                levels = names(COMPARTMENT_COLS))
partition <- partition[order(partition$compartment), ]
partition$ymax <- cumsum(partition$orthogroups)
partition$ymin <- partition$ymax - partition$orthogroups
partition$mid  <- (partition$ymax + partition$ymin) / 2
partition$lab  <- sprintf("%.1f%%", partition$percent)

pD <- ggplot(partition) +
  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = 2.6, xmax = 4,
                fill = compartment), colour = "white", linewidth = 0.6) +
  geom_text(aes(x = 3.3, y = mid, label = lab),
            family = "serif", size = 2.7,
            colour = ifelse(partition$compartment %in% c("core", "soft-core"),
                            "white", INK)) +
  annotate("text", x = 0, y = 0, family = "serif", fontface = "bold",
           size = 3.4, colour = INK,
           label = paste0(n_og, "\northogroups")) +
  scale_fill_manual(
    values = COMPARTMENT_COLS, name = "Compartment",
    labels = sprintf("%s (n = %d)", partition$compartment,
                     partition$orthogroups)) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  theme_cpx_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# ---- E: accessory-genome PCA ------------------------------------------------
pca$species <- factor(pca$species, levels = sp_levels)
pcv1 <- pca$PC1_var[1]; pcv2 <- pca$PC2_var[1]

pE <- ggplot(pca, aes(PC1, PC2, colour = species)) +
  geom_point(size = 1.5, alpha = 0.78) +
  scale_colour_manual(values = SPECIES_COLS, name = "Species") +
  labs(x = sprintf("PC1 (%.1f%%)", pcv1),
       y = sprintf("PC2 (%.1f%%)", pcv2)) +
  theme_cpx() +
  legend_inside(0.02, 0.98, just = c(0, 1)) +
  theme(legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = NA),
        legend.key.height = unit(10, "pt"))

# ---- assemble ---------------------------------------------------------------
fig1 <- (pA | pB) / (pC | pD | pE) +
  plot_layout(heights = c(1.25, 1), widths = c(1, 1)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "C", "D", "E"))) &
  theme(plot.tag = element_text(size = 14, face = "bold", family = "serif",
                                colour = INK))

save_fig(fig1, file.path(P$figures, "Figure_1"), width = 13.5, height = 10.5)

message("\nFigure 1 caption numbers:")
message(sprintf("  panel D: %s",
                paste(sprintf("%s n=%d", partition$compartment,
                              partition$orthogroups), collapse = ", ")))
message(sprintf("  panel E: PC1 = %.1f%%, PC2 = %.1f%%", pcv1, pcv2))
message(sprintf("  panel B: $n = %.3f \\times N^{%.4f}$", k, gam))
