# Stage 26 — Figure 5B forest plot, drawn in R.
#
#   export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
#   Rscript R/figure5b.R
#
# Row values live in 01_inputs/forest_rows.csv, which ships with this stage.
# Edit that file and re-run to redraw; nothing is recomputed here. Leave the
# lo and hi cells empty for a row with no bootstrap interval.
#
# Five rows have no interval, because stage 25 saved aggregate counts rather
# than the per-gene table a bootstrap needs. They are drawn as open diamonds
# and the caption says so.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})
.sd <- local({ a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd() })
source(file.path(.sd, "00_theme.R"))

P <- stage_paths()
rows <- read.csv(need(file.path(P$inputs, "forest_rows.csv")),
                 stringsAsFactors = FALSE)
rows <- rows[!is.na(rows$alpha), ]

# order by analysis family, then keep the file's order within each family
group_order <- names(FOREST_COLS)[names(FOREST_COLS) %in% rows$group]
rows$group <- factor(rows$group, levels = group_order)
rows <- rows[order(rows$group, seq_len(nrow(rows))), ]
rows$has_ci <- !(is.na(rows$lo) | is.na(rows$hi))

# Explicit numeric positions rather than factor levels: the first row of the
# file is drawn at the top, with no dependence on how a discrete scale orders
# its levels.
rows$ypos <- rev(seq_len(nrow(rows)))

message(sprintf("plotting %d rows", nrow(rows)))
if (any(!rows$has_ci)) {
  message(sprintf("%d rows without a confidence interval:", sum(!rows$has_ci)))
  for (l in rows$label[!rows$has_ci]) message("  ", l)
}

rows$annot <- ifelse(
  is.na(rows$cmh_p),
  sprintf("%.3f", rows$alpha),
  sprintf("%.3f    p = %s", rows$alpha,
          formatC(rows$cmh_p, format = "g", digits = 2)))

xlo <- min(c(rows$lo, rows$alpha), na.rm = TRUE)
xhi <- max(c(rows$hi, rows$alpha), na.rm = TRUE)
span <- xhi - xlo
text_x <- xhi + span * 0.10          # annotation column, outside the panel

p <- ggplot(rows, aes(y = ypos, colour = group)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.5) +
  # geom_errorbar with orientation="y" rather than geom_errorbarh, which
  # ggplot2 4.0 deprecates; orientation has existed since 3.3, so this works
  # on both.
  geom_errorbar(data = dplyr::filter(rows, has_ci),
                aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.24, linewidth = 0.7) +
  geom_point(data = dplyr::filter(rows, has_ci), aes(x = alpha),
             size = 2.7, shape = 21, fill = "white", stroke = 1.1) +
  geom_point(data = dplyr::filter(rows, !has_ci), aes(x = alpha),
             size = 2.9, shape = 23, fill = "white", stroke = 1.1) +
  geom_text(aes(x = text_x, label = annot), hjust = 0,
            family = "serif", size = 2.8, colour = INK) +
  scale_colour_manual(values = FOREST_COLS, name = "Analysis",
                      drop = TRUE) +
  scale_y_continuous(breaks = rows$ypos, labels = rows$label,
                     expand = expansion(add = 0.6)) +
  coord_cartesian(xlim = c(xlo - span * 0.04, xhi + span * 0.04),
                  clip = "off") +
  labs(
    x = expression(alpha*", adaptive fraction of fixed amino-acid differences"),
    y = NULL,
    caption = paste("Open diamonds mark analyses for which no bootstrap",
                    "interval was computed.")) +
  theme_cpx() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y  = element_text(size = 8.5, hjust = 1),
        legend.position  = "top",
        legend.direction = "horizontal",
        legend.key.size  = unit(10, "pt"),
        plot.caption = element_text(size = 7.5, colour = "grey40",
                                    hjust = 0, family = "serif",
                                    margin = margin(t = 8)),
        plot.margin  = margin(8, 92, 8, 8)) +
  guides(colour = guide_legend(nrow = 1))

h <- 0.38 * nrow(rows) + 2.2
save_fig(p, file.path(P$figures, "Figure_5B"), width = 8.2, height = h)
