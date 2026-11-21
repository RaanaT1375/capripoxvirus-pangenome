# Stage 28 — selective constraint by functional class.
#
#   export PANG=/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
#   Rscript R/constraint.R
#
# Reads 02_results/gene_constraint.csv from scripts/01_constraint_structure.py
# and uses the shared theme from stage 26, so the figure matches the others.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

.sd <- local({ a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd() })

pang <- Sys.getenv("PANG", "/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis")
theme_src <- file.path(pang, "revision_2026-11", "26_Figure_Regeneration",
                       "R", "00_theme.R")
if (!file.exists(theme_src)) {
  stop("stage 26's R/00_theme.R not found at\n  ", theme_src,
       "\nUnpack stage 26 first; this script reuses its theme.")
}
source(theme_src)

stage <- file.path(pang, "revision_2026-11", "28_Positive_Results")
inp <- file.path(stage, "02_results", "gene_constraint.csv")
if (!file.exists(inp)) {
  stop("missing input: ", inp,
       "\nRun scripts/01_constraint_structure.py first.")
}

d <- read.csv(inp, stringsAsFactors = FALSE)
d <- d[d$category != "Unassigned", ]
message(sprintf("%d genes in %d classes", nrow(d),
                dplyr::n_distinct(d$category)))

# order classes by median constraint, most constrained first
ord <- d |> group_by(category) |>
  summarise(m = median(omega), n = dplyr::n(), .groups = "drop") |>
  arrange(m)
d$category <- factor(d$category, levels = ord$category)
lab <- setNames(sprintf("%s\n(n = %d)", ord$category, ord$n), ord$category)

# a palette spanning the navy of the species scale
pal <- colorRampPalette(c("#002F57", "#4A7BA7", "#9BB7D4",
                          "#F7B978", "#871746"))(nlevels(d$category))

p <- ggplot(d, aes(category, omega, fill = category)) +
  geom_hline(yintercept = 1, colour = "grey70", linetype = "dashed",
             linewidth = 0.4) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85,
               linewidth = 0.35, colour = INK) +
  geom_jitter(width = 0.14, height = 0, size = 1.5, alpha = 0.6,
              colour = INK) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = lab) +
  labs(x = NULL,
       y = expression(omega[MG94]*" (inter-species)"),
       caption = paste("Dashed line marks neutrality. Every core gene lies",
                       "below it; the classes differ in how far below.")) +
  theme_cpx() +
  theme(axis.text.x = element_text(size = 8),
        plot.caption = element_text(size = 7.5, colour = "grey40",
                                    hjust = 0, family = "serif",
                                    margin = margin(t = 8)))

save_fig(p, file.path(stage, "04_figures", "constraint_by_class"),
         width = 7.6, height = 5.0)

message("\nmedians, most constrained first:")
for (i in seq_len(nrow(ord))) {
  message(sprintf("  %-22s %.4f  (n = %d)", ord$category[i], ord$m[i],
                  ord$n[i]))
}
