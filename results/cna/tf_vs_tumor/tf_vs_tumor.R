library(tidyverse)
source("/home/jeszyman/repos/science/R/figure_schema.R")
library(ggrepel)

d <- read_tsv("tf_vs_tumor.tsv", show_col_types = FALSE)

strategy_colors <- c(s1 = "#e66101", s3a = "#5e3c99", s3b = "#b2abd2")
strategy_labels <- c(s1 = "Direct (hg38)", s3a = "Disambiguate (human only)", s3b = "Disambiguate (human + ambig)")
treatment_shapes <- c(vehicle = 16, mirdametinib = 17)

p <- ggplot(d, aes(x = tumor_vol_mm3, y = tumor_fraction,
                    color = strategy, shape = treatment)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = lib_id), size = 3.5, show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  scale_shape_manual(values = treatment_shapes) +
  labs(x = expression("Terminal tumor volume (mm"^3*")"),
       y = "ichorCNA tumor fraction",
       color = "Strategy", shape = "Treatment",
       title = "Tumor fraction vs tumor volume",
       subtitle = sprintf("WU-487 PDX terminal bleeds (n = %d) | ichorCNA default-low-input with PoN",
                          n_distinct(d$lib_id))) +
  theme_scifig() +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 12, color = "gray30"))

legend_text <- sprintf(
  "WU-487 PDX terminal bleeds (n = %d). Tumor fraction estimated by ichorCNA (default-low-input preset, with panel of normals). Three alignment strategies compared: S1 (direct hg38), S3A (disambiguate, human-only reads), S3B (disambiguate, human + ambiguous reads). Single PDX line shown; cross-line validation pending.",
  n_distinct(d$lib_id))

save_plot(p, "tf_vs_tumor", w = 10, h = 7)
