library(tidyverse)
source("/home/jeszyman/repos/science/R/figure_schema.R")

d <- read_tsv("tf_vs_tumor.tsv", show_col_types = FALSE)

strategy_colors <- c(s1 = "#e66101", s3a = "#5e3c99", s3b = "#b2abd2")
strategy_labels <- c(s1 = "S1: Direct (hg38)", s3a = "S3A: Human only", s3b = "S3B: Human + ambig")
treatment_shapes <- c(Vehicle = 16, Mirdametinib = 17)

d <- d %>% mutate(treatment = str_to_title(treatment))

legend_text <- sprintf(
  "WU-487 PDX terminal bleeds (n = %d). Tumor fraction estimated by ichorCNA (default-low-input preset, with panel of normals). Three alignment strategies compared: S1 (direct hg38 alignment, retains mouse contamination), S3A (disambiguate, human-only reads), S3B (disambiguate, human + ambiguous reads). S3A yields highest tumor fractions, consistent with removal of mouse-derived reads that dilute the signal in S1. Single PDX line shown; cross-line validation pending.",
  n_distinct(d$lib_id))

p <- ggplot(d, aes(x = tumor_vol_mm3, y = tumor_fraction,
                    color = strategy, shape = treatment)) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  scale_shape_manual(values = treatment_shapes) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.2), limits = c(0, NA)) +
  labs(x = expression("Terminal tumor volume (mm"^3*")"),
       y = "ichorCNA tumor fraction",
       color = "Strategy", shape = "Treatment",
       caption = str_wrap(legend_text, width = 105)) +
  theme_scifig() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.justification = "center",
        legend.box.just = "center",
        legend.box.margin = margin(t = 10, b = 5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(4, "pt"),
        legend.spacing.y = unit(4, "pt"),
        plot.caption = element_text(size = 10, hjust = 0, color = "gray30",
                                    lineheight = 1.3,
                                    margin = margin(t = 20)),
        plot.caption.position = "plot",
        plot.margin = margin(10, 10, 10, 10, "pt")) +
  guides(color = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 1))

save_plot(p, "tf_vs_tumor", w = 7, h = 8.5)
