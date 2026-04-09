library(tidyverse)
source("/home/jeszyman/repos/science/R/figure_schema.R")
library(patchwork)
library(ggrepel)

d <- read_tsv("tf_vs_tumor.tsv", show_col_types = FALSE)

strategy_colors <- c(s1 = "#e66101", s3a = "#5e3c99", s3b = "#b2abd2")
strategy_labels <- c(s1 = "Direct (hg38)", s3a = "Disambiguate (human only)", s3b = "Disambiguate (human + ambig)")
treatment_colors <- c(vehicle = "#1b9e77", mirdametinib = "#d95f02")

# Panel A: TF vs tumor volume, all strategies
pa <- ggplot(d, aes(x = tumor_vol_mm3, y = tumor_fraction,
                     color = strategy, shape = treatment)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = lib_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  scale_shape_manual(values = c(vehicle = 16, mirdametinib = 17)) +
  labs(x = expression("Terminal tumor volume (mm"^3*")"),
       y = "ichorCNA tumor fraction",
       color = "Strategy", shape = "Treatment",
       title = "Tumor fraction vs tumor volume, all strategies") +
  theme_scifig() +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"))

# Panel B: S3A only with loess trend
s3a <- d %>% filter(strategy == "s3a")
cor_result <- cor.test(s3a$tumor_vol_mm3, s3a$tumor_fraction, method = "spearman")

pb <- ggplot(s3a, aes(x = tumor_vol_mm3, y = tumor_fraction)) +
  geom_point(aes(color = treatment), size = 5) +
  geom_smooth(method = "loess", se = FALSE, color = "#5e3c99", linewidth = 0.8) +
  geom_text_repel(aes(label = lib_id), size = 3.5) +
  scale_color_manual(values = treatment_colors) +
  annotate("text", x = max(s3a$tumor_vol_mm3) * 0.55, y = max(s3a$tumor_fraction) * 0.95,
           label = sprintf("Spearman \u03C1 = %.2f\n%s\nn = %d",
                           cor_result$estimate, format_p(cor_result$p.value), nrow(s3a)),
           size = 4, hjust = 0) +
  labs(x = expression("Terminal tumor volume (mm"^3*")"),
       y = "ichorCNA tumor fraction (S3A)",
       color = "Treatment",
       title = "S3A tumor fraction vs tumor volume") +
  theme_scifig() +
  theme(legend.position = c(0.15, 0.9),
        legend.background = element_rect(fill = "white", color = "gray80"))

p <- pa / pb +
  plot_annotation(
    title = "Tumor Fraction vs Tumor Volume",
    subtitle = "WU-487 PDX terminal bleeds (n = 6) | ichorCNA default-low-input with PoN",
    theme = theme(plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray30"),
                  plot.background = element_rect(fill = "white", color = NA))
  )

legend_text <- sprintf(
  "WU-487 PDX terminal bleeds (n = %d). Tumor fraction estimated by ichorCNA (default-low-input preset, with panel of normals). Panel A: all three alignment strategies shown; no statistical test applied. Panel B: Spearman rank correlation, S3A strategy only; rho = %.2f, %s; n = %d. Single PDX line shown; cross-line validation pending.",
  n_distinct(d$lib_id), cor_result$estimate, format_p(cor_result$p.value), nrow(s3a))

save_plot(p, "tf_vs_tumor", w = 10, h = 12)
