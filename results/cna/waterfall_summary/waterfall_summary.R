library(tidyverse)
source("/home/jeszyman/repos/science/R/figure_schema.R")
library(patchwork)

d <- read_tsv("waterfall_summary.tsv", show_col_types = FALSE) %>%
  mutate(lib_id = factor(lib_id, levels = lib_id[order(sample_order)]))

treatment_fill <- c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")

# Panel A: S3A TF waterfall
pa <- ggplot(d, aes(x = lib_id, y = s3a_tf, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_point(aes(y = s1_tf), shape = 4, size = 4, stroke = 1.5, color = "#e66101") +
  geom_text(aes(label = sprintf("%.2f", s3a_tf)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = treatment_fill) +
  labs(x = NULL, y = "Tumor fraction",
       fill = "Treatment",
       title = "Tumor fraction (S3A with PoN, sorted by TF)",
       subtitle = "Bars = S3A | X marks = S1 (contaminated)") +
  theme_scifig() +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = "white", color = "gray80"))

# Panel B: Tumor volume
pb <- ggplot(d, aes(x = lib_id, y = tumor_vol_mm3, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = treatment_fill) +
  labs(x = NULL, y = expression("Tumor vol (mm"^3*")"),
       title = "Terminal tumor volume") +
  theme_scifig() +
  guides(fill = "none")

# Panel C: Human read percentage
pc <- ggplot(d, aes(x = lib_id, y = human_pct, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = treatment_fill) +
  labs(x = NULL, y = "Human reads (%)",
       title = "Human read fraction (disambiguate)") +
  theme_scifig() +
  guides(fill = "none")

# Panel D: cfDNA concentration
pd <- ggplot(d, aes(x = lib_id, y = cfdna_conc_pg_ul, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = treatment_fill) +
  labs(x = NULL, y = "cfDNA (pg/\u00B5L)",
       title = "cfDNA concentration") +
  theme_scifig() +
  guides(fill = "none")

p <- pa / (pb | pc | pd) +
  plot_layout(heights = c(0.5, 0.5)) +
  plot_annotation(
    title = "Integrated Sample Summary",
    subtitle = "Samples sorted by S3A tumor fraction | WU-487 terminal bleeds (n = 6)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

legend_text <- "WU-487 PDX terminal bleeds (n = 6). Tumor fraction estimated by ichorCNA (default-low-input preset, with panel of normals). S3A = disambiguate (human-only reads); S1 = direct hg38 alignment (mouse contamination present). Samples sorted by descending S3A tumor fraction."

save_plot(p, "waterfall_summary", w = 14, h = 10)
