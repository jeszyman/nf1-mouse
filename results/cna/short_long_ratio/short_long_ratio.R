library(tidyverse)
source("/home/jeszyman/repos/science/R/figure_schema.R")
library(patchwork)

d <- read_tsv("short_long_ratio.tsv", show_col_types = FALSE) %>%
  mutate(strategy_label = case_when(
    strategy == "s1"  ~ "S1: Direct",
    strategy == "s3a" ~ "S3A: Human only",
    strategy == "s3b" ~ "S3B: Human+Ambig"
  ))

strategy_fill <- c("S1: Direct" = "#e66101",
                    "S3A: Human only" = "#5e3c99",
                    "S3B: Human+Ambig" = "#b2abd2")

# Panel A: Short/long ratio by sample and strategy
pa <- ggplot(d, aes(x = lib_id, y = sl_ratio, fill = strategy_label)) +
  geom_col(position = position_dodge(0.8), alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = strategy_fill) +
  labs(x = NULL, y = "Short/long ratio (<150 / >250 bp)", fill = "Strategy",
       title = "Fragment short-to-long ratio") +
  theme_scifig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Short/long ratio vs TF (S3A only)
s3a <- d %>% filter(strategy == "s3a")

pb <- ggplot(s3a, aes(x = tumor_fraction, y = sl_ratio)) +
  geom_point(aes(fill = treatment), size = 5, shape = 21, stroke = 0.8) +
  geom_text(aes(label = lib_id), vjust = -1.2, size = 3.5) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "ichorCNA tumor fraction (S3A)", y = "Short/long ratio",
       fill = "Treatment",
       title = "Short/long ratio vs tumor fraction (S3A)") +
  theme_scifig()

p <- pa + pb +
  plot_annotation(
    title = "Fragment Size Distribution: Short-to-Long Ratio",
    subtitle = "Short (<150 bp) fragments enriched in tumor cfDNA (n = 6 samples)",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray30"),
                  plot.background = element_rect(fill = "white", color = NA))
  )

legend_text <- "Short/long fragment ratio computed from 1% BAM subsample (50M read subsets). Short: <150 bp; long: >250 bp. Higher ratio indicates greater tumor-derived cfDNA enrichment. S3A tumor fraction from ichorCNA (default-low-input, with PoN)."

save_plot(p, "short_long_ratio", w = 14, h = 7)
