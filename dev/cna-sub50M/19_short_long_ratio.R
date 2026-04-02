## 19_short_long_ratio.R — Short-to-long fragment ratio analysis
## Short fragments (<150 bp) are enriched in tumor cfDNA vs long (>250 bp)
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

frag <- read_tsv("data/fragment_lengths.tsv", show_col_types = FALSE)

# Compute short-to-long ratio per sample x strategy
sl_ratio <- frag %>%
  mutate(size_class = case_when(
    insert_size < 150 ~ "short",
    insert_size > 250 ~ "long",
    TRUE ~ "nucleosomal"
  )) %>%
  group_by(lib_id, strategy, size_class) %>%
  summarise(total_pairs = sum(pairs_total), .groups = "drop") %>%
  pivot_wider(names_from = size_class, values_from = total_pairs) %>%
  mutate(sl_ratio = short / long) %>%
  left_join(sample_meta %>% select(lib_id, treatment, tumor_vol_mm3), by = "lib_id") %>%
  mutate(strategy_label = case_when(
    strategy == "s1" ~ "S1: Direct",
    strategy == "s3a" ~ "S3A: Human only",
    strategy == "s3b" ~ "S3B: Human+Ambig"
  ))

# Get TF for S3A withpon default-low-input
rerun_path <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun/ichor_rerun_summary.tsv"
ichor_rerun <- read_tsv(rerun_path, show_col_types = FALSE) %>%
  mutate(lib_id = str_extract(sample, "lib\\d+"),
         strategy = case_when(
           strategy == "1" ~ "s1",
           strategy == "3a" ~ "s3a",
           strategy == "3b" ~ "s3b"
         )) %>%
  filter(preset == "default-low-input", pon == "withpon")

sl_ratio <- sl_ratio %>%
  left_join(ichor_rerun %>% select(lib_id, strategy, tumor_fraction), by = c("lib_id", "strategy"))

theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

# Panel A: Short/long ratio by sample and strategy
pa <- ggplot(sl_ratio, aes(x = lib_id, y = sl_ratio, fill = strategy_label)) +
  geom_col(position = position_dodge(0.8), alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c("S1: Direct" = "#e66101",
                                "S3A: Human only" = "#5e3c99",
                                "S3B: Human+Ambig" = "#b2abd2")) +
  labs(x = NULL, y = "Short/long ratio (<150 / >250 bp)", fill = "Strategy",
       title = "A. Fragment short-to-long ratio") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Short/long ratio vs TF (S3A only)
s3a_sl <- sl_ratio %>% filter(strategy == "s3a")
pb <- ggplot(s3a_sl, aes(x = tumor_fraction, y = sl_ratio)) +
  geom_point(aes(fill = treatment), size = 5, shape = 21, stroke = 0.8) +
  geom_text(aes(label = lib_id), vjust = -1.2, size = 3.5) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "ichorCNA tumor fraction (S3A)", y = "Short/long ratio",
       fill = "Treatment",
       title = "B. Short/long ratio vs tumor fraction (S3A)") +
  theme_pub

combined <- pa + pb +
  plot_annotation(
    title = "Fragment Size Distribution: Short-to-Long Ratio",
    subtitle = "Short (<150 bp) fragments enriched in tumor cfDNA | 1% subsample from 50M subsets",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray30"),
                  plot.background = element_rect(fill = "white", color = NA))
  )

ggsave("plots/short_long_ratio_pub.png", combined, width = 14, height = 7, dpi = 300)
cat("Saved plots/short_long_ratio_pub.png\n")
