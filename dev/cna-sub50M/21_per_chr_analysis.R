## 21_per_chr_analysis.R — Per-chromosome read depth and contamination analysis
## Uses CNA segment data to examine chromosome-level effects
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

seg_all <- read_tsv("data/cna_segments.tsv", show_col_types = FALSE)

theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    strip.text = element_text(size = 11, face = "bold")
  )

# Compute per-chromosome median logR for each sample × strategy
chr_medians <- seg_all %>%
  filter(pon == "nopon", preset == "default-low-input") %>%
  mutate(chr_num = as.numeric(str_replace(chr, "chr", ""))) %>%
  filter(!is.na(chr_num)) %>%
  # Weight by segment length
  mutate(seg_len = end - start) %>%
  group_by(lib_id, strategy, chr_num) %>%
  summarise(
    median_logR = weighted.mean(logR, seg_len),
    .groups = "drop"
  ) %>%
  left_join(sample_meta %>% select(lib_id, treatment), by = "lib_id")

# Panel A: Heatmap of per-chr logR for lib0626 (S1 vs S3A)
lib0626_chr <- chr_medians %>%
  filter(lib_id == "lib0626", strategy %in% c("s1", "s3a")) %>%
  mutate(strategy_label = ifelse(strategy == "s1", "S1: Direct", "S3A: Disambiguated"))

p_heatmap <- ggplot(lib0626_chr,
  aes(x = factor(chr_num), y = strategy_label, fill = median_logR)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", median_logR)), size = 3.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                        midpoint = 0, limits = c(-0.6, 0.6)) +
  labs(x = "Chromosome", y = NULL, fill = "Median\nlog2R",
       title = "A. Per-chromosome CNA: lib0626 (S1 vs S3A)",
       subtitle = "S1 shows near-zero values due to mouse contamination dampening") +
  theme_pub +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "right")

# Panel B: S1 vs S3A median logR scatter per chromosome (lib0626)
chr_compare <- chr_medians %>%
  filter(lib_id == "lib0626") %>%
  select(chr_num, strategy, median_logR) %>%
  pivot_wider(names_from = strategy, values_from = median_logR)

p_scatter <- ggplot(chr_compare, aes(x = s1, y = s3a)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, color = "gray80") +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_point(aes(color = factor(chr_num)), size = 4) +
  geom_text(aes(label = chr_num), vjust = -1, size = 3) +
  scale_color_viridis_d(guide = "none") +
  labs(x = "S1 median log2R", y = "S3A median log2R",
       title = "B. Chromosome-level log2R: S1 vs S3A (lib0626)",
       subtitle = "S3A amplifies true CNA signal compressed by S1 contamination") +
  theme_pub

# Panel C: Signal compression ratio (S3A/S1 absolute logR) across all high-TF samples
compression <- chr_medians %>%
  filter(lib_id %in% c("lib0626", "lib0627")) %>%
  select(lib_id, chr_num, strategy, median_logR) %>%
  pivot_wider(names_from = strategy, values_from = median_logR) %>%
  filter(abs(s3a) > 0.05) %>%  # only chr with real signal
  mutate(compression_ratio = s3a / s1)

p_compress <- ggplot(compression,
  aes(x = factor(chr_num), y = compression_ratio, fill = lib_id)) +
  geom_col(position = position_dodge(0.8), alpha = 0.8, width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c(lib0626 = "#5e3c99", lib0627 = "#b2abd2")) +
  labs(x = "Chromosome", y = "S3A / S1 log2R ratio",
       fill = "Sample",
       title = "C. CNA signal amplification by disambiguation",
       subtitle = "Ratio > 1 = S3A detects stronger signal; below dashed line = S1 overestimates") +
  theme_pub +
  theme(axis.text.x = element_text(size = 10))

combined <- p_heatmap / (p_scatter | p_compress) +
  plot_layout(heights = c(0.4, 0.6)) +
  plot_annotation(
    title = "Per-Chromosome CNA Analysis: Contamination Impact",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.background = element_rect(fill = "white", color = NA))
  )

ggsave("plots/per_chr_analysis.png", combined, width = 16, height = 12, dpi = 300)
cat("Saved plots/per_chr_analysis.png\n")
