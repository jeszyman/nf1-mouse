## 23_fragment_pub.R — Publication-quality fragment length comparison
## Focuses on the key contrast: high-TF vs low-TF samples, S1 vs S3A
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

frag <- read_tsv("data/fragment_lengths.tsv", show_col_types = FALSE) %>%
  filter(insert_size >= 50, insert_size <= 400) %>%
  left_join(sample_meta %>% select(lib_id, treatment, tumor_vol_mm3, cfdna_conc_pg_ul),
            by = "lib_id")

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    strip.text = element_text(size = 11, face = "bold")
  )

# Normalize to % per sample × strategy
frag_pct <- frag %>%
  group_by(lib_id, strategy) %>%
  mutate(pct = 100 * pairs_total / sum(pairs_total)) %>%
  ungroup()

# Classify samples by TF
frag_pct <- frag_pct %>%
  mutate(
    tf_group = ifelse(lib_id %in% c("lib0626", "lib0627"), "High TF", "Low TF"),
    strategy_label = case_when(
      strategy == "s1" ~ "S1: Direct",
      strategy == "s3a" ~ "S3A: Human only",
      strategy == "s3b" ~ "S3B: Human+Ambig"
    )
  )

# Panel A: High-TF samples — S1 vs S3A (clear nucleosomal peak after disambig)
high_tf <- frag_pct %>%
  filter(tf_group == "High TF", strategy %in% c("s1", "s3a"))

pa <- ggplot(high_tf, aes(x = insert_size, y = pct,
                            color = strategy_label, linetype = strategy_label)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_vline(xintercept = 167, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  annotate("text", x = 175, y = max(high_tf$pct) * 0.9,
           label = "167 bp", size = 3.5, color = "gray40") +
  scale_color_manual(values = c("S1: Direct" = "#e66101", "S3A: Human only" = "#5e3c99")) +
  scale_linetype_manual(values = c("S1: Direct" = "solid", "S3A: Human only" = "solid")) +
  facet_wrap(~lib_id, ncol = 2) +
  labs(x = "Insert size (bp)", y = "% of total fragments",
       color = NULL, linetype = NULL,
       title = "A. High-TF samples: cfDNA nucleosomal peak emerges after disambiguation") +
  theme_pub +
  theme(legend.position = "bottom")

# Panel B: Low-TF samples — S1 vs S3A (no nucleosomal peak difference)
low_tf <- frag_pct %>%
  filter(tf_group == "Low TF", strategy %in% c("s1", "s3a"))

pb <- ggplot(low_tf, aes(x = insert_size, y = pct,
                           color = strategy_label)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_vline(xintercept = 167, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = c("S1: Direct" = "#e66101", "S3A: Human only" = "#5e3c99")) +
  facet_wrap(~lib_id, ncol = 2) +
  labs(x = "Insert size (bp)", y = "% of total fragments",
       color = NULL,
       title = "B. Low-TF samples: no nucleosomal peak regardless of strategy") +
  theme_pub +
  theme(legend.position = "none")

combined <- pa / pb +
  plot_layout(heights = c(0.45, 0.55)) +
  plot_annotation(
    title = "Fragment Length Profiles: S1 vs S3A",
    subtitle = "% normalization per strategy | Vertical line = 167 bp mono-nucleosomal peak",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray30"),
                  plot.background = element_rect(fill = "white", color = NA))
  )

ggsave("plots/fragment_pub.png", combined, width = 12, height = 12, dpi = 300)
cat("Saved plots/fragment_pub.png\n")
