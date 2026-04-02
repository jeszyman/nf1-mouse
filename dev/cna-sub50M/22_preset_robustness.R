## 22_preset_robustness.R — ichorCNA preset robustness analysis
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

rerun_path <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun/ichor_rerun_summary.tsv"
ichor_rerun <- read_tsv(rerun_path, show_col_types = FALSE) %>%
  mutate(lib_id = str_extract(sample, "lib\\d+"),
         strategy = case_when(
           strategy == "1" ~ "s1",
           strategy == "3a" ~ "s3a",
           strategy == "3b" ~ "s3b"
         )) %>%
  left_join(sample_meta, by = "lib_id")

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

# Focus on S3A with PoN (recommended production setting)
s3a_pon <- ichor_rerun %>%
  filter(strategy == "s3a", pon == "withpon") %>%
  mutate(preset = factor(preset,
    levels = c("aggressive", "default-low-input", "permissive"),
    labels = c("Aggressive", "Default-low-input", "Permissive")))

# Panel A: TF by preset for S3A with PoN
pa <- ggplot(s3a_pon, aes(x = lib_id, y = tumor_fraction, fill = preset)) +
  geom_col(position = position_dodge(0.8), alpha = 0.85, width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Tumor fraction", fill = "Preset",
       title = "A. S3A tumor fraction by ichorCNA preset (with PoN)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Preset agreement — TF range per sample
preset_range <- s3a_pon %>%
  group_by(lib_id) %>%
  summarise(
    tf_min = min(tumor_fraction),
    tf_max = max(tumor_fraction),
    tf_range = tf_max - tf_min,
    tf_mean = mean(tumor_fraction),
    cv = sd(tumor_fraction) / mean(tumor_fraction),
    .groups = "drop"
  ) %>%
  left_join(sample_meta %>% select(lib_id, treatment), by = "lib_id")

pb <- ggplot(preset_range, aes(x = lib_id, y = cv, fill = treatment)) +
  geom_col(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "CV across presets",
       fill = "Treatment",
       title = "B. Preset sensitivity (CV of TF across 3 presets)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel C: All strategies, all presets, withpon — dot plot
all_withpon <- ichor_rerun %>%
  filter(pon == "withpon") %>%
  mutate(
    strategy_label = case_when(
      strategy == "s1" ~ "S1",
      strategy == "s3a" ~ "S3A",
      strategy == "s3b" ~ "S3B"
    ),
    preset = factor(preset,
      levels = c("aggressive", "default-low-input", "permissive"),
      labels = c("Aggr", "Default", "Perm"))
  )

pc <- ggplot(all_withpon,
  aes(x = lib_id, y = tumor_fraction,
      color = strategy_label, shape = preset)) +
  geom_point(size = 3, alpha = 0.8, position = position_dodge(0.5)) +
  scale_color_manual(values = c(S1 = "#e66101", S3A = "#5e3c99", S3B = "#b2abd2")) +
  labs(x = NULL, y = "Tumor fraction",
       color = "Strategy", shape = "Preset",
       title = "C. All conditions: strategy x preset (with PoN)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

combined <- pa / (pb | pc) +
  plot_layout(heights = c(0.45, 0.55)) +
  plot_annotation(
    title = "ichorCNA Preset Robustness Analysis",
    subtitle = "S3A with PoN | High-TF samples (lib0626, lib0627) are robust across presets",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12, color = "gray30"),
                  plot.background = element_rect(fill = "white", color = NA))
  )

ggsave("plots/preset_robustness.png", combined, width = 14, height = 12, dpi = 300)
cat("Saved plots/preset_robustness.png\n")
