## nb03_fragment_investigation.R
## Fragment length distribution analysis and peak position investigation
## Key question: Why do S3A high-TF samples peak at 145-151 bp, not 167 bp?
library(tidyverse)
library(patchwork)

# --- Paths ---
dev_dir <- "/home/jeszyman/repos/nf1-mouse/.claude/worktrees/agent-a8c2930d/dev/cna-sub50M"
plot_dir <- file.path(dev_dir, "plots")

# --- Sample metadata ---
sample_meta <- tibble(
  lib_id    = paste0("lib0", 622:627),
  mouse_id  = paste0("mou000", 1:6),
  treatment = c(rep("vehicle", 3), rep("mirdametinib", 3)),
  cfdna_conc_pg_ul = c(467, 536, 111, 1720, 2380, 387),
  tumor_vol_mm3 = c(1823, 1410, 1048, 1902, 1868, 2551)
)

theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )


# === Load fragment data ===
frag <- read_tsv(file.path(dev_dir, "data/fragment_lengths.tsv"),
                 show_col_types = FALSE) %>%
  filter(insert_size >= 50, insert_size <= 400) %>%
  left_join(sample_meta, by = "lib_id")

# Classify TF groups based on S3A TF
ichor <- read_tsv(file.path(dev_dir, "data/ichor_rerun_full.tsv"),
                  show_col_types = FALSE) %>%
  filter(preset == "default-low-input", pon == "nopon") %>%
  mutate(strategy = case_when(
    strategy == "1" ~ "s1",
    strategy == "3a" ~ "s3a",
    strategy == "3b" ~ "s3b"
  ))

s3a_tf <- ichor %>%
  filter(strategy == "s3a") %>%
  select(lib_id, s3a_tf = tumor_fraction)

frag <- frag %>%
  left_join(s3a_tf, by = "lib_id") %>%
  mutate(tf_group = ifelse(s3a_tf > 0.15, "High TF", "Low TF"))

# Normalize to % per sample x strategy
frag_pct <- frag %>%
  group_by(lib_id, strategy) %>%
  mutate(pct = 100 * pairs_total / sum(pairs_total)) %>%
  ungroup()


# === Compute peak positions ===
peaks <- frag_pct %>%
  filter(insert_size >= 100, insert_size <= 250) %>%
  group_by(lib_id, strategy, tf_group) %>%
  slice_max(pct, n = 1) %>%
  ungroup() %>%
  select(lib_id, strategy, tf_group, peak_bp = insert_size, peak_pct = pct)

cat("\n=== Peak positions (100-250 bp range) ===\n")
print(peaks %>% arrange(strategy, lib_id), n = 30)

# Total read counts per sample x strategy
counts <- frag %>%
  group_by(lib_id, strategy) %>%
  summarise(total_pairs = sum(pairs_total), .groups = "drop")

cat("\n=== Total fragment counts (1% subsample) ===\n")
print(counts %>% pivot_wider(names_from = strategy, values_from = total_pairs))


# === Figure 7: Fragment length profiles S1 vs S3A, high-TF vs low-TF ===
high_tf <- frag_pct %>%
  filter(tf_group == "High TF", strategy %in% c("s1", "s3a")) %>%
  mutate(strategy_label = ifelse(strategy == "s1", "S1: Direct", "S3A: Human only"))

low_tf <- frag_pct %>%
  filter(tf_group == "Low TF", strategy %in% c("s1", "s3a")) %>%
  mutate(strategy_label = ifelse(strategy == "s1", "S1: Direct", "S3A: Human only"))

# Get peak positions for annotation
high_peaks <- peaks %>%
  filter(tf_group == "High TF", strategy %in% c("s1", "s3a"))

pa7 <- ggplot(high_tf, aes(x = insert_size, y = pct,
                             color = strategy_label)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 167, linetype = "dotted", color = "gray50") +
  annotate("text", x = 172, y = max(high_tf$pct) * 0.95,
           label = "167 bp\n(expected)", size = 3.5, color = "gray40", hjust = 0) +
  scale_color_manual(values = c("S1: Direct" = "#e66101",
                                 "S3A: Human only" = "#5e3c99")) +
  facet_wrap(~lib_id, ncol = 2, scales = "free_y") +
  labs(x = "Insert size (bp)", y = "% of fragments",
       color = NULL,
       title = "A. High-TF samples: S3A peak shifted to 145-151 bp") +
  theme_pub +
  theme(legend.position = "bottom")

pb7 <- ggplot(low_tf, aes(x = insert_size, y = pct,
                            color = strategy_label)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 167, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = c("S1: Direct" = "#e66101",
                                 "S3A: Human only" = "#5e3c99")) +
  facet_wrap(~lib_id, ncol = 2, scales = "free_y") +
  labs(x = "Insert size (bp)", y = "% of fragments",
       color = NULL,
       title = "B. Low-TF samples: S1 and S3A overlap, peak near 165-169 bp") +
  theme_pub +
  theme(legend.position = "none")

fig7 <- pa7 / pb7 +
  plot_layout(heights = c(0.45, 0.55)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(file.path(plot_dir, "nb_fig7_fragment_profiles.png"), fig7,
       width = 12, height = 14, dpi = 300)
cat("Saved nb_fig7_fragment_profiles.png\n")


# === Figure 8: Peak position summary ===
peaks_summary <- peaks %>%
  filter(strategy %in% c("s1", "s3a")) %>%
  mutate(strategy_label = ifelse(strategy == "s1", "S1: Direct", "S3A: Human only"))

fig8 <- ggplot(peaks_summary,
  aes(x = lib_id, y = peak_bp, fill = strategy_label, label = peak_bp)) +
  geom_col(position = position_dodge(0.8), alpha = 0.85, width = 0.7) +
  geom_text(position = position_dodge(0.8), vjust = -0.5, size = 4) +
  geom_hline(yintercept = 167, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.5, y = 169, label = "167 bp (mono-nucleosomal)",
           size = 3.5, color = "gray40", hjust = 0) +
  scale_fill_manual(values = c("S1: Direct" = "#e66101",
                                "S3A: Human only" = "#5e3c99")) +
  facet_wrap(~tf_group, scales = "free_x") +
  coord_cartesian(ylim = c(130, 180)) +
  labs(x = NULL, y = "Peak insert size (bp)",
       fill = NULL,
       title = "Fragment length peak positions by strategy and TF group") +
  theme_pub +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "nb_fig8_peak_positions.png"), fig8,
       width = 10, height = 6, dpi = 300)
cat("Saved nb_fig8_peak_positions.png\n")


# === Figure 9: Low-TF investigation ===
# Ambiguous fraction vs TF
disambig_raw <- read_tsv(file.path(dev_dir, "data/disambig_percentages.tsv"),
                         show_col_types = FALSE) %>%
  select(lib_id, human_pct, mouse_pct, ambiguous_pct)

low_tf_invest <- sample_meta %>%
  left_join(disambig_raw, by = "lib_id") %>%
  left_join(s3a_tf, by = "lib_id") %>%
  mutate(tf_group = ifelse(s3a_tf > 0.15, "High TF", "Low TF"))

pa9 <- ggplot(low_tf_invest,
  aes(x = s3a_tf, y = ambiguous_pct, color = treatment)) +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(aes(label = lib_id), size = 3.5) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "S3A tumor fraction", y = "Ambiguous reads (%)",
       color = "Treatment",
       title = "A. Ambiguous read fraction vs tumor fraction") +
  theme_pub

pb9 <- ggplot(low_tf_invest,
  aes(x = s3a_tf, y = human_pct, color = treatment)) +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(aes(label = lib_id), size = 3.5) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "S3A tumor fraction", y = "Human reads (%)",
       color = "Treatment",
       title = "B. Human read fraction vs tumor fraction") +
  theme_pub +
  theme(legend.position = "none")

pc9 <- ggplot(low_tf_invest,
  aes(x = cfdna_conc_pg_ul, y = s3a_tf, color = treatment)) +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(aes(label = lib_id), size = 3.5) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "cfDNA concentration (pg/uL)", y = "S3A tumor fraction",
       color = "Treatment",
       title = "C. cfDNA concentration vs tumor fraction") +
  theme_pub +
  theme(legend.position = "none")

pd9 <- ggplot(low_tf_invest,
  aes(x = tumor_vol_mm3, y = cfdna_conc_pg_ul, color = treatment)) +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(aes(label = lib_id), size = 3.5) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = expression("Tumor volume (mm"^3*")"),
       y = "cfDNA concentration (pg/uL)",
       color = "Treatment",
       title = "D. cfDNA concentration vs tumor volume") +
  theme_pub +
  theme(legend.position = "none")

fig9 <- (pa9 | pb9) / (pc9 | pd9) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(file.path(plot_dir, "nb_fig9_lowtf_investigation.png"), fig9,
       width = 14, height = 12, dpi = 300)
cat("Saved nb_fig9_lowtf_investigation.png\n")

cat("Part 3 (fragment investigation) complete.\n")
