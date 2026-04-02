## 05_tumor_fraction.R — ichorCNA tumor fraction comparison across strategies/presets/PoN
if (!exists("sample_meta")) source("00_load_data.R")

rerun_base <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun"
rerun_summary_path <- file.path(rerun_base, "ichor_rerun_summary.tsv")

if (!file.exists(rerun_summary_path)) {
  stop("ichorCNA re-run not complete. Run 04_rerun_ichorcna.sh first.")
}

ichor_rerun <- read_tsv(rerun_summary_path, show_col_types = FALSE) %>%
  mutate(lib_id = str_extract(sample, "lib\\d+")) %>%
  left_join(sample_meta, by = "lib_id")

write_tsv(ichor_rerun, "data/ichor_rerun_full.tsv")

# --- Original (nopon, with chrX) vs re-run (nopon, autosomes only) ---
ichor_orig <- ichor_data %>%
  mutate(pon = "nopon_chrX", lib_id = str_extract(sample, "lib\\d+"))

ichor_nopon_auto <- ichor_rerun %>%
  filter(pon == "nopon") %>%
  mutate(pon = "nopon_auto")

compare_chrx <- bind_rows(
  ichor_orig %>% select(lib_id, strategy, preset, pon, tumor_fraction, ploidy),
  ichor_nopon_auto %>% select(lib_id, strategy, preset, pon, tumor_fraction, ploidy)
)

# Effect of dropping chrX
p_chrx <- ggplot(compare_chrx %>% filter(preset == "default-low-input"),
  aes(x = lib_id, y = tumor_fraction, fill = pon)) +
  geom_col(position = "dodge", alpha = 0.8) +
  facet_wrap(~strategy, labeller = labeller(strategy = strategy_labels)) +
  labs(x = NULL, y = "Tumor fraction", fill = "ChrX handling",
       title = "Effect of dropping chrX on tumor fraction (default-low-input preset)") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/chrx_effect.png", p_chrx, width = 14, height = 5, dpi = 150)

# --- PoN effect (autosomes only) ---
ichor_pon_compare <- ichor_rerun %>%
  select(lib_id, strategy, preset, pon, tumor_fraction, ploidy)

p_pon <- ggplot(ichor_pon_compare %>% filter(preset == "default-low-input"),
  aes(x = lib_id, y = tumor_fraction, fill = pon)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c(nopon = "#fee0d2", withpon = "#de2d26")) +
  facet_wrap(~strategy, labeller = labeller(strategy = strategy_labels)) +
  labs(x = NULL, y = "Tumor fraction", fill = "PoN",
       title = "Effect of panel of normals on tumor fraction (default-low-input)") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/pon_effect.png", p_pon, width = 14, height = 5, dpi = 150)

# --- Full heatmap: sample × strategy × preset × PoN ---
p_heatmap <- ggplot(ichor_rerun,
  aes(x = interaction(strategy, pon, sep = "\n"),
      y = lib_id, fill = tumor_fraction)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", tumor_fraction)), size = 2.5) +
  scale_fill_gradient2(low = "white", mid = "#fdae61", high = "#d73027",
                        midpoint = 0.3, limits = c(0, 0.7)) +
  facet_wrap(~preset) +
  labs(x = "Strategy / PoN", y = NULL, fill = "TF",
       title = "Tumor fraction: all conditions") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
ggsave("plots/tf_heatmap.png", p_heatmap, width = 16, height = 6, dpi = 150)

# --- Cross-sample consistency (CV within treatment arms) ---
consistency <- ichor_rerun %>%
  filter(preset == "default-low-input") %>%
  group_by(strategy, pon, treatment) %>%
  summarise(
    mean_tf = mean(tumor_fraction),
    sd_tf = sd(tumor_fraction),
    cv = sd(tumor_fraction) / mean(tumor_fraction),
    n = n(),
    .groups = "drop"
  )

write_tsv(consistency, "data/tf_consistency.tsv")
cat("\nCross-sample consistency (default-low-input):\n")
print(consistency)

# --- TF vs cfDNA concentration ---
p_tf_vs_cfdna <- ggplot(ichor_rerun %>% filter(preset == "default-low-input"),
  aes(x = cfdna_conc_pg_ul, y = tumor_fraction, color = strategy, shape = pon)) +
  geom_point(size = 3) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = "cfDNA concentration (pg/µL)", y = "Tumor fraction",
       title = "Tumor fraction vs cfDNA input (default-low-input)") +
  theme_white
ggsave("plots/tf_vs_cfdna.png", p_tf_vs_cfdna, width = 9, height = 5, dpi = 150)

cat("\nScript 05 complete.\n")
