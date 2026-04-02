## 08_cross_correlations.R — Cross-variable correlations and integrative analysis
if (!exists("sample_meta")) source("00_load_data.R")

rerun_base <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun"
rerun_summary_path <- file.path(rerun_base, "ichor_rerun_summary.tsv")

ichor_rerun <- read_tsv(rerun_summary_path, show_col_types = FALSE) %>%
  mutate(lib_id = str_extract(sample, "lib\\d+")) %>%
  left_join(sample_meta, by = "lib_id")

# Focus on default-low-input, nopon as primary comparison
ichor_primary <- ichor_rerun %>%
  filter(preset == "default-low-input", pon == "nopon")

# --- TF vs cfDNA concentration (all strategies) ---
p1 <- ggplot(ichor_primary, aes(x = cfdna_conc_pg_ul, y = tumor_fraction,
                                 color = strategy, label = lib_id)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = "cfDNA concentration (pg/µL)", y = "Tumor fraction",
       title = "Tumor fraction vs cfDNA input") +
  theme_white
ggsave("plots/cross_tf_vs_cfdna.png", p1, width = 9, height = 5, dpi = 150)

# --- TF vs tumor volume ---
p2 <- ggplot(ichor_primary, aes(x = tumor_vol_mm3, y = tumor_fraction,
                                 color = strategy, label = lib_id)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = "Tumor volume (mm³)", y = "Tumor fraction",
       title = "Tumor fraction vs tumor volume") +
  theme_white
ggsave("plots/cross_tf_vs_tumor.png", p2, width = 9, height = 5, dpi = 150)

# --- Strategy agreement: S1 vs S3A tumor fraction ---
tf_wide <- ichor_primary %>%
  select(lib_id, strategy, tumor_fraction) %>%
  pivot_wider(names_from = strategy, values_from = tumor_fraction)

if ("s1" %in% names(tf_wide) && "s3a" %in% names(tf_wide)) {
  p3 <- ggplot(tf_wide, aes(x = s1, y = s3a, label = lib_id)) +
    geom_point(size = 4, color = "#5e3c99") +
    geom_text(vjust = -1, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = "TF: Direct (S1)", y = "TF: Disambiguate human-only (S3A)",
         title = "Strategy agreement: S1 vs S3A tumor fraction") +
    coord_equal() +
    theme_white
  ggsave("plots/cross_s1_vs_s3a.png", p3, width = 7, height = 7, dpi = 150)

  # Correlation
  cor_test <- cor.test(tf_wide$s1, tf_wide$s3a, method = "spearman")
  cat(sprintf("S1 vs S3A Spearman: rho=%.3f, p=%.4f\n", cor_test$estimate, cor_test$p.value))
}

# --- S3A vs S3B (effect of including ambiguous reads) ---
if ("s3a" %in% names(tf_wide) && "s3b" %in% names(tf_wide)) {
  p4 <- ggplot(tf_wide, aes(x = s3a, y = s3b, label = lib_id)) +
    geom_point(size = 4, color = "#b2abd2") +
    geom_text(vjust = -1, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = "TF: S3A (human only)", y = "TF: S3B (human + ambiguous)",
         title = "Effect of including ambiguous reads on TF") +
    coord_equal() +
    theme_white
  ggsave("plots/cross_s3a_vs_s3b.png", p4, width = 7, height = 7, dpi = 150)
}

# --- Composite summary table ---
composite <- ichor_primary %>%
  select(lib_id, strategy, tumor_fraction) %>%
  pivot_wider(names_from = strategy, values_from = tumor_fraction,
              names_prefix = "tf_") %>%
  left_join(sample_meta, by = "lib_id")

write_tsv(composite, "data/composite_summary.tsv")
cat("\nComposite summary:\n")
print(composite %>% select(lib_id, treatment, cfdna_conc_pg_ul, tumor_vol_mm3,
                            starts_with("tf_")))

cat("\nScript 08 complete.\n")
