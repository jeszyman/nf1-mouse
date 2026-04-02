## 16_strategy_concordance.R — Strategy concordance matrix and bin-level agreement
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

rerun_path <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun/ichor_rerun_summary.tsv"
ichor_rerun <- read_tsv(rerun_path, show_col_types = FALSE) %>%
  mutate(lib_id = str_extract(sample, "lib\\d+"))

# Default-low-input preset, both PoN
primary <- ichor_rerun %>%
  filter(preset == "default-low-input") %>%
  select(lib_id, strategy, pon, tumor_fraction) %>%
  mutate(strategy = case_when(
    strategy == "1" ~ "s1",
    strategy == "3a" ~ "s3a",
    strategy == "3b" ~ "s3b"
  ))

# Panel A: Strategy × PoN concordance heatmap with correlation
strategies_all <- expand_grid(
  s1 = c("s1"), s2 = c("s1", "s3a", "s3b"),
  pon = c("nopon", "withpon")
)

# Build pairwise TF comparison
tf_wide <- primary %>%
  unite("condition", strategy, pon, sep = "_") %>%
  select(lib_id, condition, tumor_fraction) %>%
  pivot_wider(names_from = condition, values_from = tumor_fraction)

# Compute all pairwise Spearman correlations
conditions <- setdiff(names(tf_wide), "lib_id")
cor_mat <- matrix(NA, length(conditions), length(conditions),
                  dimnames = list(conditions, conditions))
for (i in seq_along(conditions)) {
  for (j in seq_along(conditions)) {
    cor_mat[i, j] <- cor(tf_wide[[conditions[i]]], tf_wide[[conditions[j]]],
                         method = "spearman")
  }
}

cor_df <- as.data.frame(as.table(cor_mat)) %>%
  rename(condition1 = Var1, condition2 = Var2, rho = Freq)

p_cormat <- ggplot(cor_df, aes(x = condition1, y = condition2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                        midpoint = 0.5, limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman\nrho",
       title = "Strategy-PoN concordance matrix (default-low-input)") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 13, face = "bold"))

# Panel B: Scatter — S1 vs S3A, nopon vs withpon
scatter_data <- primary %>%
  filter(pon == "withpon") %>%
  select(lib_id, strategy, tumor_fraction) %>%
  pivot_wider(names_from = strategy, values_from = tumor_fraction) %>%
  left_join(sample_meta %>% select(lib_id, treatment, tumor_vol_mm3), by = "lib_id")

p_scatter <- ggplot(scatter_data, aes(x = s1, y = s3a)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = treatment, size = tumor_vol_mm3), alpha = 0.8) +
  geom_text(aes(label = lib_id), vjust = -1.2, size = 3.5) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  scale_size_continuous(range = c(3, 8), name = expression("Tumor vol (mm"^3*")")) +
  labs(x = "Tumor fraction: S1 (direct)", y = "Tumor fraction: S3A (disambiguate)",
       color = "Treatment",
       title = "S1 vs S3A tumor fraction (with PoN)") +
  theme_white +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 13, face = "bold"))

combined <- p_cormat + p_scatter +
  plot_annotation(
    title = "Strategy Concordance Analysis",
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )

ggsave("plots/strategy_concordance.png", combined, width = 16, height = 7, dpi = 200)
cat("Saved plots/strategy_concordance.png\n")
