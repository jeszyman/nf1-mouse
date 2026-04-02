## 15_tf_vs_tumor_corrected.R — TF vs tumor volume with corrected volumes + regression
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)
library(ggrepel)

rerun_path <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun/ichor_rerun_summary.tsv"
ichor_rerun <- read_tsv(rerun_path, show_col_types = FALSE) %>%
  mutate(lib_id = str_extract(sample, "lib\\d+"),
         strategy = paste0("s", str_replace(as.character(strategy), "([0-9])([ab])?", "\\1\\2"))) %>%
  left_join(sample_meta, by = "lib_id")

# Fix strategy labels
ichor_rerun <- ichor_rerun %>%
  mutate(strategy = case_when(
    strategy == "s1" ~ "s1",
    strategy == "s3a" ~ "s3a",
    strategy == "s3b" ~ "s3b",
    TRUE ~ strategy
  ))

# Focus: default-low-input, withpon (production-recommended)
primary <- ichor_rerun %>%
  filter(preset == "default-low-input", pon == "withpon")

# Panel A: TF vs tumor volume, colored by strategy
p_tf_vol <- ggplot(primary,
  aes(x = tumor_vol_mm3, y = tumor_fraction,
      color = strategy, shape = treatment)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = lib_id), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  scale_shape_manual(values = c(vehicle = 16, mirdametinib = 17)) +
  labs(x = expression("Terminal tumor volume (mm"^3*")"),
       y = "ichorCNA tumor fraction",
       color = "Strategy", shape = "Treatment",
       title = "A. Tumor fraction vs tumor volume (default-low-input, with PoN)") +
  theme_white +
  theme(legend.position = "right",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 13, face = "bold"))

# Panel B: S3A only with linear regression
s3a_data <- primary %>% filter(strategy == "s3a")
cor_result <- cor.test(s3a_data$tumor_vol_mm3, s3a_data$tumor_fraction, method = "spearman")

p_s3a_reg <- ggplot(s3a_data,
  aes(x = tumor_vol_mm3, y = tumor_fraction)) +
  geom_point(aes(color = treatment), size = 5) +
  geom_smooth(method = "lm", se = TRUE, color = "#5e3c99", alpha = 0.2) +
  geom_text_repel(aes(label = lib_id), size = 3.5) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  annotate("text", x = 1200, y = 0.6,
           label = sprintf("Spearman rho = %.2f\np = %.3f",
                           cor_result$estimate, cor_result$p.value),
           size = 4, hjust = 0) +
  labs(x = expression("Terminal tumor volume (mm"^3*")"),
       y = "ichorCNA tumor fraction (S3A)",
       color = "Treatment",
       title = "B. S3A: Disambiguated TF correlates with tumor burden") +
  theme_white +
  theme(legend.position = "right",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 13, face = "bold"))

combined <- p_tf_vol / p_s3a_reg +
  plot_annotation(
    title = "Tumor Fraction vs Tumor Volume (Corrected Metadata)",
    subtitle = "Previous analysis used body weight (23-27 g) instead of tumor volume (1,048-2,551 mm\u00B3)",
    theme = theme(plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "firebrick"))
  )

ggsave("plots/tf_vs_tumor_corrected.png", combined, width = 10, height = 12, dpi = 200)

cat("Spearman correlation (S3A TF vs tumor volume):\n")
cat(sprintf("  rho = %.3f, p = %.4f\n", cor_result$estimate, cor_result$p.value))
cat("Saved plots/tf_vs_tumor_corrected.png\n")
