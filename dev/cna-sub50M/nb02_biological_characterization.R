## nb02_biological_characterization.R
## Part 2: Biological characterization using S3A (the winning method)
## TF vs tumor biology, CNA profiles, fragment analysis, low-TF investigation
library(tidyverse)
library(patchwork)
library(ggrepel)

# --- Paths ---
dev_dir <- "/home/jeszyman/repos/nf1-mouse/.claude/worktrees/agent-a8c2930d/dev/cna-sub50M"
plot_dir <- file.path(dev_dir, "plots")

# --- Sample metadata ---
sample_meta <- tibble(
  lib_id    = paste0("lib0", 622:627),
  mouse_id  = paste0("mou000", 1:6),
  treatment = c(rep("vehicle", 3), rep("mirdametinib", 3)),
  cfdna_conc_pg_ul = c(467, 536, 111, 1720, 2380, 387),
  tumor_vol_mm3 = c(1823, 1410, 1048, 1902, 1868, 2551),
  days_post_implant = c(63, 63, 63, 72, 72, 72)
)

strategy_colors <- c(s1 = "#e66101", s3a = "#5e3c99", s3b = "#b2abd2")

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


# === Figure 4: TF vs tumor volume (S3A only) ===
ichor <- read_tsv(file.path(dev_dir, "data/ichor_rerun_full.tsv"),
                  show_col_types = FALSE) %>%
  filter(preset == "default-low-input", pon == "nopon") %>%
  mutate(strategy = case_when(
    strategy == "1" ~ "s1",
    strategy == "3a" ~ "s3a",
    strategy == "3b" ~ "s3b"
  ))

s3a_data <- ichor %>%
  filter(strategy == "s3a") %>%
  select(lib_id, tumor_fraction) %>%
  left_join(sample_meta, by = "lib_id") %>%
  mutate(tf_group = ifelse(tumor_fraction > 0.15, "High TF", "Low TF"))

cor_result <- cor.test(s3a_data$tumor_vol_mm3, s3a_data$tumor_fraction,
                        method = "spearman")

fig4 <- ggplot(s3a_data,
  aes(x = tumor_vol_mm3, y = tumor_fraction)) +
  geom_point(aes(color = treatment, shape = tf_group), size = 6) +
  geom_text_repel(aes(label = sprintf("%s\n(%s pg/uL)", lib_id, cfdna_conc_pg_ul)),
                  size = 3.5, max.overlaps = 10) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  scale_shape_manual(values = c("High TF" = 17, "Low TF" = 16)) +
  annotate("text", x = 1200, y = 0.55,
           label = sprintf("Spearman \u03C1 = %.2f\np = %.3f (n = 6)",
                           cor_result$estimate, cor_result$p.value),
           size = 4.5, hjust = 0) +
  labs(x = expression("Terminal tumor volume (mm"^3*")"),
       y = "S3A tumor fraction (ichorCNA)",
       color = "Treatment", shape = "TF group",
       title = "Tumor fraction vs tumor volume (S3A, default-low-input, no PoN)") +
  theme_pub +
  theme(legend.position = "right")

ggsave(file.path(plot_dir, "nb_fig4_tf_vs_tumor.png"), fig4,
       width = 10, height = 7, dpi = 300)
cat("Saved nb_fig4_tf_vs_tumor.png\n")


# === Figure 5: Waterfall summary ===
s3a_sorted <- s3a_data %>% arrange(desc(tumor_fraction))
sample_order <- s3a_sorted$lib_id

disambig <- read_tsv(file.path(dev_dir, "data/disambig_percentages.tsv"),
                     show_col_types = FALSE)

summary_data <- sample_meta %>%
  left_join(disambig %>% select(lib_id, human_pct, ambiguous_pct), by = "lib_id") %>%
  left_join(s3a_data %>% select(lib_id, s3a_tf = tumor_fraction, tf_group),
            by = "lib_id") %>%
  mutate(lib_id = factor(lib_id, levels = sample_order))

s1_data <- ichor %>%
  filter(strategy == "s1") %>%
  select(lib_id, s1_tf = tumor_fraction) %>%
  mutate(lib_id = factor(lib_id, levels = sample_order))

pa5 <- ggplot(summary_data, aes(x = lib_id, y = s3a_tf, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_point(data = s1_data, aes(x = lib_id, y = s1_tf),
             inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.5,
             color = "#e66101") +
  geom_text(aes(label = sprintf("%.0f%%", s3a_tf * 100)), vjust = -0.5, size = 4.5) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "Tumor fraction",
       fill = "Treatment",
       title = "A. S3A tumor fraction (bars) vs S1 (X marks)") +
  theme_pub +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"))

pb5 <- ggplot(summary_data, aes(x = lib_id, y = tumor_vol_mm3, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = expression("Tumor vol (mm"^3*")"),
       title = "B. Terminal tumor volume") +
  theme_pub + guides(fill = "none")

pc5 <- ggplot(summary_data, aes(x = lib_id, y = cfdna_conc_pg_ul, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "cfDNA (pg/uL)",
       title = "C. cfDNA concentration") +
  theme_pub + guides(fill = "none")

pd5 <- ggplot(summary_data, aes(x = lib_id, y = human_pct, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "Human reads (%)",
       title = "D. Human read fraction") +
  theme_pub + guides(fill = "none")

fig5 <- pa5 / (pb5 | pc5 | pd5) +
  plot_layout(heights = c(0.5, 0.5)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(file.path(plot_dir, "nb_fig5_waterfall.png"), fig5,
       width = 14, height = 10, dpi = 300)
cat("Saved nb_fig5_waterfall.png\n")


# === Figure 6: Per-chromosome CNA analysis (lib0626, lib0627) ===
seg_all <- read_tsv(file.path(dev_dir, "data/cna_segments.tsv"),
                    show_col_types = FALSE)

chr_medians <- seg_all %>%
  filter(pon == "nopon", preset == "default-low-input",
         strategy == "s3a") %>%
  mutate(chr_num = as.numeric(str_replace(chr, "chr", ""))) %>%
  filter(!is.na(chr_num)) %>%
  mutate(seg_len = end - start) %>%
  group_by(lib_id, chr_num) %>%
  summarise(median_logR = weighted.mean(logR, seg_len, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(sample_meta %>% select(lib_id, treatment), by = "lib_id") %>%
  left_join(s3a_data %>% select(lib_id, tf_group), by = "lib_id")

fig6 <- ggplot(chr_medians,
  aes(x = factor(chr_num), y = lib_id, fill = median_logR)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", median_logR)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                        midpoint = 0, limits = c(-0.8, 0.8)) +
  labs(x = "Chromosome", y = NULL, fill = "Median\nlog2R",
       title = "Per-chromosome CNA profile (S3A, all samples)") +
  theme_pub +
  theme(axis.text.y = element_text(size = 12))

ggsave(file.path(plot_dir, "nb_fig6_perchr_heatmap.png"), fig6,
       width = 16, height = 5, dpi = 300)
cat("Saved nb_fig6_perchr_heatmap.png\n")

cat("Part 2 script complete.\n")
