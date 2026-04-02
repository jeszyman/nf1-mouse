## 24_waterfall_summary.R — Integrated sample summary with TF waterfall
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
  filter(preset == "default-low-input", pon == "withpon")

disambig_pct <- read_tsv("data/disambig_percentages.tsv", show_col_types = FALSE)

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

# Sort samples by S3A TF (descending)
tf_s3a <- ichor_rerun %>%
  filter(strategy == "s3a") %>%
  arrange(desc(tumor_fraction)) %>%
  pull(lib_id)

# Build integrated summary
summary_data <- sample_meta %>%
  left_join(disambig_pct %>% select(lib_id, human_pct), by = "lib_id") %>%
  left_join(
    ichor_rerun %>%
      filter(strategy == "s3a") %>%
      select(lib_id, s3a_tf = tumor_fraction),
    by = "lib_id"
  ) %>%
  left_join(
    ichor_rerun %>%
      filter(strategy == "s1") %>%
      select(lib_id, s1_tf = tumor_fraction),
    by = "lib_id"
  ) %>%
  mutate(lib_id = factor(lib_id, levels = tf_s3a))

# Panel A: S3A TF waterfall
pa <- ggplot(summary_data, aes(x = lib_id, y = s3a_tf, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_point(aes(y = s1_tf), shape = 4, size = 4, stroke = 1.5, color = "#e66101") +
  geom_text(aes(label = sprintf("%.2f", s3a_tf)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "Tumor fraction",
       fill = "Treatment",
       title = "A. Tumor fraction (S3A with PoN, sorted by TF)",
       subtitle = "Bars = S3A | X marks = S1 (contaminated)") +
  theme_pub +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = "white", color = "gray80"))

# Panel B: Tumor volume
pb <- ggplot(summary_data, aes(x = lib_id, y = tumor_vol_mm3, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = expression("Tumor vol (mm"^3*")"), fill = NULL,
       title = "B. Terminal tumor volume") +
  theme_pub +
  guides(fill = "none")

# Panel C: Human read percentage
pc <- ggplot(summary_data, aes(x = lib_id, y = human_pct, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "Human reads (%)", fill = NULL,
       title = "C. Human read fraction (disambiguate)") +
  theme_pub +
  guides(fill = "none")

# Panel D: cfDNA concentration
pd <- ggplot(summary_data, aes(x = lib_id, y = cfdna_conc_pg_ul, fill = treatment)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "cfDNA (pg/uL)", fill = NULL,
       title = "D. cfDNA concentration") +
  theme_pub +
  guides(fill = "none")

combined <- pa / (pb | pc | pd) +
  plot_layout(heights = c(0.5, 0.5)) +
  plot_annotation(
    title = "Integrated Sample Summary",
    subtitle = "Samples sorted by S3A tumor fraction | WU-487 terminal bleeds",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave("plots/waterfall_summary.png", combined, width = 14, height = 10, dpi = 300)
cat("Saved plots/waterfall_summary.png\n")
