## 18_improved_figures.R — Publication-quality versions of key figures
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)
library(ggrepel)

rerun_path <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun/ichor_rerun_summary.tsv"
ichor_rerun <- read_tsv(rerun_path, show_col_types = FALSE) %>%
  mutate(lib_id = str_extract(sample, "lib\\d+"),
         strategy = case_when(
           strategy == "1" ~ "s1",
           strategy == "3a" ~ "s3a",
           strategy == "3b" ~ "s3b"
         )) %>%
  left_join(sample_meta, by = "lib_id")

# Publication theme
theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    strip.text = element_text(size = 12, face = "bold")
  )

# ===== FIGURE 1: Species classification summary =====
disambig_pct <- read_tsv("data/disambig_percentages.tsv", show_col_types = FALSE)

fig1_data <- disambig_pct %>%
  select(lib_id, human_pct, mouse_pct, ambiguous_pct) %>%
  pivot_longer(-lib_id, names_to = "species", values_to = "pct") %>%
  mutate(
    species = factor(species,
      levels = c("mouse_pct", "ambiguous_pct", "human_pct"),
      labels = c("Mouse", "Ambiguous", "Human")),
    lib_label = lib_id
  ) %>%
  left_join(sample_meta %>% select(lib_id, treatment), by = "lib_id")

fig1 <- ggplot(fig1_data, aes(x = lib_id, y = pct, fill = species)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(Human = "#5e3c99", Ambiguous = "#b2abd2", Mouse = "#e66101")) +
  labs(x = NULL, y = "Read pairs (%)", fill = "Species assignment",
       title = "Disambiguate species classification",
       subtitle = "50M read-pair subsets | WU-487 terminal bleeds") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/fig1_species_classification.png", fig1, width = 8, height = 6, dpi = 300)

# ===== FIGURE 2: S1 contamination audit (lib0626) =====
overlap <- read_tsv("data/lib0626_read_overlap.tsv", show_col_types = FALSE)

fig2_data <- tibble(
  category = c("Human\n(confirmed)", "Ambiguous", "Mouse\n(contamination)"),
  reads = c(overlap$in_human, overlap$in_ambig, overlap$in_mouse),
  pct = c(overlap$pct_human, overlap$pct_ambig, overlap$pct_mouse)
) %>%
  mutate(category = factor(category,
    levels = c("Mouse\n(contamination)", "Ambiguous", "Human\n(confirmed)")))

fig2 <- ggplot(fig2_data, aes(x = "S1 MAPQ\u226520\nreads", y = reads / 1e6, fill = category)) +
  geom_col(alpha = 0.85, width = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%\n(%.1fM)", pct, reads / 1e6)),
            position = position_stack(vjust = 0.5),
            size = 4.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(
    "Human\n(confirmed)" = "#5e3c99",
    "Ambiguous" = "#b2abd2",
    "Mouse\n(contamination)" = "#e66101")) +
  labs(x = NULL, y = "Unique reads (millions)", fill = NULL,
       title = "True species of S1 MAPQ\u226520 reads",
       subtitle = "lib0626 | 68.6% of high-quality reads are mouse contamination") +
  theme_pub +
  guides(fill = "none") +
  theme(axis.text.x = element_text(size = 14))

ggsave("plots/fig2_s1_contamination.png", fig2, width = 5, height = 7, dpi = 300)

# ===== FIGURE 3: TF comparison across strategies =====
primary <- ichor_rerun %>%
  filter(preset == "default-low-input", pon == "withpon")

fig3_data <- primary %>%
  mutate(strategy_label = case_when(
    strategy == "s1" ~ "S1: Direct",
    strategy == "s3a" ~ "S3A: Human only",
    strategy == "s3b" ~ "S3B: Human+Ambig"
  )) %>%
  mutate(strategy_label = factor(strategy_label,
    levels = c("S1: Direct", "S3B: Human+Ambig", "S3A: Human only")))

fig3 <- ggplot(fig3_data, aes(x = lib_id, y = tumor_fraction, fill = strategy_label)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(
    "S1: Direct" = "#e66101",
    "S3B: Human+Ambig" = "#b2abd2",
    "S3A: Human only" = "#5e3c99")) +
  geom_hline(yintercept = 0.1, linetype = "dotted", color = "gray50") +
  annotate("text", x = 0.5, y = 0.105, label = "TF = 0.10", hjust = 0, size = 3, color = "gray50") +
  labs(x = NULL, y = "ichorCNA tumor fraction", fill = "Strategy",
       title = "Tumor fraction by strategy",
       subtitle = "Default-low-input preset | with PoN | autosomes only") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("plots/fig3_tf_comparison.png", fig3, width = 10, height = 7, dpi = 300)

# ===== FIGURE 4: TF vs tumor volume (S3A only, corrected) =====
s3a_data <- primary %>% filter(strategy == "s3a")
cor_result <- cor.test(s3a_data$tumor_vol_mm3, s3a_data$tumor_fraction, method = "spearman")

fig4 <- ggplot(s3a_data, aes(x = tumor_vol_mm3, y = tumor_fraction)) +
  geom_smooth(method = "lm", se = TRUE, color = "#5e3c99", alpha = 0.15, linewidth = 0.8) +
  geom_point(aes(fill = treatment), size = 5, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(aes(label = lib_id), size = 4, box.padding = 0.5) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  annotate("text", x = 1100, y = 0.62,
           label = sprintf("Spearman \u03C1 = %.2f\np = %.2f (n = 6)",
                           cor_result$estimate, cor_result$p.value),
           size = 4.5, hjust = 0, fontface = "italic") +
  labs(x = expression(bold("Terminal tumor volume (mm"^3*")")),
       y = expression(bold("ichorCNA tumor fraction (S3A)")),
       fill = "Treatment",
       title = "Tumor fraction vs tumor volume after disambiguation",
       subtitle = "S3A with PoN | Corrected from body weight to actual tumor volume") +
  theme_pub +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"))

ggsave("plots/fig4_tf_vs_tumor.png", fig4, width = 8, height = 7, dpi = 300)

# ===== FIGURE 5: Comprehensive multi-panel summary =====
# Panel A: Species fractions (stacked bar)
pa <- fig1 + labs(title = "A. Species classification") +
  theme(plot.title = element_text(size = 13))

# Panel B: TF by strategy
pb <- fig3 + labs(title = "B. Tumor fraction by strategy") +
  theme(plot.title = element_text(size = 13))

# Panel C: TF vs tumor vol
pc <- fig4 + labs(title = "C. TF vs tumor volume (S3A)") +
  theme(plot.title = element_text(size = 13))

# Panel D: S1 contamination
pd <- fig2 + labs(title = "D. S1 contamination (lib0626)") +
  theme(plot.title = element_text(size = 13))

summary_fig <- (pa | pb) / (pd | pc) +
  plot_annotation(
    title = "PDX Read Handling: Sub-50M Cohort Summary",
    subtitle = "6 WU-487 terminal bleeds | BISCUIT aligner | ichorCNA default-low-input with PoN",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 13, color = "gray30"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave("plots/fig5_summary_multipanel.png", summary_fig, width = 18, height = 14, dpi = 300)
cat("All improved figures saved.\n")
