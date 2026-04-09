## nb01_methods_comparison.R
## Part 1: S1 vs S3A/S3B read assignment methods comparison
## Generates figures for the methods comparison section of the notebook
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
  tumor_vol_mm3 = c(1823, 1410, 1048, 1902, 1868, 2551)
)

strategy_colors <- c(s1 = "#e66101", s3a = "#5e3c99", s3b = "#b2abd2")
strategy_labels <- c(
  s1 = "S1: Direct (hg38)",
  s3a = "S3A: Disambiguate (human only)",
  s3b = "S3B: Disambiguate (human + ambig)"
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

# === Figure 1: Read assignment and alignment summary ===
disambig <- read_tsv(file.path(dev_dir, "data/disambig_percentages.tsv"),
                     show_col_types = FALSE)

align <- read_tsv(file.path(dev_dir, "data/alignment_metrics.tsv"),
                  show_col_types = FALSE)

# Panel A: Disambiguate read classification
disambig_long <- disambig %>%
  select(lib_id, treatment, human_pct, mouse_pct, ambiguous_pct) %>%
  pivot_longer(cols = c(human_pct, mouse_pct, ambiguous_pct),
               names_to = "category", values_to = "pct") %>%
  mutate(category = factor(category,
    levels = c("mouse_pct", "ambiguous_pct", "human_pct"),
    labels = c("Mouse", "Ambiguous", "Human")))

pa <- ggplot(disambig_long,
  aes(x = lib_id, y = pct, fill = category)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(Mouse = "#999999", Ambiguous = "#f0ad4e",
                                Human = "#5e3c99")) +
  geom_text(data = disambig,
    aes(x = lib_id, y = 102,
        label = sprintf("%.1f%%", human_pct)),
    inherit.aes = FALSE, size = 3.5, fontface = "bold", color = "#5e3c99") +
  labs(x = NULL, y = "% of total reads", fill = "Classification",
       title = "A. Disambiguate read classification") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Mapped read pairs by strategy
align_pairs <- align %>%
  mutate(mapped_pairs = mapped / 2,
         strategy_label = strategy_labels[strategy]) %>%
  mutate(strategy_label = factor(strategy_label, levels = strategy_labels))

pb <- ggplot(align_pairs,
  aes(x = lib_id, y = mapped_pairs / 1e6, fill = strategy)) +
  geom_col(position = position_dodge(0.8), alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = NULL, y = "Mapped read pairs (millions)", fill = "Strategy",
       title = "B. Mapped read pairs by strategy") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

fig1 <- pa / pb +
  plot_layout(heights = c(0.45, 0.55)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(file.path(plot_dir, "nb_fig1_read_assignment.png"), fig1,
       width = 10, height = 10, dpi = 300)
cat("Saved nb_fig1_read_assignment.png\n")


# === Figure 2: TF comparison across strategies ===
ichor <- read_tsv(file.path(dev_dir, "data/ichor_rerun_full.tsv"),
                  show_col_types = FALSE) %>%
  filter(preset == "default-low-input", pon == "nopon") %>%
  mutate(strategy = case_when(
    strategy == "1" ~ "s1",
    strategy == "3a" ~ "s3a",
    strategy == "3b" ~ "s3b"
  ))

# Wide format for paired comparison
ichor_wide <- ichor %>%
  select(lib_id, strategy, tumor_fraction) %>%
  pivot_wider(names_from = strategy, values_from = tumor_fraction) %>%
  left_join(sample_meta, by = "lib_id")

pa2 <- ggplot(ichor %>% left_join(sample_meta, by = "lib_id"),
  aes(x = lib_id, y = tumor_fraction, fill = strategy)) +
  geom_col(position = position_dodge(0.8), alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = NULL, y = "Tumor fraction (ichorCNA)",
       fill = "Strategy",
       title = "A. ichorCNA tumor fraction by read assignment strategy") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

pb2 <- ggplot(ichor_wide, aes(x = s1, y = s3a)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(aes(color = treatment), size = 5) +
  geom_text_repel(aes(label = lib_id), size = 4) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "S1 tumor fraction", y = "S3A tumor fraction",
       color = "Treatment",
       title = "B. S1 vs S3A: disambiguation reveals suppressed TF") +
  theme_pub

fig2 <- pa2 + pb2 +
  plot_layout(widths = c(0.55, 0.45)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(file.path(plot_dir, "nb_fig2_tf_comparison.png"), fig2,
       width = 14, height = 6, dpi = 300)
cat("Saved nb_fig2_tf_comparison.png\n")


# === Figure 3: CNA profiles S1 vs S3A (lib0626 and lib0627) ===
seg_all <- read_tsv(file.path(dev_dir, "data/cna_segments.tsv"),
                    show_col_types = FALSE)

event_colors <- c(
  NEUT = "#b0b0b0", AMP = "#d73027", GAIN = "#fc8d59",
  HETD = "#4575b4", HOMD = "#313695", HLAMP = "#a50026"
)

for (lib in c("lib0626", "lib0627")) {
  seg_lib <- seg_all %>%
    filter(lib_id == lib, pon == "nopon", preset == "default-low-input",
           strategy %in% c("s1", "s3a")) %>%
    mutate(chr_num = as.numeric(str_replace(chr, "chr", ""))) %>%
    filter(!is.na(chr_num))

  s1_tf <- ichor$tumor_fraction[ichor$lib_id == lib & ichor$strategy == "s1"]
  s3a_tf <- ichor$tumor_fraction[ichor$lib_id == lib & ichor$strategy == "s3a"]
  trt <- sample_meta$treatment[sample_meta$lib_id == lib]
  vol <- sample_meta$tumor_vol_mm3[sample_meta$lib_id == lib]

  seg_lib <- seg_lib %>%
    mutate(strategy_label = case_when(
      strategy == "s1" ~ sprintf("S1: Direct (TF = %.1f%%)", s1_tf * 100),
      strategy == "s3a" ~ sprintf("S3A: Disambiguated (TF = %.1f%%)", s3a_tf * 100)
    ))

  p <- ggplot(seg_lib, aes(x = start / 1e6, y = logR)) +
    geom_hline(yintercept = 0, color = "gray80", linewidth = 0.3) +
    geom_hline(yintercept = c(-0.3, 0.3), linetype = "dotted",
               color = "gray70", linewidth = 0.3) +
    geom_segment(aes(xend = end / 1e6, yend = logR, color = event),
                 linewidth = 2.5, alpha = 0.85) +
    scale_color_manual(values = event_colors, na.value = "gray60") +
    facet_grid(strategy_label ~ chr_num, scales = "free_x", space = "free_x") +
    scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1, 0, 1)) +
    labs(x = NULL, y = "log2 ratio", color = "CNA event",
         title = sprintf("%s (%s, tumor vol: %s mm\u00B3)",
                         lib, trt, format(vol, big.mark = ","))) +
    theme_pub +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 10, angle = 0, hjust = 0),
      panel.spacing.x = unit(0.05, "lines"),
      legend.position = "bottom"
    )

  ggsave(file.path(plot_dir, sprintf("nb_fig3_cna_%s.png", lib)), p,
         width = 18, height = 6, dpi = 300)
  cat(sprintf("Saved nb_fig3_cna_%s.png\n", lib))
}

cat("Part 1 script complete.\n")
