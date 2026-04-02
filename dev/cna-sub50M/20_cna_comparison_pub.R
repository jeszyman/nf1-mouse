## 20_cna_comparison_pub.R — Publication-quality S1 vs S3A CNA comparison
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

seg_all <- read_tsv("data/cna_segments.tsv", show_col_types = FALSE)

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  )

# Chromosome sizes for consistent x-axis
chr_order <- 1:22

# lib0626 and lib0627: paired S1 vs S3A comparison
for (lib in c("lib0626", "lib0627")) {
  seg_lib <- seg_all %>%
    filter(lib_id == lib, pon == "nopon", preset == "default-low-input",
           strategy %in% c("s1", "s3a")) %>%
    mutate(chr_num = as.numeric(str_replace(chr, "chr", ""))) %>%
    filter(!is.na(chr_num)) %>%
    mutate(strategy_label = ifelse(strategy == "s1",
      "S1: Direct alignment (contaminated)",
      "S3A: After disambiguation (clean)"))

  trt <- sample_meta$treatment[sample_meta$lib_id == lib]
  vol <- sample_meta$tumor_vol_mm3[sample_meta$lib_id == lib]

  # Get TF values
  rerun_path <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun/ichor_rerun_summary.tsv"
  ichor <- read_tsv(rerun_path, show_col_types = FALSE) %>%
    filter(str_detect(sample, lib), preset == "default-low-input", pon == "nopon")

  s1_tf <- ichor$tumor_fraction[ichor$strategy == "1"]
  s3a_tf <- ichor$tumor_fraction[ichor$strategy == "3a"]

  seg_lib <- seg_lib %>%
    mutate(strategy_label = case_when(
      strategy == "s1" ~ sprintf("S1: Direct (TF = %.2f)", s1_tf),
      strategy == "s3a" ~ sprintf("S3A: Disambiguated (TF = %.2f)", s3a_tf)
    ))

  # Color by CNA event
  event_colors <- c(
    NEUT = "#b0b0b0", AMP = "#d73027", GAIN = "#fc8d59",
    HETD = "#4575b4", HOMD = "#313695", HLAMP = "#a50026"
  )

  p <- ggplot(seg_lib, aes(x = start / 1e6, y = logR)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray80", linewidth = 0.3) +
    geom_hline(yintercept = c(-0.3, 0.3), linetype = "dotted", color = "gray70", linewidth = 0.3) +
    geom_segment(aes(xend = end / 1e6, yend = logR, color = event),
                 linewidth = 2.5, alpha = 0.85) +
    scale_color_manual(values = event_colors, na.value = "gray60") +
    facet_grid(strategy_label ~ chr_num, scales = "free_x", space = "free_x",
               switch = "y") +
    scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1, -0.5, 0, 0.5, 1)) +
    labs(x = NULL, y = "log2 ratio", color = "CNA event",
         title = sprintf("%s | %s | tumor vol: %s mm\u00B3",
                         lib, trt, format(vol, big.mark = ","))) +
    theme_pub +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 10, angle = 0, hjust = 0),
      panel.spacing.x = unit(0.05, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      legend.position = "bottom",
      legend.key.size = unit(0.8, "lines")
    )

  ggsave(sprintf("plots/cna_pub_%s.png", lib), p,
         width = 18, height = 6, dpi = 300)
  cat(sprintf("Saved plots/cna_pub_%s.png\n", lib))
}
