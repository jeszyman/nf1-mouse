## 17_cna_sidebyside.R — Side-by-side CNA profiles: S1 vs S3A for high-TF samples
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

seg_all <- read_tsv("data/cna_segments.tsv", show_col_types = FALSE)

# Focus on lib0626 and lib0627 (high TF after disambiguation)
for (lib in c("lib0626", "lib0627")) {
  seg_lib <- seg_all %>%
    filter(lib_id == lib, pon == "nopon", preset == "default-low-input") %>%
    mutate(chr_num = as.numeric(str_replace(chr, "chr", ""))) %>%
    filter(!is.na(chr_num)) %>%
    mutate(strategy_label = case_when(
      strategy == "s1" ~ "S1: Direct alignment",
      strategy == "s3a" ~ "S3A: Disambiguate (human only)",
      strategy == "s3b" ~ "S3B: Disambiguate (human + ambig)"
    ))

  if (nrow(seg_lib) == 0) next

  trt <- sample_meta$treatment[sample_meta$lib_id == lib]
  vol <- sample_meta$tumor_vol_mm3[sample_meta$lib_id == lib]

  p <- ggplot(seg_lib, aes(x = start / 1e6, y = logR)) +
    geom_segment(aes(xend = end / 1e6, yend = logR, color = event),
                 linewidth = 2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(
      values = c(NEUT = "gray50", AMP = "#d73027", GAIN = "#fc8d59",
                 HETD = "#4575b4", HOMD = "#313695", HLAMP = "#a50026"),
      na.value = "gray70"
    ) +
    facet_grid(strategy_label ~ chr_num, scales = "free_x", space = "free_x") +
    labs(x = "Position (Mb)", y = "log2 ratio", color = "CNA event",
         title = sprintf("%s (%s, tumor vol: %s mm\u00B3)", lib, trt, format(vol, big.mark = ",")),
         subtitle = "ichorCNA default-low-input, no PoN, autosomes only") +
    theme_white +
    theme(axis.text.x = element_blank(),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 9),
          panel.spacing.x = unit(0.1, "lines"),
          panel.spacing.y = unit(0.3, "lines"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "bottom")

  ggsave(sprintf("plots/cna_sidebyside_%s.png", lib), p,
         width = 18, height = 8, dpi = 200)
  cat(sprintf("Saved plots/cna_sidebyside_%s.png\n", lib))
}
