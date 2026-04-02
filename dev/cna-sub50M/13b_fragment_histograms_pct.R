## 13b_fragment_histograms_pct.R — Fragment length histograms normalized to % of total
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

frag <- read_tsv("data/fragment_lengths.tsv", show_col_types = FALSE) %>%
  filter(insert_size >= 50, insert_size <= 400) %>%
  mutate(strategy = factor(strategy, levels = c("s1", "s3a", "s3b"),
                           labels = c("S1: Direct", "S3A: Human only", "S3B: Human+Ambig")))

# Normalize each sample×strategy to percent of its own total
frag_pct <- frag %>%
  group_by(lib_id, strategy) %>%
  mutate(pct = 100 * pairs_total / sum(pairs_total)) %>%
  ungroup()

panels <- list()
for (lib in libs) {
  d <- frag_pct %>% filter(lib_id == lib)
  trt <- sample_meta$treatment[sample_meta$lib_id == lib]
  conc <- sample_meta$cfdna_conc_pg_ul[sample_meta$lib_id == lib]

  panels[[lib]] <- ggplot(d, aes(x = insert_size, y = pct, color = strategy)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    scale_color_manual(values = c("S1: Direct" = "#e66101",
                                   "S3A: Human only" = "#5e3c99",
                                   "S3B: Human+Ambig" = "#b2abd2")) +
    labs(x = "Insert size (bp)", y = "% of total",
         title = sprintf("%s (%s, %d pg/µL)", lib, trt, conc)) +
    theme_white +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.position = "none")
}

legend_plot <- ggplot(frag_pct %>% filter(lib_id == libs[1]),
                      aes(x = insert_size, y = pct, color = strategy)) +
  geom_line() +
  scale_color_manual(values = c("S1: Direct" = "#e66101",
                                 "S3A: Human only" = "#5e3c99",
                                 "S3B: Human+Ambig" = "#b2abd2")) +
  theme_white +
  theme(legend.text = element_text(size = 12))
shared_legend <- cowplot::get_legend(legend_plot)

combined <- (wrap_plots(panels, ncol = 2) | shared_legend) +
  plot_layout(widths = c(6, 1)) +
  plot_annotation(
    title = "Fragment Length Distributions (% of total per strategy)",
    subtitle = "Each strategy normalized to its own total | 1% subsample | WU-487 terminal bleeds",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

ggsave("plots/fragment_histograms_pct.png", combined, width = 14, height = 12, dpi = 150)
cat("Figure saved: plots/fragment_histograms_pct.png\n")
