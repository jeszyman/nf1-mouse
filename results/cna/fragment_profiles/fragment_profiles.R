library(tidyverse)
source("/home/jeszyman/repos/science/R/figure_schema.R")
library(patchwork)

d <- read_tsv("fragment_profiles.tsv", show_col_types = FALSE) %>%
  filter(insert_size >= 50, insert_size <= 400, strategy %in% c("s1", "s3a")) %>%
  mutate(treatment = str_to_title(treatment))

d_pct <- d %>%
  complete(insert_size = 50:400, nesting(lib_id, strategy, treatment),
           fill = list(pairs_total = 0)) %>%
  group_by(lib_id, strategy) %>%
  mutate(pct = 100 * pairs_total / sum(pairs_total)) %>%
  ungroup() %>%
  mutate(
    tf_group = ifelse(lib_id %in% c("lib0626", "lib0627"), "High TF", "Low TF"),
    strategy_label = case_when(
      strategy == "s1"  ~ "S1: Direct (hg38)",
      strategy == "s3a" ~ "S3A: Human only"
    )
  )

strategy_colors <- c("S1: Direct (hg38)" = "#e66101", "S3A: Human only" = "#5e3c99")
shared_y <- scale_y_continuous(breaks = seq(0, 6, by = 2), limits = c(0, 6.5))

base_line_theme <- list(
  theme_scifig(),
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 13, face = "bold"))
)

pa <- ggplot(filter(d_pct, tf_group == "High TF"),
             aes(x = insert_size, y = pct, color = strategy_label)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_vline(xintercept = 167, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  annotate("text", x = 200, y = 6, label = "167 bp", size = 3.5, color = "gray40") +
  scale_color_manual(values = strategy_colors) +
  shared_y +
  facet_wrap(~lib_id, ncol = 2) +
  labs(x = NULL, y = "% of total fragments", color = NULL,
       subtitle = "A. High-TF samples: nucleosomal peak emerges after disambiguation") +
  base_line_theme

pb <- ggplot(filter(d_pct, tf_group == "Low TF", lib_id %in% c("lib0622", "lib0625")),
             aes(x = insert_size, y = pct, color = strategy_label)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_vline(xintercept = 167, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  annotate("text", x = 200, y = 6, label = "167 bp", size = 3.5, color = "gray40") +
  scale_color_manual(values = strategy_colors) +
  shared_y +
  facet_wrap(~lib_id, ncol = 2) +
  labs(x = "Insert size (bp)", y = "% of total fragments", color = NULL,
       subtitle = "B. Low-TF samples: no nucleosomal peak regardless of strategy") +
  base_line_theme +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.key.width = unit(1.8, "cm"),
        legend.box.margin = margin(t = 5))

legend_text <- sprintf(
  "Mouse-read disambiguation reveals a cfDNA nucleosomal peak in high-tumor-fraction PDX plasma. Fragment length profiles from WU-487 terminal bleeds (n = %d), 50M-read subsamples. (A) High-TF samples show a nucleosomal peak at ~145-167 bp after mouse-read removal (S3A), masked in S1 by genomic DNA at 167 bp. (B) Low-TF samples show no nucleosomal peak regardless of strategy. Dotted line marks 167 bp.",
  n_distinct(d$lib_id))

p <- pa / pb +
  plot_layout(heights = c(0.5, 0.5))

# Wrap as single ggplot with caption via wrap_plots
final <- wrap_plots(p) +
  plot_annotation(
    caption = str_wrap(legend_text, width = round(10 * 11.5)),
    theme = theme(
      plot.caption = element_text(size = 12, hjust = 0, color = "gray30",
                                  lineheight = 1.3,
                                  margin = margin(t = 20)),
      plot.caption.position = "plot",
      plot.margin = margin(10, 10, 10, 10, "pt"),
      plot.background = element_rect(fill = "white", color = NA))
  )

save_plot(final, "fragment_profiles", w = 10, h = 10)
