library(tidyverse)
source("/home/jeszyman/repos/science/R/figure_schema.R")

d <- read_tsv("tf_vs_tumor.tsv", show_col_types = FALSE) %>%
  filter(strategy %in% c("s1", "s3a")) %>%
  mutate(
    strategy_label = case_when(
      strategy == "s1"  ~ "S1: Direct (hg38)",
      strategy == "s3a" ~ "S3A: Human Only"
    ),
    treatment = str_to_title(treatment)
  )

strategy_colors <- c("S1: Direct (hg38)" = "#e66101", "S3A: Human Only" = "#5e3c99")

wt <- wilcox.test(
  d$tumor_fraction[d$strategy == "s1"],
  d$tumor_fraction[d$strategy == "s3a"],
  paired = TRUE
)

legend_text <- sprintf(
  "Mouse-read disambiguation reveals higher ctDNA tumor fraction in high-TF PDX plasma. Paired ichorCNA tumor fraction estimates from S1 (direct hg38 alignment) vs S3A (mouse-read removal) across %d WU-487 terminal bleeds. Lines connect paired samples. S3A unmasks tumor fraction suppressed by mouse-read contamination in S1, with the largest gains in high-TF samples (paired Wilcoxon %s, n = %d).",
  n_distinct(d$lib_id), format_p(wt$p.value), n_distinct(d$lib_id))

p <- ggplot(d, aes(x = strategy_label, y = tumor_fraction)) +
  geom_line(aes(group = lib_id), color = "grey60", linewidth = 0.4) +
  geom_point(aes(fill = strategy_label), size = 3, shape = 21,
             color = "black", stroke = 0.7) +
  annotate("segment", x = 1, xend = 2, y = 0.72, yend = 0.72, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = 0.70, yend = 0.72, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = 0.70, yend = 0.72, linewidth = 0.4) +
  annotate("text", x = 1.5, y = 0.74, label = format_p(wt$p.value), size = 3.5) +
  scale_fill_manual(values = strategy_colors) +
  scale_y_continuous(breaks = seq(0, 0.7, by = 0.1), limits = c(0, 0.78)) +
  labs(x = NULL,
       y = "Tumor fraction (ichorCNA)",
       caption = str_wrap(legend_text, width = round(5 * 15))) +
  theme_scifig() +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)),
        plot.caption = element_text(size = 10, hjust = 0, color = "gray30",
                                    lineheight = 1.3, margin = margin(t = 20)),
        plot.caption.position = "plot",
        plot.margin = margin(10, 10, 10, 10, "pt"))

save_plot(p, "tf_vs_tumor_fig", w = 5, h = 5.5)
