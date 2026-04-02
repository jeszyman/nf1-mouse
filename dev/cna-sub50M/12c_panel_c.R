## 12c_panel_c.R — Panel C: paired stacked bars showing S1 vs disambiguate overlap
if (!exists("sample_meta")) source("00_load_data.R")

overlap <- read_tsv("data/lib0626_read_overlap.tsv", show_col_types = FALSE)

# Bar 1: S1 MAPQ>=20 reads, colored by disambiguate classification
s1_total <- overlap$s1_q20
bar1 <- tibble(
  bar = "S1 MAPQ≥20",
  category = c("Human (confirmed)", "Ambiguous", "Mouse (contamination)"),
  reads = c(overlap$in_human, overlap$in_ambig, overlap$in_mouse)
) %>% mutate(pct = 100 * reads / s1_total)

# Bar 2: Disambiguate human reads, colored by S1 MAPQ status
disambig_human_total <- 10909485  # from disambig summary
in_s1_q20 <- overlap$in_human  # 2,733,031
not_in_s1_q20 <- disambig_human_total - in_s1_q20

bar2 <- tibble(
  bar = "Disambiguate human",
  category = c("Also S1 MAPQ≥20", "S1 MAPQ<20 (missed by S1)"),
  reads = c(in_s1_q20, not_in_s1_q20)
) %>% mutate(pct = 100 * reads / disambig_human_total)

# Combine
paired <- bind_rows(bar1, bar2) %>%
  mutate(
    bar = factor(bar, levels = c("S1 MAPQ≥20", "Disambiguate human")),
    category = factor(category, levels = c(
      "Mouse (contamination)", "Ambiguous", "Human (confirmed)",
      "S1 MAPQ<20 (missed by S1)", "Also S1 MAPQ≥20"
    ))
  )

# Label positions
paired_labels <- paired %>%
  arrange(bar, desc(category)) %>%
  group_by(bar) %>%
  mutate(ymax = cumsum(reads), ymin = ymax - reads, label_y = (ymin + ymax) / 2) %>%
  ungroup()

panel_c <- ggplot(paired, aes(x = bar, y = reads / 1e6, fill = category)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(data = paired_labels,
            aes(y = label_y / 1e6, label = sprintf("%.1f%%\n(%sM)", pct, format(round(reads/1e6,1), nsmall=1))),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(
    "Human (confirmed)" = "#5e3c99",
    "Ambiguous" = "#b2abd2",
    "Mouse (contamination)" = "#e66101",
    "Also S1 MAPQ≥20" = "#5e3c99",
    "S1 MAPQ<20 (missed by S1)" = "#fee0d2"
  )) +
  labs(x = NULL, y = "Unique reads (millions)", fill = NULL,
       title = "C. Read overlap: S1 vs Disambiguate (lib0626)",
       subtitle = sprintf("Jaccard index = 0.14 | Only %.0f%% of S1's 'good' reads are truly human",
                          100 * overlap$in_human / s1_total)) +
  theme_white +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text.x = element_text(size = 12))

ggsave("plots/panel_c_overlap.png", panel_c, width = 8, height = 7, dpi = 150)
cat("Saved plots/panel_c_overlap.png\n")
