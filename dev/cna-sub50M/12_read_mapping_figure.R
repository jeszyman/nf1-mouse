## 12_read_mapping_figure.R — Comprehensive read mapping figure
## All values in read pairs for consistency with disambiguate output
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

# Convert flagstat (reads) to read pairs
flagstat_pairs <- flagstat_full %>%
  mutate(mapped_pairs = mapped / 2) %>%
  select(lib_id, strategy, mapped_pairs, treatment)

# Panel A (now first): Total human read pairs per strategy
strategy_data <- flagstat_pairs %>%
  mutate(strategy = factor(strategy, levels = c("s1", "s3b", "s3a"),
                           labels = c("S1: Direct", "S3B: Human+Ambig", "S3A: Human only")))

panel_a <- ggplot(strategy_data, aes(x = lib_id, y = mapped_pairs / 1e6, fill = strategy)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("S1: Direct" = "#e66101",
                                "S3B: Human+Ambig" = "#b2abd2",
                                "S3A: Human only" = "#5e3c99")) +
  labs(x = NULL, y = "Mapped read pairs (millions)", fill = "Strategy",
       title = "A. Human read pairs retained per strategy") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        legend.position = "right",
        plot.title = element_text(size = 13, face = "bold"))

# Panel B: Species assignment stacked bar with % labels inside
read_accounting <- disambig_data %>%
  select(lib_id, total_read_pairs, human_assigned, mouse_assigned, ambiguous) %>%
  left_join(sample_meta %>% select(lib_id, treatment, cfdna_conc_pg_ul), by = "lib_id")

panel_b_data <- read_accounting %>%
  select(lib_id, human_assigned, mouse_assigned, ambiguous) %>%
  pivot_longer(-lib_id, names_to = "category", values_to = "read_pairs") %>%
  mutate(category = factor(category,
    levels = c("mouse_assigned", "ambiguous", "human_assigned"),
    labels = c("Mouse", "Ambiguous", "Human"))) %>%
  left_join(read_accounting %>% select(lib_id, total_read_pairs), by = "lib_id") %>%
  mutate(pct = 100 * read_pairs / total_read_pairs)

# Compute label positions (midpoint of each segment in the stack)
panel_b_labels <- panel_b_data %>%
  arrange(lib_id, desc(category)) %>%
  group_by(lib_id) %>%
  mutate(ymax = cumsum(read_pairs),
         ymin = ymax - read_pairs,
         label_y = (ymin + ymax) / 2) %>%
  ungroup()

panel_b <- ggplot(panel_b_data, aes(x = lib_id, y = read_pairs / 1e6, fill = category)) +
  geom_col(alpha = 0.85) +
  geom_text(data = panel_b_labels,
            aes(x = lib_id, y = label_y / 1e6,
                label = sprintf("%s\n%.0f%%", category, pct)),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(Human = "#5e3c99", Ambiguous = "#b2abd2", Mouse = "#e66101")) +
  labs(x = NULL, y = "Read pairs (millions)", fill = NULL,
       title = "B. Disambiguate species classification (50M pairs input)") +
  theme_white +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.title = element_text(size = 13, face = "bold"))

# Panel C: S1 MAPQ>=20 true species composition (lib0626, exact read-name overlap)
# From local BAM analysis: 24.7% human, 68.6% mouse, 6.6% ambiguous
overlap <- read_tsv("data/lib0626_read_overlap.tsv", show_col_types = FALSE)

panel_c_data <- overlap %>%
  select(lib_id, in_human, in_mouse, in_ambig) %>%
  pivot_longer(-lib_id, names_to = "origin", values_to = "reads") %>%
  mutate(
    origin = factor(origin,
      levels = c("in_mouse", "in_ambig", "in_human"),
      labels = c("Mouse (contamination)", "Ambiguous", "Human (true)")),
    total = overlap$s1_q20[1],
    pct = 100 * reads / total
  )

# Label positions
panel_c_labels <- panel_c_data %>%
  arrange(lib_id, desc(origin)) %>%
  group_by(lib_id) %>%
  mutate(ymax = cumsum(reads), ymin = ymax - reads, label_y = (ymin + ymax) / 2) %>%
  ungroup()

panel_c <- ggplot(panel_c_data, aes(x = lib_id, y = reads / 1e6, fill = origin)) +
  geom_col(alpha = 0.85) +
  geom_text(data = panel_c_labels,
            aes(x = lib_id, y = label_y / 1e6,
                label = sprintf("%s\n%.1f%%", origin, pct)),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Human (true)" = "#5e3c99",
                                "Ambiguous" = "#b2abd2",
                                "Mouse (contamination)" = "#e66101")) +
  labs(x = NULL, y = "Read pairs (millions)", fill = NULL,
       title = "C. True species of S1 MAPQ≥20 reads (lib0626, exact overlap)",
       subtitle = "68.6% of S1's 'good' reads are mouse contamination") +
  theme_white +
  guides(fill = "none") +
  theme(axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "#e66101"))

# Combine all three
combined <- (panel_a | panel_b) / panel_c +
  plot_layout(heights = c(1, 0.6)) +
  plot_annotation(
    title = "PDX Read Handling: Sub-50M Cohort Read Mapping",
    subtitle = "6 WU-487 Chr8q-WT terminal bleeds (3 vehicle, 3 mirdametinib) | 50M read pairs each",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

ggsave("plots/read_mapping_comprehensive.png", combined,
       width = 14, height = 12, dpi = 150)

cat("Figure saved: plots/read_mapping_comprehensive.png\n")
