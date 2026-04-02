## 03_alignment_metrics.R — Alignment-level comparison across strategies
if (!exists("sample_meta")) source("00_load_data.R")

# --- Reads retained per strategy ---
reads_summary <- flagstat_full %>%
  select(lib_id, strategy, total, mapped, duplicates, properly_paired, treatment) %>%
  mutate(
    pct_mapped = 100 * mapped / total,
    pct_dup = 100 * duplicates / total,
    unique_molecules = mapped - duplicates,
    strategy_label = strategy_labels[strategy]
  )

write_tsv(reads_summary, "data/alignment_metrics.tsv")

# --- Total mapped reads by strategy ---
p_mapped <- ggplot(reads_summary, aes(x = lib_id, y = mapped / 1e6,
                                       fill = strategy)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = NULL, y = "Mapped reads (millions)", fill = "Strategy",
       title = "Mapped reads by strategy") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/mapped_reads.png", p_mapped, width = 10, height = 5, dpi = 150)

# --- Unique molecules (mapped - duplicates) ---
p_unique <- ggplot(reads_summary, aes(x = lib_id, y = unique_molecules / 1e6,
                                       fill = strategy)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = NULL, y = "Unique molecules (millions)", fill = "Strategy",
       title = "Unique molecules (mapped - duplicates) by strategy") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/unique_molecules.png", p_unique, width = 10, height = 5, dpi = 150)

# --- Duplication rate ---
p_dup <- ggplot(reads_summary, aes(x = lib_id, y = pct_dup, fill = strategy)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = NULL, y = "Duplication rate (%)", fill = "Strategy",
       title = "Duplication rate by strategy") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/duplication_rate.png", p_dup, width = 10, height = 5, dpi = 150)

# --- Reads retained relative to S1 ---
s1_mapped <- reads_summary %>%
  filter(strategy == "s1") %>%
  select(lib_id, s1_mapped = mapped)

retention <- reads_summary %>%
  left_join(s1_mapped, by = "lib_id") %>%
  mutate(retention_vs_s1 = 100 * mapped / s1_mapped)

p_retention <- ggplot(retention %>% filter(strategy != "s1"),
  aes(x = lib_id, y = retention_vs_s1, fill = strategy)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = NULL, y = "% of S1 mapped reads retained",
       title = "Read retention relative to direct alignment (S1)") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/retention_vs_s1.png", p_retention, width = 10, height = 5, dpi = 150)

cat("\nScript 03 complete.\n")
print(reads_summary %>% select(lib_id, strategy, mapped, unique_molecules, pct_dup))
