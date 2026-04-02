## 02_read_assignment.R — Disambiguate species assignment analysis
if (!exists("sample_meta")) source("00_load_data.R")

# --- Disambiguate fractions per sample ---
disambig_long <- disambig_data %>%
  select(lib_id, human_assigned, mouse_assigned, ambiguous) %>%
  pivot_longer(-lib_id, names_to = "category", values_to = "read_pairs") %>%
  left_join(sample_meta %>% select(lib_id, treatment, cfdna_conc_pg_ul), by = "lib_id") %>%
  mutate(category = factor(category, levels = c("human_assigned", "ambiguous", "mouse_assigned")))

p_disambig <- ggplot(disambig_long, aes(x = lib_id, y = read_pairs / 1e6, fill = category)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c(human_assigned = "#5e3c99", ambiguous = "#b2abd2",
                                mouse_assigned = "#e66101"),
                    labels = c("Human", "Ambiguous", "Mouse")) +
  labs(x = NULL, y = "Read pairs (millions)", fill = "Assignment",
       title = "Species read assignment (50M read-pair subsets)") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/disambig_fractions.png", p_disambig, width = 8, height = 5, dpi = 150)

# --- Disambiguate fractions as percentages ---
disambig_pct <- disambig_data %>%
  select(lib_id, human_pct, mouse_pct, ambiguous_pct) %>%
  left_join(sample_meta %>% select(lib_id, treatment, cfdna_conc_pg_ul), by = "lib_id")

write_tsv(disambig_pct, "data/disambig_percentages.tsv")

# --- Human fraction vs cfDNA concentration ---
p_human_vs_cfdna <- ggplot(disambig_pct,
  aes(x = cfdna_conc_pg_ul, y = human_pct, color = treatment, label = lib_id)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "cfDNA concentration (pg/µL)", y = "Human-assigned reads (%)",
       title = "Human read fraction vs cfDNA input") +
  theme_white
ggsave("plots/human_pct_vs_cfdna.png", p_human_vs_cfdna, width = 7, height = 5, dpi = 150)

cat("\nScript 02 complete. Key numbers:\n")
cat(sprintf("  Human range: %.1f–%.1f%%\n", min(disambig_pct$human_pct), max(disambig_pct$human_pct)))
cat(sprintf("  Ambiguous range: %.1f–%.1f%%\n", min(disambig_pct$ambiguous_pct), max(disambig_pct$ambiguous_pct)))
