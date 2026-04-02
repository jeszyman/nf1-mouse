## 11_comprehensive_table.R — Build comprehensive metadata table
if (!exists("sample_meta")) source("00_load_data.R")

# Disambiguate
disambig_pct <- read_tsv("data/disambig_percentages.tsv", show_col_types = FALSE) %>%
  select(lib_id, human_pct, mouse_pct, ambiguous_pct)

# Alignment: mapped reads per strategy
align <- read_tsv("data/alignment_metrics.tsv", show_col_types = FALSE) %>%
  select(lib_id, strategy, mapped) %>%
  pivot_wider(names_from = strategy, values_from = mapped, names_prefix = "mapped_")

# TF from rerun (default-low-input preset, both PoN settings)
rerun_path <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun/ichor_rerun_summary.tsv"
ichor <- read_tsv(rerun_path, show_col_types = FALSE) %>%
  filter(preset == "default-low-input") %>%
  mutate(lib_id = str_extract(sample, "lib[0-9]+"),
         col = paste0("tf_", strategy, "_", pon)) %>%
  select(lib_id, col, tumor_fraction) %>%
  pivot_wider(names_from = col, values_from = tumor_fraction)

# Bioanalyzer availability — our 6 samples are NOT in the 25 bioanalyzer runs
# Bioanalyzer runs cover different WU-487 mice (#1-2,#1-3,#2-3,#4-1,#5-1,#5-3,#5-4,#6-5,#7-1,#7-3,#7-4)
# Our mice are #2-1,#3-2,#5-2,#1-4,#3-5,#4-4 — no overlap
bio_status <- tibble(
  lib_id = paste0("lib0", 622:627),
  alt_subj = c("WU-487 2-1", "WU-487 3-2", "WU-487 5-2",
               "WU-487 1-4", "WU-487 3-5", "WU-487 4-4"),
  bioanalyzer = rep("none", 6)
)

# Join
full <- sample_meta %>%
  select(lib_id, mouse_id, treatment, cfdna_conc_pg_ul, tumor_vol_mm3) %>%
  left_join(bio_status, by = "lib_id") %>%
  left_join(disambig_pct, by = "lib_id") %>%
  left_join(align, by = "lib_id") %>%
  left_join(ichor, by = "lib_id")

write_tsv(full, "data/comprehensive_metadata.tsv")

# Convert mapped reads to read pairs for consistency
full <- full %>%
  mutate(mapped_s1 = mapped_s1 / 2,
         mapped_s3a = mapped_s3a / 2,
         mapped_s3b = mapped_s3b / 2) %>%
  rename(mapped_pairs_s1 = mapped_s1,
         mapped_pairs_s3a = mapped_s3a,
         mapped_pairs_s3b = mapped_s3b)

options(width = 200)
cat("\n=== COMPREHENSIVE SAMPLE TABLE ===\n\n")
print(as.data.frame(full))

# Export clean TSV
write_tsv(full, "data/table1.tsv")
cat("\nTable 1 exported: data/table1.tsv\n")
