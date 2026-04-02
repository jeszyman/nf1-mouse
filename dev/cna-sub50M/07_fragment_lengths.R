## 07_fragment_lengths.R — Fragment length profiles from BAMs
if (!exists("sample_meta")) source("00_load_data.R")

# --- Extract insert size distributions from BAMs using samtools stats ---
extract_insert_sizes <- function(bam_path, label, lib) {
  # Subsample ~1M reads for speed (especially over GCS mount)
  cmd <- sprintf("samtools view -s 0.01 -b '%s' | samtools stats | grep '^IS'", bam_path)
  raw <- system(cmd, intern = TRUE)
  if (length(raw) == 0) return(tibble())
  read_tsv(I(raw), col_names = FALSE, show_col_types = FALSE) %>%
    select(insert_size = X2, pairs_total = X3, inward = X4, outward = X5, other = X6) %>%
    mutate(strategy = label, lib_id = lib)
}

cat("Extracting insert size distributions from BAMs...\n")
frag_data <- list()
for (lib in libs) {
  s1_bam <- file.path(gcs_base, "strategy1/strategy1", lib,
                       paste0(lib, ".hg38.sorted.bam"))
  s3a_bam <- file.path(gcs_base, "strategy3/strategy3", lib,
                        paste0(lib, ".disambig.human.bam"))
  s3b_bam <- file.path(gcs_base, "strategy3/strategy3", lib,
                        paste0(lib, ".disambig.human_ambig.bam"))

  if (file.exists(s1_bam))  frag_data[[paste0(lib, "_s1")]]  <- extract_insert_sizes(s1_bam, "s1", lib)
  if (file.exists(s3a_bam)) frag_data[[paste0(lib, "_s3a")]] <- extract_insert_sizes(s3a_bam, "s3a", lib)
  if (file.exists(s3b_bam)) frag_data[[paste0(lib, "_s3b")]] <- extract_insert_sizes(s3b_bam, "s3b", lib)
  cat(sprintf("  %s done\n", lib))
}

frag_all <- bind_rows(frag_data) %>%
  filter(insert_size >= 50, insert_size <= 500) %>%
  left_join(sample_meta %>% select(lib_id, treatment), by = "lib_id")

write_tsv(frag_all, "data/fragment_lengths.tsv")

# --- Normalize to density per sample×strategy ---
frag_norm <- frag_all %>%
  group_by(lib_id, strategy) %>%
  mutate(density = pairs_total / sum(pairs_total)) %>%
  ungroup()

# --- Overlay by strategy per sample ---
p_frag <- ggplot(frag_norm, aes(x = insert_size, y = density, color = strategy)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = strategy_colors, labels = strategy_labels) +
  facet_wrap(~lib_id, scales = "free_y") +
  labs(x = "Insert size (bp)", y = "Density",
       title = "Fragment length distribution by strategy") +
  theme_white
ggsave("plots/fragment_lengths.png", p_frag, width = 14, height = 8, dpi = 150)

# --- Short/long ratio (100-150 / 151-220) ---
sl_ratio <- frag_all %>%
  mutate(bin = case_when(
    insert_size >= 100 & insert_size <= 150 ~ "short",
    insert_size >= 151 & insert_size <= 220 ~ "long",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(bin)) %>%
  group_by(lib_id, strategy, bin) %>%
  summarise(total = sum(pairs_total), .groups = "drop") %>%
  pivot_wider(names_from = bin, values_from = total) %>%
  mutate(sl_ratio = short / long) %>%
  left_join(sample_meta %>% select(lib_id, treatment), by = "lib_id")

write_tsv(sl_ratio, "data/short_long_ratio.tsv")

p_sl <- ggplot(sl_ratio, aes(x = lib_id, y = sl_ratio, fill = strategy)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = strategy_colors, labels = strategy_labels) +
  labs(x = NULL, y = "Short/Long ratio (100-150 / 151-220 bp)",
       title = "Fragment short/long ratio by strategy") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/short_long_ratio.png", p_sl, width = 10, height = 5, dpi = 150)

# --- Nucleosome peak position ---
peak_data <- frag_norm %>%
  filter(insert_size >= 100, insert_size <= 250) %>%
  group_by(lib_id, strategy) %>%
  slice_max(density, n = 1) %>%
  ungroup()

cat("\nNucleosome peak positions:\n")
print(peak_data %>% select(lib_id, strategy, insert_size, density))

cat("\nScript 07 complete.\n")
