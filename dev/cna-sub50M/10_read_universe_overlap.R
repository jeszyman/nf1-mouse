## 10_read_universe_overlap.R — Read universe overlap: S1 vs disambiguate categories
## Venn-style analysis with Jaccard indices showing where every read ends up
## across Strategy 1 (direct hg38) and Strategy 3 (disambiguate)
if (!exists("sample_meta")) source("00_load_data.R")

cat("Extracting read names from all BAMs (1% subsample)...\n")

get_readnames <- function(bam_path, min_mapq = 0) {
  if (!file.exists(bam_path)) return(character(0))
  cmd <- sprintf("samtools view -s 0.01 -F 4 %s '%s' | cut -f1 | sort -u",
                 if (min_mapq > 0) sprintf("-q %d", min_mapq) else "", bam_path)
  system(cmd, intern = TRUE)
}

results <- list()

for (lib in libs) {
  cat(sprintf("  %s: extracting...\n", lib))

  s1_bam <- file.path(gcs_base, "strategy1/strategy1", lib,
                       paste0(lib, ".hg38.sorted.bam"))
  human_bam <- file.path(gcs_base, "strategy3/strategy3", lib,
                          paste0(lib, ".disambig.human.bam"))
  mouse_bam <- file.path(gcs_base, "strategy3/strategy3", lib,
                          paste0(lib, ".disambig.mouse.bam"))
  ambig_bam <- file.path(gcs_base, "strategy3/strategy3", lib,
                          paste0(lib, ".disambig.ambiguous.bam"))

  # S1: all mapped reads (any MAPQ) and MAPQ>=20 subset
  s1_all <- get_readnames(s1_bam, 0)
  s1_q20 <- get_readnames(s1_bam, 20)

  # Disambiguate categories
  d_human <- get_readnames(human_bam, 0)
  d_mouse <- get_readnames(mouse_bam, 0)
  d_ambig <- get_readnames(ambig_bam, 0)

  # Compute overlaps
  s1q20_in_human <- length(intersect(s1_q20, d_human))
  s1q20_in_mouse <- length(intersect(s1_q20, d_mouse))
  s1q20_in_ambig <- length(intersect(s1_q20, d_ambig))
  s1q20_only     <- length(setdiff(s1_q20, c(d_human, d_mouse, d_ambig)))

  human_not_in_s1q20 <- length(setdiff(d_human, s1_q20))
  ambig_not_in_s1q20 <- length(setdiff(d_ambig, s1_q20))

  # Jaccard: S1_q20 vs disambig_human
  jaccard_s1_human <- length(intersect(s1_q20, d_human)) /
    length(union(s1_q20, d_human))

  # Jaccard: S1_q20 vs disambig_human+ambig (S3B)
  s3b_reads <- union(d_human, d_ambig)
  jaccard_s1_s3b <- length(intersect(s1_q20, s3b_reads)) /
    length(union(s1_q20, s3b_reads))

  results[[lib]] <- tibble(
    lib_id = lib,
    s1_q20_total = length(s1_q20),
    s1q20_in_human = s1q20_in_human,
    s1q20_in_mouse = s1q20_in_mouse,
    s1q20_in_ambig = s1q20_in_ambig,
    s1q20_unclassified = s1q20_only,
    human_not_in_s1q20 = human_not_in_s1q20,
    ambig_not_in_s1q20 = ambig_not_in_s1q20,
    pct_s1q20_is_mouse = 100 * s1q20_in_mouse / length(s1_q20),
    pct_s1q20_is_ambig = 100 * s1q20_in_ambig / length(s1_q20),
    jaccard_s1_human = jaccard_s1_human,
    jaccard_s1_s3b = jaccard_s1_s3b
  )

  cat(sprintf("    S1 MAPQ>=20: %d reads (1%% sample)\n", length(s1_q20)))
  cat(sprintf("    → in disambig human: %d (%.1f%%)\n", s1q20_in_human,
              100 * s1q20_in_human / length(s1_q20)))
  cat(sprintf("    → in disambig MOUSE: %d (%.1f%%) ← contamination\n", s1q20_in_mouse,
              100 * s1q20_in_mouse / length(s1_q20)))
  cat(sprintf("    → in disambig ambig: %d (%.1f%%)\n", s1q20_in_ambig,
              100 * s1q20_in_ambig / length(s1_q20)))
  cat(sprintf("    Jaccard(S1, S3A): %.3f | Jaccard(S1, S3B): %.3f\n",
              jaccard_s1_human, jaccard_s1_s3b))
}

overlap_data <- bind_rows(results) %>%
  left_join(sample_meta %>% select(lib_id, treatment, cfdna_conc_pg_ul), by = "lib_id")

write_tsv(overlap_data, "data/read_universe_overlap.tsv")

# --- Stacked bar: where do S1 MAPQ>=20 reads go in disambiguate? ---
overlap_long <- overlap_data %>%
  select(lib_id, s1q20_in_human, s1q20_in_mouse, s1q20_in_ambig, s1q20_unclassified) %>%
  pivot_longer(-lib_id, names_to = "category", values_to = "count") %>%
  mutate(category = factor(category,
    levels = c("s1q20_in_human", "s1q20_in_ambig", "s1q20_in_mouse", "s1q20_unclassified"),
    labels = c("Human (confirmed)", "Ambiguous", "MOUSE (contamination)", "Unclassified")))

p_contam <- ggplot(overlap_long, aes(x = lib_id, y = count, fill = category)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("Human (confirmed)" = "#5e3c99",
                                "Ambiguous" = "#b2abd2",
                                "MOUSE (contamination)" = "#e66101",
                                "Unclassified" = "gray70")) +
  labs(x = NULL, y = "Read count (1% subsample)", fill = NULL,
       title = "Where do S1 MAPQ≥20 reads actually come from?",
       subtitle = "Orange = mouse reads that pass S1's quality filter (contamination)") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/s1_read_classification.png", p_contam, width = 10, height = 6, dpi = 150)

# --- Jaccard indices ---
p_jaccard <- ggplot(overlap_data, aes(x = lib_id)) +
  geom_col(aes(y = jaccard_s1_s3b), fill = "#b2abd2", alpha = 0.8, width = 0.4,
           position = position_nudge(x = -0.2)) +
  geom_col(aes(y = jaccard_s1_human), fill = "#5e3c99", alpha = 0.8, width = 0.4,
           position = position_nudge(x = 0.2)) +
  labs(x = NULL, y = "Jaccard index",
       title = "Read set similarity: S1 (MAPQ≥20) vs disambiguate",
       subtitle = "Purple = vs S3A (human only), Light = vs S3B (human + ambig)") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/jaccard_indices.png", p_jaccard, width = 8, height = 5, dpi = 150)

cat("\nScript 10 complete.\n")
