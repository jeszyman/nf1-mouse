## 09_s1_contamination_audit.R — How many mouse reads pass S1's MAPQ filter?
## Key question: reads disambiguate calls "mouse" — what MAPQ do they have in the hg38 alignment?
## If many have MAPQ>=20, they're contamination in Strategy 1.
if (!exists("sample_meta")) source("00_load_data.R")

cat("Extracting MAPQ of disambiguate-mouse reads from hg38 alignments...\n")
cat("(Using 1% subsample for speed)\n\n")

# The Strategy 3 directory has the hg38 BAM (symlinked from S1)
# and the disambiguate mouse BAM. We need to cross-reference:
# get read names from mouse BAM, then check their MAPQ in the hg38 BAM.
#
# Faster approach: the mouse BAM was sorted from the hg38 alignment —
# reads in the mouse BAM retain their hg38 MAPQ. Let's just check the
# MAPQ distribution directly in the mouse BAM.

mouse_mapq <- list()
for (lib in libs) {
  mouse_bam <- file.path(gcs_base, "strategy3/strategy3", lib,
                          paste0(lib, ".disambig.mouse.bam"))
  if (!file.exists(mouse_bam)) {
    cat(sprintf("  %s: mouse BAM not found, skipping\n", lib))
    next
  }

  # 1% subsample MAPQ from mouse BAM (these are reads disambiguate called mouse)
  cmd <- sprintf("samtools view -s 0.01 -F 4 '%s' | cut -f5", mouse_bam)
  mapq_vals <- as.integer(system(cmd, intern = TRUE))
  mouse_mapq[[lib]] <- tibble(
    lib_id = lib,
    mapq = mapq_vals,
    category = "disambig_mouse"
  )
  n_total <- length(mapq_vals)
  n_q20 <- sum(mapq_vals >= 20)
  n_q30 <- sum(mapq_vals >= 30)
  cat(sprintf("  %s: %d mouse reads sampled, %d (%.1f%%) MAPQ>=20, %d (%.1f%%) MAPQ>=30\n",
              lib, n_total, n_q20, 100*n_q20/n_total, n_q30, 100*n_q30/n_total))
}

mouse_mapq_all <- bind_rows(mouse_mapq) %>%
  left_join(sample_meta %>% select(lib_id, treatment, cfdna_conc_pg_ul), by = "lib_id")

write_tsv(mouse_mapq_all, "data/mouse_reads_hg38_mapq.tsv")

# Scale up from 1% subsample to full counts using disambiguate totals
contamination_est <- mouse_mapq_all %>%
  group_by(lib_id) %>%
  summarise(
    sampled_total = n(),
    sampled_q20 = sum(mapq >= 20),
    sampled_q30 = sum(mapq >= 30),
    pct_q20 = 100 * sampled_q20 / sampled_total,
    pct_q30 = 100 * sampled_q30 / sampled_total,
    .groups = "drop"
  ) %>%
  left_join(disambig_data %>% select(lib_id, mouse_assigned), by = "lib_id") %>%
  mutate(
    est_mouse_q20 = round(mouse_assigned * pct_q20 / 100),
    est_mouse_q30 = round(mouse_assigned * pct_q30 / 100)
  ) %>%
  left_join(sample_meta %>% select(lib_id, treatment), by = "lib_id")

write_tsv(contamination_est, "data/s1_contamination_estimate.tsv")

cat("\nEstimated mouse contamination in S1 (reads passing MAPQ filter):\n")
print(as.data.frame(contamination_est %>%
  select(lib_id, treatment, mouse_assigned, pct_q20, est_mouse_q20, pct_q30, est_mouse_q30)))

# Plot: MAPQ distribution of disambiguate-mouse reads
p_mouse_mapq <- ggplot(mouse_mapq_all, aes(x = mapq)) +
  geom_histogram(binwidth = 1, fill = "#e66101", alpha = 0.7) +
  geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "darkred") +
  facet_wrap(~lib_id, scales = "free_y") +
  labs(x = "MAPQ (in hg38 alignment)", y = "Count (1% subsample)",
       title = "MAPQ of disambiguate-mouse reads in hg38 alignment",
       subtitle = "Reads right of dashed lines are mouse contamination in S1") +
  theme_white
ggsave("plots/mouse_reads_hg38_mapq.png", p_mouse_mapq, width = 14, height = 8, dpi = 150)

cat("\nScript 09 complete.\n")
