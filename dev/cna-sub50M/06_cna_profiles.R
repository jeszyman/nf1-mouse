## 06_cna_profiles.R — CNA profiles at genome-wide and MPNST-specific loci
if (!exists("sample_meta")) source("00_load_data.R")

rerun_base <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun"

# --- Load segment files from re-run (nopon, default-low-input as primary) ---
load_seg <- function(base_dir, pon_mode, preset) {
  seg_files <- list.files(
    file.path(base_dir, pon_mode),
    pattern = "\\.cna\\.seg$", recursive = TRUE, full.names = TRUE
  )
  # Filter to preset
  seg_files <- seg_files[grepl(preset, seg_files)]

  map_dfr(seg_files, function(f) {
    d <- read_tsv(f, show_col_types = FALSE)
    sid <- str_extract(basename(f), "lib\\d+_s\\d+[ab]?")
    logr_col <- grep("logR$", names(d), value = TRUE)[1]
    cn_col <- grep("copy.number$", names(d), value = TRUE)[1]
    event_col <- grep("event$", names(d), value = TRUE)[1]
    d %>%
      rename(logR = all_of(logr_col), copy_number = all_of(cn_col),
             event = all_of(event_col)) %>%
      select(chr, start, end, logR, copy_number, event) %>%
      mutate(
        sample_id = sid,
        lib_id = str_extract(sid, "lib\\d+"),
        strategy = str_extract(sid, "s\\d+[ab]?$"),
        pon = pon_mode,
        preset = preset
      )
  })
}

seg_nopon <- load_seg(rerun_base, "nopon", "default-low-input")
seg_withpon <- load_seg(rerun_base, "withpon", "default-low-input")
seg_all <- bind_rows(seg_nopon, seg_withpon)

if (nrow(seg_all) == 0) {
  cat("No segment files found — ichorCNA re-run may not be complete.\n")
  quit(save = "no")
}

write_tsv(seg_all, "data/cna_segments.tsv")

# --- Genome-wide CNA plot per sample (overlay strategies) ---
for (lib in libs) {
  seg_lib <- seg_all %>%
    filter(lib_id == lib, pon == "nopon") %>%
    mutate(chr_num = as.numeric(str_replace(chr, "chr", ""))) %>%
    filter(!is.na(chr_num))

  if (nrow(seg_lib) == 0) next

  p <- ggplot(seg_lib, aes(x = start / 1e6, y = logR, color = strategy)) +
    geom_segment(aes(xend = end / 1e6, yend = logR), linewidth = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = strategy_colors, labels = strategy_labels) +
    facet_grid(. ~ chr_num, scales = "free_x", space = "free_x") +
    labs(x = "Position (Mb)", y = "log2 ratio", color = "Strategy",
         title = sprintf("CNA profile: %s (%s)", lib,
                         sample_meta$treatment[sample_meta$lib_id == lib])) +
    theme_white +
    theme(axis.text.x = element_blank(), strip.text = element_text(size = 7),
          panel.spacing = unit(0.1, "lines"))
  ggsave(sprintf("plots/cna_genomewide_%s.png", lib), p,
         width = 18, height = 4, dpi = 150)
}

# --- MPNST loci zoom ---
# Extract CNA values at key loci for all samples/strategies
locus_cna <- map_dfr(1:nrow(mpnst_loci), function(i) {
  locus <- mpnst_loci[i, ]
  seg_all %>%
    filter(chr %in% c(locus$chr, locus$ichor_chr),
           pon == "nopon") %>%
    filter(start / 1e6 <= locus$end_mb, end / 1e6 >= locus$start_mb) %>%
    mutate(gene = locus$gene, expected = locus$expected_mpnst)
})

if (nrow(locus_cna) > 0) {
  p_loci <- ggplot(locus_cna, aes(x = gene, y = logR, color = strategy)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = strategy_colors, labels = strategy_labels) +
    facet_wrap(~lib_id) +
    labs(x = NULL, y = "log2 ratio at locus", color = "Strategy",
         title = "CNA at MPNST-relevant loci (all samples, nopon default-low-input)") +
    theme_white +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("plots/cna_mpnst_loci.png", p_loci, width = 14, height = 8, dpi = 150)
  write_tsv(locus_cna %>% select(lib_id, strategy, gene, chr, start, end, logR),
            "data/cna_mpnst_loci.tsv")
}

# --- Chr8q check (should be neutral for WU-487 Chr8q-WT) ---
chr8q_cna <- seg_all %>%
  filter(chr == "chr8" | chr == "8", pon == "nopon") %>%
  filter(start / 1e6 >= 50, end / 1e6 <= 145) %>%
  mutate(region = "Chr8q (50-145 Mb)")

if (nrow(chr8q_cna) > 0) {
  p_chr8q <- ggplot(chr8q_cna, aes(x = start / 1e6, y = logR, color = strategy)) +
    geom_segment(aes(xend = end / 1e6, yend = logR), linewidth = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = strategy_colors, labels = strategy_labels) +
    facet_wrap(~lib_id) +
    labs(x = "Chr8 position (Mb)", y = "log2 ratio",
         title = "Chr8q region (WU-487 is WT — should be neutral)") +
    theme_white
  ggsave("plots/cna_chr8q.png", p_chr8q, width = 14, height = 8, dpi = 150)
}

cat("\nScript 06 complete.\n")
