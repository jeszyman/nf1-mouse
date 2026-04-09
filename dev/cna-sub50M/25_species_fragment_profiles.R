## 25_species_fragment_profiles.R — Fragment length by species classification
## Shows insert size distributions for human, mouse, and ambiguous reads separately
library(tidyverse)
library(patchwork)

theme_white <- theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Extract insert sizes from each disambiguate category BAM (lib0626 local)
bam_dir <- "/mnt/data/projects/nf1-mouse/emseq/cna-sub50M-rerun/bams/lib0626"

extract_is <- function(bam_path, label) {
  cmd <- sprintf("conda run -n biotools samtools view -f 2 -s 0.005 '%s' | awk '{s=($9<0?-$9:$9); if(s>=50 && s<=500) print s}'", bam_path)
  sizes <- as.integer(system(cmd, intern = TRUE))
  tibble(insert_size = sizes, category = label)
}

cat("Extracting insert sizes by species category (lib0626, 0.5% subsample)...\n")
d_human <- extract_is(file.path(bam_dir, "lib0626.disambig.human.bam"), "Human")
cat(sprintf("  Human: %d fragments\n", nrow(d_human)))
d_mouse <- extract_is(file.path(bam_dir, "lib0626.disambig.mouse.bam"), "Mouse")
cat(sprintf("  Mouse: %d fragments\n", nrow(d_mouse)))
d_ambig <- extract_is(file.path(bam_dir, "lib0626.disambig.ambiguous.bam"), "Ambiguous")
cat(sprintf("  Ambiguous: %d fragments\n", nrow(d_ambig)))
d_s1 <- extract_is(file.path(bam_dir, "lib0626.hg38.sorted.bam"), "S1 (all)")
cat(sprintf("  S1: %d fragments\n", nrow(d_s1)))

frag_species <- bind_rows(d_human, d_mouse, d_ambig, d_s1)
write_tsv(frag_species, "data/species_fragment_profiles.tsv")

# Normalize each category to its own total
frag_pct <- frag_species %>%
  group_by(category) %>%
  count(insert_size) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = c("S1 (all)", "Mouse", "Ambiguous", "Human")))

# Panel A: Overlay all categories
p_overlay <- ggplot(frag_pct, aes(x = insert_size, y = pct, color = category)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  scale_color_manual(values = c("S1 (all)" = "#e66101", "Mouse" = "#d95f02",
                                 "Ambiguous" = "#b2abd2", "Human" = "#5e3c99")) +
  geom_vline(xintercept = 167, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 145, linetype = "dotted", color = "#5e3c99", linewidth = 0.5) +
  annotate("text", x = 170, y = max(frag_pct$pct) * 0.95, label = "167 bp", size = 4, hjust = 0) +
  annotate("text", x = 135, y = max(frag_pct$pct) * 0.85, label = "145 bp", size = 4, hjust = 1, color = "#5e3c99") +
  labs(x = "Insert size (bp)", y = "% of total per category",
       title = "lib0626: Fragment length by species classification",
       subtitle = sprintf("Human: %s fragments | Mouse: %s | Ambiguous: %s | S1: %s (0.5%% subsample)",
                          format(nrow(d_human), big.mark=","),
                          format(nrow(d_mouse), big.mark=","),
                          format(nrow(d_ambig), big.mark=","),
                          format(nrow(d_s1), big.mark=","))) +
  theme_white +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))

# Panel B: Faceted for clearer individual shapes
p_faceted <- ggplot(frag_pct, aes(x = insert_size, y = pct, fill = category)) +
  geom_area(alpha = 0.6) +
  geom_vline(xintercept = 167, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 145, linetype = "dotted", color = "#5e3c99") +
  scale_fill_manual(values = c("S1 (all)" = "#e66101", "Mouse" = "#d95f02",
                                "Ambiguous" = "#b2abd2", "Human" = "#5e3c99")) +
  facet_wrap(~category, ncol = 1, scales = "free_y") +
  labs(x = "Insert size (bp)", y = "% of total",
       title = "Fragment length profiles by disambiguate category (lib0626)") +
  guides(fill = "none") +
  theme_white +
  theme(strip.text = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13))

combined <- p_overlay / p_faceted + plot_layout(heights = c(1, 1.5))
ggsave("plots/species_fragment_profiles.png", combined, width = 12, height = 14, dpi = 150)

# Summary stats
summary_stats <- frag_species %>%
  group_by(category) %>%
  summarise(
    n = n(),
    median = median(insert_size),
    peak = {
      h <- hist(insert_size, breaks = seq(50, 500, 1), plot = FALSE)
      h$mids[which.max(h$counts)]
    },
    pct_sub150 = 100 * mean(insert_size < 150),
    pct_150_180 = 100 * mean(insert_size >= 150 & insert_size <= 180),
    pct_167 = 100 * mean(insert_size >= 165 & insert_size <= 170),
    .groups = "drop"
  )

cat("\nFragment length summary by species:\n")
print(as.data.frame(summary_stats))

cat("\nScript 25 complete.\n")
