library(tidyverse)

segments <- read_tsv("dev/cna-sub50M/data/cna_segments.tsv", show_col_types = FALSE)

# MPNST-associated loci — green highlight bands + gene diamonds
# From Szymanski 2021 Fig 3: 1q, 7p, 8q, 9p/9q, 17q
loci_bands <- tibble(
  label = c("7p", "8q", "9q", "17q"),
  chr_num = c(7, 8, 9, 17),
  start_mb = c(1, 50, 40, 25),
  end_mb   = c(75, 145, 141, 83)
)

# Only label genes that overlap actual CNA events in this sample
gene_marks <- tibble(
  gene = c("SMARCA2", "CDKN2A/B", "NF1"),
  chr_num = c(9, 9, 17),
  pos_mb = c(2, 22, 31.2)
)

seg_627 <- segments %>%
  filter(lib_id == "lib0627", strategy == "s3a", pon == "nopon",
         preset == "default-low-input") %>%
  mutate(chr_num = as.numeric(chr)) %>%
  filter(!is.na(chr_num), chr_num <= 22)

# Point colors by event
event_colors <- c(GAIN = "#D7191C", HETD = "#2C7BB6", NEUT = "gray70")

p <- ggplot() +
  # Green highlight bands for MPNST loci
  geom_rect(data = loci_bands,
            aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
            fill = "#B8E186", alpha = 0.3, inherit.aes = FALSE) +
  # CNA segments as points
  geom_point(data = seg_627,
             aes(x = (start + end) / 2e6, y = logR, color = event),
             size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  # Gene diamond markers
  geom_point(data = gene_marks,
             aes(x = pos_mb, y = -2.2),
             shape = 23, size = 3, fill = "#4DAF4A", color = "#4DAF4A") +
  geom_text(data = gene_marks,
            aes(x = pos_mb, y = -2.5, label = gene),
            size = 3.5, fontface = "italic", angle = 30, hjust = 1, vjust = 1) +
  # Locus band labels at top
  geom_text(data = loci_bands,
            aes(x = (start_mb + end_mb) / 2, y = 2.1, label = label),
            size = 4, fontface = "bold", color = "gray30") +
  scale_color_manual(values = event_colors,
                     labels = c(GAIN = "Copy number gain",
                                HETD = "Copy number loss",
                                NEUT = "Copy number neutral")) +
  facet_grid(. ~ chr_num, scales = "free_x", space = "free_x",
             labeller = labeller(chr_num = function(x) x)) +
  coord_cartesian(ylim = c(-2.8, 2.3), clip = "off") +
  scale_y_continuous(breaks = c(-2, -1, 0, 0.57, 1, 2),
                     labels = c("-2.00", "-1.00", "0.00", "0.57", "1.00", "2.00")) +
  labs(x = "Chromosome", y = "Copy Number (log2 of ratio)",
       color = NULL,
       title = "Genome-wide CNA from cfDNA — lib0627 (WU-487, TF = 65%)") +
  theme_minimal(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0.15, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    plot.margin = margin(10, 10, 30, 10)
  )

ggsave("dev/pdx-read-handling/plots/ctf_panel_c.png",
       p, width = 18, height = 6, dpi = 200, bg = "white")

cat("Saved to dev/pdx-read-handling/plots/ctf_panel_c.png\n")
