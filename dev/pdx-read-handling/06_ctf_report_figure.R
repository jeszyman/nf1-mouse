library(tidyverse)
library(patchwork)

theme_report <- theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size = 14, face = "bold"))

# --- Data ---
disambig  <- read_tsv("dev/cna-sub50M/data/disambig_percentages.tsv", show_col_types = FALSE)
summary_d <- read_tsv("dev/cna-sub50M/data/composite_summary.tsv", show_col_types = FALSE)
segments  <- read_tsv("dev/cna-sub50M/data/cna_segments.tsv", show_col_types = FALSE)

lib_order <- summary_d %>% arrange(tf_3a) %>% pull(lib_id)

# --- Panel A: Species classification stacked bar ---
disambig_long <- disambig %>%
  select(lib_id, human_pct, mouse_pct, ambiguous_pct) %>%
  pivot_longer(-lib_id, names_to = "species", values_to = "pct") %>%
  mutate(
    species = factor(
      recode(species, human_pct = "Human", mouse_pct = "Mouse", ambiguous_pct = "Ambiguous"),
      levels = c("Mouse", "Ambiguous", "Human")
    ),
    lib_id = factor(lib_id, levels = lib_order)
  )

pA <- ggplot(disambig_long, aes(x = lib_id, y = pct, fill = species)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(Mouse = "#5e3c99", Ambiguous = "#b2abd2", Human = "#e66101")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Read pairs (%)", fill = "Species\nassignment",
       title = "A. Species classification") +
  theme_report +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

# --- Panel B: TF vs tumor volume (no treatment coloring) ---
pB <- ggplot(summary_d, aes(x = tumor_vol_mm3, y = tf_3a)) +
  geom_smooth(method = "lm", se = TRUE, color = "gray50", fill = "gray85", linewidth = 0.8) +
  geom_point(size = 4, color = "#2166AC") +
  geom_text(aes(label = lib_id), nudge_y = 0.03, size = 3.5) +
  annotate("text", x = min(summary_d$tumor_vol_mm3) + 50, y = 0.6,
           label = "rho == 0.54~~(n==6)", parse = TRUE, size = 4.5, hjust = 0) +
  labs(x = expression("Terminal tumor volume (mm"^3*")"),
       y = "ichorCNA tumor fraction",
       title = "B. Tumor fraction vs tumor volume") +
  theme_report

# --- Panel C: Genome-wide CNA (lib0627) — improved version ---
loci_bands <- tibble(
  label = c("7p", "8q", "9q", "17q"),
  chr_num = c(7, 8, 9, 17),
  start_mb = c(1, 50, 40, 25),
  end_mb   = c(75, 145, 141, 83)
)

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

event_colors <- c(GAIN = "#D7191C", HETD = "#2C7BB6", NEUT = "gray70")

pC <- ggplot() +
  geom_rect(data = loci_bands,
            aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf),
            fill = "#B8E186", alpha = 0.3, inherit.aes = FALSE) +
  geom_point(data = seg_627,
             aes(x = (start + end) / 2e6, y = logR, color = event),
             size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_point(data = gene_marks,
             aes(x = pos_mb, y = -2.2),
             shape = 23, size = 3, fill = "#4DAF4A", color = "#4DAF4A") +
  geom_text(data = gene_marks,
            aes(x = pos_mb, y = -2.5, label = gene),
            size = 3.5, fontface = "italic", angle = 30, hjust = 1, vjust = 1) +
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
  labs(x = "Chromosome", y = "Copy Number (log2 of ratio)", color = NULL,
       title = "C. Genome-wide CNA from cfDNA (lib0627, TF = 65%)") +
  theme_minimal(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.15, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 10, 30, 10)
  )

# --- Combine: A+B on top, C spanning bottom ---
combined <- (pA | pB) / pC +
  plot_layout(heights = c(1, 1.2))

ggsave("dev/pdx-read-handling/plots/ctf_progress_figure.png",
       combined, width = 16, height = 12, dpi = 200, bg = "white")

cat("Saved to dev/pdx-read-handling/plots/ctf_progress_figure.png\n")
