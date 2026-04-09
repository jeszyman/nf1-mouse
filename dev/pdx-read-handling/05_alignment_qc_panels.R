library(tidyverse)
library(patchwork)

scale_label <- function(x) {
  ifelse(x >= 1e6, paste0(round(x / 1e6, 1), "M"),
  ifelse(x >= 1e3, paste0(round(x / 1e3, 1), "K"),
  as.character(x)))
}

# --- Parse flagstat files ---
parse_flagstat <- function(path) {
  lines <- readLines(path)
  vals <- as.numeric(str_extract(lines, "^\\d+"))
  tibble(
    total = vals[1], primary = vals[2], secondary = vals[3],
    supplementary = vals[4], primary_mapped = vals[8],
    primary_unmapped = vals[2] - vals[8]
  )
}

qc_dir <- "/tmp/alignment_qc"
libs <- c("lib0622", "lib0623", "lib0624", "lib0625", "lib0626", "lib0627")
refs <- c("hg38", "mm10")

# --- Build flagstat summary ---
flagstat_df <- expand_grid(lib_id = libs, ref = refs) %>%
  mutate(data = map2(lib_id, ref, ~parse_flagstat(
    file.path(qc_dir, paste0(.x, "_", .y, "_flagstat.txt"))
  ))) %>%
  unnest(data)

cat_colors <- c(
  "Primary"       = "#2166AC",
  "Secondary"     = "#FDB863",
  "Supplementary" = "#999999",
  "Unmapped"      = "#B2182B"
)

make_bar_panel <- function(df_row) {
  lib <- df_row$lib_id[1]
  ref <- df_row$ref[1]
  total <- df_row$total[1]

  cat_df <- tibble(
    category = factor(
      c("Primary", "Secondary", "Supplementary", "Unmapped"),
      levels = rev(c("Primary", "Secondary", "Supplementary", "Unmapped"))
    ),
    count = c(df_row$primary_mapped, df_row$secondary, df_row$supplementary, df_row$primary_unmapped)
  ) %>%
    mutate(pct = round(count / total * 100, 1),
           label = paste0(scale_label(count), " (", pct, "%)"))

  ggplot(cat_df, aes(x = count, y = category, fill = category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = label), hjust = -0.05, size = 4) +
    scale_fill_manual(values = cat_colors) +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.4))) +
    labs(x = "Read count", y = NULL, title = paste(lib, "\u2014", ref)) +
    theme_minimal(base_size = 13) +
    theme(axis.text.y = element_text(size = 12), legend.position = "none",
          plot.title = element_text(size = 13, face = "bold"))
}

# --- Build MAPQ panels ---
make_mapq_panel <- function(lib, ref) {
  mapq <- read_tsv(file.path(qc_dir, paste0(lib, "_", ref, "_mapq.txt")),
                   col_names = "mapq", col_types = "i")

  ggplot(mapq, aes(x = mapq)) +
    geom_histogram(binwidth = 5, fill = "#2166AC", color = NA, alpha = 0.8) +
    scale_x_continuous(breaks = seq(0, 60, 10)) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
    labs(x = "MAPQ", y = "Read count", title = paste(lib, "\u2014", ref)) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "bold"))
}

# --- Assemble: rows = libs, cols = hg38 | mm10, panels = bars + mapq ---
bar_panels <- list()
mapq_panels <- list()

for (lib in libs) {
  for (ref in refs) {
    key <- paste(lib, ref, sep = "_")
    row <- flagstat_df %>% filter(lib_id == lib, ref == !!ref)
    bar_panels[[key]] <- make_bar_panel(row)
    mapq_panels[[key]] <- make_mapq_panel(lib, ref)
  }
}

# Layout: 6 rows (libs) x 4 cols (hg38 bar, hg38 mapq, mm10 bar, mm10 mapq)
rows <- list()
for (lib in libs) {
  rows[[lib]] <- bar_panels[[paste(lib, "hg38", sep = "_")]] +
                 mapq_panels[[paste(lib, "hg38", sep = "_")]] +
                 bar_panels[[paste(lib, "mm10", sep = "_")]] +
                 mapq_panels[[paste(lib, "mm10", sep = "_")]]
}

combined <- rows[[1]] / rows[[2]] / rows[[3]] / rows[[4]] / rows[[5]] / rows[[6]] +
  plot_annotation(
    title = "Alignment QC — hg38 vs mm10",
    subtitle = "Left: alignment categories | Right: MAPQ distribution (primary reads)",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 13, color = "gray40")
    )
  )

ggsave("dev/pdx-read-handling/plots/alignment_qc_6lib.png",
       combined, width = 24, height = 30, dpi = 150, bg = "white")

cat("Saved to dev/pdx-read-handling/plots/alignment_qc_6lib.png\n")
