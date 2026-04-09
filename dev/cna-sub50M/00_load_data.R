## 00_load_data.R — Load and structure all data for sub50M CNA analysis
## Sources: ichorCNA outputs, disambiguate summaries, flagstats, metadata
library(tidyverse)
library(readxl)

# Paths
gcs_base <- "/mnt/gcs/jeszyman/projects/nf1-mouse/emseq/pdx-read-handling-sub50M"
repo     <- "/home/jeszyman/repos/nf1-mouse"
ichor_dir <- file.path(gcs_base, "ichorcna/ichorcna")
bio_dir  <- "/mnt/gcs/jeszyman/projects/nf1-mouse/inputs/bioanalyzer"

# --- Sample metadata ---
libs <- paste0("lib0", 622:627)
sample_meta <- tibble(
  lib_id    = libs,
  sample_id = paste0("smp0", 587:592),
  mouse_id  = paste0("mou000", 1:6),
  pdx_line  = "WU-487",
  chr8q     = "WT",
  treatment = c(rep("vehicle", 3), rep("mirdametinib", 3)),
  bleed_type = "terminal",
  collection_date = c(rep("2024-08-20", 3), rep("2024-08-30", 3)),
  cfdna_conc_pg_ul = c(467, 536, 111, 1720, 2380, 387),
  body_weight_g = c(23, 26, 26, 23, 21, 27),
  tumor_vol_mm3 = c(1823, 1410, 1048, 1902, 1868, 2551)
)

# --- Strategy labels ---
strategies <- c("s1", "s3a", "s3b")
strategy_labels <- c(
  s1  = "Direct (hg38)",
  s3a = "Disambiguate (human only)",
  s3b = "Disambiguate (human + ambig)"
)
strategy_colors <- c(
  s1  = "#e66101",
  s3a = "#5e3c99",
  s3b = "#b2abd2"
)

# --- ichorCNA presets ---
presets <- c("default-low-input", "aggressive", "permissive")

# --- Parse ichorCNA summary ---
ichor_data <- read_tsv(
  file.path(gcs_base, "results_summary/results_summary/ichor_summary.tsv"),
  show_col_types = FALSE
) %>%
  mutate(lib_id = str_extract(sample, "lib\\d+"),
         strategy = str_extract(sample, "s\\d+[ab]?$"))

# --- Parse disambiguate summaries ---
disambig_data <- map_dfr(libs, function(lib) {
  path <- file.path(gcs_base, "results_summary/results_summary",
                    paste0(lib, ".disambig.summary.txt"))
  if (file.exists(path)) {
    read_tsv(path, col_names = c("metric", "value"), show_col_types = FALSE) %>%
      pivot_wider(names_from = metric, values_from = value) %>%
      mutate(lib_id = lib)
  } else {
    tibble(lib_id = lib)
  }
})

# --- Parse flagstats ---
parse_flagstat <- function(path) {
  if (!file.exists(path)) return(tibble())
  lines <- readLines(path)
  vals <- as.numeric(str_extract(lines, "^\\d+"))
  tibble(
    total = vals[1], primary = vals[2], secondary = vals[3],
    supplementary = vals[4], duplicates = vals[5],
    mapped = vals[7], primary_mapped = vals[8],
    properly_paired = vals[10], singletons = vals[12]
  )
}

flagstat_data <- bind_rows(
  # Strategy 1
  map_dfr(libs, function(lib) {
    parse_flagstat(file.path(gcs_base, "strategy1/strategy1", lib,
                             paste0(lib, ".hg38.flagstat.txt"))) %>%
      mutate(lib_id = lib, strategy = "s1")
  }),
  # Strategy 3A
  map_dfr(libs, function(lib) {
    parse_flagstat(file.path(gcs_base, "strategy3/strategy3", lib,
                             paste0(lib, ".disambig.human.flagstat.txt"))) %>%
      mutate(lib_id = lib, strategy = "s3a")
  }),
  # Strategy 3B
  map_dfr(libs, function(lib) {
    parse_flagstat(file.path(gcs_base, "strategy3/strategy3", lib,
                             paste0(lib, ".disambig.human_ambig.flagstat.txt"))) %>%
      mutate(lib_id = lib, strategy = "s3b")
  })
)

# --- Join metadata ---
ichor_full <- ichor_data %>%
  left_join(sample_meta, by = "lib_id")

flagstat_full <- flagstat_data %>%
  left_join(sample_meta, by = "lib_id")

# --- Plot theme ---
theme_white <- theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# --- MPNST CNA loci of interest ---
mpnst_loci <- tibble(
  gene = c("CDKN2A", "CDKN2B", "TP53", "NF1", "SUZ12", "EED", "EZH2"),
  chr  = c("chr9",   "chr9",   "chr17", "chr17", "chr7", "chr11", "chr7"),
  ichor_chr = c("9",  "9",     "17",    "17",    "7",    "11",    "7"),
  start_mb = c(21.9, 22.0, 7.6, 31.1, 159.1, 85.9, 148.8),
  end_mb   = c(22.0, 22.1, 7.7, 31.4, 159.3, 86.1, 149.0),
  role = c("tumor suppressor", "tumor suppressor", "tumor suppressor",
           "tumor suppressor", "PRC2", "PRC2", "PRC2"),
  expected_mpnst = c("loss", "loss", "loss/LOH", "loss", "loss", "loss", "loss")
)

# Chr8q region of interest
chr8q_region <- tibble(
  name = "Chr8q gain region",
  chr = "chr8", ichor_chr = "8",
  start_mb = 50, end_mb = 145,
  note = "Defining feature of Chr8q-gain PDX lines; WU-487 is WT"
)

cat("Data loaded successfully.\n")
cat(sprintf("  %d samples, %d ichorCNA results, %d flagstat entries\n",
            nrow(sample_meta), nrow(ichor_data), nrow(flagstat_data)))
cat(sprintf("  %d disambiguate summaries\n", nrow(disambig_data)))
