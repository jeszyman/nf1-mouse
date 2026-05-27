#!/usr/bin/env Rscript
# 04_chr8_plots.R -- Phase 5: chr8 zoomed log2 ratio plots + summary TSV.
#
# Reads per-bin logR and segment calls from ichorCNA *.cna.seg, plots
# chr8 with the hg38 centromere band annotated and the p/q arms
# labeled, and writes a per-sample summary TSV (TF, ploidy, median
# logR for chr8p/chr8q, q-minus-p delta, mapped reads, dup rate).
#
# Pure base R -- system R 4.6.0 has a SETLENGTH ABI break that prevents
# loading ggplot2/data.table/rlang.so from /usr/local/lib/R/site-library;
# the R_libs/ tree has only a patched data.table for ichorCNA itself.
# Base R is sufficient for these per-chromosome plots.
#
# Usage (from dev/chr8-amp-detect/):
#   sudo -n -u jupyter /usr/bin/Rscript 04_chr8_plots.R            # nopon
#   sudo -n -u jupyter /usr/bin/Rscript 04_chr8_plots.R hdulp-pon
#   sudo -n -u jupyter /usr/bin/Rscript 04_chr8_plots.R both
#
# Writes to /mnt/data/projects/nf1-mouse/chr8-amp-detect/plots/.

args <- commandArgs(trailingOnly = TRUE)
arg <- if (length(args) >= 1) args[1] else "nopon"
configs <- if (arg == "both") c("nopon", "hdulp-pon") else arg
stopifnot(all(configs %in% c("nopon", "hdulp-pon")))

BASE_DIR   <- "/mnt/data/projects/nf1-mouse/chr8-amp-detect"
ICHOR_ROOT <- file.path(BASE_DIR, "ichorcna")
BAM_DIR    <- file.path(BASE_DIR, "bams")
PLOTS_DIR  <- file.path(BASE_DIR, "plots")
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

## hg38 chr8 coordinates (UCSC gap.txt centromere acen)
CHR8_LEN     <- 145138636
CENTRO_START <-  44033744
CENTRO_END   <-  45877265
CENTRO_MID   <- (CENTRO_START + CENTRO_END) / 2

## Panel layout: AMP top row, WT bottom row, left-to-right by library
samples <- data.frame(
  library_id   = c("lib0644", "lib0645", "lib0626", "lib0627"),
  sample_id    = c("smp0927", "smp0928", "smp0626", "smp0627"),
  pdx_line     = c("MN2",     "MN2",     "WU-487",  "WU-487"),
  chr8q_status = c("AMP",     "AMP",     "WT",      "WT"),
  stringsAsFactors = FALSE
)

## ichorCNA event labels seen in cna.seg / seg.txt
event_colors <- c(
  HOMD  = "#1f4f8b",
  HETD  = "#2ca02c",
  NEUT  = "grey55",
  GAIN  = "#d62728",
  AMP   = "#d62728",
  HLAMP = "#7f0000"
)
default_color <- "grey55"

ichor_subdir <- function(lib, config) {
  file.path(ICHOR_ROOT, sprintf("%s_%s_default-low-input", lib, config))
}
ichor_id <- function(lib, config) sprintf("%s_%s_dli", lib, config)

read_params <- function(lib, config) {
  f <- file.path(ichor_subdir(lib, config),
                 paste0(ichor_id(lib, config), ".params.txt"))
  txt <- readLines(f)
  num <- function(prefix) {
    ln <- grep(paste0("^", prefix, ":"), txt, value = TRUE)[1]
    as.numeric(sub(paste0("^", prefix, ":\\s*"), "", ln))
  }
  list(tf = num("Tumor Fraction"), ploidy = num("Ploidy"))
}

read_chr8_bins <- function(lib, config) {
  f <- file.path(ichor_subdir(lib, config),
                 paste0(ichor_id(lib, config), ".cna.seg"))
  d <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[d$chr == "chr8", ]
  logr  <- grep("\\.logR$", names(d), value = TRUE)[1]
  evt   <- grep("\\.event$", names(d), value = TRUE)[1]
  cn    <- grep("\\.copy\\.number$", names(d), value = TRUE)[1]
  data.frame(
    start   = d$start,
    end     = d$end,
    mid     = (d$start + d$end) / 2,
    logR    = suppressWarnings(as.numeric(d[[logr]])),
    event   = d[[evt]],
    cn      = as.integer(d[[cn]]),
    stringsAsFactors = FALSE
  )
}

read_chr8_segs <- function(lib, config) {
  f <- file.path(ichor_subdir(lib, config),
                 paste0(ichor_id(lib, config), ".seg.txt"))
  d <- read.delim(f, stringsAsFactors = FALSE)
  d[d$chrom == "chr8", ]
}

read_flagstat <- function(lib) {
  fd <- file.path(BAM_DIR, sprintf("%s.dedup.flagstat.txt",  lib))
  fs <- file.path(BAM_DIR, sprintf("%s.sorted.flagstat.txt", lib))
  if (!file.exists(fd) || !file.exists(fs))
    return(list(mapped = NA_integer_, dup_pct = NA_real_))
  primary_mapped <- function(p) {
    ln <- grep("primary mapped", readLines(p), value = TRUE)[1]
    as.numeric(sub("^([0-9]+).*", "\\1", ln))
  }
  pre  <- primary_mapped(fs)
  post <- primary_mapped(fd)
  list(mapped = post, dup_pct = 100 * (1 - post / pre))
}

plot_chr8_panel <- function(lib, line, status, bins, segs, parm) {
  rng <- range(c(bins$logR, segs$seg.median.logR, 0), na.rm = TRUE)
  ylim_top <- max(2, ceiling(rng[2]))
  ylim_bot <- min(-2, floor(rng[1]))
  cols <- event_colors[bins$event]
  cols[is.na(cols)] <- default_color
  title <- sprintf("%s  %s (chr8q %s)   TF=%.3f  ploidy=%.2f",
                   lib, line, status, parm$tf, parm$ploidy)
  plot(bins$mid / 1e6, bins$logR,
       xlim = c(0, CHR8_LEN / 1e6),
       ylim = c(ylim_bot, ylim_top),
       pch = 20, cex = 0.55, col = cols,
       xlab = "chr8 position (Mb)",
       ylab = "log2 ratio",
       main = title, cex.main = 0.9, las = 1)
  abline(h = 0, col = "grey85")
  ## centromere band
  rect(CENTRO_START / 1e6, ylim_bot, CENTRO_END / 1e6, ylim_top,
       col = adjustcolor("grey70", alpha.f = 0.30), border = NA)
  abline(v = CENTRO_MID / 1e6, col = "grey40", lty = 2)
  ## arm labels in upper margin
  axis(3, at = CENTRO_START / 2e6, labels = "p",
       tick = FALSE, line = -1.0, cex.axis = 0.85, col.axis = "grey30")
  axis(3, at = (CENTRO_END + CHR8_LEN) / 2e6, labels = "q",
       tick = FALSE, line = -1.0, cex.axis = 0.85, col.axis = "grey30")
  ## segment median lines
  if (nrow(segs) > 0) {
    seg_cols <- event_colors[segs$call]
    seg_cols[is.na(seg_cols)] <- "black"
    segments(x0 = segs$start / 1e6, x1 = segs$end / 1e6,
             y0 = segs$seg.median.logR, y1 = segs$seg.median.logR,
             col = seg_cols, lwd = 3)
  }
}

run_one_config <- function(config) {
  png_path <- file.path(PLOTS_DIR, sprintf("chr8_zoomed_%s.png", config))
  png(png_path, width = 1700, height = 1150, res = 150)
  par(mfrow = c(2, 2), mar = c(4.2, 4.3, 2.6, 1.0), oma = c(0, 0, 2.4, 0))

  summary_rows <- vector("list", nrow(samples))
  for (i in seq_len(nrow(samples))) {
    lib  <- samples$library_id[i]
    line <- samples$pdx_line[i]
    stat <- samples$chr8q_status[i]
    bins <- read_chr8_bins(lib, config)
    segs <- read_chr8_segs(lib, config)
    parm <- read_params(lib, config)
    fs   <- read_flagstat(lib)
    plot_chr8_panel(lib, line, stat, bins, segs, parm)

    p_bins <- bins[bins$end   <= CENTRO_START, ]
    q_bins <- bins[bins$start >= CENTRO_END,   ]
    med_p <- median(p_bins$logR, na.rm = TRUE)
    med_q <- median(q_bins$logR, na.rm = TRUE)
    summary_rows[[i]] <- data.frame(
      library_id        = lib,
      sample_id         = samples$sample_id[i],
      pdx_line          = line,
      chr8q_status      = stat,
      config            = config,
      mapped_reads      = fs$mapped,
      dup_pct           = round(fs$dup_pct, 2),
      ichor_tf          = round(parm$tf, 4),
      ichor_ploidy      = round(parm$ploidy, 3),
      n_bins_chr8p      = nrow(p_bins),
      n_bins_chr8q      = nrow(q_bins),
      median_logR_chr8p = round(med_p, 4),
      median_logR_chr8q = round(med_q, 4),
      delta_q_minus_p   = round(med_q - med_p, 4),
      stringsAsFactors  = FALSE
    )
  }
  mtext(sprintf("chr8 log2 ratio  --  ichorCNA %s, default-low-input, 1 Mb bins",
                config),
        outer = TRUE, line = 0.6, cex = 1.05, font = 2)
  dev.off()
  cat("wrote:", png_path, "\n")

  summary_df <- do.call(rbind, summary_rows)
  tsv_path <- file.path(PLOTS_DIR, sprintf("summary_%s.tsv", config))
  write.table(summary_df, tsv_path, sep = "\t", quote = FALSE,
              row.names = FALSE)
  cat("wrote:", tsv_path, "\n")
  invisible(summary_df)
}

for (cfg in configs) run_one_config(cfg)
