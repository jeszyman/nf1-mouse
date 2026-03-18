# Direct vs Disambiguate PDX Read Assignment — Full Workup Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compare Strategy 1 (direct-to-human) and Strategy 3 (disambiguate) on test data, then run both on all 18 mouse samples via a new Mayo GCP VM to pick a winner from full data.

**Architecture:** Four phases — (A) transfer Strategy 1/3 content from old 3-strategy heading into new bioinfo-dev, close out old heading, (B) evaluate/harden run scripts for production (update reference paths to NCBI hg38 + UCSC mm10), (C) provision Mayo GCP VM with right specs, (D) download references, pull FASTQs, and run all 18 samples on VM. ICRG (Strategy 2) is bookmarked in its own analysis heading.

**Tech Stack:** R (biotools conda env), bash, BISCUIT, pysam, GCP (vm-spec-wizard, service catalog), gsutil

---

## Reference Standardization (CRITICAL)

The test runs used Ensembl references (numeric chromosome naming: `1`, `2`, ...). Production will use:
- **Human**: NCBI GRCh38 no-alt analysis set (NO decoys) — `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz` from NCBI `seqs_for_alignment_pipelines.ucsc_ids/`. 195 contigs (25 primary + 169 unlocalized/unplaced + EBV), chr-prefixed, no ALTs, no decoys. Needs download + BISCUIT indexing.
- **Mouse**: UCSC mm10 from `hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz` — 66 contigs (22 primary + 22 random + 22 unplaced), chr-prefixed, no ALTs, no patches, no strain-specific contigs. Needs download + BISCUIT indexing.

Symmetric references: neither has decoys or ALT contigs. Eliminates reference asymmetry bias in disambiguation.

Reference download, BISCUIT indexing, and all alignment runs happen on the VM (Phase D). No local re-run needed.

---

## Key Data Context (from prior test runs, Ensembl refs — for orientation only)

Prior comparison.tsv from test runs (3M read-pair subsets of lib0625/lib0626, Ensembl refs):

| sample  | strategy          | total_reads | mapped   | MAPQ30  | human_reads |
|---------|-------------------|-------------|----------|---------|-------------|
| lib0625 | strategy1_direct  | 7,273,869   | 7,032,278| 753,156 | 753,156     |
| lib0625 | strategy3_disambig| 7,273,869   | 7,032,278| 113,339 | 1,118,997   |
| lib0626 | strategy1_direct  | 7,250,678   | 7,123,450| 890,806 | 890,806     |
| lib0626 | strategy3_disambig| 7,250,678   | 7,123,450| 225,330 | 1,344,539   |

Disambiguate summaries (read-pair level):
- lib0625: 505K human (16.8%), 2.14M mouse (71.4%), 352K ambiguous (11.7%)
- lib0626: 615K human (20.5%), 2.02M mouse (67.4%), 364K ambiguous (12.1%)

Critical observation: Strategy 1 MAPQ30 (753K) > Strategy 3 MAPQ30 (113K), but Strategy 3 total human reads (1.12M) >> Strategy 1 (753K). The discrepancy is because disambiguate assigns reads to human by NM-tag comparison, but those reads retain their original MAPQ from the hg38 alignment — many are MAPQ<30 because they have alternative alignments within the human genome. This is expected and means MAPQ filtering applies differently in the two strategies.

---

## Phase A: Org Restructuring + Bioinfo-dev Setup (jeff-beast)

### Task 1: Transfer content from old heading and create dev directory

Transfer Strategy 1 and Strategy 3 content from `* PDX read handling evaluation` into the new bioinfo-dev heading. Close out the old heading. The new heading focuses exclusively on direct-vs-disambiguate comparison.

**Files:**
- Modify: `nf1-mouse.org` (restructure headings)
- Create: `dev/direct-vs-disambiguate/00_load_data.R`

- [ ] **Step 1: Create org heading and dev directory (org-edit skill)**

- [ ] **Step 2: Create 00_load_data.R

**Files:**
- Create: `dev/direct-vs-disambiguate/00_load_data.R`

- [ ] **Step 1: Create directory and 00_load_data.R**

```r
# 00_load_data.R — Load comparison data from Strategy 1 and Strategy 3
library(tidyverse)

base_dir <- "/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling"

# Comparison TSV
comparison <- read_tsv(file.path(base_dir, "comparison.tsv")) %>%
  filter(strategy != "strategy2_icrg")

# Disambiguate summaries
disambig_summary <- bind_rows(
  read_tsv(file.path(base_dir, "strategy3/lib0625/lib0625.disambig.summary.txt"),
           col_names = c("metric", "value")) %>% mutate(sample = "lib0625"),
  read_tsv(file.path(base_dir, "strategy3/lib0626/lib0626.disambig.summary.txt"),
           col_names = c("metric", "value")) %>% mutate(sample = "lib0626")
)

# Flagstat parser
parse_flagstat <- function(path) {
  lines <- readLines(path)
  vals <- as.numeric(str_extract(lines, "^\\d+"))
  tibble(
    total = vals[1], primary = vals[2], secondary = vals[3],
    supplementary = vals[4], duplicates = vals[5],
    mapped = vals[7], primary_mapped = vals[8],
    properly_paired = vals[10], singletons = vals[12]
  )
}

# Load flagstats for both strategies
flagstats <- bind_rows(
  parse_flagstat(file.path(base_dir, "strategy1/lib0625/lib0625.hg38.flagstat.txt")) %>%
    mutate(sample = "lib0625", strategy = "direct"),
  parse_flagstat(file.path(base_dir, "strategy1/lib0626/lib0626.hg38.flagstat.txt")) %>%
    mutate(sample = "lib0626", strategy = "direct"),
  parse_flagstat(file.path(base_dir, "strategy3/lib0625/lib0625.disambig.human.flagstat.txt")) %>%
    mutate(sample = "lib0625", strategy = "disambiguate_human"),
  parse_flagstat(file.path(base_dir, "strategy3/lib0626/lib0626.disambig.human.flagstat.txt")) %>%
    mutate(sample = "lib0626", strategy = "disambiguate_human"),
  parse_flagstat(file.path(base_dir, "strategy3/lib0625/lib0625.hg38.flagstat.txt")) %>%
    mutate(sample = "lib0625", strategy = "disambiguate_hg38_full"),
  parse_flagstat(file.path(base_dir, "strategy3/lib0626/lib0626.hg38.flagstat.txt")) %>%
    mutate(sample = "lib0626", strategy = "disambiguate_hg38_full")
)

# Plot theme
theme_white <- theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

dir.create("data", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
```

- [ ] **Step 2: Run to verify data loads**

Run: `cd /home/jeszyman/repos/nf1-mouse/dev/direct-vs-disambiguate && conda run -n biotools Rscript 00_load_data.R 2>&1`

### Task 2: MAPQ distribution comparison

**Files:**
- Create: `dev/direct-vs-disambiguate/01_mapq_comparison.R`

This is the core comparison. For each sample, extract per-read MAPQ from:
- Strategy 1: the full hg38 BAM (all reads, human and mouse together)
- Strategy 3: the disambiguate human-only BAM (only NM-assigned human reads)

- [ ] **Step 1: Write MAPQ extraction and comparison script**

```r
# 01_mapq_comparison.R — MAPQ distributions: direct vs disambiguate human BAMs
source("00_load_data.R")

# Extract MAPQ distributions from BAMs using samtools
extract_mapq <- function(bam_path, label) {
  cmd <- sprintf("samtools view -F 4 '%s' | cut -f5", bam_path)
  mapq <- as.integer(system(cmd, intern = TRUE))
  tibble(mapq = mapq, source = label)
}

samples <- c("lib0625", "lib0626")
mapq_data <- list()

for (s in samples) {
  # Strategy 1: all mapped reads from direct hg38 alignment
  s1_bam <- file.path(base_dir, sprintf("strategy1/%s/%s.hg38.sorted.bam", s, s))
  # Strategy 3: human-assigned reads only
  s3_bam <- file.path(base_dir, sprintf("strategy3/%s/%s.disambig.human.bam", s, s))

  mapq_data[[paste0(s, "_direct")]] <- extract_mapq(s1_bam, "direct") %>% mutate(sample = s)
  mapq_data[[paste0(s, "_disambig")]] <- extract_mapq(s3_bam, "disambig_human") %>% mutate(sample = s)
}

d_mapq <- bind_rows(mapq_data)
write_tsv(d_mapq, "data/mapq_distributions.tsv")

# Plot: MAPQ histograms faceted by sample and strategy
p1 <- ggplot(d_mapq, aes(x = mapq, fill = source)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c(direct = "#e66101", disambig_human = "#5e3c99")) +
  facet_grid(sample ~ source, scales = "free_y") +
  labs(x = "MAPQ", y = "Read count",
       title = "MAPQ distribution: Direct vs Disambiguate human reads") +
  theme_white
ggsave("plots/mapq_comparison.png", p1, width = 12, height = 8, dpi = 150)

# Summary table
mapq_summary <- d_mapq %>%
  group_by(sample, source) %>%
  summarise(
    n_reads = n(),
    pct_mapq0 = 100 * mean(mapq == 0),
    pct_mapq30plus = 100 * mean(mapq >= 30),
    pct_mapq60 = 100 * mean(mapq == 60),
    median_mapq = median(mapq),
    .groups = "drop"
  )
write_tsv(mapq_summary, "data/mapq_summary.tsv")
print(mapq_summary)
```

- [ ] **Step 2: Run MAPQ comparison**

Run: `cd /home/jeszyman/repos/nf1-mouse/dev/direct-vs-disambiguate && conda run -n biotools Rscript 01_mapq_comparison.R 2>&1`

- [ ] **Step 3: Review plots and summary, interpret results**

### Task 3: Read-level overlap between strategies

**Files:**
- Create: `dev/direct-vs-disambiguate/02_read_overlap.R`

Key question: Of the MAPQ>=30 reads in Strategy 1, how many are also classified as human by Strategy 3? Of the human reads Strategy 3 finds, how many would pass MAPQ>=30 in Strategy 1?

- [ ] **Step 1: Write read overlap analysis**

```r
# 02_read_overlap.R — Read-level overlap between strategies
source("00_load_data.R")

# For each sample: extract read names from both strategies
get_qnames <- function(bam_path, min_mapq = 0) {
  cmd <- sprintf("samtools view -F 4 -q %d '%s' | cut -f1 | sort -u", min_mapq, bam_path)
  system(cmd, intern = TRUE)
}

for (s in c("lib0625", "lib0626")) {
  cat(sprintf("\n=== %s ===\n", s))

  s1_bam <- file.path(base_dir, sprintf("strategy1/%s/%s.hg38.sorted.bam", s, s))
  s3_human <- file.path(base_dir, sprintf("strategy3/%s/%s.disambig.human.bam", s, s))

  # Strategy 1: MAPQ>=30 read names (the "human" set for direct alignment)
  s1_q30 <- get_qnames(s1_bam, min_mapq = 30)
  # Strategy 3: all human-assigned read names (any MAPQ)
  s3_all <- get_qnames(s3_human)
  # Strategy 3: MAPQ>=30 human reads
  s3_q30 <- get_qnames(s3_human, min_mapq = 30)

  # Overlaps
  overlap_s1q30_s3all <- length(intersect(s1_q30, s3_all))
  s1_only <- length(setdiff(s1_q30, s3_all))
  s3_only <- length(setdiff(s3_all, s1_q30))

  cat(sprintf("Strategy 1 MAPQ>=30: %d reads\n", length(s1_q30)))
  cat(sprintf("Strategy 3 human (any MAPQ): %d reads\n", length(s3_all)))
  cat(sprintf("Strategy 3 human MAPQ>=30: %d reads\n", length(s3_q30)))
  cat(sprintf("Overlap (S1 Q30 ∩ S3 human): %d (%.1f%% of S1)\n",
              overlap_s1q30_s3all, 100 * overlap_s1q30_s3all / length(s1_q30)))
  cat(sprintf("S1 Q30 only (not in S3 human): %d — these may be mouse contamination\n", s1_only))
  cat(sprintf("S3 human only (not in S1 Q30): %d — these are low-MAPQ human or missed by S1\n", s3_only))
}
```

- [ ] **Step 2: Run overlap analysis**

Run: `cd /home/jeszyman/repos/nf1-mouse/dev/direct-vs-disambiguate && conda run -n biotools Rscript 02_read_overlap.R 2>&1`

Note: This step will be slow (~5-10 min per sample) because it extracts and sorts millions of read names. Consider whether samtools+awk is faster than R for this.

- [ ] **Step 3: Interpret overlap results**

Key interpretive questions:
- If most S1 MAPQ>=30 reads are also in S3 human: disambiguate is a superset (adds true human reads that S1 misses)
- If S1 has reads NOT in S3: those are potential mouse contamination (mouse reads that map well to human by chance)
- The ratio of S1-only to S3-only informs which strategy is more conservative vs more sensitive

### Task 4: Per-chromosome read distribution

**Files:**
- Create: `dev/direct-vs-disambiguate/03_chrom_distribution.R`

- [ ] **Step 1: Write chromosome distribution comparison**

```r
# 03_chrom_distribution.R — Per-chromosome read counts for both strategies
source("00_load_data.R")

parse_idxstats <- function(path, label, sample) {
  read_tsv(path, col_names = c("chrom", "length", "mapped", "unmapped"),
           col_types = "ciii") %>%
    filter(chrom != "*", mapped > 0) %>%
    mutate(source = label, sample = sample)
}

idx_data <- bind_rows(
  parse_idxstats(file.path(base_dir, "strategy1/lib0625/lib0625.hg38.idxstats.txt"),
                 "direct", "lib0625"),
  parse_idxstats(file.path(base_dir, "strategy1/lib0626/lib0626.hg38.idxstats.txt"),
                 "direct", "lib0626")
)

# For Strategy 3, we need to generate idxstats from the disambiguate human BAM
for (s in c("lib0625", "lib0626")) {
  bam <- file.path(base_dir, sprintf("strategy3/%s/%s.disambig.human.bam", s, s))
  idx_out <- file.path(base_dir, sprintf("strategy3/%s/%s.disambig.human.idxstats.txt", s, s))
  if (!file.exists(idx_out)) {
    system2("samtools", c("idxstats", bam), stdout = idx_out)
  }
  idx_data <- bind_rows(idx_data,
    parse_idxstats(idx_out, "disambig_human", s))
}

# Focus on standard chromosomes
std_chroms <- c(as.character(1:22), "X", "Y")
idx_std <- idx_data %>%
  filter(chrom %in% std_chroms) %>%
  mutate(chrom = factor(chrom, levels = std_chroms))

p <- ggplot(idx_std, aes(x = chrom, y = mapped, fill = source)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c(direct = "#e66101", disambig_human = "#5e3c99")) +
  facet_wrap(~sample, ncol = 1, scales = "free_y") +
  labs(x = "Chromosome", y = "Mapped reads",
       title = "Per-chromosome reads: Direct vs Disambiguate") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/chrom_distribution.png", p, width = 14, height = 8, dpi = 150)

# Ratio plot: disambiguate / direct per chromosome
ratio_data <- idx_std %>%
  select(chrom, mapped, source, sample) %>%
  pivot_wider(names_from = source, values_from = mapped) %>%
  mutate(ratio = disambig_human / direct)

p_ratio <- ggplot(ratio_data, aes(x = chrom, y = ratio)) +
  geom_col(fill = "#5e3c99", alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~sample, ncol = 1) +
  labs(x = "Chromosome", y = "Disambiguate / Direct read ratio",
       title = "Per-chromosome enrichment: Disambiguate vs Direct") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/chrom_ratio.png", p_ratio, width = 14, height = 6, dpi = 150)
```

- [ ] **Step 2: Run chromosome distribution analysis**

Run: `cd /home/jeszyman/repos/nf1-mouse/dev/direct-vs-disambiguate && conda run -n biotools Rscript 03_chrom_distribution.R 2>&1`

- [ ] **Step 3: Interpret chromosome results**

If the ratio is uniform across chromosomes: strategies differ by a global sensitivity threshold.
If specific chromosomes are enriched/depleted: those regions have more/less mouse orthology.

### Task 5: Commit Phase A results

- [ ] **Step 1: Commit dev directory and analysis outputs**

```bash
git add dev/direct-vs-disambiguate/
git commit -m "Add direct-vs-disambiguate bioinfo-dev: MAPQ, overlap, chrom analyses"
```

---

## Phase B: Script Evaluation for Production

### Task 6: Evaluate run scripts for full-FASTQ readiness

**Files:**
- Review: `dev/pdx-read-handling/01_strategy1_direct.sh`
- Review: `dev/pdx-read-handling/03_strategy3_disambiguate.sh`
- Review: `dev/pdx-read-handling/disambiguate.py`

Current scripts were tested on 3M read-pair subsets (~6 GB FASTQs). Full samples are ~50-60 GB each (18 samples total).

- [ ] **Step 1: Identify production scaling issues**

Known issues to evaluate:
1. **disambiguate.py memory**: Loads ALL read names into RAM as Python sets. At 3M reads, ~500 MB. Full samples (~200M reads each?) could need 30+ GB just for read-name sets. Need to estimate actual read counts from full FASTQs.
2. **BISCUIT alignment RAM**: ~18-21 GB per process (from prior OOM crash). Cannot parallelize BISCUIT on a machine with <42 GB RAM without sequential runs.
3. **Disk space**: Each strategy produces sorted BAMs. Strategy 3 produces 5 BAMs per sample (hg38, mm10, human, mouse, ambiguous). Need to estimate total disk for 18 samples x 2 strategies.
4. **FASTQ glob pattern**: Scripts use `ls *_R1_001.fastq.gz` — verify this matches the GCS bucket naming convention. The data map says: `{num}_nf1_lib{num}_S{lane}_R{1,2}_001.fastq.gz`.
5. **Script hardcoded paths**: All scripts use `/mnt/data/projects/nf1-mouse/` — need to parameterize or document the VM disk layout.

- [ ] **Step 2: Estimate full sample sizes**

Check GCS for actual FASTQ sizes:
```bash
gsutil ls -l gs://chaudhuri-lab-bucket1/jeszyman/data/nf1/Project_JohnShern_CS039619_112EMseqlib_061225/Flowcell_22WFNNLT4/Sample_95_nf1_lib95/ 2>&1
```

- [ ] **Step 3: Decision — fix disambiguate.py memory or accept RAM requirement**

Options:
a) Accept high RAM: provision VM with 64+ GB → simplest, no code changes
b) Streaming approach: rewrite disambiguate.py to stream-merge name-sorted BAMs instead of loading all read names → lower RAM but more complex
c) Chunk approach: process reads in batches → medium complexity

Recommendation: Option (a) for now — RAM is cheap on GCP, and 18 samples is a one-time run. Document the memory estimate.

- [ ] **Step 4: Decide on script changes (if any)**

Evaluate whether to:
- Add `set -x` for logging
- Add timing (`time` or `date` stamps)
- Write per-sample log files
- Parameterize base directory for VM portability

### Task 7: Estimate compute requirements

- [ ] **Step 1: Calculate from test run data**

From the test run (3M read pairs, ~5.5 hrs total for all 3 strategies):
- Strategy 1 (direct BISCUIT): ~1 hr per sample at 3M reads
- Strategy 3 (disambiguate): ~2.5 hrs per sample at 3M reads (two alignments + python)
- Full samples are ~20-30x larger → estimate 20-30 hrs per sample for Strategy 3

For 18 samples running both strategies sequentially:
- Strategy 1: 18 x ~20 hrs = ~360 hrs (15 days sequential)
- Strategy 3: 18 x ~60 hrs = ~1080 hrs (45 days sequential)
- With 2 parallel streams: ~23 days total
- With 4 parallel streams (32 cores, 84 GB): ~12 days total

Decision needed: how aggressive on parallelization? More cores = more cost but faster.

- [ ] **Step 2: Write VM spec recommendation**

Based on estimates:

| Component | Spec | Rationale |
|-----------|------|-----------|
| Machine   | **n2-highmem-64** | 64 vCPUs, 512 GB RAM — run 6 BISCUIT at 8 threads (126 GB), streaming disambiguate adds <1 GB |
| Data disk | **10 TB SSD** | 18 samples x 2 strategies x ~350 GB BAMs + 1 TB FASTQs + refs + headroom |
| Zone      | **us-west4-a** | 0/3000 CPUs used, 20250/81920 GB SSD used — plenty of room |

With 6 parallel streams: ~2-3 weeks total wall time for both strategies on all 18 samples.

- [ ] **Step 3: Commit any script changes**

---

## Phase C: VM Provisioning

### Task 8: Run VM Spec Wizard

- [ ] **Step 1: Run the wizard**

```bash
pip install -e ~/repos/chaudhuri-lab/vm-spec-wizard/ 2>/dev/null
vm-spec-wizard
```

Use the spec recommendation from Task 7. Target zone: `us-west4-a` (same as existing jeff-nf1, least contention per quota table).

- [ ] **Step 2: Review generated spec file**

Spec file will be at `gs://chaudhuri-lab-bucket1/vm-specs/<name>-<date>.md`

- [ ] **Step 3: Create VM via Service Catalog**

Follow the [Service Catalog](https://console.cloud.google.com/catalog/browse?authuser=1&hl=en&project=aif-usr-p-chaudhuri-lab-83f0) "6.1 Shielded Virtual Machine" workflow. Map spec fields to form fields per `chaudhuri-lab.org` L628-635.

Gotchas:
- c3d machines cannot use "regional data disk"
- Only allowed zones: us-central1-{a,b,c,f}, us-east1-{b,c,d}, us-west4-{a,b,c}

### Task 9: Post-provision VM setup

- [ ] **Step 1: SSH into new VM**

```bash
gcloud_vm_connect <vm-name> <zone>
```

- [ ] **Step 2: Run provision-vm.sh**

This handles conda, base packages, GCS mount, etc.

- [ ] **Step 3: Install project-specific conda envs**

```bash
# On VM:
conda env create -f ~/repos/basecamp/environment.yaml
conda env create -f ~/repos/biotools/environment.yaml
# Or if envs are already on the VM image, just verify:
conda run -n biotools biscuit --version
conda run -n biotools python3 -c "import pysam; print(pysam.__version__)"
```

- [ ] **Step 4: Clone nf1-mouse repo and set up disk layout**

```bash
git clone <repo-url> ~/repos/nf1-mouse
mkdir -p /mnt/data/projects/nf1-mouse/emseq/{fastqs,pdx-read-handling/{strategy1,strategy3}}
mkdir -p /mnt/data/projects/nf1-mouse/ref/{biscuit_hg38,biscuit/mm10}
```

- [ ] **Step 5: Set up references on VM**

Either pull indexed refs from bucket or rebuild:
```bash
# Option A: gsutil rsync refs from bucket (faster if already synced)
gsutil -m rsync -r gs://chaudhuri-lab-bucket1/jeszyman/projects/nf1/ref/ /mnt/data/projects/nf1-mouse/ref/

# Option B: Run reference setup script (slower but guaranteed correct)
cd ~/repos/nf1-mouse/dev/pdx-read-handling
conda run -n emseq bash 00_setup_refs.sh
```

---

## Phase D: Production Run

### Task 10: Download references and sync to VM

- [ ] **Step 1: On VM — download NCBI hg38 no-alt analysis set and BISCUIT-index**

```bash
mkdir -p /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt
cd /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
conda run -n emseq biscuit index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
conda run -n emseq samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

- [ ] **Step 2: On VM — download UCSC mm10 and BISCUIT-index**

```bash
mkdir -p /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10
cd /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10
wget -q https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
conda run -n emseq biscuit index mm10.fa
conda run -n emseq samtools faidx mm10.fa
```

- [ ] **Step 3: Verify both references are chr-prefixed and symmetric (no decoys, no ALTs)**

```bash
head -1 /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
head -1 /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa
# Count contigs — hg38: 195, mm10: 66
grep -c "^>" /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
grep -c "^>" /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa
```

### Task 11: Pull all 18 sample FASTQs to VM

- [ ] **Step 1: Pull FASTQs from GCS bucket**

18 samples: lib95-112 on the bucket (= lib0622-0639 in our naming).

```bash
# On VM — pull all 18 PDX sample FASTQs
for nci_num in $(seq 95 112); do
  sample_dir="Sample_${nci_num}_nf1_lib${nci_num}"
  gsutil -m rsync -r \
    "gs://chaudhuri-lab-bucket1/jeszyman/data/nf1/Project_JohnShern_CS039619_112EMseqlib_061225/Flowcell_22WFNNLT4/${sample_dir}/" \
    "/mnt/data/projects/nf1-mouse/emseq/fastqs/lib$(printf '%04d' $((nci_num + 527)))/"
done
```

Note: NCI lib95 = our lib0622 (offset 527). Verify naming mapping before running.

- [ ] **Step 2: Verify FASTQ integrity (md5 check)**

```bash
for d in /mnt/data/projects/nf1-mouse/emseq/fastqs/lib0*/; do
  cd "$d"
  md5sum -c *.md5 2>&1 | tail -1
done
```

### Task 12: Run Strategy 1 on all 18 samples

- [ ] **Step 1: Run Strategy 1 (direct) in parallel batches**

```bash
# Run 2-4 at a time depending on VM cores
for sample in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627 \
              lib0628 lib0629 lib0630 lib0631 lib0632 lib0633 \
              lib0634 lib0635 lib0636 lib0637 lib0638 lib0639; do
  echo "$sample"
done | xargs -P 2 -I {} bash -c \
  'conda run -n emseq bash ~/repos/nf1-mouse/dev/pdx-read-handling/01_strategy1_direct.sh {} 8 \
   > /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/strategy1/{}.log 2>&1'
```

- [ ] **Step 2: Verify all Strategy 1 outputs exist**

```bash
for s in lib0{622..639}; do
  ls /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/strategy1/$s/*.flagstat.txt 2>/dev/null \
    || echo "MISSING: $s"
done
```

### Task 13: Run Strategy 3 on all 18 samples

- [ ] **Step 1: Run Strategy 3 (disambiguate) in parallel batches**

```bash
for sample in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627 \
              lib0628 lib0629 lib0630 lib0631 lib0632 lib0633 \
              lib0634 lib0635 lib0636 lib0637 lib0638 lib0639; do
  echo "$sample"
done | xargs -P 2 -I {} bash -c \
  'conda run -n emseq bash ~/repos/nf1-mouse/dev/pdx-read-handling/03_strategy3_disambiguate.sh {} 8 \
   > /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/strategy3/{}.log 2>&1'
```

- [ ] **Step 2: Verify all Strategy 3 outputs**

```bash
for s in lib0{622..639}; do
  ls /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/strategy3/$s/*.disambig.summary.txt 2>/dev/null \
    || echo "MISSING: $s"
done
```

### Task 14: Generate full comparison and sync back

- [ ] **Step 1: Run comparison script on all 18 samples**

Modify `04_compare.sh` to loop over all 18 samples instead of just lib0625/lib0626.

- [ ] **Step 2: gsutil_rsync results back to bucket**

```bash
gsutil_rsync /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/ \
  gs://chaudhuri-lab-bucket1/jeszyman/projects/nf1/emseq/pdx-read-handling/
```

- [ ] **Step 3: Pull results to jeff-beast for analysis**

```bash
# From jeff-beast:
gsutil_rsync gs://chaudhuri-lab-bucket1/jeszyman/projects/nf1/emseq/pdx-read-handling/ \
  /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/
```

- [ ] **Step 4: Stop VM to avoid billing**

```bash
gcloud compute instances stop <vm-name> --zone=<zone> --discard-local-ssd=false
```

---

## Decision Points Requiring User Input

1. ~~disambiguate.py memory~~ — RESOLVED: rewritten to streaming approach, O(1) RAM
2. ~~VM spec~~ — RESOLVED: n2-highmem-64, 10 TB SSD, us-west4-a
3. **Task 11 Step 1**: Verify NCI→internal lib naming mapping before bulk FASTQ pull

### Comparison Metrics (resolved)

**3-way comparison:** Strategy 1 (direct MAPQ>=30), Strategy 3A (disambiguate, discard ambiguous), Strategy 3B (disambiguate, ambiguous → human)

**Alignment-level (per sample):**
- Total human reads (raw yield)
- MAPQ distribution of human reads
- Per-chromosome read distribution

**Cross-sample consistency:**
- Human fraction across 6 WU-487 terminals (same PDX, should be similar)
- Human fraction across 12 JH-2-055 serials (should track tumor burden)

**Contamination signal:**
- Use disambiguate mouse-assigned reads as ground truth; count how many of those pass MAPQ>=30 in the direct hg38 BAM — these are confirmed mouse reads that Strategy 1 keeps as "human"

**Future TODOs (in Ideas section of org heading):**
- ichorCNA on 3A vs 3B to see if ambiguous reads affect CNA calls
- Overall CNA comparison across all 3 strategies
- LINE dPCR correlation for serial samples (where dPCR works)
