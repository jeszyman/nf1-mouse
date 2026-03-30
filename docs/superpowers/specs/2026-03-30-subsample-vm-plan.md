# Subsampled PDX Alignment on New GCP VM

**Purpose:** Run BISCUIT alignment + ichorCNA on subsampled FASTQs to produce copy number results for CTF Year 4 progress report (due 2026-04-01).

**Why subsampled:** Full-size alignment (~380M read pairs/sample) takes ~6 days/sample. Subsampling to 50M read pairs (~13%) yields ~2.5x human coverage after disambiguate — sufficient for arm-level CNA calling with ichorCNA — and finishes in hours, not weeks.

---

## Prerequisites (user handles before Claude starts)

1. **New GCP VM provisioned** with adequate resources (recommend n2-highmem-16 or larger; 6 parallel BISCUIT processes at 8 threads × 25 GB RAM = 150 GB minimum)
2. **References copied to new VM** via GCS intermediary:
   - From production VM: `gsutil -m cp -r /mnt/data/projects/nf1-mouse/ref/biscuit/ gs://<bucket>/nf1-mouse-refs/`
   - On new VM: `gsutil -m cp -r gs://<bucket>/nf1-mouse-refs/ /mnt/data/projects/nf1-mouse/ref/biscuit/`
   - Required refs:
     - `/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/` (hg38 no-alt + BISCUIT index files + .fai)
     - `/mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/` (mm10 + BISCUIT index files + .fai)
3. **Repo cloned** on new VM: `git clone <nf1-mouse repo URL>`
4. **GCS mount or gsutil access** to FASTQ bucket (FASTQs at `gs://chaudhuri-lab-bucket1/jeszyman/data/nf1/Project_JohnShern_CS039619_112EMseqlib_061225/Flowcell_22WFNNLT4/`)

---

## Samples

6 WU-487 terminal bleeds (Chr8q-WT):

| Library | GCS Sample Dir | Mouse | cfDNA (pg/uL) |
|---------|---------------|-------|---------------|
| lib0622 | Sample_95_nf1_lib95 | mou0001 | 467 |
| lib0623 | Sample_96_nf1_lib96 | mou0002 | 536 |
| lib0624 | Sample_97_nf1_lib97 | mou0003 | 111 |
| lib0625 | Sample_98_nf1_lib98 | mou0004 | 1720 |
| lib0626 | Sample_99_nf1_lib99 | mou0005 | 2380 |
| lib0627 | Sample_100_nf1_lib100 | mou0006 | 387 |

FASTQ naming convention: `{num}_nf1_lib{num}_S{lane}_R{1,2}_001.fastq.gz`

---

## Phase 1: Environment Setup

```bash
# Create conda env from repo yaml
conda env create -f config/nf1-mouse-conda-env.yaml
# Env contains: biscuit, samtools, pysam

# Verify references exist and are indexed
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/*.bis.amb
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/*.bis.amb
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/*.fai
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/*.fai

# Install seqkit for subsampling
conda install -n nf1-mouse -c bioconda seqkit

# Install ichorCNA
# ichorCNA is an R package — install via R
conda install -n nf1-mouse -c bioconda -c conda-forge r-base bioconductor-hmmcopy r-optparse
conda run -n nf1-mouse R -e 'if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("HMMcopy")'
# Install ichorCNA from GitHub
conda run -n nf1-mouse R -e 'install.packages("devtools", repos="https://cran.r-project.org"); devtools::install_github("broadinstitute/ichorCNA")'

# Also need readCounter (HMMcopy suite) for generating wig files
conda install -n nf1-mouse -c bioconda hmmcopy
```

Create directory structure:
```bash
mkdir -p /mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M/{lib0622,lib0623,lib0624,lib0625,lib0626,lib0627}
mkdir -p /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/{strategy1,strategy3,comparison,ichorcna}
```

---

## Phase 2: FASTQ Subsample

Pull FASTQs from GCS and subsample to 50M read pairs (100M individual reads per file) using `seqkit head`. Process R1 and R2 identically to preserve pair synchronization.

```bash
GCS_BASE="gs://chaudhuri-lab-bucket1/jeszyman/data/nf1/Project_JohnShern_CS039619_112EMseqlib_061225/Flowcell_22WFNNLT4"
LOCAL_FULL="/mnt/data/projects/nf1-mouse/emseq/fastqs"
LOCAL_SUB="/mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M"

declare -A SAMPLE_MAP=(
  [lib0622]="Sample_95_nf1_lib95"
  [lib0623]="Sample_96_nf1_lib96"
  [lib0624]="Sample_97_nf1_lib97"
  [lib0625]="Sample_98_nf1_lib98"
  [lib0626]="Sample_99_nf1_lib99"
  [lib0627]="Sample_100_nf1_lib100"
)

for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  SAMPLE_DIR="${SAMPLE_MAP[$LIB]}"
  mkdir -p "${LOCAL_FULL}/${LIB}" "${LOCAL_SUB}/${LIB}"

  # Pull full FASTQs from GCS (if not already local)
  gsutil -m cp "${GCS_BASE}/${SAMPLE_DIR}/*_R1_001.fastq.gz" "${LOCAL_FULL}/${LIB}/"
  gsutil -m cp "${GCS_BASE}/${SAMPLE_DIR}/*_R2_001.fastq.gz" "${LOCAL_FULL}/${LIB}/"

  # Subsample: take first 50M reads (= 50M read pairs since R1 and R2 are synchronized)
  R1=$(ls "${LOCAL_FULL}/${LIB}"/*_R1_001.fastq.gz)
  R2=$(ls "${LOCAL_FULL}/${LIB}"/*_R2_001.fastq.gz)

  conda run -n nf1-mouse seqkit head -n 50000000 "$R1" -o "${LOCAL_SUB}/${LIB}/$(basename $R1)"
  conda run -n nf1-mouse seqkit head -n 50000000 "$R2" -o "${LOCAL_SUB}/${LIB}/$(basename $R2)"

  # Verify pair count matches
  R1_COUNT=$(conda run -n nf1-mouse seqkit stats -T "${LOCAL_SUB}/${LIB}/$(basename $R1)" | tail -1 | cut -f4)
  R2_COUNT=$(conda run -n nf1-mouse seqkit stats -T "${LOCAL_SUB}/${LIB}/$(basename $R2)" | tail -1 | cut -f4)
  echo "${LIB}: R1=${R1_COUNT} R2=${R2_COUNT}"

  # Delete full FASTQs after subsampling to save disk
  rm "${LOCAL_FULL}/${LIB}"/*_R1_001.fastq.gz "${LOCAL_FULL}/${LIB}"/*_R2_001.fastq.gz
done
```

**Note on `seqkit head`:** Takes the first N records from the FASTQ. Since R1 and R2 are written in lockstep by the sequencer, taking the first 50M from each preserves pairing. No need for random sampling — systematic bias from positional subsampling is negligible for alignment/CNA purposes.

---

## Phase 3: Strategy 1 — Direct hg38 Alignment

Align subsampled reads to hg38 with BISCUIT. Includes the header truncation workaround (`samtools view -bS -t ref.fai`).

```bash
HG38_REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
HG38_FAI="${HG38_REF}.fai"
THREADS=8

for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
(
  FASTQ_DIR="/mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M/${LIB}"
  OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy1/${LIB}"
  mkdir -p "$OUTDIR"

  R1=$(ls "${FASTQ_DIR}"/*_R1_001.fastq.gz)
  R2=$(ls "${FASTQ_DIR}"/*_R2_001.fastq.gz)

  conda run -n nf1-mouse bash -c "
    biscuit align -@ ${THREADS} '${HG38_REF}' '${R1}' '${R2}' \
      | samtools view -bS -t '${HG38_FAI}' \
      | samtools sort -@ ${THREADS} -o '${OUTDIR}/${LIB}.hg38.sorted.bam'

    samtools index '${OUTDIR}/${LIB}.hg38.sorted.bam'
    samtools flagstat '${OUTDIR}/${LIB}.hg38.sorted.bam' > '${OUTDIR}/${LIB}.hg38.flagstat.txt'
    samtools idxstats '${OUTDIR}/${LIB}.hg38.sorted.bam' > '${OUTDIR}/${LIB}.hg38.idxstats.txt'
  "
  echo "Strategy 1 complete: ${LIB}"
) &
done
wait
echo "All Strategy 1 alignments complete."
```

**Checkpoint:** After lib0622 completes, inspect the BAM header and flagstat before waiting for the rest:
```bash
samtools view -H /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy1/lib0622/lib0622.hg38.sorted.bam | head -20
cat /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy1/lib0622/lib0622.hg38.flagstat.txt
```

---

## Phase 4: Strategy 3 — Disambiguate

Strategy 3 reuses the hg38 BAMs from Strategy 1. Only mm10 alignment is new. Then run disambiguate.py to classify reads.

### 4a: mm10 alignment

```bash
MM10_REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa"
MM10_FAI="${MM10_REF}.fai"
THREADS=8

for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
(
  FASTQ_DIR="/mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M/${LIB}"
  OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"
  mkdir -p "$OUTDIR"

  R1=$(ls "${FASTQ_DIR}"/*_R1_001.fastq.gz)
  R2=$(ls "${FASTQ_DIR}"/*_R2_001.fastq.gz)

  # Symlink the Strategy 1 hg38 BAM into Strategy 3 directory
  ln -sf "../../strategy1/${LIB}/${LIB}.hg38.sorted.bam" "${OUTDIR}/${LIB}.hg38.sorted.bam"
  ln -sf "../../strategy1/${LIB}/${LIB}.hg38.sorted.bam.bai" "${OUTDIR}/${LIB}.hg38.sorted.bam.bai"

  # Align to mm10 (with header truncation workaround)
  conda run -n nf1-mouse bash -c "
    biscuit align -@ ${THREADS} '${MM10_REF}' '${R1}' '${R2}' \
      | samtools view -bS -t '${MM10_FAI}' \
      | samtools sort -@ ${THREADS} -o '${OUTDIR}/${LIB}.mm10.sorted.bam'

    samtools index '${OUTDIR}/${LIB}.mm10.sorted.bam'
  "
  echo "mm10 alignment complete: ${LIB}"
) &
done
wait
echo "All mm10 alignments complete."
```

### 4b: Disambiguate

```bash
REPO_DIR="$HOME/repos/nf1-mouse"

for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
(
  OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"

  conda run -n nf1-mouse python3 "${REPO_DIR}/dev/pdx-read-handling/disambiguate.py" \
    --human-bam "${OUTDIR}/${LIB}.hg38.sorted.bam" \
    --mouse-bam "${OUTDIR}/${LIB}.mm10.sorted.bam" \
    --out-prefix "${OUTDIR}/${LIB}" \
    --threads 8

  conda run -n nf1-mouse bash -c "
    samtools flagstat '${OUTDIR}/${LIB}.hg38.sorted.bam' > '${OUTDIR}/${LIB}.hg38.flagstat.txt'
    samtools flagstat '${OUTDIR}/${LIB}.mm10.sorted.bam' > '${OUTDIR}/${LIB}.mm10.flagstat.txt'
    samtools flagstat '${OUTDIR}/${LIB}.disambig.human.bam' > '${OUTDIR}/${LIB}.disambig.human.flagstat.txt'
  "
  echo "Disambiguate complete: ${LIB}"
) &
done
wait
echo "All disambiguations complete."
```

---

## Phase 5: Comparison Metrics

Adapted from `04_compare.sh`. Collects flagstat and disambiguate summary stats across both strategies.

```bash
BASEDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M"
OUTFILE="${BASEDIR}/comparison/comparison.tsv"
mkdir -p "${BASEDIR}/comparison"

echo -e "sample\tstrategy\ttotal_read_pairs\tmapped_reads\tmapq30_reads\thuman_reads\tambiguous_reads\tambiguous_pct" > "$OUTFILE"

for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  TOTAL_PAIRS=50000000  # Known from subsampling

  # Strategy 1
  S1_BAM="$BASEDIR/strategy1/${LIB}/${LIB}.hg38.sorted.bam"
  if [ -f "$S1_BAM" ]; then
    MAPPED=$(conda run -n nf1-mouse samtools view -c -F 4 "$S1_BAM")
    MAPQ30=$(conda run -n nf1-mouse samtools view -c -F 4 -q 30 "$S1_BAM")
    echo -e "${LIB}\tstrategy1_direct\t${TOTAL_PAIRS}\t${MAPPED}\t${MAPQ30}\t${MAPQ30}\t0\t0.0" >> "$OUTFILE"
  fi

  # Strategy 3
  S3_SUMMARY="$BASEDIR/strategy3/${LIB}/${LIB}.disambig.summary.txt"
  S3_BAM="$BASEDIR/strategy3/${LIB}/${LIB}.disambig.human.bam"
  if [ -f "$S3_SUMMARY" ]; then
    HUMAN=$(awk -F'\t' '$1=="human_assigned"{print $2}' "$S3_SUMMARY")
    AMBIG=$(awk -F'\t' '$1=="ambiguous"{print $2}' "$S3_SUMMARY")
    AMBIG_PCT=$(awk -F'\t' '$1=="ambiguous_pct"{print $2}' "$S3_SUMMARY")
    MAPQ30=$(conda run -n nf1-mouse samtools view -c -F 4 -q 30 "$S3_BAM")
    echo -e "${LIB}\tstrategy3_disambig\t${TOTAL_PAIRS}\t${MAPPED}\t${MAPQ30}\t${HUMAN}\t${AMBIG}\t${AMBIG_PCT}" >> "$OUTFILE"
  fi
done

echo "Comparison written to $OUTFILE"
column -t "$OUTFILE"
```

---

## Phase 6: ichorCNA

### 6a: Generate read count wig files

ichorCNA requires read depth in fixed windows. Use HMMcopy's `readCounter` at 1 Mb resolution.

```bash
BASEDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M"
ICHOR_DIR="${BASEDIR}/ichorcna"
mkdir -p "$ICHOR_DIR"

# Generate wig files for both strategies
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  # Strategy 1 BAM
  S1_BAM="$BASEDIR/strategy1/${LIB}/${LIB}.hg38.sorted.bam"
  conda run -n nf1-mouse readCounter --window 1000000 --quality 20 \
    --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" \
    "$S1_BAM" > "${ICHOR_DIR}/${LIB}.strategy1.wig"

  # Strategy 3 disambiguated human BAM
  S3_BAM="$BASEDIR/strategy3/${LIB}/${LIB}.disambig.human.bam"
  if [ -f "$S3_BAM" ]; then
    conda run -n nf1-mouse readCounter --window 1000000 --quality 20 \
      --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" \
      "$S3_BAM" > "${ICHOR_DIR}/${LIB}.strategy3.wig"
  fi
done
```

### 6b: Run ichorCNA

**Checkpoint:** Run on lib0622 first to verify ichorCNA works without a panel of normals. If the output looks reasonable, batch the rest.

```bash
ICHOR_DIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/ichorcna"

# ichorCNA R script location (from installed package)
ICHORCNA_RSCRIPT=$(conda run -n nf1-mouse R -e 'cat(system.file("scripts", "runIchorCNA.R", package="ichorCNA"))' 2>/dev/null)

for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  for STRAT in strategy1 strategy3; do
    WIG="${ICHOR_DIR}/${LIB}.${STRAT}.wig"
    OUTDIR_ICHOR="${ICHOR_DIR}/${LIB}_${STRAT}"
    mkdir -p "$OUTDIR_ICHOR"

    if [ -f "$WIG" ]; then
      conda run -n nf1-mouse Rscript "$ICHORCNA_RSCRIPT" \
        --id "${LIB}_${STRAT}" \
        --WIG "$WIG" \
        --outDir "$OUTDIR_ICHOR" \
        --includeHOMD FALSE \
        --chrs "c(1:22, 'X')" \
        --chrTrain "c(1:22)" \
        --estimateNormal TRUE \
        --estimatePloidy TRUE \
        --estimateScPrevalence FALSE \
        --txnE 0.9999 \
        --txnStrength 10000 \
        --genomeBuild "hg38" \
        --genomeStyle "UCSC" \
        2>&1 | tee "${OUTDIR_ICHOR}/ichorCNA.log"
    fi
  done
done
```

**Note:** Running without `--normalPanel` (no PoN). ichorCNA falls back to GC/mappability bias correction. Results may be noisier but should show arm-level events clearly at ~2.5x coverage. 81 healthy human EM-seq controls exist on GCS for future PoN generation.

**GC/mappability files:** ichorCNA requires pre-computed GC and mappability wig files for the chosen bin size. These ship with the ichorCNA package. Locate them with:
```bash
conda run -n nf1-mouse R -e 'cat(system.file("extdata", package="ichorCNA"))'
```
Then pass `--gcWig` and `--mapWig` for hg38 at 1000000 bin size (e.g., `gc_hg38_1000kb.wig`, `map_hg38_1000kb.wig`). If the package doesn't include hg38 wigs, download from the ichorCNA GitHub `inst/extdata/` directory.

---

## Phase 7: Results Collection

Gather all outputs into a single results directory for easy review.

```bash
BASEDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M"
RESULTS="${BASEDIR}/results_summary"
mkdir -p "$RESULTS"

# Copy comparison table
cp "${BASEDIR}/comparison/comparison.tsv" "$RESULTS/"

# Copy all flagstat files
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  cp "${BASEDIR}/strategy1/${LIB}/${LIB}.hg38.flagstat.txt" "$RESULTS/${LIB}.s1.flagstat.txt"
  cp "${BASEDIR}/strategy3/${LIB}/${LIB}.disambig.summary.txt" "$RESULTS/${LIB}.s3.disambig.summary.txt" 2>/dev/null
  cp "${BASEDIR}/strategy3/${LIB}/${LIB}.disambig.human.flagstat.txt" "$RESULTS/${LIB}.s3.human.flagstat.txt" 2>/dev/null
done

# Copy ichorCNA plots and segment files
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  for STRAT in strategy1 strategy3; do
    ICHOR_OUT="${BASEDIR}/ichorcna/${LIB}_${STRAT}"
    if [ -d "$ICHOR_OUT" ]; then
      cp "${ICHOR_OUT}"/*.pdf "$RESULTS/" 2>/dev/null
      cp "${ICHOR_OUT}"/*.seg "$RESULTS/" 2>/dev/null
      cp "${ICHOR_OUT}"/*.params.txt "$RESULTS/" 2>/dev/null
    fi
  done
done

echo "Results collected in $RESULTS"
ls -la "$RESULTS"
```

---

## Key things for Claude on the new VM to know

1. **Header truncation workaround is mandatory.** BISCUIT 1.8.0 truncates `@SQ` lines. Always pipe through `samtools view -bS -t ref.fai` between `biscuit align` and `samtools sort`. Applies to both hg38 AND mm10 alignments.

2. **Strategy 3 reuses Strategy 1 hg38 BAMs.** Do NOT re-align to hg38 — symlink from `strategy1/` into `strategy3/` dirs.

3. **disambiguate.py** is at `dev/pdx-read-handling/disambiguate.py` in the repo. It requires pysam (in the nf1-mouse conda env). It is O(1) RAM — no memory concerns even with large BAMs.

4. **ichorCNA checkpoint:** Run lib0622 through ichorCNA first. If the CNA profile and tumor fraction estimate look reasonable (expect diploid-ish with possible CDKN2A deletion on 9p), proceed with the rest. If it fails or produces degenerate output, troubleshoot before batching.

5. **Disk cleanup:** Delete full-size FASTQs after subsampling. Subsampled FASTQs are ~3-5 GB each, BAMs ~5-10 GB each. Total disk for 6 samples through both strategies: ~200-300 GB.

6. **Conda env yaml** is at `config/nf1-mouse-conda-env.yaml` (biscuit, samtools, pysam). You'll need to add seqkit, R, HMMcopy, and ichorCNA on top.

7. **These are Chr8q-WT samples.** Don't expect Chr8q gain. Look for other MPNST-associated CNAs (CDKN2A/9p deletion, SUZ12/17q loss, etc.) as proof that the pipeline detects real tumor biology.
