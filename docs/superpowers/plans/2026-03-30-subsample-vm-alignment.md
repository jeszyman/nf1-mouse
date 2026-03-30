# Subsampled VM Alignment Implementation Plan

> **For agentic workers:** Read this plan top-to-bottom and execute tasks sequentially. Each task has verification steps — do not proceed to the next task until verification passes. This is a pipeline execution plan on a fresh GCP VM, not a software development task.

**Goal:** Produce BISCUIT alignment + ichorCNA copy number results on 50M-read-pair subsamples of 6 WU-487 terminal bleed EM-seq libraries, for the CTF Year 4 progress report.

**Architecture:** Subsample FASTQs from GCS, align to hg38 (Strategy 1) and mm10 (Strategy 3), disambiguate to separate human/mouse reads, run ichorCNA on human BAMs. All work under `/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/`.

**Tech Stack:** BISCUIT 1.8.0, samtools, pysam, seqkit, ichorCNA, HMMcopy readCounter, conda

---

## Prerequisites Checklist

Before starting any task, verify these are done (user handles these):

- [ ] GCP VM is running with >= 150 GB RAM, >= 500 GB disk on `/mnt/data/`
- [ ] References exist at `/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/` and `/mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/`
- [ ] Repo cloned (need `dev/pdx-read-handling/disambiguate.py` and `config/nf1-mouse-conda-env.yaml`)
- [ ] GCS access configured (`gsutil` works for `gs://chaudhuri-lab-bucket1/`)

---

### Task 1: Create conda environment and install all tools

**Files:**
- Read: `config/nf1-mouse-conda-env.yaml`

- [ ] **Step 1: Create base conda env**

```bash
conda env create -f config/nf1-mouse-conda-env.yaml
```

Expected: env `nf1-mouse` created with biscuit, samtools, pysam.

- [ ] **Step 2: Install seqkit for FASTQ subsampling**

```bash
conda install -n nf1-mouse -c bioconda -c conda-forge seqkit -y
```

- [ ] **Step 3: Install R, HMMcopy, and ichorCNA dependencies**

```bash
conda install -n nf1-mouse -c bioconda -c conda-forge r-base r-optparse hmmcopy -y
```

- [ ] **Step 4: Install ichorCNA R package**

```bash
conda run -n nf1-mouse R -e '
if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cran.r-project.org")
BiocManager::install("HMMcopy", update=FALSE, ask=FALSE)
install.packages("devtools", repos="https://cran.r-project.org")
devtools::install_github("broadinstitute/ichorCNA")
'
```

- [ ] **Step 5: Verify all tools are available**

```bash
conda run -n nf1-mouse biscuit --version
conda run -n nf1-mouse samtools --version | head -1
conda run -n nf1-mouse seqkit version
conda run -n nf1-mouse readCounter 2>&1 | head -1
conda run -n nf1-mouse R -e 'library(ichorCNA); cat("ichorCNA loaded\n")'
```

Expected: all commands return version info or success messages, no errors.

- [ ] **Step 6: Verify references exist and are indexed**

```bash
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/*.bis.amb
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa.fai
ls /mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/*.bis.amb
```

Expected: all files exist. If `.fai` is missing, generate with `samtools faidx`. If `.bis.amb` is missing, generate with `biscuit index`.

- [ ] **Step 7: Locate ichorCNA GC/mappability wig files**

```bash
EXTDATA=$(conda run -n nf1-mouse R -e 'cat(system.file("extdata", package="ichorCNA"))' 2>/dev/null)
echo "ichorCNA extdata: $EXTDATA"
ls "$EXTDATA"/*.wig 2>/dev/null || echo "No wig files found in package"
```

If hg38 1Mb wig files are not found (`gc_hg38_1000kb.wig`, `map_hg38_1000kb.wig`), download them:
```bash
EXTDATA=$(conda run -n nf1-mouse R -e 'cat(system.file("extdata", package="ichorCNA"))' 2>/dev/null)
# Check ichorCNA GitHub for the correct download URLs for hg38 1Mb wigs
# They are typically at: https://raw.githubusercontent.com/broadinstitute/ichorCNA/master/inst/extdata/
```

Record the full paths to the GC and map wig files — they are needed in Task 6.

---

### Task 2: Create directory structure and subsample FASTQs

- [ ] **Step 1: Create output directories**

```bash
mkdir -p /mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M/{lib0622,lib0623,lib0624,lib0625,lib0626,lib0627}
mkdir -p /mnt/data/projects/nf1-mouse/emseq/fastqs/{lib0622,lib0623,lib0624,lib0625,lib0626,lib0627}
mkdir -p /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/{strategy1,strategy3,comparison,ichorcna}
```

- [ ] **Step 2: Pull FASTQs from GCS and subsample — lib0622 first as test**

```bash
GCS_BASE="gs://chaudhuri-lab-bucket1/jeszyman/data/nf1/Project_JohnShern_CS039619_112EMseqlib_061225/Flowcell_22WFNNLT4"
LOCAL_FULL="/mnt/data/projects/nf1-mouse/emseq/fastqs"
LOCAL_SUB="/mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M"

LIB="lib0622"
SAMPLE_DIR="Sample_95_nf1_lib95"

gsutil -m cp "${GCS_BASE}/${SAMPLE_DIR}/"*_R1_001.fastq.gz "${LOCAL_FULL}/${LIB}/"
gsutil -m cp "${GCS_BASE}/${SAMPLE_DIR}/"*_R2_001.fastq.gz "${LOCAL_FULL}/${LIB}/"

R1=$(ls "${LOCAL_FULL}/${LIB}"/*_R1_001.fastq.gz)
R2=$(ls "${LOCAL_FULL}/${LIB}"/*_R2_001.fastq.gz)

conda run -n nf1-mouse seqkit head -n 50000000 "$R1" -o "${LOCAL_SUB}/${LIB}/$(basename $R1)"
conda run -n nf1-mouse seqkit head -n 50000000 "$R2" -o "${LOCAL_SUB}/${LIB}/$(basename $R2)"
```

- [ ] **Step 3: Verify lib0622 subsample pair counts match**

```bash
conda run -n nf1-mouse seqkit stats -T "${LOCAL_SUB}/lib0622/"*_R1_001.fastq.gz | tail -1
conda run -n nf1-mouse seqkit stats -T "${LOCAL_SUB}/lib0622/"*_R2_001.fastq.gz | tail -1
```

Expected: both show 50,000,000 reads (or fewer if the original has fewer than 50M — unlikely at ~380M pairs). R1 and R2 counts must match exactly.

- [ ] **Step 4: Delete full-size lib0622 FASTQs to save disk**

```bash
rm "${LOCAL_FULL}/lib0622/"*_R1_001.fastq.gz "${LOCAL_FULL}/lib0622/"*_R2_001.fastq.gz
```

- [ ] **Step 5: Subsample remaining 5 libraries**

Process each one sequentially (GCS download is the bottleneck, not CPU):

```bash
GCS_BASE="gs://chaudhuri-lab-bucket1/jeszyman/data/nf1/Project_JohnShern_CS039619_112EMseqlib_061225/Flowcell_22WFNNLT4"
LOCAL_FULL="/mnt/data/projects/nf1-mouse/emseq/fastqs"
LOCAL_SUB="/mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M"

declare -A SAMPLE_MAP=(
  [lib0623]="Sample_96_nf1_lib96"
  [lib0624]="Sample_97_nf1_lib97"
  [lib0625]="Sample_98_nf1_lib98"
  [lib0626]="Sample_99_nf1_lib99"
  [lib0627]="Sample_100_nf1_lib100"
)

for LIB in lib0623 lib0624 lib0625 lib0626 lib0627; do
  SAMPLE_DIR="${SAMPLE_MAP[$LIB]}"

  gsutil -m cp "${GCS_BASE}/${SAMPLE_DIR}/"*_R1_001.fastq.gz "${LOCAL_FULL}/${LIB}/"
  gsutil -m cp "${GCS_BASE}/${SAMPLE_DIR}/"*_R2_001.fastq.gz "${LOCAL_FULL}/${LIB}/"

  R1=$(ls "${LOCAL_FULL}/${LIB}"/*_R1_001.fastq.gz)
  R2=$(ls "${LOCAL_FULL}/${LIB}"/*_R2_001.fastq.gz)

  conda run -n nf1-mouse seqkit head -n 50000000 "$R1" -o "${LOCAL_SUB}/${LIB}/$(basename $R1)"
  conda run -n nf1-mouse seqkit head -n 50000000 "$R2" -o "${LOCAL_SUB}/${LIB}/$(basename $R2)"

  # Verify
  R1_COUNT=$(conda run -n nf1-mouse seqkit stats -T "${LOCAL_SUB}/${LIB}/$(basename $R1)" | tail -1 | cut -f4)
  R2_COUNT=$(conda run -n nf1-mouse seqkit stats -T "${LOCAL_SUB}/${LIB}/$(basename $R2)" | tail -1 | cut -f4)
  echo "${LIB}: R1=${R1_COUNT} R2=${R2_COUNT}"

  # Clean up full FASTQs
  rm "${LOCAL_FULL}/${LIB}"/*_R1_001.fastq.gz "${LOCAL_FULL}/${LIB}"/*_R2_001.fastq.gz
done
```

- [ ] **Step 6: Final verification — all 6 subsampled libraries present**

```bash
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  echo -n "${LIB}: "
  ls /mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M/${LIB}/*_R1_001.fastq.gz 2>/dev/null && echo "OK" || echo "MISSING"
done
```

Expected: all 6 show "OK".

---

### Task 3: Strategy 1 — BISCUIT alignment to hg38

**Critical:** The `samtools view -bS -t ref.fai` step between `biscuit align` and `samtools sort` is mandatory. BISCUIT 1.8.0 truncates @SQ header lines for references with long contig names. This workaround injects correct headers from the .fai file.

- [ ] **Step 1: Align lib0622 to hg38 as test**

```bash
HG38_REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
HG38_FAI="${HG38_REF}.fai"
THREADS=8
LIB="lib0622"
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
```

- [ ] **Step 2: Verify lib0622 BAM header and flagstat**

```bash
conda run -n nf1-mouse samtools view -H /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy1/lib0622/lib0622.hg38.sorted.bam | head -20
cat /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy1/lib0622/lib0622.hg38.flagstat.txt
```

Expected:
- Header shows full `@SQ` lines with `SN:` and `LN:` tags (not truncated)
- Flagstat shows ~100M reads total (50M pairs), some percentage mapped
- Stop and troubleshoot if header is truncated or mapped % is unexpectedly low (<50%)

- [ ] **Step 3: Align remaining 5 libraries in parallel**

```bash
HG38_REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
HG38_FAI="${HG38_REF}.fai"
THREADS=8

for LIB in lib0623 lib0624 lib0625 lib0626 lib0627; do
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

**Note on parallelism:** 5 parallel × 8 threads × ~25 GB RAM = 40 cores + 125 GB RAM. Adjust if VM is smaller — reduce to 3 parallel if under 100 GB RAM.

- [ ] **Step 4: Verify all 6 Strategy 1 BAMs**

```bash
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  BAM="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy1/${LIB}/${LIB}.hg38.sorted.bam"
  if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
    TOTAL=$(grep "in total" "/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy1/${LIB}/${LIB}.hg38.flagstat.txt" | awk '{print $1}')
    echo "${LIB}: OK (${TOTAL} reads)"
  else
    echo "${LIB}: MISSING or not indexed"
  fi
done
```

Expected: all 6 show ~100M reads total. Do not proceed until all are present.

---

### Task 4: Strategy 3 — mm10 alignment

Strategy 3 reuses the hg38 BAMs from Strategy 1 (via symlink). This task only does the mm10 alignment.

**Critical:** Apply the same `samtools view -bS -t ref.fai` header workaround for mm10 — the bug may affect mm10's contig names too.

- [ ] **Step 1: Align lib0622 to mm10 as test**

```bash
MM10_REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa"
MM10_FAI="${MM10_REF}.fai"
THREADS=8
LIB="lib0622"
FASTQ_DIR="/mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M/${LIB}"
OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"
mkdir -p "$OUTDIR"

R1=$(ls "${FASTQ_DIR}"/*_R1_001.fastq.gz)
R2=$(ls "${FASTQ_DIR}"/*_R2_001.fastq.gz)

# Symlink Strategy 1 hg38 BAM
ln -sf "../../strategy1/${LIB}/${LIB}.hg38.sorted.bam" "${OUTDIR}/${LIB}.hg38.sorted.bam"
ln -sf "../../strategy1/${LIB}/${LIB}.hg38.sorted.bam.bai" "${OUTDIR}/${LIB}.hg38.sorted.bam.bai"

conda run -n nf1-mouse bash -c "
  biscuit align -@ ${THREADS} '${MM10_REF}' '${R1}' '${R2}' \
    | samtools view -bS -t '${MM10_FAI}' \
    | samtools sort -@ ${THREADS} -o '${OUTDIR}/${LIB}.mm10.sorted.bam'

  samtools index '${OUTDIR}/${LIB}.mm10.sorted.bam'
"
```

- [ ] **Step 2: Verify lib0622 mm10 BAM header**

```bash
conda run -n nf1-mouse samtools view -H /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/lib0622/lib0622.mm10.sorted.bam | head -10
conda run -n nf1-mouse samtools flagstat /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/lib0622/lib0622.mm10.sorted.bam
```

Expected: header shows full `@SQ` lines. Flagstat shows ~100M reads, mapped % should be substantial (mouse host reads map well to mm10).

- [ ] **Step 3: Align remaining 5 libraries to mm10 in parallel**

```bash
MM10_REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa"
MM10_FAI="${MM10_REF}.fai"
THREADS=8

for LIB in lib0623 lib0624 lib0625 lib0626 lib0627; do
(
  FASTQ_DIR="/mnt/data/projects/nf1-mouse/emseq/fastqs_sub50M/${LIB}"
  OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"
  mkdir -p "$OUTDIR"

  R1=$(ls "${FASTQ_DIR}"/*_R1_001.fastq.gz)
  R2=$(ls "${FASTQ_DIR}"/*_R2_001.fastq.gz)

  ln -sf "../../strategy1/${LIB}/${LIB}.hg38.sorted.bam" "${OUTDIR}/${LIB}.hg38.sorted.bam"
  ln -sf "../../strategy1/${LIB}/${LIB}.hg38.sorted.bam.bai" "${OUTDIR}/${LIB}.hg38.sorted.bam.bai"

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

- [ ] **Step 4: Verify all 6 mm10 BAMs and hg38 symlinks**

```bash
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"
  MM10_OK="NO"
  HG38_OK="NO"
  [ -f "${OUTDIR}/${LIB}.mm10.sorted.bam" ] && [ -f "${OUTDIR}/${LIB}.mm10.sorted.bam.bai" ] && MM10_OK="OK"
  [ -L "${OUTDIR}/${LIB}.hg38.sorted.bam" ] && [ -f "${OUTDIR}/${LIB}.hg38.sorted.bam" ] && HG38_OK="OK"
  echo "${LIB}: mm10=${MM10_OK} hg38_symlink=${HG38_OK}"
done
```

Expected: all show `mm10=OK hg38_symlink=OK`.

---

### Task 5: Disambiguate — classify reads as human/mouse/ambiguous

`disambiguate.py` name-sorts both BAMs, walks in lockstep, compares NM tags (true mismatches — BISCUIT puts conversion mismatches in ZC), assigns each read pair to the species with fewer mismatches. Ties go to ambiguous.

- [ ] **Step 1: Run disambiguate on lib0622 as test**

```bash
REPO_DIR="$HOME/repos/nf1-mouse"
LIB="lib0622"
OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"

conda run -n nf1-mouse python3 "${REPO_DIR}/dev/pdx-read-handling/disambiguate.py" \
  --human-bam "${OUTDIR}/${LIB}.hg38.sorted.bam" \
  --mouse-bam "${OUTDIR}/${LIB}.mm10.sorted.bam" \
  --out-prefix "${OUTDIR}/${LIB}" \
  --threads 8
```

- [ ] **Step 2: Verify lib0622 disambiguate output**

```bash
LIB="lib0622"
OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"

cat "${OUTDIR}/${LIB}.disambig.summary.txt"
conda run -n nf1-mouse samtools flagstat "${OUTDIR}/${LIB}.disambig.human.bam"
```

Expected:
- Summary shows total ~50M read pairs, with human, mouse, and ambiguous counts
- Human BAM exists, is coordinate-sorted, and is indexed
- Human fraction should be >0 (these are PDX with human tumor — expect 30-70% human depending on tumor burden)
- If human_assigned is 0 or near-0, stop and investigate — likely a BAM header or symlink issue

- [ ] **Step 3: Run disambiguate on remaining 5 libraries in parallel**

```bash
REPO_DIR="$HOME/repos/nf1-mouse"

for LIB in lib0623 lib0624 lib0625 lib0626 lib0627; do
(
  OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"

  conda run -n nf1-mouse python3 "${REPO_DIR}/dev/pdx-read-handling/disambiguate.py" \
    --human-bam "${OUTDIR}/${LIB}.hg38.sorted.bam" \
    --mouse-bam "${OUTDIR}/${LIB}.mm10.sorted.bam" \
    --out-prefix "${OUTDIR}/${LIB}" \
    --threads 8

  echo "Disambiguate complete: ${LIB}"
) &
done
wait
echo "All disambiguations complete."
```

- [ ] **Step 4: Generate flagstat for all Strategy 3 BAMs**

```bash
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}"
  conda run -n nf1-mouse bash -c "
    samtools flagstat '${OUTDIR}/${LIB}.hg38.sorted.bam' > '${OUTDIR}/${LIB}.hg38.flagstat.txt'
    samtools flagstat '${OUTDIR}/${LIB}.mm10.sorted.bam' > '${OUTDIR}/${LIB}.mm10.flagstat.txt'
    samtools flagstat '${OUTDIR}/${LIB}.disambig.human.bam' > '${OUTDIR}/${LIB}.disambig.human.flagstat.txt'
  "
done
```

- [ ] **Step 5: Verify all 6 disambiguate summaries**

```bash
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  SUMMARY="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/strategy3/${LIB}/${LIB}.disambig.summary.txt"
  if [ -f "$SUMMARY" ]; then
    HUMAN=$(awk -F'\t' '$1=="human_pct"{print $2}' "$SUMMARY")
    MOUSE=$(awk -F'\t' '$1=="mouse_pct"{print $2}' "$SUMMARY")
    AMBIG=$(awk -F'\t' '$1=="ambiguous_pct"{print $2}' "$SUMMARY")
    echo "${LIB}: human=${HUMAN}% mouse=${MOUSE}% ambig=${AMBIG}%"
  else
    echo "${LIB}: MISSING summary"
  fi
done
```

Expected: all 6 present with non-zero human percentage.

---

### Task 6: Comparison metrics

- [ ] **Step 1: Generate comparison TSV across both strategies**

```bash
BASEDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M"
OUTFILE="${BASEDIR}/comparison/comparison.tsv"
mkdir -p "${BASEDIR}/comparison"

echo -e "sample\tstrategy\ttotal_read_pairs\tmapped_reads\tmapq30_reads\thuman_reads\tambiguous_reads\tambiguous_pct" > "$OUTFILE"

for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  TOTAL_PAIRS=50000000

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

- [ ] **Step 2: Review comparison table**

```bash
column -t /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/comparison/comparison.tsv
```

Expected: 12 rows (6 samples × 2 strategies). Strategy 3 should show fewer human reads than Strategy 1 (disambiguate removes mouse-origin reads that Strategy 1 keeps). Ambiguous % for Strategy 3 should be ~10-15% based on prior ICRG analysis.

---

### Task 7: ichorCNA — copy number profiling

- [ ] **Step 1: Generate readCounter wig files for all BAMs**

```bash
BASEDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M"
ICHOR_DIR="${BASEDIR}/ichorcna"
mkdir -p "$ICHOR_DIR"

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

- [ ] **Step 2: Verify wig files are non-empty**

```bash
ICHOR_DIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/ichorcna"
for WIG in ${ICHOR_DIR}/*.wig; do
  LINES=$(wc -l < "$WIG")
  echo "$(basename $WIG): ${LINES} lines"
done
```

Expected: each wig file has ~2,800-3,000 lines (one header line + ~2,800 1Mb bins across autosomes+chrX).

- [ ] **Step 3: Run ichorCNA on lib0622 Strategy 1 as checkpoint**

Replace `<GC_WIG>` and `<MAP_WIG>` with the actual paths found in Task 1 Step 7.

```bash
ICHOR_DIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/ichorcna"
ICHORCNA_RSCRIPT=$(conda run -n nf1-mouse R -e 'cat(system.file("scripts", "runIchorCNA.R", package="ichorCNA"))' 2>/dev/null)
EXTDATA=$(conda run -n nf1-mouse R -e 'cat(system.file("extdata", package="ichorCNA"))' 2>/dev/null)

LIB="lib0622"
STRAT="strategy1"
WIG="${ICHOR_DIR}/${LIB}.${STRAT}.wig"
OUTDIR_ICHOR="${ICHOR_DIR}/${LIB}_${STRAT}"
mkdir -p "$OUTDIR_ICHOR"

conda run -n nf1-mouse Rscript "$ICHORCNA_RSCRIPT" \
  --id "${LIB}_${STRAT}" \
  --WIG "$WIG" \
  --gcWig "${EXTDATA}/gc_hg38_1000kb.wig" \
  --mapWig "${EXTDATA}/map_hg38_1000kb.wig" \
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
```

- [ ] **Step 4: Verify lib0622 ichorCNA output**

```bash
OUTDIR_ICHOR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/ichorcna/lib0622_strategy1"
ls -la "$OUTDIR_ICHOR"
cat "${OUTDIR_ICHOR}"/*.params.txt 2>/dev/null
```

Expected:
- Output directory contains: `.seg` file, `.params.txt`, `.pdf` (genome-wide CNA plot), `.RData`
- params.txt shows estimated tumor fraction and ploidy
- If ichorCNA crashes or produces empty output, check the log for errors. Common issues: wrong wig format, missing GC/map wig files, chromosome name mismatch (chr-prefix vs no prefix)

**Stop here and review the CNA plot before batching.** These are WU-487 (Chr8q-WT). Expect flat Chr8q. May see CDKN2A/9p deletion, SUZ12/17q loss, or other MPNST-associated events. If the plot shows all diploid with no events, that's also acceptable — it means the sample has low tumor fraction or the tumor is genuinely near-diploid.

- [ ] **Step 5: Run ichorCNA on all remaining samples and strategies**

```bash
ICHOR_DIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/ichorcna"
ICHORCNA_RSCRIPT=$(conda run -n nf1-mouse R -e 'cat(system.file("scripts", "runIchorCNA.R", package="ichorCNA"))' 2>/dev/null)
EXTDATA=$(conda run -n nf1-mouse R -e 'cat(system.file("extdata", package="ichorCNA"))' 2>/dev/null)

for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  for STRAT in strategy1 strategy3; do
    WIG="${ICHOR_DIR}/${LIB}.${STRAT}.wig"
    OUTDIR_ICHOR="${ICHOR_DIR}/${LIB}_${STRAT}"

    # Skip if already done (lib0622_strategy1 from checkpoint)
    if [ -f "${OUTDIR_ICHOR}/"*.params.txt 2>/dev/null ]; then
      echo "Skipping ${LIB}_${STRAT} (already done)"
      continue
    fi

    mkdir -p "$OUTDIR_ICHOR"

    if [ -f "$WIG" ]; then
      conda run -n nf1-mouse Rscript "$ICHORCNA_RSCRIPT" \
        --id "${LIB}_${STRAT}" \
        --WIG "$WIG" \
        --gcWig "${EXTDATA}/gc_hg38_1000kb.wig" \
        --mapWig "${EXTDATA}/map_hg38_1000kb.wig" \
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

- [ ] **Step 6: Verify all ichorCNA runs completed**

```bash
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  for STRAT in strategy1 strategy3; do
    OUTDIR_ICHOR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/ichorcna/${LIB}_${STRAT}"
    if ls "${OUTDIR_ICHOR}"/*.params.txt &>/dev/null; then
      TF=$(grep "Tumor Fraction" "${OUTDIR_ICHOR}"/*.params.txt 2>/dev/null | awk '{print $NF}')
      echo "${LIB}_${STRAT}: OK (TF=${TF})"
    else
      echo "${LIB}_${STRAT}: FAILED or MISSING"
    fi
  done
done
```

Expected: all 12 runs (6 samples × 2 strategies) show "OK" with a tumor fraction estimate.

---

### Task 8: Collect results

- [ ] **Step 1: Gather all outputs into results_summary directory**

```bash
BASEDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M"
RESULTS="${BASEDIR}/results_summary"
mkdir -p "$RESULTS"

# Comparison table
cp "${BASEDIR}/comparison/comparison.tsv" "$RESULTS/"

# Flagstat and disambig summaries
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  cp "${BASEDIR}/strategy1/${LIB}/${LIB}.hg38.flagstat.txt" "$RESULTS/${LIB}.s1.flagstat.txt"
  cp "${BASEDIR}/strategy3/${LIB}/${LIB}.disambig.summary.txt" "$RESULTS/${LIB}.s3.disambig.summary.txt" 2>/dev/null
  cp "${BASEDIR}/strategy3/${LIB}/${LIB}.disambig.human.flagstat.txt" "$RESULTS/${LIB}.s3.human.flagstat.txt" 2>/dev/null
done

# ichorCNA plots, segments, params
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

- [ ] **Step 2: Print final summary**

```bash
echo "=== COMPARISON TABLE ==="
column -t /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/results_summary/comparison.tsv
echo ""
echo "=== TUMOR FRACTIONS ==="
for LIB in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627; do
  for STRAT in strategy1 strategy3; do
    PARAMS="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/results_summary/${LIB}_${STRAT}"*.params.txt
    if [ -f $PARAMS ]; then
      TF=$(grep "Tumor Fraction" $PARAMS 2>/dev/null | awk '{print $NF}')
      PLOIDY=$(grep "Ploidy" $PARAMS 2>/dev/null | awk '{print $NF}')
      echo "${LIB}_${STRAT}: TF=${TF} Ploidy=${PLOIDY}"
    fi
  done
done
echo ""
echo "=== FILES ==="
ls /mnt/data/projects/nf1-mouse/emseq/pdx-read-handling-sub50M/results_summary/
```

Expected: comparison TSV with 12 rows, tumor fraction estimates for all 12 runs, PDF CNA plots, seg files, and params files all present in `results_summary/`.

---

## Troubleshooting Reference

| Problem | Likely cause | Fix |
|---------|-------------|-----|
| BAM header shows truncated @SQ lines | Missing `samtools view -bS -t ref.fai` step | Re-run alignment with the workaround pipe |
| disambiguate.py shows 0 human reads | Symlink broken or hg38 BAM missing | Check `ls -la` on symlinks in strategy3 dir |
| ichorCNA crashes with chromosome name error | UCSC vs Ensembl chr naming mismatch | Check if BAM uses `chr1` vs `1`; match `--genomeStyle` |
| readCounter produces empty wig | BAM not indexed or wrong chromosome names | Run `samtools index` first; check `samtools idxstats` |
| Memory pressure during parallel alignment | Too many parallel processes | Reduce from 5 to 3 parallel; each uses ~25 GB |
| ichorCNA shows TF=0 for all samples | Low tumor burden or PoN needed | Acceptable for progress report; note in results |
