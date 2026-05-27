#!/usr/bin/env bash
# Phase C: Serial bwa-meth alignment of equal-depth-subsampled libraries
# at 96 threads per sample. Reuses recipe from 02_bwameth_align.sh but:
#   - reads inputs from chr8-amp-detect/${TAG}/subsampled/
#   - writes BAMs to chr8-amp-detect/${TAG}/bams/
#   - serial (avoids prior 4-way VM crash)
#   - deletes sorted.bam after dedup verified, to keep peak disk under 1.5 TB
# Env: nf1-mouse conda env (bwameth + samtools)
# Run as: sudo -n -u jupyter nohup bash 02b_bwameth_align_serial.sh > /mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/logs/00_driver.log 2>&1 &
set -euo pipefail

TAG_FILE=/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/logs/00_target_tag.txt
if [[ ! -s "$TAG_FILE" ]]; then
  echo "ERROR: $TAG_FILE missing. Run 01b_subsample_equal.sh first." >&2
  exit 1
fi
TARGET_TAG=$(cat "$TAG_FILE")

REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
BASE="/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M"
SUB="${BASE}/subsampled"
OUTDIR="${BASE}/bams"
LOGDIR="${BASE}/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

THREADS_ALIGN=96
THREADS_IO=16

source /opt/conda/etc/profile.d/conda.sh
conda activate nf1-mouse

LIBS=(lib0626 lib0627 lib0644 lib0645)

for LIB in "${LIBS[@]}"; do
  R1="${SUB}/${LIB}_R1.${TARGET_TAG}.fq.gz"
  R2="${SUB}/${LIB}_R2.${TARGET_TAG}.fq.gz"
  SORTED="${OUTDIR}/${LIB}.hg38.sorted.bam"
  DEDUP="${OUTDIR}/${LIB}.hg38.dedup.bam"
  LOG="${LOGDIR}/${LIB}.align.log"

  if [[ -s "$DEDUP" && -s "${DEDUP}.bai" ]]; then
    echo "[$(date)] ${LIB}: dedup BAM exists, skipping"
    continue
  fi

  if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "[$(date)] ${LIB}: subsampled FASTQs missing ($R1) — skipping" | tee -a "$LOG"
    continue
  fi

  {
    echo "[$(date)] ${LIB}: bwa-meth align start (${THREADS_ALIGN} threads, tag=${TARGET_TAG})"
    bwameth.py --threads "$THREADS_ALIGN" --reference "$REF" "$R1" "$R2" \
      | samtools fixmate -m -u -@ "$THREADS_IO" - - \
      | samtools sort -@ "$THREADS_IO" -m 4G -o "$SORTED" -
    samtools index -@ "$THREADS_IO" "$SORTED"

    echo "[$(date)] ${LIB}: markdup"
    samtools markdup -r -@ "$THREADS_IO" "$SORTED" "$DEDUP"
    samtools index -@ "$THREADS_IO" "$DEDUP"

    samtools flagstat -@ "$THREADS_IO" "$DEDUP" > "${OUTDIR}/${LIB}.dedup.flagstat.txt"
    samtools idxstats "$DEDUP" > "${OUTDIR}/${LIB}.dedup.idxstats.txt"
    samtools flagstat -@ "$THREADS_IO" "$SORTED" > "${OUTDIR}/${LIB}.sorted.flagstat.txt"

    # Free disk: drop sorted.bam now that dedup is verified
    if [[ -s "$DEDUP" && -s "${DEDUP}.bai" ]]; then
      rm -f "$SORTED" "${SORTED}.bai"
      echo "[$(date)] ${LIB}: removed sorted.bam"
    fi

    echo "[$(date)] ${LIB}: done"
  } >> "$LOG" 2>&1
done

echo "[$(date)] all libs done"
