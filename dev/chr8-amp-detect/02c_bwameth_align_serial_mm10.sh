#!/usr/bin/env bash
set -euo pipefail
# Phase C-mm10: serial bwa-meth alignment of the same 450M-subsampled FASTQs
# to mm10. Pairs with 02b (hg38 alignment) for downstream PDX disambiguation.
# Inputs: chr8-amp-detect/600M/subsampled/<lib>_R{1,2}.450M.fq.gz
# Outputs: chr8-amp-detect/600M/bams_mm10/<lib>.mm10.dedup.bam (+ .bai, flagstat, idxstats)
TAG_FILE=/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/logs/00_target_tag.txt
[[ -s "$TAG_FILE" ]] || { echo "ERROR: $TAG_FILE missing" >&2; exit 1; }
TARGET_TAG=$(cat "$TAG_FILE")

REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa"
BASE="/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M"
SUB="${BASE}/subsampled"
OUTDIR="${BASE}/bams_mm10"
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
  SORTED="${OUTDIR}/${LIB}.mm10.sorted.bam"
  DEDUP="${OUTDIR}/${LIB}.mm10.dedup.bam"
  LOG="${LOGDIR}/${LIB}.mm10.align.log"

  if [[ -s "$DEDUP" && -s "${DEDUP}.bai" ]]; then
    echo "[$(date)] ${LIB}: mm10 dedup BAM exists, skipping"
    continue
  fi
  [[ -s "$R1" && -s "$R2" ]] || { echo "[$(date)] ${LIB}: subsampled FASTQs missing" >&2; continue; }

  {
    echo "[$(date)] ${LIB}: mm10 bwa-meth align start (${THREADS_ALIGN} threads, tag=${TARGET_TAG})"
    bwameth.py --threads "$THREADS_ALIGN" --reference "$REF" "$R1" "$R2" \
      | samtools fixmate -m -u -@ "$THREADS_IO" - - \
      | samtools sort -@ "$THREADS_IO" -m 4G -o "$SORTED" -
    samtools index -@ "$THREADS_IO" "$SORTED"

    echo "[$(date)] ${LIB}: markdup"
    samtools markdup -r -@ "$THREADS_IO" "$SORTED" "$DEDUP"
    samtools index -@ "$THREADS_IO" "$DEDUP"

    samtools flagstat -@ "$THREADS_IO" "$DEDUP" > "${OUTDIR}/${LIB}.mm10.dedup.flagstat.txt"
    samtools idxstats "$DEDUP" > "${OUTDIR}/${LIB}.mm10.dedup.idxstats.txt"

    if [[ -s "$DEDUP" && -s "${DEDUP}.bai" ]]; then
      rm -f "$SORTED" "${SORTED}.bai"
      echo "[$(date)] ${LIB}: removed sorted.bam"
    fi
    echo "[$(date)] ${LIB}: done"
  } >> "$LOG" 2>&1
done
echo "[$(date)] all libs done"
