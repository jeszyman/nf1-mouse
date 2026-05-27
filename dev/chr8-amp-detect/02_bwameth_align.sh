#!/usr/bin/env bash
set -euo pipefail
# Phase 3: bwa-meth human-only direct alignment of one subsampled sample.
# Env: nf1-mouse (bwameth.py + samtools)
# Usage: 02_bwameth_align.sh <LIB> [THREADS]
#   LIB     library id, e.g. lib0626 (subsampled fastqs must already exist)
#   THREADS bwa-meth alignment threads (default 32)

LIB=$1
THREADS=${2:-32}
SORT_THREADS=12

REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
BASE="/mnt/data/projects/nf1-mouse/chr8-amp-detect"
SUB="${BASE}/subsampled"
OUTDIR="${BASE}/bams"
mkdir -p "$OUTDIR"

R1="${SUB}/${LIB}_R1.50M.fq.gz"
R2="${SUB}/${LIB}_R2.50M.fq.gz"
SORTED="${OUTDIR}/${LIB}.hg38.sorted.bam"
DEDUP="${OUTDIR}/${LIB}.hg38.dedup.bam"

if [[ -s "$DEDUP" && -s "${DEDUP}.bai" ]]; then
  echo "[$(date)] ${LIB}: dedup BAM already exists, skipping"
  exit 0
fi

echo "[$(date)] ${LIB}: bwa-meth align start (${THREADS} threads)"

# bwa-meth emits mate-adjacent records; fixmate -m tags them for markdup,
# then coordinate-sort to the QC sorted.bam.
bwameth.py --threads "$THREADS" --reference "$REF" "$R1" "$R2" \
  | samtools fixmate -m -u -@ "$SORT_THREADS" - - \
  | samtools sort -@ "$SORT_THREADS" -o "$SORTED" -
samtools index -@ "$SORT_THREADS" "$SORTED"

echo "[$(date)] ${LIB}: markdup"
samtools markdup -r -@ "$SORT_THREADS" "$SORTED" "$DEDUP"
samtools index -@ "$SORT_THREADS" "$DEDUP"

samtools flagstat -@ "$SORT_THREADS" "$DEDUP" > "${OUTDIR}/${LIB}.dedup.flagstat.txt"
samtools idxstats "$DEDUP" > "${OUTDIR}/${LIB}.dedup.idxstats.txt"
samtools flagstat -@ "$SORT_THREADS" "$SORTED" > "${OUTDIR}/${LIB}.sorted.flagstat.txt"

echo "[$(date)] ${LIB}: done"
