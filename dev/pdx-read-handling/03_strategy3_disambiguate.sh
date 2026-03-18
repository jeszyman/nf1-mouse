#!/usr/bin/env bash
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/nf1-mouse/nf1-mouse.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-18 08:24:29
# ============================================================

set -euo pipefail
# Strategy 3: Disambiguate via dual BISCUIT alignment
# Env: emseq (biscuit + samtools)
# BISCUIT NM tag = true mismatches only (conversion in ZC tag)

SAMPLE=$1
THREADS=${2:-8}
HG38_REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
MM10_REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ucsc_mm10/mm10.fa"
FASTQ_DIR="/mnt/data/projects/nf1-mouse/emseq/fastqs/${SAMPLE}"
OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/strategy3/${SAMPLE}"
mkdir -p "$OUTDIR"

R1=$(ls "${FASTQ_DIR}"/*_R1_001.fastq.gz)
R2=$(ls "${FASTQ_DIR}"/*_R2_001.fastq.gz)

# Align to human
biscuit align -@ "$THREADS" "$HG38_REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "${OUTDIR}/${SAMPLE}.hg38.sorted.bam"
samtools index "${OUTDIR}/${SAMPLE}.hg38.sorted.bam"

# Align to mouse
biscuit align -@ "$THREADS" "$MM10_REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "${OUTDIR}/${SAMPLE}.mm10.sorted.bam"
samtools index "${OUTDIR}/${SAMPLE}.mm10.sorted.bam"

# Disambiguate: compare NM tags, assign to species with fewer mismatches
conda run -n biotools python3 /home/jeszyman/repos/nf1-mouse/dev/pdx-read-handling/disambiguate.py \
  --human-bam "${OUTDIR}/${SAMPLE}.hg38.sorted.bam" \
  --mouse-bam "${OUTDIR}/${SAMPLE}.mm10.sorted.bam" \
  --out-prefix "${OUTDIR}/${SAMPLE}" \
  --threads "$THREADS"

samtools flagstat "${OUTDIR}/${SAMPLE}.hg38.sorted.bam" > "${OUTDIR}/${SAMPLE}.hg38.flagstat.txt"
samtools flagstat "${OUTDIR}/${SAMPLE}.mm10.sorted.bam" > "${OUTDIR}/${SAMPLE}.mm10.flagstat.txt"
samtools flagstat "${OUTDIR}/${SAMPLE}.disambig.human.bam" > "${OUTDIR}/${SAMPLE}.disambig.human.flagstat.txt"
