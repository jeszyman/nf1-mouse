#!/usr/bin/env bash
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/nf1-mouse/nf1-mouse.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-24 11:06:46
# ============================================================

set -euo pipefail
# Strategy 1: Direct BISCUIT alignment to hg38
# Env: emseq (biscuit + samtools)

SAMPLE=$1
THREADS=${2:-8}
REF="/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
FASTQ_DIR="/mnt/data/projects/nf1-mouse/emseq/fastqs/${SAMPLE}"
OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/strategy1/${SAMPLE}"
mkdir -p "$OUTDIR"

R1=$(ls "${FASTQ_DIR}"/*_R1_001.fastq.gz)
R2=$(ls "${FASTQ_DIR}"/*_R2_001.fastq.gz)

biscuit align -@ "$THREADS" "$REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "${OUTDIR}/${SAMPLE}.hg38.sorted.bam"

samtools index "${OUTDIR}/${SAMPLE}.hg38.sorted.bam"
samtools flagstat "${OUTDIR}/${SAMPLE}.hg38.sorted.bam" > "${OUTDIR}/${SAMPLE}.hg38.flagstat.txt"
samtools idxstats "${OUTDIR}/${SAMPLE}.hg38.sorted.bam" > "${OUTDIR}/${SAMPLE}.hg38.idxstats.txt"
