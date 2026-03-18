#!/usr/bin/env bash
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/nf1-mouse/nf1-mouse.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-18 06:45:51
# ============================================================

set -euo pipefail
# Strategy 2: ICRG combined hg38+mm10 reference
# Env: emseq (biscuit + samtools)

SAMPLE=$1
THREADS=${2:-8}
REF="/mnt/data/projects/nf1-mouse/ref/biscuit/icrg_hg38_mm10/icrg_hg38_mm10.fa"
FASTQ_DIR="/mnt/data/projects/nf1-mouse/emseq/fastqs/${SAMPLE}"
OUTDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/strategy2/${SAMPLE}"
mkdir -p "$OUTDIR"

R1=$(ls "${FASTQ_DIR}"/*_R1_001.fastq.gz)
R2=$(ls "${FASTQ_DIR}"/*_R2_001.fastq.gz)

# Align to combined reference
biscuit align -@ "$THREADS" "$REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "${OUTDIR}/${SAMPLE}.icrg.sorted.bam"

samtools index "${OUTDIR}/${SAMPLE}.icrg.sorted.bam"

# Extract human chromosomes (those without mm10_ prefix)
samtools view -@ "$THREADS" -b "${OUTDIR}/${SAMPLE}.icrg.sorted.bam" \
  $(samtools idxstats "${OUTDIR}/${SAMPLE}.icrg.sorted.bam" | cut -f1 | grep -v "^mm10_" | grep -v "^\*" | tr '\n' ' ') \
  > "${OUTDIR}/${SAMPLE}.icrg.human_only.sorted.bam"

samtools index "${OUTDIR}/${SAMPLE}.icrg.human_only.sorted.bam"
samtools flagstat "${OUTDIR}/${SAMPLE}.icrg.sorted.bam" > "${OUTDIR}/${SAMPLE}.icrg.flagstat.txt"
samtools flagstat "${OUTDIR}/${SAMPLE}.icrg.human_only.sorted.bam" > "${OUTDIR}/${SAMPLE}.icrg.human_only.flagstat.txt"
samtools idxstats "${OUTDIR}/${SAMPLE}.icrg.sorted.bam" > "${OUTDIR}/${SAMPLE}.icrg.idxstats.txt"
