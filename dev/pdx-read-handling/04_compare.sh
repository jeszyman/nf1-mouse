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
# Compare strategies across both samples
# Collect: total reads, mapped reads, MAPQ>=30, human-assigned reads

BASEDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling"
OUTFILE="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling/comparison.tsv"

echo -e "sample\tstrategy\ttotal_reads\tmapped_reads\tmapq30_reads\thuman_reads" > "$OUTFILE"

for SAMPLE in lib0625 lib0626; do
  # Strategy 1
  S1_BAM="$BASEDIR/strategy1/${SAMPLE}/${SAMPLE}.hg38.sorted.bam"
  if [ -f "$S1_BAM" ]; then
    TOTAL=$(samtools view -c "$S1_BAM")
    MAPPED=$(samtools view -c -F 4 "$S1_BAM")
    MAPQ30=$(samtools view -c -F 4 -q 30 "$S1_BAM")
    echo -e "${SAMPLE}\tstrategy1_direct\t${TOTAL}\t${MAPPED}\t${MAPQ30}\t${MAPQ30}" >> "$OUTFILE"
  fi

  # Strategy 2
  S2_BAM="$BASEDIR/strategy2/${SAMPLE}/${SAMPLE}.icrg.human_only.sorted.bam"
  S2_FULL="$BASEDIR/strategy2/${SAMPLE}/${SAMPLE}.icrg.sorted.bam"
  if [ -f "$S2_BAM" ]; then
    TOTAL=$(samtools view -c "$S2_FULL")
    MAPPED=$(samtools view -c -F 4 "$S2_FULL")
    HUMAN=$(samtools view -c "$S2_BAM")
    MAPQ30=$(samtools view -c -F 4 -q 30 "$S2_BAM")
    echo -e "${SAMPLE}\tstrategy2_icrg\t${TOTAL}\t${MAPPED}\t${MAPQ30}\t${HUMAN}" >> "$OUTFILE"
  fi

  # Strategy 3
  S3_BAM="$BASEDIR/strategy3/${SAMPLE}/${SAMPLE}.disambig.human.bam"
  S3_HG38="$BASEDIR/strategy3/${SAMPLE}/${SAMPLE}.hg38.sorted.bam"
  if [ -f "$S3_BAM" ]; then
    TOTAL=$(samtools view -c "$S3_HG38")
    MAPPED=$(samtools view -c -F 4 "$S3_HG38")
    HUMAN=$(samtools view -c "$S3_BAM")
    MAPQ30=$(samtools view -c -F 4 -q 30 "$S3_BAM")
    echo -e "${SAMPLE}\tstrategy3_disambig\t${TOTAL}\t${MAPPED}\t${MAPQ30}\t${HUMAN}" >> "$OUTFILE"
  fi
done

echo "Comparison written to $OUTFILE"
cat "$OUTFILE" | column -t
