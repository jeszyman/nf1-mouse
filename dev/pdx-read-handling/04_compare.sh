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
# Compare strategies across all samples
# Collect: total read pairs (from FASTQ), mapped, MAPQ>=30, human-assigned, ambiguous rate

BASEDIR="/mnt/data/projects/nf1-mouse/emseq/pdx-read-handling"
OUTFILE="${BASEDIR}/comparison.tsv"

echo -e "sample\tstrategy\ttotal_read_pairs\tmapped_reads\tmapq30_reads\thuman_reads\tambiguous_reads\tambiguous_pct" > "$OUTFILE"

for SAMPLE in lib0622 lib0623 lib0624 lib0625 lib0626 lib0627 \
              lib0628 lib0629 lib0630 lib0631 lib0632 lib0633 \
              lib0634 lib0635 lib0636 lib0637 lib0638 lib0639; do

  # Count total read pairs from FASTQ (shared denominator for both strategies)
  FASTQ_DIR="/mnt/data/projects/nf1-mouse/emseq/fastqs/${SAMPLE}"
  if [ -d "$FASTQ_DIR" ]; then
    TOTAL_PAIRS=$(zcat "${FASTQ_DIR}"/*_R1_001.fastq.gz | awk 'NR%4==1' | wc -l)
  else
    continue
  fi

  # Strategy 1
  S1_BAM="$BASEDIR/strategy1/${SAMPLE}/${SAMPLE}.hg38.sorted.bam"
  if [ -f "$S1_BAM" ]; then
    MAPPED=$(samtools view -c -F 4 "$S1_BAM")
    MAPQ30=$(samtools view -c -F 4 -q 30 "$S1_BAM")
    echo -e "${SAMPLE}\tstrategy1_direct\t${TOTAL_PAIRS}\t${MAPPED}\t${MAPQ30}\t${MAPQ30}\t0\t0.0" >> "$OUTFILE"
  fi

  # Strategy 3
  S3_SUMMARY="$BASEDIR/strategy3/${SAMPLE}/${SAMPLE}.disambig.summary.txt"
  S3_BAM="$BASEDIR/strategy3/${SAMPLE}/${SAMPLE}.disambig.human.bam"
  S3_HG38="$BASEDIR/strategy3/${SAMPLE}/${SAMPLE}.hg38.sorted.bam"
  if [ -f "$S3_SUMMARY" ]; then
    MAPPED=$(samtools view -c -F 4 "$S3_HG38")
    HUMAN=$(awk -F'\t' '$1=="human_assigned"{print $2}' "$S3_SUMMARY")
    AMBIG=$(awk -F'\t' '$1=="ambiguous"{print $2}' "$S3_SUMMARY")
    AMBIG_PCT=$(awk -F'\t' '$1=="ambiguous_pct"{print $2}' "$S3_SUMMARY")
    MAPQ30=$(samtools view -c -F 4 -q 30 "$S3_BAM")
    echo -e "${SAMPLE}\tstrategy3_disambig\t${TOTAL_PAIRS}\t${MAPPED}\t${MAPQ30}\t${HUMAN}\t${AMBIG}\t${AMBIG_PCT}" >> "$OUTFILE"
  fi
done

echo "Comparison written to $OUTFILE"
cat "$OUTFILE" | column -t
