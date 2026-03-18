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
# Setup references on local disk at /mnt/data/projects/nf1-mouse/ref/
# Env: emseq (bwameth.py + biscuit)
# hg38 refs already symlinked from breast project; this script handles mm10 + ICRG
REFDIR="/mnt/data/projects/nf1-mouse/ref"

# -- mm10: download Ensembl GRCm38 primary assembly if not present
MM10_DIR="$REFDIR/biscuit/mm10"
mkdir -p "$MM10_DIR"
if [ ! -f "$MM10_DIR/mm10.fa" ]; then
  echo "Downloading mm10 (Ensembl GRCm38 release 102)..."
  wget -q -O "$MM10_DIR/mm10.fa.gz" \
    "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
  gunzip "$MM10_DIR/mm10.fa.gz"
fi

# -- Index mm10 for BISCUIT (Strategy 3)
if [ ! -f "$MM10_DIR/mm10.fa.bis.amb" ]; then
  echo "Indexing mm10 for BISCUIT..."
  biscuit index "$MM10_DIR/mm10.fa"
fi
samtools faidx "$MM10_DIR/mm10.fa" 2>/dev/null || true

# -- ICRG combined reference: hg38 + mm10 with mouse chroms prefixed mm10_ (Strategy 2)
ICRG_DIR="$REFDIR/bwa_meth/icrg_hg38_mm10"
HG38_FA="$REFDIR/bwa_meth/ncbi_decoy_hg38/ncbi_decoy_hg38.fa"
mkdir -p "$ICRG_DIR"
if [ ! -f "$ICRG_DIR/icrg_hg38_mm10.fa" ]; then
  echo "Building ICRG combined reference (hg38 + mm10_ prefixed mouse)..."
  cp "$HG38_FA" "$ICRG_DIR/icrg_hg38_mm10.fa"
  sed 's/^>/>mm10_/' "$MM10_DIR/mm10.fa" >> "$ICRG_DIR/icrg_hg38_mm10.fa"
fi
if [ ! -f "$ICRG_DIR/icrg_hg38_mm10.fa.bwameth.c2t" ]; then
  echo "Indexing ICRG combined reference for bwa-meth..."
  bwameth.py index "$ICRG_DIR/icrg_hg38_mm10.fa"
fi
samtools faidx "$ICRG_DIR/icrg_hg38_mm10.fa" 2>/dev/null || true

echo "Reference setup complete."
echo "  hg38 bwa-meth: $REFDIR/bwa_meth/ncbi_decoy_hg38/"
echo "  hg38 biscuit:  $REFDIR/biscuit_hg38/"
echo "  mm10 biscuit:  $MM10_DIR/"
echo "  ICRG bwa-meth: $ICRG_DIR/"
