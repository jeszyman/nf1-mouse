#!/usr/bin/env bash
set -euo pipefail
# Phase 4: HMMcopy readCounter -> 1Mb wig, then patched ichorCNA (no-PoN,
# default-low-input, autosomes only) on a single dedup BAM.
# Usage: sudo -n -u jupyter bash 03_ichorcna.sh <LIB>
#   LIB    library id, e.g. lib0626 (dedup BAM must already exist)
#
# Why "sudo -n -u jupyter": every input under /mnt/data/projects/nf1-mouse/
# is owned by the jupyter group with drwxrws--- (no other access). The
# patched ichorCNA package + R deps live in R_libs/ which is jupyter-only
# readable as well. Run the whole script as jupyter so plain file tests
# and writes work uniformly.

LIB=$1

BAMDIR="/mnt/data/projects/nf1-mouse/chr8-amp-detect/bams"
OUTBASE="/mnt/data/projects/nf1-mouse/chr8-amp-detect/ichorcna"
WIGDIR="${OUTBASE}/wig"
OUTDIR="${OUTBASE}/${LIB}_nopon_default-low-input"

READCOUNTER="/opt/conda/envs/hmmcopy-bin/bin/readCounter"
SAMTOOLS="/opt/conda/envs/nf1-mouse/bin/samtools"
R_LIBS_PATH="/mnt/data/projects/nf1-mouse/R_libs"
# runIchorCNA.R is not shipped inside the installed R package — it lives in
# the patched fork repo's scripts/ dir. The jupyter user can read this path
# because the repo dir is g+rx world-traversable.
RUNICHOR="/home/ext_szymanski_jeffrey_mayo_edu/repos/ichorCNA-patched/scripts/runIchorCNA.R"
EXTDATA="${R_LIBS_PATH}/ichorCNA/extdata"

BAM="${BAMDIR}/${LIB}.hg38.dedup.bam"
WIG="${WIGDIR}/${LIB}.1Mb.wig"

if [[ ! -s "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM" >&2
  exit 1
fi

mkdir -p "$WIGDIR" "$OUTDIR"

# Detect autosome contig names from the BAM header (UCSC 'chr1' vs NCBI '1')
CHRS=$("$SAMTOOLS" view -H "$BAM" \
  | awk -F'\t' '/^@SQ/{for(i=1;i<=NF;i++) if($i~/^SN:/){sub("SN:","",$i); print $i}}' \
  | grep -E '^(chr)?([1-9]|1[0-9]|2[0-2])$' \
  | paste -sd ',')
echo "[$(date)] ${LIB}: autosome contigs -> ${CHRS}"

# 1) Bin-count to 1Mb wig (autosomes only, q>=20)
if [[ ! -s "$WIG" ]]; then
  echo "[$(date)] ${LIB}: readCounter -> ${WIG}"
  "$READCOUNTER" \
    --window 1000000 \
    --quality 20 \
    --chromosome "$CHRS" \
    "$BAM" > "$WIG"
else
  echo "[$(date)] ${LIB}: wig already exists, skipping readCounter"
fi

# 2) ichorCNA (patched fork) with no-PoN, default-low-input, autosomes only
# Build the chrs R-expression from CHRS (e.g. "chr1,chr2,..." -> "c('chr1','chr2',...)")
CHRS_R="c('$(echo "$CHRS" | sed "s/,/','/g")')"

echo "[$(date)] ${LIB}: ichorCNA -> ${OUTDIR}"
R_LIBS="$R_LIBS_PATH" /usr/bin/Rscript "$RUNICHOR" \
  --id "${LIB}_nopon_dli" \
  --WIG "$WIG" \
  --gcWig "${EXTDATA}/gc_hg38_1000kb.wig" \
  --mapWig "${EXTDATA}/map_hg38_1000kb.wig" \
  --centromere "${EXTDATA}/GRCh38.GCA_000001405.2_centromere_acen.txt" \
  --outDir "$OUTDIR" \
  --includeHOMD FALSE \
  --chrs "$CHRS_R" \
  --chrTrain "$CHRS_R" \
  --estimateNormal TRUE \
  --estimatePloidy TRUE \
  --estimateScPrevalence FALSE \
  --txnE 0.9999 \
  --txnStrength 10000 \
  --genomeBuild "hg38" \
  --genomeStyle "UCSC" \
  2>&1 | tee "${OUTDIR}/ichorCNA.log"

echo "[$(date)] ${LIB}: done -> ${OUTDIR}"
