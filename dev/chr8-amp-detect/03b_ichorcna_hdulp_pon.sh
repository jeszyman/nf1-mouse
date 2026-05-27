#!/usr/bin/env bash
set -euo pipefail
# Phase 4 (with PoN): same as 03_ichorcna.sh but adds --normalPanel pointing
# at ichorCNA's bundled HD-ULP hg38 1Mb panel of normals. Output goes to a
# separate *_hdulp-pon_* directory so the no-PoN runs are preserved.
# Usage: sudo -n -u jupyter bash 03b_ichorcna_hdulp_pon.sh <LIB>

LIB=$1

BAMDIR="/mnt/data/projects/nf1-mouse/chr8-amp-detect/bams"
OUTBASE="/mnt/data/projects/nf1-mouse/chr8-amp-detect/ichorcna"
WIGDIR="${OUTBASE}/wig"
OUTDIR="${OUTBASE}/${LIB}_hdulp-pon_default-low-input"

READCOUNTER="/opt/conda/envs/hmmcopy-bin/bin/readCounter"
SAMTOOLS="/opt/conda/envs/nf1-mouse/bin/samtools"
R_LIBS_PATH="/mnt/data/projects/nf1-mouse/R_libs"
RUNICHOR="/home/ext_szymanski_jeffrey_mayo_edu/repos/ichorCNA-patched/scripts/runIchorCNA.R"
EXTDATA="${R_LIBS_PATH}/ichorCNA/extdata"
PON="${EXTDATA}/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds"

BAM="${BAMDIR}/${LIB}.hg38.dedup.bam"
WIG="${WIGDIR}/${LIB}.1Mb.wig"

if [[ ! -s "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM" >&2
  exit 1
fi
if [[ ! -s "$PON" ]]; then
  echo "ERROR: PoN not found: $PON" >&2
  exit 1
fi

mkdir -p "$WIGDIR" "$OUTDIR"

CHRS=$("$SAMTOOLS" view -H "$BAM" \
  | awk -F'\t' '/^@SQ/{for(i=1;i<=NF;i++) if($i~/^SN:/){sub("SN:","",$i); print $i}}' \
  | grep -E '^(chr)?([1-9]|1[0-9]|2[0-2])$' \
  | paste -sd ',')
echo "[$(date)] ${LIB}: autosome contigs -> ${CHRS}"

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

CHRS_R="c('$(echo "$CHRS" | sed "s/,/','/g")')"

echo "[$(date)] ${LIB}: ichorCNA (+HD-ULP PoN) -> ${OUTDIR}"
R_LIBS="$R_LIBS_PATH" /usr/bin/Rscript "$RUNICHOR" \
  --id "${LIB}_hdulp-pon_dli" \
  --WIG "$WIG" \
  --gcWig "${EXTDATA}/gc_hg38_1000kb.wig" \
  --mapWig "${EXTDATA}/map_hg38_1000kb.wig" \
  --centromere "${EXTDATA}/GRCh38.GCA_000001405.2_centromere_acen.txt" \
  --normalPanel "$PON" \
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
