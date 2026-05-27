#!/usr/bin/env bash
# Phase D: thin wrapper around 03_ichorcna.sh that points at the 600M BAMs
# and writes outputs under chr8-amp-detect/600M/ichorcna/.
# Usage: sudo -n -u jupyter bash 03b_ichorcna_600M.sh <LIB>
set -euo pipefail

LIB=$1
BAMDIR="/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/bams"
OUTBASE="/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/ichorcna"
WIGDIR="${OUTBASE}/wig"
OUTDIR="${OUTBASE}/${LIB}_nopon_default-low-input"

READCOUNTER="/opt/conda/envs/hmmcopy-bin/bin/readCounter"
SAMTOOLS="/opt/conda/envs/nf1-mouse/bin/samtools"
R_LIBS_PATH="/mnt/data/projects/nf1-mouse/R_libs"
RUNICHOR="/home/ext_szymanski_jeffrey_mayo_edu/repos/ichorCNA-patched/scripts/runIchorCNA.R"
EXTDATA="${R_LIBS_PATH}/ichorCNA/extdata"

BAM="${BAMDIR}/${LIB}.hg38.dedup.bam"
WIG="${WIGDIR}/${LIB}.1Mb.wig"

if [[ ! -s "$BAM" ]]; then echo "ERROR: BAM not found: $BAM" >&2; exit 1; fi
mkdir -p "$WIGDIR" "$OUTDIR"

CHRS=$("$SAMTOOLS" view -H "$BAM" \
  | awk -F'\t' '/^@SQ/{for(i=1;i<=NF;i++) if($i~/^SN:/){sub("SN:","",$i); print $i}}' \
  | grep -E '^(chr)?([1-9]|1[0-9]|2[0-2])$' | paste -sd ',')
echo "[$(date)] ${LIB}: autosomes -> ${CHRS}"

if [[ ! -s "$WIG" ]]; then
  "$READCOUNTER" --window 1000000 --quality 20 --chromosome "$CHRS" "$BAM" > "$WIG"
fi

CHRS_R="c('$(echo "$CHRS" | sed "s/,/','/g")')"
R_LIBS="$R_LIBS_PATH" /usr/bin/Rscript "$RUNICHOR" \
  --id "${LIB}_nopon_dli_600M" --WIG "$WIG" \
  --gcWig "${EXTDATA}/gc_hg38_1000kb.wig" \
  --mapWig "${EXTDATA}/map_hg38_1000kb.wig" \
  --centromere "${EXTDATA}/GRCh38.GCA_000001405.2_centromere_acen.txt" \
  --outDir "$OUTDIR" --includeHOMD FALSE \
  --chrs "$CHRS_R" --chrTrain "$CHRS_R" \
  --estimateNormal TRUE --estimatePloidy TRUE --estimateScPrevalence FALSE \
  --txnE 0.9999 --txnStrength 10000 \
  --genomeBuild "hg38" --genomeStyle "UCSC" \
  2>&1 | tee "${OUTDIR}/ichorCNA.log"
echo "[$(date)] ${LIB}: ichorCNA done -> ${OUTDIR}"
