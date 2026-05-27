#!/usr/bin/env bash
set -euo pipefail
# Resume bwa-meth alignments serially after the May 18 VM crash.
# Runs ONE library at a time with 32 align threads; ~90 min each on this VM.
# Env: nf1-mouse

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALIGN="${SCRIPT_DIR}/02_bwameth_align.sh"

LIBS=(lib0626 lib0627 lib0644 lib0645)
THREADS=32

LOGDIR="/mnt/data/projects/nf1-mouse/chr8-amp-detect/logs"
mkdir -p "$LOGDIR"

source /opt/conda/etc/profile.d/conda.sh
conda activate nf1-mouse

for LIB in "${LIBS[@]}"; do
  LOG="${LOGDIR}/align_${LIB}.resume.log"
  echo "[$(date)] === starting ${LIB} (serial, ${THREADS} threads) ===" | tee -a "$LOG"
  bash "$ALIGN" "$LIB" "$THREADS" >>"$LOG" 2>&1
  echo "[$(date)] === finished ${LIB} ===" | tee -a "$LOG"
done

echo "[$(date)] all libraries done"
