#!/usr/bin/env bash
# Master driver for the chr8-amp-detect 600M re-run.
# Chains:
#   1. wait for lib0626/0627 gsutil rsync FASTQs to settle
#   2. run 01b_subsample_equal.sh (counts + equal-depth seqtk)
#   3. run 02b_bwameth_align_serial.sh (serial 96-thread alignment, all 4 libs)
#   4. for each lib, run 03b_ichorcna_600M.sh
# Each phase is idempotent / resumable per the underlying scripts.
# Phase 5 (plots) + 6 (org write-up) are not run here — kicked off
# manually after a human looks at ichorCNA output.
# Usage: sudo -n -u jupyter setsid nohup bash 00_master_600M.sh \
#          > /mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/logs/00_master.log 2>&1 &
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

LOG=/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/logs
mkdir -p "$LOG"

# 1) Wait for lib0626 + lib0627 FASTQ pairs to exist with non-zero size.
echo "[$(date)] master: waiting for lib0626/0627 fastqs to land"
while :; do
  F1=/mnt/data/projects/nf1-mouse/chr8-amp-detect/fastqs/lib0626/99_nf1_lib99_S5_R1_001.fastq.gz
  F2=/mnt/data/projects/nf1-mouse/chr8-amp-detect/fastqs/lib0626/99_nf1_lib99_S5_R2_001.fastq.gz
  F3=/mnt/data/projects/nf1-mouse/chr8-amp-detect/fastqs/lib0627/100_nf1_lib100_S7_R1_001.fastq.gz
  F4=/mnt/data/projects/nf1-mouse/chr8-amp-detect/fastqs/lib0627/100_nf1_lib100_S7_R2_001.fastq.gz
  # All 4 must exist, be > 1GB (sanity), AND the rsync process not running.
  if [[ -s "$F1" && -s "$F2" && -s "$F3" && -s "$F4" \
        && $(stat -c %s "$F1") -gt 1073741824 \
        && $(stat -c %s "$F3") -gt 1073741824 ]] \
     && ! pgrep -f "gsutil.*Sample_99_nf1_lib99\|gsutil.*Sample_100_nf1_lib100" >/dev/null; then
    echo "[$(date)] master: lib0626/0627 fastqs settled"
    break
  fi
  sleep 60
done

# 2) Equal-depth subsample
echo "[$(date)] master: phase B (subsample)"
bash "${SCRIPTDIR}/01b_subsample_equal.sh"

# 3) Serial bwa-meth alignment
echo "[$(date)] master: phase C (alignment)"
bash "${SCRIPTDIR}/02b_bwameth_align_serial.sh"

# 4) ichorCNA per lib
echo "[$(date)] master: phase D (ichorCNA)"
for LIB in lib0626 lib0627 lib0644 lib0645; do
  bash "${SCRIPTDIR}/03b_ichorcna_600M.sh" "$LIB"
done

echo "[$(date)] master: ALL DONE. Ready for phase E (plots) + phase F (org write-up)."
