#!/usr/bin/env bash
# Equal-depth subsample for the 4 chr8-amp-detect libraries.
# - Counts R1 reads in each raw FASTQ (parallel)
# - Picks TARGET = floor(min(counts) / 50_000_000) * 50_000_000
# - seqtk subsamples R1 and R2 of each lib with same -s42 seed
# Env: nf1-mouse conda env (provides seqtk)
# Run as: sudo -n -u jupyter bash 01b_subsample_equal.sh
set -euo pipefail

BASE="/mnt/data/projects/nf1-mouse/chr8-amp-detect"
OUTDIR_ROOT="${BASE}/600M"      # actual target tag is decided dynamically below
LOGDIR="${OUTDIR_ROOT}/logs"
mkdir -p "$LOGDIR"

# Source conda + activate env so seqtk is on PATH
source /opt/conda/etc/profile.d/conda.sh
conda activate nf1-mouse

declare -A R1=(
  [lib0626]="${BASE}/fastqs/lib0626/$(basename $(ls ${BASE}/fastqs/lib0626/*_R1_001.fastq.gz 2>/dev/null | head -1))"
  [lib0627]="${BASE}/fastqs/lib0627/$(basename $(ls ${BASE}/fastqs/lib0627/*_R1_001.fastq.gz 2>/dev/null | head -1))"
  [lib0644]="${BASE}/fastqs/lib0644/5_MN2_4_5_2_12_26_S5_R1_001.fastq.gz"
  [lib0645]="${BASE}/fastqs/lib0645/6_MN2_3_4_2_12_26_S6_R1_001.fastq.gz"
)
declare -A R2
for k in "${!R1[@]}"; do
  R2[$k]="${R1[$k]/_R1_/_R2_}"
done

LIBS=(lib0626 lib0627 lib0644 lib0645)

# ---- Step 1: count R1 reads per lib, in parallel ----
echo "[$(date)] counting R1 reads (parallel)"
declare -A COUNT_PID
for L in "${LIBS[@]}"; do
  ( zcat "${R1[$L]}" | awk 'NR%4==1' | wc -l > "${LOGDIR}/${L}.count.txt"; \
    echo "[$(date)] ${L}: $(cat ${LOGDIR}/${L}.count.txt) reads" ) &
  COUNT_PID[$L]=$!
done
for L in "${LIBS[@]}"; do wait "${COUNT_PID[$L]}"; done

# ---- Step 2: pick TARGET ----
MIN=$(for L in "${LIBS[@]}"; do cat "${LOGDIR}/${L}.count.txt"; done | sort -n | head -1)
TARGET=$(( (MIN / 50000000) * 50000000 ))
TARGET_M=$(( TARGET / 1000000 ))
TARGET_TAG="${TARGET_M}M"

echo "[$(date)] MIN=${MIN}  TARGET=${TARGET}  TAG=${TARGET_TAG}"
for L in "${LIBS[@]}"; do
  printf "  %s: %s\n" "$L" "$(cat ${LOGDIR}/${L}.count.txt)"
done | tee -a "${LOGDIR}/00_target.txt"
echo "${TARGET}" > "${LOGDIR}/00_target_value.txt"
echo "${TARGET_TAG}" > "${LOGDIR}/00_target_tag.txt"

SUBDIR="${OUTDIR_ROOT}/subsampled"
mkdir -p "$SUBDIR"

# ---- Step 3: seqtk subsample, STREAMING Bernoulli (fraction < 1) ----
# seqtk sample <FRAC> (float < 1) uses streaming Bernoulli sampling: each read
# kept with probability FRAC, constant memory. Same seed on R1/R2 = identical
# decisions = synced pairs. Reservoir mode (seqtk sample <INT>) OOM'd at 450M
# (~135 GB RAM/process) and crashed the 220-GB VM twice on 2026-05-23.
#
# We oversample by 0.5% so all libs end up with >= TARGET reads. The exact
# per-lib output count varies by ~sqrt(N) reads — fine for ichorCNA bin coverage.
for L in "${LIBS[@]}"; do
  O1="${SUBDIR}/${L}_R1.${TARGET_TAG}.fq.gz"
  O2="${SUBDIR}/${L}_R2.${TARGET_TAG}.fq.gz"
  if [[ -s "$O1" && -s "$O2" ]]; then
    echo "[$(date)] ${L}: subsample outputs exist, skipping"
    continue
  fi
  COUNT=$(cat "${LOGDIR}/${L}.count.txt")
  # FRAC = min(0.999, TARGET * 1.005 / COUNT), 4-decimal precision.
  FRAC=$(awk -v t="$TARGET" -v c="$COUNT" 'BEGIN{f=t*1.005/c; if(f>0.999)f=0.999; printf "%.4f", f}')
  echo "[$(date)] ${L}: seqtk sample -s42 ${FRAC}  (streaming, COUNT=${COUNT}, TARGET~${TARGET})"
  # R1 + R2 in parallel — streaming mode uses negligible RAM so this is safe.
  seqtk sample -s42 "${R1[$L]}" "$FRAC" | gzip > "$O1" &
  P1=$!
  seqtk sample -s42 "${R2[$L]}" "$FRAC" | gzip > "$O2" &
  P2=$!
  wait $P1 $P2
  echo "[$(date)] ${L}: done"
done

# ---- Step 4: verify ----
echo "[$(date)] verification"
for L in "${LIBS[@]}"; do
  O1="${SUBDIR}/${L}_R1.${TARGET_TAG}.fq.gz"
  O2="${SUBDIR}/${L}_R2.${TARGET_TAG}.fq.gz"
  C1=$(zcat "$O1" | awk 'NR%4==1' | wc -l)
  C2=$(zcat "$O2" | awk 'NR%4==1' | wc -l)
  echo "  ${L}: R1=${C1}  R2=${C2}  (target=${TARGET})"
done

echo "[$(date)] all done. TARGET_TAG=${TARGET_TAG}"
