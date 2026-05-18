#!/usr/bin/env bash
# Rebuild results/qc/fastqc-summary.tsv from per-lib row files.
#
# Reads canonical per-lib QC row files at:
#   ${NF1_MOUSE_DATA_ROOT}/emseq/qc/<stage>/<lib>/qc-<stage>.tsv
# where <stage> is currently just "raw-fastq" and will grow (trimmed-fastq,
# aligned-hg38, disambig-human, markdup, ...) as the pipeline extends.
#
# NF1_MOUSE_DATA_ROOT defaults to /mnt/data/projects/nf1-mouse on hosts with
# a local working store (beast, big-boy); set it to /mnt/gcs/jeszyman/projects/nf1-mouse
# to rebuild from the bucket mirror.

set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"
DATA_ROOT=${NF1_MOUSE_DATA_ROOT:-/mnt/data/projects/nf1-mouse}
QC_ROOT="${DATA_ROOT}/emseq/qc"
OUT_TSV="${REPO_ROOT}/results/qc/fastqc-summary.tsv"

HEADER="library_id	nci_alt_id	processing_stage	r1_n_reads	r2_n_reads	r1_mean_len	r2_mean_len	r1_pct_gc	r2_pct_gc	r1_mean_q	r2_mean_q	pct_dup	naive_depth_x"

mkdir -p "$(dirname "$OUT_TSV")"

if [[ ! -d "$QC_ROOT" ]]; then
  echo "ERROR: $QC_ROOT not present on this host (set NF1_MOUSE_DATA_ROOT)" >&2
  exit 1
fi

# Stages to collect, in order of pipeline progression. Add new stages as they exist.
STAGES=(raw-fastq)

{
  echo "$HEADER"
  for STAGE in "${STAGES[@]}"; do
    STAGE_DIR="${QC_ROOT}/${STAGE}"
    [[ -d "$STAGE_DIR" ]] || continue
    # cat all per-lib row files, sorted by lib_id
    find "$STAGE_DIR" -maxdepth 2 -name "qc-${STAGE}.tsv" -print0 \
      | sort -z \
      | xargs -0 cat
  done
} > "$OUT_TSV"

N=$(( $(wc -l < "$OUT_TSV") - 1 ))
echo "wrote $OUT_TSV (header + $N row files from $QC_ROOT)" >&2
