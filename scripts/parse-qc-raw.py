#!/usr/bin/env python3
"""Parse raw-FASTQ QC outputs for a single library; write a per-lib TSV row file.

Combines FastQC (per-read metrics) and fastp (pair-level duplication).

Usage:
    parse-qc-raw.py <lib_id> <qc_dir> [<nci_alt_id>]
    e.g. parse-qc-raw.py lib0622 /mnt/data/.../qc/raw-fastq/lib0622/ nf1_lib95

If <nci_alt_id> is omitted, the parser tries to recover it from the FASTQ
filename inside the FastQC zip (regex `nf1_lib\\d+`). After the
data_link_raw_fastq rule renames inputs to <lib>.raw_R{1,2}.fastq.gz, the
caller (snakemake rule) must pass nci_alt_id explicitly.

Writes <qc_dir>/qc-raw-fastq.tsv — one row, no header. This is the canonical
per-lib row; rebuild-summary.sh concatenates all such row files into the
repo-committed results/qc/fastqc-summary.tsv (Option B: per-lib row file +
concat view).

Expects in <qc_dir>:
    *_R1*_fastqc.zip
    *_R2*_fastqc.zip
    fastp.json

Schema (13 columns, long-format on processing_stage):
    library_id  nci_alt_id  processing_stage  r1_n_reads  r2_n_reads
    r1_mean_len  r2_mean_len  r1_pct_gc  r2_pct_gc  r1_mean_q  r2_mean_q
    pct_dup  naive_depth_x
"""
from __future__ import annotations
import json
import re
import sys
import zipfile
from pathlib import Path


HEADER = [
    "library_id", "nci_alt_id", "processing_stage",
    "r1_n_reads", "r2_n_reads",
    "r1_mean_len", "r2_mean_len",
    "r1_pct_gc", "r2_pct_gc",
    "r1_mean_q", "r2_mean_q",
    "pct_dup",
    "naive_depth_x",
]

HG38_BP = 3_100_000_000  # for naive depth (assumes 100% on-target)
PROCESSING_STAGE = "raw-fastq"


def parse_fastqc_zip(zpath: Path) -> dict:
    """Extract per-read metrics from a single FastQC zip."""
    with zipfile.ZipFile(zpath) as zf:
        data_name = next(n for n in zf.namelist() if n.endswith("/fastqc_data.txt"))
        text = zf.read(data_name).decode()

    basic = {}
    in_basic = False
    for line in text.splitlines():
        if line.startswith(">>Basic Statistics"):
            in_basic = True
            continue
        if in_basic and line.startswith(">>END_MODULE"):
            break
        if in_basic and "\t" in line and not line.startswith("#"):
            k, v = line.split("\t", 1)
            basic[k] = v

    filename = basic.get("Filename", "")
    m = re.search(r"(nf1_lib\d+)", filename)
    nci = m.group(1) if m else ""

    n_reads = int(basic.get("Total Sequences", 0))
    pct_gc = float(basic.get("%GC", 0))
    seq_len = basic.get("Sequence length", "0")
    if "-" in seq_len:
        lo, hi = (int(x) for x in seq_len.split("-"))
        mean_len = (lo + hi) / 2.0
    else:
        mean_len = float(seq_len)

    # Per-base sequence quality: mean across positions
    quals = []
    in_qual = False
    for line in text.splitlines():
        if line.startswith(">>Per base sequence quality"):
            in_qual = True
            continue
        if in_qual and line.startswith(">>END_MODULE"):
            break
        if in_qual and "\t" in line and not line.startswith("#"):
            cols = line.split("\t")
            try:
                quals.append(float(cols[1]))
            except (ValueError, IndexError):
                pass
    mean_q = round(sum(quals) / len(quals), 2) if quals else float("nan")

    return {
        "nci": nci,
        "n_reads": n_reads,
        "mean_len": mean_len,
        "pct_gc": pct_gc,
        "mean_q": mean_q,
    }


def parse_fastp_json(jpath: Path) -> float:
    """Return pair-level duplication rate as a percentage (0-100)."""
    with open(jpath) as fh:
        d = json.load(fh)
    # fastp duplication.rate is a fraction [0, 1]
    return round(100.0 * d["duplication"]["rate"], 2)


def main(argv: list[str]) -> int:
    if len(argv) not in (3, 4):
        print("usage: parse-qc-raw.py <lib_id> <qc_dir> [<nci_alt_id>]", file=sys.stderr)
        return 2
    lib_id, qc_dir = argv[1], Path(argv[2])
    nci_arg = argv[3] if len(argv) == 4 else ""

    zips = sorted(qc_dir.glob("*_fastqc.zip"))
    if len(zips) != 2:
        print(f"ERROR: expected 2 fastqc.zip in {qc_dir}, found {len(zips)}", file=sys.stderr)
        return 1
    r1_zip = next(z for z in zips if "_R1" in z.name)
    r2_zip = next(z for z in zips if "_R2" in z.name)
    r1 = parse_fastqc_zip(r1_zip)
    r2 = parse_fastqc_zip(r2_zip)

    fastp_json = qc_dir / "fastp.json"
    if not fastp_json.exists():
        print(f"ERROR: missing fastp.json in {qc_dir}", file=sys.stderr)
        return 1
    pct_dup = parse_fastp_json(fastp_json)

    total_bases = r1["n_reads"] * r1["mean_len"] + r2["n_reads"] * r2["mean_len"]
    naive_depth_x = round(total_bases / HG38_BP, 2)

    nci_alt_id = nci_arg or r1["nci"]

    row = [
        lib_id, nci_alt_id, PROCESSING_STAGE,
        r1["n_reads"], r2["n_reads"],
        r1["mean_len"], r2["mean_len"],
        r1["pct_gc"], r2["pct_gc"],
        r1["mean_q"], r2["mean_q"],
        pct_dup,
        naive_depth_x,
    ]
    out_path = qc_dir / f"qc-{PROCESSING_STAGE}.tsv"
    out_path.write_text("\t".join(str(x) for x in row) + "\n")
    print(f"wrote {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
