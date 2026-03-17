#!/usr/bin/env python3
# Conda env: biotools (has pysam 0.23.3)
# Run: conda run -n biotools python3 disambiguate.py --human-bam ... --mouse-bam ... --out-prefix ...
"""Disambiguate PDX reads by comparing NM tags from dual BISCUIT alignments.

For each read pair present in both BAMs, compare the sum of NM tags
(true mismatches only — BISCUIT separates conversion mismatches into ZC).
Assign to the species with fewer mismatches. Ties go to ambiguous.

Reads unmapped in both species go to ambiguous.
Reads mapped in only one species go to that species.
"""

import argparse
import pysam
from collections import defaultdict


def get_nm_by_qname(bam_path):
    """Extract NM tag sum per read pair (read name -> total NM)."""
    nm = defaultdict(int)
    mapped = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            mapped.add(read.query_name)
            try:
                nm[read.query_name] += read.get_tag("NM")
            except KeyError:
                nm[read.query_name] += 0
    return nm, mapped


def main():
    parser = argparse.ArgumentParser(description="Disambiguate PDX reads")
    parser.add_argument("--human-bam", required=True)
    parser.add_argument("--mouse-bam", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    print("Reading human BAM NM tags...")
    human_nm, human_mapped = get_nm_by_qname(args.human_bam)
    print(f"  {len(human_mapped)} mapped read pairs")

    print("Reading mouse BAM NM tags...")
    mouse_nm, mouse_mapped = get_nm_by_qname(args.mouse_bam)
    print(f"  {len(mouse_mapped)} mapped read pairs")

    all_reads = human_mapped | mouse_mapped
    human_only = set()
    mouse_only = set()
    ambiguous = set()

    for qname in all_reads:
        in_human = qname in human_mapped
        in_mouse = qname in mouse_mapped

        if in_human and not in_mouse:
            human_only.add(qname)
        elif in_mouse and not in_human:
            mouse_only.add(qname)
        else:
            h_nm = human_nm[qname]
            m_nm = mouse_nm[qname]
            if h_nm < m_nm:
                human_only.add(qname)
            elif m_nm < h_nm:
                mouse_only.add(qname)
            else:
                ambiguous.add(qname)

    print(f"Classification: {len(human_only)} human, {len(mouse_only)} mouse, {len(ambiguous)} ambiguous")

    # Write human-assigned reads from the human BAM
    out_human = f"{args.out_prefix}.disambig.human.bam"
    out_mouse = f"{args.out_prefix}.disambig.mouse.bam"
    out_ambig = f"{args.out_prefix}.disambig.ambiguous.bam"

    with pysam.AlignmentFile(args.human_bam, "rb") as inbam:
        with pysam.AlignmentFile(out_human, "wb", header=inbam.header) as outbam:
            for read in inbam:
                if read.query_name in human_only:
                    outbam.write(read)

    with pysam.AlignmentFile(args.mouse_bam, "rb") as inbam:
        with pysam.AlignmentFile(out_mouse, "wb", header=inbam.header) as outbam:
            for read in inbam:
                if read.query_name in mouse_only:
                    outbam.write(read)

    with pysam.AlignmentFile(args.human_bam, "rb") as inbam:
        with pysam.AlignmentFile(out_ambig, "wb", header=inbam.header) as outbam:
            for read in inbam:
                if read.query_name in ambiguous:
                    outbam.write(read)

    pysam.sort("-@", str(args.threads), "-o", out_human + ".tmp", out_human)
    import os
    os.rename(out_human + ".tmp", out_human)
    pysam.index(out_human)

    # Write summary
    with open(f"{args.out_prefix}.disambig.summary.txt", "w") as f:
        f.write(f"total_read_pairs\t{len(all_reads)}\n")
        f.write(f"human_assigned\t{len(human_only)}\n")
        f.write(f"mouse_assigned\t{len(mouse_only)}\n")
        f.write(f"ambiguous\t{len(ambiguous)}\n")
        f.write(f"human_pct\t{100*len(human_only)/len(all_reads):.1f}\n")
        f.write(f"mouse_pct\t{100*len(mouse_only)/len(all_reads):.1f}\n")
        f.write(f"ambiguous_pct\t{100*len(ambiguous)/len(all_reads):.1f}\n")

    print(f"Wrote {out_human}, {out_mouse}, {out_ambig}")


if __name__ == "__main__":
    main()
