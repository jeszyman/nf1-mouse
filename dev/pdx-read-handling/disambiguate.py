#!/usr/bin/env python3
# Conda env: biotools (has pysam 0.23.3)
# Run: conda run -n biotools python3 disambiguate.py --human-bam ... --mouse-bam ... --out-prefix ...
"""Disambiguate PDX reads by comparing NM tags from dual BISCUIT alignments.

Streaming approach: name-sort both BAMs, walk in lockstep, compare NM per
read pair, write output immediately. RAM is O(1) regardless of input size.

For each read pair present in both BAMs, compare the sum of NM tags
(true mismatches only — BISCUIT separates conversion mismatches into ZC).
Assign to the species with fewer mismatches. Ties go to ambiguous.

Reads mapped in only one species go to that species.
"""

import argparse
import os
import subprocess
import sys
from itertools import groupby

import pysam


def name_sort_bam(bam_path, threads):
    """Name-sort a BAM file, return path to sorted file."""
    out = bam_path.replace(".sorted.bam", ".namesorted.bam")
    if out == bam_path:
        out = bam_path + ".namesorted.bam"
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(threads), "-o", out, bam_path],
        check=True,
    )
    return out


def nm_sum_for_group(reads):
    """Sum NM tags across all alignments in a read-pair group. Returns (nm_sum, read_list, has_mapped)."""
    nm = 0
    mapped = False
    read_list = []
    for r in reads:
        read_list.append(r)
        if not r.is_unmapped:
            mapped = True
            try:
                nm += r.get_tag("NM")
            except KeyError:
                pass
    return nm, read_list, mapped


def main():
    parser = argparse.ArgumentParser(description="Disambiguate PDX reads (streaming)")
    parser.add_argument("--human-bam", required=True, help="Coordinate-sorted human BAM")
    parser.add_argument("--mouse-bam", required=True, help="Coordinate-sorted mouse BAM")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    # Name-sort both BAMs
    print("Name-sorting human BAM...")
    human_ns = name_sort_bam(args.human_bam, args.threads)
    print("Name-sorting mouse BAM...")
    mouse_ns = name_sort_bam(args.mouse_bam, args.threads)

    # Open name-sorted BAMs and output files
    h_bam = pysam.AlignmentFile(human_ns, "rb", check_sq=False)
    m_bam = pysam.AlignmentFile(mouse_ns, "rb", check_sq=False)

    out_human_path = f"{args.out_prefix}.disambig.human.bam"
    out_mouse_path = f"{args.out_prefix}.disambig.mouse.bam"
    out_ambig_path = f"{args.out_prefix}.disambig.ambiguous.bam"

    out_human = pysam.AlignmentFile(out_human_path, "wb", header=h_bam.header)
    out_mouse = pysam.AlignmentFile(out_mouse_path, "wb", header=m_bam.header)
    out_ambig = pysam.AlignmentFile(out_ambig_path, "wb", header=h_bam.header)

    # Group reads by query name in each BAM
    h_groups = groupby(h_bam, key=lambda r: r.query_name)
    m_groups = groupby(m_bam, key=lambda r: r.query_name)

    counts = {"human": 0, "mouse": 0, "ambiguous": 0, "total": 0}

    # Merge-walk: advance whichever iterator has the smaller query name
    h_name, h_reads = None, None
    m_name, m_reads = None, None

    def next_h():
        try:
            name, grp = next(h_groups)
            return name, list(grp)
        except StopIteration:
            return None, None

    def next_m():
        try:
            name, grp = next(m_groups)
            return name, list(grp)
        except StopIteration:
            return None, None

    h_name, h_reads = next_h()
    m_name, m_reads = next_m()

    while h_name is not None or m_name is not None:
        if h_name is not None and (m_name is None or h_name < m_name):
            # Human only
            h_nm, h_list, h_mapped = nm_sum_for_group(iter(h_reads))
            if h_mapped:
                counts["human"] += 1
                for r in h_list:
                    out_human.write(r)
            counts["total"] += 1
            h_name, h_reads = next_h()

        elif m_name is not None and (h_name is None or m_name < h_name):
            # Mouse only
            m_nm, m_list, m_mapped = nm_sum_for_group(iter(m_reads))
            if m_mapped:
                counts["mouse"] += 1
                for r in m_list:
                    out_mouse.write(r)
            counts["total"] += 1
            m_name, m_reads = next_m()

        else:
            # Same read name in both — compare NM
            h_nm, h_list, h_mapped = nm_sum_for_group(iter(h_reads))
            m_nm, m_list, m_mapped = nm_sum_for_group(iter(m_reads))
            counts["total"] += 1

            if h_mapped and not m_mapped:
                counts["human"] += 1
                for r in h_list:
                    out_human.write(r)
            elif m_mapped and not h_mapped:
                counts["mouse"] += 1
                for r in m_list:
                    out_mouse.write(r)
            elif h_mapped and m_mapped:
                if h_nm < m_nm:
                    counts["human"] += 1
                    for r in h_list:
                        out_human.write(r)
                elif m_nm < h_nm:
                    counts["mouse"] += 1
                    for r in m_list:
                        out_mouse.write(r)
                else:
                    counts["ambiguous"] += 1
                    for r in h_list:
                        out_ambig.write(r)
            else:
                counts["ambiguous"] += 1

            h_name, h_reads = next_h()
            m_name, m_reads = next_m()

        if counts["total"] % 5_000_000 == 0 and counts["total"] > 0:
            print(f"  Processed {counts['total']:,} read pairs...", file=sys.stderr)

    out_human.close()
    out_mouse.close()
    out_ambig.close()
    h_bam.close()
    m_bam.close()

    print(f"Classification: {counts['human']:,} human, {counts['mouse']:,} mouse, {counts['ambiguous']:,} ambiguous")

    # Sort and index human output
    print("Coordinate-sorting human output...")
    pysam.sort("-@", str(args.threads), "-o", out_human_path + ".tmp", out_human_path)
    os.rename(out_human_path + ".tmp", out_human_path)
    pysam.index(out_human_path)

    # Clean up name-sorted intermediates
    os.remove(human_ns)
    os.remove(mouse_ns)

    # Write summary
    with open(f"{args.out_prefix}.disambig.summary.txt", "w") as f:
        f.write(f"total_read_pairs\t{counts['total']}\n")
        f.write(f"human_assigned\t{counts['human']}\n")
        f.write(f"mouse_assigned\t{counts['mouse']}\n")
        f.write(f"ambiguous\t{counts['ambiguous']}\n")
        if counts["total"] > 0:
            f.write(f"human_pct\t{100*counts['human']/counts['total']:.1f}\n")
            f.write(f"mouse_pct\t{100*counts['mouse']/counts['total']:.1f}\n")
            f.write(f"ambiguous_pct\t{100*counts['ambiguous']/counts['total']:.1f}\n")

    print(f"Wrote {out_human_path}, {out_mouse_path}, {out_ambig_path}")


if __name__ == "__main__":
    main()
