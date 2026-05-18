#########1#########2#########3#########4#########5#########6#########7#########8
#
# nf1-mouse data setup module.
#
# Stage 0 rules: bridge immutable raw delivery (under D_INPUTS/fastqs/) to
# lib-named symlinks under D_EMSEQ/fastqs/ that downstream rules consume.
# Required wrapper-provided variables:
#   D_INPUTS, D_EMSEQ, D_LOGS, CONDA_NF1_MOUSE
#
#########1#########2#########3#########4#########5#########6#########7#########8

from glob import glob


def _raw_fastq(library_id, read):
    """Resolve the immutable raw FASTQ path for one (library_id, read)."""
    hits = sorted(glob(f"{D_INPUTS}/fastqs/{library_id}/*_R{read}_001.fastq.gz"))
    if not hits:
        raise FileNotFoundError(
            f"data_link_raw_fastq: no R{read} FASTQ for {library_id} "
            f"under {D_INPUTS}/fastqs/{library_id}/"
        )
    if len(hits) > 1:
        raise ValueError(
            f"data_link_raw_fastq: ambiguous R{read} FASTQ for {library_id}: {hits}"
        )
    return hits[0]


# ── Raw FASTQ symlinks ─────────────────────────────────────────────────────
# Decouple raw delivery filenames (e.g. 95_nf1_lib95_S55_R1_001.fastq.gz) from
# downstream rule wildcards. Convention: <lib>.raw_R{1,2}.fastq.gz are
# symlinks to the immutable raw; *.processed_R{1,2}.fastq.gz are real files
# produced by downstream rules (trim/filter/etc.).
rule data_link_raw_fastq:
    message: "data: link raw FASTQ pair for {wildcards.library_id}"
    conda: CONDA_NF1_MOUSE
    input:
        r1 = lambda wc: _raw_fastq(wc.library_id, 1),
        r2 = lambda wc: _raw_fastq(wc.library_id, 2),
    log:
        cmd = f"{D_LOGS}/{{library_id}}_data_link_raw_fastq.log",
    output:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R2.fastq.gz",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")"
        exec &>> "{log.cmd}"
        echo "[data_link_raw_fastq] $(date) lib={wildcards.library_id}"
        mkdir -p "$(dirname "{output.r1}")"
        ln -srf "{input.r1}" "{output.r1}"
        ln -srf "{input.r2}" "{output.r2}"
        ls -l "{output.r1}" "{output.r2}"
        """
