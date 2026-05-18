#########1#########2#########3#########4#########5#########6#########7#########8
#
# nf1-mouse QC module.
#
# Raw-FASTQ QC pass: FastQC per read + fastp pair-level stats + per-lib
# parsed row + repo-committed summary view. Long-format on processing_stage
# (raw-fastq) — future stages (trimmed-fastq, aligned-hg38, …) follow the
# same pattern.
#
# Required wrapper-provided variables:
#   D_EMSEQ, D_LOGS, D_BENCHMARK, R_NF1_MOUSE, CONDA_NF1_MOUSE, LIBRARIES
#
#########1#########2#########3#########4#########5#########6#########7#########8


# ── Raw-FASTQ FastQC ───────────────────────────────────────────────────────
rule qc_fastqc:
    message: "qc: FastQC raw FASTQ for {wildcards.library_id}"
    conda: CONDA_NF1_MOUSE
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R2.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_qc_fastqc.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_qc_fastqc.tsv"
    threads: 4
    output:
        r1_zip = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/{{library_id}}.raw_R1_fastqc.zip",
        r1_html = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/{{library_id}}.raw_R1_fastqc.html",
        r2_zip = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/{{library_id}}.raw_R2_fastqc.zip",
        r2_html = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/{{library_id}}.raw_R2_fastqc.html",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")" "$(dirname "{output.r1_zip}")"
        exec &>> "{log.cmd}"
        echo "[qc_fastqc] $(date) lib={wildcards.library_id} threads={threads}"
        fastqc --quiet --threads {threads} \
          --outdir "$(dirname "{output.r1_zip}")" \
          "{input.r1}" "{input.r2}"
        """


# ── Raw-FASTQ fastp (stats only) ───────────────────────────────────────────
# Omitting --out1/--out2 puts fastp in stats-only mode: only JSON + HTML
# reports are written, no filtered FASTQs. Pair-level duplication.rate is the
# headline metric we consume downstream.
rule qc_fastp:
    message: "qc: fastp stats-only on raw FASTQ for {wildcards.library_id}"
    conda: CONDA_NF1_MOUSE
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R2.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_qc_fastp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_qc_fastp.tsv"
    threads: 4
    output:
        json = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/fastp.json",
        html = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/fastp.html",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")" "$(dirname "{output.json}")"
        exec &>> "{log.cmd}"
        echo "[qc_fastp] $(date) lib={wildcards.library_id} threads={threads}"
        fastp \
          --in1 "{input.r1}" --in2 "{input.r2}" \
          --json "{output.json}" --html "{output.html}" \
          --disable_adapter_trimming \
          --disable_quality_filtering \
          --disable_length_filtering \
          --thread {threads}
        """


# ── Per-lib row file (canonical QC record at processing_stage=raw-fastq) ───
rule qc_parse_raw:
    message: "qc: parse raw-FASTQ QC outputs → per-lib row file for {wildcards.library_id}"
    conda: CONDA_NF1_MOUSE
    input:
        r1_zip = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/{{library_id}}.raw_R1_fastqc.zip",
        r2_zip = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/{{library_id}}.raw_R2_fastqc.zip",
        fastp_json = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/fastp.json",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_qc_parse_raw.log",
    params:
        qc_dir = lambda wc: f"{D_EMSEQ}/qc/raw-fastq/{wc.library_id}",
        script = f"{R_NF1_MOUSE}/scripts/parse-qc-raw.py",
        nci_alt_id = lambda wc: NCI_ALT_IDS[wc.library_id],
    output:
        row = f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/qc-raw-fastq.tsv",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")"
        exec &>> "{log.cmd}"
        echo "[qc_parse_raw] $(date) lib={wildcards.library_id} nci={params.nci_alt_id}"
        python3 "{params.script}" \
          "{wildcards.library_id}" "{params.qc_dir}" "{params.nci_alt_id}"
        """


# ── Repo summary view (aggregates per-lib rows across all stages) ──────────
rule qc_rebuild_summary:
    message: "qc: rebuild repo summary view from per-lib row files"
    conda: CONDA_NF1_MOUSE
    input:
        expand(f"{D_EMSEQ}/qc/raw-fastq/{{library_id}}/qc-raw-fastq.tsv",
               library_id=LIBRARIES),
    log:
        cmd = f"{D_LOGS}/qc_rebuild_summary.log",
    params:
        script = f"{R_NF1_MOUSE}/scripts/rebuild-summary.sh",
        data_root = D_DATA,
    output:
        tsv = f"{R_NF1_MOUSE}/results/qc/fastqc-summary.tsv",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")" "$(dirname "{output.tsv}")"
        exec &>> "{log.cmd}"
        echo "[qc_rebuild_summary] $(date) n_libs={input}" | head -c 200
        echo
        NF1_MOUSE_DATA_ROOT="{params.data_root}" bash "{params.script}"
        """
