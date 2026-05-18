#########1#########2#########3#########4#########5#########6#########7#########8
#
# nf1-mouse alignment module.
#
# One rule per aligner (bwameth, biscuit, bismark). Each rule is parameterized
# by {library_id} and {ref_name}; produces a sorted+indexed BAM with paired
# flagstat + idxstats sidecars. Outputs land at:
#   <D_EMSEQ>/align/<align_method>/<ref_name>/<lib>.<ref_name>.<align_method>.sorted.bam
#
# Strategy 1 (direct-to-hg38) and Strategy 3 (dual-align hg38 + mm10, then
# disambiguate downstream) are *not* separate modes here — they're just
# different rule_all targets that consume the same per-aligner-per-ref BAMs.
#
# Required wrapper-provided variables:
#   D_EMSEQ, D_REF, D_LOGS, D_BENCHMARK, D_DATA, CONDA_NF1_MOUSE
#
#########1#########2#########3#########4#########5#########6#########7#########8


# ── bwa-meth alignment ─────────────────────────────────────────────────────
rule align_bwameth:
    message: "align: bwameth {wildcards.library_id} → {wildcards.ref_name}"
    conda: CONDA_NF1_MOUSE
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R2.fastq.gz",
        fa = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa",
        c2t = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}_align_bwameth.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}_align_bwameth.tsv"
    params:
        sort_tmp = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{wc.ref_name}.bwameth",
    threads: 8
    output:
        bam = f"{D_EMSEQ}/align/bwameth/{{ref_name}}/{{library_id}}.{{ref_name}}.bwameth.sorted.bam",
        bai = f"{D_EMSEQ}/align/bwameth/{{ref_name}}/{{library_id}}.{{ref_name}}.bwameth.sorted.bam.bai",
        flag = f"{D_EMSEQ}/align/bwameth/{{ref_name}}/{{library_id}}.{{ref_name}}.bwameth.flagstat.txt",
        idx  = f"{D_EMSEQ}/align/bwameth/{{ref_name}}/{{library_id}}.{{ref_name}}.bwameth.idxstats.txt",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")" "$(dirname "{output.bam}")" "$(dirname "{params.sort_tmp}")"
        exec &>> "{log.cmd}"
        echo "[align_bwameth] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"
        bwameth.py --threads {threads} --reference "{input.fa}" "{input.r1}" "{input.r2}" \
          | samtools view -@ {threads} -u -F 4 - \
          | samtools sort -@ {threads} -T "{params.sort_tmp}" -o "{output.bam}"
        samtools index "{output.bam}"
        samtools flagstat "{output.bam}" > "{output.flag}"
        samtools idxstats "{output.bam}" > "{output.idx}"
        """


# ── BISCUIT alignment ──────────────────────────────────────────────────────
rule align_biscuit:
    message: "align: biscuit {wildcards.library_id} → {wildcards.ref_name}"
    conda: CONDA_NF1_MOUSE
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R2.fastq.gz",
        fa = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa",
        fai = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.fai",
        bis_amb = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.bis.amb",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}_align_biscuit.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}_align_biscuit.tsv"
    params:
        sort_tmp = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{wc.ref_name}.biscuit",
    threads: 8
    output:
        bam = f"{D_EMSEQ}/align/biscuit/{{ref_name}}/{{library_id}}.{{ref_name}}.biscuit.sorted.bam",
        bai = f"{D_EMSEQ}/align/biscuit/{{ref_name}}/{{library_id}}.{{ref_name}}.biscuit.sorted.bam.bai",
        flag = f"{D_EMSEQ}/align/biscuit/{{ref_name}}/{{library_id}}.{{ref_name}}.biscuit.flagstat.txt",
        idx  = f"{D_EMSEQ}/align/biscuit/{{ref_name}}/{{library_id}}.{{ref_name}}.biscuit.idxstats.txt",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")" "$(dirname "{output.bam}")" "$(dirname "{params.sort_tmp}")"
        exec &>> "{log.cmd}"
        echo "[align_biscuit] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"
        # samtools view -t REF.fai overrides BISCUIT's occasionally-truncated @SQ header (1.8.0 bug).
        biscuit align -@ {threads} "{input.fa}" "{input.r1}" "{input.r2}" \
          | samtools view -uS -t "{input.fai}" - \
          | samtools sort -@ {threads} -T "{params.sort_tmp}" -o "{output.bam}"
        samtools index "{output.bam}"
        samtools flagstat "{output.bam}" > "{output.flag}"
        samtools idxstats "{output.bam}" > "{output.idx}"
        """


# ── Bismark alignment ──────────────────────────────────────────────────────
rule align_bismark:
    message: "align: bismark {wildcards.library_id} → {wildcards.ref_name}"
    conda: CONDA_NF1_MOUSE
    input:
        r1 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_EMSEQ}/fastqs/{{library_id}}/{{library_id}}.raw_R2.fastq.gz",
        fa = f"{D_REF}/bismark/{{ref_name}}/{{ref_name}}.fa",
        ct = f"{D_REF}/bismark/{{ref_name}}/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2",
    log:
        cmd = f"{D_LOGS}/{{library_id}}.{{ref_name}}_align_bismark.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}.{{ref_name}}_align_bismark.tsv"
    params:
        ref_dir = lambda wc: f"{D_REF}/bismark/{wc.ref_name}",
        out_dir = lambda wc: f"{D_EMSEQ}/align/bismark/{wc.ref_name}",
        tmp_dir = lambda wc: f"{D_DATA}/tmp/{wc.library_id}.{wc.ref_name}.bismark",
        # Bismark names its output <R1_basename_no_.fastq.gz>_bismark_bt2_pe.bam
        raw_bam = lambda wc: f"{D_EMSEQ}/align/bismark/{wc.ref_name}/{wc.library_id}.raw_R1_bismark_bt2_pe.bam",
    threads: 8
    output:
        bam = f"{D_EMSEQ}/align/bismark/{{ref_name}}/{{library_id}}.{{ref_name}}.bismark.sorted.bam",
        bai = f"{D_EMSEQ}/align/bismark/{{ref_name}}/{{library_id}}.{{ref_name}}.bismark.sorted.bam.bai",
        flag = f"{D_EMSEQ}/align/bismark/{{ref_name}}/{{library_id}}.{{ref_name}}.bismark.flagstat.txt",
        idx  = f"{D_EMSEQ}/align/bismark/{{ref_name}}/{{library_id}}.{{ref_name}}.bismark.idxstats.txt",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")" "{params.out_dir}" "{params.tmp_dir}"
        exec &>> "{log.cmd}"
        echo "[align_bismark] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"
        bismark --genome "{params.ref_dir}" \
          -1 "{input.r1}" -2 "{input.r2}" \
          -o "{params.out_dir}" \
          --temp_dir "{params.tmp_dir}" \
          --parallel {threads}
        samtools sort -@ {threads} -T "{params.tmp_dir}/sort" -o "{output.bam}" "{params.raw_bam}"
        rm -f "{params.raw_bam}"
        samtools index "{output.bam}"
        samtools flagstat "{output.bam}" > "{output.flag}"
        samtools idxstats "{output.bam}" > "{output.idx}"
        """
