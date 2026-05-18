#########1#########2#########3#########4#########5#########6#########7#########8
#
# nf1-mouse reference module.
#
# Builds per-aligner indexes for each named reference. FASTA convention:
#   <D_REF>/<align_method>/<ref_name>/<ref_name>.fa
# Each aligner's index files live alongside the FASTA in that subdir.
#
# Source FASTAs: ref_stage_fasta symlinks the source-delivery FASTA (path
# from config.references.<ref>.source_path, relative to D_DATA) into the
# canonical per-aligner location so all aligners share a uniform ref-name
# basename.
#
# Required wrapper-provided variables:
#   D_DATA, D_REF, D_LOGS, D_BENCHMARK, CONDA_NF1_MOUSE, REF_META
#
#########1#########2#########3#########4#########5#########6#########7#########8


# ── Stage source FASTA into per-aligner canonical path ────────────────────
rule ref_stage_fasta:
    message: "ref: stage source FASTA → {output.fa}"
    conda: CONDA_NF1_MOUSE
    input:
        src = lambda wc: f"{D_DATA}/{REF_META[wc.ref_name]['source_path']}",
    log:
        cmd = f"{D_LOGS}/{{align_method}}.{{ref_name}}_ref_stage_fasta.log",
    output:
        fa = f"{D_REF}/{{align_method}}/{{ref_name}}/{{ref_name}}.fa",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")" "$(dirname "{output.fa}")"
        exec &>> "{log.cmd}"
        echo "[ref_stage_fasta] $(date) aligner={wildcards.align_method} ref={wildcards.ref_name}"
        ln -sf "{input.src}" "{output.fa}"
        ls -l "{output.fa}"
        """


# ── BWA-meth index ─────────────────────────────────────────────────────────
# bwameth.py index runs `bwa index` (v1 format) on the C-to-T converted FASTA.
# bwa-meth's align step uses `bwa-mem2 mem`, which requires the v2 index files
# (.0123, .bwt.2bit.64) — so we follow with `bwa-mem2 index` on the same .c2t
# reference.
rule ref_bwameth_index:
    message: "ref: bwameth.py + bwa-mem2 index {wildcards.ref_name}"
    conda: CONDA_NF1_MOUSE
    input:
        fa = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_ref_bwameth_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_ref_bwameth_index.tsv"
    threads: 8
    output:
        c2t     = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t",
        amb     = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t.amb",
        ann     = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t.ann",
        pac     = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t.pac",
        bm2_idx = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t.0123",
        bm2_bwt = f"{D_REF}/bwameth/{{ref_name}}/{{ref_name}}.fa.bwameth.c2t.bwt.2bit.64",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")"
        exec &>> "{log.cmd}"
        echo "[ref_bwameth_index] $(date) ref={wildcards.ref_name} threads={threads}"
        echo "  step 1/2: bwameth.py index ({wildcards.ref_name})"
        bwameth.py index "{input.fa}"
        echo "  step 2/2: bwa-mem2 index on .bwameth.c2t"
        bwa-mem2 index "{input.fa}.bwameth.c2t"
        """


# ── BISCUIT index ──────────────────────────────────────────────────────────
rule ref_biscuit_index:
    message: "ref: biscuit index {wildcards.ref_name}"
    conda: CONDA_NF1_MOUSE
    input:
        fa = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_ref_biscuit_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_ref_biscuit_index.tsv"
    threads: 8
    output:
        # BISCUIT writes .bis.{amb,ann,pac}, .par.{bwt,sa}, .dau.{bwt,sa}
        bis_amb = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.bis.amb",
        bis_ann = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.bis.ann",
        bis_pac = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.bis.pac",
        par_bwt = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.par.bwt",
        par_sa  = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.par.sa",
        dau_bwt = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.dau.bwt",
        dau_sa  = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.dau.sa",
        fai     = f"{D_REF}/biscuit/{{ref_name}}/{{ref_name}}.fa.fai",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")"
        exec &>> "{log.cmd}"
        echo "[ref_biscuit_index] $(date) ref={wildcards.ref_name} threads={threads}"
        biscuit index "{input.fa}"
        samtools faidx "{input.fa}"
        """


# ── Bismark genome prep ────────────────────────────────────────────────────
# Bismark expects a *directory* containing the FASTA; index goes into a
# Bisulfite_Genome/ subdir there.
rule ref_bismark_index:
    message: "ref: bismark_genome_preparation {wildcards.ref_name}"
    conda: CONDA_NF1_MOUSE
    input:
        fa = f"{D_REF}/bismark/{{ref_name}}/{{ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_ref_bismark_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_ref_bismark_index.tsv"
    params:
        ref_dir = lambda wc: f"{D_REF}/bismark/{wc.ref_name}",
    threads: 8
    output:
        ct_marker = f"{D_REF}/bismark/{{ref_name}}/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2",
        ga_marker = f"{D_REF}/bismark/{{ref_name}}/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2",
    shell:
        """
        mkdir -p "$(dirname "{log.cmd}")"
        exec &>> "{log.cmd}"
        echo "[ref_bismark_index] $(date) ref={wildcards.ref_name} threads={threads}"
        bismark_genome_preparation --bowtie2 --parallel {threads} "{params.ref_dir}"
        """
