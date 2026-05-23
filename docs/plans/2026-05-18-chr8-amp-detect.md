# Chr8q amp detection feasibility test — Hirbe/Wang JCI revision

## Context

The Hirbe/Wang JCI manuscript ("Dual inhibition of FAK and MEK is a
therapeutic vulnerability in chromosome 8 gain MPNST", 201331-JCI-RG-1,
Guangfeng Wang first author) is in major revision. Reviewer A comment
#4 asked whether Chr8 gain can be serially monitored in mouse blood —
this is the direct motivation for our PDX cfDNA Chr8q work.

Project status at the time of this plan:

- 18 existing EM-seq libraries are all **Chr8q-WT** (6 WU-487 terminal,
  12 JH-2-055 serial). Gap identified [2026-03-18].
- A 14-sample shipment of MN2 (Chr8q-gain) and other Chr8q-gain
  plasma went to NCI on [2026-03-30]. The NCI flowcell
  `Flowcell_235NHYLT3` (Project_JohnShern_CS041084_14EMseqlib_042826)
  was transferred to GCS at
  `gs://chaudhuri-lab-bucket1/jeszyman/data/Flowcell_235NHYLT3/`
  on [2026-05-12].
- The sub-50M cohort analysis (`dev/cna-sub50M`, 6 WU-487 terminals)
  established that on Chr8q-WT samples, direct human-only alignment
  (Strategy 1) compresses tumor fraction 3-5x vs disambiguation
  (Strategy 3A), but two high-TF samples (lib0626 S3A TF=0.46, lib0627
  S3A TF=0.65) showed *some* CNA signal even on S1. Chr8q is a broad
  ~95 Mb event; the open question is whether it's broad enough to be
  visible in direct alignment of high-TF Chr8q-gain samples.

This plan tests the simplest fast pipeline (bwa-meth human-only direct
alignment, no disambiguation, 50M read pairs) on 2 Chr8q-AMP samples
and 2 Chr8q-WT positive-control samples, with the question:
**"Does chr8q gain show up at all under the fastest possible pipeline
on the best available samples?"**

If yes: simple/fast pipeline is sufficient for the manuscript
deliverable. If no: fall back to disambiguation (Strategy 3A) on the
same samples in a follow-on session.

## Readout

Primary readout is visual: per-sample chr8 log2 ratio coverage plot.
If MN2 (AMP) samples show a visible step-up in log2 ratio at the
p-to-q transition around the centromere (~chr8:45 Mb) and WU-487
samples do not, the chr8q amp is detectable under the simple pipeline
and answers the reviewer.

Secondary readout: ichorCNA no-PoN tumor fraction estimate per sample
(does not gate the answer; reported alongside).

ichorCNA configuration is **no-PoN, default-low-input preset,
autosomes only**. No with-PoN comparison — the prior BISCUIT PoN is
not validated for bwa-meth alignments and a single configuration keeps
the experiment focused.

## Samples (n=4)

### AMP arm (Chr8q-gain)

| library | isolation | sample  | mouse                  | line | treatment              | date       | bioanalyzer cfDNA conc (pg/uL) |
|---------+-----------+---------+------------------------+------+------------------------+------------+--------------------------------|
| lib0644 | iso0023   | smp0927 | mou0141 ("MN2 4-5")    | MN2  | mirdametinib_volasertib| 2026-02-12 | 251.43                         |
| lib0645 | iso0024   | smp0928 | mou0142 ("MN2 3-4")    | MN2  | mirdametinib_volasertib| 2026-02-12 | 183.58                         |

Top 2 MN2 samples of the 14 NCI shipment samples by bioanalyzer
cfDNA concentration (column `bioanalyzer_conc_pg_ul` in
`cfdna_isolations`; the 167bp peak column is sparsely populated and
was not used for ranking). Both from the combo-treatment arm, both
from the [2026-02-12] timepoint. Lib prep at NCI on flowcell
`235NHYLT3`. Library IDs assigned during Phase 0 metadata gap-fill
(lib0644 = flowcell Sample 5 / iso0023; lib0645 = flowcell Sample 6
/ iso0024).

The 14-sample NCI batch comprises:
- 4 JH-2-055 (Chr8q-WT) re-sequenced serial bleeds, lib0640–lib0643
  (Samples 1–4: mice #4-4 and #2-2, day 5 and day 19 timepoints)
- 10 MN2 (Chr8q-gain) cohort-2 samples, lib0644–lib0653
  (Samples 5–14: mice MN2 4-5, 3-4, 3-3, 3-2, 4-2 across two
  timepoints 2/12/26 and 2/19/26)

Phase 0 also confirmed there is no untreated terminal-bleed
Chr8q-gain backup line on this flowcell — all gain samples are MN2
cohort-2 combo or vehicle arms. So n=1 gain line (MN2) for this
feasibility test.

### WT arm (Chr8q-WT, positive controls for high TF)

| library | sample  | mouse   | line   | treatment    | bleed   | known S3A TF |
|---------+---------+---------+--------+--------------+---------+--------------|
| lib0626 | smp0626 | mou0005 | WU-487 | mirdametinib | terminal| 0.46         |
| lib0627 | smp0627 | mou0006 | WU-487 | mirdametinib | terminal| 0.65         |

Highest-TF samples from prior `dev/cna-sub50M` analysis. Lib prep at
WashU. Note that lib0626 here uses the same library ID as smp0626 in
the cna-sub50M cohort — re-aligning the original FASTQ.

Library prep batch difference (NCI vs WashU) and flowcell difference
are known confounds — unavoidable but acknowledged.

## Pipeline (all 4 samples, bwa-meth human-only direct alignment)

### Phase 0 — Data model gap-fill (DONE 2026-05-18)

Executed on jeff-beast. Outputs:

1. Listed Flowcell_235NHYLT3 contents on GCS. 14 sample directories
   named `Sample_<position>_<line>_<animal>_<m>_<d>_<yy>`. Each
   contains `<position>_<line>_<animal>_<m>_<d>_<yy>_S<N>_R{1,2}_001.fastq.gz`.
2. Built position → isolation → sample → mouse → mouse-tag mapping
   from the cfdna_isolations sheet (positions 1-14 match isolations
   iso0019-iso0032 by shipment_position) and from the sample
   directory names. Verified Sample 5 = iso0023 = smp0927 (MN2 4-5)
   and Sample 6 = iso0024 = smp0928 (MN2 3-4). Flowcell has 4
   JH-2-055 + 10 MN2; no other Chr8q-gain line on this flowcell.
3. Appended 14 rows to `data/metadata.xlsx::sequencing_libraries`
   (lib0640 → lib0653 mapped to iso0021, iso0019, iso0022, iso0020,
   iso0023, iso0024, iso0025, iso0026, iso0027, iso0028, iso0029,
   iso0030, iso0031, iso0032 in flowcell-sample order). Existing
   schema preserved (library_id, cfdna_iso_id, prep_kit,
   sequencing_platform, read_length, lib_conc_ng_ul); prep_kit set
   to `nebmeth_emseq` matching existing NCI-prepped libs; remaining
   fields blank (lib conc and read length not in the manifest yet).
4. dPCR check: `data/sources/nf1_murine_dpcr_batch2_3_21_25.pdf`
   exists from a 2026-03-21 dPCR batch. The current xlsx has no
   structured dPCR column for the NCI MN2 isolations; PDF not parsed
   in Phase 0. Deferred — not blocking the alignment phases.
5. `frictionless validate data/metadata_schema.yaml --type package`
   reports all 6 tables VALID after the row additions.
6. Committed `metadata: add 14 sequencing_libraries rows for NCI
   Flowcell_235NHYLT3 (lib0640-lib0653)`.

Caveat noted during Phase 0: the bioanalyzer numbers cited above
(251.43, 183.58) live in column `bioanalyzer_conc_pg_ul` (total
cfDNA), not `bioanalyzer_167bp_peak_pg_ul` (which is only populated
for a single sample). The MN2 ranking is unchanged.

### Phase 1 — VM and reference setup (DONE 2026-05-19)

1. Provision GCE VM (240-300 vCPU, ~250 GB RAM, ~4 TB local disk for
   FASTQs + BAMs + ichorCNA intermediates). VM type TBD; n2-highmem
   or c2-standard-60s used previously.
2. Bootstrap conda env `nf1-mouse` from `config/nf1-mouse-conda-env.yaml`
   (already has `bwameth` and `biscuit` per the env yaml).
3. Pull hg38 no-alt reference (`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`)
   from the bucket — same reference used for BISCUIT to keep contig
   set consistent.
4. Build bwa-meth index: `bwameth.py index <ref.fa>` — one-time,
   runs in background while Phase 2 proceeds.

Done on `big-boy` VM (the jupyter-owned GCP image). Reference at
`/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`
with bwa-meth index alongside.

### Phase 2 — FASTQ subsample (DONE 2026-05-19)

1. `gsutil rsync` 4 sample FASTQ pairs to VM local disk:
   - 2 MN2 from `Flowcell_235NHYLT3/` (lib IDs from Phase 0)
   - 2 WU-487 (lib0626, lib0627) from existing bucket location
2. `seqtk sample -s42 <R1>.fastq.gz 50000000 | gzip > <R1>.50M.fq.gz`
   and same for R2 with same seed — preserves pair sync.
3. Verify each subsampled fastq has 50,000,000 reads
   (`zcat | awk 'NR%4==1' | wc -l`) and R1/R2 have matched IDs.

Outputs at `/mnt/data/projects/nf1-mouse/chr8-amp-detect/subsampled/`
(`lib0626`, `lib0627`, `lib0644`, `lib0645` × R1/R2 × 50M).

### Phase 3 — bwa-meth alignment (DONE 2026-05-20, serial after VM crash)

For each of 4 samples in parallel:

```
bwameth.py --threads 16 --reference hg38 R1.50M.fq.gz R2.50M.fq.gz \
  | samtools sort -@ 8 -o sample.hg38.sorted.bam -
samtools index sample.hg38.sorted.bam
samtools markdup -r sample.hg38.sorted.bam sample.hg38.dedup.bam
samtools index sample.hg38.dedup.bam
samtools flagstat sample.hg38.dedup.bam > sample.flagstat.txt
samtools idxstats sample.hg38.dedup.bam > sample.idxstats.txt
```

QC outputs: alignment rate, MAPQ distribution, duplication rate,
per-chromosome read counts. Flag any sample with <10% reads mapping
to hg38 or >40% duplication.

Implemented in `dev/chr8-amp-detect/02_bwameth_align.sh` (per-sample)
and `02_resume_serial.sh` (serial driver, after a 4-way-parallel
attempt crashed the VM on 2026-05-18). All 4 dedup BAMs produced at
`/mnt/data/projects/nf1-mouse/chr8-amp-detect/bams/<lib>.hg38.dedup.bam`;
primary-mapped reads 8.5M–11.5M per sample (≥16% of 50M pairs as
specified). Files are owned by the `jupyter` group — use
`sudo -n -u jupyter` for any write into that subtree.

### Phase 4 — ichorCNA (no-PoN)

1. HMMcopy `readCounter` → wig (1 Mb bins, autosomes only) for each
   of the 4 dedup BAMs.
2. ichorCNA with `default-low-input` preset, autosomes only, **no
   PoN**. 4 ichorCNA runs total.

**Tooling install — DONE 2026-05-22:**

- Patched ichorCNA fork installed: `/mnt/data/projects/nf1-mouse/R_libs/ichorCNA`
  (v0.3.2 from `https://github.com/jeszyman/ichorCNA-patched`). All
  required R deps colocated in the same `R_libs/` tree (HMMcopy,
  GenomicRanges, GenomeInfoDb, plyr, optparse, Rcpp, and a **local
  `data.table`** that bypasses the broken system-R 4.6.0 data.table
  with its `SETLENGTH` ABI mismatch). The patched fork was required
  because upstream ichorCNA fails to build under R 4.6.0 due to that
  same data.table ABI issue during byte-compile/lazy-load.
- `readCounter` binary installed: `/opt/conda/envs/hmmcopy-bin/bin/readCounter`
  (bioconda `hmmcopy` package; also provides `mapCounter`, `gcCounter`).
- Verified load: `sudo -n -u jupyter R_LIBS=/mnt/data/projects/nf1-mouse/R_libs /usr/bin/Rscript -e 'suppressMessages(library(ichorCNA))'`
  reports `ichorCNA loaded OK, version: 0.3.2`.

**Invocation conventions (must use):**

- Run any R that touches ichorCNA as the `jupyter` user with
  `R_LIBS=/mnt/data/projects/nf1-mouse/R_libs` exported. The `R_libs/`
  directory is `drwxrws---` jupyter-only; system R is `/usr/bin/Rscript`
  (4.6.0). `R_LIBS` (not `R_LIBS_USER`) is what the prior session
  found honored.
- `readCounter` must be invoked via its absolute path
  (`/opt/conda/envs/hmmcopy-bin/bin/readCounter`) since it lives in
  a dedicated conda env and is not on the default `PATH`.
- The runIchorCNA driver script is staged at `/tmp/runIchorCNA.R`
  (the ichorCNA package's bundled scripts/runIchorCNA.R, copied for
  convenience).

### Phase 5 — Chr8 coverage plots and TF table

R script (`04_chr8_plots.R`) reads ichorCNA per-bin log2 outputs:

1. Per sample, extract per-bin log2 ratios on chr8.
2. Plot per-sample chr8 log2 ratio along chr8 position, with the
   centromere (~chr8:45 Mb) and the p/q boundary annotated. One panel
   per sample, 4 panels total, AMP on top and WT on bottom for visual
   comparison.
3. Also plot the genome-wide log2 ratio profile per sample (4 panels)
   for context — same style as `dev/cna-sub50M/plots/cna_pub_*.png`.
4. Write a small TSV summary: sample, library_id, line, chr8q_status,
   mapped_reads, dup_pct, ichor_tf, median_log2_chr8p,
   median_log2_chr8q. Median chr8p and chr8q values are descriptive,
   not pass/fail.

Pass/fail is read visually from the chr8 coverage plots: AMP samples
should show a visible step-up in log2 ratio at the p→q transition;
WT samples should not.

### Phase 6 — Org write-up

New top-level section in `nf1-mouse.org` (analogous to
"Sub-50M cohort CNA analysis"), e.g. "Chr8q amp detection feasibility":

- Background (Hirbe revision, reviewer comment, prior sub-50M result
  on Chr8q-WT)
- Methods (4 samples, bwa-meth direct, 50M, ichorCNA no-PoN)
- Results (chr8 per-sample coverage plots, ichorCNA TF table,
  genome-wide log2 profiles)
- Result summary in sci-write voice
- Limitations (n=2/arm, NCI vs WashU prep, bwa-meth vs BISCUIT, no
  pre-sequencing TF orthogonal validation, no-PoN baseline noise)

All bash and R blocks tangle to `dev/chr8-amp-detect/`.

## File layout

```
dev/chr8-amp-detect/
  00_setup_refs.sh            # bwameth index on hg38 no-alt
  01_subsample.sh             # seqtk sample 50M pairs
  02_bwameth_align.sh         # per-sample bwa-meth + sort + dedup
  03_ichorcna.sh              # per-sample readCounter + ichorCNA
  04_chr8_plots.R             # chr8 coverage plots, genome plots, TF table
  plots/
  README.org                  # short scope note pointing to nf1-mouse.org
```

Match the dev/cna-sub50M and dev/pdx-read-handling patterns (numbered
bash + R, plots/, exploratory not literate).

## Critical files to modify

- `data/metadata.xlsx` — 14 new `sequencing_libraries` rows (Phase 0)
- `data/metadata_schema.yaml` — only if a new column is required for
  the NCI flowcell rows (e.g., `sequencer_facility`); otherwise
  unchanged
- `nf1-mouse.org` — new "Chr8q amp detection feasibility" section
  (use org-edit skill before structural edits) with tangled blocks
  feeding `dev/chr8-amp-detect/`
- `config/nf1-mouse-conda-env.yaml` — verify `bwameth`, `seqtk`,
  `ichorCNA`, `hmmcopy_utils` are all present; tangle if changed

## Existing utilities to reuse

- `dev/pdx-read-handling/01_strategy1_direct.sh` — template for the
  per-sample alignment script structure (env activation, FASTQ
  discovery, samtools flagstat/idxstats pattern); swap BISCUIT→bwameth
- `dev/cna-sub50M/03_alignment_metrics.R` / `06_cna_profiles.R` —
  templates for log2-ratio aggregation and chr8 plotting; adapt
  rather than re-write
- ichorCNA invocation parameters from `dev/cna-sub50M/` — already
  validated for autosomes + default-low-input
- `data/metadata_schema.yaml` `sequencing_libraries` schema — match
  exactly for the new 14 NCI rows
- Per-sample chr8 coverage plot is a thin adaptation of
  `dev/cna-sub50M/plots/cna_pub_lib0626.png` — subset to chr8, mark
  centromere/q-boundary; ~30 lines of new R

## Compute budget (rough)

- bwameth.py index hg38: 2-4 hours (one-time, can run while subsampling)
- bwameth.py align, 50M pairs, 16 threads: ~30-60 min per sample
  (bwa-mem-based, ~5x BISCUIT throughput). 4 samples in parallel on
  240-vCPU VM: ~1 hour wall.
- HMMcopy readCounter + ichorCNA: minutes per sample
- gsutil rsync: ~20 GB of subsampled FASTQ (faster than pulling full
  ~71 GB libraries); typically <30 min
- Total wall time once VM is up: half a day

## Verification

1. Phase 0 acceptance: `frictionless validate
   data/metadata_schema.yaml --type package` reports all 6 tables
   valid; both MN2 target isolations have lib rows.
2. Phase 3 QC: per-sample human-mapped read count ≥ 8M (≥16% of 50M
   pairs, consistent with prior 18-22% PDX human fraction); dup rate
   reasonable (<25% for these libraries).
3. Phase 4: ichorCNA converges (no fit failures), TF estimates
   reported.
4. Phase 5: chr8 coverage plots produced for all 4 samples; visual
   inspection records whether AMP samples show a p→q step-up and WT
   samples are flat. Result summary states the answer to the
   reviewer's question in one sentence.
5. Org section exports cleanly to HTML via the running-results
   convention (verify ToC, figures resolve).

## Out of scope (deferred to follow-on if Phase 5 says "not detectable")

- Strategy 3A disambiguation of the same 4 samples
- Serial-timepoint analysis (combo timepoint 1 vs timepoint 2)
- All 14 NCI samples (only top 2 by bioanalyzer here)
- A 2nd Chr8q-gain PDX line beyond MN2 (unless found in Phase 0)
- PoN re-build from JH-2-055 normals
- Methylation analysis on these BAMs (this is CNA-only)
