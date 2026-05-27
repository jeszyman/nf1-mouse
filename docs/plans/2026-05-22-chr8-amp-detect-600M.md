# Chr8q amp detection — expand to equal-depth ~600M-pair re-run

## Context

Local Claude sessions on the `big-boy` GCE VM keep getting stuck mid-execution
of the chr8-amp-detect plan
(`~/repos/nf1-mouse/docs/plans/2026-05-18-chr8-amp-detect.md`). That plan
ran the 4-sample feasibility test at 50M read pairs and produced dedup BAMs
for all 4 libraries but never reached Phase 4 (ichorCNA) or Phase 5 (chr8
plots). At 50M, signal was equivocal; the Hirbe/Wang JCI revision needs
the strongest possible answer to Reviewer A's question, so the user wants
to push depth as high as the data allows.

Goal: re-run the same 4 samples (lib0626, lib0627, lib0644, lib0645) at
the maximum **equal** read-pair depth available across all four libraries,
rounded down to a friendly round number (e.g. if smallest = 612M pairs,
use 600M). Keep the existing 50M outputs intact for depth comparison.
Phases 4–6 of the original plan then run on the new BAMs.

Why equal depth: ichorCNA's per-bin coverage and the visual chr8 log2-ratio
readout are both depth-sensitive. WT-vs-AMP visual contrast is only
interpretable if every sample sees the same number of reads attempting
alignment.

## Decisions locked in this session

- Target depth: **floor of the smallest library's available read-pair count, rounded down to nearest 50M.** Computed dynamically by `01_subsample.sh` after counting all four raw FASTQs.
- Parallelism: **serial, 96 `bwameth.py` threads per sample.** Prior 4-way parallel at 50M crashed the VM; serial 96-thread is the safe high-throughput choice on a 224-vCPU box.
- Output layout: **new subdir** `chr8-amp-detect/600M/` (use the actual chosen N in the directory name; assume `600M` as placeholder for the plan). 50M outputs untouched.
- FASTQ access: **`gsutil` only** for lib0626/0627. No reading from the `/mnt/gcs/` mount.
- lib0644/0645 raw FASTQs already on disk under `chr8-amp-detect/fastqs/lib0644/` and `.../lib0645/` (NCI Flowcell_235NHYLT3) — reuse, don't re-download.

## Source paths (verified from repo)

- `lib0626` → `gs://chaudhuri-lab-bucket1/jeszyman/data/nf1/Project_JohnShern_CS039619_112EMseqlib_061225/Flowcell_22WFNNLT4/Sample_99_nf1_lib99/` (mapping `lib0626=nf1_lib99` per `results/qc/fastqc-summary.tsv` and `nf1-mouse.org` Run-1 mapping `lib0622-lib0639 = nf1_lib95-lib112`)
- `lib0627` → `gs://chaudhuri-lab-bucket1/jeszyman/data/nf1/Project_JohnShern_CS039619_112EMseqlib_061225/Flowcell_22WFNNLT4/Sample_100_nf1_lib100/`
- `lib0644` → already on VM: `/mnt/data/projects/nf1-mouse/chr8-amp-detect/fastqs/lib0644/5_MN2_4_5_2_12_26_S5_R{1,2}_001.fastq.gz` (also `gs://chaudhuri-lab-bucket1/jeszyman/data/Flowcell_235NHYLT3/Sample_5_MN2_4_5_2_12_26/`)
- `lib0645` → already on VM: `/mnt/data/projects/nf1-mouse/chr8-amp-detect/fastqs/lib0645/6_MN2_3_4_2_12_26_S6_R{1,2}_001.fastq.gz`

## Compute / disk budget

- `big-boy`: 224 vCPU, 220 GB RAM, 2.5 TB free on `/mnt/data` (verified this session).
- Raw FASTQs on disk for the NCI libs are ~46–49 GB per mate (≈700M+ pairs each). lib0626/lib0627 likely similar based on the prior 18-sample EM-seq batch. Expected target after equal-depth flooring: 600M pairs each.
- 600M-pair `bwa-meth` BAM per sample ≈ 100–110 GB sorted + 100 GB dedup. 4 × 2 ≈ 800 GB intermediates. Comfortably fits in 2.5 TB free, with `sorted.bam` deleted after dedup to keep peak below 1.5 TB.
- Wall time at 96 threads, serial, 600M pairs: ~3–6 h per sample for `bwameth.py` + sort + dedup → ~12–24 h total alignment phase.
- Phase 4 (readCounter + ichorCNA) is minutes per sample.
- Phase 5 (chr8 plots + summary) is seconds.

## File layout — new artifacts

```
~/repos/nf1-mouse/dev/chr8-amp-detect/
  01b_subsample_equal.sh        # NEW — count + equal-depth subsample
  02b_bwameth_align_serial.sh   # NEW — serial driver for 4 libs at 600M
  03b_ichorcna_pon.sh           # already exists; reuse, point at new BAMs
  04b_chr8_plots_600M.R         # NEW — plots for the 600M depth, plus
                                #       depth-comparison overlay vs 50M
  README.org                    # update with the 600M run section

/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/
  subsampled/                   # equal-depth fq.gz pairs
  bams/                         # *.hg38.dedup.bam(.bai), flagstat, idxstats
  ichorcna/                     # per-sample ichorCNA dirs
  plots/                        # chr8 per-sample + genome-wide PNGs
  logs/                         # per-phase logs
```

## Plan

### Phase A — Raw FASTQ staging (lib0626/0627 only)

1. `gsutil -m rsync -r` the two `Sample_99_nf1_lib99/` and `Sample_100_nf1_lib100/` directories to `/mnt/data/projects/nf1-mouse/chr8-amp-detect/fastqs/lib0626/` and `.../lib0627/`. Use the same `sudo -n -u jupyter` execution pattern as Phase 0 so the files inherit the `jupyter` group.
2. Verify R1+R2 present and non-zero per library.

### Phase B — Equal-depth selection and seqtk subsample

`dev/chr8-amp-detect/01b_subsample_equal.sh`:

1. For each of {lib0626, lib0627, lib0644, lib0645}, count R1 reads via `zcat <R1>.fastq.gz | awk 'NR%4==1' | wc -l`. Run the four counts in parallel (each is single-core; 4× ~30 min on these file sizes).
2. `MIN = min(counts)`; `TARGET = floor(MIN / 50_000_000) * 50_000_000`. Echo all four counts and the chosen TARGET to the log.
3. For each lib, `seqtk sample -s42 <R1> $TARGET | gzip > <R1>.${TARGET_TAG}.fq.gz` and same R2 with same seed. `TARGET_TAG` like `600M`.
4. Verify: each output R1/R2 has exactly TARGET reads and matched IDs (sample first 4 read names from R1 + R2, confirm they share `/1` vs `/2` or matching base IDs).

Outputs: `/mnt/data/projects/nf1-mouse/chr8-amp-detect/600M/subsampled/<lib>_R{1,2}.${TARGET_TAG}.fq.gz`.

### Phase C — Serial bwa-meth alignment

`dev/chr8-amp-detect/02b_bwameth_align_serial.sh` — driver that calls a per-library inner function (or `02_bwameth_align.sh` adapted with new in/out paths and `THREADS=96`). For each lib in order:

```
bwameth.py --threads 96 --reference $REF $R1 $R2 \
  | samtools fixmate -m -u -@ 16 - - \
  | samtools sort -@ 16 -m 4G -o $SORTED -
samtools index -@ 16 $SORTED
samtools markdup -r -@ 16 $SORTED $DEDUP
samtools index -@ 16 $DEDUP
samtools flagstat -@ 16 $DEDUP > $LIB.dedup.flagstat.txt
samtools idxstats $DEDUP > $LIB.dedup.idxstats.txt
rm -f $SORTED $SORTED.bai      # free ~100 GB after dedup verified
```

- Reference, conda env, and `samtools` invocation match `02_bwameth_align.sh`.
- `bwameth.py` at 96 threads is the long pole; sort/markdup at 16 keeps RAM safe.
- Per-lib log files to `chr8-amp-detect/600M/logs/<lib>.align.log` with `set -o pipefail` and `date`-stamped phase markers, so a stuck or crashed session can be resumed at the next lib.
- Resume logic: if `$DEDUP` and `$DEDUP.bai` both exist non-empty, skip the lib (mirrors the existing script).

### Phase D — ichorCNA (no-PoN, default-low-input, autosomes only)

Reuse `dev/chr8-amp-detect/03_ichorcna.sh` (or its `03b` variant), pointed at the new `600M/bams/` directory:

1. `/opt/conda/envs/hmmcopy-bin/bin/readCounter -w 1000000 -q 20 -c chr1,chr2,...,chr22 <dedup.bam> > <lib>.1Mb.wig`.
2. Run patched ichorCNA as `jupyter` with `R_LIBS=/mnt/data/projects/nf1-mouse/R_libs`, default-low-input preset, no PoN, autosomes only — same invocation conventions documented in Phase 4 of the upstream plan.
3. Outputs into `chr8-amp-detect/600M/ichorcna/<lib>/`.

### Phase E — chr8 coverage plots and depth-comparison summary

`dev/chr8-amp-detect/04b_chr8_plots_600M.R`:

1. Per sample, parse ichorCNA per-bin `*.cna.seg` (or `*.bin.tsv`) for chr8, plot log2 ratio vs chr8 position with centromere (~chr8:45 Mb) and p/q boundary annotated. 4-panel layout: AMP samples top row, WT positive controls bottom row.
2. Genome-wide log2 ratio profile per sample — same style as `dev/cna-sub50M/plots/cna_pub_*.png`.
3. Depth-comparison overlay: for each sample, overlay the 50M and 600M chr8 log2 traces. This is the deliverable that directly shows whether deeper sequencing resolves the equivocal Chr8q signal seen at 50M.
4. TSV summary: lib, line, chr8q_status, depth (50M or 600M), mapped_reads, dup_pct, ichor_tf, median_log2_chr8p, median_log2_chr8q.

Adapt from `dev/cna-sub50M/03_alignment_metrics.R` and `06_cna_profiles.R` rather than writing from scratch (per the existing plan's "Existing utilities to reuse" section).

### Phase F — Org write-up (`nf1-mouse.org`)

Add a new subsection under the existing "Chr8q amp detection feasibility" heading (created in Phase 6 of the upstream plan), titled "Equal-depth ~600M-pair re-run":

- Motivation (50M signal equivocal; reviewer needs strongest answer)
- Equal-depth target chosen, with the four raw read counts and the floor-rounded N
- Methods delta vs 50M (same pipeline, just deeper)
- Result panels: per-sample chr8 plots, genome-wide plots, 50M/600M overlay, TSV summary
- Result summary in sci-write voice answering the Reviewer A question definitively at this depth
- Limitations carried forward from the 50M plan (n=2/arm, NCI vs WashU prep, no PoN, bwa-meth vs BISCUIT)

Use `org-edit` skill before structural edits. Tangle so the new bash and R blocks write to `dev/chr8-amp-detect/`.

## Critical files to modify

- `dev/chr8-amp-detect/01b_subsample_equal.sh` (new) — equal-depth seqtk
- `dev/chr8-amp-detect/02b_bwameth_align_serial.sh` (new) — serial 96-thread driver
- `dev/chr8-amp-detect/04b_chr8_plots_600M.R` (new) — chr8 plots + 50M/600M overlay
- `dev/chr8-amp-detect/README.org` — note the 600M sibling run
- `nf1-mouse.org` — new "Equal-depth ~600M-pair re-run" subsection

## Existing utilities to reuse (from `~/repos/nf1-mouse`)

- `dev/chr8-amp-detect/02_bwameth_align.sh` — alignment recipe (`bwameth.py | fixmate | sort | markdup` + flagstat/idxstats). Wrap in serial driver with new paths + 96 threads.
- `dev/chr8-amp-detect/02_resume_serial.sh` — serial-driver structure
- `dev/chr8-amp-detect/03_ichorcna.sh` and `03b_ichorcna_hdulp_pon.sh` — ichorCNA invocation already wired with the `R_LIBS`/`jupyter`/`readCounter`-absolute-path conventions
- `dev/cna-sub50M/03_alignment_metrics.R`, `06_cna_profiles.R`, `17_cna_sidebyside.R` — log2-ratio plotting templates; the 50M/600M overlay is closest to `17_cna_sidebyside.R`
- ichorCNA patched fork + R libs at `/mnt/data/projects/nf1-mouse/R_libs/` (already installed per Phase 4 of upstream plan, dated 2026-05-22)
- `readCounter` at `/opt/conda/envs/hmmcopy-bin/bin/readCounter`
- Reference + bwa-meth index at `/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/`

## Verification

1. Phase A: `gsutil ls` confirms lib0626/0627 in `chr8-amp-detect/fastqs/`; both R1 and R2 non-zero.
2. Phase B: for each of 4 libs, R1 read count == R2 read count == chosen TARGET; first/last read IDs paired.
3. Phase C: per-sample `flagstat` shows primary-mapped reads in expected band (≥16% of TARGET reads, scaled from the 50M observation of 8.5–11.5M mapped per 50M pairs); dedup rate plausible (<25%).
4. Phase D: ichorCNA exits 0 for all 4; `<lib>.params.txt` reports a converged TF.
5. Phase E: 4 chr8 panels, 4 genome panels, 4 depth-overlay panels, 1 TSV — all produced. Visual inspection records whether MN2 samples now show a clear p→q step.
6. Phase F: org section exports to HTML via the project's running-results convention; all figure links resolve.

## Out of scope (carry over from upstream plan)

- Strategy 3A disambiguation
- Serial-timepoint analysis
- All 14 NCI samples (only top 2 by bioanalyzer)
- A 2nd Chr8q-gain PDX line
- PoN re-build
- Methylation analysis on these BAMs

## Follow-on notes (post-2026-05-26)

- **No-PoN + HD-ULP-PoN ichorCNA on the 450M bwa-meth hg38-only BAMs produced shared-bias artifacts**: all 4 samples (2 WT + 2 AMP) showed similar genome-wide CNA patterns with implausibly high TFs (43-55%). MN2 samples called as chr8 HETD (deletion) instead of expected chr8q gain. Diagnosed as mouse-host cfDNA cross-mapping into the human reference under Strategy 1 (direct human-only alignment). HD-ULP PoN partially removed the artifact but did not recover the chr8q gain signal.
- **In-progress remedy**: Strategy 3A disambiguation via dual bwa-meth alignment (hg38 + mm10) followed by NM-tag comparison. mm10 bwa-meth index build kicked off 2026-05-26. Alignment driver staged at =dev/chr8-amp-detect/02c_bwameth_align_serial_mm10.sh=. bwa-meth-aware disambiguator at =dev/pdx-read-handling/disambiguate_bwameth.py= (functionally identical to BISCUIT version; docstring documents caveats).
- **Deferred: hg38 re-alignment of human-classified reads after disambiguation.** Common downstream step after PDX read parsing. Skipped for the immediate chr8q-amp ichorCNA question because (a) ichorCNA bins by read count and tolerates broken-pair flags, and (b) it would add ~40h compute. Add later as =02d_bwameth_realign_human.sh= if pursuing methylation calling or fragmentomics on these libraries (EM-seq data — eventually wanted), or if reviewers ask for clean paired-end disambiguated BAMs.

## TODO: Cross-classification analysis (deferred)

For each of the 4 libs, compute the per-read cross-classification table comparing the straight bwa-meth hg38-only verdict (any read mapped to hg38) against the post-disambiguation verdict (human / mouse / ambiguous). Smoke-test on 50k chr1 reads from lib0626 [2026-05-26] showed ~60% of reads that bwa-meth assigned to hg38 are actually mouse, ~9% ambiguous, ~38% real human. The per-sample full-BAM numbers will quantify mouse-contamination rate per library, and the per-bin / per-chromosome breakdown will reveal which genomic regions concentrate the cross-mapping (likely pericentromeric, ribosomal, and other high-homology zones). Useful for:

- Methods-section reporting of contamination magnitude
- Identifying regions that would benefit from extra masking even after disambiguation
- Comparing host-contamination rates across sample lines (WU-487 vs MN2 — host:graft ratios likely differ)
- Sanity-checking that the disambiguators bwa-meth NM scoring behaves consistently across libraries

## TODO: Cross-classification analysis (deferred)

For each of the 4 libs, compute the per-read cross-classification table comparing the straight bwa-meth hg38-only verdict (any read mapped to hg38) against the post-disambiguation verdict (human / mouse / ambiguous). Smoke-test on 50k chr1 reads from lib0626 [2026-05-26] showed ~60% of reads that bwa-meth assigned to hg38 are actually mouse, ~9% ambiguous, ~38% real human. The per-sample full-BAM numbers will quantify mouse-contamination rate per library, and the per-bin / per-chromosome breakdown will reveal which genomic regions concentrate the cross-mapping (likely pericentromeric, ribosomal, and other high-homology zones). Useful for:

- Methods-section reporting of contamination magnitude
- Identifying regions that would benefit from extra masking even after disambiguation
- Comparing host-contamination rates across sample lines (WU-487 vs MN2 — host:graft ratios likely differ)
- Sanity-checking that the disambiguators bwa-meth NM scoring behaves consistently across libraries

Script template: cross_classify.sh used for the smoke test; adapt to operate on full BAMs and emit per-chromosome contamination rates.
