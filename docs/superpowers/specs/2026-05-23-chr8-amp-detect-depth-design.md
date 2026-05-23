# Chr8q amp detection — depth-increase phase (design)

## Context

Follow-on to the 50M-pair feasibility test in
`docs/plans/2026-05-18-chr8-amp-detect.md` (Phases 0–5 complete). That
run aligned 50M-pair subsamples of 4 libraries (2 MN2 chr8q-AMP, 2
WU-487 chr8q-WT) with bwa-meth direct-to-hg38 and ran ichorCNA in two
configs (no-PoN, HD-ULP-PoN). Outcome: ichorCNA did NOT call chr8q as
GAIN in either MN2 sample under either config. Per-bin
`delta_q_minus_p` was 0.21/0.23 for MN2 (AMP) vs 0.16/0.10 (no-PoN) or
0.19/0.05 (HD-ULP-PoN) for WU-487 (WT) — MN2 signal sits within the WT
noise band. Visual chr8q step-up was not convincing. See
`dev/chr8-amp-detect/plots/chr8_zoomed_nopon.png` and
`chr8_zoomed_hdulp-pon.png` and the two `summary_*.tsv` files.

Diagnosis: per-bin sampling variance, not disambiguation. The
ultimately correct fix is a mouse-plasma PoN built from naive PDX
mouse cfDNA EM-seq libraries. The near-term fix is to throw more
reads at the same pipeline. All 4 libraries have substantial unused
depth on disk:

| library | source | R1 size | est. pair count | currently used |
|---------|--------|---------|-----------------|----------------|
| lib0626 | WashU  | 25 GB   | ~595M           | 50M (~8%)     |
| lib0627 | WashU  | 27 GB   | ~625M           | 50M (~8%)     |
| lib0644 | NCI    | 46 GB   | ~1.3–1.5B       | 50M (~3%)     |
| lib0645 | NCI    | 49 GB   | ~1.3–1.5B       | 50M (~3%)     |

No resequencing required.

## Goal & readout

Test whether the chr8q amp in MN2 PDX cfDNA becomes cleanly
detectable when the same bwa-meth direct pipeline is re-run at ~12×
deeper sequencing depth.

- Primary readout: visual step-up in per-bin log2 ratio at the chr8
  p→q transition in MN2 (lib0644, lib0645), absent in WU-487 WT
  (lib0626, lib0627).
- Secondary: ichorCNA actually calls chr8q as a `GAIN` segment in
  MN2 in at least one PoN config (it didn't at 50M).
- Quantitative: `delta_q_minus_p` MN2-vs-WT separation should grow
  roughly 3× the 50M ratio (noise scales as 1/√N, expected ~3.46×
  noise drop from 12× depth).

Answers Hirbe/Wang JCI reviewer A #4 (chr8 gain serially monitorable
in mouse blood).

## Design decisions (locked)

1. **Intervention**: depth increase only. No disambiguation (S3A),
   no aligner change. Single-variable A/B vs the 50M no-PoN and
   HD-ULP-PoN results. Disambiguation is the deferred fallback if
   depth alone fails.
2. **Target depth**: 600M pairs each, matched. lib0626 likely uses
   essentially all ~595M reads available (no subsample cushion);
   lib0644/lib0645 still have ~2× headroom. Acceptable.
3. **Aligner**: bwa-meth → hg38 no-alt (unchanged from 50M run).
   Reference + index already at
   `/mnt/data/projects/nf1-mouse/ref/biscuit/ncbi_hg38_noalt/`.
4. **ichorCNA configs**: both no-PoN and HD-ULP-PoN, both
   default-low-input preset, autosomes only. Matches 50M run
   structure for clean comparison.
5. **Subsampling**: `seqtk sample -s42 <R1> 600000000` paired with
   same seed for R2. Same seed as 50M run; the 50M and 600M subsamples
   are independent draws from the same source for samples where the
   source > 600M, but lib0626 (~595M) will be effectively the entire
   library.
6. **Compute strategy**: serial alignment (4 samples × ~2–3 hr each
   = ~10–12 hr wall). The 4-way parallel attempt crashed the VM at
   50M; do not retry parallelism at 600M.
7. **Storage**: ~200 GB subsampled fastqs (retained) + ~250–300 GB
   BAMs. Fits in 4 TB VM disk. Outputs go to a parallel tree
   `_600M` so 50M results remain intact for side-by-side figures.
8. **Plan file**: new `docs/plans/2026-05-23-chr8-amp-detect-depth.md`,
   not an amendment to the 50M plan. The 50M plan is closed except
   for Phase 6 (org write-up), which will fold both 50M and 600M
   results into one `nf1-mouse.org` section.
9. **Deliverable figure**: 2-row × 4-column comparison panel per
   PoN config (rows = depth: 50M, 600M; cols = lib0644, lib0645,
   lib0626, lib0627). Two PNGs total. Plus a unified `summary.tsv`
   with both depths.

## Pipeline

All scripts live in `dev/chr8-amp-detect/` and reuse the existing
files where possible:

| step | script | reuse status |
|------|--------|--------------|
| subsample | NEW `01_subsample_600M.sh` | new wrapper around `seqtk sample` |
| align | `02_bwameth_align.sh` | reuse unchanged; called per lib with new output dir |
| serial driver | `02_resume_serial.sh` | reuse pattern (new copy `02_serial_600M.sh`) |
| ichor no-PoN | `03_ichorcna.sh` | reuse unchanged; per-lib BAM input |
| ichor HD-ULP-PoN | `03b_ichorcna_hdulp_pon.sh` | reuse unchanged |
| plots | EXTEND `04_chr8_plots.R` | extend to produce a 50M-vs-600M comparison panel |

Operational conventions (locked from prior phases):

- Run as `sudo -n -u jupyter` for anything touching
  `/mnt/data/projects/nf1-mouse/` (jupyter-owned tree, `drwxrws---`).
- ichorCNA driver invoked via `/usr/bin/Rscript` with
  `R_LIBS=/mnt/data/projects/nf1-mouse/R_libs` exported.
- `readCounter` at `/opt/conda/envs/hmmcopy-bin/bin/readCounter`
  (absolute path; not on default PATH).
- Long alignment runs MUST use `nohup` per
  [[feedback_nohup]] — tool timeouts kill background tasks otherwise.

## File layout

```
dev/chr8-amp-detect/
  01_subsample.sh                 # existing 50M (unchanged)
  01_subsample_600M.sh            # NEW: 600M subsample
  02_bwameth_align.sh             # existing (unchanged)
  02_resume_serial.sh             # existing (unchanged)
  02_serial_600M.sh               # NEW: serial driver pointing at 600M tree
  03_ichorcna.sh                  # existing (unchanged)
  03b_ichorcna_hdulp_pon.sh       # existing (unchanged)
  04_chr8_plots.R                 # EXTEND: add comparison-panel mode
  plots/
    chr8_zoomed_nopon.png         # existing 50M
    chr8_zoomed_hdulp-pon.png     # existing 50M
    chr8_compare_50M_600M_nopon.png   # NEW
    chr8_compare_50M_600M_hdulp-pon.png # NEW
    summary_combined.tsv          # NEW: 50M + 600M rows
```

On the VM:

```
/mnt/data/projects/nf1-mouse/chr8-amp-detect/
  subsampled/                     # existing 50M fastqs
  subsampled_600M/                # NEW
  bams/                           # existing 50M BAMs
  bams_600M/                      # NEW
  ichorcna/                       # existing 50M ichor outputs
  ichorcna_600M/                  # NEW
```

## Phases (preview)

1. **Re-subsample to 600M** — seqtk sample -s42 R1/R2 to
   `subsampled_600M/`. Verify pair counts.
2. **bwa-meth align (serial)** — 4 samples sequentially, ~10–12 hr
   wall, under nohup. QC: per-sample human-mapped ≥ 80M, dup rate
   reasonable.
3. **ichorCNA both PoN configs** — 8 runs total (4 libs × 2 configs).
   readCounter once per BAM, ichorCNA twice per BAM.
4. **Comparison figure + combined summary TSV** — extend the R
   script; write the 2×4 panel PNGs and the combined TSV.
5. **Org write-up** — fold both 50M and 600M into the existing
   pending Phase 6 of the 50M plan; one section in `nf1-mouse.org`
   under "Chr8q amp detection feasibility".

## Verification

1. Per-sample human-mapped reads ≥ 80M (target was ≥10× the 50M
   floor of 8M).
2. ichorCNA converges for all 8 runs; TF estimates reported.
3. In at least one PoN config, both MN2 samples have a chr8q segment
   called as `GAIN` (or higher); both WT samples remain `NEUT` on
   chr8q.
4. `delta_q_minus_p` separation between MN2 and WT grows ≥ 2×
   relative to the 50M run (allowing some slack vs the theoretical
   3× from 1/√N).
5. Comparison PNGs render the chr8q step-up visibly in MN2 at 600M
   and not at 50M.

## Out of scope

- Disambiguation (Strategy 3A) — deferred fallback if depth alone fails.
- Mouse-plasma PoN — the ultimately correct fix; separate future
  effort (requires naive PDX mouse cfDNA EM-seq).
- All 14 NCI samples — only the original top-2 MN2 by bioanalyzer
  here.
- A 2nd chr8q-gain PDX line beyond MN2.
- Methylation analysis on these BAMs (CNA-only).

## Open details (handle in implementation plan)

- Whether to delete the 50M subsampled fastqs after 600M completes
  (saves ~20 GB, but they're cheap to keep).
- Whether the comparison figure should add a delta_q_minus_p
  annotation per panel, or keep that only in the summary TSV.
- Whether to also produce genome-wide quicklook PNGs for the 600M
  runs in addition to chr8-zoomed (the 50M run produced both).
