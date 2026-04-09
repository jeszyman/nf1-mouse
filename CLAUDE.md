# NF1 Mouse — K08 Aim 2

## Project Identity
- K08 Aim 2: cfDNA monitoring of tumor biology and therapy response in NF1 MPNST PDX models
- 6 PDX lines (JH-2-002, WU-356, JH-2-055, WU-487, WU-225, MN2), 140 mice, 5 treatment arms
- Exploratory analysis repo — not a pipeline/biopipe

## Code Conventions
- Primary org file: `nf1-mouse.org` — notes, analysis code, and literate docs live here
- Exploratory analysis notebooks go in `dev/<analysis-name>/` (R scripts, not org-babel)
- `data/metadata.xlsx` is the single source of truth for metadata (6 relational sheets)
- Schema changes must update both org YAML blocks and the xlsx sheets
- Conda env: `nf1-mouse` (yaml at `config/nf1-mouse-conda-env.yaml`)

## Key Files
- `nf1-mouse.org` — primary org file
- `data/metadata.xlsx` — metadata workbook (mice, pdx_lines, samples, cfdna_isolations, sequencing_libraries, tumor_measurements)
- `config/common.yaml` — project config
- `dev/` — exploratory analysis notebooks (cna-sub50M, icrg-analysis, pdx-read-handling)
- `resources/nf1-mouse.bib` — project bibliography

## Data
- Heavy data at `/mnt/data/projects/nf1-mouse/` — repo stays light
- GCS is source of truth for raw inputs
- Do not commit large data files to the repo
