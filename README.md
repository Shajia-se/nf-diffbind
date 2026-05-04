# nf-diffbind

Nextflow DSL2 module for DiffBind differential binding analysis.

## Input Modes (Priority Order)

1. `--samplesheet` (existing DiffBind sample sheet)
2. `--samples_master` auto-generation

## Mode 1: Explicit `samplesheet`

Required columns:
- `SampleID`
- `Condition`
- `Replicate`
- `bamReads`
- `Peaks`
- `PeakCaller`

Optional:
- `bamControl`
- `Batch`

## Mode 2: Auto from `samples_master`

Required columns:
- `sample_id`
- `condition`

Optional columns used:
- `replicate`
- `library_type`
- `is_control`
- `use_for_diffbind`
- `enabled`

Auto behavior:
- keeps enabled, non-control, chip rows with `use_for_diffbind=true`
- resolves BAM from `${chipfilter_output}/${sample_id}*.nomulti.bam`
- resolves peak from `${macs3_output}/strict_q0.01/${sample_id}_peaks.${diffbind_peak_ext}` by default
- override profile with `--diffbind_macs3_profile`
- writes an internal generated CSV and runs DiffBind with it

## Run

Explicit sheet:
```bash
nextflow run main.nf -profile hpc --samplesheet samplesheet.csv
```

Auto from `samples_master`:
```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output \
  --macs3_output /path/to/nf-macs3/macs3_output
```

Default MACS3 profile used by DiffBind: `strict_q0.01`

Resume:
```bash
nextflow run main.nf -profile hpc -resume
```

## Key Outputs

Output directory: `${project_folder}/${diffbind_output}`

Core DiffBind:
- `01_general_QC.pdf`
- `sample_info.after_count.tsv`
- `libsizes.tsv`, `libsizes.xlsx`
- `contrasts.tsv`
- `significant.<contrast>.tsv/.xlsx`
- `all_peaks.<contrast>.tsv/.xlsx/.bed`
- `db_after_count_and_normalize.RData`
- `diffbind_session.RData`

Added overlap/summary outputs:
- `peak_universe_upset_input.tsv`
- `peak_universe_condition_sizes.tsv`
- `peak_universe_pairwise_overlap.tsv`
- `peak_universe_venn.pdf` (optional, if `VennDiagram` is installed)
- `condition_unique_up.<contrast>.bed`
- `condition_unique_down.<contrast>.bed`
- `diffbind_summary.tsv`, `diffbind_summary.xlsx`

## Parameters

- `diff_fdr` (default `0.05`): FDR threshold for unique peak export
- `diff_lfc` (default `1`): absolute Fold threshold for unique peak export
- `diffbind_peak_ext` (default `narrowPeak`): peak suffix for auto mode

## Note

If `Condition` and `Batch` are perfectly confounded, differential results are exploratory.
