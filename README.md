# nf-diffbind

Nextflow DSL2 module for DiffBind differential binding analysis.

## Input

Required `samplesheet.csv` columns:

- `SampleID`
- `Condition`
- `Replicate`
- `bamReads`
- `Peaks`
- `PeakCaller`

Optional:

- `bamControl`
- `Batch` (used for confounding warning)

## Run

```bash
nextflow run main.nf -profile hpc --samplesheet samplesheet.csv
```

Resume:

```bash
nextflow run main.nf -profile hpc --samplesheet samplesheet.csv -resume
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

## Note

If `Condition` and `Batch` are perfectly confounded, differential results are exploratory.
