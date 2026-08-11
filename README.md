# nf-diffbind

`nf-diffbind` runs DiffBind differential binding analysis from filtered BAMs and MACS3 peaks.

This module is optional after `nf-macs3`. Use it when there are enough replicates per condition for differential binding analysis.

## Input Modes

1. Explicit DiffBind samplesheet with `--samplesheet`
2. Auto-generate from `samples_master`

## Auto Mode

Required:

```bash
--samples_master /path/to/samples_master.csv
--chipfilter_output /path/to/chipfilter_output
--macs3_output /path/to/macs3_output
```

Required `samples_master` columns:

```csv
sample_id,condition
```

Useful optional columns:

```csv
replicate,library_type,is_control,use_for_diffbind,enabled
```

Auto mode uses:

```text
${chipfilter_output}/${sample_id}*.nomulti.bam
${macs3_output}/${diffbind_macs3_profile}/${sample_id}_peaks.${diffbind_peak_ext}
```

Default MACS3 profile:

```text
strict_q0.01
```

## Explicit Samplesheet

Required columns:

```csv
SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller
```

Optional:

```csv
bamControl,Batch
```

## Recommended Run

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --macs3_output /path/to/macs3_output \
  --project_folder /path/to/output_project \
  --diffbind_output diffbind_output
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --macs3_output /path/to/macs3_output \
  --project_folder /path/to/output_project
```

## Key Outputs

Output directory:

```text
${project_folder}/${diffbind_output}
```

Important outputs:

- `01_general_QC.pdf`
- `sample_info.after_count.tsv`
- `libsizes.tsv` / `libsizes.xlsx`
- `contrasts.tsv`
- `significant.<contrast>.tsv/.xlsx`
- `all_peaks.<contrast>.tsv/.xlsx/.bed`
- `condition_unique_up.<contrast>.bed`
- `condition_unique_down.<contrast>.bed`
- `diffbind_summary.tsv/.xlsx`
- `diffbind_session.RData`

## Key Parameters

| Parameter | Meaning | Default |
| --- | --- | --- |
| `samplesheet` | Existing DiffBind samplesheet | `null` |
| `samples_master` | Metadata CSV for auto mode | `null` |
| `chipfilter_output` | Folder with `*.nomulti.bam` | `null` |
| `macs3_output` | MACS3 output folder | `null` |
| `diffbind_macs3_profile` | MACS3 profile for peaks | `strict_q0.01` |
| `diff_fdr` | FDR cutoff for exported significant peaks | `0.05` |
| `diff_lfc` | Absolute fold-change cutoff | `1` |
| `mail_user` | HPC notification email | `molendo.hpc@gmail.com` |

## Notes

- Differential results need biological replication. If groups lack replicates, the module writes QC/summary outputs and skips contrast analysis.
- If `Condition` and `Batch` are perfectly confounded, results should be treated as exploratory.
- Actual execution should be tested on an environment with Nextflow available.
