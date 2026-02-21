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
- `Batch` (used for confounding warning only)

Example:

```csv
SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller
WT_GAR0968,WT,1,/path/WT_rep1.clean.bam,/path/WT_rep1_peaks.narrowPeak,macs
WT_GAR0979,WT,2,/path/WT_rep2.clean.bam,/path/WT_rep2_peaks.narrowPeak,macs
TG_GAR1585,TG,1,/path/TG_rep1.clean.bam,/path/TG_rep1_peaks.narrowPeak,macs
TG_GAR1586,TG,2,/path/TG_rep2.clean.bam,/path/TG_rep2_peaks.narrowPeak,macs
```

## Run

```bash
nextflow run main.nf -profile hpc --samplesheet samplesheet.csv
```

Resume:

```bash
nextflow run main.nf -profile hpc --samplesheet samplesheet.csv -resume
```

## Output

Output directory: `${project_folder}/${diffbind_output}`

Common files:

- `01_general_QC.pdf`
- `sample_info.after_count.tsv`
- `libsizes.tsv`, `libsizes.xlsx`
- `consensus_peaks.tsv`, `consensus_peaks.xlsx`
- `db_after_count_and_normalize.RData`
- `diffbind_session.RData`

If differential contrasts are possible (>=2 replicates per group):

- `02_<group1>.vs.<group2>.pdf`
- `significant.<group1>.vs.<group2>.tsv/.xlsx`
- `all_peaks.<group1>.vs.<group2>.tsv/.xlsx/.bed`
- `normalized_expression_of_significant_peaks.<group1>.vs.<group2>.tsv/.xlsx`

## Note

If `Condition` and `Batch` are perfectly confounded, DiffBind cannot separate batch from biology. In that case, results should be treated as exploratory.
