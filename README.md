# nf-diffbind

A clean and modular DiffBind pipeline for differential binding analysis of ChIP-seq data.

---

## 📌 Input

### Required:

- `samplesheet.csv`  
  Must contain:

  | SampleID | Condition | Replicate | bamReads | bamControl | Peaks | PeakCaller |
  |----------|-----------|-----------|----------|------------|--------|------------|

Example:

SampleID,Condition,Replicate,bamReads,bamControl,Peaks,PeakCaller
WT_1,WT,1,/path/bam/WT.Rep1.md.bam,/path/input/WT_input.Rep1.md.bam,/path/peaks/WT.Rep1_peaks.xls,macs
WT_2,WT,2,/path/bam/WT.Rep2.md.bam,/path/input/WT_input.Rep2.md.bam,/path/peaks/WT.Rep2_peaks.xls,macs
TG_1,TG,1,/path/bam/TG.Rep1.md.bam,/path/input/TG_input.Rep1.md.bam,/path/peaks/TG.Rep1_peaks.xls,macs
TG_2,TG,2,/path/bam/TG.Rep2.md.bam,/path/input/TG_input.Rep2.md.bam,/path/peaks/TG.Rep2_peaks.xls,macs


---

## 🚀 Run
```
nextflow run main.nf -profile hpc --samplesheet samplesheet.csv --outdir diffbind_output
```

---

## 📤 Output structure

diffbind_output/  
01_QC.pdf  
02_PCA.pdf  
significant_WT_vs_TG.xlsx  
all_peaks_WT_vs_TG.bed  
consensus_peaks.tsv  
db_after_count.RData  
db_normalized.RData  


## 🧠 Notes

- DiffBind analysis method = DESeq2    
- Condition grouping uses column `Condition`    
- Requires ≥ 2 replicates per group  
- Compatible with nf-macs3 output  
