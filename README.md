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

WT_1,WT,1,/path/WT.Rep1.bam,/path/WT.input.bam,/path/WT.Rep1_peaks.xls,macs
TG_1,TG,1,/path/TG.Rep1.bam,/path/TG.input.bam,/path/TG.Rep1_peaks.xls,macs

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
