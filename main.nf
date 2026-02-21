#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def diffbind_output = params.diffbind_output ?: 'diffbind_output'

process diffbind_run {
  tag "diffbind"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${diffbind_output}", mode: 'copy', overwrite: true

  input:
    path samplesheet

  output:
    path "*"

  script:
  """
  set -euo pipefail

  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp

  cat > diffbind_run.R << 'RS'
  suppressPackageStartupMessages({
    library(DiffBind)
    library(openxlsx)
  })

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: diffbind_run.R samplesheet.csv", call. = FALSE)

  samplesheet <- normalizePath(args[1])
  message("Using sample sheet: ", samplesheet)

  samples <- read.csv(samplesheet, stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c("SampleID","Condition","Replicate","bamReads","Peaks","PeakCaller")
  miss <- setdiff(required_cols, colnames(samples))
  if (length(miss) > 0) {
    stop("Missing required columns in sample sheet: ", paste(miss, collapse = ", "))
  }

  if (!("bamControl" %in% colnames(samples))) samples[["bamControl"]] <- NA

  for (i in seq_len(nrow(samples))) {
    b <- samples[["bamReads"]][i]
    p <- samples[["Peaks"]][i]
    if (!file.exists(b)) stop("bamReads not found: ", b)
    if (!file.exists(p)) stop("Peaks not found: ", p)
  }

  message("Peak counts per sample")
  samples[["nPeaks"]] <- vapply(samples[["Peaks"]], function(p) length(readLines(p, warn = FALSE)), integer(1))
  print(samples[, c("SampleID","Condition","Replicate","nPeaks")])

  samples <- subset(samples, nPeaks > 0)
  if (nrow(samples) < 2) stop("Less than 2 samples with non-empty peaks.")

  write.csv(samples, "samples.filtered.csv", row.names = FALSE)

  # Basic confounding warning (if Batch exists)
  if ("Batch" %in% colnames(samples)) {
    tab <- table(samples[["Condition"]], samples[["Batch"]])
    if (all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)) {
      message("WARNING: Condition and Batch are perfectly confounded; differential results are exploratory.")
    }
  }

  db <- dba(sampleSheet = samples)

  pdf("01_general_QC.pdf")
  plot(db)

  db <- dba.count(db, bParallel = FALSE)
  plot(db)

  info <- dba.show(db)
  write.table(as.data.frame(info), "sample_info.after_count.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  if (all(c("Reads","FRiP","ID") %in% colnames(info))) {
    qc <- data.frame(
      SampleID = info[, "ID"],
      Reads = info[, "Reads"],
      FRiP = info[, "FRiP"],
      PeakReads = round(info[, "Reads"] * info[, "FRiP"])
    )
    write.table(qc, "libsizes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    openxlsx::write.xlsx(qc, "libsizes.xlsx", rowNames = FALSE)
  }

  try(dba.plotPCA(db, attributes = DBA_CONDITION, label = DBA_ID), silent = TRUE)

  db <- dba.normalize(db)
  plot(db)
  dev.off()

  save(db, file = "db_after_count_and_normalize.RData")

  has_reps <- sum(table(samples[["Condition"]]) > 1) > 1
  if (!has_reps) {
    message("Not enough replicates in more than one group; skipping differential contrast analysis.")
  } else {
    db <- dba.contrast(db, categories = DBA_CONDITION, minMembers = 2)
    db <- dba.analyze(db, method = DBA_DESEQ2)

    comparisons <- dba.show(db, bContrasts = TRUE)
    write.table(as.data.frame(comparisons), "contrasts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

    if (nrow(comparisons) > 0) {
      for (comp in seq_len(nrow(comparisons))) {
        comp_name <- paste0(comparisons[["Group"]][comp], ".vs.", comparisons[["Group2"]][comp])
        pdf(paste0("02_", comp_name, ".pdf"))

        sig <- dba.report(db, contrast = comp)

        if (!is.null(sig) && length(sig@ranges) > 1) {
          try(dba.plotPCA(db, contrast = comp, method = DBA_DESEQ2, attributes = DBA_CONDITION, label = DBA_ID), silent = TRUE)
          try(plot(db, contrast = comp), silent = TRUE)
          try(dba.plotVolcano(db, contrast = comp), silent = TRUE)
          try(dba.plotBox(db, contrast = comp), silent = TRUE)

          try({
            corvals <- dba.plotHeatmap(db, correlations = FALSE, scale = "row", contrast = comp)
            write.table(corvals,
                        paste0("normalized_expression_of_significant_peaks.", comp_name, ".tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            openxlsx::write.xlsx(as.data.frame(corvals),
                                 paste0("normalized_expression_of_significant_peaks.", comp_name, ".xlsx"),
                                 rowNames = FALSE)
          }, silent = TRUE)

          sig_df <- as.data.frame(sig)
          write.table(sig_df, paste0("significant.", comp_name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
          openxlsx::write.xlsx(sig_df, paste0("significant.", comp_name, ".xlsx"), rowNames = FALSE)
        }

        try(dba.plotMA(db, method = DBA_DESEQ2, contrast = comp), silent = TRUE)
        try(dba.plotMA(db, bXY = TRUE, contrast = comp), silent = TRUE)

        all_peaks <- as.data.frame(dba.report(db, th = 1, contrast = comp))
        write.table(all_peaks, paste0("all_peaks.", comp_name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        openxlsx::write.xlsx(all_peaks, paste0("all_peaks.", comp_name, ".xlsx"), rowNames = FALSE)

        if (ncol(all_peaks) >= 3) {
          write.table(all_peaks[, 1:3],
                      paste0("all_peaks.", comp_name, ".bed"),
                      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        }

        dev.off()
      }
    }
  }

  consensus <- as.data.frame(db[["binding"]])
  write.table(consensus, "consensus_peaks.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  openxlsx::write.xlsx(consensus, "consensus_peaks.xlsx", rowNames = FALSE)

  save.image("diffbind_session.RData")
  sessionInfo()
  RS

  Rscript diffbind_run.R ${samplesheet}
  """
}

workflow {
  if (!params.samplesheet) {
    exit 1, "ERROR: Please provide --samplesheet"
  }

  Channel
    .fromPath(params.samplesheet, checkIfExists: true)
    .ifEmpty { exit 1, "ERROR: samplesheet not found: ${params.samplesheet}" }
    .set { ch_samplesheet }

  diffbind_run(ch_samplesheet)
}
