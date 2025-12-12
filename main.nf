process diffbind_run {

  tag "DiffBind"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${diffbind_output}", mode: 'copy'

  input:
    path samplesheet

  output:
    path "*"

  script:
  """
  set -euo pipefail

  # 确保 R 有可写的临时目录
  export TMPDIR=\$PWD

  cat > diffbind_run.R << 'EOF'
  #!/usr/bin/env Rscript

  suppressPackageStartupMessages({
    library(DiffBind)
    library(openxlsx)
  })

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: diffbind_run.R samplesheet.csv", call. = FALSE)
  }

  samplesheet <- normalizePath(args[1])
  message(">>> Using sampleSheet: ", samplesheet)

  outdir <- getwd()
  message(">>> Working directory: ", outdir)

  samples <- read.csv(samplesheet, stringsAsFactors = FALSE)

  required_cols <- c("SampleID","Condition","Replicate","bamReads","Peaks","PeakCaller")
  if (!all(required_cols %in% colnames(samples))) {
    stop(paste0("Sample sheet must contain columns: ",
                paste(required_cols, collapse = ", ")))
  }

  if (!"bamControl" %in% colnames(samples)) {
    samples[["bamControl"]] <- NA
  }

  message(">>> Checking peak files...")
  samples[["nPeaks"]] <- NA_integer_

  for (i in seq_len(nrow(samples))) {
    peak_file <- samples[["Peaks"]][i]
    if (!file.exists(peak_file)) {
      warning("Peak file not found: ", peak_file)
      samples[["nPeaks"]][i] <- 0L
    } else {
      n <- length(readLines(peak_file, warn = FALSE))
      samples[["nPeaks"]][i] <- n
    }
  }

  message("Peak counts per sample:")
  print(samples[, c("SampleID","Condition","nPeaks")])

  samples <- subset(samples, nPeaks > 0)
  if (nrow(samples) < 2) {
    stop("Less than 2 samples with peaks > 0. DiffBind will not run.")
  }

  write.csv(samples, "samples_filtered.csv", row.names = FALSE)

  message(">>> Building DiffBind object...")
  db <- dba(sampleSheet = samples)

  pdf("01_general_QC.pdf")
  plot(db)

  db <- dba.count(db, bParallel = FALSE)
  plot(db)

  info <- dba.show(db)
  print(info)

  if ("FRiP" %in% colnames(info)) {
    libsizes <- cbind(
      LibReads  = info[, "Reads"],
      FRiP      = info[, "FRiP"],
      PeakReads = round(info[, "Reads"] * info[, "FRiP"])
    )
    rownames(libsizes) <- info[, "ID"]
    write.xlsx(as.data.frame(libsizes), "libsizes.xlsx", row.names = TRUE)
  }

  try(dba.plotPCA(db, attributes = DBA_CONDITION, label = DBA_ID), silent = TRUE)

  db <- dba.normalize(db)
  plot(db)

  dev.off()

  save(db, file = "db_after_count_and_normalize.RData")

  group_col <- "Condition"

  if (sum(table(samples[[group_col]]) > 1) > 1) {

    message(">>> Setting up contrasts by ", group_col, " ...")

    db <- dba.contrast(db, categories = DBA_CONDITION, minMembers = 2)
    db <- dba.analyze(db, method = DBA_DESEQ2)

    comparisons <- dba.show(db, bContrasts = TRUE)
    message(">>> Contrasts:")
    print(comparisons)

    for (comp in seq_len(nrow(comparisons))) {

      comp_name <- paste0(comparisons$Group[comp], ".vs.", comparisons$Group2[comp])
      message(">>> Processing contrast: ", comp_name)

      pdf(paste0("02_", comp_name, ".pdf"))

      sig <- dba.report(db, contrast = comp)

      if (!is.null(sig) && length(sig@ranges) > 1) {

        try(dba.plotPCA(db, contrast = comp, method = DBA_DESEQ2,
                        attributes = DBA_CONDITION, label = DBA_ID),
            silent = TRUE)

        plot(db, contrast = comp)

        try(dba.plotVolcano(db, contrast = comp), silent = TRUE)
        try(dba.plotBox(db, contrast = comp), silent = TRUE)

        try({
          corvals <- dba.plotHeatmap(db, correlations = FALSE,
                                     scale = "row", contrast = comp)
          write.table(corvals,
                      paste0("normalized_expression_of_significant_peaks.", comp_name, ".tsv"),
                      sep = "\t", quote = FALSE, row.names = FALSE)
          write.xlsx(as.data.frame(corvals),
                     paste0("normalized_expression_of_significant_peaks.", comp_name, ".xlsx"),
                     row.names = FALSE)
        }, silent = TRUE)

        sig_df <- as.data.frame(sig)
        write.table(sig_df,
                    paste0("significant.", comp_name, ".tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        write.xlsx(sig_df,
                   paste0("significant.", comp_name, ".xlsx"),
                   row.names = FALSE)

      } else {
        message("No significant peaks in ", comp_name)
      }

      try(dba.plotMA(db, method = DBA_DESEQ2, contrast = comp), silent = TRUE)
      try(dba.plotMA(db, bXY = TRUE, contrast = comp), silent = TRUE)

      all_peaks <- as.data.frame(dba.report(db, th = 1, contrast = comp))
      write.table(all_peaks,
                  paste0("all_peaks.", comp_name, ".tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      write.xlsx(all_peaks,
                 paste0("all_peaks.", comp_name, ".xlsx"))
      write.table(all_peaks[, 1:3],
                  paste0("all_peaks.", comp_name, ".bed"),
                  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

      dev.off()
    }

  } else {
    message(">>> Not enough replicates in more than one group for differential analysis.")
  }

  message(">>> Writing consensus binding matrix...")

  consensus <- as.data.frame(db[["binding"]])
  write.xlsx(consensus, "consensus_peaks.xlsx", row.names = FALSE)
  write.table(consensus, "consensus_peaks.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  save.image("diffbind_session.RData")

  message(">>> Done.")
  sessionInfo()
  EOF

  Rscript diffbind_run.R ${samplesheet}
  """
}
