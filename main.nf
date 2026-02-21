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

  fdr_cutoff <- as.numeric("${params.diff_fdr ?: 0.05}")
  lfc_cutoff <- as.numeric("${params.diff_lfc ?: 1}")

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

  ## Peak-universe overlap by condition from input peak files
  read_peak_ids <- function(path) {
    D <- tryCatch(read.delim(path, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(D) || ncol(D) < 3) return(character(0))
    paste(D[[1]], D[[2]], D[[3]], sep = ":")
  }

  cond_levels <- unique(samples[["Condition"]])
  peak_sets <- list()
  for (cond in cond_levels) {
    rows <- samples[samples[["Condition"]] == cond, , drop = FALSE]
    ids <- unique(unlist(lapply(rows[["Peaks"]], read_peak_ids), use.names = FALSE))
    peak_sets[[as.character(cond)]] <- ids
  }

  if (length(peak_sets) >= 2) {
    universe <- unique(unlist(peak_sets, use.names = FALSE))
    mem <- data.frame(peak_id = universe, stringsAsFactors = FALSE)
    for (nm in names(peak_sets)) {
      mem[[nm]] <- as.integer(universe %in% peak_sets[[nm]])
    }
    write.table(mem, "peak_universe_upset_input.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

    overlap_stats <- data.frame(
      condition = names(peak_sets),
      n_peaks = vapply(peak_sets, length, integer(1)),
      stringsAsFactors = FALSE
    )
    write.table(overlap_stats, "peak_universe_condition_sizes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

    # Pairwise overlap matrix
    conds <- names(peak_sets)
    M <- matrix(0, nrow = length(conds), ncol = length(conds), dimnames = list(conds, conds))
    for (i in seq_along(conds)) {
      for (j in seq_along(conds)) {
        M[i, j] <- length(intersect(peak_sets[[conds[i]]], peak_sets[[conds[j]]]))
      }
    }
    write.table(M, "peak_universe_pairwise_overlap.tsv", sep = "\t", quote = FALSE, col.names = NA)

    # Optional venn if exactly two groups and package is available
    if (length(peak_sets) == 2 && requireNamespace("VennDiagram", quietly = TRUE)) {
      grDevices::pdf("peak_universe_venn.pdf", width = 6, height = 6)
      grid::grid.newpage()
      VennDiagram::draw.pairwise.venn(
        area1 = length(peak_sets[[1]]),
        area2 = length(peak_sets[[2]]),
        cross.area = length(intersect(peak_sets[[1]], peak_sets[[2]])),
        category = names(peak_sets),
        fill = c("#4C78A8", "#F58518"),
        alpha = c(0.5, 0.5),
        cex = 1.2,
        cat.cex = 1.1
      )
      grDevices::dev.off()
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
  summary_rows <- list()

  if (!has_reps) {
    message("Not enough replicates in more than one group; skipping differential contrast analysis.")
  } else {
    db <- dba.contrast(db, categories = DBA_CONDITION, minMembers = 2)
    db <- dba.analyze(db, method = DBA_DESEQ2)

    comparisons <- dba.show(db, bContrasts = TRUE)
    write.table(as.data.frame(comparisons), "contrasts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

    if (nrow(comparisons) > 0) {
      for (comp in seq_len(nrow(comparisons))) {
        raw_comp_name <- paste0(comparisons[["Group"]][comp], ".vs.", comparisons[["Group2"]][comp])
        comp_name <- gsub("[^A-Za-z0-9._-]", "_", raw_comp_name)
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

        n_all <- nrow(all_peaks)
        n_sig <- NA_integer_
        n_up <- NA_integer_
        n_down <- NA_integer_

        if (n_all > 0 && all(c("FDR", "Fold", "seqnames", "start", "end") %in% colnames(all_peaks))) {
          sig2 <- subset(all_peaks, !is.na(FDR) & FDR <= fdr_cutoff & !is.na(Fold) & abs(Fold) >= lfc_cutoff)
          up <- subset(sig2, Fold >= lfc_cutoff)
          down <- subset(sig2, Fold <= -lfc_cutoff)

          n_sig <- nrow(sig2)
          n_up <- nrow(up)
          n_down <- nrow(down)

          if (n_up > 0) {
            write.table(up[, c("seqnames", "start", "end")],
                        paste0("condition_unique_up.", comp_name, ".bed"),
                        sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
          }
          if (n_down > 0) {
            write.table(down[, c("seqnames", "start", "end")],
                        paste0("condition_unique_down.", comp_name, ".bed"),
                        sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
          }
        }

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          contrast = comp_name,
          n_all = n_all,
          n_sig_fdr_lfc = n_sig,
          n_up = n_up,
          n_down = n_down,
          stringsAsFactors = FALSE
        )

        dev.off()
      }
    }
  }

  if (length(summary_rows) > 0) {
    S <- do.call(rbind, summary_rows)
    write.table(S, "diffbind_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    openxlsx::write.xlsx(S, "diffbind_summary.xlsx", rowNames = FALSE)
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
