# Fig 6I: reproduce genotype-associated TF motif scatter from archived results.
# Root: CH_publication_release. Original scripts, workbooks and figures are read-only.
# Workbook: ../pipeline/resource/motif_associations_genotypes.xlsx
# No genotype DA refitting; motif GC values use the formula in 5_clone_analysis_final.R.

reproduce_Fig6I_from_workbook <- function(obj = NULL, root = getwd(),
                                         workbook = NULL, output_root = NULL) {
  needed <- c("openxlsx", "ggplot2", "ggrepel", "scales", "Matrix", "Signac", "SeuratObject")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Packages unavailable: ", paste(missing, collapse = ", "), call. = FALSE)

  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome")))
    stop("Working directory must be CH_publication_release.", call. = FALSE)
  if (is.null(workbook))
    workbook <- file.path(root, "pipeline/resource/motif_associations_genotypes.xlsx")
  if (!file.exists(workbook)) stop("Missing historical workbook: ", workbook, call. = FALSE)
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  historic <- normalizePath(file.path(root, "multiome/figures"), mustWork = TRUE)
  if (identical(output_root, historic) || startsWith(output_root, paste0(historic, "/")))
    stop("Refusing to write inside historical figure directory.", call. = FALSE)

  src <- file.path(output_root, "source_data/Fig6I_motif_genotype_associations.tsv.gz")
  fig <- file.path(output_root, "figures/Fig6I_motif_genotypes_REVIEW.pdf")
  existing <- c(src, fig)[file.exists(c(src, fig))]
  if (length(existing)) stop("Existing publication files protected:\n",
                             paste(existing, collapse = "\n"), call. = FALSE)

  sheets <- openxlsx::getSheetNames(workbook)
  if (!all(c("DNMT3A", "TET2") %in% sheets))
    stop("Expected DNMT3A and TET2 sheets.", call. = FALSE)
  fields <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  read_sheet <- function(s) {
    dat <- openxlsx::read.xlsx(workbook, sheet = s, colNames = TRUE)
    if (!all(fields %in% colnames(dat))) stop("Missing numeric columns in ", s)
    motif <- trimws(as.character(dat[[1]])) # unnamed leading Excel column
    if (anyNA(motif) || any(!nzchar(motif)) || anyDuplicated(motif))
      stop("Missing or repeated motif names in ", s,
           "; identifiers must be resolved before a safe genotype join.")
    if (!all(vapply(dat[, fields, drop = FALSE], is.numeric, logical(1))))
      stop("Non-numeric statistical results in ", s)
    out <- dat[, fields, drop = FALSE]
    rownames(out) <- motif
    out
  }
  dn <- read_sheet("DNMT3A")
  te <- read_sheet("TET2")
  if (!setequal(rownames(dn), rownames(te)))
    stop("Motif sets differ between DNMT3A and TET2 sheets.", call. = FALSE)
  motifs <- rownames(dn)
  te <- te[motifs, , drop = FALSE]

  if (is.null(obj)) {
    checkpoint <- file.path(root, "multiome/data/processed/multiome_seurat.rds")
    message("Loading Multiome checkpoint for historical motif GC-content: ", checkpoint)
    obj <- readRDS(checkpoint)
  } else {
    message("Using supplied Multiome Seurat object; no checkpoint reload.")
  }
  if (!inherits(obj, "Seurat") || !"ATAC" %in% SeuratObject::Assays(obj))
    stop("Expected a Multiome Seurat object with an ATAC assay.")
  atac <- obj[["ATAC"]]
  mm <- atac@motifs@data
  motif_names <- unlist(atac@motifs@motif.names)
  if (is.null(colnames(mm)) || !all(colnames(mm) %in% names(motif_names)))
    stop("Motif matrix columns do not align to motif annotations.")
  named_motifs <- unname(as.character(motif_names[colnames(mm)]))
  if (anyDuplicated(named_motifs))
    stop("Motif annotation contains duplicated motif names; manual ID mapping required.")
  mf <- atac@meta.features
  if (!"GC.percent" %in% colnames(mf)) stop("ATAC peak GC.percent is unavailable.")
  if (!identical(rownames(mm), rownames(mf))) {
    if (!all(rownames(mm) %in% rownames(mf))) stop("Motif/peak annotation rows do not match.")
    mf <- mf[rownames(mm), , drop = FALSE]
  }
  gc <- as.numeric(mf$GC.percent)
  if (any(!is.finite(gc))) stop("Missing/nonfinite peak GC.percent values.")

  # Exact historical motif.gc formula, 5_clone_analysis_final.R lines 64-70:
  # remove chrM/X/Y from numerator; denominator remains all peaks.
  contig <- sub("-.*", "", rownames(mm))
  keep <- !(contig %in% c("chrM", "chrX", "chrY"))
  numerator <- as.numeric(Matrix::crossprod(mm[keep, , drop = FALSE], gc[keep]))
  denominator <- as.numeric(Matrix::colSums(mm))
  motif_gc <- numerator / denominator
  names(motif_gc) <- named_motifs
  if (any(!is.finite(motif_gc)) || !all(motifs %in% names(motif_gc)))
    stop("Cannot align finite historical motif GC contents to workbook names.")

  out <- data.frame(motif = motifs, motif_gc_percent = unname(motif_gc[motifs]),
                    stringsAsFactors = FALSE)
  for (nm in fields) {
    out[[paste0("DNMT3A_", nm)]] <- as.numeric(dn[[nm]])
    out[[paste0("TET2_", nm)]] <- as.numeric(te[[nm]])
  }
  stopifnot(nrow(out) == length(motifs), !anyDuplicated(out$motif),
            all(is.finite(out$DNMT3A_log2FoldChange)),
            all(is.finite(out$TET2_log2FoldChange)))

  # Historical plot: x=DNMT3A logFC, y=TET2 logFC, color=motif mean GC.
  gc_values <- as.numeric(motif_gc)
  p <- ggplot2::ggplot(out, ggplot2::aes(
    x = DNMT3A_log2FoldChange, y = TET2_log2FoldChange,
    label = motif, colour = motif_gc_percent)) +
    ggplot2::geom_point(size = 1.35) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggrepel::geom_text_repel(size = 3.1, max.overlaps = 10, seed = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_gradientn(
      colours = c("blue", "grey70", "darkred"),
      values = scales::rescale(c(min(gc_values), mean(gc_values), max(gc_values))),
      name = "motif GC") +
    ggplot2::labs(x = "logFC, DNMT3A", y = "logFC, TET2", title = "DNMT3A vs TET2") +
    ggplot2::theme(text = ggplot2::element_text(size = 13))

  dir.create(dirname(src), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(fig), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(src, "wt")
  tryCatch(utils::write.table(out, con, sep = "\t", quote = FALSE,
                              row.names = FALSE, col.names = TRUE, na = "NA"),
           finally = close(con))
  check <- utils::read.delim(src, check.names = FALSE)
  if (!identical(out$motif, check$motif) ||
      !isTRUE(all.equal(out$DNMT3A_log2FoldChange, check$DNMT3A_log2FoldChange,
                        tolerance = 1e-12)) ||
      !isTRUE(all.equal(out$TET2_log2FoldChange, check$TET2_log2FoldChange,
                        tolerance = 1e-12)) ||
      !isTRUE(all.equal(out$motif_gc_percent, check$motif_gc_percent,
                        tolerance = 1e-12)))
    stop("Source Data round-trip check failed.")
  grDevices::pdf(fig, width = 8, height = 7, useDingbats = FALSE)
  tryCatch(print(p), finally = grDevices::dev.off())
  summary <- data.frame(panel = "Fig6I", motifs = nrow(out),
                        source_data_roundtrip = "PASS", visual_comparison = "REVIEW")
  print(summary, row.names = FALSE)
  cat("Source Data:", src, "\nReview figure:", fig, "\n")
  invisible(summary)
}
