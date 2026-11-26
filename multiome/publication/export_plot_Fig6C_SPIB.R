# Figure 6C: SPIB SCT expression and chromVAR motif accessibility.
# Run from CH_publication_release (the repository root).
# Historical plotting recipe: multiome/scripts/2_multiome_processing1.R,
# lines 171-196. The plot uses sct_SPIB and MA0081.2 on wnn.umap.
#
# First pass, with the already-loaded checkpoint:
#   source("multiome/publication/export_plot_Fig6C_SPIB.R")
#   Fig6C_result <- reproduce_Fig6C_SPIB(obj = obj)
#
# For a later standalone runner (do NOT execute yet):
#   Fig6C_result <- reproduce_Fig6C_SPIB()
#
# Does not normalize, recompute chromVAR, modify the Seurat object,
# or write to historical scripts/, figures/, or data/ directories.
# Refuses to overwrite existing publication outputs.

reproduce_Fig6C_SPIB <- function(obj = NULL, root = getwd(),
                                 output_root = NULL) {
  needed <- c("Seurat", "SeuratObject", "Signac", "ggplot2")
  unavailable <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Unavailable package(s): ",
                                paste(unavailable, collapse = ", "), call. = FALSE)

  root <- normalizePath(root, mustWork = TRUE)
  input <- file.path(root, "multiome/data/processed/multiome_seurat.rds")
  if (!file.exists(input)) stop("Missing checkpoint; run from CH_publication_release: ", input)
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  historical <- normalizePath(file.path(root, "multiome/figures"), mustWork = TRUE)
  if (identical(output_root, historical) ||
      startsWith(output_root, paste0(historical, "/"))) {
    stop("Refusing to write inside historical multiome/figures.")
  }

  src_dir <- file.path(output_root, "source_data")
  fig_dir <- file.path(output_root, "figures")
  src_file <- file.path(src_dir, "Fig6C_SPIB_RNA_chromVAR.tsv.gz")
  fig_file <- file.path(fig_dir, "Fig6C_SPIB_reproduced_REVIEW.pdf")
  if (any(file.exists(c(src_file, fig_file)))) {
    stop("Existing outputs will not be overwritten:\n",
         paste(c(src_file, fig_file)[file.exists(c(src_file, fig_file))],
               collapse = "\n"), call. = FALSE)
  }
  if (is.null(obj)) {
    message("Loading historical Seurat checkpoint: ", input)
    obj <- readRDS(input)
  } else {
    message("Using supplied Seurat object; checkpoint not reloaded.")
  }
  if (!inherits(obj, "Seurat")) stop("obj is not a Seurat object.")
  if (ncol(obj) != 39787L) stop("Unexpected number of cells: ", ncol(obj))
  for (assay in c("SCT", "chromvar")) {
    if (!assay %in% SeuratObject::Assays(obj)) stop("Missing assay: ", assay)
  }
  if (!"wnn.umap" %in% names(obj@reductions)) stop("Missing wnn.umap reduction.")

  motif_id <- Signac::ConvertMotifID(obj, name = "SPIB")
  if (length(motif_id) != 1L || !identical(unname(motif_id), "MA0081.2")) {
    stop("Unexpected SPIB motif: ", paste(motif_id, collapse = ", "))
  }

  xy <- SeuratObject::Embeddings(obj, reduction = "wnn.umap")
  ids <- rownames(xy)
  if (!is.matrix(xy) || ncol(xy) != 2L || length(ids) != 39787L ||
      anyNA(ids) || anyDuplicated(ids) || any(!is.finite(xy))) {
    stop("Unexpected WNN embedding.")
  }
  if (!setequal(ids, colnames(obj))) stop("WNN cells differ from Seurat cells.")

  # Read one row per assay; do not materialize a dense 39,787-cell assay matrix.
  sct_data <- SeuratObject::GetAssayData(obj, assay = "SCT", layer = "data")
  chromvar_data <- SeuratObject::GetAssayData(obj, assay = "chromvar", layer = "data")
  if (!"SPIB" %in% rownames(sct_data)) stop("SPIB absent from SCT data layer.")
  if (!motif_id %in% rownames(chromvar_data)) stop("SPIB motif absent from chromvar data layer.")
  rna <- as.numeric(sct_data["SPIB", ids])
  motif <- as.numeric(chromvar_data[motif_id, ids])
  rm(sct_data, chromvar_data)
  if (length(rna) != length(ids) || length(motif) != length(ids) ||
      any(!is.finite(rna[!is.na(rna)])) ||
      any(!is.finite(motif[!is.na(motif)])) ||
      !any(is.finite(rna)) || !any(is.finite(motif))) {
    stop("Invalid RNA or chromVAR values.")
  }

  # Reuse the previously validated embedding export and its cell metadata.
  wnn_file <- file.path(root, "multiome/publication/source_data",
                        "Fig6B_FigS9B_WNN_UMAP.tsv.gz")
  if (!file.exists(wnn_file)) stop("Missing validated embedding Source Data: ", wnn_file)
  con <- gzfile(wnn_file, "rt")
  wnn <- tryCatch(utils::read.delim(con, check.names = FALSE,
                                    stringsAsFactors = FALSE),
                  finally = close(con))
  req <- c("cell_id", "UMAP_1", "UMAP_2", "sample", "transferred_cluster")
  if (!identical(names(wnn), req) || nrow(wnn) != length(ids) ||
      anyDuplicated(wnn$cell_id) || !setequal(wnn$cell_id, ids)) {
    stop("The previously validated WNN table has unexpected cells or fields.")
  }
  ix <- match(ids, wnn$cell_id)
  if (!isTRUE(all.equal(unname(xy[, 1]), wnn$UMAP_1[ix], tolerance = 1e-14)) ||
      !isTRUE(all.equal(unname(xy[, 2]), wnn$UMAP_2[ix], tolerance = 1e-14))) {
    stop("WNN coordinates disagree with the earlier validated Source Data.")
  }
  metadata <- obj[[]]
  metadata <- metadata[ids, , drop = FALSE]
  sample_now <- as.character(metadata$dataset)
  cluster_now <- as.character(metadata$projection)
  sample_old <- wnn$sample[ix]
  cluster_old <- as.character(wnn$transferred_cluster[ix])
  same_cluster <- (is.na(cluster_now) & is.na(cluster_old)) |
    (!is.na(cluster_now) & !is.na(cluster_old) & cluster_now == cluster_old)
  if (!identical(sample_now, sample_old) || !all(same_cluster)) {
    stop("Cell sample/cluster metadata differ from validated WNN Source Data.")
  }

  # Export RAW chromVAR scores, and also the display cutoff applied by FeaturePlot.
  # The historical plotting script sets min.cutoff=0 for the chromVAR panel.
  tab <- data.frame(
    cell_id = ids,
    UMAP_1 = unname(xy[, 1]),
    UMAP_2 = unname(xy[, 2]),
    sample = sample_now,
    transferred_cluster = cluster_now,
    SPIB_SCT_data = rna,
    SPIB_chromVAR_raw = motif,
    SPIB_chromVAR_display_min0 = pmax(motif, 0),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  for (directory in c(src_dir, fig_dir)) {
    if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
      stop("Cannot create: ", directory)
    }
  }
  tmp_tsv <- tempfile(pattern = ".fig6c-", tmpdir = src_dir, fileext = ".tsv.gz")
  tmp_pdf <- tempfile(pattern = ".fig6c-", tmpdir = fig_dir, fileext = ".pdf")
  on.exit(unlink(c(tmp_tsv, tmp_pdf)), add = TRUE)
  to_write <- tab
  num_cols <- vapply(to_write, is.numeric, logical(1))
  to_write[num_cols] <- lapply(to_write[num_cols], function(z) {
    ans <- rep(NA_character_, length(z))
    good <- !is.na(z)
    ans[good] <- formatC(z[good], format = "g", digits = 17)
    ans
  })
  con <- gzfile(tmp_tsv, "wt")
  tryCatch(utils::write.table(to_write, con, sep = "\t", quote = FALSE,
                              row.names = FALSE, na = "NA"),
           finally = close(con))
  con <- gzfile(tmp_tsv, "rt")
  back <- tryCatch(utils::read.delim(con, stringsAsFactors = FALSE,
                                     check.names = FALSE, na.strings = "NA"),
                   finally = close(con))
  if (!identical(back$cell_id, tab$cell_id) ||
      !identical(back$sample, tab$sample) ||
      !identical(as.character(back$transferred_cluster), tab$transferred_cluster)) {
    stop("Failed Source Data identifier/annotation read-back check.")
  }
  compare_numeric <- function(a, b) {
    identical(is.na(a), is.na(b)) &&
      isTRUE(all.equal(a[!is.na(a)], b[!is.na(b)], tolerance = 1e-14))
  }
  if (!all(vapply(names(tab)[num_cols], function(nm) {
    compare_numeric(tab[[nm]], back[[nm]])
  }, logical(1)))) stop("Failed Source Data numerical read-back check.")

  # Recreate the historical FeaturePlot calls directly from the saved checkpoint.
  # Sourcing the historical 2_multiome_processing1.R would run unrelated analyses.
  gene_plot <- Seurat::FeaturePlot(
    obj, features = "sct_SPIB", order = TRUE, reduction = "wnn.umap"
  )
  motif_plot <- Seurat::FeaturePlot(
    obj, features = motif_id, min.cutoff = 0,
    cols = c("lightgrey", "darkred"), order = TRUE,
    reduction = "wnn.umap"
  )
  grDevices::pdf(tmp_pdf, width = 8, height = 4, useDingbats = FALSE)
  tryCatch({
    if (requireNamespace("patchwork", quietly = TRUE)) {
      print(patchwork::wrap_plots(gene_plot, motif_plot, ncol = 2))
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(gene_plot, motif_plot, ncol = 2)
    } else {
      stop("Either patchwork or gridExtra is needed to assemble the two panels.")
    }
  }, finally = grDevices::dev.off())
  if (!file.exists(tmp_pdf) || file.info(tmp_pdf)$size < 5000L) {
    stop("PDF creation failed or resulted in an unexpectedly small file.")
  }

  if (any(file.exists(c(src_file, fig_file)))) stop("Output files appeared during processing.")
  if (!file.rename(tmp_tsv, src_file)) stop("Could not finalize Source Data.")
  if (!file.rename(tmp_pdf, fig_file)) stop("Could not finalize figure PDF.")

  summary <- data.frame(
    panel = "Fig6C", cells = length(ids),
    rna_feature = "sct_SPIB", motif_id = motif_id,
    rna_na = sum(is.na(rna)), motif_na = sum(is.na(motif)),
    motif_negative = sum(motif < 0, na.rm = TRUE),
    rna_min = min(rna, na.rm = TRUE), rna_max = max(rna, na.rm = TRUE),
    motif_min = min(motif, na.rm = TRUE),
    motif_max = max(motif, na.rm = TRUE),
    source_data_roundtrip = "PASS", previous_wnn_match = "PASS",
    graphic_written = file.info(fig_file)$size > 5000L,
    visual_agreement_with_artwork = "REVIEW"
  )
  print(summary, row.names = FALSE)
  message("Source Data: ", src_file)
  message("Figure for review: ", fig_file)
  invisible(summary)
}
