# Multiome publication: Fig 6B, Fig S9B and Fig S9C.
# Run from CH_publication_release, NOT multiome/.
#
# First pass, reusing the checkpoint already loaded in this session:
#   source("multiome/publication/export_multiome_batch1_embeddings.R")
#   embedding_summary <- export_multiome_batch1_embeddings(obj = obj)
#
# Standalone entry for the eventual publication runner:
#   source("multiome/publication/export_multiome_batch1_embeddings.R")
#   embedding_summary <- export_multiome_batch1_embeddings()
#
# Sourcing defines the function only. The explicit call runs it.
# No checkpoint updates, normalization, integration, clustering or plotting.
# No changes to historical data/, scripts/ or figures/.
# Refuses to replace existing outputs.
#
# Historical sources inspected:
#   multiome/scripts/1_rna_integration.R: saves projection.clusters in
#     data/processed/clusters_annotation.rds.
#   multiome/scripts/2_multiome_processing1.R: plots dataset / projection;
#     preserves unassigned projection cells (grey80).
# The expected counts below are from the user's checkpoint inspection, NOT
# an independent numerical reference extracted from publication artwork.

export_multiome_batch1_embeddings <- function(obj = NULL, root = getwd(),
                                               output_root = NULL) {
  for (pkg in c("Seurat", "SeuratObject", "Signac")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Required package cannot load: ", pkg, call. = FALSE)
    }
  }
  root <- normalizePath(root, mustWork = TRUE)
  input_file <- file.path(root, "multiome/data/processed/multiome_seurat.rds")
  annotation_file <- file.path(root, "multiome/data/processed/clusters_annotation.rds")
  for (f in c(input_file, annotation_file)) {
    if (!file.exists(f)) {
      stop("Missing input: ", f, "\nRun from CH_publication_release.", call. = FALSE)
    }
  }
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  source_dir <- file.path(output_root, "source_data")
  validation_dir <- file.path(output_root, "validation")
  output_files <- c(
    "wnn.umap" = "Fig6B_FigS9B_WNN_UMAP.tsv.gz",
    "umap.rna" = "FigS9C_RNA_UMAP.tsv.gz",
    "umap.atac" = "FigS9C_ATAC_UMAP.tsv.gz"
  )
  panel_names <- c("Fig6B/FigS9B", "FigS9C_RNA", "FigS9C_ATAC")
  final_files <- file.path(source_dir, unname(output_files))
  report_file <- file.path(validation_dir, "Batch1_embeddings_export_summary.tsv")
  existing <- c(final_files, report_file)[file.exists(c(final_files, report_file))]
  if (length(existing)) {
    stop("Existing outputs retained; no files were overwritten:\n",
         paste(existing, collapse = "\n"), call. = FALSE)
  }
  if (is.null(obj)) {
    message("Reading historical checkpoint: ", input_file)
    obj <- readRDS(input_file)
    object_source <- "read_from_checkpoint"
  } else {
    message("Using explicitly supplied Seurat object; checkpoint is not reloaded.")
    object_source <- "explicit_object_argument"
  }
  if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.")
  meta <- obj[[]]
  ids <- colnames(obj)
  required <- c("dataset", "projection")
  if (!all(required %in% names(meta))) stop("Missing dataset/projection metadata.")
  if (length(ids) != 39787L || is.null(ids) || anyNA(ids) || any(ids == "") ||
      anyDuplicated(ids) || is.null(rownames(meta)) || anyDuplicated(rownames(meta)) ||
      !setequal(ids, rownames(meta))) stop("Unexpected or inconsistent checkpoint cell IDs.")
  meta <- meta[ids, , drop = FALSE]
  sample <- as.character(meta$dataset)
  cluster <- as.character(meta$projection)
  expected_samples <- c(C1 = 9690L, C2 = 6840L, D = 8799L, T = 14458L)
  expected_clusters <- c("1" = 4020L, "2" = 2115L, "3" = 1213L, "4" = 8760L,
                         "5" = 4728L, "6" = 7382L, "7" = 1820L, "8" = 1848L,
                         "9" = 348L, "10" = 1220L)
  if (anyNA(sample) || !setequal(unique(sample), names(expected_samples)) ||
      !identical(as.integer(table(sample)[names(expected_samples)]),
                 unname(expected_samples))) stop("Sample counts differ from inspected checkpoint.")
  if (!setequal(unique(cluster[!is.na(cluster)]), names(expected_clusters)) ||
      !identical(as.integer(table(cluster)[names(expected_clusters)]),
                 unname(expected_clusters)) || sum(is.na(cluster)) != 6333L) {
    stop("Transferred-cluster counts differ from inspected checkpoint.")
  }

  # Compare saved annotation labels by cell ID, never by assumed row order.
  annotation <- readRDS(annotation_file)
  if (!is.list(annotation) || is.null(annotation$projection.clusters)) {
    stop("Annotation file lacks projection.clusters; no exports written.")
  }
  reference <- annotation$projection.clusters
  ref_ids <- names(reference)
  if (is.null(ref_ids) || anyNA(ref_ids) || any(ref_ids == "") ||
      anyDuplicated(ref_ids) || !all(ids %in% ref_ids)) {
    stop("Annotation cell IDs cannot be joined unambiguously; no exports written.")
  }
  expected <- as.character(reference[match(ids, ref_ids)])
  same <- (is.na(cluster) & is.na(expected)) |
    (!is.na(cluster) & !is.na(expected) & cluster == expected)
  mismatch <- sum(!same)
  if (mismatch != 0L) {
    stop(mismatch, " annotation mismatches; no exports written. Example cell IDs: ",
         paste(head(ids[!same], 8L), collapse = ", "), call. = FALSE)
  }
  message("PASS: 39,787 cells; sample and cluster counts match the inspected checkpoint.")
  message("PASS: saved annotation agrees for all cells, including 6,333 unassigned cells.")

  # Only copy the small embeddings and metadata, not assay matrices.
  tables <- vector("list", length(output_files))
  report <- vector("list", length(output_files))
  for (i in seq_along(output_files)) {
    reduction <- names(output_files)[i]
    if (!reduction %in% names(obj@reductions)) stop("Missing reduction: ", reduction)
    xy <- SeuratObject::Embeddings(obj, reduction = reduction)
    xy_ids <- rownames(xy)
    if (!is.matrix(xy) || ncol(xy) != 2L || nrow(xy) != length(ids) ||
        is.null(xy_ids) || anyNA(xy_ids) || anyDuplicated(xy_ids) ||
        !setequal(xy_ids, ids) || any(!is.finite(xy))) {
      stop("Invalid cell IDs or coordinates in ", reduction, "; no exports written.")
    }
    mi <- match(xy_ids, ids)
    # Preserve checkpoint embedding order (also useful for reproducing point order).
    tables[[i]] <- data.frame(
      cell_id = xy_ids, UMAP_1 = unname(xy[, 1]), UMAP_2 = unname(xy[, 2]),
      sample = sample[mi], transferred_cluster = cluster[mi],
      stringsAsFactors = FALSE, row.names = NULL
    )
    report[[i]] <- data.frame(
      panel = panel_names[i], reduction = reduction, n_cells = nrow(xy),
      n_annotated = sum(!is.na(cluster)), n_unassigned = sum(is.na(cluster)),
      n_samples = length(unique(sample)), annotation_mismatches = mismatch,
      all_coordinates_finite = TRUE, coordinate_roundtrip_max_abs_difference = NA_real_,
      export_matches_checkpoint = FALSE, object_source = object_source,
      status = "PENDING", stringsAsFactors = FALSE
    )
  }

  for (d in c(source_dir, validation_dir)) {
    if (!dir.exists(d) && !dir.create(d, recursive = TRUE)) stop("Cannot create: ", d)
  }
  temporary <- vapply(seq_along(final_files), function(i)
    tempfile(pattern = ".embedding-export-", tmpdir = source_dir, fileext = ".tsv.gz"),
    character(1))
  report_tmp <- tempfile(pattern = ".embedding-report-", tmpdir = validation_dir)
  on.exit(unlink(c(temporary, report_tmp)), add = TRUE)
  write_gz <- function(tab, file) {
    con <- gzfile(file, open = "wt")
    on.exit(close(con))
    write.table(tab, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  read_gz <- function(file) {
    con <- gzfile(file, open = "rt")
    on.exit(close(con))
    read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
               colClasses = "character", na.strings = "NA")
  }
  for (i in seq_along(tables)) {
    tab <- tables[[i]]
    text <- tab
    # Explicit 17-significant-digit formatting avoids avoidable TSV precision loss.
    text$UMAP_1 <- trimws(formatC(tab$UMAP_1, format = "g", digits = 17, decimal.mark = "."))
    text$UMAP_2 <- trimws(formatC(tab$UMAP_2, format = "g", digits = 17, decimal.mark = "."))
    write_gz(text, temporary[i])
    back <- read_gz(temporary[i])
    same_meta <- identical(back$cell_id, tab$cell_id) &&
      identical(back$sample, tab$sample) &&
      identical(back$transferred_cluster, tab$transferred_cluster)
    original_xy <- c(tab$UMAP_1, tab$UMAP_2)
    read_xy <- c(as.numeric(back$UMAP_1), as.numeric(back$UMAP_2))
    same_xy <- length(read_xy) == length(original_xy) && all(is.finite(read_xy)) &&
      all(abs(read_xy - original_xy) <= 1e-14 * pmax(1, abs(original_xy)))
    if (!same_meta || !same_xy) stop("TSV read-back check failed: ", output_files[i])
    report[[i]]$coordinate_roundtrip_max_abs_difference <- max(abs(read_xy - original_xy))
    report[[i]]$export_matches_checkpoint <- TRUE
    report[[i]]$status <- "PASS"
  }
  summary <- do.call(rbind, report)
  write.table(summary, report_tmp, sep = "\t", quote = FALSE, row.names = FALSE)
  # Publish only after all three temporary exports have passed read-back checks.
  targets <- c(final_files, report_file)
  if (any(file.exists(targets))) stop("An output appeared during export; refusing to replace it.")
  for (i in seq_along(temporary)) {
    if (!file.rename(temporary[i], final_files[i])) stop("Cannot finalize: ", final_files[i])
  }
  if (!file.rename(report_tmp, report_file)) stop("Cannot finalize summary: ", report_file)
  print(summary, row.names = FALSE)
  message("Three Source Data tables exported and checked. No figures or original data were modified.")
  invisible(summary)
}
