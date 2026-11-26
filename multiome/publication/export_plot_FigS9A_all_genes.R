# S9A alternative: gene-detection fractions WITHOUT the historical DE gene list.
# Run from CH_publication_release. Sourcing this file only defines the function.
#
#   source("multiome/publication/export_plot_FigS9A_all_genes.R")
#   S9A_result <- reproduce_FigS9A_all_genes()
#
# To reuse an explicitly loaded integration object instead of reloading:
#   S9A_result <- reproduce_FigS9A_all_genes(
#       integration = s9ad_checkpoint$rna.combined)
#
# Historical basis: multiome/scripts/1_rna_integration.R, lines 219-229.
# CHANGE: all rows of the joined RNA counts are exported, rather than gg
# (the gene intersection derived from the 60 DE comparisons).
# The source table retains zero-detection genes. Only genes with positive
# detection fractions in BOTH groups can be shown on the logarithmic axes.
# Cells with missing projected annotations are NOT excluded.
# No PrepSCTFindMarkers(), FindMarkers(), normalization, or UMAP computation.
# No original checkpoint, historical figure, or existing approved panel is edited.
#
# Layer access: https://satijalab.github.io/seurat-object/reference/Layers.html

reproduce_FigS9A_all_genes <- function(integration = NULL, root = getwd(),
                                      output_root = NULL,
                                      chunk_cells = 1000L,
                                      overwrite = FALSE) {
  needed <- c("Seurat", "SeuratObject", "Matrix", "ggplot2")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Packages unavailable: ", paste(missing, collapse = ", "),
                            ". No packages were installed.", call. = FALSE)
  if (length(chunk_cells) != 1L || !is.numeric(chunk_cells) ||
      !is.finite(chunk_cells) || chunk_cells < 1 || chunk_cells > 100000L ||
      chunk_cells != floor(chunk_cells)) {
    stop("chunk_cells must be an integer between 1 and 100000.")
  }
  if (length(overwrite) != 1L || is.na(overwrite) || !is.logical(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }
  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome"))) {
    stop("Run from CH_publication_release, not from multiome/.")
  }
  input <- file.path(root, "multiome/data/processed/rna_integration.RData")
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  src_dir <- file.path(output_root, "source_data")
  fig_dir <- file.path(output_root, "figures")
  protected <- file.path(root, c("multiome/data", "multiome/figures", "multiome/scripts"))
  for (p in c(src_dir, fig_dir)) {
    for (b in protected) {
      if (identical(p, b) || startsWith(p, paste0(b, "/"))) {
        stop("Refusing to write into historical directories: ", p)
      }
    }
  }
  src_file <- file.path(src_dir, "FigS9A_gene_detection_all_genes.tsv.gz")
  fig_file <- file.path(fig_dir, "FigS9A_gene_detection_all_genes_REVIEW.pdf")
  targets <- c(src_file, fig_file)
  if (!overwrite && any(file.exists(targets))) {
    stop("All-gene outputs already exist; they were retained:\n",
         paste(targets[file.exists(targets)], collapse = "\n"),
         "\nUse a separate output_root, or overwrite=TRUE only for this derivative.",
         call. = FALSE)
  }

  if (is.null(integration)) {
    if (!file.exists(input)) stop("Missing checkpoint: ", input)
    message("Loading rna_integration.RData (large checkpoint; no DE will run).")
    e <- new.env(parent = emptyenv())
    loaded <- load(input, envir = e)
    if (!"rna.combined" %in% loaded) stop("Checkpoint lacks rna.combined.")
    integration <- get("rna.combined", envir = e, inherits = FALSE)
    # Release the other objects loaded from the RData, not any caller objects.
    rm(list = ls(e, all.names = TRUE), envir = e)
    rm(e)
    invisible(gc(verbose = FALSE))
    object_source <- "rna_integration.RData"
  } else {
    message("Using explicitly supplied integration object; checkpoint not reloaded.")
    object_source <- "explicit_integration_argument"
  }
  if (!inherits(integration, "Seurat")) stop("integration must be a Seurat object.")
  if (!"RNA" %in% SeuratObject::Assays(integration)) stop("RNA assay is missing.")
  layers <- SeuratObject::Layers(integration[["RNA"]], search = NA)
  count_layers <- grep("^counts($|\\.)", layers, value = TRUE)
  if (!identical(count_layers, "counts")) {
    stop("Expected the saved, joined RNA counts layer. Found: ",
         paste(count_layers, collapse = ", "), ". No layers were joined or changed.")
  }
  counts <- SeuratObject::LayerData(integration, assay = "RNA", layer = "counts")
  meta <- integration[[]]
  # Only counts and small metadata are needed from here onward.
  rm(integration)
  invisible(gc(verbose = FALSE))
  if (!inherits(counts, "sparseMatrix")) {
    stop("Expected sparse RNA counts; refusing to create a large dense matrix.")
  }
  genes <- rownames(counts)
  cells <- colnames(counts)
  if (is.null(genes) || is.null(cells) || !length(genes) || !length(cells) ||
      anyNA(genes) || anyNA(cells) || anyDuplicated(genes) || anyDuplicated(cells) ||
      any(!nzchar(genes)) || any(!nzchar(cells)) ||
      any(grepl("[\t\r\n]", genes))) stop("Invalid gene or cell identifiers.")
  if (!"dataset" %in% names(meta) || anyDuplicated(rownames(meta))) {
    stop("Missing dataset annotations or duplicated metadata cell IDs.")
  }
  ix <- match(cells, rownames(meta))
  if (anyNA(ix) || !setequal(cells, rownames(meta))) {
    stop("RNA counts and metadata cell IDs do not agree.")
  }
  dataset <- as.character(meta$dataset[ix])
  expected <- c(exp1 = 51282L, exp2 = 48479L, C1 = 11004L,
                C2 = 7354L, D = 9413L, T = 16267L)
  if (anyNA(dataset) || !setequal(unique(dataset), names(expected))) {
    stop("Unexpected dataset labels. Use rna_integration.RData, not multiome_seurat.rds.")
  }
  observed <- as.integer(table(factor(dataset, levels = names(expected))))
  if (!identical(observed, as.integer(expected))) {
    stop("Dataset counts do not match the inspected 143,799-cell integration checkpoint.")
  }
  old_cells <- which(dataset %in% c("exp1", "exp2"))
  new_cells <- which(dataset %in% c("C1", "C2", "D", "T"))
  n_old <- length(old_cells)
  n_new <- length(new_cells)
  message("Cells: ", n_old, " scRNA-seq; ", n_new,
          " Multiome RNA. Counting detection for all ", length(genes), " RNA genes.")

  count_detected <- function(indices) {
    total <- numeric(nrow(counts))
    for (start in seq.int(1L, length(indices), by = as.integer(chunk_cells))) {
      take <- indices[seq.int(start, min(length(indices), start + chunk_cells - 1L))]
      block <- counts[, take, drop = FALSE]
      if (anyNA(block) || !is.finite(max(block)) || min(block) < 0) {
        stop("RNA counts contain missing, nonfinite, or negative values.")
      }
      total <- total + as.numeric(Matrix::rowSums(block > 0))
    }
    total
  }
  # No dense coercion and no full-size binary copy of the gene-by-cell matrix.
  detected_old <- count_detected(old_cells)
  detected_new <- count_detected(new_cells)
  rm(counts, meta)
  invisible(gc(verbose = FALSE))
  stopifnot(all(detected_old >= 0 & detected_old <= n_old),
            all(detected_new >= 0 & detected_new <= n_new))
  tab <- data.frame(
    gene = genes,
    scRNA_detected_cells = as.integer(detected_old),
    scRNA_total_cells = n_old,
    scRNA_detected_fraction = detected_old / n_old,
    multiome_detected_cells = as.integer(detected_new),
    multiome_total_cells = n_new,
    multiome_detected_fraction = detected_new / n_new,
    plotted_on_log_axes = detected_old > 0 & detected_new > 0,
    stringsAsFactors = FALSE
  )
  if (!any(tab$plotted_on_log_axes)) stop("No genes have positive detection in both groups.")

  for (d in c(src_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  tmp_src <- tempfile(".S9A_source_", tmpdir = src_dir, fileext = ".tsv.gz")
  tmp_pdf <- tempfile(".S9A_plot_", tmpdir = fig_dir, fileext = ".pdf")
  on.exit(unlink(c(tmp_src, tmp_pdf)), add = TRUE)
  con <- gzfile(tmp_src, "wt")
  tryCatch(utils::write.table(tab, con, sep = "\t", row.names = FALSE,
                             quote = FALSE, na = "NA"), finally = close(con))
  con <- gzfile(tmp_src, "rt")
  saved <- tryCatch(utils::read.delim(con, check.names = FALSE,
                                     stringsAsFactors = FALSE,
                                     colClasses = c("character", rep("numeric", 6), "logical")),
                    finally = close(con))
  # Test the Source Data round-trip, not equivalence to the historical DE selection.
  if (!isTRUE(all.equal(tab, saved, tolerance = 1e-14, check.attributes = FALSE))) {
    stop("Source Data read-back does not match the calculated detection values.")
  }
  plot_data <- saved[saved$plotted_on_log_axes, , drop = FALSE]
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = scRNA_detected_fraction, y = multiome_detected_fraction)) +
    ggplot2::geom_point(size = 1, alpha = 0.05) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "scRNA-seq", y = "multiome", title = "mRNA coverage",
                  subtitle = "All genes; no DE-based selection") +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   plot.subtitle = ggplot2::element_text(size = 9))
  ggplot2::ggsave(tmp_pdf, plot = p, device = grDevices::pdf,
                  width = 4, height = 4, units = "in", useDingbats = FALSE)
  sizes <- file.info(c(tmp_src, tmp_pdf))$size
  if (anyNA(sizes) || any(sizes <= 0)) stop("An output was empty or missing.")
  if (!overwrite && any(file.exists(targets))) stop("Outputs appeared during calculation; retained them.")
  if (!file.rename(tmp_src, src_file)) stop("Could not finalize Source Data: ", src_file)
  if (!file.rename(tmp_pdf, fig_file)) stop("Could not finalize figure: ", fig_file)

  summary <- data.frame(
    panel = "FigS9A_all_genes_alternative", scRNA_cells = n_old,
    multiome_cells = n_new, genes_exported = nrow(saved),
    genes_plotted = nrow(plot_data),
    genes_zero_in_one_or_both = sum(!saved$plotted_on_log_axes),
    DE_gene_selection = "NOT_APPLIED", source_data_roundtrip = "PASS",
    historical_equivalence = "NOT_CLAIMED", object_source = object_source,
    stringsAsFactors = FALSE
  )
  print(summary, row.names = FALSE)
  message("Source Data: ", src_file, "\nFigure: ", fig_file)
  message("All-gene derivative generated. Original files were not modified.")
  invisible(list(summary = summary, source_data = src_file, figure = fig_file))
}
