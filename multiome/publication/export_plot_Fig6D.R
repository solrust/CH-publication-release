# Fig 6D: clone overlays on the validated WNN embedding.
# Run with CH_publication_release as the working directory.
#
# source("multiome/publication/export_plot_Fig6D.R")
# Fig6D_result <- reproduce_Fig6D()
#
# Reads only the existing WNN Source Data and final cell metadata.
# No Seurat object, clone calling, or embedding reconstruction is required.
# Outputs are placed under multiome/publication/ and never overwrite existing files.

reproduce_Fig6D <- function(root = getwd(), output_root = NULL,
                                wnn_source_file = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Required packages: ggplot2 and gridExtra. No outputs were written.",
         call. = FALSE)
  }

  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome"))) {
    stop("Run from CH_publication_release (not from multiome/).", call. = FALSE)
  }
  # For an independent publication run, use the WNN table generated in that
  # run, not an approved first-pass Source Data file.
  if (is.null(wnn_source_file)) {
    wnn_source_file <- file.path(root, "multiome/publication/source_data",
                                  "Fig6B_FigS9B_WNN_UMAP.tsv.gz")
  }
  wnn_file <- wnn_source_file
  metadata_file <- file.path(root, "multiome/data/processed/metainfo.txt")
  if (any(!file.exists(c(wnn_file, metadata_file)))) {
    stop("Missing input(s):\n",
         paste(c(wnn_file, metadata_file)[!file.exists(c(wnn_file, metadata_file))],
               collapse = "\n"), call. = FALSE)
  }

  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  protected <- normalizePath(file.path(root, "multiome/figures"), mustWork = TRUE)
  if (identical(output_root, protected) ||
      startsWith(output_root, paste0(protected, "/"))) {
    stop("Refusing to write to the historical figures directory.", call. = FALSE)
  }
  out_data <- file.path(output_root, "source_data", "Fig6D_clone_umap_cells.tsv.gz")
  out_figure <- file.path(output_root, "figures", "Fig6D_reproduced_REVIEW.pdf")
  if (any(file.exists(c(out_data, out_figure)))) {
    stop("Existing output(s) retained; nothing overwritten:\n",
         paste(c(out_data, out_figure)[file.exists(c(out_data, out_figure))],
               collapse = "\n"), call. = FALSE)
  }

  # Exact IDs identified from final Multiome metadata and historical Fig 6D.
  examples <- data.frame(
    panel = c("myeloid biased", "lymphoid biased", "HSC biased"),
    clone_id = c(
      "C1_AAAATGCCAACAGCACACAGTCGAATTCGTTGAGAGGTCA",
      "D_CTGGTGATCTCAAATTACGGGTGAGTTTGTCAACAGTTGA",
      "T_TGGATGCGCCCATGATACAAGGGACTAGGTGTGGAGAGCG"
    ),
    sample = c("C1", "D", "T"),
    size = c(110L, 818L, 117L),
    source_flag = c("is_clone_110_C1", "is_clone_818_D", "is_clone_117_T"),
    display_name = c("clone 1   AAA...", "clone 7   CTG...", "clone 11   TG..."),
    color = c("#1737FF", "#B52AFF", "#FF9F00"),
    stringsAsFactors = FALSE
  )

  con <- gzfile(wnn_file, "rt")
  wnn <- tryCatch(
    utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE),
    finally = close(con)
  )
  meta <- utils::read.delim(metadata_file, row.names = 1,
                            check.names = FALSE, stringsAsFactors = FALSE)
  required_wnn <- c("cell_id", "UMAP_1", "UMAP_2", "sample", "transferred_cluster")
  required_meta <- c("sample", "cloneid", "projection", "wnnUMAP_1", "wnnUMAP_2")
  if (!all(required_wnn %in% names(wnn)) ||
      !all(required_meta %in% names(meta))) {
    stop("Missing expected columns in WNN export or metainfo.txt.", call. = FALSE)
  }
  if (nrow(wnn) != 39787L || nrow(meta) != 39787L ||
      anyNA(wnn$cell_id) || anyDuplicated(wnn$cell_id) ||
      anyDuplicated(rownames(meta)) ||
      !setequal(wnn$cell_id, rownames(meta))) {
    stop("Unexpected cell counts or ambiguous cell identifiers.", call. = FALSE)
  }

  # Join by cell identifier, never by assumed row order.
  meta <- meta[match(wnn$cell_id, rownames(meta)), , drop = FALSE]
  if (!identical(as.character(wnn$sample), as.character(meta$sample))) {
    stop("Sample labels differ between the validated WNN table and metainfo.txt.",
         call. = FALSE)
  }
  annotations_match <- (is.na(wnn$transferred_cluster) & is.na(meta$projection)) |
    (!is.na(wnn$transferred_cluster) & !is.na(meta$projection) &
       as.character(wnn$transferred_cluster) == as.character(meta$projection))
  if (any(!annotations_match)) {
    stop("Transferred cluster labels differ between WNN and final metadata.",
         call. = FALSE)
  }
  xy <- cbind(as.numeric(wnn$UMAP_1), as.numeric(wnn$UMAP_2))
  historical_xy <- cbind(as.numeric(meta$wnnUMAP_1), as.numeric(meta$wnnUMAP_2))
  if (any(!is.finite(xy)) || any(!is.finite(historical_xy))) {
    stop("Nonfinite WNN coordinates.", call. = FALSE)
  }
  max_coord_diff <- max(abs(xy - historical_xy))
  if (max_coord_diff > 1e-5) {
    stop("The saved WNN Source Data differs from final metadata (max difference ",
         signif(max_coord_diff, 4), ").", call. = FALSE)
  }

  dat <- data.frame(
    cell_id = wnn$cell_id,
    sample = as.character(wnn$sample),
    cloneid = as.character(meta$cloneid),
    projection = as.character(wnn$transferred_cluster),
    wnnUMAP_1 = xy[, 1],
    wnnUMAP_2 = xy[, 2],
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(examples))) {
    ex <- examples[i, ]
    is_highlighted <- !is.na(dat$cloneid) & dat$cloneid == ex$clone_id
    if (sum(is_highlighted) != ex$size ||
        !all(dat$sample[is_highlighted] == ex$sample)) {
      stop("Clone count/sample check failed: ", ex$clone_id,
           " (observed ", sum(is_highlighted), ", expected ", ex$size, ").",
           call. = FALSE)
    }
    dat[[ex$source_flag]] <- as.integer(is_highlighted)
  }

  # Draw historical-style grey background with colored open clone circles.
  # The WNN coordinate system is shared with the validated Fig 6B table.
  x_limits <- range(dat$wnnUMAP_1) + c(-0.6, 0.6)
  y_limits <- range(dat$wnnUMAP_2) + c(-0.6, 0.6)
  plots <- lapply(seq_len(nrow(examples)), function(i) {
    ex <- examples[i, ]
    highlighted <- dat[dat[[ex$source_flag]] == 1L, , drop = FALSE]
    ggplot2::ggplot(dat, ggplot2::aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
      ggplot2::geom_point(color = "#BDBDBD", size = 0.18, shape = 16) +
      ggplot2::geom_point(data = highlighted, shape = 21, fill = "white",
                          color = ex$color, size = 0.95, stroke = 0.65) +
      ggplot2::coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE) +
      ggplot2::labs(title = ex$panel,
                    subtitle = paste0(ex$display_name, "\nsize = ", ex$size)) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, face = "plain", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, hjust = 0),
        plot.margin = ggplot2::margin(4, 4, 4, 4)
      )
  })

  for (d in c(dirname(out_data), dirname(out_figure))) {
    if (!dir.exists(d) && !dir.create(d, recursive = TRUE)) {
      stop("Could not create output folder: ", d)
    }
  }
  temp_data <- tempfile(".Fig6D-data-", tmpdir = dirname(out_data), fileext = ".tsv.gz")
  temp_figure <- tempfile(".Fig6D-figure-", tmpdir = dirname(out_figure), fileext = ".pdf")
  on.exit(unlink(c(temp_data, temp_figure)), add = TRUE)

  # Write with full coordinate precision, then verify by reading back.
  to_write <- dat
  for (field in c("wnnUMAP_1", "wnnUMAP_2")) {
    to_write[[field]] <- formatC(dat[[field]], digits = 17L, format = "g")
  }
  write_table <- function(x, file) {
    con <- gzfile(file, "wt")
    on.exit(close(con))
    utils::write.table(x, con, sep = "\t", row.names = FALSE,
                       quote = FALSE, na = "NA")
  }
  write_table(to_write, temp_data)
  back <- utils::read.delim(gzfile(temp_data, "rt"), check.names = FALSE,
                            stringsAsFactors = FALSE)
  same_text <- identical(back$cell_id, dat$cell_id) &&
    identical(back$sample, dat$sample) &&
    identical(back$cloneid, dat$cloneid) &&
    identical(as.character(back$projection), as.character(dat$projection))
  same_flags <- all(vapply(examples$source_flag, function(f)
    identical(as.integer(back[[f]]), dat[[f]]), logical(1)))
  max_readback_diff <- max(abs(as.matrix(back[, c("wnnUMAP_1", "wnnUMAP_2")]) - xy))
  if (!same_text || !same_flags || !is.finite(max_readback_diff) ||
      max_readback_diff > 1e-12) {
    stop("Fig 6D Source Data read-back verification failed.", call. = FALSE)
  }

  write_pdf <- function(file, grobs) {
    grDevices::pdf(file, width = 11.4, height = 3.9, useDingbats = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    gridExtra::grid.arrange(grobs = grobs, ncol = 3)
  }
  write_pdf(temp_figure, plots)
  if (!is.finite(file.info(temp_figure)$size) || file.info(temp_figure)$size <= 0) {
    stop("Fig 6D PDF is empty; no final files written.", call. = FALSE)
  }
  if (any(file.exists(c(out_data, out_figure)))) {
    stop("An output was created by another process; refusing to overwrite it.",
         call. = FALSE)
  }
  if (!file.rename(temp_data, out_data)) stop("Unable to finalize ", out_data)
  if (!file.rename(temp_figure, out_figure)) stop("Unable to finalize ", out_figure)

  summary <- data.frame(
    panel = examples$panel,
    sample = examples$sample,
    clone_size = examples$size,
    highlighted_cells = vapply(examples$source_flag, function(f)
      sum(dat[[f]]), integer(1)),
    n_background_cells = nrow(dat),
    max_WNN_difference_vs_metainfo = max_coord_diff,
    source_data_roundtrip_max_abs_difference = max_readback_diff,
    stringsAsFactors = FALSE
  )
  print(summary, row.names = FALSE)
  message("Fig 6D Source Data: ", out_data)
  message("Fig 6D figure for comparison: ", out_figure)
  message("First-pass plot created; visual agreement with publication is still REVIEW.")
  invisible(summary)
}
