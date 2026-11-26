# Fig 6B, Fig S9B, Fig S9C: first-pass graphics from exported Source Data.
# Working directory: CH_publication_release.
# source("multiome/publication/plot_multiome_batch1_embeddings.R")
# plot_summary <- plot_multiome_batch1_embeddings()
#
# No Seurat checkpoint is loaded; no analysis, coordinates, or annotations change.
# The function writes ONLY to the publication output directory. Original figures,
# source tables, historical scripts, and the Seurat object are never modified.
# Existing output files are protected unless overwrite=TRUE is explicitly set.
#
# Provenance:
# - Historical recipe: multiome/scripts/2_multiome_processing1.R, lines 11-115.
# - Cluster palette: exact RGB label-box fills read from the supplied original
#   multiome_integration.png and cross-checked in the RNA/ATAC PNGs. These are
#   the displayed colors of old.scrna$clust.final.cols, not a newly chosen palette.
# - Sample palette: exact RGB values of the C1/C2/D/T legend markers in that PNG.
# - Cluster display names: supplied Figure_6_V3 artwork; IDs in tables unchanged.
# - Point sizing/stroke and shuffle seed follow Seurat 5.3.1 visualization.R.
#   https://raw.githubusercontent.com/satijalab/seurat/v5.3.1/R/visualization.R
#
# PANEL MAPPING NOTE:
# The supplied S9B artwork visually matches the sample half of
# multiome_integration_rna.png, NOT multiome_integration.png (WNN).
# Consequently the S9B REVIEW plot below uses the exported RNA embedding.
# The three-page review PDF includes both RNA and WNN sample views for approval.
# The existing WNN TSV's FigS9B filename tag is not treated as proof of mapping.
#
# Reproducibility: fixed seeds for point order and label layout; RNG state restored.
# Label placement/font rasterization need not match historical pixels exactly.
# First-pass visual approval is still required; the report does NOT claim it.
# For the eventual clean rerun, pass a separate output_root. No rerun is run here.

plot_multiome_batch1_embeddings <- function(root = getwd(), source_dir = NULL,
                                            output_root = NULL, overwrite = FALSE) {
  pkgs <- c("ggplot2", "ggrepel", "cowplot", "gridExtra", "png")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Packages cannot load: ", paste(missing, collapse = ", "),
                            ". No files were changed.", call. = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome"))) {
    stop("Run from CH_publication_release, not multiome/.", call. = FALSE)
  }
  if (is.null(source_dir)) source_dir <- file.path(root, "multiome/publication/source_data")
  source_dir <- normalizePath(source_dir, mustWork = TRUE)
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  protected <- normalizePath(file.path(root, "multiome/figures"), mustWork = TRUE)
  if (identical(output_root, dirname(protected)) || identical(output_root, protected) ||
      startsWith(output_root, paste0(protected, "/"))) {
    stop("Refusing to write into the historical figure directory.", call. = FALSE)
  }
  figure_dir <- file.path(output_root, "figures")
  validation_dir <- file.path(output_root, "validation")
  tables <- c(wnn = "Fig6B_FigS9B_WNN_UMAP.tsv.gz",
              rna = "FigS9C_RNA_UMAP.tsv.gz", atac = "FigS9C_ATAC_UMAP.tsv.gz")
  originals <- file.path(protected, c("multiome_integration.png",
                                      "multiome_integration_rna.png",
                                      "multiome_integration_atac.png"))
  files_needed <- c(file.path(source_dir, tables), originals)
  if (any(!file.exists(files_needed))) stop("Missing files:\n",
    paste(files_needed[!file.exists(files_needed)], collapse = "\n"), call. = FALSE)
  pdf_names <- c("Fig6B_reproduced_REVIEW.pdf", "FigS9B_reproduced_REVIEW.pdf",
                 "FigS9C_reproduced_REVIEW.pdf", "Batch1_embeddings_visual_review.pdf")
  targets <- c(file.path(figure_dir, pdf_names),
               file.path(validation_dir, "Batch1_embeddings_plot_summary.tsv"))
  if (!overwrite && any(file.exists(targets))) stop(
    "Existing outputs retained:\n", paste(targets[file.exists(targets)], collapse = "\n"),
    "\nUse a separate output_root or explicitly set overwrite=TRUE.", call. = FALSE)

  # Set and restore all relevant RNG settings, even on error.
  previous_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")

  read_export <- function(file) {
    con <- gzfile(file, "rt"); on.exit(close(con))
    x <- utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
                          colClasses = "character", na.strings = "NA")
    required <- c("cell_id", "UMAP_1", "UMAP_2", "sample", "transferred_cluster")
    if (!identical(names(x), required)) stop("Unexpected columns in ", basename(file))
    x$UMAP_1 <- as.numeric(x$UMAP_1); x$UMAP_2 <- as.numeric(x$UMAP_2)
    if (nrow(x) != 39787L || anyNA(x$cell_id) || any(x$cell_id == "") ||
        anyDuplicated(x$cell_id) || any(!is.finite(x$UMAP_1)) ||
        any(!is.finite(x$UMAP_2)) || anyNA(x$sample)) stop("Invalid export: ", file)
    expected_sample <- c(C1 = 9690L, C2 = 6840L, D = 8799L, T = 14458L)
    expected_cluster <- c("1" = 4020L, "2" = 2115L, "3" = 1213L, "4" = 8760L,
       "5" = 4728L, "6" = 7382L, "7" = 1820L, "8" = 1848L, "9" = 348L, "10" = 1220L)
    if (!setequal(unique(x$sample), names(expected_sample)) ||
        !identical(as.integer(table(x$sample)[names(expected_sample)]), unname(expected_sample)) ||
        !setequal(unique(stats::na.omit(x$transferred_cluster)), names(expected_cluster)) ||
        !identical(as.integer(table(x$transferred_cluster)[names(expected_cluster)]),
                   unname(expected_cluster)) || sum(is.na(x$transferred_cluster)) != 6333L) {
      stop("Sample/cluster counts differ from the validated exports: ", file)
    }
    x
  }
  dat <- lapply(file.path(source_dir, tables), read_export); names(dat) <- names(tables)
  ref <- dat$wnn
  for (x in dat) {
    ix <- match(ref$cell_id, x$cell_id)
    if (anyNA(ix) || !identical(ref$sample, x$sample[ix]) ||
        !identical(ref$transferred_cluster, x$transferred_cluster[ix]))
      stop("Embeddings disagree on cell IDs, samples or annotations.")
  }

  cluster_colors <- c("1" = "#A6CEE3", "2" = "#1F78B4", "3" = "#B2DF8A",
    "4" = "#33A02C", "5" = "#FB9A99", "6" = "#E31A1C", "7" = "#FDBF6F",
    "8" = "#FF7F00", "9" = "#CAB2D6", "10" = "#6A3D9A", Unassigned = "#CCCCCC")
  sample_colors <- c(C1 = "#F8766D", C2 = "#7CAE00", D = "#00BFC4", T = "#C77CFF")
  cluster_labels <- c("1 - HSC", "2 - LMPP-1", "3 - MLP", "4 - Pro-B", "5 - pDC",
     "6 - LMPP-2", "7 - early GMP", "8 - late GMP", "9 - EMPP/MkP", "10 - EryP")

  make_plot <- function(x, grouping, title, show_cluster_legend = FALSE) {
    ids <- unique(x$transferred_cluster[!is.na(x$transferred_cluster)])
    set.seed(1L)
    index <- sample.int(nrow(x))
    d <- x[index, , drop = FALSE]
    groups <- if (grouping == "cluster") d$transferred_cluster else d$sample
    groups[is.na(groups)] <- "Unassigned" # display-only; table remains unchanged
    palette <- if (grouping == "cluster") cluster_colors else sample_colors
    d$plot_group <- factor(groups, levels = names(palette))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = plot_group)) +
      ggplot2::geom_point(size = min(1583 / nrow(d), 1), stroke = 0.600075815011372,
                         shape = 16, alpha = 1) +
      ggplot2::scale_color_manual(values = palette, drop = FALSE,
         breaks = if (grouping == "cluster") as.character(1:10) else names(sample_colors),
         labels = if (grouping == "cluster") cluster_labels else names(sample_colors)) +
      cowplot::theme_cowplot(font_size = 12, font_family = "sans") +
      ggplot2::labs(x = "UMAP 1", y = "UMAP 2", title = title,
                   color = if (grouping == "cluster") "Clusters from RNA" else NULL) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14,
                                                      face = "plain")) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
    if (grouping == "cluster") {
      centers <- do.call(rbind, lapply(ids, function(id) {
        ind <- !is.na(x$transferred_cluster) & x$transferred_cluster == id
        data.frame(UMAP_1 = mean(x$UMAP_1[ind]), UMAP_2 = mean(x$UMAP_2[ind]), label = id)
      }))
      p <- p + ggrepel::geom_label_repel(data = centers,
        mapping = ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
        inherit.aes = FALSE, fill = unname(cluster_colors[centers$label]),
        color = "black", size = 3, family = "sans", seed = 1L,
        max.time = Inf, max.iter = 10000, max.overlaps = Inf, show.legend = FALSE) +
        ggplot2::theme(legend.position = if (show_cluster_legend) "right" else "none")
    } else {
      # Use existing ggplot2 convention, including versions before 3.5.0.
      p <- p + ggplot2::theme(legend.position = c(0.12, 0.82),
          legend.background = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank())
    }
    # Validate the actual point layer, not just file sizes.
    built <- ggplot2::ggplot_build(p)$data[[1]]
    expected_color <- unname(palette[as.character(d$plot_group)])
    if (nrow(built) != nrow(d) || anyNA(built$colour) ||
        !identical(as.numeric(built$x), d$UMAP_1) ||
        !identical(as.numeric(built$y), d$UMAP_2) ||
        !identical(toupper(as.character(built$colour)), toupper(expected_color))) {
      stop("Plot-layer check failed: ", title)
    }
    list(plot = p, n_points = nrow(built), max_coord_difference = 0,
         n_unassigned = sum(is.na(d$transferred_cluster)))
  }
  plots <- list(
    wnn_cluster = make_plot(dat$wnn, "cluster", "RNA-ATAC coembedding", TRUE),
    wnn_sample = make_plot(dat$wnn, "sample", "Samples: WNN embedding"),
    rna_cluster = make_plot(dat$rna, "cluster", "RNA-only embedding"),
    rna_sample = make_plot(dat$rna, "sample", "Sample distribution: RNA embedding"),
    atac_cluster = make_plot(dat$atac, "cluster", "ATAC-only embedding"),
    atac_sample = make_plot(dat$atac, "sample", "Samples: ATAC embedding")
  )
  message("PASS: all six plot layers retain 39,787 points, unchanged coordinates and mapped colors.")
  message("S9B mapping for review: RNA-only sample view; WNN alternative included in review PDF.")

  for (d in c(figure_dir, validation_dir)) {
    if (!dir.exists(d) && !dir.create(d, recursive = TRUE)) stop("Cannot create ", d)
  }
  temps <- c(vapply(seq_along(pdf_names), function(i)
    tempfile(".embedding-plot-", tmpdir = figure_dir, fileext = ".pdf"), character(1)),
    tempfile(".embedding-summary-", tmpdir = validation_dir, fileext = ".tsv"))
  on.exit(unlink(temps), add = TRUE)
  ggplot2::ggsave(temps[1], plot = plots$wnn_cluster$plot, device = "pdf",
                 width = 6, height = 4.2, units = "in", useDingbats = FALSE)
  ggplot2::ggsave(temps[2], plot = plots$rna_sample$plot, device = "pdf",
                 width = 4.2, height = 4.2, units = "in", useDingbats = FALSE)
  s9c <- gridExtra::arrangeGrob(plots$atac_cluster$plot, plots$rna_cluster$plot, nrow = 1)
  ggplot2::ggsave(temps[3], plot = s9c, device = "pdf",
                 width = 8.4, height = 4.2, units = "in", useDingbats = FALSE)

  write_review <- function(file) {
    grDevices::pdf(file, width = 9, height = 9.2, onefile = TRUE, useDingbats = FALSE)
    device_id <- grDevices::dev.cur()
    on.exit({
      if (device_id %in% grDevices::dev.list()) grDevices::dev.off(device_id)
    }, add = TRUE)
    keys <- list(c("wnn_cluster", "wnn_sample"), c("rna_cluster", "rna_sample"),
                 c("atac_cluster", "atac_sample"))
    headings <- c("WNN: Fig 6B and alternative sample view",
                  "RNA: Fig S9C (right) and candidate Fig S9B",
                  "ATAC: Fig S9C (left) and sample view")
    for (i in seq_along(originals)) {
      image <- png::readPNG(originals[i])
      # Keep a 2:1 physical viewport, matching original 600 x 300 images.
      old <- grid::grobTree(grid::rasterGrob(image, interpolate = FALSE),
                            vp = grid::viewport(width = grid::unit(8, "in"),
                                                height = grid::unit(4, "in")))
      # Drop the added publication legend in the direct recipe comparison only.
      left <- plots[[keys[[i]][1]]]$plot + ggplot2::theme(legend.position = "none")
      right <- plots[[keys[[i]][2]]]$plot
      new <- gridExtra::arrangeGrob(left, right, nrow = 1)
      sheet <- gridExtra::arrangeGrob(
        grid::textGrob(paste0(headings[i], "\nTOP: original R output | BOTTOM: regenerated from TSV"),
                       gp = grid::gpar(fontsize = 12)),
        old, new, ncol = 1, heights = c(0.6, 4.1, 4.1))
      grid::grid.newpage(); grid::grid.draw(sheet)
    }
  }
  write_review(temps[4])

  # Check completed PDF structure without introducing a new PDF package.
  check_pdf <- function(file) {
    size <- file.info(file)$size
    if (is.na(size) || size < 1000) stop("Empty/truncated PDF: ", file)
    con <- base::file(file, "rb"); on.exit(close(con))
    if (!identical(rawToChar(readBin(con, "raw", 5L)), "%PDF-")) stop("Invalid PDF header.")
    seek(con, where = max(0, size - 256), origin = "start")
    tail <- rawToChar(readBin(con, "raw", n = 256L))
    if (!grepl("%%EOF", tail, fixed = TRUE)) stop("PDF was not finalized: ", file)
    TRUE
  }
  invisible(vapply(temps[1:4], check_pdf, logical(1)))
  summary <- data.frame(
    view = names(plots), n_cells = vapply(plots, `[[`, numeric(1), "n_points"),
    n_unassigned = vapply(plots, `[[`, numeric(1), "n_unassigned"),
    coordinate_max_abs_difference = 0, plot_layer_check = "PASS",
    original_artwork_comparison = "PENDING_VISUAL_REVIEW",
    sample_S9B_mapping = "RNA-only inferred from supplied artwork; review required",
    stringsAsFactors = FALSE, row.names = NULL)
  utils::write.table(summary, temps[5], sep = "\t", quote = FALSE, row.names = FALSE)
  if (!overwrite && any(file.exists(targets))) stop("An output appeared; refusing to overwrite.")
  for (i in seq_along(targets)) {
    if (!file.rename(temps[i], targets[i])) stop("Could not finalize ", targets[i])
  }
  print(summary[, 1:6], row.names = FALSE)
  message("Wrote 3 panel PDFs, one 3-page comparison PDF, and one plot-check summary.")
  message("Original data and figures were untouched. Visual agreement is not yet approved.")
  invisible(summary)
}
