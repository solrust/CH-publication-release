# Fig 6G: recreate the published motif/clone-phenotype dot plot from its
# historical workbook (multiome/scripts/5_clone_analysis_final.R, lines 366-406).
# Run with CH_publication_release as the working directory.
# No Seurat object, model fitting, or historical file modifications.

reproduce_Fig6G_from_workbook <- function(
    root = getwd(),
    workbook = NULL,
    output_root = NULL) {

  for (pkg in c("openxlsx", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Missing package: ", pkg, call. = FALSE)
  }
  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome")))
    stop("Run from CH_publication_release, not multiome/", call. = FALSE)
  if (is.null(workbook)) {
    workbook <- file.path(root, "pipeline", "resource",
                          "motif_associations.xlsx")
  }
  if (!file.exists(workbook))
    stop("Historical workbook not found: ", workbook, call. = FALSE)
  if (is.null(output_root))
    output_root <- file.path(root, "multiome", "publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  historical_figures <- normalizePath(file.path(root, "multiome", "figures"),
                                      mustWork = TRUE)
  if (identical(output_root, historical_figures) ||
      startsWith(output_root, paste0(historical_figures, "/"))) {
    stop("Refusing to write to historical figures directory", call. = FALSE)
  }

  source_path <- file.path(output_root, "source_data",
                           "Fig6G_motif_clone_phenotype_associations.tsv.gz")
  figure_path <- file.path(output_root, "figures",
                           "Fig6G_motif_phenotypes_REVIEW.pdf")
  if (any(file.exists(c(source_path, figure_path)))) {
    stop("Existing publication output(s) protected:\n",
         paste(c(source_path, figure_path)[file.exists(c(source_path, figure_path))],
               collapse = "\n"), call. = FALSE)
  }

  # Same eight names and ordering as pheno.names.print in original script.
  phenotype_order <- c("HSC bias", "lymphoid bias", "dendritic bias",
                       "myeloid bias", "erythroid bias", "mk bias",
                       "growth rate", "clone size")
  sheets <- openxlsx::getSheetNames(workbook)
  missing_sheets <- setdiff(phenotype_order, sheets)
  if (length(missing_sheets))
    stop("Workbook missing phenotype sheet(s): ",
         paste(missing_sheets, collapse = ", "), call. = FALSE)

  cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  by_phenotype <- setNames(lapply(phenotype_order, function(pheno) {
    dat <- openxlsx::read.xlsx(workbook, sheet = pheno, colNames = TRUE)
    # The first Excel column contains motif names; its header is blank.
    missing <- setdiff(cols, names(dat))
    if (length(missing))
      stop("Missing columns in sheet '", pheno, "': ",
           paste(missing, collapse = ", "), call. = FALSE)
    ids <- trimws(as.character(dat[[1]]))
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop("Missing/duplicated motif identifiers in sheet: ", pheno,
           call. = FALSE)
    result <- data.frame(motif = ids, phenotype = pheno,
                         stringsAsFactors = FALSE)
    for (nm in cols) {
      if (!is.numeric(dat[[nm]]))
        stop("Non-numeric column ", nm, " in sheet: ", pheno,
             call. = FALSE)
      result[[nm]] <- dat[[nm]]
    }
    result
  }), phenotype_order)

  motif_order <- by_phenotype[[1]]$motif
  for (pheno in phenotype_order) {
    ids <- by_phenotype[[pheno]]$motif
    if (!setequal(ids, motif_order)) {
      stop("Motif sets differ among sheets (first mismatch: ", pheno,
           "). Review workbook before plotting.", call. = FALSE)
    }
    by_phenotype[[pheno]] <- by_phenotype[[pheno]][
      match(motif_order, ids), , drop = FALSE]
  }

  full <- do.call(rbind, by_phenotype)
  rownames(full) <- NULL
  to_matrix <- function(field) {
    x <- do.call(cbind, lapply(by_phenotype, function(d) d[[field]]))
    rownames(x) <- motif_order
    colnames(x) <- phenotype_order
    x
  }
  fc <- to_matrix("log2FoldChange")
  pv <- to_matrix("pvalue")
  padj <- to_matrix("padj")

  # Original plot: select motifs by min adjusted P across ALL eight
  # phenotypes (including mk bias), then omit mk bias from display.
  selected <- apply(padj, 1, function(x) any(is.finite(x) & x < 0.2))
  display_phenotypes <- setdiff(phenotype_order, "mk bias")
  full$selected_by_any_FDR_0.2 <- unname(selected[full$motif])
  full$plotted_in_Fig6G <- full$selected_by_any_FDR_0.2 &
    full$phenotype %in% display_phenotypes
  if (sum(selected) < 2L) stop("Too few selected motifs to cluster", call. = FALSE)

  plot_fc <- fc[selected, display_phenotypes, drop = FALSE]
  if (any(!is.finite(plot_fc)))
    stop("Selected motif effects contain nonfinite values; review workbook",
         call. = FALSE)
  motif_levels <- rownames(plot_fc)[stats::hclust(stats::dist(plot_fc))$order]
  phenotype_levels <- colnames(plot_fc)[
    stats::hclust(stats::dist(t(plot_fc)))$order]

  plot_dat <- full[full$plotted_in_Fig6G, , drop = FALSE]
  plot_dat$motif <- factor(plot_dat$motif, levels = motif_levels)
  plot_dat$phenotype <- factor(plot_dat$phenotype, levels = phenotype_levels)
  plot_dat$direction <- ifelse(plot_dat$log2FoldChange > 0, "up", "down")
  # Display only: avoid infinite alpha for P values rounded to zero.
  plot_dat$minus_log10_p <- -log10(pmax(plot_dat$pvalue,
                                        .Machine$double.xmin))
  if (any(!is.finite(plot_dat$minus_log10_p)))
    stop("Selected plotting P values are missing/nonfinite", call. = FALSE)
  plot_dat$p_below_0.01 <- plot_dat$pvalue < 0.01

  p <- ggplot2::ggplot(plot_dat) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(
      x = motif, y = phenotype, size = abs(log2FoldChange),
      stroke = p_below_0.01, fill = direction, alpha = minus_log10_p
    ), colour = "black", shape = 21) +
    ggplot2::scale_size_continuous(range = c(0, 8)) +
    ggplot2::scale_fill_manual(values = c(up = "red", down = "blue")) +
    ggplot2::labs(x = "TF motif", y = "HSC phenotype") +
    ggplot2::guides(size = "none", alpha = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 20)
    )

  dir.create(dirname(source_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(source_path, open = "wt")
  tryCatch(utils::write.table(full, con, sep = "\t", quote = FALSE,
                              row.names = FALSE, na = "NA"),
           finally = close(con))
  reloaded <- utils::read.delim(gzfile(source_path), check.names = FALSE)
  if (!identical(full$motif, reloaded$motif) ||
      !identical(full$phenotype, reloaded$phenotype) ||
      !isTRUE(all.equal(full$log2FoldChange, reloaded$log2FoldChange,
                        tolerance = 1e-12, check.attributes = FALSE)) ||
      !isTRUE(all.equal(full$pvalue, reloaded$pvalue,
                        tolerance = 1e-12, check.attributes = FALSE))) {
    stop("Source Data roundtrip verification failed", call. = FALSE)
  }

  grDevices::pdf(figure_path, width = 10, height = 5,
                 useDingbats = FALSE)
  tryCatch(print(p), finally = grDevices::dev.off())

  summary <- data.frame(
    total_motifs = length(motif_order),
    phenotypes_in_workbook = length(phenotype_order),
    motifs_selected_FDR_0.2 = sum(selected),
    phenotypes_displayed = length(display_phenotypes),
    points_plotted = nrow(plot_dat),
    source_data_roundtrip = "PASS",
    visual_agreement = "REVIEW"
  )
  print(summary, row.names = FALSE)
  cat("Source Data:", source_path, "\n")
  cat("Figure:", figure_path, "\n")
  invisible(summary)
}
