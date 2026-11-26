# Supplementary Fig. 10A and 10F - plot historical motif-association volcanoes.
# Run with CH_publication_release as working directory.
# Reuses the already-validated numeric exports for Fig 6G and Fig 6I.
# Does not recompute DA/motif associations, load Seurat, or overwrite historical outputs.
# Sources: multiome/scripts/5_clone_analysis_final.R, lines 350-363 and 479-491.

plot_FigS10AF_from_existing_source_data <- function(root = getwd(),
                                                     output_root = NULL,
                                                     max_labels_per_panel = 40L) {
  needed <- c("ggplot2", "ggrepel", "gridExtra", "scales")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing plotting package(s): ", paste(absent, collapse = ", "), call. = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome", "publication")))
    stop("Run from CH_publication_release.", call. = FALSE)
  if (is.null(output_root)) output_root <- file.path(root, "multiome", "publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  historic <- normalizePath(file.path(root, "multiome", "figures"), mustWork = TRUE)
  if (identical(output_root, historic) || startsWith(output_root, paste0(historic, "/")))
    stop("Refusing to write under the historical figures directory.", call. = FALSE)

  source_dir <- file.path(root, "multiome", "publication", "source_data")
  src_a <- file.path(source_dir, "Fig6G_motif_clone_phenotype_associations.tsv.gz")
  src_f <- file.path(source_dir, "Fig6I_motif_genotype_associations.tsv.gz")
  if (!all(file.exists(c(src_a, src_f)))) stop("Missing validated upstream Source Data:\n",
    paste(c(src_a, src_f)[!file.exists(c(src_a, src_f))], collapse = "\n"), call. = FALSE)

  output_a <- file.path(output_root, "figures", "FigS10A_motif_phenotype_volcano_REVIEW.pdf")
  output_f <- file.path(output_root, "figures", "FigS10F_motif_genotype_volcano_REVIEW.pdf")
  if (any(file.exists(c(output_a, output_f)))) stop("Existing figure output(s) protected:\n",
      paste(c(output_a, output_f)[file.exists(c(output_a, output_f))], collapse = "\n"), call. = FALSE)

  max_labels_per_panel <- as.integer(max_labels_per_panel)
  if (is.na(max_labels_per_panel) || max_labels_per_panel < 0L)
    stop("max_labels_per_panel must be a nonnegative integer.", call. = FALSE)

  # === S10A: seven displayed clone phenotypes; mk bias remains in Fig 6G Source Data. ===
  a <- utils::read.delim(src_a, check.names = FALSE, stringsAsFactors = FALSE)
  required_a <- c("motif", "phenotype", "log2FoldChange", "pvalue", "padj")
  if (!all(required_a %in% names(a))) stop("Fig 6G Source Data lacks expected fields.", call. = FALSE)
  if (anyNA(a$motif) || any(!nzchar(a$motif)) || anyDuplicated(a[, c("motif", "phenotype")]))
    stop("Missing or repeated phenotype/motif identifiers in Fig 6G Source Data.", call. = FALSE)
  all_phenos <- c("HSC bias", "lymphoid bias", "dendritic bias", "myeloid bias",
                  "erythroid bias", "mk bias", "growth rate", "clone size")
  if (!setequal(unique(a$phenotype), all_phenos))
    stop("Unexpected phenotype set in Fig 6G Source Data.", call. = FALSE)
  display_a <- setdiff(all_phenos, "mk bias")
  n_motifs_a <- length(unique(a$motif))
  if (nrow(a) != n_motifs_a * length(all_phenos))
    stop("Incomplete motif-by-phenotype grid in Fig 6G Source Data.", call. = FALSE)

  choose_labels <- function(d, limit) {
    eligible <- which(is.finite(d$pvalue) & is.finite(d$log2FoldChange) & d$pvalue < 0.05)
    if (length(eligible) && limit > 0L) {
      eligible <- eligible[order(d$pvalue[eligible])]
      d$label <- ""
      d$label[head(eligible, limit)] <- d$motif[head(eligible, limit)]
    } else d$label <- ""
    d
  }
  p_a <- lapply(display_a, function(pheno) {
    d <- a[a$phenotype == pheno, , drop = FALSE]
    d <- d[is.finite(d$log2FoldChange) & is.finite(d$pvalue) & d$pvalue >= 0, , drop = FALSE]
    d$p_plot <- pmax(d$pvalue, .Machine$double.xmin)
    d <- choose_labels(d, max_labels_per_panel)
    ggplot2::ggplot(d, ggplot2::aes(x = log2FoldChange, y = p_plot)) +
      ggplot2::geom_point(size = 0.65, color = "black", alpha = 0.85) +
      ggrepel::geom_text_repel(ggplot2::aes(label = label),
                               size = 2.0, seed = 61010, max.overlaps = Inf,
                               box.padding = 0.12, point.padding = 0.05,
                               min.segment.length = 0, na.rm = TRUE) +
      ggplot2::scale_y_log10() +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::labs(title = pheno, x = "logFC", y = "p-value") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5))
  })

  # === S10F: the already exported Fig 6I table includes both genotype results
  # and motif GC content. No Seurat checkpoint reload is necessary. ===
  f <- utils::read.delim(src_f, check.names = FALSE, stringsAsFactors = FALSE)
  required_f <- c("motif", "motif_gc_percent", "DNMT3A_log2FoldChange",
                  "DNMT3A_pvalue", "DNMT3A_padj", "TET2_log2FoldChange",
                  "TET2_pvalue", "TET2_padj")
  if (!all(required_f %in% names(f))) stop("Fig 6I Source Data lacks expected fields.", call. = FALSE)
  if (anyNA(f$motif) || any(!nzchar(f$motif)) || anyDuplicated(f$motif) ||
      any(!is.finite(f$motif_gc_percent)))
    stop("Invalid motif identifiers/GC metadata in Fig 6I Source Data.", call. = FALSE)
  gc <- f$motif_gc_percent
  if (diff(range(gc)) <= 0) stop("Invalid GC range", call. = FALSE)
  gc_breaks <- scales::rescale(c(min(gc), mean(gc), max(gc)))

  p_f <- lapply(c("DNMT3A", "TET2"), function(geno) {
    d <- data.frame(motif = f$motif,
                    motif_gc_percent = gc,
                    log2FoldChange = f[[paste0(geno, "_log2FoldChange")]],
                    pvalue = f[[paste0(geno, "_pvalue")]],
                    padj = f[[paste0(geno, "_padj")]],
                    stringsAsFactors = FALSE)
    d <- d[is.finite(d$log2FoldChange) & is.finite(d$pvalue) & d$pvalue >= 0, , drop = FALSE]
    d$p_plot <- pmax(d$pvalue, .Machine$double.xmin)
    d <- choose_labels(d, max_labels_per_panel)
    nsig <- sum(!is.na(d$padj) & d$padj <= 0.2)
    nonsig <- d$pvalue[is.finite(d$padj) & d$padj > 0.2 & is.finite(d$pvalue) & d$pvalue > 0]
    threshold <- if (nsig > 1L && length(nonsig)) min(nonsig) else NA_real_
    p <- ggplot2::ggplot(d, ggplot2::aes(
      x = log2FoldChange, y = p_plot, colour = motif_gc_percent)) +
      ggplot2::geom_point(size = 0.8, alpha = 0.85) +
      ggrepel::geom_text_repel(ggplot2::aes(label = label),
                               size = 2.1, seed = 61010, max.overlaps = Inf,
                               box.padding = 0.12, point.padding = 0.05,
                               min.segment.length = 0, na.rm = TRUE) +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_gradientn(
        colours = c("blue", "grey70", "darkred"), values = gc_breaks,
        limits = range(gc), name = "Motif GC") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = paste0(geno, " vs WT TF motif association"),
                    x = "logFC", y = "p-value")
    if (is.finite(threshold))
      p <- p + ggplot2::geom_hline(yintercept = threshold, linetype = 2)
    p
  })

  dir.create(dirname(output_a), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output_a, width = 16, height = 9, useDingbats = FALSE)
  tryCatch({
    gridExtra::grid.arrange(grobs = p_a, layout_matrix = rbind(c(1, 2, 3, NA), c(4, 5, 6, 7)))
  }, finally = grDevices::dev.off())

  grDevices::pdf(output_f, width = 12, height = 6, useDingbats = FALSE)
  tryCatch(gridExtra::grid.arrange(grobs = p_f, ncol = 2),
           finally = grDevices::dev.off())

  result <- data.frame(
    panel = c("FigS10A", "FigS10F"),
    n_panels = c(length(p_a), length(p_f)),
    n_motifs = c(n_motifs_a, nrow(f)),
    source_data_status = "REUSED_VALIDATED",
    artwork_comparison = "REVIEW",
    stringsAsFactors = FALSE
  )
  print(result, row.names = FALSE)
  cat("Fig S10A PDF:", output_a, "\n")
  cat("Fig S10F PDF:", output_f, "\n")
  cat("Shared S10A numerical source:", src_a, "\n")
  cat("Shared S10F numerical source:", src_f, "\n")
  invisible(result)
}
