#!/usr/bin/env Rscript
# CH manuscript Figure 3D/3E: canonical render from accepted-panel pathway selections.
# Requires the canonical TSVs from export_Fig3DE.R, NOT mutant_de.RData.
# Run from the root of CH_publication_release.

root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "RNA_analysis", "scripts", "6_mutant_analysis.R"))) {
  stop("Run this script from the CH_publication_release repository root.")
}

for (package in c("ggplot2", "ggrepel", "gridExtra", "RColorBrewer")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package not installed: ", package)
  }
}

indir <- file.path(root, "RNA_analysis", "publication", "source_data")
outdir <- file.path(root, "RNA_analysis", "publication", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 3_sc_processing.R defined clust.final.cols using brewer.pal(10, "Paired"),
# with final cell state 1 = HSC. Avoid loading 7.4 GiB scrna.RData just for color.
hsc_color <- unname(RColorBrewer::brewer.pal(10, "Paired")[[1]])
pcut <- 1e-3  # Original 6_mutant_analysis.R labeling threshold (raw DE p-value)

for (genotype in c("DNMT3A", "TET2")) {
  panel <- if (genotype == "DNMT3A") "Fig3D" else "Fig3E"
  genes_file <- file.path(indir, paste0(panel, "_genes.tsv"))
  pathways_file <- file.path(indir, paste0(panel, "_pathways.tsv"))
  if (!file.exists(genes_file) || !file.exists(pathways_file)) {
    stop("Source Data TSVs missing for ", panel,
         ". Run export_Fig3DE.R first.")
  }
  genes <- utils::read.delim(genes_file, check.names = FALSE,
                             stringsAsFactors = FALSE, na.strings = "NA")
  pathways <- utils::read.delim(pathways_file, check.names = FALSE,
                                stringsAsFactors = FALSE, na.strings = "NA")
  stopifnot(all(c("gene_id", "log2FoldChange", "pvalue",
                  "volcano_label_candidate") %in% names(genes)))
  stopifnot(all(c("pathway", "NES", "pval", "plot_order",
                  "signed_neglog10_raw_pvalue") %in% names(pathways)))

  # Keep the historical mapping: raw DESeq2 p-value, log10 axis, shrunken LFC.
  valid_p <- genes$pvalue[is.finite(genes$pvalue) & genes$pvalue > 0]
  if (!length(valid_p)) stop("No positive DE p-values in ", panel)
  x_annotation <- max(genes$log2FoldChange, na.rm = TRUE) * 0.6

  # Scientific quantities and most plotting styles follow 6_mutant_analysis.R,
  # lines 27-83. Unlike the historical code, draw the cluster annotation ONCE.
  # Text label positions can vary slightly across ggrepel/ggplot2 versions.
  p1 <- ggplot2::ggplot(
    genes,
    ggplot2::aes(x = log2FoldChange, y = pvalue,
                 label = volcano_label_candidate,
                 alpha = I(ifelse(is.finite(pvalue) & pvalue <= pcut, 1, 0.1)))
  ) +
    ggplot2::geom_point(size = 2) +
    ggrepel::geom_text_repel(size = 4, max.overlaps = 10, seed = 12345) +
    ggplot2::annotate("point", x = x_annotation, y = 0.4, size = 20,
                      shape = 21, fill = hsc_color, alpha = 1) +
    ggplot2::annotate("text", x = x_annotation, y = 0.4,
                      label = "1", size = 10) +
    ggplot2::labs(x = paste(genotype, "/ WT, log"), y = "p-value") +
    ggplot2::scale_y_log10(limits = c(min(valid_p), 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20))

  if (!setequal(pathways$plot_order, seq_len(nrow(pathways)))) {
    stop("Invalid pathway plot order in ", panel)
  }
  ordered_levels <- pathways$pathway[order(pathways$plot_order)]
  pathways$pathway <- factor(pathways$pathway, levels = ordered_levels)
  pathways$direction_sign <- factor(sign(pathways$NES))

  p2 <- ggplot2::ggplot(
    pathways,
    ggplot2::aes(x = pathway, y = signed_neglog10_raw_pvalue,
                 fill = direction_sign)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(),
                      color = grDevices::adjustcolor("black", 1)) +
    ggplot2::scale_fill_manual(values = c(
      "-1" = grDevices::adjustcolor("blue", 0.2),
      "1"  = grDevices::adjustcolor("red", 0.2)
    )) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "p-value, log") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   axis.text.x = ggplot2::element_text(size = 20),
                   axis.text.y = ggplot2::element_text(size = 10),
                   legend.position = "none")

  pdf_file <- file.path(outdir, paste0(panel, "_reproduced.pdf"))
  grDevices::pdf(pdf_file, width = 12, height = 6, useDingbats = FALSE)
  tryCatch(
    gridExtra::grid.arrange(gridExtra::arrangeGrob(p1),
                            gridExtra::arrangeGrob(p2), widths = c(1, 1)),
    finally = grDevices::dev.off()
  )
  message("Wrote: ", pdf_file)
}

capture.output(utils::sessionInfo(),
               file = file.path(outdir, "Fig3DE_plot_sessionInfo.txt"))
message("Canonical Figure 3D/E PDFs generated from accepted-artwork Source Data; review before replacing submitted figures.")
