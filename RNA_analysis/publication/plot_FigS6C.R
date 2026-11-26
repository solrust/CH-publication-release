#!/usr/bin/env Rscript
# Recreate Supplementary Figure 6C from the exported TSV files alone.
# Does NOT load any large RData object or rerun differential expression.
# Run from CH_publication_release root.

root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "RNA_analysis", "scripts", "5_mutant_DE.R"))) {
  stop("Please run from the CH_publication_release repository root.")
}
for (pkg in c("ggplot2", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
}

indir <- file.path(root, "RNA_analysis", "publication", "source_data")
outdir <- file.path(root, "RNA_analysis", "publication", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
input <- file.path(indir, "FigS6C_DE_concordance.tsv.gz")
corr_input <- file.path(indir, "FigS6C_correlations.tsv")
if (!file.exists(input) || !file.exists(corr_input)) {
  stop("Missing FigS6C Source Data. Run export_FigS6C.R first.")
}

con <- gzfile(input, open = "rt")
x <- tryCatch(
  utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
                    na.strings = "NA"),
  finally = close(con)
)
corr <- utils::read.delim(corr_input, stringsAsFactors = FALSE)
needed <- c("cluster", "gene_id", "log2FoldChange_DNMT3A",
            "log2FoldChange_TET2", "plot_group", "plot_label_candidate")
if (!all(needed %in% names(x))) stop("Missing columns in FigS6C source table")
if (!all(c("cluster", "displayed_r", "label_x", "label_y") %in% names(corr))) {
  stop("Missing columns in FigS6C correlation table")
}
if (!setequal(unique(x$cluster), 1:10) ||
    !identical(sort(as.integer(corr$cluster)), 1:10)) {
  stop("Expected ten states in both FigS6C tables")
}

x$cluster <- factor(x$cluster, levels = 1:10)
corr$cluster <- factor(corr$cluster, levels = 1:10)
corr$r_label <- format(corr$displayed_r, nsmall = 2, trim = TRUE)
# Only text-bearing rows are needed for ggrepel: faster and less memory use.
labels <- x[!is.na(x$plot_label_candidate), , drop = FALSE]

# Match the original manual color assignment in 3_sc_processing.R:
# col.genotype <- adjustcolor(c('black', 'red', 'green'), 0.5)
# A genuine unclassified point is grey.  Legend omits the NA class.
colors <- c(
  "both" = "black",
  "DNMT3A" = grDevices::adjustcolor("red", alpha.f = 0.5),
  "TET2" = grDevices::adjustcolor("green", alpha.f = 0.5)
)
# Preserve the historical NA/unclassified genes as a grey background layer.
# In v1, na.translate = FALSE silently removed these points.
x$draw_group <- ifelse(is.na(x$plot_group), "background", as.character(x$plot_group))
x$draw_group <- factor(x$draw_group, levels = c("background", "both", "DNMT3A", "TET2"))

p <- ggplot2::ggplot() +
  # Draw all other genes first, beneath the highlighted differential genes.
  ggplot2::geom_point(
    data = x[x$draw_group == "background", , drop = FALSE],
    ggplot2::aes(x = log2FoldChange_DNMT3A,
                 y = log2FoldChange_TET2),
    color = "grey75", size = 0.7, alpha = 0.7, na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = x[x$draw_group != "background", , drop = FALSE],
    ggplot2::aes(x = log2FoldChange_DNMT3A,
                 y = log2FoldChange_TET2,
                 color = draw_group),
    size = 1, alpha = 0.75, na.rm = TRUE
  ) +
  ggplot2::geom_text(
    data = corr,
    ggplot2::aes(x = label_x, y = label_y, label = r_label),
    size = 10
  ) +
  ggrepel::geom_text_repel(
    data = labels,
    ggplot2::aes(x = log2FoldChange_DNMT3A,
                 y = log2FoldChange_TET2,
                 label = plot_label_candidate),
    seed = 12345, max.overlaps = 10
  ) +
  ggplot2::scale_color_manual(
    values = colors, breaks = c("both", "DNMT3A", "TET2"),
    drop = FALSE
  ) +
  ggplot2::facet_wrap(~cluster, nrow = 2, scales = "free") +
  ggplot2::labs(x = "DNMT3A, logFC", y = "TET2, logFC", color = "DE") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 20),
    legend.text = ggplot2::element_text(size = 20),
    axis.title = ggplot2::element_text(size = 30)
  )

outfile <- file.path(outdir, "FigS6C_reproduced.pdf")
grDevices::pdf(outfile, width = 20, height = 8, useDingbats = FALSE)
tryCatch(
  print(p),  # Necessary when a ggplot is sent to a PDF device from source()
  finally = grDevices::dev.off()
)
capture.output(utils::sessionInfo(),
               file = file.path(outdir, "FigS6C_sessionInfo.txt"))
message("Wrote: ", outfile)
message("Review against Fig. S6C in figure_S6_V2.ai before marking as validated.")
