# indels_calling/publication/run_publication.R
# CH publication reproduction: Supplementary Fig. S5G (TET2 and DNMT3A)
# Execute with getwd() == CH_publication_release.
# Historical algorithm adapted from indels_calling/scripts/1_processing.R.
# Original data, original script, and approved PDFs are never modified.
#
# IMPORTANT HISTORICAL DEFINITIONS:
#   Numerator = distinct (CB, UB) pairs with an I/D anywhere in an unspliced
#               alignment that spans the ENTIRE CRISPR target window.
#   Denominator = distinct UB sequences (not distinct CB,UB pairs) among all
#                 unspliced alignments that span the window.
#   These quantities are UMIs, although the historical axes say "reads".
# A CB+UB denominator is also exported for an explicit technical audit.
# The counting rules are intentionally NOT changed here.

.ic_specs <- data.frame(
  gene = c("TET2", "DNMT3A"),
  experiment = c(127L, 123L),
  chr = c("chr4", "chr2"),
  broad_start = c(105146934L, 25225874L),
  broad_end = c(105276332L, 25344483L),
  target_start = c(105234020L, 25247596L),
  target_end = c(105234039L, 25247612L),
  stringsAsFactors = FALSE
)

.ic_bam_path <- function(root, experiment, run) {
  file.path(root, "..", sprintf("human_larry_exp%d", experiment),
            sprintf("RNA%d_augment", run), "outs", "possorted_genome_bam.bam")
}

indels_preflight <- function(root = getwd(), mode = c("recount", "plot"),
                             checkpoint_path = NULL) {
  mode <- match.arg(mode)
  root <- normalizePath(root, mustWork = TRUE)
  metadata <- file.path(root, "RNA_analysis/data/processed/hashtag_meta.rds")
  if (!file.exists(metadata)) stop("Missing sample metadata: ", metadata, call. = FALSE)

  if (mode == "recount") {
    inputs <- do.call(rbind, lapply(seq_len(nrow(.ic_specs)), function(i) {
      s <- .ic_specs[i, ]
      do.call(rbind, lapply(1:4, function(run) {
        bam <- .ic_bam_path(root, s$experiment, run)
        data.frame(gene = s$gene, run = run,
                   BAM = file.exists(bam),
                   BAI = file.exists(paste0(bam, ".bai")) ||
                         file.exists(sub("\\.bam$", ".bai", bam)),
                   path = bam, stringsAsFactors = FALSE)
      }))
    }))
    print(inputs, row.names = FALSE)
    if (any(!inputs$BAM | !inputs$BAI)) {
      stop("Missing BAM or index; see preflight output.", call. = FALSE)
    }
  } else {
    if (is.null(checkpoint_path) || !file.exists(checkpoint_path)) {
      stop("Plot mode requires a valid checkpoint_path to an existing Source Data TSV.gz.",
           call. = FALSE)
    }
  }

  pkgs <- if (mode == "recount") {
    c("GenomicRanges", "IRanges", "S4Vectors", "Rsamtools",
      "GenomicAlignments", "data.table", "ggplot2")
  } else {
    "ggplot2"
  }
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  print(data.frame(package = pkgs, available = unname(ok)), row.names = FALSE)
  if (any(!ok)) {
    stop("Required package(s) missing: ", paste(pkgs[!ok], collapse = ", "),
         ". Install into the historical R/Bioconductor environment first.", call. = FALSE)
  }
  invisible(TRUE)
}

.ic_metadata_genotype <- function(meta, gene, experiment, run) {
  d <- meta[as.character(meta$exp) == as.character(experiment) &
            as.character(meta$RNA) == as.character(run), , drop = FALSE]
  g <- unique(as.character(d$genotype))
  g <- g[!is.na(g) & nzchar(g)]
  if (length(g) != 1L || !g %in% c(gene, "WT")) {
    stop("Unexpected genotype assignment: ", gene, " RNA", run,
         " (", paste(g, collapse = ", "), ")", call. = FALSE)
  }
  g
}

.ic_count_one <- function(root, spec, run, meta, control_label) {
  gene <- spec$gene
  geno <- .ic_metadata_genotype(meta, gene, spec$experiment, run)
  bam <- .ic_bam_path(root, spec$experiment, run)
  cat(sprintf("\n=== %s / RNA%d / %s ===\n", gene, run, geno))
  cat("Reading indexed BAM interval:", spec$chr,
      spec$broad_start, "-", spec$broad_end, "\n")

  # The same region, flags, grouping, UB/CB tags, and CIGAR rules as the
  # historical indels_calling/scripts/1_processing.R.
  targets <- GenomicRanges::GRanges(
    seqnames = spec$chr,
    ranges = IRanges::IRanges(start = spec$broad_start, end = spec$broad_end),
    strand = "-"
  )
  param <- Rsamtools::ScanBamParam(
    which = targets,
    tag = c("UB", "CB"),
    what = c("cigar", "pos", "qname"),
    flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE)
  )
  grouped <- GenomicAlignments::readGAlignmentsList(bam, param = param)
  flat <- unlist(grouped, use.names = FALSE)
  rm(grouped)
  if (!length(flat)) stop("No alignments in BAM interval: ", bam, call. = FALSE)
  tags <- S4Vectors::mcols(flat)
  if (is.null(tags$CB) || is.null(tags$UB)) {
    stop("Missing CB or UB tag from extracted alignments: ", bam, call. = FALSE)
  }

  dt <- data.table::data.table(
    Pos = as.integer(IRanges::start(flat)),
    Cigar = as.character(GenomicAlignments::cigar(flat)),
    CellBC = as.character(tags$CB),
    UB = as.character(tags$UB)
  )
  rm(flat, tags)
  dt <- dt[!is.na(UB) & UB != "" & !is.na(CellBC) & CellBC != ""]
  dt[, IsSplicedBridge := grepl("N", Cigar)]
  dt[, HasIndel := grepl("I", Cigar) | grepl("D", Cigar)]
  exonic <- dt[IsSplicedBridge == FALSE]
  rm(dt)

  exonic[, End := Pos + GenomicAlignments::cigarWidthAlongReferenceSpace(Cigar) - 1L]
  robust <- exonic[Pos <= spec$target_start & End >= spec$target_end]
  rm(exonic)
  if (!nrow(robust)) stop("No non-spliced reads span target site: ", gene,
                            " RNA", run, call. = FALSE)

  # Preserve original numerator and denominator *exactly*.
  edited <- data.table::uniqueN(robust[HasIndel == TRUE, .(CellBC, UB)])
  total_ub_only <- data.table::uniqueN(robust$UB)
  # Additional audit, not used for the published figure:
  total_cell_ub <- data.table::uniqueN(robust[, .(CellBC, UB)])
  spanning_reads <- nrow(robust)
  reads_with_any_indel <- sum(robust$HasIndel)
  if (edited > total_ub_only) {
    warning("Historical numerator exceeds historical denominator for ",
            gene, " RNA", run, "; this may reflect CB/UB counting asymmetry.")
  }
  cat(sprintf("historical edited=%d; historical total=%d; CB+UB total=%d; spanning alignments=%d\n",
              edited, total_ub_only, total_cell_ub, spanning_reads))

  data.frame(
    gene = gene,
    experiment = as.integer(spec$experiment),
    RNA_run = as.integer(run),
    metadata_genotype = geno,
    display_genotype = if (geno == "WT") control_label else geno,
    chromosome = spec$chr,
    target_start = as.integer(spec$target_start),
    target_end = as.integer(spec$target_end),
    historical_edited_CB_UB_any_indel = as.integer(edited),
    historical_total_UB_only = as.integer(total_ub_only),
    audited_total_CB_UB = as.integer(total_cell_ub),
    spanning_alignments = as.integer(spanning_reads),
    spanning_alignments_with_any_indel = as.integer(reads_with_any_indel),
    historical_edited_fraction = edited / total_ub_only,
    audited_CB_UB_fraction = edited / total_cell_ub,
    BAM_relative_path = file.path("..", sprintf("human_larry_exp%d", spec$experiment),
                                  sprintf("RNA%d_augment", run), "outs",
                                  "possorted_genome_bam.bam"),
    stringsAsFactors = FALSE
  )
}

.ic_write_tsv_gz <- function(df, path) {
  con <- gzfile(path, open = "wt")
  tryCatch(
    utils::write.table(df, file = con, sep = "\t", row.names = FALSE,
                       col.names = TRUE, quote = FALSE, na = "NA"),
    finally = close(con)
  )
  invisible(path)
}

.ic_plot <- function(df, gene) {
  s <- .ic_specs[.ic_specs$gene == gene, ]
  dat <- df[df$gene == gene, , drop = FALSE]
  ctrl <- unique(dat$display_genotype[dat$metadata_genotype == "WT"])
  if (nrow(dat) != 4L || length(ctrl) != 1L) {
    stop("Expected four observations and one WT display label for ", gene,
         call. = FALSE)
  }
  dat$display_genotype <- factor(dat$display_genotype,
                                 levels = c(gene, ctrl))
  palette <- stats::setNames(c("#F8766D", "#00BFC4"), c(gene, ctrl))
  title <- sprintf("%s: %s:%s-%s", gene, s$chr,
                   format(s$target_start, big.mark = ",", scientific = FALSE, trim = TRUE),
                   format(s$target_end, big.mark = ",", scientific = FALSE, trim = TRUE))
  ggplot2::ggplot(dat, ggplot2::aes(x = historical_total_UB_only,
                                   y = historical_edited_CB_UB_any_indel,
                                   fill = display_genotype)) +
    ggplot2::geom_abline(intercept = 0, slope = 1,
                         linetype = "dashed", size = 0.6) +
    ggplot2::geom_point(shape = 21, size = 3, colour = "black") +
    ggplot2::scale_fill_manual(values = palette, breaks = c(gene, ctrl),
                               drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq(0, 250, 50)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 250, 50)) +
    ggplot2::coord_cartesian(xlim = c(0, 250), ylim = c(0, 250),
                             expand = FALSE) +
    ggplot2::labs(x = "total reads", y = "edited reads",
                  title = title, fill = "genotype") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      legend.position = c(0.03, 0.97),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill = "white", colour = "grey60"),
      panel.grid.minor = ggplot2::element_line(colour = "grey94")
    )
}

.ic_write_pdf <- function(df, path) {
  ps <- lapply(c("TET2", "DNMT3A"), function(gene) .ic_plot(df, gene))
  grDevices::pdf(path, width = 10.4, height = 4.8, useDingbats = FALSE)
  tryCatch({
    grid::grid.newpage()
    print(ps[[1]], vp = grid::viewport(x = 0.25, y = 0.5,
                                       width = 0.48, height = 0.97))
    print(ps[[2]], vp = grid::viewport(x = 0.75, y = 0.5,
                                       width = 0.48, height = 0.97))
  }, finally = grDevices::dev.off())
  invisible(path)
}

# mode="recount": calculate all 8 rows independently from indexed BAMs.
# mode="plot": use an existing Source Data TSV.gz to re-create only the PDF.
# Output never overwrites existing Source Data or PDF files.
run_indels_publication <- function(
    mode = c("recount", "plot"),
    root = getwd(),
    output_root = file.path(root, "indels_calling", "publication"),
    checkpoint_path = NULL,
    control_label = "WT") {
  mode <- match.arg(mode)
  if (!is.character(control_label) || length(control_label) != 1L ||
      is.na(control_label) || !nzchar(control_label)) {
    stop("control_label must be a single nonempty label.", call. = FALSE)
  }
  root <- normalizePath(root, mustWork = TRUE)
  indels_preflight(root = root, mode = mode, checkpoint_path = checkpoint_path)

  source_path <- file.path(output_root, "source_data",
                           "FigS5G_target_site_editing_counts.tsv.gz")
  figure_path <- file.path(output_root, "figures",
                           "FigS5G_target_site_editing_REVIEW.pdf")
  if (file.exists(source_path) || file.exists(figure_path)) {
    stop("Target Source Data or PDF exists; choose a fresh output_root. Files preserved.",
         call. = FALSE)
  }

  if (mode == "recount") {
    meta <- readRDS(file.path(root, "RNA_analysis/data/processed/hashtag_meta.rds"))
    if (!is.data.frame(meta) || !all(c("exp", "RNA", "genotype") %in% names(meta))) {
      stop("Unexpected hashtag_meta.rds structure.", call. = FALSE)
    }
    rows <- list()
    for (i in seq_len(nrow(.ic_specs))) {
      spec <- .ic_specs[i, ]
      for (run in 1:4) {
        rows[[length(rows) + 1L]] <- .ic_count_one(root, spec, run, meta, control_label)
      }
    }
    data <- do.call(rbind, rows)
  } else {
    data <- utils::read.delim(gzfile(checkpoint_path), stringsAsFactors = FALSE,
                              check.names = FALSE)
    required_cols <- c("gene", "RNA_run", "metadata_genotype", "display_genotype",
                       "historical_total_UB_only", "historical_edited_CB_UB_any_indel")
    if (!all(required_cols %in% names(data))) {
      stop("Input checkpoint is missing expected columns.", call. = FALSE)
    }
  }
  if (nrow(data) != 8L ||
      !identical(as.character(data$gene), rep(c("TET2", "DNMT3A"), each = 4L)) ||
      !identical(as.integer(data$RNA_run), rep(1:4, times = 2L))) {
    stop("Expected exactly 8 rows ordered by TET2, DNMT3A and RNA1-4.",
         call. = FALSE)
  }

  dir.create(dirname(source_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
  .ic_write_tsv_gz(data, source_path)
  .ic_write_pdf(data, figure_path)
  cat("\nPASS: Fig S5G; 8 sample-level rows; one two-panel PDF.\n")
  cat("Source Data:", source_path, "\nPDF:", figure_path, "\n")
  invisible(list(source_data = source_path, figure = figure_path, counts = data))
}
