# Reproduce Fig 6H HOXA genotype coverage from historical DA results and fragments.
# Working directory: CH_publication_release (parent of multiome/).
# Historical checkpoints, resource workbook and multiome/figures/ are read-only.
# Does NOT rerun genotype differential accessibility or alter stored fragment files.

.f6h_get_tracks <- function(p) {
  # Signac CoveragePlot returns a nested patchwork of ggplots. Its coverage
  # ggplot stores the exact normalized, smoothed, downsampled points that are
  # rendered. Recursively find that layer without recalculating the signal.
  out <- list()
  walk <- function(q, depth = 0L) {
    if (depth > 8L || !inherits(q, "ggplot")) return(invisible(NULL))
    d <- tryCatch(q$data, error = function(e) NULL)
    if (is.data.frame(d) &&
        all(c("position", "coverage", "group") %in% names(d))) {
      out[[length(out) + 1L]] <<- d
    }
    children <- tryCatch(q$patches$plots, error = function(e) NULL)
    if (length(children)) for (child in children) walk(child, depth + 1L)
    invisible(NULL)
  }
  walk(p)
  if (!length(out)) stop(
    "Could not locate normalized coverage data within Signac plot. ",
    "No output files were written; inspect its ggplot/patchwork structure."
  )
  # There should be just one coverage track with three genotype facets.
  counts <- vapply(out, function(d) length(unique(as.character(d$group))), integer(1L))
  d <- out[[which.max(counts)]]
  if (!setequal(as.character(unique(d$group)), c("WT", "TET2", "DNMT3A"))) {
    stop("Extracted coverage plot has unexpected genotype groups: ",
         paste(unique(as.character(d$group)), collapse = ", "))
  }
  d
}

reproduce_Fig6H_HOXA <- function(obj = NULL, root = getwd(), output_root = NULL) {
  required <- c("SeuratObject", "Seurat", "Signac", "ggplot2", "openxlsx")
  absent <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Packages unavailable: ", paste(absent, collapse = ", "))

  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome")))
    stop("Run from CH_publication_release, not from multiome/.")
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  historical_figdir <- normalizePath(file.path(root, "multiome/figures"), mustWork = TRUE)
  if (identical(output_root, historical_figdir) ||
      startsWith(output_root, paste0(historical_figdir, "/"))) {
    stop("Refusing to write in historical multiome/figures.")
  }
  outfile <- file.path(output_root, "source_data/Fig6H_HOXA_genotype_coverage.tsv.gz")
  figfile <- file.path(output_root, "figures/Fig6H_HOXA_genotype_REVIEW.pdf")
  existing <- c(outfile, figfile)[file.exists(c(outfile, figfile))]
  if (length(existing)) stop("Existing publication outputs protected:\n",
                             paste(existing, collapse = "\n"))

  # Historical Fig 6H is the FIRST page of DNMT3A_coverage_plot.pdf.
  # Its original script takes the most significant positively changing peak
  # from padj-sorted results and extends 10,000 bp in both directions.
  workbook <- file.path(root, "..", "pipeline/resource/DApeaks_genotypes.xlsx")
  if (!file.exists(workbook)) stop("Historical workbook missing: ", workbook)
  da <- openxlsx::read.xlsx(workbook, sheet = "DNMT3A", colNames = TRUE)
  needed <- c("log2FoldChange", "padj")
  if (!all(needed %in% names(da))) stop("Expected DNMT3A results columns not found")
  ids <- as.character(da[[1L]])
  if (anyNA(ids) || anyDuplicated(ids)) stop("Missing or duplicate peak identifiers")
  positive <- which(is.finite(da$log2FoldChange) & da$log2FoldChange > 0)
  ranked <- positive[order(da$padj[positive], na.last = TRUE)]
  if (!length(ranked)) stop("No positive differential ATAC peaks")
  peak <- ids[ranked[1L]]
  pieces <- strsplit(peak, "-", fixed = TRUE)[[1L]]
  if (length(pieces) != 3L || pieces[1L] != "chr7" ||
      anyNA(suppressWarnings(as.numeric(pieces[2:3]))) ||
      as.numeric(pieces[2L]) < 27110000 ||
      as.numeric(pieces[3L]) > 27150000) {
    stop("First positive DNMT3A peak is not at the published HOXA locus: ", peak,
         ". Stop rather than select a peak by guesswork.")
  }
  cat("Historical first positively changing DNMT3A peak:", peak, "\n")

  if (is.null(obj)) {
    input <- file.path(root, "multiome/data/processed/multiome_seurat.rds")
    message("Loading historical Multiome Seurat checkpoint: ", input)
    obj <- readRDS(input)
  } else {
    message("Using supplied Seurat object; checkpoint not reloaded.")
  }
  if (!inherits(obj, "Seurat") || !"ATAC" %in% SeuratObject::Assays(obj))
    stop("A Seurat object with ATAC assay is required")

  meta_file <- file.path(root, "multiome/data/processed/metainfo.txt")
  meta <- read.delim(meta_file, row.names = 1L, check.names = FALSE)
  cells <- colnames(obj)
  if (length(cells) != 39787L || !setequal(cells, rownames(meta))) {
    stop("Seurat cells do not match the archived 39,787-cell metainfo table")
  }
  geno <- as.character(meta[cells, "genotype"])
  sample <- as.character(meta[cells, "sample"])
  valid <- c("WT", "TET2", "DNMT3A")
  if (any(!is.na(geno) & !geno %in% valid))
    stop("Unrecognized genotype in historical metainfo")
  print(table(sample, factor(geno, levels = valid), useNA = "ifany"))
  if (any(!valid %in% geno)) stop("Not all three genotype groups are represented")
  # Work only on this in-memory object; do not save or alter historical RDS.
  obj@meta.data$Fig6H_genotype <- factor(geno, levels = valid)

  # Relocate the four original Signac fragment objects on IRIS, verifying MD5.
  atac <- obj[["ATAC"]]
  frags <- Signac::Fragments(atac)
  if (length(frags) != 4L) stop("Expected four historical fragment objects")
  original_paths <- vapply(frags, function(f) slot(f, "path"), character(1L))
  samples <- sub(".*/multiome_larry_(C1|C2|D|T)/outs/atac_fragments\\.tsv\\.gz$",
                 "\\1", original_paths)
  if (!setequal(samples, c("C1", "C2", "D", "T")) || anyDuplicated(samples))
    stop("Cannot map the historical four fragments to samples")
  new_paths <- file.path(root, "..", "larry_atac", paste0("multiome_larry_", samples),
                         "outs/atac_fragments.tsv.gz")
  for (i in seq_along(frags)) {
    if (!file.exists(new_paths[i]) ||
        !file.exists(paste0(new_paths[i], ".tbi"))) {
      stop("Missing fragment or Tabix index: ", new_paths[i])
    }
    if (!file.exists(original_paths[i])) {
      message("Verifying relocated fragment: ", samples[i])
      frags[[i]] <- Signac::UpdatePath(frags[[i]], new.path = new_paths[i], verbose = FALSE)
    }
  }
  atac <- Signac::`Fragments<-`(atac, value = NULL)
  atac <- Signac::`Fragments<-`(atac, value = frags)
  obj[["ATAC"]] <- atac

  # The published panel combines historical C1_WT/C2_WT into WT and
  # C2_TET2/T_TET2 into TET2, retaining D_DNMT3A as DNMT3A.
  message("Building Signac coverage across the HOXA locus (one region only).")
  plot <- Signac::CoveragePlot(
    object = obj, assay = "ATAC", region = peak,
    group.by = "Fig6H_genotype", idents = valid,
    extend.upstream = 10000, extend.downstream = 10000,
    annotation = TRUE, peaks = TRUE, links = FALSE,
    window = 100
  )
  dat <- .f6h_get_tracks(plot)
  numeric_cols <- c("position", "coverage")
  if (!all(vapply(dat[numeric_cols], is.numeric, logical(1L))) ||
      any(!is.finite(dat$position)) || any(!is.finite(dat$coverage))) {
    stop("Nonfinite or nonnumeric plot values")
  }
  out <- data.frame(
    panel = "Fig6H", chromosome = "chr7",
    genomic_position = as.numeric(dat$position),
    genotype = as.character(dat$group),
    normalized_signal_plotted = as.numeric(dat$coverage),
    selected_DNMT3A_peak = peak,
    stringsAsFactors = FALSE
  )
  if ("norm.value" %in% names(dat))
    out$normalized_signal_before_smoothing <- as.numeric(dat$norm.value)
  out <- out[order(match(out$genotype, valid), out$genomic_position), , drop = FALSE]
  if (anyDuplicated(out[c("genotype", "genomic_position")]))
    stop("Duplicated plotted positions per genotype; inspect patchwork plot")
  if (nrow(out) < 300L) stop("Unexpectedly few coverage points")

  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(figfile), recursive = TRUE, showWarnings = FALSE)
  # Write exact numbers used in the rendered Signac track (not an approximation
  # from peak counts or a post-hoc re-normalization).
  con <- gzfile(outfile, "wt")
  tryCatch(utils::write.table(out, con, sep = "\t", row.names = FALSE,
                              quote = FALSE, na = "NA"),
           finally = close(con))
  readback <- utils::read.delim(gzfile(outfile), check.names = FALSE)
  if (nrow(readback) != nrow(out) ||
      !identical(as.character(readback$genotype), out$genotype) ||
      max(abs(readback$normalized_signal_plotted - out$normalized_signal_plotted),
          na.rm = TRUE) > 1e-10) {
    unlink(outfile)
    stop("Source Data round-trip check failed")
  }

  grDevices::pdf(figfile, width = 7.2, height = 5.5, useDingbats = FALSE)
  tryCatch(print(plot), finally = grDevices::dev.off())

  cat("Selected HOXA peak:", peak, "\n")
  cat("Genotype cells:\n")
  print(table(factor(geno, levels = valid)))
  cat("Plotted coverage points:", nrow(out), "\n")
  cat("Source Data:", outfile, "\n")
  cat("Review figure:", figfile, "\n")
  cat("Source Data round-trip: PASS; historical artwork comparison: REVIEW\n")
  invisible(list(peak = peak, genotype_counts = table(factor(geno, levels = valid)),
                 points = nrow(out), source_data = outfile, figure = figfile))
}
