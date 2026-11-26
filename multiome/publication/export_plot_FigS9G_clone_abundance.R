# Supplementary Figure 9G: per-clone abundance in scRNA-seq and Multiome.
# Input files are read-only; designed for the CH_publication_release working directory.
# Does not load the large multiome Seurat checkpoint or run clone calling.
# This is a first-pass RECONSTRUCTION of the plotted quantities; historical
# plotting code from RNA_clone_analysis has not yet been confirmed.
# Multiple Multiome clone records may map to one RNA clone (per historical lookup).
# The source table reports multiplicity and does not collapse these records.

reproduce_FigS9G_clone_abundance <- function(
    root = getwd(),
    output_root = NULL,
    min_clones = 20L,
    require_size_consistency = TRUE
) {
  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome")) ||
      !dir.exists(file.path(root, "RNA_clone_analysis"))) {
    stop("Run from CH_publication_release, which contains multiome/ and RNA_clone_analysis/.", call. = FALSE)
  }
  paths <- file.path(root, "RNA_clone_analysis/data/processed",
                     c("clone_multiome_info.RData", "clone_info.RData"))
  if (any(!file.exists(paths))) {
    stop("Missing clone checkpoint(s):\n", paste(paths[!file.exists(paths)], collapse = "\n"),
         call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not available in this R session; no files have been modified.", call. = FALSE)
  }
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  original_figures <- normalizePath(file.path(root, "multiome/figures"), mustWork = TRUE)
  if (identical(output_root, original_figures) ||
      startsWith(output_root, paste0(original_figures, "/"))) {
    stop("Refusing to write into original multiome/figures/.", call. = FALSE)
  }
  out_data <- file.path(output_root, "source_data/FigS9G_clone_abundance_RNA_vs_Multiome.tsv.gz")
  out_pdf <- file.path(output_root, "figures/FigS9G_clone_abundance_REVIEW.pdf")
  if (any(file.exists(c(out_data, out_pdf)))) {
    stop("Existing outputs will not be overwritten:\n",
         paste(c(out_data, out_pdf)[file.exists(c(out_data, out_pdf))], collapse = "\n"),
         call. = FALSE)
  }

  m <- new.env(parent = baseenv())
  r <- new.env(parent = baseenv())
  load(paths[1], envir = m)
  load(paths[2], envir = r)
  cat("=== Checkpoint contents ===\n")
  cat("Multiome: ", paste(ls(m), collapse = ", "), "\n", sep = "")
  cat("RNA: ", paste(ls(r), collapse = ", "), "\n", sep = "")
  if (!all(c("clone.info.multi", "cell.clone.multi") %in% ls(m)) ||
      !"clone.info" %in% ls(r)) {
    stop("Expected clone.info.multi, cell.clone.multi, or clone.info absent. Inspect archive.",
         call. = FALSE)
  }
  multi <- as.data.frame(m$clone.info.multi)
  cellmap <- as.data.frame(m$cell.clone.multi)
  rna <- as.data.frame(r$clone.info)
  cat("clone.info.multi: ", nrow(multi), " rows; columns: ",
      paste(names(multi), collapse = ", "), "\n", sep = "")
  cat("cell.clone.multi: ", nrow(cellmap), " rows; columns: ",
      paste(names(cellmap), collapse = ", "), "\n", sep = "")
  cat("clone.info: ", nrow(rna), " rows; columns: ",
      paste(names(rna), collapse = ", "), "\n", sep = "")
  required <- list(multi = c("larry", "size"),
                   cellmap = c("cloneid", "sample"), rna = "size")
  for (name in names(required)) {
    actual <- switch(name, multi = names(multi), cellmap = names(cellmap), rna = names(rna))
    missing <- setdiff(required[[name]], actual)
    if (length(missing)) stop("Missing ", name, " columns: ", paste(missing, collapse = ", "),
                              call. = FALSE)
  }
  if (is.null(rownames(multi)) || is.null(rownames(rna)) ||
      anyDuplicated(rownames(multi)) || anyDuplicated(rownames(rna))) {
    stop("Clone identifier rownames missing or duplicated; cannot join safely.", call. = FALSE)
  }
  clone_id <- rownames(multi)
  larry_id <- as.character(multi$larry)
  # Historical multiome/scripts/3_mutant_analysis.R uses
  # clone.info[clone.info.multi$larry, ], allowing multiple Multiome clone
  # records to map to the same RNA clone. Do NOT require unique larry IDs.
  valid_larry <- !is.na(larry_id) & nzchar(larry_id)
  multiplicity <- integer(length(larry_id))
  multiplicity[valid_larry] <- as.integer(table(larry_id[valid_larry])[larry_id[valid_larry]])
  duplicates <- valid_larry & multiplicity > 1L
  cat("\n=== Historical LARRY matching ===\n")
  cat("Multiome records with repeated RNA LARRY ID: ", sum(duplicates),
      "; distinct repeated IDs: ", length(unique(larry_id[duplicates])), "\n", sep = "")
  if (any(duplicates)) {
    tab <- data.frame(larry_id = larry_id[duplicates],
                      multiome_clone_id = clone_id[duplicates],
                      multiome_sample = as.character(multi$sample[duplicates]),
                      multiome_size = multi$size[duplicates],
                      stringsAsFactors = FALSE)
    cat("Examples (review to confirm many-to-one identity):\n")
    print(utils::head(tab[order(tab$larry_id), , drop = FALSE], 12), row.names = FALSE)
  }
  to_number <- function(x) suppressWarnings(as.numeric(as.character(x)))
  size_multi <- to_number(multi$size)
  # R's match() explicitly allows repeated query IDs and returns the same
  # RNA clone record for each occurrence, as the historical row lookup does.
  idx <- match(larry_id, rownames(rna))
  size_rna <- to_number(rna$size[idx])
  if (!any(!is.na(idx))) {
    stop("No Multiome LARRY IDs match rownames(clone.info); no files written.", call. = FALSE)
  }
  if ("sample" %in% names(rna)) {
    cat("RNA sample IDs for matched records (first 10): ",
        paste(utils::head(unique(as.character(rna$sample[idx[!is.na(idx)]])), 10), collapse = ", "),
        "\n", sep = "")
  }

  cmid <- as.character(cellmap$cloneid)
  cmsample <- as.character(cellmap$sample)
  valid_cm <- !is.na(cmid) & nzchar(cmid) & cmid %in% clone_id
  mapped <- data.frame(clone_id = cmid[valid_cm], sample = cmsample[valid_cm],
                       stringsAsFactors = FALSE)
  by_id <- split(mapped$sample, mapped$clone_id)
  sample_from_cells <- vapply(clone_id, function(id) {
    vals <- unique(by_id[[id]])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) == 1L) vals else NA_character_
  }, character(1))
  multi_cell_count <- tabulate(match(mapped$clone_id, clone_id), nbins = length(clone_id))
  sample <- sample_from_cells
  if ("sample" %in% names(multi)) {
    provided <- as.character(multi$sample)
    use_provided <- is.na(sample) | !nzchar(sample)
    sample[use_provided] <- provided[use_provided]
    disagree <- !is.na(sample_from_cells) & !is.na(provided) &
      nzchar(provided) & provided != sample_from_cells
    if (any(disagree)) {
      stop("Contradictory sample labels between clone.info.multi and cell.clone.multi: ",
           sum(disagree), " clones.", call. = FALSE)
    }
  }
  sample[is.na(sample) | !sample %in% c("C1", "C2", "D", "T")] <- "NA"

  same_size <- is.finite(size_multi) & size_multi == multi_cell_count
  has_cells <- multi_cell_count > 0
  if (!any(has_cells)) {
    stop("No clone IDs from cell.clone.multi match clone.info.multi; no files written.",
         call. = FALSE)
  }
  cat("\n=== Identity and count checks ===\n")
  cat("LARRY IDs matched to scRNA clones: ", sum(!is.na(idx)), "/", length(idx), "\n", sep = "")
  cat("Distinct matched RNA clones: ", length(unique(larry_id[!is.na(idx)])), "\n", sep = "")
  cat("Multiome clones observed in cell map: ", sum(has_cells), "/", length(clone_id), "\n", sep = "")
  cat("Multiome saved size equals cell-map count: ",
      sum(same_size & has_cells), "/", sum(has_cells), "\n", sep = "")
  cat("Unknown sample labels: ", sum(sample == "NA"), "\n", sep = "")
  if (isTRUE(require_size_consistency) &&
      sum(has_cells) > 0 && sum(same_size & has_cells) / sum(has_cells) < 0.95) {
    cat("Example size mismatches:\n")
    print(utils::head(data.frame(clone_id, recorded = size_multi,
                                 from_cells = multi_cell_count)[has_cells & !same_size, ], 10),
          row.names = FALSE)
    stop("Stored Multiome 'size' does not agree with per-clone cell counts for >=5% of clones. ",
         "No files written; confirm the historical definition before plotting.", call. = FALSE)
  }

  df <- data.frame(
    clone_id = clone_id, larry_id = larry_id,
    larry_mapping_multiplicity = multiplicity,
    sample = sample, multiome_size = size_multi, RNA_size = size_rna,
    multiome_cells_from_mapping = multi_cell_count,
    multiome_size_matches_mapping = same_size,
    stringsAsFactors = FALSE
  )
  good <- is.finite(df$multiome_size) & df$multiome_size > 0 &
    is.finite(df$RNA_size) & df$RNA_size > 0
  df$plotted_on_log_axes <- good
  cat("Clones with positive sizes in both assays: ", sum(good), "\n", sep = "")
  if (sum(good) < min_clones) {
    stop("Too few paired clone abundances to reproduce S9G. Check LARRY-ID mapping; no files written.",
         call. = FALSE)
  }
  cat("Plotted clone counts by sample:\n")
  print(table(df$sample[good]))
  cat("Range RNA: ", paste(range(df$RNA_size[good]), collapse = " - "), "\n", sep = "")
  cat("Range Multiome: ", paste(range(df$multiome_size[good]), collapse = " - "), "\n", sep = "")

  # The artwork shows raw clone counts on log10 axes with a diagonal y=x.
  # No normalization, regression or filtering of clone size is introduced.
  vis <- df[good, , drop = FALSE]
  vis$sample <- factor(vis$sample, levels = c("C1", "C2", "D", "T", "NA"))
  colors <- c(C1 = "#F8766D", C2 = "#7CAE00", D = "#00BFC4",
              T = "#C77CFF", `NA` = "#7F7F7F")
  p <- ggplot2::ggplot(vis, ggplot2::aes(x = multiome_size, y = RNA_size,
                                        colour = sample)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         colour = "black", linewidth = 0.5) +
    ggplot2::geom_point(size = 1.1, alpha = 0.75) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = colors, drop = FALSE) +
    ggplot2::labs(x = "multiome", y = "RNA", colour = NULL,
                  title = "Single clone abundance") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = c(0.85, 0.22),
                   panel.grid.minor = ggplot2::element_blank())

  dir.create(dirname(out_data), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(out_data, open = "wt")
  utils::write.table(df, con, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "NA")
  close(con)
  con_read <- gzfile(out_data, open = "rt")
  readback <- tryCatch(utils::read.delim(con_read, check.names = FALSE,
                                        stringsAsFactors = FALSE),
                       finally = close(con_read))
  if (nrow(readback) != nrow(df) ||
      !isTRUE(all.equal(readback$multiome_size, df$multiome_size,
                        check.attributes = FALSE)) ||
      !isTRUE(all.equal(readback$RNA_size, df$RNA_size,
                        check.attributes = FALSE))) {
    unlink(out_data)
    stop("Source Data round-trip check failed; CSV removed.", call. = FALSE)
  }
  grDevices::pdf(out_pdf, width = 5.1, height = 4.6, useDingbats = FALSE)
  tryCatch(print(p), finally = grDevices::dev.off())
  cat("\nSource Data: ", out_data, "\n", sep = "")
  cat("Figure for historical comparison: ", out_pdf, "\n", sep = "")
  cat("STATUS: REVIEW - must compare with Supplementary Figure 9G artwork, ",
      "and verify clone-size definition against the original clone-processing code.\n", sep = "")
  invisible(list(n_clones = nrow(df), n_plotted = sum(good),
                 n_matched = sum(!is.na(idx)),
                 n_size_consistent = sum(same_size & has_cells),
                 source_data = out_data, figure = out_pdf, status = "REVIEW"))
}
