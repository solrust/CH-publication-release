# CH publication release -- MYC_analysis: Suppl. Fig. 7G-H and 10D
# Historical calculation: MYC_analysis/scripts/MYC_phenotypes.R
# Run from CH_publication_release; inputs and historical figures are read-only.
# S7H/S10D are lightweight; S7G loads RNA_analysis/data/processed/scrna.RData (~7 GiB).

# Load the Matrix version used for CH publication reproduction.
if (!"Matrix" %in% loadedNamespaces()) {
  library(Matrix, lib.loc = Sys.getenv("CH_R_LIB", unset = "/packages"))
}

stopifnot(
  getNamespaceVersion("Matrix") >= package_version("1.6-4")
)

.myc_paths <- function(root, output_root) {
  list(root = root,
       output_root = output_root,
       assoc = file.path(root, "Association_analysis/data/processed/use/association_summary.RData"),
       hsc = file.path(root, "Association_analysis/data/processed/hsc/hsc_axes.RData"),
       rna = file.path(root, "RNA_analysis/data/processed/scrna.RData"),
       motif = file.path(root, "pipeline/resource/motif_associations.xlsx"),
       reference = file.path(root, "MYC_analysis/figures"))
}

myc_preflight <- function(panels = c("S7G", "S7H", "S10D"), root = getwd()) {
  allowed <- c("S7G", "S7H", "S10D")
  stopifnot(length(panels) > 0, !anyDuplicated(panels), all(panels %in% allowed))
  root <- normalizePath(root, mustWork = TRUE)
  p <- .myc_paths(root, file.path(root, "MYC_analysis/publication"))
  required <- c(association = p$assoc,
                if ("S7G" %in% panels) c(HSC_axes = p$hsc, scRNA = p$rna),
                if ("S10D" %in% panels) c(motif_workbook = p$motif))
  tab <- data.frame(input = names(required),
                    present = file.exists(required),
                    GiB = round(ifelse(file.exists(required), file.info(required)$size / 1024^3, NA_real_), 3),
                    path = unname(required), row.names = NULL)
  print(tab, row.names = FALSE)
  pkgs <- c("ggplot2", if ("S10D" %in% panels) "readxl")
  pkgok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  print(data.frame(package = names(pkgok), loadable = unname(pkgok)), row.names = FALSE)
  if (any(!tab$present)) stop("Required input(s) missing; see table above.", call. = FALSE)
  if (any(!pkgok)) stop("Required plotting package(s) unavailable; see table above.", call. = FALSE)
  invisible(tab)
}

.myc_load <- function(path) {
  e <- new.env(parent = baseenv())
  # Loading into an ordinary R environment avoids attaching saved object names globally.
  # Keep saved objects private without masking global analysis variables.
  load(path, envir = e)
  e
}

.myc_assoc <- function(p) {
  e <- .myc_load(p$assoc)
  if (!exists("gsea.clone", e, inherits = FALSE) ||
      !exists("test.basic.list", e, inherits = FALSE))
    stop("Association checkpoint lacks gsea.clone or test.basic.list.", call. = FALSE)
  groups <- e$gsea.clone$gsea.clone
  phenotypes <- names(e$test.basic.list)
  if (!length(groups) || length(groups) != length(phenotypes) ||
      anyNA(phenotypes) || any(!nzchar(phenotypes)))
    stop("Cannot align GSEA results to historical phenotype names.", call. = FALSE)
  list(groups = groups, phenotypes = phenotypes)
}

.myc_gsea <- function(assoc) {
  paths <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2")
  rows <- list()
  for (path in paths) for (i in seq_along(assoc$groups)) {
    tbl <- as.data.frame(assoc$groups[[i]])
    if (!all(c("pathway", "pval", "NES") %in% names(tbl)))
      stop("GSEA table missing pathway, pval or NES.", call. = FALSE)
    idx <- which(grepl(path, as.character(tbl$pathway), fixed = TRUE))
    if (length(idx) != 1L) stop("Expected one ", path, " row for ", assoc$phenotypes[i],
                                "; found ", length(idx), ".", call. = FALSE)
    # Keep the numerical statistics; don't serialize list-valued leadingEdge into TSV.
    cols <- intersect(c("pathway", "pval", "padj", "log2err", "ES", "NES", "size"), names(tbl))
    x <- tbl[idx, cols, drop = FALSE]
    x$phenotype <- assoc$phenotypes[i]
    rows[[length(rows) + 1L]] <- x
  }
  x <- do.call(rbind, rows)
  rownames(x) <- NULL
  x$pval <- as.numeric(x$pval)
  x$NES <- as.numeric(x$NES)
  if (any(!is.finite(x$pval)) || any(x$pval < 0 | x$pval > 1) || any(!is.finite(x$NES)))
    stop("Invalid MYC GSEA statistics.", call. = FALSE)
  x$logP <- pmin(-log10(pmax(x$pval, .Machine$double.xmin)), 10)
  x
}

.myc_motif <- function(p) {
  sheets <- readxl::excel_sheets(p$motif)
  rows <- lapply(sheets, function(sheet) {
    tab <- as.data.frame(readxl::read_excel(p$motif, sheet = sheet))
    required <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
    if (!all(required %in% names(tab))) stop("Missing statistical columns in sheet ", sheet)
    hit <- which(as.character(tab[[1]]) == "MYC")
    if (length(hit) != 1L) stop("Expected one MYC motif in sheet ", sheet,
                                 "; found ", length(hit), call. = FALSE)
    x <- tab[hit, required, drop = FALSE]
    x$phenotype <- sheet
    x$motif <- "MYC motif"
    x
  })
  x <- do.call(rbind, rows)
  rownames(x) <- NULL
  x$plotted <- x$phenotype != "mk bias"
  x$minus_log_p <- -log(pmax(as.numeric(x$pvalue), .Machine$double.xmin))
  if (any(!is.finite(x$log2FoldChange)) || any(!is.finite(x$pvalue)) ||
      any(x$pvalue < 0 | x$pvalue > 1)) stop("Invalid MYC motif statistics.")
  x
}

.myc_hsc <- function(p, assoc) {
  message("Loading HSC axes checkpoint: ", p$hsc)
  hsc <- .myc_load(p$hsc)
  if (!exists("prr", hsc, inherits = FALSE)) stop("HSC axes lack prr.")
  pcs <- hsc$prr$x
  rotation <- hsc$prr$rotation
  if (!is.matrix(pcs) || !is.matrix(rotation) || ncol(pcs) < 40L)
    stop("Unexpected HSC PCA checkpoint structure.")
  cls <- rownames(pcs)
  if (is.null(cls) || anyDuplicated(cls)) stop("HSC PCA cell IDs absent or repeated.")
  message("Loading large scRNA checkpoint (S7G only): ", p$rna)
  rna <- .myc_load(p$rna)
  if (!all(c("emb", "p2counts") %in% ls(rna, all.names = TRUE)))
    stop("scRNA checkpoint lacks emb or p2counts.")
  missing_emb <- setdiff(cls, rownames(rna$emb))
  missing_counts <- setdiff(cls, rownames(rna$p2counts))

  cat("S7G HSC cell alignment:",
      length(cls), "HSC cells;",
      length(missing_emb), "missing from UMAP;",
      length(missing_counts), "missing from counts\n")

  if (length(missing_emb) || length(missing_counts)) {
    stop("S7G cell-ID mismatch. UMAP missing: ",
         paste(head(missing_emb, 5), collapse = ", "),
         "; counts missing: ",
         paste(head(missing_counts, 5), collapse = ", "),
         call. = FALSE)
  }
  if (!"MYC" %in% colnames(rna$p2counts)) stop("No MYC in RNA expression features.")
  # Follow the historical formula: regress MYC expression on 40 HSC PCs.
  y <- as.numeric(rna$p2counts[cls, "MYC"])
  names(y) <- cls
  fit <- stats::lm(y ~ pcs[, seq_len(40L), drop = FALSE])
  score_myc <- as.numeric(stats::fitted(fit))
  # Follow the historical selection: V1 leading edge from first phenotype result.
  first <- as.data.frame(assoc$groups[[1]])
  idx <- which(grepl("HALLMARK_MYC_TARGETS_V1", as.character(first$pathway), fixed = TRUE))
  if (length(idx) != 1L || !"leadingEdge" %in% names(first))
    stop("Cannot recover historical MYC V1 leading-edge genes.")
  targets <- as.character(unlist(first$leadingEdge[idx], use.names = FALSE))
  targets <- intersect(rownames(rotation), targets)
  if (!length(targets)) stop("No MYC leading-edge genes represented in HSC PCA loadings.")
  target_weight <- colMeans(rotation[targets, , drop = FALSE])
  score_targets <- as.numeric(pcs[cls, , drop = FALSE] %*% target_weight)
  coords <- rna$emb[cls, , drop = FALSE]
  if (ncol(coords) < 2L) stop("HSC RNA embedding must have >=2 columns.")
  out <- data.frame(cell_id = cls,
                    UMAP_1 = as.numeric(coords[, 1]),
                    UMAP_2 = as.numeric(coords[, 2]),
                    MYC_expression_score = score_myc,
                    MYC_targets_score = score_targets,
                    stringsAsFactors = FALSE)
  if (any(!is.finite(as.matrix(out[, -1])))) stop("Nonfinite MYC scores or coordinates.")
  attr(out, "n_leading_edge_genes") <- length(targets)
  out
}

.myc_write_table <- function(df, path) {
  con <- gzfile(path, open = "wt")
  tryCatch(utils::write.table(df, con, sep = "\t", quote = FALSE,
                              row.names = FALSE, na = "NA"),
           finally = close(con))
  back <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(identical(dim(back), dim(df)), identical(names(back), names(df)))
  for (nm in names(df)) {
    a <- df[[nm]]; b <- back[[nm]]
    same <- if (is.numeric(a)) isTRUE(all.equal(a, b, tolerance = 1e-10,
                                               check.attributes = FALSE)) else identical(a, b)
    if (!same) stop("Source Data round-trip mismatch: ", nm, call. = FALSE)
  }
  invisible(TRUE)
}

.myc_save_pdf <- function(plot, path, width = 7, height = 7) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  tryCatch(print(plot), finally = grDevices::dev.off())
  if (!file.exists(path) || file.info(path)$size < 1000)
    stop("PDF is missing or empty: ", path, call. = FALSE)
}

.myc_figure <- function(panel, data, path) {
  if (panel == "S7H") {
    gg <- ggplot2::ggplot(data, ggplot2::aes(x = pathway, y = phenotype,
                                            size = logP, fill = NES)) +
      ggplot2::geom_point(shape = 21, color = "black", stroke = .5) +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      ggplot2::scale_size_continuous(range = c(0,20), limits = c(0,10), breaks = c(1,2,5,10)) +
      ggplot2::scale_x_discrete(labels = function(x) sub(".*(HALLMARK_MYC_TARGETS_V[12]).*", "\\1", x)) +
      ggplot2::theme_minimal() + ggplot2::labs(y = "clonal phenotype", x = NULL,
                                              size = "logP", fill = "NES") +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = 20),
                     axis.text.y = ggplot2::element_text(size = 15),
                     axis.text.x = ggplot2::element_text(size = 15, angle = 30,
                                                         vjust = 1, hjust = 1))
    .myc_save_pdf(gg, path)
  } else if (panel == "S10D") {
    d <- data[data$plotted, , drop = FALSE]
    gg <- ggplot2::ggplot(d, ggplot2::aes(y = phenotype, x = motif,
                                         size = minus_log_p, fill = log2FoldChange)) +
      ggplot2::geom_point(shape = 21, color = "black", stroke = .5) +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      ggplot2::scale_size_continuous(range = c(0,20), limits = c(0,10), breaks = c(1,2,5,10)) +
      ggplot2::theme_minimal() + ggplot2::labs(y = "clonal phenotype", x = NULL,
                                              size = "-logP", fill = "log2FoldChange") +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = 20),
                     axis.text.y = ggplot2::element_text(size = 15),
                     axis.text.x = ggplot2::element_text(size = 15, angle = 30,
                                                         vjust = 1, hjust = 1))
    .myc_save_pdf(gg, path)
  } else if (panel == "S7G") {
    grDevices::pdf(path, width = 6, height = 3, useDingbats = FALSE)
    tryCatch({
      oldpar <- par(mfrow = c(1, 2), mar = rep(2,4))
      for (j in seq_along(c("MYC_expression_score", "MYC_targets_score"))) {
        nm <- c("MYC_expression_score", "MYC_targets_score")[j]
        scor <- data[[nm]]
        bins <- as.integer(cut(rank(scor), 10))
        cols <- grDevices::colorRampPalette(c("blue", "grey90", "red"))(10)[bins]
        plot(data$UMAP_1, data$UMAP_2, pch = 19, cex = .5,
             col = grDevices::adjustcolor(cols, alpha.f = .5),
             xaxt = "n", yaxt = "n", xlim = c(-2,5), ylim = c(2,10),
             main = c("MYC", "MYC targets")[j], font.main = 1, cex.main = 2,
             xlab = "", ylab = "")
      }
      par(oldpar)
    }, finally = grDevices::dev.off())
    if (!file.exists(path) || file.info(path)$size < 1000)
      stop("S7G PDF is missing or empty.", call. = FALSE)
  } else stop("Unknown panel.")
}

run_myc_publication <- function(panels = c("S7H", "S10D", "S7G"),
                                root = getwd(), output_root = NULL) {
  root <- normalizePath(root, mustWork = TRUE)
  if (is.null(output_root)) output_root <- file.path(root, "MYC_analysis/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  ref <- normalizePath(file.path(root, "MYC_analysis/figures"), mustWork = TRUE)
  if (identical(output_root, ref) || startsWith(output_root, paste0(ref, "/")))
    stop("Refusing to write to historical MYC_analysis/figures.", call. = FALSE)
  myc_preflight(panels, root)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_root, "source_data"), showWarnings = FALSE)
  dir.create(file.path(output_root, "figures"), showWarnings = FALSE)
  p <- .myc_paths(root, output_root)
  filenames <- list(S7G = c("FigS7G_MYC_HSC_scores.tsv.gz", "FigS7G_MYC_HSC_REVIEW.pdf"),
                    S7H = c("FigS7H_MYC_targets_GSEA.tsv.gz", "FigS7H_MYC_targets_GSEA_REVIEW.pdf"),
                    S10D = c("FigS10D_MYC_motif_associations.tsv.gz", "FigS10D_MYC_motif_REVIEW.pdf"))
  # Prevent replacing approved outputs or partially completed panels.
  for (panel in panels) {
    file_dest <- c(file.path(output_root, "source_data", filenames[[panel]][1]),
                   file.path(output_root, "figures", filenames[[panel]][2]))
    if (any(file.exists(file_dest))) stop("Existing outputs retained for ", panel,
                                         ":\n", paste(file_dest[file.exists(file_dest)], collapse="\n"),
                                         call. = FALSE)
  }
  assoc <- NULL
  results <- list()
  for (panel in panels) {
    cat("\n=== Fig ", panel, " ===\n", sep = "")
    staged <- tempfile(pattern = paste0(".myc_", panel, "_"), tmpdir = output_root)
    dir.create(staged)
    tryCatch({
      if (is.null(assoc) && panel %in% c("S7H", "S7G")) assoc <- .myc_assoc(p)
      data <- switch(panel, S7H = .myc_gsea(assoc),
                     S10D = .myc_motif(p),
                     S7G = .myc_hsc(p, assoc))
      data_file <- file.path(staged, filenames[[panel]][1])
      fig_file <- file.path(staged, filenames[[panel]][2])
      .myc_write_table(data, data_file)
      .myc_figure(panel, data, fig_file)
      dst_data <- file.path(output_root, "source_data", filenames[[panel]][1])
      dst_fig <- file.path(output_root, "figures", filenames[[panel]][2])
      if (!file.rename(data_file, dst_data)) stop("Could not commit MYC Source Data file.")
      if (!file.rename(fig_file, dst_fig)) stop("Could not commit MYC figure PDF.")
      cat("PASS: ", panel, "; ", nrow(data), " numerical rows; 1 PDF.\n", sep="")
      results[[panel]] <- data.frame(panel = panel, n_rows = nrow(data),
                                     source_data = dst_data, figure = dst_fig,
                                     status = "EXPORTED; VISUAL REVIEW PENDING",
                                     stringsAsFactors = FALSE)
    }, finally = unlink(staged, recursive = TRUE))
  }
  invisible(do.call(rbind, results))
}
