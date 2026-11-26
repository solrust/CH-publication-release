# Fig S9F: Corces purified HSPC ATAC versus RNA-annotated Multiome pseudobulks.
# Run from CH_publication_release. Historical data and figures remain unchanged.
# Source: multiome/scripts/7_rev2_corces_comparison.R
# The Corces signature RDS is an archived historical intermediate.
# A compact recount checkpoint is retained so subsequent plotting does not rescan fragments.

.s9f_population_order <- c(
  "HSC", "LMPP-1", "LMPP-2", "MLP", "early GMP", "late GMP",
  "EMPP/MkP", "EryP", "pDC", "Pro-B"
)
.s9f_reference_order <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP")
.s9f_projection_labels <- c(
  `1`="HSC", `2`="LMPP-1", `3`="MLP", `4`="Pro-B", `5`="pDC",
  `6`="LMPP-2", `7`="early GMP", `8`="late GMP", `9`="EMPP/MkP", `10`="EryP"
)

# Values transcribed from the original multiome/figures/corces_atac_benchmark/
# corces_multiome_accessibility_correlation.pdf. These are rounded to two decimals;
# they are a validation reference, NOT the new numerical Source Data.
.s9f_published_rounded <- matrix(c(
   0.65, 0.65, 0.69, 0.35, 0.56, 0.12, 0.23,
   0.49, 0.50, 0.75, 0.30, 0.66, 0.09, 0.29,
   0.54, 0.56, 0.74, 0.35, 0.67, 0.11, 0.27,
   0.45, 0.46, 0.72, 0.28, 0.61, 0.09, 0.29,
   0.34, 0.36, 0.60, 0.25, 0.78, 0.02, 0.19,
   0.24, 0.26, 0.49, 0.22, 0.78,-0.01, 0.13,
   0.25, 0.33, 0.35, 0.56, 0.39, 0.19, 0.09,
   0.16, 0.24, 0.24, 0.41, 0.21, 0.70, 0.08,
   0.23, 0.25, 0.55, 0.13, 0.74,-0.01, 0.23,
   0.35, 0.36, 0.63, 0.16, 0.59, 0.04, 0.53
), nrow=10L, byrow=TRUE,
 dimnames=list(.s9f_population_order, .s9f_reference_order))

inspect_FigS9F_Corces <- function(root=getwd(), obj=NULL) {
  root <- normalizePath(root, mustWork=TRUE)
  sigfile <- file.path(root, "multiome/data/processed/corces_atac_benchmark/corces_hspc_hg38_signatures.rds")
  seuratfile <- file.path(root, "multiome/data/processed/multiome_seurat.rds")
  metafile <- file.path(root, "multiome/data/processed/corces_atac_benchmark/multiome_atac_pseudobulk_metadata.csv")
  cat("=== S9F saved inputs ===\n")
  for (f in c(sigfile, seuratfile, metafile)) {
    cat(if (file.exists(f)) "PRESENT " else "MISSING ",
        if (file.exists(f)) round(file.info(f)$size/1024^2, 1) else NA_real_,
        "MiB ", f, "\n", sep="")
  }
  if (!file.exists(sigfile)) stop("Missing historical Corces signatures RDS.")
  ref <- readRDS(sigfile)
  cat("\nReference object fields:", paste(names(ref), collapse=", "), "\n")
  required <- c("signatures", "ranges", "mean_reference_profiles", "reference_metadata")
  stopifnot(all(required %in% names(ref)))
  cat("Reference union regions:", length(ref$ranges), "\n")
  cat("Mean reference profile dimensions:", paste(dim(ref$mean_reference_profiles), collapse=" x "), "\n")
  cat("Reference populations:", paste(colnames(ref$mean_reference_profiles), collapse=", "), "\n")
  if (file.exists(metafile)) {
    meta <- read.csv(metafile, check.names=FALSE, row.names=1)
    cat("Saved Multiome pseudobulk groups:", nrow(meta), "\n")
    cat("Saved populations:", paste(unique(as.character(meta$population)), collapse=", "), "\n")
  }
  cat("\n=== Packages (availability only; nothing installed) ===\n")
  pkgs <- c("Seurat", "SeuratObject", "Signac", "Matrix", "edgeR", "pheatmap", "gridExtra")
  for (p in pkgs) {
    where <- find.package(p, quiet=TRUE)
    cat(p, if (length(where)) paste0("PRESENT at ", where[1L]) else "MISSING", "\n")
  }
  cat("\nNote: recounting requires accessible, indexed Multiome fragment files.\n")
  if (!is.null(obj)) {
    if (!requireNamespace("Signac", quietly=TRUE)) stop("Signac cannot load; cannot check fragments")
    frags <- Signac::Fragments(obj[["ATAC"]])
    paths <- vapply(frags, function(f) slot(f, "path"), character(1L))
    cat("Fragment files:", length(paths), "\n")
    for (p in paths) cat(if (file.exists(p)) "OK " else "MISSING ", p, "\n")
  } else {
    cat("Fragment paths not inspected: pass an already-loaded multiome obj if available.\n")
  }
  invisible(TRUE)
}

reproduce_FigS9F_Corces <- function(obj=NULL, root=getwd(), output_root=NULL,
                                    overwrite=FALSE, recount_fragments=FALSE) {
  root <- normalizePath(root, mustWork=TRUE)
  sigfile <- file.path(root, "multiome/data/processed/corces_atac_benchmark/corces_hspc_hg38_signatures.rds")
  inputfile <- file.path(root, "multiome/data/processed/multiome_seurat.rds")
  historical_figures <- normalizePath(file.path(root, "multiome/figures"), mustWork=TRUE)
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork=FALSE)
  if (identical(output_root, historical_figures) ||
      startsWith(output_root, paste0(historical_figures, "/"))) {
    stop("Refusing to write inside historical multiome/figures.")
  }
  output_data <- file.path(output_root, "source_data", "FigS9F_Corces_Multiome_Pearson.tsv.gz")
  output_pdf <- file.path(output_root, "figures", "FigS9F_Corces_ATAC_benchmark_REVIEW.pdf")
  pbfile <- file.path(root, "multiome/publication/intermediate/FigS9F_Multiome_Corces_pseudobulk.rds")
  targets <- c(output_data, output_pdf)
  if (!overwrite && any(file.exists(targets))) {
    stop("Existing outputs will not be overwritten:\n",
         paste(targets[file.exists(targets)], collapse="\n"), call.=FALSE)
  }
  if (!file.exists(sigfile)) stop("Missing Corces signatures checkpoint: ", sigfile)
  for (pkg in c("Matrix", "edgeR", "pheatmap", "gridExtra")) {
    if (!requireNamespace(pkg, quietly=TRUE)) stop("Required package unavailable: ", pkg)
  }
  ref <- readRDS(sigfile)
  stopifnot(all(c("ranges", "mean_reference_profiles") %in% names(ref)))
  reference_ids <- names(ref$ranges)
  if (is.null(reference_ids) || anyDuplicated(reference_ids)) {
    stop("Historical Corces reference region IDs are absent or duplicated")
  }
  if (!all(.s9f_reference_order %in% colnames(ref$mean_reference_profiles))) {
    stop("The Corces reference lacks one of the seven published HSPC types")
  }
  reference_profiles <- ref$mean_reference_profiles[reference_ids, .s9f_reference_order, drop=FALSE]
  if (anyNA(reference_profiles)) stop("Corces reference profiles have missing region IDs")

  if (file.exists(pbfile) && !recount_fragments) {
    message("Using previously saved, compact Multiome pseudobulk checkpoint: ", pbfile)
    pb <- readRDS(pbfile)
  } else {
    # First execution must recount fragments. The historical script did not save
    # the Multiome pseudobulk matrix, only its metadata.
    for (pkg in c("Seurat", "SeuratObject", "Signac")) {
      if (!requireNamespace(pkg, quietly=TRUE)) stop("Required for fragment recount: ", pkg)
    }
    if (is.null(obj)) {
      message("Loading historical Multiome Seurat checkpoint: ", inputfile)
      obj <- readRDS(inputfile)
    } else {
      message("Using caller's already-loaded Multiome Seurat object")
    }
    stopifnot(inherits(obj, "Seurat"))
    required_meta <- c("projection", "dataset")
    if (!all(required_meta %in% colnames(obj@meta.data))) {
      stop("Missing required Seurat metadata: ",
           paste(setdiff(required_meta, colnames(obj@meta.data)), collapse=", "))
    }
    fragments <- Signac::Fragments(obj[["ATAC"]])
    if (!length(fragments)) stop("No fragments attached to the saved ATAC assay")
    fragpaths <- vapply(fragments, function(f) slot(f, "path"), character(1L))
    if (any(!file.exists(fragpaths))) {
      # Relocate historical Signac Fragment objects in memory only. UpdatePath
      # checks the stored fragment/index MD5 hashes; mismatches stop the run.
      # Do not change or save the historical Seurat checkpoint.
      samples <- sub(
        ".*/multiome_larry_(C1|C2|D|T)/outs/atac_fragments\\.tsv\\.gz$",
        "\\1", fragpaths
      )
      if (!all(samples %in% c("C1", "C2", "D", "T")) ||
          anyDuplicated(samples) || length(samples) != 4L) {
        stop("Cannot unambiguously map historical fragment paths to four samples:\n",
             paste(fragpaths, collapse="\n"))
      }
      relocated <- file.path(root, "..", "larry_atac",
                             paste0("multiome_larry_", samples),
                             "outs", "atac_fragments.tsv.gz")
      missing <- !file.exists(relocated) |
        !(file.exists(paste0(relocated, ".tbi")) |
          file.exists(paste0(relocated, ".csi")))
      if (any(missing)) stop("Fragment or index missing at IRIS location:\n",
                             paste(relocated[missing], collapse="\n"))
      for (i in seq_along(fragments)) {
        if (!file.exists(fragpaths[i])) {
          message("Relocating fragment for ", samples[i], ": ", relocated[i])
          fragments[[i]] <- Signac::UpdatePath(
            object = fragments[[i]], new.path = relocated[i]
          )
        }
      }
      fragpaths <- vapply(fragments, function(f) slot(f, "path"), character(1L))
    }
    if (any(!file.exists(fragpaths))) {
      stop("Fragment path(s) still missing:\n",
           paste(fragpaths[!file.exists(fragpaths)], collapse="\n"))
    }
    message("Recounting Multiome fragments over ", length(ref$ranges),
            " historical reference regions. This is the potentially expensive step.")
    fm <- Signac::FeatureMatrix(fragments=fragments, features=ref$ranges,
                                 cells=colnames(obj), process_n=2000, verbose=TRUE)
    rownames(fm) <- reference_ids
    meta <- obj@meta.data[colnames(fm), c("dataset", "projection"), drop=FALSE]
    pop <- unname(.s9f_projection_labels[as.character(meta$projection)])
    missing_map <- is.na(pop)
    pop[missing_map] <- as.character(meta$projection[missing_map])
    meta$population <- pop
    keep <- !is.na(meta$population) & meta$population %in% .s9f_population_order &
            !is.na(meta$dataset)
    fm <- fm[, keep, drop=FALSE]
    meta <- meta[keep, , drop=FALSE]
    groups <- interaction(meta$dataset, meta$population, sep="__", drop=TRUE)
    sizes <- table(groups)
    keep <- as.character(groups) %in% names(sizes)[sizes >= 30L]
    fm <- fm[, keep, drop=FALSE]
    meta <- meta[keep, , drop=FALSE]
    groups <- droplevels(groups[keep])
    if (!ncol(fm) || !nlevels(groups)) stop("No pseudobulks with at least 30 cells")
    design <- Matrix::sparseMatrix(
      i=seq_along(groups), j=as.integer(groups), x=rep(1, length(groups)),
      dims=c(length(groups), nlevels(groups)),
      dimnames=list(NULL, levels(groups)))
    pbcounts <- fm %*% design
    first <- match(levels(groups), as.character(groups))
    pbmeta <- meta[first, c("dataset", "population"), drop=FALSE]
    rownames(pbmeta) <- levels(groups)
    pbmeta$n_cells <- as.integer(table(groups)[rownames(pbmeta)])
    pb <- list(counts=pbcounts, metadata=pbmeta)
    dir.create(dirname(pbfile), recursive=TRUE, showWarnings=FALSE)
    if (!file.exists(pbfile)) {
      saveRDS(pb, pbfile, compress=TRUE)
      message("Saved compact pseudobulk checkpoint: ", pbfile)
    } else if (recount_fragments) {
      message("Existing pseudobulk checkpoint preserved: ", pbfile,
              "; this run uses the new in-memory recount.")
    }
    rm(fm, design, pbcounts)
  }
  stopifnot(all(c("counts", "metadata") %in% names(pb)))
  counts <- pb$counts
  meta <- pb$metadata
  if (!identical(colnames(counts), rownames(meta))) stop("Pseudobulk metadata is misaligned")
  if (!all(reference_ids %in% rownames(counts))) stop("Pseudobulk counts lack Corces reference regions")
  if (any(Matrix::colSums(counts) <= 0)) stop("At least one pseudobulk has no ATAC fragments")
  counts <- counts[reference_ids, , drop=FALSE]
  present <- .s9f_population_order %in% as.character(meta$population)
  if (!all(present)) stop("Missing published Multiome population(s): ",
                         paste(.s9f_population_order[!present], collapse=", "))

  dge <- edgeR::DGEList(counts=as.matrix(counts))
  dge <- edgeR::calcNormFactors(dge, method="TMM")
  logcpm <- edgeR::cpm(dge, log=TRUE, prior.count=2)
  population_profiles <- vapply(.s9f_population_order, function(pop) {
    rowMeans(logcpm[, as.character(meta$population)==pop, drop=FALSE])
  }, numeric(nrow(logcpm)))
  rownames(population_profiles) <- rownames(logcpm)
  shared <- intersect(rownames(reference_profiles), rownames(population_profiles))
  shared <- shared[Matrix::rowSums(counts[shared, , drop=FALSE]) > 0]
  if (length(shared) < 100) stop("Too few overlapping accessible reference regions")

  raw_cor <- stats::cor(x=reference_profiles[shared, , drop=FALSE],
                        y=population_profiles[shared, , drop=FALSE],
                        method="pearson", use="pairwise.complete.obs")
  raw_cor <- t(raw_cor)  # rows: Multiome; columns: Corces immunophenotypes
  if (!identical(dimnames(raw_cor), dimnames(.s9f_published_rounded))) {
    stop("Unexpected correlation matrix ordering or labels")
  }
  if (any(!is.finite(raw_cor))) stop("Nonfinite correlations; investigate pseudobulk inputs")
  normalized <- apply(raw_cor, 2, function(v) {
    rng <- range(v)
    if (diff(rng)==0) rep(0.5, length(v)) else (v-rng[1])/diff(rng)
  })
  dimnames(normalized) <- dimnames(raw_cor)
  top_row <- apply(raw_cor, 2, which.max)
  display <- matrix(sprintf("%.2f", raw_cor), nrow=nrow(raw_cor),
                    dimnames=dimnames(raw_cor))
  for (j in seq_len(ncol(raw_cor))) display[top_row[j], j] <- paste0(display[top_row[j], j], " *")
  expected <- .s9f_published_rounded
  difference <- abs(raw_cor-expected)
  matched <- difference <= 0.00501
  cat("\n=== S9F numerical comparison with historical PDF (2-decimal values) ===\n")
  cat("Shared accessible regions:", length(shared), "(historical: 34973)\n")
  cat("Published values matched within rounding:", sum(matched), "/", length(matched), "\n")
  cat("Maximum absolute difference from displayed historical value:",
      signif(max(difference),4), "\n")
  result_status <- if (length(shared)==34973L && all(matched)) "PASS" else "REVIEW"
  cat("First-pass status:", result_status, "\n")

  long <- expand.grid(multiome_population=rownames(raw_cor),
                      reference_population=colnames(raw_cor),
                      KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  idx <- cbind(match(long$multiome_population, rownames(raw_cor)),
               match(long$reference_population, colnames(raw_cor)))
  long$pearson_r <- raw_cor[idx]
  long$color_normalized_within_reference <- normalized[idx]
  long$best_match_in_reference <- idx[,1L] == top_row[idx[,2L]]
  long$shared_reference_regions <- length(shared)

  dir.create(dirname(output_data), recursive=TRUE, showWarnings=FALSE)
  dir.create(dirname(output_pdf), recursive=TRUE, showWarnings=FALSE)
  tmp_data <- paste0(output_data, ".tmp")
  tmp_pdf <- tempfile(pattern = ".FigS9F_plot_",
                      tmpdir = dirname(output_pdf), fileext = ".pdf")
  on.exit(unlink(c(tmp_data, tmp_pdf)), add=TRUE)
  con <- gzfile(tmp_data, "wt")
  tryCatch(utils::write.table(long, con, sep="\t", row.names=FALSE,
                              col.names=TRUE, quote=FALSE, na="NA"),
           finally=close(con))

  grDevices::pdf(tmp_pdf, width=5, height=5)
  pdf_device <- grDevices::dev.cur()
  tryCatch({
    ph <- pheatmap::pheatmap(normalized, cluster_rows=FALSE, cluster_cols=FALSE,
            legend_breaks=c(0,0.5,1), legend_labels=c("Low","Mid","High"),
            display_numbers=display, number_color="grey25", border_color=NA,
            angle_col=90,
            main=paste0("Multiome versus purified HSPC ATAC\n",
                        length(shared), " reference regions\n"), silent=TRUE)
    combined <- gridExtra::arrangeGrob(ph$gtable,
      left=grid::textGrob("Multiome populations", rot=90,
                          gp=grid::gpar(fontsize=16, fontface="bold")),
      bottom=grid::textGrob("HSPC immunophenotypes",
                            gp=grid::gpar(fontsize=16, fontface="bold")))
    grid::grid.newpage()
    grid::grid.draw(combined)
  }, finally={
    if (pdf_device %in% grDevices::dev.list())
      grDevices::dev.off(which = pdf_device)
  })
  if (!file.exists(tmp_pdf) || file.info(tmp_pdf)$size < 1000) {
    stop("Plot PDF is missing or empty; nothing replaced")
  }
  if (overwrite) unlink(targets[file.exists(targets)])
  if (!file.rename(tmp_data, output_data)) stop("Failed to finalize Source Data")
  if (!file.rename(tmp_pdf, output_pdf)) stop("Failed to finalize figure")
  message("Source Data: ", output_data)
  message("Figure: ", output_pdf)
  invisible(list(status=result_status, n_reference_regions=length(shared),
                 n_rounded_values_matched=sum(matched),
                 source_data=output_data, figure=output_pdf, pseudobulk_checkpoint=pbfile))
}
