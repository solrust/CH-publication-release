# Multiome: canonical publication entry point (fresh, non-resumable run).
# Working directory: CH_publication_release (repository root).
#
# SOURCE DATA and PDF generation from existing processed scientific inputs.
# No historical data/figure files or approved first-pass outputs are overwritten.
# Every run requires an EMPTY NEW output directory. Deliberately no resume flag.
#
# Example (after restarting RStudio in the historical /packages environment):
#   source("multiome/publication/run_publication.R")
#   result <- run_multiome_publication(mode = "checkpoint")
#
# mode="checkpoint": regenerate 12 of 15 tables from upstream checkpoints,
#   workbooks, or metadata; copy 3 validated numerical tables (Fig6F, Fig6H,
#   FigS9A) as documented numerical checkpoints rather than repeating their
#   expensive upstream calculations. Regenerate 17 PDFs from those tables.
# mode="full": additionally independently re-estimate Fig6F glmmTMB models
#   from saved peak modules, recount Fig6H from fragments, and recompute S9A
#   detection fractions from the integrated RNA object. Still reuses archived
#   S9E motif results, S9F pseudobulk counts and Fig6J split predictions.
#
# Historical producer code and original figures remain under multiome/scripts/
# and multiome/figures/. Figure-specific publication functions remain in this
# folder for scientific provenance; do not run them manually for production.

run_multiome_publication <- function(
    root = getwd(),
    output_root = NULL,
    mode = c("checkpoint", "full"),
    tolerance = 1e-8) {
  mode <- match.arg(mode)
  root <- normalizePath(root, mustWork = TRUE)
  pub <- file.path(root, "multiome", "publication")
  if (!dir.exists(pub) || !dir.exists(file.path(root, "RNA_clone_analysis")))
    stop("Working directory must be CH_publication_release", call. = FALSE)
  if (is.null(output_root)) output_root <- file.path(pub, "production_clean")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  # The production destination must be a new child of publication/, not pub itself.
  # Since the historical multiome/{figures,data} are siblings, they cannot overlap.
  if (!startsWith(output_root, paste0(pub, "/")) ||
      identical(output_root, pub))
    stop("output_root must be a NEW subdirectory inside multiome/publication.", call. = FALSE)
  if (dir.exists(output_root) &&
      length(list.files(output_root, all.files = TRUE, no.. = TRUE)))
    stop("Refusing nonempty output directory: ", output_root, "\n",
         "A clean run must start in a NEW empty output directory. ",
         "Do not delete an existing production run until you have inspected it.",
         call. = FALSE)
  for (dir in c("source_data", "figures", "intermediate", "validation"))
    dir.create(file.path(output_root, dir), showWarnings = FALSE, recursive = TRUE)

  # A named private namespace avoids polluting the interactive R session.
  env <- new.env(parent = globalenv())
  load_script <- function(filename) {
    f <- file.path(pub, filename)
    if (!file.exists(f)) stop("Missing publication script: ", f, call. = FALSE)
    sys.source(f, envir = env)
    invisible(NULL)
  }
  if ("Matrix" %in% loadedNamespaces()) {
    if (getNamespaceVersion("Matrix") < package_version("1.6-4") ||
        normalizePath(getNamespaceInfo("Matrix", "path")) !=
          normalizePath(file.path(Sys.getenv("CH_R_LIB", unset = "/packages"), "Matrix")))
      stop("Restart R and load Matrix 1.6-4 from /packages before this workflow.")
  } else {
    if (!dir.exists(file.path(Sys.getenv("CH_R_LIB", unset = "/packages"), "Matrix"))) stop("Missing Matrix in CH_R_LIB (default /packages)")
    suppressPackageStartupMessages(library(Matrix, lib.loc = Sys.getenv("CH_R_LIB", unset = "/packages")))
  }
  obj_cache <- NULL
  get_obj <- function() {
    if (is.null(obj_cache)) {
      if (!requireNamespace("Seurat", quietly = TRUE) ||
          !requireNamespace("Signac", quietly = TRUE))
        stop("Seurat/Signac cannot load in the current R environment")
      message("Loading Multiome Seurat checkpoint once")
      obj_cache <<- readRDS(file.path(root,
                         "multiome/data/processed/multiome_seurat.rds"))
    }
    obj_cache
  }
  copy_checkpoint <- function(filename) {
    from <- file.path(pub, "intermediate", filename)
    to <- file.path(output_root, "intermediate", filename)
    if (!file.exists(from)) stop("Missing intermediate checkpoint: ", from)
    if (!file.copy(from, to, overwrite = FALSE))
      stop("Cannot stage checkpoint: ", filename)
  }
  # Provenance is explicit for each saved table; copied numerical checkpoints
  # are distinguished from independent re-exports in the final report.
  source_provenance <- c(
    "Fig6B_FigS9B_WNN_UMAP.tsv.gz" = "SEURAT_CHECKPOINT",
    "FigS9C_RNA_UMAP.tsv.gz" = "SEURAT_CHECKPOINT",
    "FigS9C_ATAC_UMAP.tsv.gz" = "SEURAT_CHECKPOINT",
    "Fig6C_SPIB_RNA_chromVAR.tsv.gz" = "SEURAT_CHECKPOINT",
    "Fig6D_clone_umap_cells.tsv.gz" = "CLONE_METADATA_PLUS_WNN_TABLE",
    "Fig6F_lineage_module_associations.tsv.gz" = if (mode == "checkpoint") "CACHED_SOURCE_DATA_NOT_REFIT" else "MODULE_CHECKPOINT_GLMMTMB_REFIT",
    "Fig6G_motif_clone_phenotype_associations.tsv.gz" = "HISTORICAL_WORKBOOK",
    "Fig6H_HOXA_genotype_coverage.tsv.gz" = if (mode == "checkpoint") "CACHED_SOURCE_DATA_NOT_REFIT" else "FRAGMENT_RECOUNT",
    "Fig6I_motif_genotype_associations.tsv.gz" = "HISTORICAL_WORKBOOK_PLUS_SEURAT_GC",
    "Fig6J_heritability_splits.tsv.gz" = "SAVED_7200_SPLIT_RESULTS",
    "Fig6J_heritability_summary.tsv" = "SAVED_7200_SPLIT_RESULTS",
    "FigS9A_gene_detection_all_genes.tsv.gz" = if (mode == "checkpoint") "CACHED_SOURCE_DATA_NOT_REFIT" else "RNA_INTEGRATION_COUNTS_RECOMPUTE",
    "FigS9E_cluster_motif_enrichment.tsv.gz" = "SAVED_MOTIF_ENRICHMENT_RESULTS",
    "FigS9F_Corces_Multiome_Pearson.tsv.gz" = "SAVED_PSEUDOBULK_COUNTS",
    "FigS9G_clone_abundance_RNA_vs_Multiome.tsv.gz" = "RNA_MULTIOME_CLONE_CHECKPOINTS"
  )
  record <- function(filename, origin) {
    source_provenance[filename] <<- origin
    invisible(NULL)
  }
  make_stage <- function(label, expr, outputs = character()) {
    cat("\n=== ", label, " ===\n", sep = "")
    # This function never skips. Output_root was checked empty at entry.
    force(expr)
    # For each stage, ensure that it produced every promised output.
    if (length(outputs)) {
      paths <- file.path(output_root, outputs)
      good <- file.exists(paths) & !is.na(file.info(paths)$size) &
        file.info(paths)$size >= 100L
      if (!all(good)) stop("Incomplete stage ", label, ": ",
        paste(paths[!good], collapse = ", "), call. = FALSE)
    }
    invisible(NULL)
  }
  # 1) Global RNA/ATAC/WNN coordinates and existing annotation.
  make_stage("Fig6B + FigS9B/C embeddings", {
    load_script("export_multiome_batch1_embeddings.R")
    env$export_multiome_batch1_embeddings(obj = get_obj(), root = root,
                                             output_root = output_root)
    for (nm in c("Fig6B_FigS9B_WNN_UMAP.tsv.gz", "FigS9C_RNA_UMAP.tsv.gz",
                 "FigS9C_ATAC_UMAP.tsv.gz")) record(nm, "SEURAT_CHECKPOINT")
    load_script("plot_multiome_batch1_embeddings.R")
    env$plot_multiome_batch1_embeddings(root = root,
        source_dir = file.path(output_root, "source_data"), output_root = output_root)
  }, outputs = c("source_data/Fig6B_FigS9B_WNN_UMAP.tsv.gz",
                 "source_data/FigS9C_RNA_UMAP.tsv.gz",
                 "source_data/FigS9C_ATAC_UMAP.tsv.gz",
                 "figures/Fig6B_reproduced_REVIEW.pdf",
                 "figures/FigS9B_reproduced_REVIEW.pdf",
                 "figures/FigS9C_reproduced_REVIEW.pdf",
                 "figures/Batch1_embeddings_visual_review.pdf",
                 "validation/Batch1_embeddings_plot_summary.tsv"))
  # 2) SPIB expression and chromVAR scores.
  make_stage("Fig6C", {
    load_script("export_plot_Fig6C_SPIB.R")
    env$reproduce_Fig6C_SPIB(obj = get_obj(), root = root,
                            output_root = output_root)
    record("Fig6C_SPIB_RNA_chromVAR.tsv.gz", "SEURAT_CHECKPOINT")
  }, outputs = c("source_data/Fig6C_SPIB_RNA_chromVAR.tsv.gz",
                 "figures/Fig6C_SPIB_reproduced_REVIEW.pdf"))
  # 3) Clone overlays rejoin metadata to the existing approved WNN data.
  make_stage("Fig6D", {
    load_script("export_plot_Fig6D.R")
    env$reproduce_Fig6D(root = root, output_root = output_root,
      wnn_source_file = file.path(output_root, "source_data",
        "Fig6B_FigS9B_WNN_UMAP.tsv.gz"))
    record("Fig6D_clone_umap_cells.tsv.gz", "CLONE_METADATA_PLUS_WNN_TABLE")
  }, outputs = c("source_data/Fig6D_clone_umap_cells.tsv.gz",
                 "figures/Fig6D_reproduced_REVIEW.pdf"))
  # 4) Motif statistics from original exported historical workbooks.
  make_stage("Fig6G + Fig6I", {
    load_script("export_plot_Fig6G_from_workbook.R")
    env$reproduce_Fig6G_from_workbook(root = root, output_root = output_root)
    record("Fig6G_motif_clone_phenotype_associations.tsv.gz", "HISTORICAL_WORKBOOK")
    load_script("export_plot_Fig6I_from_workbook.R")
    env$reproduce_Fig6I_from_workbook(obj = get_obj(), root = root,
                                      output_root = output_root)
    record("Fig6I_motif_genotype_associations.tsv.gz", "HISTORICAL_WORKBOOK_PLUS_SEURAT_GC")
  }, outputs = c("source_data/Fig6G_motif_clone_phenotype_associations.tsv.gz",
                 "figures/Fig6G_motif_phenotypes_REVIEW.pdf",
                 "source_data/Fig6I_motif_genotype_associations.tsv.gz",
                 "figures/Fig6I_motif_genotypes_REVIEW.pdf"))
  # 5) Matched scRNA/ATAC clone sizes.
  make_stage("FigS9G", {
    load_script("export_plot_FigS9G_clone_abundance.R")
    env$reproduce_FigS9G_clone_abundance(root = root, output_root = output_root)
    record("FigS9G_clone_abundance_RNA_vs_Multiome.tsv.gz", "RNA_MULTIOME_CLONE_CHECKPOINTS")
  }, outputs = c("source_data/FigS9G_clone_abundance_RNA_vs_Multiome.tsv.gz",
                 "figures/FigS9G_clone_abundance_REVIEW.pdf"))
  # 6) Corces ATAC validation from saved pseudobulk counts; no fragment recount.
  make_stage("FigS9F", {
    load_script("export_plot_FigS9F_Corces.R")
    env$reproduce_FigS9F_Corces(root = root, output_root = output_root,
                                 recount_fragments = FALSE)
    record("FigS9F_Corces_Multiome_Pearson.tsv.gz", "SAVED_PSEUDOBULK_COUNTS")
  }, outputs = c("source_data/FigS9F_Corces_Multiome_Pearson.tsv.gz",
                 "figures/FigS9F_Corces_ATAC_benchmark_REVIEW.pdf"))
  # 7) Heritability outputs from the archived 7200 per-split Spearman values.
  make_stage("Fig6J", {
    load_script("reproduce_Fig6J_heritability.R")
    if (!exists("Fig6J_finalize_saved", envir = env, inherits = FALSE))
      stop("Fig6J_finalize_saved() missing from the canonical Fig6J script")
    copy_checkpoint("Fig6J_split_results_checkpoint.rds")
    env$Fig6J_finalize_saved(root = root, output_root = output_root)
    record("Fig6J_heritability_splits.tsv.gz", "SAVED_7200_SPLIT_RESULTS")
    record("Fig6J_heritability_summary.tsv", "SAVED_7200_SPLIT_RESULTS")
  }, outputs = c("source_data/Fig6J_heritability_splits.tsv.gz",
                 "source_data/Fig6J_heritability_summary.tsv",
                 "figures/Fig6J_WT_vs_TET2_REVIEW.pdf",
                 "figures/Fig6J_WT_vs_DNMT3A_REVIEW.pdf"))
  # 8) Re-express archived Signac motif-enrichment results numerically.
  make_stage("FigS9E", {
    cp <- readRDS(file.path(pub, "intermediate/FigS9E_motif_enrichment_checkpoint.rds"))
    if (!identical(cp$clusters, as.character(1:10)) || length(cp$results) != 10L)
      stop("FigS9E motif checkpoint incomplete")
    tabs <- lapply(cp$results, function(x) {
      y <- x$enriched_motifs
      if (is.null(y) || !nrow(y)) return(NULL)
      y <- as.data.frame(y)
      if (!"cluster" %in% names(y)) y$cluster <- x$cluster
      if (!"motif_id" %in% names(y)) y$motif_id <- rownames(y)
      y$selected_peaks <- length(x$significant_peak_ids)
      y$n_tested_peaks <- x$n_tested_peaks
      y$random_seed <- x$seed
      y
    })
    tabs <- Filter(Negate(is.null), tabs)
    fields <- unique(unlist(lapply(tabs, names)))
    tabs <- lapply(tabs, function(y) {
      for (field in setdiff(fields, names(y))) y[[field]] <- NA
      y[, fields, drop = FALSE]
    })
    full <- do.call(rbind, tabs)
    rownames(full) <- NULL
    destination <- file.path(output_root, "source_data/FigS9E_cluster_motif_enrichment.tsv.gz")
    con <- gzfile(destination, "wt")
    tryCatch(utils::write.table(full, con, sep = "\t", row.names = FALSE,
                                col.names = TRUE, quote = FALSE, na = "NA"),
             finally = close(con))
    record("FigS9E_cluster_motif_enrichment.tsv.gz", "SAVED_MOTIF_ENRICHMENT_RESULTS")
  }, outputs = c("source_data/FigS9E_cluster_motif_enrichment.tsv.gz"))
  # 9) The three large first-pass calculations: reuse approved numerical Source
  # Data in routine mode, or explicitly repeat their upstream calculations.
  cached <- c("Fig6F_lineage_module_associations.tsv.gz",
              "Fig6H_HOXA_genotype_coverage.tsv.gz",
              "FigS9A_gene_detection_all_genes.tsv.gz")
  if (mode == "checkpoint") {
    make_stage("Fig6F/6H/S9A: approved numerical checkpoints", {
      for (name in cached) {
        from <- file.path(pub, "source_data", name)
        to <- file.path(output_root, "source_data", name)
        if (!file.exists(from) || !file.copy(from, to, overwrite = FALSE))
          stop("Required cached Source Data missing or cannot be staged: ", name)
        record(name, "CACHED_SOURCE_DATA_NOT_REFIT")
      }
    }, outputs = paste0("source_data/", cached))
  } else {
    make_stage("Fig6F: glmmTMB refit from saved lineage modules", {
      load_script("reproduce_Fig6F_lineage_modules.R")
      copy_checkpoint("Fig6F_lineage_modules_checkpoint.rds")
      env$Fig6F_fit_plot(obj = get_obj(), root = root, output_root = output_root)
      record(cached[1L], "MODULE_CHECKPOINT_GLMMTMB_REFIT")
    }, outputs = c("source_data/Fig6F_lineage_module_associations.tsv.gz",
                   "figures/Fig6F_lineage_modules_REVIEW.pdf"))
    make_stage("Fig6H: Signac HOXA coverage from fragments", {
      load_script("export_plot_Fig6H_HOXA.R")
      env$reproduce_Fig6H_HOXA(obj = get_obj(), root = root, output_root = output_root)
      record(cached[2L], "FRAGMENT_RECOUNT")
    }, outputs = c("source_data/Fig6H_HOXA_genotype_coverage.tsv.gz",
                   "figures/Fig6H_HOXA_genotype_REVIEW.pdf"))
    obj_cache <- NULL
    invisible(gc(verbose = FALSE))
    make_stage("FigS9A: all-gene detection from RNA integration", {
      load_script("export_plot_FigS9A_all_genes.R")
      env$reproduce_FigS9A_all_genes(root = root, output_root = output_root)
      record(cached[3L], "RNA_INTEGRATION_COUNTS_RECOMPUTE")
    }, outputs = c("source_data/FigS9A_gene_detection_all_genes.tsv.gz",
                   "figures/FigS9A_gene_detection_all_genes_REVIEW.pdf"))
  }
  obj_cache <- NULL
  invisible(gc(verbose = FALSE))

  # Check all 15 panel-level tables against the approved first pass.
  # The check is separate from numerical provenance; copied tables are marked.
  read_table <- function(filename) {
    con <- if (grepl("\\.gz$", filename)) gzfile(filename, "rt") else file(filename, "rt")
    on.exit(close(con), add = TRUE)
    utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
                      quote = "", comment.char = "", na.strings = "NA")
  }
  compare <- function(a, b) {
    if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) return(FALSE)
    for (i in seq_along(a)) {
      if (is.numeric(a[[i]]) && is.numeric(b[[i]])) {
        if (!isTRUE(all.equal(a[[i]], b[[i]], tolerance = tolerance,
                              check.attributes = FALSE))) return(FALSE)
      } else if (!identical(as.character(a[[i]]), as.character(b[[i]]))) return(FALSE)
    }
    TRUE
  }
  expected <- c(
    "Fig6B_FigS9B_WNN_UMAP.tsv.gz", "FigS9C_RNA_UMAP.tsv.gz",
    "FigS9C_ATAC_UMAP.tsv.gz", "Fig6C_SPIB_RNA_chromVAR.tsv.gz",
    "Fig6D_clone_umap_cells.tsv.gz", "Fig6F_lineage_module_associations.tsv.gz",
    "Fig6G_motif_clone_phenotype_associations.tsv.gz",
    "Fig6H_HOXA_genotype_coverage.tsv.gz",
    "Fig6I_motif_genotype_associations.tsv.gz",
    "Fig6J_heritability_splits.tsv.gz", "Fig6J_heritability_summary.tsv",
    "FigS9A_gene_detection_all_genes.tsv.gz",
    "FigS9E_cluster_motif_enrichment.tsv.gz",
    "FigS9F_Corces_Multiome_Pearson.tsv.gz",
    "FigS9G_clone_abundance_RNA_vs_Multiome.tsv.gz")
  stopifnot(setequal(expected, names(source_provenance)))
  data_report <- do.call(rbind, lapply(expected, function(nm) {
    target <- file.path(output_root, "source_data", nm)
    baseline <- file.path(pub, "source_data", nm)
    if (!file.exists(target) || !file.exists(baseline))
      stop("Missing generated or approved Source Data: ", nm)
    x <- read_table(target); y <- read_table(baseline)
    if (!nrow(x)) stop("Empty generated Source Data table: ", nm)
    data.frame(file = nm, n_rows = nrow(x),
               provenance = unname(source_provenance[nm]),
               numeric_comparison = if (compare(x, y)) "PASS" else "REVIEW",
               stringsAsFactors = FALSE)
  }))
  if (any(data_report$numeric_comparison != "PASS"))
    stop("A regenerated Source Data table differs from approved results; inspect the output folder.")

  # Render PDFs FROM THE NEWLY WRITTEN TABLES, not from the approved tables.
  # PDFs created by an earlier stage of THIS fresh run are retained; the
  # renderer creates any remaining publication panels from the new tables.
  cat("\n=== Render PDFs from regenerated Source Data ===\n")
  pdf_report <- plot_multiome_publication_PDFs(root = root,
                  source_dir = file.path(output_root, "source_data"),
                  output_dir = file.path(output_root, "figures"), overwrite = FALSE)
  expected_pdf <- unique(pdf_report$pdf)
  if (length(expected_pdf) != 17L) stop("Expected 17 PDF components; got ", length(expected_pdf))
  for (nm in expected_pdf) {
    path <- file.path(output_root, "figures", nm)
    if (!file.exists(path) || file.info(path)$size < 1000L)
      stop("Invalid generated PDF: ", nm)
  }
  # One production report records provenance for BOTH numerical tables and PDFs.
  report_dir <- file.path(output_root, "validation")
  final_report <- rbind(
    data.frame(type = "source_data", file = data_report$file,
               n_rows = data_report$n_rows, provenance = data_report$provenance,
               status = data_report$numeric_comparison, stringsAsFactors = FALSE),
    data.frame(type = "pdf", file = pdf_report$pdf, n_rows = NA_integer_,
               provenance = pdf_report$note,
               status = pdf_report$status, stringsAsFactors = FALSE)
  )
  utils::write.table(final_report, file.path(report_dir, "production_report.tsv"),
       quote = FALSE, sep = "\t", row.names = FALSE)
  cat("\n=== Reproduction complete ===\n")
  cat(nrow(data_report), "numeric Source Data tables;", length(expected_pdf), "PDF components.\n")
  cat("Level:", mode, "; approved data and original figures unchanged.\n")
  cat("PROVENANCE NOTE: checkpoint mode copies 3 approved numerical tables ",
      "(Fig6F, Fig6H, FigS9A), never silently treats them as recomputed.\n", sep="")
  cat("Output:", output_root, "\n")
  cat("Excluded by agreement: Fig S9D; Fig S10B/C/D/E/G (delegated/deferred).\n")
  invisible(list(output_dir = output_root, source_data = data_report, pdf = pdf_report))
}


# Internal PDF renderer: draws 17 components from the Source Data just produced.
plot_multiome_publication_PDFs <- function(
    root = getwd(),
    output_dir = NULL,
    source_dir = NULL,
    overwrite = FALSE) {
  root <- normalizePath(root, mustWork = TRUE)
  pub <- file.path(root, "multiome", "publication")
  original <- file.path(root, "multiome", "figures")
  approved <- file.path(pub, "figures")
  src <- if (is.null(source_dir)) file.path(pub, "source_data") else normalizePath(source_dir, mustWork = TRUE)
  if (!dir.exists(src) || !dir.exists(original)) {
    stop("Run from CH_publication_release with existing Multiome Source Data.", call. = FALSE)
  }
  if (is.null(output_dir)) output_dir <- file.path(approved, "production")
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  old_dir <- normalizePath(original, mustWork = TRUE)
  ref_dir <- normalizePath(approved, mustWork = TRUE)
  if (identical(output_dir, old_dir) || startsWith(output_dir, paste0(old_dir, "/")) ||
      identical(output_dir, ref_dir)) {
    stop("Refusing to overwrite historical or approved first-pass figure folders.", call. = FALSE)
  }
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE))
    stop("Cannot create output folder: ", output_dir)
  for (p in c("ggplot2", "ggrepel", "gridExtra", "scales", "pheatmap")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing plotting package: ", p)
  }
  read_data <- function(filename, required) {
    f <- file.path(src, filename)
    if (!file.exists(f)) stop("Missing approved Source Data: ", f, call. = FALSE)
    con <- if (grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "rt")
    on.exit(close(con), add = TRUE)
    x <- utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
                           na.strings = "NA", quote = "", comment.char = "")
    missing <- setdiff(required, names(x))
    if (length(missing)) stop(filename, " is missing: ", paste(missing, collapse = ", "))
    if (!nrow(x)) stop(filename, " has no rows")
    x
  }
  written <- list()
  check_pdf <- function(f) {
    if (!file.exists(f) || is.na(file.info(f)$size) || file.info(f)$size < 1000)
      return(FALSE)
    con <- file(f, "rb")
    on.exit(close(con), add = TRUE)
    identical(rawToChar(readBin(con, "raw", n = 5L)), "%PDF-")
  }
  save_pdf <- function(name, draw, width, height, note = "FROM_SOURCE_DATA_TABLES") {
    target <- file.path(output_dir, name)
    if (file.exists(target) && !isTRUE(overwrite)) {
      if (!check_pdf(target)) stop("Existing PDF is invalid: ", target)
      cat("FROM EARLIER STAGE OF THIS RUN: ", name, "\n", sep = "")
      written[[length(written) + 1L]] <<- data.frame(pdf = name,
        status = "GENERATED_IN_STAGE", note = note)
      return(invisible(target))
    }
    tmp <- tempfile(".multiome_render_", tmpdir = output_dir, fileext = ".pdf")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    grDevices::pdf(tmp, width = width, height = height, useDingbats = FALSE, onefile = TRUE)
    dev <- grDevices::dev.cur()
    err <- tryCatch({
      draw()
      NULL
    }, error = identity)
    if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev)
    if (inherits(err, "error")) stop("PDF plot error in ", name, ": ", conditionMessage(err))
    if (!check_pdf(tmp)) stop("Empty/invalid regenerated PDF: ", name)
    if (file.exists(target)) {
      if (!isTRUE(overwrite)) stop("Target appeared during rendering: ", target)
      unlink(target)
    }
    if (!file.rename(tmp, target)) stop("Cannot finalize PDF: ", target)
    cat("PASS  ", name, " (", format(file.info(target)$size, big.mark = ","), " bytes)\n", sep = "")
    written[[length(written) + 1L]] <<- data.frame(pdf = name, status = "GENERATED_IN_RENDERER", note = note)
    invisible(target)
  }
  single <- function(name, plot, width, height, note = "FROM_SOURCE_DATA_TABLES") {
    save_pdf(name, function() print(plot), width, height, note)
  }
  joined <- function(name, plots, ncol, width, height, note = "FROM_SOURCE_DATA_TABLES") {
    save_pdf(name, function() gridExtra::grid.arrange(grobs = plots, ncol = ncol),
             width, height, note)
  }
  colors_cluster <- c(`1`="#A6CEE3", `2`="#1F78B4", `3`="#B2DF8A",
                      `4`="#33A02C", `5`="#FB9A99", `6`="#E31A1C",
                      `7`="#FDBF6F", `8`="#FF7F00", `9`="#CAB2D6",
                      `10`="#6A3D9A", Unassigned="#CCCCCC")
  colors_sample <- c(C1="#F8766D", C2="#7CAE00", D="#00BFC4", T="#C77CFF", `NA`="#7F7F7F")
  # 1) WNN, RNA, ATAC embeddings. S9B is the RNA-only SAMPLE embedding.
  embedding <- function(filename) read_data(filename,
        c("cell_id", "UMAP_1", "UMAP_2", "sample", "transferred_cluster"))
  wnn <- embedding("Fig6B_FigS9B_WNN_UMAP.tsv.gz")
  rna <- embedding("FigS9C_RNA_UMAP.tsv.gz")
  atac <- embedding("FigS9C_ATAC_UMAP.tsv.gz")
  if (any(vapply(list(wnn, rna, atac), nrow, integer(1)) != 39787L) ||
      anyDuplicated(wnn$cell_id) ||
      !identical(wnn$cell_id, rna$cell_id) ||
      !identical(wnn$cell_id, atac$cell_id)) stop("Embedding cell IDs or row counts disagree")
  embed_plot <- function(d, group = c("cluster", "sample"), title) {
    group <- match.arg(group)
    d$shown <- if (group == "cluster") as.character(d$transferred_cluster) else as.character(d$sample)
    d$shown[is.na(d$shown)] <- "Unassigned"
    pal <- if (group == "cluster") colors_cluster else colors_sample
    d$shown <- factor(d$shown, levels = names(pal))
    p <- ggplot2::ggplot(d, ggplot2::aes(UMAP_1, UMAP_2, colour = shown)) +
      ggplot2::geom_point(size = 0.15, stroke = 0, alpha = 1) +
      ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
      ggplot2::coord_fixed() + ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = NULL)
    if (group == "cluster") {
      centers <- do.call(rbind, lapply(as.character(1:10), function(id) {
        y <- d[as.character(d$shown) == id, , drop = FALSE]
        data.frame(UMAP_1 = mean(y$UMAP_1), UMAP_2 = mean(y$UMAP_2), label = id)
      }))
      p <- p + ggrepel::geom_label_repel(data = centers,
                  ggplot2::aes(UMAP_1, UMAP_2, label = label),
                  inherit.aes = FALSE, size = 3, seed = 1L, max.overlaps = Inf)
    }
    p
  }
  single("Fig6B_reproduced_REVIEW.pdf", embed_plot(wnn, "cluster", "RNA-ATAC coembedding"), 6, 4.2)
  single("FigS9B_reproduced_REVIEW.pdf", embed_plot(rna, "sample", "Sample distribution: RNA embedding"), 4.2, 4.2)
  joined("FigS9C_reproduced_REVIEW.pdf",
         list(embed_plot(atac, "cluster", "ATAC-only embedding"),
              embed_plot(rna, "cluster", "RNA-only embedding")), 2, 8.4, 4.2)
  rm(wnn, rna, atac)

  # 2) Fig 6C: SPIB expression and chromVAR scores using saved WNN coordinates.
  d <- read_data("Fig6C_SPIB_RNA_chromVAR.tsv.gz",
                 c("UMAP_1", "UMAP_2", "SPIB_SCT_data", "SPIB_chromVAR_raw"))
  spib <- function(z, col, title, low, high) {
    z$color_val <- if (col == "SPIB_chromVAR_raw") pmax(z[[col]], 0) else z[[col]]
    d <- z[order(z$color_val), , drop = FALSE]
    ggplot2::ggplot(d, ggplot2::aes(UMAP_1, UMAP_2, colour = color_val)) +
      ggplot2::geom_point(size = 0.12) +
      ggplot2::scale_colour_gradient(low = low, high = high) +
      ggplot2::coord_fixed() + ggplot2::theme_bw(base_size = 10) +
      ggplot2::labs(title = title, x = "wnnUMAP_1", y = "wnnUMAP_2", colour = NULL)
  }
  joined("Fig6C_SPIB_reproduced_REVIEW.pdf",
         list(spib(d, "SPIB_SCT_data", "sct_SPIB", "lightgrey", "#2166AC"),
              spib(d, "SPIB_chromVAR_raw", "chromvar_MA0081.2", "lightgrey", "darkred")),
         2, 8, 4)
  rm(d)

  # 3) Fig 6D: published clone highlights; full WNN background retained.
  d <- read_data("Fig6D_clone_umap_cells.tsv.gz",
                 c("wnnUMAP_1", "wnnUMAP_2", "is_clone_110_C1", "is_clone_818_D", "is_clone_117_T"))
  cl <- list(list(field="is_clone_110_C1", name="myeloid biased", color="#1737FF", n=110L),
             list(field="is_clone_818_D", name="lymphoid biased", color="#B52AFF", n=818L),
             list(field="is_clone_117_T", name="HSC biased", color="#FF9F00", n=117L))
  cp <- lapply(cl, function(ex) {
    take <- d[d[[ex$field]] == 1L, , drop = FALSE]
    if (nrow(take) != ex$n) stop("Fig 6D clone size changed: ", ex$field)
    ggplot2::ggplot(d, ggplot2::aes(wnnUMAP_1, wnnUMAP_2)) +
      ggplot2::geom_point(color = "#BDBDBD", size = 0.18) +
      ggplot2::geom_point(data = take, colour = ex$color, fill = "white",
                          shape = 21, size = 0.95, stroke = 0.65) +
      ggplot2::coord_fixed() + ggplot2::theme_void(base_size = 10) +
      ggplot2::labs(title = ex$name, subtitle = paste0("n = ", ex$n))
  })
  joined("Fig6D_reproduced_REVIEW.pdf", cp, 3, 11.4, 3.9)
  rm(d, cp)

  # 4) Fig 6F: validated GLMM association coefficients (no model refitting).
  d <- read_data("Fig6F_lineage_module_associations.tsv.gz",
                 c("phenotype", "module", "estimate", "pvalue", "plotted_in_Fig6F"))
  d <- d[!is.na(d$plotted_in_Fig6F) & d$plotted_in_Fig6F, , drop = FALSE]
  if (nrow(d) != 25L) stop("Expected 25 published Fig 6F associations")
  d$phenotype <- factor(d$phenotype,
    levels = c("HSC bias", "myeloid bias", "erythroid bias", "dendritic bias", "lymphoid bias"))
  d$module <- factor(d$module, levels = c("erythroid", "myeloid", "Mk", "lymphoid", "dendritic"))
  d$direction <- ifelse(d$estimate > 0, "up", "down")
  d$logp <- -log10(pmax(d$pvalue, .Machine$double.xmin))
  p <- ggplot2::ggplot(d, ggplot2::aes(module, phenotype, size = abs(estimate),
                stroke = pvalue < 0.05, fill = direction, alpha = logp)) +
    ggplot2::geom_point(shape = 21, colour = "black") +
    ggplot2::scale_fill_manual(values = c(up = "red", down = "blue")) +
    ggplot2::guides(size = "none", alpha = "none") + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = "Lineage-defining peaks", y = "HSC phenotypes") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom", legend.title = ggplot2::element_blank())
  single("Fig6F_lineage_modules_REVIEW.pdf", p, 5, 4.2)
  rm(d, p)

  # 5) Fig 6G: historical motif associations, selection FDR < 0.2 in any phenotype.
  d <- read_data("Fig6G_motif_clone_phenotype_associations.tsv.gz",
                 c("motif", "phenotype", "log2FoldChange", "pvalue", "padj", "selected_by_any_FDR_0.2"))
  phenos <- c("HSC bias", "lymphoid bias", "dendritic bias", "myeloid bias",
              "erythroid bias", "growth rate", "clone size")
  d <- d[d$selected_by_any_FDR_0.2 & d$phenotype %in% phenos, , drop = FALSE]
  if (length(unique(d$motif)) != 27L || nrow(d) != 189L) stop("Unexpected Fig6G selected motifs")
  fx <- stats::reshape(d[, c("motif", "phenotype", "log2FoldChange")],
           idvar = "motif", timevar = "phenotype", direction = "wide")
  vals <- as.matrix(fx[, -1, drop = FALSE]); rownames(vals) <- fx$motif
  d$motif <- factor(d$motif, levels = rownames(vals)[stats::hclust(stats::dist(vals))$order])
  d$phenotype <- factor(d$phenotype, levels = phenos)
  d$direction <- ifelse(d$log2FoldChange > 0, "up", "down")
  p <- ggplot2::ggplot(d, ggplot2::aes(motif, phenotype, size = abs(log2FoldChange),
              fill = direction, alpha = -log10(pmax(pvalue, .Machine$double.xmin)))) +
    ggplot2::geom_point(shape = 21, colour = "black") +
    ggplot2::scale_fill_manual(values = c(up = "red", down = "blue")) +
    ggplot2::guides(size = "none", alpha = "none") + ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "TF motif", y = "HSC phenotype") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom", legend.title = ggplot2::element_blank())
  single("Fig6G_motif_phenotypes_REVIEW.pdf", p, 10, 5)
  rm(d, p, fx, vals)

  # 6) Fig 6H: coverage profiles as stored in the approved numerical export.
  # Gene-model coordinates are not in this table: do not fabricate a gene track.
  d <- read_data("Fig6H_HOXA_genotype_coverage.tsv.gz",
                 c("genomic_position", "genotype", "normalized_signal_plotted", "selected_DNMT3A_peak"))
  if (length(unique(d$selected_DNMT3A_peak)) != 1L ||
      !setequal(unique(d$genotype), c("WT", "TET2", "DNMT3A"))) stop("Fig 6H metadata mismatch")
  peak <- strsplit(unique(d$selected_DNMT3A_peak), "-", fixed = TRUE)[[1]]
  if (length(peak) != 3L || peak[1] != "chr7") stop("Invalid archived HOXA peak ID")
  d$genotype <- factor(d$genotype, levels = c("WT", "TET2", "DNMT3A"))
  p <- ggplot2::ggplot(d, ggplot2::aes(genomic_position, normalized_signal_plotted,
                                      group = genotype, color = genotype)) +
    ggplot2::annotate("rect", xmin = as.numeric(peak[2]), xmax = as.numeric(peak[3]),
                      ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.10) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::facet_grid(genotype ~ ., scales = "fixed") +
    ggplot2::scale_color_manual(values = c(WT="#4DAF4A", TET2="#80A543", DNMT3A="#E41A1C")) +
    ggplot2::guides(color = "none") + ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "chr7 position (bp)", y = "Normalized signal",
                  title = "HOXA locus", subtitle = "Signal-only replot; gene models retained in approved artwork")
  single("Fig6H_HOXA_genotype_REVIEW.pdf", p, 7.2, 5.5,
         note = "COVERAGE_FROM_SOURCE_DATA;FULL_GENE_ANNOTATION_NOT_REPLOTTED")
  rm(d, p)

  # 7) Fig 6I: motif effects by genotype, GC-content coloring.
  d <- read_data("Fig6I_motif_genotype_associations.tsv.gz",
                 c("motif", "motif_gc_percent", "DNMT3A_log2FoldChange", "TET2_log2FoldChange"))
  p <- ggplot2::ggplot(d, ggplot2::aes(DNMT3A_log2FoldChange, TET2_log2FoldChange,
                                      color = motif_gc_percent, label = motif)) +
    ggplot2::geom_point(size = 1.35) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggrepel::geom_text_repel(seed = 1, size = 3.1, max.overlaps = 10) +
    ggplot2::scale_color_gradientn(colours = c("blue", "grey70", "darkred"),
                                   values = scales::rescale(c(min(d$motif_gc_percent),
                                         mean(d$motif_gc_percent), max(d$motif_gc_percent)))) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(x = "logFC, DNMT3A", y = "logFC, TET2", title = "DNMT3A vs TET2", color = "motif GC")
  single("Fig6I_motif_genotypes_REVIEW.pdf", p, 8, 7)
  rm(d, p)

  # 8) Fig 6J: identical plot helper used in the validated prediction script.
  d <- read_data("Fig6J_heritability_summary.tsv",
                 c("modality", "genotype", "phenotype_index", "phenotype_label", "mean_spearman", "sd_spearman"))
  if (nrow(d) != 24L) stop("Fig 6J summary must contain 24 combinations")
  f6j <- new.env(parent = globalenv())
  sys.source(file.path(pub, "reproduce_Fig6J_heritability.R"), envir = f6j)
  single("Fig6J_WT_vs_TET2_REVIEW.pdf", f6j$.f6j_plot(d, "TET2"), 7, 6)
  single("Fig6J_WT_vs_DNMT3A_REVIEW.pdf", f6j$.f6j_plot(d, "DNMT3A"), 7, 6)
  rm(d, f6j)

  # 9) Supplementary Fig 9A all-gene approximation accepted for publication.
  d <- read_data("FigS9A_gene_detection_all_genes.tsv.gz",
                 c("scRNA_detected_fraction", "multiome_detected_fraction", "plotted_on_log_axes"))
  d <- d[!is.na(d$plotted_on_log_axes) & d$plotted_on_log_axes, , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(scRNA_detected_fraction, multiome_detected_fraction)) +
    ggplot2::geom_point(size = 1, alpha = 0.05) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::labs(x = "scRNA-seq", y = "multiome", title = "mRNA coverage",
                  subtitle = "All genes; no DE-based selection")
  single("FigS9A_gene_detection_all_genes_REVIEW.pdf", p, 4, 4)
  rm(d, p)

  # 10) Supplementary Fig 9E: motif enrichment per cluster (saved test results).
  d <- read_data("FigS9E_cluster_motif_enrichment.tsv.gz",
                 c("cluster", "motif.name", "fold.enrichment", "pvalue", "p.adjust"))
  labels <- c(`1`="HSC", `3`="MLP", `4`="Pro-B", `5`="pDC", `7`="early GMP",
              `8`="late GMP", `9`="EMPP/MkP", `10`="EryP")
  plots <- lapply(names(labels), function(id) {
    x <- d[as.character(d$cluster) == id, , drop = FALSE]
    x <- x[is.finite(x$fold.enrichment) & x$fold.enrichment > 0 &
            is.finite(x$pvalue) & x$pvalue > 0, , drop = FALSE]
    if (!nrow(x)) stop("Fig S9E missing motif rows for cluster ", id)
    x$label <- ifelse(is.finite(x$p.adjust) & x$p.adjust < 0.01,
                      as.character(x$motif.name), NA_character_)
    ggplot2::ggplot(x, ggplot2::aes(fold.enrichment, pvalue, label = label)) +
      ggplot2::geom_point(size = 0.6, color = colors_cluster[id]) +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggrepel::geom_text_repel(seed = 1, size = 2.1, max.overlaps = 25,
                               max.time = 0.4, na.rm = TRUE) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(title = paste0("cluster ", id, " - ", labels[id]),
                    x = "motif fold enrichment", y = "p-value")
  })
  joined("FigS9E_cluster_motif_enrichment_REVIEW.pdf", plots, 4, 15, 8.2)
  rm(d, plots)

  # 11) Supplementary Fig 9F: full-precision correlation and column-scaled colors.
  d <- read_data("FigS9F_Corces_Multiome_Pearson.tsv.gz",
                 c("multiome_population", "reference_population", "pearson_r",
                   "color_normalized_within_reference", "best_match_in_reference"))
  pops <- c("HSC", "LMPP-1", "LMPP-2", "MLP", "early GMP", "late GMP",
            "EMPP/MkP", "EryP", "pDC", "Pro-B")
  refs <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP")
  if (nrow(d) != 70L || anyDuplicated(d[, c("multiome_population", "reference_population")]))
    stop("Fig S9F requires exactly 70 unique correlations")
  ind <- cbind(match(d$multiome_population, pops), match(d$reference_population, refs))
  if (anyNA(ind)) stop("Unknown Corces/Multiome population label")
  val <- matrix(NA_real_, nrow = 10L, ncol = 7L, dimnames = list(pops, refs))
  raw <- val
  best <- matrix(FALSE, nrow = 10L, ncol = 7L, dimnames = list(pops, refs))
  val[ind] <- d$color_normalized_within_reference
  raw[ind] <- d$pearson_r
  best[ind] <- as.logical(d$best_match_in_reference)
  if (anyNA(val) || anyNA(raw)) stop("Incomplete S9F heatmap")
  display <- matrix(sprintf("%.2f", raw), nrow = nrow(raw), dimnames = dimnames(raw))
  display[best] <- paste0(display[best], " *")
  ph <- pheatmap::pheatmap(val, cluster_rows = FALSE, cluster_cols = FALSE,
        legend_breaks = c(0, 0.5, 1), legend_labels = c("Low", "Mid", "High"),
        display_numbers = display, number_color = "grey25", border_color = NA,
        angle_col = 90, silent = TRUE,
        main = "Multiome versus purified HSPC ATAC\n34973 reference regions")
  save_pdf("FigS9F_Corces_ATAC_benchmark_REVIEW.pdf", function() {
    grid::grid.newpage()
    grid::grid.draw(gridExtra::arrangeGrob(ph$gtable,
      left = grid::textGrob("Multiome populations", rot = 90),
      bottom = grid::textGrob("HSPC immunophenotypes")))
  }, 5, 5)
  rm(d, ph)

  # 12) Supplementary Fig 9G: published clone-abundance scatter.
  d <- read_data("FigS9G_clone_abundance_RNA_vs_Multiome.tsv.gz",
                 c("sample", "RNA_size", "multiome_size", "plotted_on_log_axes"))
  d <- d[!is.na(d$plotted_on_log_axes) & d$plotted_on_log_axes, , drop = FALSE]
  if (nrow(d) != 485L) stop("S9G matched clone count is not 485")
  d$sample <- factor(as.character(d$sample), levels = names(colors_sample))
  p <- ggplot2::ggplot(d, ggplot2::aes(multiome_size, RNA_size, color = sample)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    ggplot2::geom_point(size = 1.1, alpha = 0.75) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(values = colors_sample, drop = FALSE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = "multiome", y = "RNA", color = NULL, title = "Single clone abundance")
  single("FigS9G_clone_abundance_REVIEW.pdf", p, 5.1, 4.6)
  rm(d, p)

  # 13) Supplementary Fig 10A/F from shared Fig 6G/6I Source Data.
  a <- read_data("Fig6G_motif_clone_phenotype_associations.tsv.gz",
                 c("motif", "phenotype", "log2FoldChange", "pvalue"))
  choose_labels <- function(x, limit = 40L) {
    x$label <- ""
    ok <- which(is.finite(x$pvalue) & x$pvalue < 0.05 & is.finite(x$log2FoldChange))
    if (length(ok)) {
      choose <- head(ok[order(x$pvalue[ok])], limit)
      x$label[choose] <- x$motif[choose]
    }
    x
  }
  plots <- lapply(phenos, function(phenotype) {
    x <- a[a$phenotype == phenotype & is.finite(a$pvalue) & a$pvalue >= 0, , drop = FALSE]
    x <- choose_labels(x)
    x$p_plot <- pmax(x$pvalue, .Machine$double.xmin)
    ggplot2::ggplot(x, ggplot2::aes(log2FoldChange, p_plot)) +
      ggplot2::geom_point(size = 0.65) +
      ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 2,
                                seed = 61010, max.overlaps = Inf, na.rm = TRUE) +
      ggplot2::scale_y_log10() + ggplot2::theme_bw(base_size = 10) +
      ggplot2::labs(title = phenotype, x = "logFC", y = "p-value")
  })
  save_pdf("FigS10A_motif_phenotype_volcano_REVIEW.pdf", function() {
      gridExtra::grid.arrange(grobs = plots,
         layout_matrix = rbind(c(1, 2, 3, NA), c(4, 5, 6, 7)))
    }, 16, 9)
  rm(a, plots)
  f <- read_data("Fig6I_motif_genotype_associations.tsv.gz",
                 c("motif", "motif_gc_percent", "DNMT3A_log2FoldChange", "DNMT3A_pvalue",
                   "TET2_log2FoldChange", "TET2_pvalue"))
  plots <- lapply(c("DNMT3A", "TET2"), function(geno) {
    x <- data.frame(motif = f$motif, motif_gc_percent = f$motif_gc_percent,
      log2FoldChange = f[[paste0(geno, "_log2FoldChange")]],
      pvalue = f[[paste0(geno, "_pvalue")]])
    x <- x[is.finite(x$log2FoldChange) & is.finite(x$pvalue) & x$pvalue >= 0, , drop = FALSE]
    x <- choose_labels(x)
    x$p_plot <- pmax(x$pvalue, .Machine$double.xmin)
    ggplot2::ggplot(x, ggplot2::aes(log2FoldChange, p_plot, colour = motif_gc_percent)) +
      ggplot2::geom_point(size = 0.8) +
      ggrepel::geom_text_repel(ggplot2::aes(label = label), seed = 61010,
                                 size = 2.1, max.overlaps = Inf, na.rm = TRUE) +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_gradientn(colours = c("blue", "grey70", "darkred"),
        values = scales::rescale(c(min(f$motif_gc_percent), mean(f$motif_gc_percent),
                                    max(f$motif_gc_percent))), name = "Motif GC") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = paste0(geno, " vs WT TF motif association"),
                    x = "logFC", y = "p-value")
  })
  joined("FigS10F_motif_genotype_volcano_REVIEW.pdf", plots, 2, 12, 6)
  rm(f, plots)

  out <- do.call(rbind, written)
  rownames(out) <- NULL
  cat("\n=== PDF regeneration summary ===\n")
  print(out, row.names = FALSE)
  cat("\nPDF directory: ", output_dir, "\n", sep = "")
  cat("This PDF renderer does not run analytical tests or fragment recounts.\n")
  cat("Not included: S9D, S10B/C/D/E/G (deferred or other projects); Fig 6A/E are artwork.\n")
  invisible(out)
}
