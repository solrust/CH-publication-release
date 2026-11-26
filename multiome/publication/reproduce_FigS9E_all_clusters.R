# Fig S9E: cluster-specific ATAC peak and motif-enrichment reproduction.
# Run from CH_publication_release, using the historical Multiome R environment.
# Based on multiome/scripts/2_multiome_processing1.R. No historical files modified.
# Reuses completed cluster-1 pilot; saves one resumable intermediate checkpoint,
# one full numerical motif table, and one 8-panel review PDF.

reproduce_FigS9E_all_clusters <- function(obj = NULL, root = getwd(), seed = 1L) {
  required <- c("Seurat", "SeuratObject", "Signac", "ggplot2", "ggrepel", "gridExtra")
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Packages unavailable: ", paste(unavailable, collapse = ", "))

  root <- normalizePath(root, mustWork = TRUE)
  project <- file.path(root, "multiome")
  pub <- file.path(project, "publication")
  if (!dir.exists(pub)) stop("Run from CH_publication_release; missing multiome/publication/.")
  historical <- normalizePath(file.path(project, "figures"), mustWork = TRUE)
  stopifnot(!identical(pub, historical))

  checkpoint_file <- file.path(pub, "intermediate", "FigS9E_motif_enrichment_checkpoint.rds")
  pilot_file <- file.path(pub, "intermediate", "FigS9E_cluster1_pilot.rds")
  source_file <- file.path(pub, "source_data", "FigS9E_cluster_motif_enrichment.tsv.gz")
  figure_file <- file.path(pub, "figures", "FigS9E_cluster_motif_enrichment_REVIEW.pdf")
  if (any(file.exists(c(source_file, figure_file)))) {
    stop("Final outputs already exist; leave approved results unchanged:\n",
         paste(c(source_file, figure_file)[file.exists(c(source_file, figure_file))], collapse = "\n"))
  }
  if (!file.exists(pilot_file)) stop("Missing completed cluster-1 pilot: ", pilot_file)
  if (is.null(obj)) {
    message("Loading historical Seurat object once.")
    obj <- readRDS(file.path(project, "data", "processed", "multiome_seurat.rds"))
  } else {
    message("Using supplied Seurat object; checkpoint not reloaded.")
  }
  if (!inherits(obj, "Seurat") || !"ATAC" %in% SeuratObject::Assays(obj)) {
    stop("Expected Seurat object with ATAC assay.")
  }
  if (!all(c("projection", "nCount_ATAC") %in% colnames(obj@meta.data))) {
    stop("Missing expected ATAC/projection metadata.")
  }
  if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
  obj <- SeuratObject::`DefaultAssay<-`(obj, value = "ATAC")
  obj <- SeuratObject::`Idents<-`(obj, value = "projection")
  meta <- obj[["ATAC"]]@meta.features
  if (!"GC.percent" %in% colnames(meta)) stop("Missing GC.percent peak metadata.")

  seed <- as.integer(seed)
  if (length(seed) != 1 || is.na(seed)) stop("seed must be one integer.")
  cluster_ids <- as.character(1:10)
  shown_ids <- as.character(c(1, 3, 4, 5, 7, 8, 9, 10))
  cluster_names <- c("1" = "HSC", "3" = "MLP", "4" = "Pro-B", "5" = "pDC",
                     "7" = "early GMP", "8" = "late GMP", "9" = "EMPP/MkP", "10" = "EryP")
  colors <- c("1" = "#A6CEE3", "2" = "#1F78B4", "3" = "#B2DF8A",
              "4" = "#33A02C", "5" = "#FB9A99", "6" = "#E31A1C",
              "7" = "#FDBF6F", "8" = "#FF7F00", "9" = "#CAB2D6", "10" = "#6A3D9A")

  if (file.exists(checkpoint_file)) {
    state <- readRDS(checkpoint_file)
    if (!identical(state$seed, seed) ||
        !identical(state$clusters, cluster_ids) ||
        !identical(state$n_cells, ncol(obj))) {
      stop("Existing checkpoint is from different input/settings; inspect before resuming.")
    }
    message("Resuming: ", length(state$results), " cluster(s) already saved.")
  } else {
    pilot <- readRDS(pilot_file)
    if (!identical(as.character(pilot$cluster), "1") ||
        !identical(as.integer(pilot$seed), seed) ||
        !all(c("enriched_motifs", "significant_peak_ids", "gc_background_ids") %in% names(pilot))) {
      stop("Cluster-1 pilot does not match this workflow/seed.")
    }
    enriched1 <- pilot$enriched_motifs
    if (!all(c("motif.name", "fold.enrichment", "pvalue", "p.adjust") %in% names(enriched1))) {
      stop("Unexpected motif columns in cluster-1 pilot.")
    }
    state <- list(
      method = "Seurat Wilcoxon cluster vs rest; Signac GC-matched FindMotifs",
      historical_script = "multiome/scripts/2_multiome_processing1.R",
      seed = seed, n_cells = ncol(obj), clusters = cluster_ids,
      FDR_threshold = 0.01, log2FC_threshold = 2,
      GC_background_n = 50000,
      results = list("1" = list(
        cluster = "1", seed = seed,
        n_tested_peaks = pilot$n_tested_peaks,
        significant_peak_ids = pilot$significant_peak_ids,
        gc_background_ids = pilot$gc_background_ids,
        enriched_motifs = enriched1,
        elapsed_DA_seconds = pilot$elapsed_da_seconds,
        elapsed_motif_seconds = pilot$elapsed_motif_seconds
      ))
    )
  }
  dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
  save_checkpoint <- function(x) {
    f <- tempfile(pattern = ".s9e_", tmpdir = dirname(checkpoint_file), fileext = ".rds")
    on.exit(if (file.exists(f)) unlink(f), add = TRUE)
    saveRDS(x, f, compress = TRUE)
    if (!file.rename(f, checkpoint_file)) stop("Could not update checkpoint: ", checkpoint_file)
  }
  if (!file.exists(checkpoint_file)) save_checkpoint(state)

  for (cl in cluster_ids) {
    if (cl %in% names(state$results)) {
      message("Cluster ", cl, ": reusing saved result.")
      next
    }
    message("Cluster ", cl, ": Wilcoxon differential accessibility.")
    dt <- system.time({
      markers <- Seurat::FindMarkers(
        object = obj, ident.1 = cl, test.use = "wilcox", latent.vars = "nCount_ATAC"
      )
    })
    if (!all(c("p_val_adj", "avg_log2FC") %in% colnames(markers))) {
      stop("Unexpected differential-accessibility columns for cluster ", cl)
    }
    selected <- rownames(markers)[
      is.finite(markers$p_val_adj) & markers$p_val_adj < 0.01 &
      is.finite(markers$avg_log2FC) & markers$avg_log2FC > 2
    ]
    message("  ", nrow(markers), " tested peaks; ", length(selected), " selected.")
    result <- list(
      cluster = cl, seed = seed + as.integer(cl) - 1L,
      n_tested_peaks = nrow(markers),
      significant_peak_ids = selected,
      selected_peak_DA = markers[selected, c("avg_log2FC", "p_val_adj"), drop = FALSE],
      gc_background_ids = character(), enriched_motifs = NULL,
      elapsed_DA_seconds = as.numeric(dt["elapsed"]), elapsed_motif_seconds = 0
    )
    rm(markers)
    if (length(selected) > 0) {
      if (!all(selected %in% rownames(meta))) {
        stop("Selected peaks absent from ATAC metadata for cluster ", cl)
      }
      accessible <- intersect(Signac::AccessiblePeaks(obj, idents = cl), rownames(meta))
      if (!length(accessible)) stop("No accessible background peaks for cluster ", cl)
      set.seed(result$seed)
      mt <- system.time({
        bg <- Signac::MatchRegionStats(
          meta.feature = meta[accessible, , drop = FALSE],
          query.feature = meta[selected, , drop = FALSE],
          features.match = "GC.percent", n = 50000
        )
        enriched <- Signac::FindMotifs(obj, features = selected, background = bg)
      })
      req <- c("motif.name", "fold.enrichment", "pvalue", "p.adjust")
      if (!all(req %in% colnames(enriched))) {
        stop("Unexpected motif table columns for cluster ", cl)
      }
      enriched$motif_id <- rownames(enriched)
      enriched$cluster <- cl
      result$gc_background_ids <- bg
      result$enriched_motifs <- enriched
      result$elapsed_motif_seconds <- as.numeric(mt["elapsed"])
      rm(bg, enriched)
    }
    state$results[[cl]] <- result
    save_checkpoint(state)  # Durable progress after every cluster; restartable.
    message("  Saved cluster ", cl, ".")
    gc(verbose = FALSE)
  }

  # Full numerical motif table for all clusters with selected peaks.
  tables <- lapply(state$results, function(x) {
    tab <- x$enriched_motifs
    if (is.null(tab) || !nrow(tab)) return(NULL)
    tab <- as.data.frame(tab)
    if (!"cluster" %in% names(tab)) tab$cluster <- x$cluster
    if (!"motif_id" %in% names(tab)) tab$motif_id <- rownames(tab)
    tab$selected_peaks <- length(x$significant_peak_ids)
    tab$n_tested_peaks <- x$n_tested_peaks
    tab$random_seed <- x$seed
    tab
  })
  tables <- Filter(Negate(is.null), tables)
  if (!length(tables)) stop("No motif enrichments were generated.")
  union_cols <- unique(unlist(lapply(tables, names)))
  tables <- lapply(tables, function(x) {
    for (nm in setdiff(union_cols, names(x))) x[[nm]] <- NA
    x[, union_cols, drop = FALSE]
  })
  combined <- do.call(rbind, tables)
  rownames(combined) <- NULL

  plot_one <- function(cl) {
    x <- state$results[[cl]]$enriched_motifs
    if (is.null(x) || !nrow(x)) {
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::labs(title = paste("cluster", cl, "- no selected peaks")))
    }
    x <- as.data.frame(x)
    x$label <- ifelse(is.finite(x$p.adjust) & x$p.adjust < 0.01,
                      as.character(x$motif.name), NA_character_)
    x <- x[is.finite(x$fold.enrichment) & x$fold.enrichment > 0 &
             is.finite(x$pvalue) & x$pvalue > 0, , drop = FALSE]
    ggplot2::ggplot(x, ggplot2::aes(x = fold.enrichment, y = pvalue, label = label)) +
      ggplot2::geom_point(color = unname(colors[cl]), size = 0.6) +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggrepel::geom_text_repel(na.rm = TRUE, max.overlaps = 25, seed = seed,
                              size = 2.1, max.time = 0.4) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(title = paste0("cluster ", cl, " - ", cluster_names[cl]),
                    x = "motif fold enrichment", y = "p-value")
  }
  grobs <- lapply(shown_ids, plot_one)
  dir.create(dirname(source_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
  tmp_tsv <- tempfile(pattern = ".S9E_motifs_", tmpdir = dirname(source_file), fileext = ".tsv.gz")
  tmp_pdf <- tempfile(pattern = ".S9E_figure_", tmpdir = dirname(figure_file), fileext = ".pdf")
  on.exit({if (file.exists(tmp_tsv)) unlink(tmp_tsv)
           if (file.exists(tmp_pdf)) unlink(tmp_pdf)}, add = TRUE)
  con <- gzfile(tmp_tsv, open = "wt")
  tryCatch(utils::write.table(combined, file = con, sep = "\t", row.names = FALSE,
                              col.names = TRUE, quote = FALSE, na = "NA"),
           finally = close(con))
  reread <- utils::read.delim(tmp_tsv, check.names = FALSE)
  if (nrow(reread) != nrow(combined) ||
      !identical(as.character(reread$motif_id), as.character(combined$motif_id)) ||
      !identical(as.character(reread$cluster), as.character(combined$cluster))) {
    stop("S9E Source Data round-trip validation failed.")
  }
  for (field in c("fold.enrichment", "pvalue", "p.adjust")) {
    original <- suppressWarnings(as.numeric(combined[[field]]))
    restored <- suppressWarnings(as.numeric(reread[[field]]))
    if (!isTRUE(all.equal(original, restored, tolerance = 1e-12,
                          check.attributes = FALSE))) {
      stop("S9E numeric round-trip failed for ", field)
    }
  }
  grDevices::pdf(tmp_pdf, width = 15, height = 8.2, onefile = TRUE)
  tryCatch(gridExtra::grid.arrange(grobs = grobs, ncol = 4, nrow = 2),
           finally = grDevices::dev.off())
  if (file.info(tmp_pdf)$size <= 1000) stop("S9E review PDF appears empty.")
  if (!file.rename(tmp_tsv, source_file)) stop("Cannot finalize Source Data: ", source_file)
  if (!file.rename(tmp_pdf, figure_file)) stop("Cannot finalize review PDF: ", figure_file)

  summary <- do.call(rbind, lapply(cluster_ids, function(cl) {
    x <- state$results[[cl]]
    data.frame(cluster = cl, tested_peaks = x$n_tested_peaks,
               selected_peaks = length(x$significant_peak_ids),
               motifs = if (is.null(x$enriched_motifs)) 0L else nrow(x$enriched_motifs),
               DA_seconds = round(x$elapsed_DA_seconds, 1),
               motifs_seconds = round(x$elapsed_motif_seconds, 1))
  }))
  print(summary, row.names = FALSE)
  message("Source Data: ", source_file)
  message("Review figure: ", figure_file)
  message("Resumable intermediate: ", checkpoint_file)
  invisible(summary)
}
