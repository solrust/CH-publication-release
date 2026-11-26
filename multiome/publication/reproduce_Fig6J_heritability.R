# Fig 6J: historical HSC-to-clone-phenotype prediction (RNA PCA vs ATAC LSI).
# Working directory: CH_publication_release (NOT multiome/).
# Historical basis: multiome/scripts/6_heritability_multiome.R.
# Does NOT alter historical RDS/RData inputs or multiome/figures/.
# Outputs are staged under multiome/publication/; no upstream recomputation.
# Stage A: Fig6J_prepare() once; Stage B: Fig6J_pilot(); Stage C: Fig6J_run().

.f6j_paths <- function(root = getwd(), output_root = NULL) {
  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome")) ||
      !dir.exists(file.path(root, "RNA_clone_analysis"))) {
    stop("Run from the CH_publication_release root.", call. = FALSE)
  }
  if (is.null(output_root)) output_root <- file.path(root, "multiome/publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  original <- normalizePath(file.path(root, "multiome/figures"), mustWork = TRUE)
  if (identical(output_root, original) ||
      startsWith(output_root, paste0(original, "/"))) {
    stop("Refusing to write into original multiome/figures/.", call. = FALSE)
  }
  list(root = root, out = output_root,
       input = file.path(output_root, "intermediate/Fig6J_compact_HSC_inputs.rds"),
       progress = file.path(output_root, "intermediate/Fig6J_split_results_checkpoint.rds"),
       splits = file.path(output_root, "source_data/Fig6J_heritability_splits.tsv.gz"),
       summary = file.path(output_root, "source_data/Fig6J_heritability_summary.tsv"),
       fig_tet2 = file.path(output_root, "figures/Fig6J_WT_vs_TET2_REVIEW.pdf"),
       fig_dnmt3a = file.path(output_root, "figures/Fig6J_WT_vs_DNMT3A_REVIEW.pdf"))
}

.f6j_read_objects <- function(file, requested) {
  e <- new.env(parent = emptyenv())
  # load() does not evaluate the saved objects; their classes can load lazily.
  got <- load(file, envir = e)
  missing <- setdiff(requested, got)
  if (length(missing)) stop("Missing objects in ", basename(file), ": ",
                            paste(missing, collapse = ", "), call. = FALSE)
  mget(requested, envir = e, inherits = FALSE)
}

Fig6J_prepare <- function(obj = NULL, root = getwd(), output_root = NULL) {
  p <- .f6j_paths(root, output_root)
  if (file.exists(p$input)) {
    cat("Compact Fig 6J input checkpoint already exists; not overwriting:\n", p$input, "\n")
    return(invisible(readRDS(p$input)))
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("SeuratObject cannot load in the current R environment.", call. = FALSE)
  if (is.null(obj)) {
    file <- file.path(p$root, "multiome/data/processed/multiome_seurat.rds")
    message("Loading historical Multiome checkpoint once: ", file)
    obj <- readRDS(file)
  } else message("Reusing supplied Multiome checkpoint; no reload.")
  if (!inherits(obj, "Seurat") ||
      !all(c("pca", "lsi") %in% names(obj@reductions)) ||
      !"projection" %in% colnames(obj@meta.data))
    stop("Multiome checkpoint is missing PCA, LSI or projected cell states.", call. = FALSE)

  rna_file <- file.path(p$root, "RNA_clone_analysis/data/processed/clone_info.RData")
  multi_file <- file.path(p$root, "RNA_clone_analysis/data/processed/clone_multiome_info.RData")
  r <- .f6j_read_objects(rna_file, c("clone.info", "clone.states2"))
  m <- .f6j_read_objects(multi_file, c("cell.clone.multi", "clone.info.multi"))
  ri <- r$clone.info
  states <- r$clone.states2
  ci <- m$clone.info.multi
  cellmap <- m$cell.clone.multi
  if (ncol(states) < 10L || !all(c("4", "8", "9", "10") %in% colnames(states)) ||
      !identical(as.character(colnames(states)[c(4, 8, 9, 10)]),
                 c("4", "8", "9", "10"))) {
    stop("Historical phenotype indices 4/8/9/10 do not match clone.states2 columns. Inspect before proceeding.",
         call. = FALSE)
  }
  if (!all(c("cloneid") %in% names(cellmap)) ||
      !"larry" %in% names(ci) || !"genotype" %in% names(ri))
    stop("Historical clone object columns are not as expected.", call. = FALSE)
  if (anyDuplicated(rownames(cellmap))) stop("Repeated cell barcodes in cell.clone.multi.")

  ids <- rownames(cellmap)
  cloneid <- as.character(cellmap$cloneid)
  larry <- as.character(ci$larry[match(cloneid, rownames(ci))])
  state_idx <- match(larry, rownames(states))
  genotype <- as.character(ri$genotype[match(larry, rownames(ri))])
  scores <- as.data.frame(states[state_idx, c(4, 8, 9, 10), drop = FALSE])
  names(scores) <- c("pheno4", "pheno8", "pheno9", "pheno10")
  for (nm in names(scores)) scores[[nm]] <- suppressWarnings(as.numeric(scores[[nm]]))

  # Exact historical cell-name conversion: drop suffix after first '-'.
  seurat_cells <- colnames(obj)
  trimmed <- sub("-.*$", "", seurat_cells)
  if (anyDuplicated(trimmed)) {
    stop("Stripping Seurat cell-name suffix creates duplicate cell IDs; ",
         "historical mapping needs inspection before safely proceeding.", call. = FALSE)
  }
  projected <- as.character(obj$projection)
  names(projected) <- trimmed
  keep <- !is.na(genotype) & !is.na(scores$pheno9) &
          !is.na(projected[ids]) & projected[ids] == "1" &
          genotype %in% c("WT", "TET2", "DNMT3A")
  selected <- ids[keep]
  if (!length(selected)) stop("No matching HSC cells with genotype/phenotypes found.")

  md <- data.frame(cell_id = selected, cloneid = cloneid[keep],
                   larry = larry[keep], genotype = genotype[keep],
                   projected = projected[selected], scores[keep, , drop = FALSE],
                   row.names = NULL, stringsAsFactors = FALSE)
  if (any(!is.finite(as.matrix(md[, c("pheno4", "pheno8", "pheno9", "pheno10")])))) {
    stop("Nonfinite phenotype values in historical HSC selection; ",
         "inspect before reproducing the original lm/cor calls.", call. = FALSE)
  }

  extract_reduction <- function(reduction) {
    z <- obj@reductions[[reduction]]@cell.embeddings
    z_ids <- sub("-.*$", "", rownames(z))
    if (anyDuplicated(z_ids)) stop("Repeated cell IDs in ", reduction, " after trimming.")
    ind <- match(selected, z_ids)
    if (anyNA(ind)) stop("Missing HSC cells in ", reduction, " embedding.")
    z <- z[ind, , drop = FALSE]       # all saved dimensions, as in historical script
    rownames(z) <- selected
    if (ncol(z) != 50L || any(!is.finite(z)))
      stop("Expected 50 finite dimensions for ", reduction)
    z
  }
  inputs <- list(
    metadata = md,
    reductions = list(RNA = extract_reduction("pca"), ATAC = extract_reduction("lsi")),
    historical_parameters = list(cluster_state = 1L, phenotype_indices = c(4L,8L,9L,10L),
                                 train_fraction = 0.9, repetitions = 300L,
                                 model = "lm(y ~ .)", metric = "Spearman correlation"),
    source_files = c("RNA_clone_analysis/data/processed/clone_info.RData",
                     "RNA_clone_analysis/data/processed/clone_multiome_info.RData",
                     "multiome/data/processed/multiome_seurat.rds")
  )
  cat("=== Historical Fig 6J compact input ===\n")
  cat("Clone-mapped HSC cells retained:", nrow(md), "\n")
  print(table(md$genotype))
  cat("RNA PCs:", ncol(inputs$reductions$RNA),
      "; ATAC LSI dimensions:", ncol(inputs$reductions$ATAC), "\n")
  for (g in c("WT", "TET2", "DNMT3A")) {
    n <- sum(md$genotype == g)
    if (floor(n * 0.9) < 52L || (n - floor(n * 0.9)) < 3L)
      stop("Insufficient HSC cells for original 50-predictor model in genotype ", g,
           " (n=", n, ").")
  }
  dir.create(dirname(p$input), recursive = TRUE, showWarnings = FALSE)
  saveRDS(inputs, p$input, compress = "xz")
  message("Saved compact HSC inputs: ", p$input)
  invisible(inputs)
}

.f6j_one_group <- function(inputs, modality, genotype, phenotype,
                           n_rep, seed) {
  md <- inputs$metadata
  red <- inputs$reductions[[modality]]
  pool <- which(md$genotype == genotype)    # HSC restriction applied in prepare()
  y <- md[[paste0("pheno", phenotype)]]
  n <- length(pool)
  n_train <- floor(0.9 * n)
  if (n_train <= ncol(red) + 1L || n - n_train < 3L)
    stop("Insufficient cells for ", modality, " / ", genotype, " / ", phenotype)
  set.seed(seed)
  rho <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    train <- pool[sample.int(n, size = n_train, replace = FALSE)]
    test <- setdiff(pool, train)
    train_df <- data.frame(y = y[train], red[train, , drop = FALSE], check.names = TRUE)
    fit <- stats::lm(y ~ ., data = train_df)
    test_df <- data.frame(red[test, , drop = FALSE], check.names = TRUE)
    pred <- stats::predict(fit, newdata = test_df)
    rho[i] <- suppressWarnings(stats::cor(y[test], pred, method = "spearman"))
    if (!is.finite(rho[i])) stop("Nonfinite prediction correlation in ",
                                modality, "/", genotype, "/", phenotype,
                                ", repetition ", i, call. = FALSE)
  }
  data.frame(modality = modality, genotype = genotype,
             phenotype_index = phenotype,
             phenotype_label = c(`4` = "Lymph", `8` = "GMP", `9` = "MkP", `10` = "Ery")[[as.character(phenotype)]],
             repetition = seq_len(n_rep), seed = seed,
             n_cells = n, n_train = n_train, n_test = n - n_train,
             spearman = rho, stringsAsFactors = FALSE)
}

Fig6J_pilot <- function(n_rep = 5L, seed = 610061L, root = getwd(), output_root = NULL) {
  p <- .f6j_paths(root, output_root)
  if (!file.exists(p$input)) stop("Run Fig6J_prepare() first.")
  inp <- readRDS(p$input)
  n_rep <- as.integer(n_rep)
  if (n_rep < 2L) stop("Pilot requires at least two repetitions.")
  groups <- expand.grid(modality = c("RNA", "ATAC"),
                        genotype = c("WT", "TET2", "DNMT3A"),
                        phenotype = c(4L, 8L, 9L, 10L), stringsAsFactors = FALSE)
  pilot <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    x <- groups[i, ]
    .f6j_one_group(inp, x$modality, x$genotype, x$phenotype,
                   n_rep = n_rep, seed = seed + 101L * i)
  }))
  ag <- aggregate(spearman ~ modality + genotype + phenotype_label,
                  data = pilot, FUN = mean)
  cat("=== Fig 6J pilot means (NOT final Source Data) ===\n")
  print(ag, row.names = FALSE)
  cat("Pilot repetitions per combination:", n_rep,
      "; all 24 combinations completed. No files written.\n")
  invisible(pilot)
}

.f6j_plot <- function(sumtab, mutant) {
  wt <- subset(sumtab, genotype == "WT")
  mu <- subset(sumtab, genotype == mutant)
  names(wt)[names(wt) == "mean_spearman"] <- "WT_mean"
  names(wt)[names(wt) == "sd_spearman"] <- "WT_sd"
  names(mu)[names(mu) == "mean_spearman"] <- "mutant_mean"
  names(mu)[names(mu) == "sd_spearman"] <- "mutant_sd"
  tab <- merge(wt[, c("modality", "phenotype_index", "phenotype_label", "WT_mean", "WT_sd")],
               mu[, c("modality", "phenotype_index", "mutant_mean", "mutant_sd")],
               by = c("modality", "phenotype_index"), sort = FALSE)
  if (nrow(tab) != 8L) stop("Expected 8 plotted RNA/ATAC phenotype points for ", mutant)
  tab$modality <- factor(tab$modality, levels = c("ATAC", "RNA"))
  ggplot2::ggplot(tab, ggplot2::aes(x = WT_mean, y = mutant_mean, color = modality)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey40", linewidth = 1.2) +
    ggplot2::geom_segment(ggplot2::aes(x = WT_mean - WT_sd, xend = WT_mean + WT_sd,
                                      yend = mutant_mean), linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(y = mutant_mean - mutant_sd,
                                      yend = mutant_mean + mutant_sd,
                                      xend = WT_mean), linewidth = 0.4) +
    ggplot2::geom_point(size = 3.2) +
    ggrepel::geom_text_repel(ggplot2::aes(label = phenotype_label),
                            seed = 610061, show.legend = FALSE, size = 4.2,
                            max.overlaps = Inf, min.segment.length = 0) +
    ggplot2::scale_color_manual(values = c(ATAC = "#F8766D", RNA = "#00BFC4")) +
    ggplot2::coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6)) +
    ggplot2::labs(x = "WT HSC to Phenotype, Spearman",
                  y = paste0(mutant, " HSC to Phenotype, Spearman"),
                  title = "HSC to Phenotype heritability", color = "modality") +
    ggplot2::theme_bw(base_size = 13)
}

Fig6J_run <- function(n_rep = 300L, seed = 610061L,
                      root = getwd(), output_root = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggrepel", quietly = TRUE))
    stop("ggplot2 and ggrepel are required for the published-style plots.", call. = FALSE)
  p <- .f6j_paths(root, output_root)
  if (!file.exists(p$input)) stop("Run Fig6J_prepare() first.")
  targets <- c(p$splits, p$summary, p$fig_tet2, p$fig_dnmt3a)
  if (any(file.exists(targets)))
    stop("Final Fig6J outputs already exist; no overwrite:\n",
         paste(targets[file.exists(targets)], collapse = "\n"), call. = FALSE)
  n_rep <- as.integer(n_rep)
  seed <- as.integer(seed)
  if (n_rep != 300L) stop("Final publication run must use the historical 300 repetitions.")
  inp <- readRDS(p$input)
  groups <- expand.grid(modality = c("RNA", "ATAC"),
                        genotype = c("WT", "TET2", "DNMT3A"),
                        phenotype = c(4L, 8L, 9L, 10L), stringsAsFactors = FALSE)
  if (file.exists(p$progress)) {
    cp <- readRDS(p$progress)
    if (!identical(cp$seed, seed) || !identical(cp$n_rep, n_rep))
      stop("An existing progress checkpoint uses different seed/repetition settings.")
  } else cp <- list(seed = seed, n_rep = n_rep, result = list())

  for (i in seq_len(nrow(groups))) {
    x <- groups[i, ]
    key <- paste(x$modality, x$genotype, x$phenotype, sep = "_")
    if (!is.null(cp$result[[key]])) {
      message("Reusing completed group: ", key)
      next
    }
    message("Fitting ", key, " (300 splits; 90/10) ...")
    res <- .f6j_one_group(inp, x$modality, x$genotype, x$phenotype,
                          n_rep = n_rep, seed = seed + 101L * i)
    cp$result[[key]] <- res
    dir.create(dirname(p$progress), recursive = TRUE, showWarnings = FALSE)
    saveRDS(cp, p$progress)
    cat("  mean=", round(mean(res$spearman), 4),
        " sd=", round(stats::sd(res$spearman), 4), "\n", sep = "")
  }
  # Keep the computation and finalization separate. This is also the
  # recovery route if the source-data export failed after all fits finished.
  Fig6J_finalize_saved(root = root, output_root = output_root)
}

# Safely finish Fig 6J using the 24 already-computed groups.
# This does NOT fit models or load the 7 GB Seurat checkpoint.
Fig6J_finalize_saved <- function(root = getwd(), output_root = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggplot2 and ggrepel are needed for Fig 6J plots.", call. = FALSE)
  }
  p <- .f6j_paths(root, output_root)
  if (!file.exists(p$progress)) stop("Missing completed split checkpoint: ", p$progress)
  cp <- readRDS(p$progress)
  if (!identical(as.integer(cp$n_rep), 300L)) {
    stop("Checkpoint does not contain the required 300 repetitions.")
  }
  groups <- expand.grid(modality = c("RNA", "ATAC"),
                        genotype = c("WT", "TET2", "DNMT3A"),
                        phenotype = c(4L, 8L, 9L, 10L),
                        stringsAsFactors = FALSE)
  keys <- paste(groups$modality, groups$genotype, groups$phenotype, sep = "_")
  if (length(cp$result) != length(keys) ||
      !setequal(names(cp$result), keys)) {
    stop("The progress checkpoint does not contain exactly the 24 model groups.")
  }
  splits <- do.call(rbind, cp$result[keys])
  rownames(splits) <- NULL
  required <- c("modality", "genotype", "phenotype_index", "phenotype_label",
                "repetition", "n_cells", "n_train", "n_test", "spearman")
  if (!all(required %in% names(splits)) || nrow(splits) != 7200L ||
      any(!is.finite(splits$spearman))) {
    stop("Checkpoint split results are incomplete or invalid.")
  }
  for (key in keys) {
    d <- cp$result[[key]]
    if (nrow(d) != 300L || !identical(as.integer(d$repetition), 1:300)) {
      stop("Invalid split count/index for ", key)
    }
  }
  inds <- split(seq_len(nrow(splits)),
                interaction(splits$modality, splits$genotype,
                            splits$phenotype_index, drop = TRUE))
  sumtab <- do.call(rbind, lapply(inds, function(ii) {
    d <- splits[ii, , drop = FALSE]
    data.frame(modality = d$modality[1], genotype = d$genotype[1],
               phenotype_index = d$phenotype_index[1],
               phenotype_label = d$phenotype_label[1],
               n_splits = nrow(d), n_cells = d$n_cells[1],
               n_train = d$n_train[1], n_test = d$n_test[1],
               mean_spearman = mean(d$spearman),
               sd_spearman = stats::sd(d$spearman),
               stringsAsFactors = FALSE)
  }))
  rownames(sumtab) <- NULL
  if (nrow(sumtab) != 24L || any(sumtab$n_splits != 300L)) {
    stop("Incomplete Fig 6J summary.")
  }

  # Explicitly close the compressed connection before reading it back.
  # Write to a temporary path first; only replace an incomplete output
  # after its replacement has passed the complete round-trip check.
  dir.create(dirname(p$splits), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(p$fig_tet2), recursive = TRUE, showWarnings = FALSE)
  validate_splits <- function(path) {
    check <- tryCatch({
      con <- gzfile(path, "rt")
      tryCatch(utils::read.delim(con, check.names = FALSE),
               finally = close(con))
    }, error = function(e) NULL)
    if (is.null(check) || nrow(check) != nrow(splits)) return(FALSE)
    for (field in c("modality", "genotype", "phenotype_index", "repetition")) {
      if (!identical(as.character(check[[field]]), as.character(splits[[field]])))
        return(FALSE)
    }
    isTRUE(all.equal(as.numeric(check$spearman), as.numeric(splits$spearman),
                     tolerance = 1e-12))
  }
  if (file.exists(p$splits) && validate_splits(p$splits)) {
    message("Existing compressed Source Data already valid; retaining it.")
  } else {
    temp <- paste0(p$splits, ".tmp-", Sys.getpid())
    on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
    con <- gzfile(temp, "wt")
    tryCatch(utils::write.table(splits, con, sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "NA"),
             finally = close(con))
    if (!validate_splits(temp)) stop("New gzip Source Data did not pass validation.")
    if (!file.rename(temp, p$splits)) {
      if (!file.copy(temp, p$splits, overwrite = TRUE)) {
        stop("Could not replace incomplete gzip Source Data; validated temporary file retained.")
      }
      unlink(temp)
    }
    message("Repaired and validated compressed split Source Data.")
  }

  if (file.exists(p$summary)) {
    old <- utils::read.delim(p$summary, check.names = FALSE)
    if (nrow(old) != nrow(sumtab) ||
        !identical(as.character(old$modality), as.character(sumtab$modality)) ||
        !isTRUE(all.equal(old$mean_spearman, sumtab$mean_spearman,
                          tolerance = 1e-12)) ||
        !isTRUE(all.equal(old$sd_spearman, sumtab$sd_spearman,
                          tolerance = 1e-12))) {
      stop("An existing Fig 6J summary differs from the checkpoint; inspect before overwriting.")
    }
    message("Existing 24-row summary agrees with checkpoint; retaining it.")
  } else {
    utils::write.table(sumtab, p$summary, sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
  }

  for (mut in c("TET2", "DNMT3A")) {
    fig <- if (mut == "TET2") p$fig_tet2 else p$fig_dnmt3a
    if (file.exists(fig)) {
      message("Existing figure retained: ", fig)
      next
    }
    grDevices::pdf(fig, width = 8, height = 8, useDingbats = FALSE)
    tryCatch(print(.f6j_plot(sumtab, mut)), finally = grDevices::dev.off())
    message("Generated: ", fig)
  }
  cat("=== Fig 6J FINALIZED from saved results ===\n")
  cat("Model groups:", nrow(sumtab), "; split correlations:", nrow(splits),
      "; numeric round-trip: PASS\n")
  cat("Split Source Data:", p$splits, "\nSummary:", p$summary,
      "\nFigures:\n", p$fig_tet2, "\n", p$fig_dnmt3a,
      "\nHistorical visual agreement: REVIEW\n")
  invisible(sumtab)
}
