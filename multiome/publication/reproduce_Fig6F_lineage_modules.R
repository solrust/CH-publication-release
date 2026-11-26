# Fig 6F: publication-only regeneration from archived Multiome/clone checkpoints.
# Source: multiome/scripts/3_mutant_analysis.R (ATAC pseudobulks) and
#         multiome/scripts/5_clone_analysis_final.R (lines 17-56, 103-215).
# WORKING DIRECTORY: CH_publication_release (repository root).
# Original scripts, figures, and data/processed/ are never modified.
#
# Usage in RStudio (with sufficient memory):
#   source("multiome/publication/reproduce_Fig6F_lineage_modules.R")
#   Fig6F_preflight(obj = obj)       # cheap metadata checks
#   Fig6F_make_modules(obj = obj)    # seven DESeq2 contrasts; resumable
#   # After glmmTMB is available in the historical R environment:
#   Fig6F_fit_plot(obj = obj)        # eight phenotypes; final Source Data + PDF
#
# The checkpoint under publication/intermediate is rewritten on successful
# progress (only this new working checkpoint). Final exports refuse overwrite.

.f6f_paths <- function(root = getwd(), output_root = NULL) {
  root <- normalizePath(root, mustWork = TRUE)
  if (!dir.exists(file.path(root, "multiome", "scripts")))
    stop("Working directory must be CH_publication_release.", call. = FALSE)
  if (is.null(output_root)) output_root <- file.path(root, "multiome", "publication")
  output_root <- normalizePath(output_root, mustWork = FALSE)
  historical <- normalizePath(file.path(root, "multiome", "figures"), mustWork = TRUE)
  if (identical(output_root, historical) ||
      startsWith(output_root, paste0(historical, "/")))
    stop("Refusing to write into historical multiome/figures.", call. = FALSE)
  list(root = root, output = output_root,
       checkpoint = file.path(output_root, "intermediate", "Fig6F_lineage_modules_checkpoint.rds"),
       source = file.path(output_root, "source_data", "Fig6F_lineage_module_associations.tsv.gz"),
       figure = file.path(output_root, "figures", "Fig6F_lineage_modules_REVIEW.pdf"))
}

.f6f_required <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "),
                            ". No installation was attempted.", call. = FALSE)
}

.f6f_save <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".Fig6F-in-progress-", tmpdir = dirname(file))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = TRUE)
  if (!file.rename(tmp, file))
    stop("Cannot update working checkpoint: ", file, call. = FALSE)
}

.f6f_load_obj <- function(obj, p) {
  if (!is.null(obj)) {
    if (!inherits(obj, "Seurat")) stop("obj must be a Seurat object.")
    message("Using supplied Seurat object; no reload.")
    return(obj)
  }
  .f6f_required(c("SeuratObject"))
  f <- file.path(p$root, "multiome", "data", "processed", "multiome_seurat.rds")
  if (!file.exists(f)) stop("Missing checkpoint: ", f)
  message("Loading historical Multiome object: ", f)
  readRDS(f)
}

.f6f_load_meta <- function(p, cells) {
  e <- new.env(parent = emptyenv())
  f1 <- file.path(p$root, "RNA_clone_analysis", "data", "processed", "clone_info.RData")
  f2 <- file.path(p$root, "RNA_clone_analysis", "data", "processed", "clone_multiome_info.RData")
  fh <- file.path(p$root, "RNA_analysis", "data", "processed", "hashtag_meta.rds")
  for (f in c(f1, f2, fh)) if (!file.exists(f)) stop("Missing input: ", f)
  load(f1, envir = e)
  load(f2, envir = e)
  require_fields <- c("clone.info", "clone.states2", "clone.info.multi", "cell.clone.multi")
  for (n in require_fields) if (!exists(n, envir = e, inherits = FALSE))
    stop("Input RData missing historical object: ", n)
  hash <- readRDS(fh)
  if (!"exp" %in% names(hash) || is.null(rownames(hash)))
    stop("hashtag_meta.rds lacks rownames or exp column.")
  # Historical 3_mutant_analysis.R appends '-1' to cell.clone.multi rownames.
  cc <- e$cell.clone.multi
  original <- rownames(cc)
  if (is.null(original) || anyDuplicated(original)) stop("Invalid original cell.clone.multi rownames.")
  suffixed <- paste0(original, "-1")
  n_direct <- sum(original %in% cells)
  n_suffix <- sum(suffixed %in% cells)
  mapped <- if (n_suffix > n_direct) suffixed else original
  if (anyDuplicated(mapped) || sum(mapped %in% cells) < 1000)
    stop("Cannot align clone cell IDs to Multiome ATAC cells. Direct matches: ",
         n_direct, "; '-1' matches: ", n_suffix)
  rownames(cc) <- mapped
  mi <- e$clone.info.multi
  rn <- e$clone.info
  for (n in c("larry", "size", "sample")) if (!n %in% names(mi))
    stop("clone.info.multi missing: ", n)
  if (!all(c("id", "genotype") %in% names(rn))) stop("clone.info missing id or genotype.")
  matched_rna <- rn[as.character(mi$larry), c("id", "genotype"), drop = FALSE]
  rownames(matched_rna) <- rownames(mi)
  cc$id <- as.character(matched_rna[as.character(cc$cloneid), "id"])
  cc$larry <- as.character(mi[as.character(cc$cloneid), "larry"])
  cc$sample <- as.character(cc$sample)
  list(e = e, hash = hash, cell = cc, mi = mi, matched_rna = matched_rna,
       n_direct = n_direct, n_suffix = n_suffix)
}

.f6f_atac <- function(obj) {
  if (!"ATAC" %in% names(obj@assays)) stop("No ATAC assay.")
  mat <- SeuratObject::GetAssayData(obj, assay = "ATAC", layer = "counts")
  if (is.null(colnames(mat)) || is.null(rownames(mat)) || anyDuplicated(colnames(mat)))
    stop("ATAC count matrix lacks unique cell or peak identifiers.")
  if (is.null(obj@meta.data$projection)) stop("Missing transferred projection annotations.")
  list(counts = mat, cells = colnames(mat),
       projection = as.character(obj@meta.data[colnames(mat), "projection"]))
}

Fig6F_preflight <- function(obj = NULL, root = getwd()) {
  .f6f_required(c("SeuratObject", "Matrix", "DESeq2"))
  p <- .f6f_paths(root)
  obj <- .f6f_load_obj(obj, p)
  a <- .f6f_atac(obj)
  m <- .f6f_load_meta(p, a$cells)
  tab <- table(a$projection, useNA = "ifany")
  cat("=== Fig6F historical input inspection ===\n")
  cat("ATAC peaks:", nrow(a$counts), "  cells:", ncol(a$counts), "\n")
  cat("Clone cells matching after historical suffix:", m$n_suffix,
      "(without suffix:", m$n_direct, ")\n")
  cat("Matched RNA sample IDs:", sum(!is.na(m$cell$id)), "\n")
  cat("Unique IDs:", paste(sort(unique(na.omit(m$cell$id))), collapse = ", "), "\n")
  cat("Annotation counts:\n"); print(tab)
  cat("Phenotype columns available:",
      paste(intersect(c("1", "4", "5", "8", "9", "10"), colnames(m$e$clone.states2)), collapse = ", "), "\n")
  invisible(data.frame(n_peaks = nrow(a$counts), n_cells = ncol(a$counts),
                       clone_cells_matched = sum(rownames(m$cell) %in% a$cells),
                       rna_id_available = sum(!is.na(m$cell$id))))
}

.f6f_pseudobulks <- function(counts, projection, cc, clusters) {
  # Reproduces the count aggregation in 3_mutant_analysis.R, but uses
  # sparse matrix multiplication to avoid repeatedly copying peak matrices.
  available <- intersect(colnames(counts), rownames(cc))
  sample_ids <- unique(as.character(cc[available, "id"]))
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  if (length(sample_ids) < 2) stop("Too few matched RNA sample IDs.")
  out <- setNames(vector("list", length(clusters)), as.character(clusters))
  for (cluster in clusters) {
    keep <- which(!is.na(projection) & projection == as.character(cluster) &
                    colnames(counts) %in% available)
    keys <- colnames(counts)[keep]
    ids <- as.character(cc[keys, "id"])
    valid <- !is.na(ids) & ids %in% sample_ids
    keep <- keep[valid]; ids <- ids[valid]
    if (!length(keep)) stop("No assigned cells in projection cluster ", cluster)
    present <- unique(ids)
    group <- Matrix::sparseMatrix(i = seq_along(ids), j = match(ids, present),
                                  x = 1, dims = c(length(ids), length(present)),
                                  dimnames = list(NULL, present))
    pb <- counts[, keep, drop = FALSE] %*% group
    out[[as.character(cluster)]] <- as.matrix(pb)
    cat("Cluster", cluster, ":", length(keep), "assigned cells,",
        length(present), "RNA sample groups\n")
  }
  # Historical Reduce('+', clust.atac.list[c(1,2,6)]) adds columns by position.
  # Require identical sample groups rather than silently misalign IDs.
  sets <- lapply(out, colnames)
  reference <- sets[[1]]
  if (!all(vapply(sets, function(x) setequal(x, reference), logical(1))))
    stop("Historical pseudobulk sample groups differ among cell states. ",
         "STOP: do not sum columns positionally; inspect ID coverage.")
  out <- lapply(out, function(x) x[, reference, drop = FALSE])
  out
}

Fig6F_make_modules <- function(obj = NULL, root = getwd(), output_root = NULL) {
  .f6f_required(c("SeuratObject", "Matrix", "DESeq2"))
  p <- .f6f_paths(root, output_root)
  if (file.exists(p$checkpoint)) {
    cp <- readRDS(p$checkpoint)
    if (isTRUE(cp$modules_complete)) {
      message("Fig6F lineage modules already computed; reusing checkpoint: ", p$checkpoint)
      return(invisible(cp$module_sizes))
    }
  } else cp <- list(da = list(), modules_complete = FALSE)
  obj <- .f6f_load_obj(obj, p)
  a <- .f6f_atac(obj)
  m <- .f6f_load_meta(p, a$cells)
  if (!is.null(cp$n_peaks) &&
      (!identical(cp$n_peaks, nrow(a$counts)) || !identical(cp$n_cells, ncol(a$counts))))
    stop("Checkpoint was created from different ATAC dimensions; inspect before resuming.")
  cp$n_peaks <- nrow(a$counts); cp$n_cells <- ncol(a$counts)
  cluster_deriv <- c("1", "4", "8", "5", "9", "10", "3")
  cluster_prog <- c("1", "2", "6")
  clusters <- unique(c(cluster_deriv, cluster_prog))
  pb <- .f6f_pseudobulks(a$counts, a$projection, m$cell, clusters)
  for (cl in cluster_deriv) {
    if (!is.null(cp$da[[cl]])) {
      message("Lineage contrast ", cl, ": reusing checkpoint.")
      next
    }
    message("Lineage contrast ", cl, ": DESeq2, cluster vs pooled 1/2/6.")
    back <- Reduce(`+`, pb[cluster_prog])
    target <- pb[[cl]]
    if (!identical(colnames(target), colnames(back))) stop("Sample-group order mismatch.")
    joined <- cbind(target, back)
    sample_ids <- rep(colnames(target), 2)
    if (any(!sample_ids %in% rownames(m$hash)))
      stop("Some pseudobulk RNA IDs not found in hashtag_meta.rds: ",
           paste(setdiff(sample_ids, rownames(m$hash)), collapse = ", "))
    inf <- m$hash[sample_ids, , drop = FALSE]
    inf$cond <- factor(rep(c("clust", "back"), each = ncol(target)),
                       levels = c("back", "clust"))
    inf$exp <- factor(inf$exp)
    colnames(joined) <- rownames(inf)
    if (!all(is.finite(joined)) || any(joined < 0)) stop("Invalid pseudobulk counts.")
    if (any(colSums(joined) == 0)) stop("Zero-total sample in contrast ", cl)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(joined),
                                          colData = inf, design = ~ cond + exp)
    dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = FALSE)
    coef <- "cond_clust_vs_back"
    if (!coef %in% DESeq2::resultsNames(dds))
      stop("Unexpected coefficient names: ", paste(DESeq2::resultsNames(dds), collapse = ", "))
    rs <- DESeq2::lfcShrink(dds, coef = coef, type = "normal")
    cp$da[[cl]] <- data.frame(log2FoldChange = as.numeric(rs$log2FoldChange),
                               padj = as.numeric(rs$padj),
                               pvalue = as.numeric(rs$pvalue),
                               row.names = rownames(rs))
    rm(dds, rs); gc(verbose = FALSE)
    .f6f_save(cp, p$checkpoint)
    cat("Saved lineage cluster", cl, "\n")
  }
  # Historical threshold and mutually exclusive maximum-effect assignment.
  peak_order <- rownames(cp$da[[cluster_deriv[1]]])
  if (!all(vapply(cluster_deriv, function(k)
      identical(rownames(cp$da[[k]]), peak_order), logical(1))))
    stop("Peak order differs across lineage contrasts; cannot align logFC by position.")
  logfc <- do.call(cbind, lapply(cluster_deriv, function(k) cp$da[[k]]$log2FoldChange))
  rownames(logfc) <- rownames(cp$da[[cluster_deriv[1]]])
  colnames(logfc) <- cluster_deriv
  max_lfc <- apply(logfc, 1, max) # identical NA semantics to historical code
  peak_sets <- setNames(lapply(cluster_deriv, function(k) {
    d <- cp$da[[k]]
    which_keep <- which(!is.na(d$padj) & d$padj < 0.1 &
                          is.finite(d$log2FoldChange) & d$log2FoldChange >= 1 &
                          d$log2FoldChange >= max_lfc)
    intersect(rownames(d)[which_keep], rownames(a$counts))
  }), cluster_deriv)
  peak_sets <- peak_sets[lengths(peak_sets) > 0L]
  required <- c("4", "8", "5", "9", "10")
  if (!all(required %in% names(peak_sets)))
    stop("Missing published lineage modules: ",
         paste(setdiff(required, names(peak_sets)), collapse = ", "),
         ". Check DE contrasts before moving to GLM.")
  mod_counts <- do.call(rbind, lapply(peak_sets, function(ids)
    Matrix::colSums(a$counts[ids, , drop = FALSE])))
  rownames(mod_counts) <- names(peak_sets)
  colnames(mod_counts) <- a$cells
  cp$peak_sets <- peak_sets
  cp$module_counts <- mod_counts
  cp$module_sizes <- data.frame(cluster = names(peak_sets),
                                n_peaks = lengths(peak_sets), row.names = NULL)
  cp$modules_complete <- TRUE
  .f6f_save(cp, p$checkpoint)
  cat("=== Published lineage peak modules ready ===\n")
  print(cp$module_sizes, row.names = FALSE)
  cat("Checkpoint:", p$checkpoint, "\n")
  invisible(cp$module_sizes)
}

.f6f_phenotype_meta <- function(p, m) {
  edna <- new.env(parent = emptyenv())
  f <- file.path(p$root, "DNA_analysis", "data", "processed", "DNA_clones.RData")
  if (!file.exists(f)) stop("Missing DNA clone rate file: ", f)
  load(f, envir = edna)
  if (!exists("rate", envir = edna, inherits = FALSE))
    stop("DNA_clones.RData has no rate object.")
  clons <- as.character(m$mi$larry)
  raw <- m$mi
  rate <- edna$rate[clons]
  r <- m$e$clone.info
  # In the historical script this also appends a second column named 'size'
  # (RNA-relative clone size), while the original Multiome 'size' column is
  # used by cells_partition's size.cut >= 50. Preserve that column order.
  summary.size <- aggregate(r[clons, "size"], by = list(r[clons, "id"]), sum)
  rownames(summary.size) <- as.character(summary.size[, 1])
  denom <- summary.size[as.character(r[clons, "id"]), 2]
  pheno <- cbind(raw,
                 m$matched_rna[, c("id", "genotype"), drop = FALSE],
                 rate = as.numeric(rate),
                 size = log(as.numeric(r[clons, "size"]) / denom),
                 m$e$clone.states2[clons, c("1", "4", "5", "8", "9", "10"), drop = FALSE])
  rownames(pheno) <- rownames(raw)
  if (anyDuplicated(rownames(pheno))) stop("Duplicate Multiome clone IDs.")
  pheno
}

.f6f_grouped <- function(module_counts, library_sizes, projection, cc,
                         high_cells, low_cells, cluster_ids = c("1", "2", "6")) {
  out <- list()
  for (condition in c("1", "2")) {
    chosen <- if (condition == "1") low_cells else high_cells
    for (cluster in cluster_ids) {
      eligible <- intersect(chosen, colnames(module_counts)[!is.na(projection) &
                                   projection == cluster])
      eligible <- intersect(eligible, rownames(cc))
      sample_ids <- unique(as.character(cc[eligible, "id"]))
      sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
      for (sample in sample_ids) {
        ids <- eligible[as.character(cc[eligible, "id"]) == sample]
        y <- rowSums(module_counts[, ids, drop = FALSE])
        size <- sum(library_sizes[ids])
        if (size <= 0) next
        out[[length(out) + 1]] <- list(clone = condition, clust = cluster,
                                       batch = sample, counts = y,
                                       size = size, n_cells = length(ids))
      }
    }
  }
  if (!length(out)) stop("No positive-size pseudobulk groups.")
  frame <- data.frame(clone = vapply(out, `[[`, character(1), "clone"),
                      clust = vapply(out, `[[`, character(1), "clust"),
                      batch = vapply(out, `[[`, character(1), "batch"),
                      size = vapply(out, `[[`, numeric(1), "size"),
                      n_cells = vapply(out, `[[`, integer(1), "n_cells"))
  frame$clone <- factor(frame$clone, levels = c("1", "2"))
  frame$clust <- factor(frame$clust, levels = cluster_ids)
  ymat <- do.call(cbind, lapply(out, `[[`, "counts"))
  rownames(ymat) <- rownames(module_counts)
  list(info = frame, counts = ymat)
}

Fig6F_fit_plot <- function(obj = NULL, root = getwd(), output_root = NULL) {
  .f6f_required(c("SeuratObject", "Matrix", "glmmTMB", "ggplot2"))
  p <- .f6f_paths(root, output_root)
  if (file.exists(p$source) || file.exists(p$figure))
    stop("Approved/previous outputs protected. Existing file(s):\n",
         paste(c(p$source, p$figure)[file.exists(c(p$source, p$figure))],
               collapse = "\n"), call. = FALSE)
  if (!file.exists(p$checkpoint)) stop("Run Fig6F_make_modules first.")
  cp <- readRDS(p$checkpoint)
  if (!isTRUE(cp$modules_complete)) stop("Lineage peak modules are incomplete.")
  obj <- .f6f_load_obj(obj, p)
  a <- .f6f_atac(obj)
  m <- .f6f_load_meta(p, a$cells)
  if (!identical(colnames(cp$module_counts), a$cells))
    stop("Saved module counts do not match the current object's cell order.")
  phenotype <- .f6f_phenotype_meta(p, m)
  lib_sizes <- Matrix::colSums(a$counts)
  names(lib_sizes) <- a$cells
  pnames <- c("1" = "HSC bias", "4" = "lymphoid bias",
              "5" = "dendritic bias", "8" = "myeloid bias",
              "10" = "erythroid bias", "9" = "mk bias",
              "rate" = "growth rate", "size" = "clone size")
  fate <- c("1" = "HSC", "4" = "lymphoid", "8" = "myeloid",
            "5" = "dendritic", "9" = "Mk", "10" = "erythroid", "3" = "cluster3")
  results <- list()
  for (ph in names(pnames)) {
    cat("Phenotype:", pnames[[ph]], "\n")
    scores <- as.numeric(phenotype[[ph]])
    eligible <- is.finite(scores) & as.numeric(phenotype$size) >= 50
    if (!any(eligible)) stop("No eligible clones for phenotype: ", ph)
    med <- stats::median(scores[eligible], na.rm = TRUE)
    clone_hi <- as.character(phenotype$larry[eligible & scores >= med])
    clone_lo <- as.character(phenotype$larry[eligible & scores < med])
    cc <- m$cell
    high_cells <- rownames(cc)[!is.na(cc$larry) & cc$larry %in% clone_hi]
    low_cells <- rownames(cc)[!is.na(cc$larry) & cc$larry %in% clone_lo]
    sums <- .f6f_grouped(cp$module_counts, lib_sizes, a$projection,
                         cc, high_cells, low_cells)
    info <- sums$info
    if (length(unique(info$clone)) < 2 || length(unique(info$clust)) < 2)
      stop("Cannot fit both conditions and progenitor clusters for: ", ph)
    for (fate_id in rownames(sums$counts)) {
      dat <- info
      dat$y <- as.numeric(sums$counts[fate_id, ])
      fit <- glmmTMB::glmmTMB(y ~ clone + clust + offset(log(size)),
                               data = dat, family = glmmTMB::nbinom2(link = "log"))
      co <- summary(fit)$coefficients$cond
      if (!"clone2" %in% rownames(co))
        stop("GLM did not produce expected high-vs-low clone coefficient for ", ph)
      z <- co["clone2", ]
      results[[length(results) + 1]] <- data.frame(
        phenotype = pnames[[ph]], module = unname(fate[[fate_id]]),
        module_cluster = fate_id, estimate = unname(z[1]),
        std_error = unname(z[2]), z_value = unname(z[3]),
        pvalue = unname(z[4]), n_pseudobulks = nrow(dat),
        n_cells_in_groups = sum(dat$n_cells), stringsAsFactors = FALSE)
    }
  }
  full <- do.call(rbind, results)
  rownames(full) <- NULL
  # Published Figure 6F has 5 lineage modules x 5 phenotype rows.
  published_modules <- c("erythroid", "myeloid", "Mk", "lymphoid", "dendritic")
  published_phenotypes <- c("HSC bias", "myeloid bias", "erythroid bias",
                            "dendritic bias", "lymphoid bias")
  full$plotted_in_Fig6F <- full$module %in% published_modules &
                           full$phenotype %in% published_phenotypes
  subset <- full[full$plotted_in_Fig6F, , drop = FALSE]
  if (nrow(subset) != 25L || anyDuplicated(subset[c("module", "phenotype")]))
    stop("Expected 5 x 5 published module/phenotype combinations; inspect results.")
  subset$module <- factor(subset$module, levels = published_modules)
  subset$phenotype <- factor(subset$phenotype, levels = published_phenotypes)
  subset$direction <- ifelse(subset$estimate > 0, "up", "down")
  subset$negative_log10_p <- -log10(pmax(subset$pvalue, .Machine$double.xmin))
  subset$significant <- subset$pvalue < 0.05
  plot <- ggplot2::ggplot(subset,
         ggplot2::aes(module, phenotype, size = abs(estimate),
                      stroke = significant, fill = direction,
                      alpha = negative_log10_p)) +
    ggplot2::geom_point(colour = "black", shape = 21) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = c(up = "red", down = "blue")) +
    ggplot2::guides(size = "none", alpha = "none") +
    ggplot2::labs(x = "Lineage-defining peaks", y = "HSC phenotypes") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom", legend.title = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 12))
  dir.create(dirname(p$source), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(p$figure), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(p$source, "wt")
  tryCatch(utils::write.table(full, con, sep = "\t", row.names = FALSE,
                              col.names = TRUE, quote = FALSE, na = "NA"),
           finally = close(con))
  check <- utils::read.delim(p$source, check.names = FALSE)
  if (!identical(as.character(full$phenotype), check$phenotype) ||
      !isTRUE(all.equal(full$estimate, check$estimate, tolerance = 1e-12)) ||
      !isTRUE(all.equal(full$pvalue, check$pvalue, tolerance = 1e-12)))
    stop("Source Data round-trip validation failed.")
  grDevices::pdf(p$figure, width = 5, height = 4.2, useDingbats = FALSE)
  tryCatch(print(plot), finally = grDevices::dev.off())
  cat("Fig6F Source Data:", p$source, "\n")
  cat("Fig6F review figure:", p$figure, "\n")
  cat("Total module/phenotype results:", nrow(full), "; plotted:", nrow(subset), "\n")
  invisible(full)
}
