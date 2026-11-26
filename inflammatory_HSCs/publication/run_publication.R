# CH revision -- inflammatory_HSCs (Supplementary Fig. 10B and 10C)
# Run from CH_publication_release. Historical inputs / figures are read-only.
# Two separately executable stages because the scRNA and Multiome Seurat checkpoints are large.

.ihsc_required <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "),
                            ". No installation attempted.", call. = FALSE)
}
.ihsc_root <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  if (!file.exists(file.path(root, "inflammatory_HSCs/scripts/analysis.R")))
    stop("Use CH_publication_release as working directory (root=...).", call. = FALSE)
  root
}
.ihsc_path <- function(root, path) file.path(root, path)
.ihsc_expected <- function(root) {
  c(
    "inflammatory_HSCs/Zeng2026_HSCiM/hsc_gene_set.rds",
    "inflammatory_HSCs/Zeng2026_HSCiM/Source_Data_Extended_Data_Fig5.xlsx",
    "Association_analysis/data/processed/hsc/hsc_axes.RData",
    "RNA_analysis/data/processed/scrna.RData",
    "RNA_clone_analysis/data/processed/clone_info.RData",
    "DNA_analysis/data/processed/DNA_clones.RData",
    "multiome/data/processed/multiome_seurat.rds"
  )
}

ihsc_preflight <- function(root = getwd()) {
  root <- .ihsc_root(root)
  paths <- .ihsc_expected(root)
  full <- file.path(root, paths)
  ans <- data.frame(input = paths, exists = file.exists(full),
                    GiB = round(ifelse(file.exists(full), file.info(full)$size / 1024^3, NA_real_), 3),
                    stringsAsFactors = FALSE)
  cat("=== Historical inputs ===\n")
  print(ans, row.names = FALSE)
  packages <- c("ggplot2", "glmnet", "pROC", "Matrix", "readxl", "ggrepel", "SeuratObject", "Signac")
  cat("\n=== R package availability (no installations) ===\n")
  print(data.frame(package = packages, available = vapply(packages, requireNamespace, logical(1), quietly = TRUE)), row.names = FALSE)
  cat("\n=== Production panel ownership ===\n",
      "S10B: scRNA discriminant + ROC + RNA HSC embedding + clone associations\n",
      "S10C: Multiome discriminant + ROC + RNA-only HSC embedding + motif comparison\n",
      "No upstream download / multiome_input.RData loading / model run in preflight.\n", sep = "")
  invisible(ans)
}

.ihsc_ensure_outputs <- function(output_root, files) {
  hist <- normalizePath(file.path(getwd(), "inflammatory_HSCs/figures"), mustWork = FALSE)
  output_root <- normalizePath(output_root, mustWork = FALSE)
  if (identical(output_root, hist) || startsWith(output_root, paste0(hist, "/")))
    stop("Refusing to write into historical inflammatory_HSCs/figures.")
  actual <- file.path(output_root, files)
  if (any(file.exists(actual)))
    stop("Output(s) already exist; first-pass reference preserved:\n",
         paste(actual[file.exists(actual)], collapse = "\n"),
         "\nUse a different output_root for an independent run.", call. = FALSE)
  dir.create(file.path(output_root, "source_data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_root, "figures"), recursive = TRUE, showWarnings = FALSE)
  output_root
}
.ihsc_tsv <- function(df, path) {
  con <- gzfile(path, "wt")
  tryCatch(
    write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"),
    finally = close(con)
  )
  chk <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(chk) != nrow(df) || !identical(names(chk), names(df)))
    stop("Source Data round-trip failed: ", path)
  invisible(path)
}
.ihsc_pdf <- function(path, width, height, plot_fun) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  tryCatch(plot_fun(), finally = grDevices::dev.off())
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0)
    stop("Empty plot file: ", path)
  invisible(path)
}
.ihsc_signature <- function(root) {
  g <- readRDS(file.path(root, "inflammatory_HSCs/Zeng2026_HSCiM/hsc_gene_set.rds"))
  stopifnot(all(c("hscI.genes", "hscII.genes") %in% names(g)))
  g
}
.ihsc_model <- function(hsc.sig.expr, gsets, n_pc, label) {
  .ihsc_required(c("glmnet", "pROC"))
  z <- scale(hsc.sig.expr, center = TRUE, scale = TRUE)
  keep <- colSums(!is.finite(z)) == 0L
  z <- z[, keep, drop = FALSE]
  cat(label, ": ", nrow(z), " HSC cells, ", ncol(z), " variable signature genes\n", sep = "")
  if (ncol(z) <= n_pc) stop(label, ": fewer than ", n_pc+1, " usable signature genes")
  pca <- stats::prcomp(z, center = FALSE, scale. = FALSE)
  pc <- seq_len(n_pc)
  gene.pc <- sweep(pca$rotation[, pc, drop = FALSE], 2, pca$sdev[pc], "*")
  class <- rep(NA_integer_, nrow(gene.pc))
  class[rownames(gene.pc) %in% gsets$hscI.genes] <- 0L
  class[rownames(gene.pc) %in% gsets$hscII.genes] <- 1L
  if (anyNA(class) || length(unique(class)) != 2L) stop(label, ": unexpected gene labels")
  set.seed(123)
  cvfit <- glmnet::cv.glmnet(gene.pc, class, family = "binomial", alpha = 0,
                              standardize = TRUE, nfolds = 10)
  fit.coef <- coef(cvfit, s = "lambda.1se")
  beta <- as.numeric(fit.coef[-1, 1]); names(beta) <- colnames(gene.pc)
  gene.score <- drop(gene.pc %*% beta)
  gene.df <- data.frame(gene = rownames(gene.pc), score = gene.score,
                        class = factor(class, levels = c(0, 1), labels = c("HSC-I", "HSC-II")))
  roc <- pROC::roc(gene.df$class, gene.df$score,
                   levels = c("HSC-I", "HSC-II"), direction = "<", quiet = TRUE)
  roc.df <- data.frame(FPR = 1 - roc$specificities, TPR = roc$sensitivities)
  cell.axis <- beta * pca$sdev[pc]
  names(cell.axis) <- colnames(pca$x[, pc, drop = FALSE])
  score <- drop(pca$x[, pc, drop = FALSE] %*% cell.axis)
  names(score) <- rownames(pca$x)
  list(genes = gene.df, roc = roc.df, score = score,
       AUC = as.numeric(pROC::auc(roc)))
}
.ihsc_discriminant_figures <- function(fit, dir_fig, prefix) {
  .ihsc_required("ggplot2")
  gg <- asNamespace("ggplot2")
  box <- gg$ggplot(fit$genes, gg$aes(x = class, y = score, fill = class)) +
    gg$geom_boxplot(outlier.shape = NA) + gg$geom_jitter(width = 0.15, size = 1) +
    gg$theme_bw() + gg$labs(x = NULL, y = "HSC-I/HSC-II discriminant score") +
    gg$theme(legend.position = "none", axis.title = gg$element_text(size = 20),
             axis.text = gg$element_text(size = 16))
  .ihsc_pdf(file.path(dir_fig, paste0(prefix, "_gene_discriminant_boxplot.pdf")), 5, 5,
            function() print(box))
  roc <- gg$ggplot(fit$roc, gg$aes(FPR, TPR)) + gg$geom_line(linewidth = 1.5) +
    gg$geom_abline(slope = 1, intercept = 0, linetype = 2) +
    gg$coord_equal() + gg$theme_bw() +
    gg$labs(x = "False positive rate", y = "True positive rate",
            title = sprintf("HSC-I vs HSC-II genes, AUC = %.2f", fit$AUC)) +
    gg$theme(axis.title = gg$element_text(size = 20), axis.text = gg$element_text(size = 16))
  .ihsc_pdf(file.path(dir_fig, paste0(prefix, "_gene_discriminant_ROC.pdf")), 5, 5,
            function() print(roc))
}
.ihsc_embedding <- function(fit, emb, out, xlim, ylim) {
  cls <- names(fit$score)
  stopifnot(all(cls %in% rownames(emb)), ncol(emb) >= 2)
  coords <- as.matrix(emb[cls, 1:2, drop = FALSE])
  pts <- data.frame(cell_id = cls, UMAP_1 = coords[, 1], UMAP_2 = coords[, 2],
                    HSCiM_score = unname(fit$score), stringsAsFactors = FALSE)
  colors <- grDevices::colorRampPalette(c("blue", "grey90", "red"))(10)[
    as.integer(cut(rank(fit$score), breaks = 10, include.lowest = TRUE))]
  .ihsc_pdf(out, 6, 6, function() {
    plot(coords, pch = 19, cex = 0.5, col = grDevices::adjustcolor(colors, 0.5),
         xlab = "", ylab = "", xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim,
         cex.main = 2, font.main = 1, main = "HSC-II / HSC-I axis")
  })
  pts
}

.ihsc_clone_correlations <- function(root, hsc_scores, fig_dir) {
  .ihsc_required(c("Matrix", "ggplot2"))
  e <- new.env(parent = globalenv())
  load(file.path(root, "RNA_clone_analysis/data/processed/clone_info.RData"), envir = e)
  load(file.path(root, "DNA_analysis/data/processed/DNA_clones.RData"), envir = e)
  stopifnot(all(c("clone.info", "clone.states2", "cell.clone") %in% ls(e)),
            "rate" %in% ls(e))
  clons <- rownames(e$clone.states2)
  clone_info <- e$clone.info[clons, , drop = FALSE]
  summary.size <- aggregate(clone_info$size, list(clone_info$id), sum)
  rownames(summary.size) <- summary.size[, 1]
  clone.size <- log(clone_info$size / summary.size[clone_info$id, 2])
  clone.pheno <- cbind("growth rate" = e$rate[clons],
                       "clone size" = clone.size,
                       e$clone.states2[clons, c(1, 4, 5, 8, 10), drop = FALSE])
  rownames(clone.pheno) <- clons
  colnames(clone.pheno) <- c("growth rate", "clone size", "HSC bias", "lymphoid bias",
                             "dendritic bias", "myeloid bias", "erythroid bias")
  shared.cells <- intersect(names(hsc_scores), rownames(e$cell.clone))
  cell.df <- data.frame(cell_id = shared.cells,
                        cloneid = as.character(e$cell.clone[shared.cells, "cloneid"]),
                        HSCiM = unname(hsc_scores[shared.cells]))
  clone.hscim <- aggregate(HSCiM ~ cloneid, data = cell.df, FUN = mean, na.rm = TRUE)
  n.hsc <- table(cell.df$cloneid)
  clone.hscim$n_HSC <- as.integer(n.hsc[clone.hscim$cloneid])
  rownames(clone.hscim) <- clone.hscim$cloneid
  clone.hscim <- clone.hscim[clone.hscim$n_HSC >= 5, , drop = FALSE]
  shared.clones <- intersect(rownames(clone.hscim), rownames(clone.pheno))
  combined <- data.frame(clone_id = shared.clones,
                         HSCiM = clone.hscim[shared.clones, "HSCiM"],
                         n_HSC = clone.hscim[shared.clones, "n_HSC"],
                         clone.pheno[shared.clones, , drop = FALSE],
                         check.names = FALSE, row.names = NULL)
  if (!nrow(combined)) stop("No shared HSC clones; check historical clone identifiers.")
  res <- do.call(rbind, lapply(colnames(clone.pheno), function(ph) {
    x <- combined$HSCiM; y <- combined[[ph]]
    valid <- is.finite(x) & is.finite(y)
    if (sum(valid) < 4) stop("Too few valid clones for ", ph)
    fit <- suppressWarnings(cor.test(x[valid], y[valid], method = "spearman", exact = FALSE))
    data.frame(phenotype = ph, n_clones = sum(valid), rho = unname(fit$estimate), pvalue = fit$p.value)
  }))
  res$padj <- p.adjust(res$pvalue, "BH")
  res <- res[order(res$pvalue), , drop = FALSE]
  res$logP <- pmin(-log10(pmax(res$pvalue, .Machine$double.xmin)), 10)
  res$program <- "HSC-II"
  gg <- asNamespace("ggplot2")
  p <- gg$ggplot(res, gg$aes(y = phenotype, x = program, size = logP, fill = rho)) +
    gg$geom_point(shape = 21, color = "black", stroke = 0.5) +
    gg$scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    gg$scale_size_continuous(range = c(0, 20), limits = c(0, 10), breaks = c(1, 2, 5, 10)) +
    gg$theme_minimal() + gg$labs(y = "clonal phenotype", x = NULL,
                                  size = "-log10(P)", fill = "Spearman rho") +
    gg$theme(axis.title.y = gg$element_text(size = 20),
             axis.text.y = gg$element_text(size = 15),
             axis.text.x = gg$element_text(size = 15, angle = 30, vjust = 1, hjust = 1))
  .ihsc_pdf(file.path(fig_dir, "FigS10B_HSCiM_clones.pdf"), 5, 5, function() print(p))
  list(clones = combined, correlation = res)
}

.ihsc_run_B <- function(root, output_root) {
  cat("=== S10B scRNA HSC-iM ===\n")
  needed <- c("source_data/FigS10B_gene_discriminant.tsv.gz",
              "source_data/FigS10B_ROC.tsv.gz",
              "source_data/FigS10B_HSC_axis_cells.tsv.gz",
              "source_data/FigS10B_clone_scores.tsv.gz",
              "source_data/FigS10B_clone_correlations.tsv.gz",
              "figures/FigS10B_gene_discriminant_boxplot.pdf",
              "figures/FigS10B_gene_discriminant_ROC.pdf",
              "figures/FigS10B_HSC_axis.pdf",
              "figures/FigS10B_HSCiM_clones.pdf")
  output_root <- .ihsc_ensure_outputs(output_root, needed)
  gsets <- .ihsc_signature(root)
  hsc <- new.env(parent = globalenv()); load(file.path(root, "Association_analysis/data/processed/hsc/hsc_axes.RData"), envir = hsc)
  rna <- new.env(parent = globalenv()); load(file.path(root, "RNA_analysis/data/processed/scrna.RData"), envir = rna)
  if (!all(c("p2counts", "emb") %in% ls(rna)) || !"prr" %in% ls(hsc))
    stop("Historical RNA checkpoint missing p2counts/emb or HSC prr.")
  cls <- rownames(hsc$prr$x)
  sig.use <- intersect(union(gsets$hscI.genes, gsets$hscII.genes), colnames(rna$p2counts))
  if (!all(cls %in% rownames(rna$p2counts))) stop("HSC cells missing from RNA counts")
  expr <- as.matrix(rna$p2counts[cls, sig.use, drop = FALSE])
  fit <- .ihsc_model(expr, gsets, n_pc = 100L, label = "S10B")
  rm(expr, hsc); gc(verbose = FALSE)
  pts <- .ihsc_embedding(fit, rna$emb, file.path(output_root, "figures/FigS10B_HSC_axis.pdf"),
                         xlim = c(-2, 5), ylim = c(2, 10))
  rm(rna); gc(verbose = FALSE)
  .ihsc_discriminant_figures(fit, file.path(output_root, "figures"), "FigS10B")
  clone <- .ihsc_clone_correlations(root, fit$score, file.path(output_root, "figures"))
  .ihsc_tsv(fit$genes, file.path(output_root, needed[1]))
  .ihsc_tsv(fit$roc, file.path(output_root, needed[2]))
  .ihsc_tsv(pts, file.path(output_root, needed[3]))
  .ihsc_tsv(clone$clones, file.path(output_root, needed[4]))
  .ihsc_tsv(clone$correlation, file.path(output_root, needed[5]))
  cat(sprintf("S10B: AUC=%.4f; %d HSC cells; %d matched clones; 5 Source Data tables; 4 PDFs\n",
              fit$AUC, nrow(pts), nrow(clone$clones)))
  invisible(list(AUC = fit$AUC, cells = nrow(pts), clones = nrow(clone$clones)))
}

.ihsc_run_C <- function(root, output_root) {
  cat("=== S10C Multiome HSC-iM ===\n")
  .ihsc_required(c("ggplot2", "glmnet", "pROC", "readxl", "ggrepel", "SeuratObject", "Signac"))
  needed <- c("source_data/FigS10C_gene_discriminant.tsv.gz",
              "source_data/FigS10C_ROC.tsv.gz",
              "source_data/FigS10C_HSC_axis_cells.tsv.gz",
              "source_data/FigS10C_chromVAR_motif_associations.tsv.gz",
              "source_data/FigS10C_published_motif_comparison.tsv.gz",
              "figures/FigS10C_gene_discriminant_boxplot.pdf",
              "figures/FigS10C_gene_discriminant_ROC.pdf",
              "figures/FigS10C_HSC_axis.pdf",
              "figures/FigS10C_motif_comparison.pdf")
  output_root <- .ihsc_ensure_outputs(output_root, needed)
  gsets <- .ihsc_signature(root)
  seurat <- readRDS(file.path(root, "multiome/data/processed/multiome_seurat.rds"))
  cls <- colnames(seurat)[!is.na(seurat$projection) & as.character(seurat$projection) == "1"]
  # Historical published S10C: SCT counts normalized by each cell's TOTAL
  # SCT counts, before signature-gene standardization in .ihsc_model().
  # Use all SCT count features in the denominator, not only signature genes.
  sct <- seurat[["SCT"]]
  sig.use <- intersect(union(gsets$hscI.genes, gsets$hscII.genes), rownames(sct@counts))
  depth <- Matrix::colSums(sct@counts[, cls, drop = FALSE])
  stopifnot(length(depth) == length(cls), all(is.finite(depth)), all(depth > 0))
  expr <- t(as.matrix(sct@counts[sig.use, cls, drop = FALSE]))
  expr <- sweep(expr, 1L, depth[rownames(expr)], "/")
  stopifnot(all(is.finite(expr)))
  fit <- .ihsc_model(expr, gsets, n_pc = 20L, label = "S10C")
  coords <- SeuratObject::Embeddings(seurat[["umap.rna"]])
  pts <- .ihsc_embedding(fit, coords, file.path(output_root, "figures/FigS10C_HSC_axis.pdf"),
                         xlim = c(-5, -1), ylim = c(-8, 0))
  .ihsc_discriminant_figures(fit, file.path(output_root, "figures"), "FigS10C")
  # Same linear-model calculation as the historical script, with per-motif FDR adjustment.
  chromvar <- as.matrix(seurat[["chromvar"]]@data[, cls, drop = FALSE])
  model <- stats::lm(t(chromvar) ~ fit$score[cls])
  coef <- stats::coef(model)[2, ]
  # summary.mlm() prefixes response names (often "Response <motif>").
  # Do not index those prefixed names directly with the motif IDs: that
  # would silently turn all P values and FDRs into NA.
  coefficient_matrix <- stats::coef(model)
  pval <- vapply(summary(model), function(x) x$coefficients[2, "Pr(>|t|)"], numeric(1))
  motif.ids <- rownames(chromvar)
  if (length(pval) != ncol(coefficient_matrix) ||
      !setequal(colnames(coefficient_matrix), motif.ids)) {
    stop("S10C motif model responses do not align with chromVAR motif IDs")
  }
  # The summary list follows the response column order of the fitted model.
  names(pval) <- colnames(coefficient_matrix)
  motif.names <- Signac::ConvertMotifID(seurat[["ATAC"]], id = motif.ids)
  names(motif.names) <- motif.ids  # historical mapping is unnamed before assignment
  motif.res <- data.frame(motif = motif.ids, TF = unname(motif.names[motif.ids]),
                          beta = as.numeric(coef[motif.ids]),
                          pval = as.numeric(pval[motif.ids]), stringsAsFactors = FALSE)
  motif.res$padj <- p.adjust(motif.res$pval, "BH")
  n_valid <- sum(is.finite(motif.res$pval))
  n_sig <- sum(motif.res$padj < 0.05, na.rm = TRUE)
  if (n_valid < 0.95 * nrow(motif.res) || n_sig == 0L) {
    stop("S10C motif P-value/FDR alignment failed; refusing to write incomplete Source Data")
  }
  cat(sprintf("S10C motif tests: %d/%d finite P values; %d motifs with FDR < 0.05\n",
              n_valid, nrow(motif.res), n_sig))
  rm(seurat, expr, chromvar, model, sct); gc(verbose = FALSE)
  reference <- as.data.frame(readxl::read_excel(
    file.path(root, "inflammatory_HSCs/Zeng2026_HSCiM/Source_Data_Extended_Data_Fig5.xlsx"),
    sheet = "ExtFig5g"))
  for (field in c("Motif_JASPAR", "Motif_Name", "FoldEnrichment_HSCiMvsHSCI"))
    if (!field %in% names(reference)) stop("Missing published motif field: ", field)
  cmp <- merge(motif.res, reference, by.x = "motif", by.y = "Motif_JASPAR")
  if (!nrow(cmp)) stop("No shared motif IDs with published reference")
  cmp <- cmp[is.finite(cmp$beta) & is.finite(cmp$FoldEnrichment_HSCiMvsHSCI) &
               cmp$FoldEnrichment_HSCiMvsHSCI > 0, , drop = FALSE]
  ct <- suppressWarnings(cor.test(cmp$beta, cmp$FoldEnrichment_HSCiMvsHSCI, method = "spearman", exact = FALSE))
  cmp$log_fold_enrichment <- log(cmp$FoldEnrichment_HSCiMvsHSCI)
  cmp$significance <- factor(ifelse(cmp$padj < 0.01, "FDR < 0.01",
                                    ifelse(cmp$padj < 0.05, "FDR < 0.05", "NS")),
                             levels = c("NS", "FDR < 0.05", "FDR < 0.01"))
  cmp$label <- ifelse(cmp$padj < 0.05, cmp$Motif_Name, NA_character_)
  gg <- asNamespace("ggplot2")
  p <- gg$ggplot(cmp, gg$aes(beta, log_fold_enrichment)) +
    gg$geom_point(gg$aes(color = significance), size = 1.5) +
    gg$scale_color_manual(values = c("NS" = "grey75", "FDR < 0.05" = "lightcoral", "FDR < 0.01" = "red")) +
    ggrepel::geom_text_repel(gg$aes(label = label), size = 4, na.rm = TRUE, max.overlaps = 15, seed = 123) +
    gg$annotate("text", x = -Inf, y = Inf,
                label = sprintf("Spearman = %.2f\nP = %.2g", unname(ct$estimate), ct$p.value),
                hjust = -0.1, vjust = 1.1, size = 5) +
    gg$theme_bw(base_size = 18) +
    gg$labs(x = "Motif association with HSC-iM score",
            y = "log(published HSC-iM motif enrichment)", color = NULL)
  .ihsc_pdf(file.path(output_root, "figures/FigS10C_motif_comparison.pdf"), 7, 7, function() print(p))
  .ihsc_tsv(fit$genes, file.path(output_root, needed[1]))
  .ihsc_tsv(fit$roc, file.path(output_root, needed[2]))
  .ihsc_tsv(pts, file.path(output_root, needed[3]))
  .ihsc_tsv(motif.res, file.path(output_root, needed[4]))
  .ihsc_tsv(cmp, file.path(output_root, needed[5]))
  cat(sprintf("S10C: AUC=%.4f; %d HSC cells; %d matched motifs; Spearman=%.3f; 5 Source Data tables; 4 PDFs\n",
              fit$AUC, nrow(pts), nrow(cmp), unname(ct$estimate)))
  invisible(list(AUC = fit$AUC, cells = nrow(pts), motifs = nrow(cmp), spearman = unname(ct$estimate)))
}

# Atomic per-panel staging: a failed first pass leaves no partial approved outputs.
.ihsc_stage_panel <- function(panel, root, output_root) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  hist <- normalizePath(file.path(root, "inflammatory_HSCs/figures"), mustWork = TRUE)
  dest <- normalizePath(output_root, mustWork = TRUE)
  if (identical(dest, hist) || startsWith(dest, paste0(hist, "/")))
    stop("Cannot write into historical figure directory")
  existing <- unlist(lapply(c("source_data", "figures"), function(dir) {
    location <- file.path(dest, dir)
    if (!dir.exists(location)) return(character())
    list.files(location, pattern = paste0("^Fig", panel, "_"), full.names = TRUE)
  }), use.names = FALSE)
  if (length(existing)) stop("Panel outputs already exist; they are protected:\n",
                             paste(existing, collapse = "\n"))
  stage <- tempfile(paste0(".stage_", panel, "_"), tmpdir = dest)
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  result <- if (panel == "S10B") .ihsc_run_B(root, stage) else .ihsc_run_C(root, stage)
  files <- list.files(stage, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  if (!length(files)) stop("No outputs in staging area")
  for (src in files) {
    relative <- substring(src, nchar(stage) + 2L)
    target <- file.path(dest, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(target) || !file.rename(src, target))
      stop("Cannot finalize staged file: ", target)
  }
  cat("COMMITTED ", panel, ": ", length(files),
      " Source Data and figure files under ", dest, "\n", sep = "")
  invisible(result)
}

run_inflammatory_publication <- function(panels = c("S10B", "S10C"),
                                         root = getwd(),
                                         output_root = NULL) {
  root <- .ihsc_root(root)
  if (is.null(output_root)) output_root <- file.path(root, "inflammatory_HSCs/publication")
  panels <- match.arg(panels, c("S10B", "S10C"), several.ok = TRUE)
  result <- list()
  for (panel in panels) result[[panel]] <- .ihsc_stage_panel(panel, root, output_root)
  cat("Finished requested panel(s): ", paste(panels, collapse = ", "), "\n", sep = "")
  invisible(result)
}
