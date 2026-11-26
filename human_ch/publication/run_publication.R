# human_ch publication reproduction: Figure 3G + Supplementary Figure 5E
# Source from CH_publication_release root in R 4.2.3 / /packages Matrix 1.6-4.
# Checkpoint export only: does NOT refit GSEA or rebuild Seurat integration.

hc_require <- function() {
  for (p in c("ggplot2", "ggrepel")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing R package: ", p)
  }
}

hc_hsc_pathways <- c(
  "dormant (vs active) HSCs", "active (vs dormant) HSCs",
  "AML stem cell dormant", "AML stem cell active",
  "quiescent  (vs activated) HSCs", "activated (vs quiescent) HSCs",
  "GMP", "LMPP"
)

hc_gsea_source <- function(g) {
  stopifnot(all(c("gsea_xeno_dnmt3a", "gsea_xeno_tet2", "gsea_human_ch", "gsea_dnmt3a_r882") %in% names(g)))
  hnames <- c("DNMT3A_env", "TET2_env", "DNMT3A_mut", "TET2_mut")
  stopifnot(all(hnames %in% names(g$gsea_human_ch)))
  groups <- c(list("hDNMT3A_R882" = g$gsea_dnmt3a_r882,
                   "DNMT3A_vs_CTL, xeno" = g$gsea_xeno_dnmt3a,
                   "TET2_vs_CTL, xeno" = g$gsea_xeno_tet2), g$gsea_human_ch[hnames])
  main <- unlist(lapply(groups, `[[`, "mainPathways"), use.names = FALSE)
  main <- as.character(main)
  recur <- names(which(table(main) > 1L))
  uni <- unique(main)
  # Exact historical pathway-selection rule in scripts/4_comparative_DE.R.
  selected <- unique(c(uni[grepl("M5", uni)], hc_hsc_pathways, recur))
  selected <- setdiff(selected, c(
    "M5930: HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "M5915: HALLMARK_APICAL_JUNCTION", "M5898: HALLMARK_DNA_REPAIR",
    "GO:0048589: developmental growth",
    "R-HSA-71291: Metabolism of amino acids and derivatives",
    "GO:0090670: RNA localization to Cajal body",
    "AML stem cell dormant", "AML stem cell active", "active (vs dormant) HSCs"
  ))
  out <- lapply(names(groups), function(nm) {
    x <- groups[[nm]]$fgsea_res
    req <- c("pathway", "NES", "pval", "padj")
    if (!all(req %in% names(x))) stop("Missing fgsea columns for ", nm)
    x <- as.data.frame(x)[, intersect(c("pathway", "NES", "pval", "padj", "size"), names(x)), drop = FALSE]
    x <- x[!is.na(x$pathway) & x$pathway %in% selected, , drop = FALSE]
    x$data <- nm
    x
  })
  tab <- do.call(rbind, out)
  rownames(tab) <- NULL
  if (!"size" %in% names(tab)) tab$size <- NA_integer_
  catg <- rep("Signaling/Other", nrow(tab))
  catg[grepl("LMPP|GMP", tab$pathway)] <- "Bias"
  catg[grepl("HSCs", tab$pathway)] <- "Dormancy"
  catg[grepl("IL2|INFLAM|INTERFERON|TNFA", tab$pathway)] <- "Inflammation"
  catg[grepl("OXIDATIVE|FATTY_ACID|Translation", tab$pathway)] <- "Metabolism"
  catg[grepl("MYC|G2M|MITOTIC", tab$pathway)] <- "Proliferation"
  tab$Category <- catg
  tab$minus_log10_FDR <- -log10(pmax(tab$padj, .Machine$double.xmin))
  tab <- tab[, c("data", "Category", "pathway", "NES", "pval", "padj", "size", "minus_log10_FDR")]
  tab <- tab[order(tab$data, tab$Category, tab$pathway), , drop = FALSE]
  rownames(tab) <- NULL
  stopifnot(!anyDuplicated(tab[c("data", "pathway")]), nrow(tab) > 0L)
  selection <- data.frame(pathway = sort(selected),
                          recurrent_main_pathway = sort(selected) %in% recur,
                          historical_M5_match = grepl("M5", sort(selected)),
                          historical_HSC_list = sort(selected) %in% hc_hsc_pathways,
                          stringsAsFactors = FALSE)
  list(genesets = tab, selection = selection)
}

hc_plot_3g <- function(df, filename) {
  d <- df
  d$Category <- factor(d$Category, levels = c("Bias", "Dormancy", "Inflammation", "Metabolism", "Proliferation", "Signaling/Other"))
  d$data <- factor(d$data, levels = sort(unique(d$data)))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = NES, y = pathway)) +
    ggplot2::geom_point(ggplot2::aes(color = NES, size = minus_log10_FDR)) +
    ggplot2::scale_color_gradientn(colors = c("darkblue", "grey90", "darkred")) +
    ggplot2::facet_grid(Category ~ data, scales = "free_y", space = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::scale_size(range = c(0.5, 12), breaks = c(1, 2, 3, 5, 10),
                        labels = c("1", "2 (p=0.01)", "3", "5 (p<1e-5)", "10")) +
    ggplot2::labs(size = "-log10(FDR)", color = "NES", title = "Pathway Comparison")
  ggplot2::ggsave(filename, p, width = 16, height = 6, units = "in")
}

hc_cluster_labels <- c(`1` = "HSC", `2` = "LMPP-1", `3` = "MLP", `4` = "Pro-B",
                       `5` = "pDC", `6` = "LMPP-2", `7` = "early GMP", `8` = "late GMP",
                       `9` = "EMPP/MkP", `10` = "EryP")

hc_s5e_source <- function(meta) {
  required <- c("Cluster", "clust.final", "umaprpca_1", "umaprpca_2", "exp", "dataset")
  stopifnot(is.data.frame(meta), all(required %in% names(meta)),
            !anyNA(meta[, c("umaprpca_1", "umaprpca_2")]),
            !is.null(rownames(meta)), !anyDuplicated(rownames(meta)))
  t <- data.frame(cell_id = rownames(meta),
                  UMAP1 = as.numeric(meta$umaprpca_1), UMAP2 = as.numeric(meta$umaprpca_2),
                  exp = as.character(meta$exp), dataset = as.character(meta$dataset),
                  human_annotation = as.character(meta$Cluster),
                  xeno_cluster = as.character(meta$clust.final),
                  stringsAsFactors = FALSE)
  t$xeno_annotation <- unname(hc_cluster_labels[t$xeno_cluster])
  t$cell_source <- ifelse(!is.na(t$dataset), "xenograft",
                          ifelse(!is.na(t$human_annotation), "human annotated", "human unannotated"))
  s <- as.data.frame(table(t$cell_source), stringsAsFactors = FALSE)
  names(s) <- c("cell_source", "n_cells")
  s <- s[s$n_cells > 0, , drop = FALSE]
  stopifnot(nrow(t) == 34073L, sum(!is.na(t$human_annotation)) > 0,
            sum(!is.na(t$xeno_cluster)) > 0)
  list(cells = t, summary = s)
}

hc_plot_s5e <- function(d, filename) {
  n <- !is.na(d$human_annotation)
  h <- d[n, , drop = FALSE]
  cen_h <- aggregate(h[, c("UMAP1", "UMAP2")], by = list(Annotation = h$human_annotation), FUN = mean)
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = d[!n, ], ggplot2::aes(UMAP1, UMAP2), size = 0.45, color = "grey80") +
    ggplot2::geom_point(data = h, ggplot2::aes(UMAP1, UMAP2, color = human_annotation), size = 0.45) +
    ggrepel::geom_text_repel(data = cen_h, ggplot2::aes(UMAP1, UMAP2, label = Annotation), size = 3.2, max.overlaps = Inf) +
    ggplot2::theme_bw() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size=7)) +
    ggplot2::labs(x = "", y = "")
  n2 <- !is.na(d$xeno_cluster)
  x <- d[n2, , drop = FALSE]
  cen_x <- aggregate(x[, c("UMAP1", "UMAP2")], by = list(cluster = x$xeno_cluster), FUN = median)
  cen_x$label <- paste0(cen_x$cluster, ": ", unname(hc_cluster_labels[cen_x$cluster]))
  # Approximate display palette. Replace with original old.scrna$clust.final.cols
  # when exact historical color matching is essential (not available in metadata).
  palette <- c(`1`="#bad7ea", `2`="#447fce", `3`="#a9df87", `4`="#55b445", `5`="#f17e92",
               `6`="#bf70bb", `7`="#ee722c", `8`="#fba440", `9`="#8769aa", `10`="#d6afd7")
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = d[!n2, ], ggplot2::aes(UMAP1, UMAP2), size = 0.45, color = "grey80") +
    ggplot2::geom_point(data = x, ggplot2::aes(UMAP1, UMAP2, color = xeno_cluster), size = 0.45) +
    ggrepel::geom_text_repel(data = cen_x, ggplot2::aes(UMAP1, UMAP2, label = label), size = 3.2, max.overlaps = Inf) +
    ggplot2::scale_color_manual(values = palette, breaks = names(hc_cluster_labels), labels = hc_cluster_labels) +
    ggplot2::theme_bw() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size=7)) +
    ggplot2::labs(x = "", y = "")
  grDevices::pdf(filename, width = 17, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1), newpage = FALSE)
  print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2), newpage = FALSE)
  invisible(filename)
}

hc_write_tsv <- function(x, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con), add=TRUE)
  utils::write.table(x, con, sep="\t", quote=TRUE, row.names=FALSE, col.names=TRUE, na="NA")
}

human_ch_export <- function(panels = c("3G", "S5E"), root = getwd(), validate_replay = TRUE) {
  hc_require()
  panels <- unique(panels)
  stopifnot(length(panels) > 0, all(panels %in% c("3G", "S5E")))
  stopifnot(dir.exists(file.path(root, "human_ch", "data")))
  pub <- file.path(root, "human_ch", "publication")
  sd <- file.path(pub, "source_data")
  figs <- file.path(pub, "figures")
  valdir <- file.path(pub, "validation")
  for (dir in c(sd, figs, valdir)) dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  output <- list()
  if ("3G" %in% panels) {
    raw <- readRDS(file.path(root, "human_ch", "data", "datasets_DE_GSEA.rds"))
    a <- hc_gsea_source(raw)
    if (validate_replay) stopifnot(isTRUE(all.equal(a, hc_gsea_source(readRDS(file.path(root, "human_ch", "data", "datasets_DE_GSEA.rds"))), check.attributes=FALSE)))
    output[["3G"]] <- a
  }
  if ("S5E" %in% panels) {
    raw <- readRDS(file.path(root, "human_ch", "data", "integrated_seurat_metadata.rds"))
    a <- hc_s5e_source(raw)
    if (validate_replay) stopifnot(isTRUE(all.equal(a, hc_s5e_source(readRDS(file.path(root, "human_ch", "data", "integrated_seurat_metadata.rds"))), check.attributes=FALSE)))
    output[["S5E"]] <- a
  }
  specs <- list()
  if ("3G" %in% panels) specs[["3G"]] <- c("source_data/Fig3G_GSEA_pathways.tsv.gz", "source_data/Fig3G_pathway_selection.tsv",
                                               "figures/Fig3G_GSEA_comparison_REVIEW.pdf")
  if ("S5E" %in% panels) specs[["S5E"]] <- c("source_data/FigS5E_integrated_UMAP_cells.tsv.gz", "source_data/FigS5E_cell_group_summary.tsv",
                                                 "figures/FigS5E_integrated_UMAP_REVIEW.pdf")
  targets <- unlist(specs, use.names=FALSE)
  existing <- file.path(pub, targets)[file.exists(file.path(pub, targets))]
  if (length(existing)) stop("Refusing to overwrite existing publication output(s): ", paste(existing, collapse=", "))
  staged <- tempfile("human_ch_stage_", tmpdir=pub)
  dir.create(staged)
  on.exit(unlink(staged, recursive=TRUE), add=TRUE)
  dir.create(file.path(staged, "source_data"))
  dir.create(file.path(staged, "figures"))
  if ("3G" %in% panels) {
    a <- output[["3G"]]
    hc_write_tsv(a$genesets, file.path(staged, specs[["3G"]][1]))
    hc_write_tsv(a$selection, file.path(staged, specs[["3G"]][2]))
    hc_plot_3g(a$genesets, file.path(staged, specs[["3G"]][3]))
  }
  if ("S5E" %in% panels) {
    a <- output[["S5E"]]
    hc_write_tsv(a$cells, file.path(staged, specs[["S5E"]][1]))
    hc_write_tsv(a$summary, file.path(staged, specs[["S5E"]][2]))
    hc_plot_s5e(a$cells, file.path(staged, specs[["S5E"]][3]))
  }
  stopifnot(all(file.exists(file.path(staged, targets))), all(file.info(file.path(staged, targets))$size > 0))
  if (!all(file.copy(file.path(staged, targets), file.path(pub, targets), overwrite=FALSE))) stop("File promotion failed; inspect publication outputs")
  v <- do.call(rbind, lapply(names(output), function(panel) {
    a <- output[[panel]]
    n <- vapply(a, nrow, integer(1))
    data.frame(panel=panel, table=names(n), n_rows=as.integer(n),
               status=if (validate_replay) "CHECKPOINT_REPLAY_PASS" else "NOT_TESTED", stringsAsFactors=FALSE)
  }))
  report <- file.path(valdir, "Fig3G_S5E_checkpoint_replay.tsv")
  if (!file.exists(report)) hc_write_tsv(v, report)
  print(v, row.names=FALSE)
  cat("Saved PDFs and Source Data. Historical artwork visual review: PENDING.\n")
  invisible(list(status="NUMERICAL_CHECKPOINT_PASS; VISUAL_REVIEW_PENDING", report=v, paths=file.path(pub, targets)))
}

# === BEGIN human_ch S5F checkpoint extension (2026-11-10) ===
# S5F: cell-state pseudobulk means of saved integration RNA scale.data layers,
# and Pearson correlation between human reference and xenograft states.
# Historical calculate_cluster_centroids() definition was unavailable; explicit
# per-gene arithmetic means over clusters with >=10 cells are reconstructed here.

hc_s5f_centroids <- function(expression_matrix, labels, min_cells=10L) {
  stopifnot(length(dim(expression_matrix)) == 2L,
            !is.null(rownames(expression_matrix)),
            !is.null(colnames(expression_matrix)))
  if (anyDuplicated(rownames(expression_matrix)) || anyDuplicated(colnames(expression_matrix)))
    stop('S5F duplicated gene/cell IDs')
  labels <- labels[colnames(expression_matrix)]
  labels <- as.character(labels)
  cell_count <- table(labels[!is.na(labels) & nzchar(labels)])
  keep <- names(cell_count[cell_count >= min_cells])
  if (!length(keep)) stop('S5F no groups with >=', min_cells, ' cells')
  means <- vapply(keep, function(k) {
    cols <- which(!is.na(labels) & labels == k)
    rowMeans(expression_matrix[, cols, drop=FALSE], na.rm=FALSE)
  }, numeric(nrow(expression_matrix)))
  rownames(means) <- rownames(expression_matrix)
  colnames(means) <- keep
  stopifnot(!anyNA(means), all(is.finite(means)))
  list(mean_expression=means,
       counts=data.frame(cluster=keep, n_cells=as.integer(cell_count[keep]),
                         stringsAsFactors=FALSE))
}

hc_s5f_source <- function(integrated) {
  if (!requireNamespace('SeuratObject', quietly=TRUE)) stop('S5F needs SeuratObject')
  # Exact two stored integrated-layer matrices used by the historical script.
  assay <- integrated@assays$RNA
  x <- tryCatch(assay$scale.data.1, error=function(e) NULL)
  y <- tryCatch(assay$scale.data.2, error=function(e) NULL)
  if (is.null(x) || is.null(y)) stop('S5F missing historical scale.data.1/scale.data.2 layers')
  meta <- integrated@meta.data
  if (!all(c('clust.final','Cluster') %in% names(meta))) stop('S5F expected metadata missing')
  stopifnot(all(colnames(x) %in% rownames(meta)), all(colnames(y) %in% rownames(meta)))
  # Orientation and assignment match 2_run_integrations_comparison.R.
  xs <- hc_s5f_centroids(x, setNames(meta$clust.final, rownames(meta)))
  ys <- hc_s5f_centroids(y, setNames(meta$Cluster, rownames(meta)))
  shared <- intersect(rownames(xs$mean_expression), rownames(ys$mean_expression))
  if (length(shared) < 100) stop('S5F too few shared scale-data genes: ', length(shared))
  xmean <- xs$mean_expression[shared,,drop=FALSE]
  ymean <- ys$mean_expression[shared,,drop=FALSE]
  cor_xy <- cor(xmean, ymean, method='pearson', use='pairwise.complete.obs')
  # Required finite pattern; do not coerce undefined correlations to zero.
  stopifnot(identical(rownames(cor_xy), colnames(xmean)),
            identical(colnames(cor_xy), colnames(ymean)))
  # Independent correlation calculation through standardized columns
  # (not a second independent integration / preprocessing).
  zsc <- function(m) {
    m <- sweep(m, 2L, colMeans(m), '-')
    sweep(m, 2L, sqrt(colSums(m*m)/(nrow(m)-1)), '/')
  }
  z1 <- zsc(xmean); z2 <- zsc(ymean)
  check <- crossprod(z1, z2)/(nrow(xmean)-1)
  finite <- is.finite(cor_xy) & is.finite(check)
  if (!identical(is.finite(cor_xy), is.finite(check)) ||
      (any(finite) && max(abs(cor_xy[finite]-check[finite])) > 1e-8))
    stop('S5F Pearson correlation internal check failed')
  grid <- expand.grid(xeno_cluster=rownames(cor_xy),
                      reference_cluster=colnames(cor_xy), stringsAsFactors=FALSE)
  grid$pearson_r <- as.vector(cor_xy)
  grid$n_shared_genes <- length(shared)
  grid$xeno_cells <- xs$counts$n_cells[match(grid$xeno_cluster, xs$counts$cluster)]
  grid$reference_cells <- ys$counts$n_cells[match(grid$reference_cluster, ys$counts$cluster)]
  grid <- grid[order(grid$xeno_cluster, grid$reference_cluster),,drop=FALSE]
  rownames(grid) <- NULL
  cell_counts <- rbind(data.frame(dataset='xenograft', xs$counts),
                       data.frame(dataset='reference', ys$counts))
  centroids <- data.frame(gene=shared,
                          setNames(as.data.frame(xmean, check.names=FALSE), paste0('xeno_', colnames(xmean))),
                          setNames(as.data.frame(ymean, check.names=FALSE), paste0('reference_', colnames(ymean))),
                          check.names=FALSE)
  list(correlations=grid, cell_counts=cell_counts, gene_centroids=centroids,
       matrix=cor_xy)
}

hc_s5f_pdf <- function(m, path) {
  if (!requireNamespace('pheatmap', quietly=TRUE)) stop('S5F needs pheatmap')
  plotmat <- t(m)
  if (exists('hc_cluster_labels', mode='numeric') || exists('hc_cluster_labels', mode='character')) {
    lab <- hc_cluster_labels[colnames(plotmat)]
    lab[is.na(lab)] <- colnames(plotmat)[is.na(lab)]
    colnames(plotmat) <- paste(names(lab), lab, sep=' - ')
  }
  grDevices::pdf(path, width=8, height=8)
  on.exit(grDevices::dev.off(), add=TRUE)
  pheatmap::pheatmap(plotmat, cluster_rows=TRUE, cluster_cols=TRUE,
                     border_color=NA, legend_breaks=c(-1,-.5,0,.5,1),
                     main=paste0('Reference versus xenograft GE cluster correlation\n',
                                 attr(m, 'n_shared_genes'), ' shared genes'),
                     angle_col=45)
  invisible(path)
}

hc_s5f_export <- function(root=getwd(), validate_replay=TRUE) {
  if (basename(normalizePath(root)) != 'CH_publication_release')
    warning('Check working directory: expected CH_publication_release root')
  input <- file.path(root, 'human_ch/data/integrated_seurat.rds')
  stopifnot(file.exists(input))
  pub <- file.path(root, 'human_ch/publication')
  dir.create(pub, recursive=TRUE, showWarnings=FALSE)
  paths <- file.path(pub, c(
    'source_data/FigS5F_cluster_correlations.tsv.gz',
    'source_data/FigS5F_cluster_cell_counts.tsv',
    'source_data/FigS5F_shared_gene_centroids.tsv.gz',
    'figures/FigS5F_reference_xenograft_correlations_REVIEW.pdf',
    'validation/S5F_checkpoint_correlation_check.tsv'))
  if (any(file.exists(paths))) stop('S5F refuses to overwrite: ', paste(paths[file.exists(paths)], collapse=', '))
  for (d in unique(dirname(paths))) dir.create(d, recursive=TRUE, showWarnings=FALSE)
  obj <- readRDS(input)
  a <- hc_s5f_source(obj)
  rm(obj); invisible(gc())
  stopifnot(nrow(a$correlations) == nrow(a$matrix)*ncol(a$matrix),
            nrow(a$gene_centroids) == unique(a$correlations$n_shared_genes))
  stage <- tempfile('hc_s5f_', tmpdir=pub)
  dir.create(stage)
  on.exit(unlink(stage, recursive=TRUE), add=TRUE)
  staged <- file.path(stage, basename(paths))
  hc_write_tsv(a$correlations, staged[1]); hc_write_tsv(a$cell_counts, staged[2]); hc_write_tsv(a$gene_centroids, staged[3])
  attr(a$matrix, 'n_shared_genes') <- nrow(a$gene_centroids)
  hc_s5f_pdf(a$matrix, staged[4])
  report <- data.frame(table=c('correlations','cell_counts','gene_centroids'),
                       rows=c(nrow(a$correlations), nrow(a$cell_counts), nrow(a$gene_centroids)),
                       validation='CHECKPOINT_CORRELATION_INTERNAL_PASS',
                       note='S5F reconstructed centroid means; historical centroid function unavailable',
                       stringsAsFactors=FALSE)
  hc_write_tsv(report, staged[5])
  stopifnot(all(file.exists(staged)), all(file.info(staged)$size > 0))
  # Refuse if any output was created concurrently. No overwriting.
  if (any(file.exists(paths))) stop('S5F output appeared during run; refusing promotion')
  for (i in seq_along(staged)) {
    if (!file.rename(staged[i], paths[i])) stop('Unable to promote ', staged[i], ' to ', paths[i])
  }
  print(report, row.names=FALSE)
  cat('S5F source data + PDF saved; HISTORICAL VISUAL REVIEW PENDING.\n')
  invisible(list(status='S5F_CHECKPOINT_INTERNAL_PASS; HISTORICAL_VISUAL_REVIEW_PENDING',
                 paths=paths, n_genes=nrow(a$gene_centroids),
                 n_xeno=nrow(a$matrix), n_reference=ncol(a$matrix)))
}
# === END human_ch S5F checkpoint extension ===

# --- Fig3C: PCA-ranked HSC-to-GMP pseudotime (historical checkpoint reconstruction) ---
# Based on human_ch/scripts/3_lineage_enrichment.R; no Seurat reintegration.
# Historical steps retained: pcaMethods::pca(nPcs=3), rank(PC1), reference
# quantile bins [lower, upper) including the final upper-exclusive endpoint,
# and sample/genotype densities normalized by HSC cluster counts.

hc_3c_dependencies <- function() {
  for (pkg in c('pcaMethods', 'ggplot2'))
    if (!requireNamespace(pkg, quietly=TRUE)) stop('Missing R package: ', pkg)
  if (!'Matrix' %in% loadedNamespaces() ||
      getNamespaceVersion('Matrix') != package_version('1.6-4') ||
      normalizePath(getNamespaceInfo('Matrix', 'path')) != normalizePath(file.path(Sys.getenv("CH_R_LIB", unset = "/packages"), "Matrix")))
    stop('Start a fresh R session and load /packages/Matrix 1.6-4 first.')
}

hc_3c_get_xeno <- function(root) {
  e <- new.env(parent=baseenv())
  load(file.path(root, 'RNA_analysis/data/processed/scrna.RData'), envir=e)
  must <- c('p2pca', 'sampleid', 'clust.final', 'emb', 'col.genotype')
  if (!all(must %in% ls(e))) stop('RNA scrna.RData missing: ', paste(setdiff(must,ls(e)), collapse=', '))
  hash <- readRDS(file.path(root, 'RNA_analysis/data/processed/hashtag_meta.rds'))
  sid <- e$sampleid
  if (is.null(names(sid)) || anyDuplicated(names(sid))) stop('Xenograft sampleid requires unique cell names')
  ix <- match(as.character(sid), rownames(hash))
  if (anyNA(ix)) {
    unmatched <- unique(as.character(sid[is.na(ix)]))
    message('Fig3C xenograft: ', sum(is.na(ix)),
            ' unmatched sample assignments (', length(unmatched),
            ' distinct values; examples: ',
            paste(utils::head(unmatched, 8L), collapse=', '),
            '). Kept for historical PCA; no genotype/experiment assigned.')
  }
  if (all(is.na(ix))) stop('Fig3C: no sample IDs match hashtag metadata; investigate ID format')
  g <- as.character(hash$genotype[ix]); names(g) <- names(sid)
  ex <- as.character(hash$exp[ix]); names(ex) <- names(sid)
  cl <- e$clust.final[names(sid)]
  list(pca=e$p2pca, emb=e$emb, sampleid=as.character(sid),
       cellid=names(sid), genotype=g, exp=ex, cluster=cl,
       palette=e$col.genotype)
}

hc_3c_get_human <- function(root) {
  obj <- readRDS(file.path(root, 'human_ch/data/human_ch_seurat.rds'))
  emb <- obj@reductions$pca@cell.embeddings
  um <- obj@reductions$umap.pca@cell.embeddings
  sample <- obj$Sample
  geno <- obj$Genotype
  cl <- obj$Cluster
  cells <- names(sample)
  if (is.null(cells)) cells <- colnames(obj)
  if (length(cells) != length(sample) || anyDuplicated(cells)) stop('Human Seurat cell identifiers invalid')
  # Historical 3_lineage_enrichment.R, human branch: restrict the initial
  # sampleid vector BEFORE PCA to these four genotypes (not just plot groups).
  eligible <- c('DNMT3A','TET2','WT-CH','WT-control')
  keep <- !is.na(geno) & as.character(geno) %in% eligible
  message('Fig3C human: ', sum(keep), ' cells in historical four-genotype universe; ',
          sum(!keep), ' excluded before PCA.')
  if (!all(c('DNMT3A','TET2','WT-CH','WT-control') %in% unique(as.character(geno[keep]))))
    stop('Fig3C human: expected genotype absent; review input labels')
  sample <- sample[keep]; geno <- geno[keep]; cl <- cl[keep]; cells <- cells[keep]
  names(sample) <- names(geno) <- names(cl) <- cells
  sid <- paste(sample, geno, sep='_'); names(sid) <- cells
  list(pca=emb, emb=um, sampleid=sid, cellid=cells,
       genotype=as.character(geno), exp=rep(NA_character_, length(cells)),
       cluster=cl,
       palette=NULL)
}

hc_3c_prepare <- function(raw, dataset, joint) {
  stopifnot(dataset %in% c('xeno', 'human'))
  ids <- raw$cellid
  for (v in c('sampleid','genotype','exp','cluster'))
    if (length(raw[[v]]) != length(ids)) stop(dataset, ' invalid length: ', v)
  if (anyDuplicated(ids)) stop('Duplicate cell IDs in ', dataset)
  sid <- raw$sampleid; names(sid) <- ids
  gn <- raw$genotype; names(gn) <- ids
  ex <- raw$exp; names(ex) <- ids
  cl <- raw$cluster; names(cl) <- ids
  states <- if (dataset == 'xeno') c(8,7,6,1) else c('HSC/MPP','LMPP','Early GMP','Late GMP')
  hsc_state <- if (dataset == 'xeno') '1' else 'HSC/MPP'
  pc <- as.matrix(raw$pca)
  if (is.null(rownames(pc)) || ncol(pc) < 3L) stop(dataset, ': missing PCA cell IDs or PCs')
  if (anyDuplicated(rownames(pc))) stop(dataset, ': duplicate PCA IDs')
  chosen <- ids[!is.na(cl) & as.character(cl) %in% as.character(states)]
  chosen <- chosen[chosen %in% rownames(pc)]
  if (length(chosen) < 20L) stop(dataset, ': too few selected pseudotime cells')
  # Historical input order is names(clustid), not PCA row order.
  pr <- pcaMethods::pca(pc[chosen, , drop=FALSE], nPcs=3)
  sc <- pr@scores[,1]
  if (length(sc) != length(chosen)) stop(dataset, ': PCA score count mismatch')
  if (is.null(names(sc))) names(sc) <- chosen
  if (!setequal(names(sc),chosen)) stop(dataset, ': PCA score identifier mismatch')
  sc <- sc[chosen]
  ranks <- rank(sc)
  # The source script uses integrated UMAP when do.joint_emb=TRUE.
  um <- as.matrix(raw$emb)
  if (is.null(rownames(um)) || ncol(um) < 2L) stop(dataset, ': no original UMAP')
  emb_ids <- rownames(joint)[rownames(joint) %in% rownames(um)]
  um <- as.matrix(joint[emb_ids, c('umaprpca_1','umaprpca_2'), drop=FALSE])
  colnames(um) <- c('UMAP1','UMAP2')
  hsc_ids <- ids[!is.na(cl) & as.character(cl) == hsc_state]
  hsc_n <- table(sid[hsc_ids])
  base <- data.frame(cell_id=chosen, score_PC1=as.numeric(sc), pseudotime_rank=as.numeric(ranks),
                     sample_id=as.character(sid[chosen]), genotype=as.character(gn[chosen]),
                     experiment=as.character(ex[chosen]), cluster=as.character(cl[chosen]),
                     stringsAsFactors=FALSE)
  if (anyNA(base$sample_id) || anyNA(base$genotype)) {
    if (dataset != 'xeno') stop(dataset, ': NA sample/genotype')
    message('Fig3C xeno PCA selection: ',
            sum(is.na(base$sample_id) | is.na(base$genotype)),
            ' cells have incomplete sample/genotype annotation; ',
            'comparison-specific experiment filtering follows historical script.')
  }
  list(cells=base, hsc_n=hsc_n, um=um, palette=raw$palette)
}

hc_3c_density <- function(prepared, dataset, mutation) {
  stopifnot(dataset %in% c('xeno','human'), mutation %in% c('DNMT3A','TET2'))
  d <- prepared$cells
  if (dataset == 'xeno') {
    exp_no <- if (mutation == 'DNMT3A') '123' else '127'
    d <- d[!is.na(d$experiment) & d$experiment == exp_no, , drop=FALSE]
    if (anyNA(d$sample_id) || anyNA(d$genotype))
      stop('Fig3C xeno ', mutation, ': retained experiment cells have missing sample/genotype IDs')
    ref <- 'WT'; step <- 0.02
  } else {
    other_mut <- if (mutation == 'DNMT3A') 'TET2' else 'DNMT3A'
    d <- d[d$genotype != other_mut, , drop=FALSE]
    ref <- 'WT-control'; step <- 0.1
  }
  if (nrow(d) == 0L || !any(d$genotype %in% ref)) stop(dataset, '/', mutation, ': no reference cells')
  q <- as.numeric(stats::quantile(d$pseudotime_rank[d$genotype %in% ref],
                                  probs=seq(0,1,by=step), names=FALSE))
  if (anyDuplicated(q)) stop(dataset, '/', mutation, ': duplicate quantile edges; review required')
  # findInterval yields the same historical [lower,upper) bins; endpoint excluded.
  bin <- findInterval(d$pseudotime_rank, q, rightmost.closed=FALSE)
  valid <- bin >= 1L & bin < length(q) & d$pseudotime_rank < q[length(q)]
  d$bin <- ifelse(valid, bin, NA_integer_)
  # Historical code stores a sample row only for a nonzero sample-bin count.
  use <- d[!is.na(d$bin), , drop=FALSE]
  sb <- stats::aggregate(rep.int(1L,nrow(use)), by=list(bin=use$bin, sid=use$sample_id), FUN=sum)
  names(sb)[3] <- 'stat'
  gmap <- unique(d[,c('sample_id','genotype')]); names(gmap) <- c('sid','genotype')
  if (anyDuplicated(gmap$sid)) stop(dataset, '/', mutation, ': sample maps to multiple genotypes')
  sb$genotype <- as.character(gmap$genotype[match(sb$sid,gmap$sid)])
  # Original cluster label: modal HSC-to-GMP cluster across all cells in bin.
  top <- vapply(sort(unique(use$bin)), function(k) {
    v <- table(use$cluster[use$bin == k]); names(v)[which.max(v)]
  }, character(1))
  sb$cluster_mode <- unname(top[as.character(sb$bin)])
  sb$hsc_count <- as.integer(prepared$hsc_n[sb$sid])
  sb$stat_norm <- sb$stat / sb$hsc_count
  gp <- stats::aggregate(stat ~ genotype + bin, data=sb, FUN=sum)
  sample_ids <- unique(sb$sid)
  denom <- data.frame(sid=sample_ids,
                      genotype=gmap$genotype[match(sample_ids,gmap$sid)],
                      hsc_count=as.integer(prepared$hsc_n[sample_ids]), stringsAsFactors=FALSE)
  gn_total <- stats::aggregate(hsc_count ~ genotype, data=denom, FUN=sum)
  gp$genotype_hsc_count <- gn_total$hsc_count[match(gp$genotype,gn_total$genotype)]
  gp$stat_norm <- gp$stat / gp$genotype_hsc_count
  cuts <- data.frame(bin=seq_len(length(q)-1L), lower=q[-length(q)], upper=q[-1L],
                     reference_genotype=ref, step=step, stringsAsFactors=FALSE)
  panel <- paste(dataset, mutation, sep='_')
  for (obj in c('d','sb','gp','denom','cuts')) {
    z <- get(obj)
    z$panel <- panel
    assign(obj,z)
  }
  d <- d[order(d$cell_id), , drop=FALSE]
  sb <- sb[order(sb$bin, sb$sid), , drop=FALSE]
  gp <- gp[order(gp$bin, gp$genotype), , drop=FALSE]
  denom <- denom[order(denom$sid), , drop=FALSE]
  # Validate sample-level/bin and genotype-level aggregation independently.
  direct <- as.data.frame(table(use$bin, use$sample_id, dnn=c('bin','sid')), stringsAsFactors=FALSE)
  direct <- direct[direct$Freq > 0, ,drop=FALSE]
  key <- paste(sb$bin, sb$sid)
  idx <- match(key,paste(direct$bin,direct$sid))
  stopifnot(!anyNA(idx), identical(as.integer(sb$stat),as.integer(direct$Freq[idx])))
  bygeno <- stats::aggregate(stat ~ bin + genotype, data=sb, FUN=sum)
  key2 <- paste(gp$bin,gp$genotype)
  idx2 <- match(key2,paste(bygeno$bin,bygeno$genotype))
  stopifnot(!anyNA(idx2), all(gp$stat == bygeno$stat[idx2]), sum(sb$stat) == nrow(use))
  list(cells=d, sample_bin=sb, genotype_bin=gp, sample_HSC=denom,
       reference_bins=cuts, background=prepared$um,
       all_ranks=setNames(prepared$cells$pseudotime_rank,prepared$cells$cell_id),
       palette=prepared$palette)
}

hc_3c_plot <- function(z, filename, dataset, mutation) {
  dc <- z$cells; sb <- z$sample_bin; gp <- z$genotype_bin
  bg <- z$background
  ids <- intersect(rownames(bg), names(z$all_ranks))
  d <- data.frame(UMAP1=as.numeric(bg[ids,1]), UMAP2=as.numeric(bg[ids,2]),
                  rank=unname(z$all_ranks[ids]))
  # Match original Fig3C script: first 10k integrated cells used as grey background.
  grey <- data.frame(UMAP1=as.numeric(bg[seq_len(min(10000L,nrow(bg))),1]),
                     UMAP2=as.numeric(bg[seq_len(min(10000L,nrow(bg))),2]))
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data=grey, ggplot2::aes(x=UMAP1,y=UMAP2),
                        size=0.05,alpha=0.1,color='grey20') +
    ggplot2::geom_point(data=d, ggplot2::aes(x=UMAP1,y=UMAP2,color=rank),size=0.05) +
    ggplot2::scale_color_gradient2(low='blue',mid='white',high='red',
                                    midpoint=mean(d$rank)) +
    ggplot2::theme_bw() + ggplot2::labs(x='',y='') +
    ggplot2::theme(legend.position='none')
  p2 <- ggplot2::ggplot(sb,ggplot2::aes(x=bin,y=stat_norm,group=sid,color=genotype)) +
    ggplot2::geom_point(alpha=0.3, na.rm=TRUE) +
    ggplot2::geom_line(alpha=0.3,na.rm=TRUE) +
    ggplot2::geom_smooth(span=0.3,linewidth=0.5,se=FALSE,na.rm=TRUE) +
    ggplot2::geom_smooth(data=gp,ggplot2::aes(x=bin,y=stat_norm,group=genotype,color=genotype),
                         inherit.aes=FALSE,span=if(dataset=='xeno') 0.3 else 1,
                         linewidth=3,se=FALSE,na.rm=TRUE) +
    ggplot2::labs(x='pseudotime',y='density to HSC') + ggplot2::theme_bw() +
    ggplot2::theme(legend.position=c(0.1,0.65),legend.justification=c(0,0))
  palette <- z$palette
  if (!is.null(palette) && length(palette)) {
    palette <- as.character(palette)
    names(palette) <- names(z$palette)
    if (dataset == 'human') names(palette)[names(palette)=='WT'] <- 'WT-CH'
    palette['WT-CH'] <- 'lightblue'
    palette['WT-control'] <- 'black'
    present <- unique(as.character(sb$genotype))
    if (all(present %in% names(palette)))
      p2 <- p2 + ggplot2::scale_color_manual(values=palette)
  }
  grDevices::pdf(filename,width=8,height=4,onefile=TRUE)
  on.exit(grDevices::dev.off(),add=TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,2)))
  print(p1,vp=grid::viewport(layout.pos.row=1,layout.pos.col=1),newpage=FALSE)
  print(p2,vp=grid::viewport(layout.pos.row=1,layout.pos.col=2),newpage=FALSE)
  invisible(filename)
}

hc_fig3c_export <- function(root=getwd(),datasets=c('xeno','human')) {
  hc_3c_dependencies()
  stopifnot(all(datasets %in% c('xeno','human')),length(datasets)>0L)
  datasets <- unique(datasets)
  pub <- file.path(root,'human_ch/publication')
  sd <- file.path(pub,'source_data')
  figs <- file.path(pub,'figures')
  val <- file.path(pub,'validation')
  for (dr in c(sd,figs,val)) dir.create(dr,recursive=TRUE,showWarnings=FALSE)
  prefix <- 'Fig3C_'
  targets <- c('source_data/Fig3C_pseudotime_cells.tsv.gz',
               'source_data/Fig3C_sample_bin_density.tsv.gz',
               'source_data/Fig3C_genotype_bin_density.tsv.gz',
               'source_data/Fig3C_HSC_sample_denominators.tsv',
               'source_data/Fig3C_reference_quantile_bins.tsv',
               'validation/Fig3C_internal_validation.tsv',
               as.vector(outer(datasets,c('DNMT3A','TET2'),
                 function(ds,gn) paste0('figures/Fig3C_',ds,'_',gn,'_REVIEW.pdf'))))
  # Full four-panel run by default. Partial runs must not collide with completed outputs.
  existing <- file.path(pub,targets)[file.exists(file.path(pub,targets))]
  if (length(existing)) stop('Refusing to overwrite: ',paste(existing,collapse=', '))
  stage <- tempfile('Fig3C_stage_',tmpdir=pub); dir.create(stage)
  on.exit(unlink(stage,recursive=TRUE),add=TRUE)
  for (sub in c('source_data','figures','validation')) dir.create(file.path(stage,sub))
  all <- list(); checks <- list(); joint <- readRDS(file.path(root,'human_ch/data/integrated_seurat_metadata.rds'))
  historical_palette <- NULL
  for (ds in datasets) {
    message('Fig3C: loading ',ds,' checkpoint; fitting historical pcaMethods PCA once')
    raw <- if(ds=='xeno') hc_3c_get_xeno(root) else hc_3c_get_human(root)
    if (ds == 'xeno') historical_palette <- raw$palette
    if (ds == 'human' && !is.null(historical_palette)) raw$palette <- historical_palette
    prep <- hc_3c_prepare(raw,ds,joint)
    rm(raw); gc(verbose=FALSE)
    for (gn in c('DNMT3A','TET2')) {
      z <- hc_3c_density(prep,ds,gn)
      panel <- paste(ds,gn,sep='_'); all[[panel]] <- z
      fn <- file.path(stage,paste0('figures/Fig3C_',ds,'_',gn,'_REVIEW.pdf'))
      hc_3c_plot(z,fn,ds,gn)
      checks[[panel]] <- data.frame(panel=panel,selected_cells=nrow(z$cells),
         valid_binned_cells=sum(!is.na(z$cells$bin)),sample_bin_rows=nrow(z$sample_bin),
         genotype_bin_rows=nrow(z$genotype_bin),
         status='CHECKPOINT_DERIVATION_INTERNAL_PASS',stringsAsFactors=FALSE)
      message('Fig3C ',panel,': ',nrow(z$cells),' cells; ',nrow(z$sample_bin),' sample bins')
      # Drop plotting-only large matrices after PDF is generated.
      all[[panel]]$background <- NULL
      all[[panel]]$all_ranks <- NULL
      all[[panel]]$palette <- NULL
    }
    rm(prep); gc(verbose=FALSE)
  }
  fields <- list(cells='Fig3C_pseudotime_cells.tsv.gz',
                 sample_bin='Fig3C_sample_bin_density.tsv.gz',
                 genotype_bin='Fig3C_genotype_bin_density.tsv.gz',
                 sample_HSC='Fig3C_HSC_sample_denominators.tsv',
                 reference_bins='Fig3C_reference_quantile_bins.tsv')
  for (field in names(fields)) {
    z <- do.call(rbind,lapply(all,`[[`,field)); rownames(z) <- NULL
    hc_write_tsv(z,file.path(stage,'source_data',fields[[field]]))
  }
  validation <- do.call(rbind,checks); rownames(validation)<-NULL
  hc_write_tsv(validation,file.path(stage,'validation/Fig3C_internal_validation.tsv'))
  if (!all(file.exists(file.path(stage,targets))) ||
      any(file.info(file.path(stage,targets))$size <= 0))
    stop('Fig3C: incomplete output; no promotion')
  for (p in targets) {
    if (!file.copy(file.path(stage,p),file.path(pub,p),overwrite=FALSE))
      stop('Fig3C output promotion failed: ',p)
  }
  print(validation,row.names=FALSE)
  cat('Fig3C: Source Data + ',2L*length(datasets),' PDFs generated; historical visual review PENDING.\n')
  invisible(list(status='CHECKPOINT_DERIVATION_INTERNAL_PASS; HISTORICAL_VISUAL_REVIEW_PENDING',
                 report=validation, paths=file.path(pub,targets)))
}

# --- Fig3C: render human PDFs from saved Source Data (palette only) ---
# Does not refit PCA, recalculate bins, or alter journal Source Data.
hc_fig3c_replot_palette <- function(root=getwd(), output_dir=NULL) {
  if (is.null(output_dir))
    output_dir <- file.path(root,'human_ch/publication/validation/Fig3C_palette_review')
  sd <- file.path(root,'human_ch/publication/source_data')
  read_sd <- function(nm) {
    f <- file.path(sd,nm)
    if (!file.exists(f)) stop('Missing approved Source Data: ',f)
    utils::read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)
  }
  cells <- read_sd('Fig3C_pseudotime_cells.tsv.gz')
  bins <- read_sd('Fig3C_sample_bin_density.tsv.gz')
  gbin <- read_sd('Fig3C_genotype_bin_density.tsv.gz')
  stopifnot(all(c('panel','cell_id','pseudotime_rank') %in% names(cells)),
            all(c('panel','bin','sid','stat_norm','genotype') %in% names(bins)),
            all(c('panel','bin','stat_norm','genotype') %in% names(gbin)))
  meta <- readRDS(file.path(root,'human_ch/data/integrated_seurat_metadata.rds'))
  stopifnot(!is.null(rownames(meta)),
            all(c('dataset','umaprpca_1','umaprpca_2') %in% colnames(meta)))
  human <- meta[is.na(meta$dataset),c('umaprpca_1','umaprpca_2'),drop=FALSE]
  if (!nrow(human) || anyNA(human)) stop('Human integrated coordinates incomplete')
  um <- as.matrix(human); colnames(um) <- c('UMAP1','UMAP2')
  ranks <- cells[grepl('^human_',cells$panel),c('cell_id','pseudotime_rank'),drop=FALSE]
  if (!nrow(ranks)) stop('No human pseudotime cell records')
  if (anyDuplicated(ranks$cell_id)) {
    dup <- duplicated(ranks$cell_id)
    first <- match(ranks$cell_id[dup],ranks$cell_id)
    if (any(ranks$pseudotime_rank[dup] != ranks$pseudotime_rank[first]))
      stop('Human cell pseudotime ranks disagree across panels')
    ranks <- ranks[!duplicated(ranks$cell_id),,drop=FALSE]
  }
  all_ranks <- setNames(ranks$pseudotime_rank,ranks$cell_id)
  if (!all(names(all_ranks) %in% rownames(um)))
    stop('Human pseudotime cells missing from integrated UMAP')
  if (dir.exists(output_dir) && length(list.files(output_dir,all.files=TRUE,no..=TRUE)))
    stop('Nonempty review destination: ',output_dir)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  stopifnot(dir.exists(output_dir))
  outs <- character()
  for (mutation in c('DNMT3A','TET2')) {
    panel <- paste0('human_',mutation)
    d <- cells[cells$panel==panel,,drop=FALSE]
    b <- bins[bins$panel==panel,,drop=FALSE]
    g <- gbin[gbin$panel==panel,,drop=FALSE]
    if (!nrow(d) || !nrow(b) || !nrow(g)) stop('No saved data for ',panel)
    if (!setequal(unique(as.character(b$genotype)),c(mutation,'WT-CH','WT-control')))
      stop('Expected three-group saved data absent for ',panel)
    z <- list(cells=d,sample_bin=b,genotype_bin=g,
              background=um,all_ranks=all_ranks,
              palette=c('DNMT3A'='red','TET2'='green',
                        'WT-CH'='lightblue','WT-control'='black'))
    fn <- file.path(output_dir,paste0('Fig3C_human_',mutation,'_PALETTE_REVIEW.pdf'))
    hc_3c_plot(z,fn,'human',mutation)
    if (!file.exists(fn) || file.info(fn)$size < 1000)
      stop('Nonempty PDF not generated: ',fn)
    outs <- c(outs,fn)
  }
  cat('PASS: two human Fig3C PDFs rendered from saved numerical Source Data; PCA not rerun.\n')
  cat('No active PDFs, Source Data, or scientific inputs overwritten.\n')
  print(outs)
  invisible(outs)
}

