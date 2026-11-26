# Publication-compatible reproduction of RNA_analysis/R/project_gset.R.
# Keeps the historical math: mean expression -> PCA -> linear fitted values,
# and first-PC loadings. Makes matrix subsetting and method dispatch explicit.
# The original exploratory helper remains unchanged.

project_gset_publication <- function(gset, mat, pr = NULL) {
  if (length(dim(mat)) != 2L || is.null(rownames(mat)) ||
      is.null(colnames(mat))) {
    stop('Expression input must have two dimensions and named rows/columns.')
  }
  genes <- intersect(as.character(gset), colnames(mat))
  if (!length(genes)) return(NULL)

  expr <- mat[, genes, drop = FALSE]
  if (length(dim(expr)) != 2L) stop('Gene subsetting did not retain 2 dimensions.')
  is_expressed <- Matrix::colSums(expr > 0) > 0
  genes <- genes[is_expressed]
  if (!length(genes)) return(NULL)
  if (length(genes) < 2L && is.null(pr)) {
    stop('Only one expressed program gene remains; original PCA requires >=2.')
  }
  expr <- mat[, genes, drop = FALSE]
  av <- Matrix::rowMeans(expr)

  if (is.null(pr)) {
    # prcomp() works on a dense subset of only the program genes, not the
    # entire expression matrix. This matches the historical PCA calculation.
    expr_small <- as.matrix(expr)
    if (any(!is.finite(expr_small))) stop('Nonfinite program expression values.')
    if (any(apply(expr_small, 2L, stats::sd) == 0)) {
      stop('A program gene has zero variance; original scaled PCA cannot run.')
    }
    pca <- stats::prcomp(expr_small, center = TRUE, scale. = TRUE)
    gscore <- pca$rotation[, 1L]
    pr <- pca$x
  } else {
    gscore <- NULL
  }

  fit <- stats::lm(av ~ pr)
  list(av = av, av.fit = fit$fitted.values,
       gscore = gscore, n_matched_expressed_genes = length(genes))
}
