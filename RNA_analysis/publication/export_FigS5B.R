# Source Data for Supplementary Fig. 5B: mitochondrial/doublet QC
# From historical p2.rds and hashtag annotations, without rerunning Pagoda2.
# Run from CH_publication_release root in a fresh 64 GB (or larger) R session.
# Key improvement over historical 3_sc_processing.R: sparse correlation,
# never as.matrix(p2$counts) for all genes.

if (!requireNamespace('Matrix', quietly=TRUE)) stop('Matrix package is required')
input_p2 <- 'RNA_analysis/data/p2objects/p2.rds'
input_hash <- 'RNA_analysis/data/processed/hashtag_info.rds'
sdir <- 'RNA_analysis/publication/source_data'
vdir <- 'RNA_analysis/publication/validation'
dir.create(sdir, recursive=TRUE, showWarnings=FALSE)
dir.create(vdir, recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(input_p2), file.exists(input_hash))

message('Loading p2.rds (large; allow 64 GB session).')
if (requireNamespace('pagoda2', quietly=TRUE)) suppressPackageStartupMessages(library(pagoda2))
p2 <- readRDS(input_p2)
mat <- p2$counts
emb <- p2$embeddings$PCA$UMAP
clust <- p2$clusters$PCA$leiden
if (!inherits(mat, 'dgCMatrix')) {
  stop('Expected p2$counts as a sparse dgCMatrix; observed ', paste(class(mat),collapse='/'),
       '. Do not convert the full matrix to dense; investigate first.')
}
stopifnot(!is.null(rownames(mat)), !is.null(colnames(mat)),
          !is.null(rownames(emb)), !is.null(names(clust)),
          identical(rownames(mat), rownames(emb)),
          setequal(rownames(emb), names(clust)))
clust <- clust[rownames(emb)]
n <- nrow(mat)
message('Cells: ', n, '; genes: ', ncol(mat), '; Leiden clusters: ', length(unique(clust)))

# Reconstruct the historical sampleid0 labels exactly as in 3_sc_processing.R.
hashtags <- readRDS(input_hash)
stopifnot(is.list(hashtags), length(hashtags)==8L)
sampleid0 <- unlist(lapply(seq_along(hashtags),function(i) {
  expid <- if (i<5L) '123' else '127'
  x <- paste(expid, hashtags[[i]], sep='_')
  names(x) <- names(hashtags[[i]])
  x[hashtags[[i]]=='empty'] <- 'empty'
  x[hashtags[[i]]=='doublets'] <- 'doublets'
  x
}))
if (!all(rownames(emb) %in% names(sampleid0))) {
  stop('Some Pagoda2 cells are missing from hashtag_info.rds. Check cell IDs.')
}
is_doublet <- sampleid0[rownames(emb)] == 'doublets'
stopifnot(!anyNA(is_doublet))

# Old script: mitochondrial genes -> PC1 -> correlations with ALL genes.
# The historical as.matrix(p2$counts) would require ~30 GB additional RAM.
mt_genes <- colnames(mat)[grepl('MT-',colnames(mat),fixed=TRUE)]
if (length(mt_genes) < 2L) stop('Fewer than 2 mitochondrial genes in p2$counts.')
message('Initial mitochondrial genes: ', length(mt_genes))
initial_pc <- stats::prcomp(as.matrix(mat[,mt_genes,drop=FALSE]),center=TRUE,scale.=TRUE)
pc <- as.numeric(initial_pc$x[,1])

# Calculate exact Pearson correlations using sparse sums and crossproducts.
# Equivalent to cor(pc, as.matrix(mat)) without creating a dense full matrix.
cs <- Matrix::colSums(mat)
css <- Matrix::colSums(mat^2)
sspc <- sum(pc^2) - sum(pc)^2/n
num <- as.numeric(Matrix::crossprod(mat,pc)) - sum(pc)*cs/n
den <- sqrt(pmax(0,sspc - cs^2/n) * sspc)
mitcor <- num / den
names(mitcor) <- colnames(mat)
mt_genes1 <- unique(c(names(mitcor)[!is.na(mitcor) & abs(mitcor)>0.4], mt_genes))
message('Selected mitochondrial-correlated genes: ', length(mt_genes1))
estimated_dense_GiB <- n*length(mt_genes1)*8/1024^3
if (estimated_dense_GiB>4) stop('Second PCA would allocate at least ',
  round(estimated_dense_GiB,1), ' GiB of dense matrix. Defer S5B and review memory first.')

# Historical mito-PC1 calculation and threshold of 11.
mito_pca <- stats::prcomp(as.matrix(mat[,mt_genes1,drop=FALSE]),center=TRUE,scale.=TRUE)
mito_score <- as.numeric(mito_pca$x[,1])
mt_cutoff <- 11
cluster_levels <- sort(unique(as.integer(as.character(clust))))
cluster_factor <- factor(as.integer(as.character(clust)),levels=cluster_levels)
med <- tapply(mito_score,cluster_factor,stats::median)
mito_flag <- mito_score >= mt_cutoff | as.character(clust) %in% names(med)[med>mt_cutoff]

sizes <- table(cluster_factor)
doubles <- table(factor(as.integer(as.character(clust))[is_doublet], levels=cluster_levels))
frac <- as.numeric(doubles)/as.numeric(sizes)
cluster_qc <- data.frame(cluster=cluster_levels,n_cells=as.integer(sizes),
     n_hashtag_doublets=as.integer(doubles), doublet_fraction=frac,
     median_mito_PC1=as.numeric(med),
     mito_cluster_flag=as.numeric(med)>mt_cutoff,
     doublet_cluster_flag=frac>=0.3)
doublet_flag <- is_doublet | as.character(clust) %in% as.character(cluster_qc$cluster[cluster_qc$doublet_cluster_flag])

cells <- data.frame(cell_id=rownames(emb), UMAP_1=as.numeric(emb[,1]),
  UMAP_2=as.numeric(emb[,2]), leiden_cluster=as.integer(as.character(clust)),
  mito_PC1=mito_score, hashtag_doublet=is_doublet, flagged_mito=mito_flag,
  flagged_doublet=doublet_flag)
stopifnot(nrow(cells)==nrow(emb), all(is.finite(cells$mito_PC1)))
write.table(cells,gzfile(file.path(sdir,'FigS5B_cell_QC.tsv.gz')),
            sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(cluster_qc,file.path(sdir,'FigS5B_cluster_QC.tsv'),
            sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(data.frame(gene=names(mitcor),mito_PC1_correlation=mitcor,
            selected_for_second_PCA=names(mitcor)%in%mt_genes1),
            gzfile(file.path(sdir,'FigS5B_mito_gene_correlations.tsv.gz')),
            sep='\t',row.names=FALSE,quote=FALSE,na='NA')
summary <- data.frame(metric=c('cells','genes','leiden_clusters','mitochondrial_genes',
  'second_PC_genes','mito_score_cutoff','doublet_fraction_cutoff','hashtag_doublets',
  'flagged_mito_cells','flagged_doublet_cells'),
  value=c(n,ncol(mat),length(cluster_levels),length(mt_genes),length(mt_genes1),
          mt_cutoff,0.3,sum(is_doublet),sum(mito_flag),sum(doublet_flag)))
write.table(summary,file.path(vdir,'FigS5B_export_summary.tsv'),
  sep='\t',row.names=FALSE,quote=FALSE)
print(summary,row.names=FALSE)
message('S5B QC tables exported; compare with accepted S5B artwork before approval.')
