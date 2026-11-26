# Fig3B: nearest-neighbour smoothing from portable cell-level Source Data.
# Unlike original O(N^2) distance loops this uses RANN::nn2 (exact kd-tree search).
# Main analysis settings preserved: query own genotype/experiment cells,
# 29 nearest neighbors (historical rank(distance) < 30), palette index
# as.integer(log2(fold_change) * 50) + 51, clipped to [1, 100].
if (!requireNamespace('RANN', quietly=TRUE)) {
  stop('Fig3B smoothing needs CRAN package RANN. Install with install.packages("RANN")')
}
sd_dir <- 'RNA_analysis/publication/source_data'
val_dir <- 'RNA_analysis/publication/validation'
a <- read.delim(gzfile(file.path(sd_dir,'Fig3A_umap_cells.tsv.gz')),
                check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('','NA'))
f <- read.delim(file.path(sd_dir,'FigS6B_state_fold_changes.tsv'), stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(a$cell_id), all(c('123','127') %in% a$experiment),
          all(c('123','127') %in% f$experiment))

results <- list()
for (expid in c('123','127')) {
  geno <- if (expid == '123') 'DNMT3A' else 'TET2'
  q <- a[!is.na(a$cell_state) & !is.na(a$experiment) &
           a$experiment == expid & a$genotype == geno,,drop=FALSE]
  stopifnot(nrow(q) >= 29L)
  fold <- f[f$experiment == expid,]
  vals <- fold$fold_change[match(q$cell_state,fold$cell_state)]
  if(anyNA(vals) || any(vals <= 0)) stop('Missing/invalid cell-state fold changes for ',expid)
  encoded <- as.integer(log2(vals)*50) + 51L
  encoded <- pmax(1L,pmin(100L,encoded))
  xy <- as.matrix(q[,c('umap_1','umap_2')])
  # For continuous coordinates, nearest 29 is equivalent to rank(y) < 30.
  # Exact distance ties can cause differences; documented in README.
  knn <- RANN::nn2(xy,query=xy,k=29,searchtype='standard')
  smoothed <- rowMeans(matrix(encoded[as.vector(knn$nn.idx)],
                               nrow=nrow(q),ncol=29L))
  q$cell_state_fold_change <- vals
  q$encoded_state_index <- encoded
  q$smoothed_index <- smoothed
  q$nearest_neighbor_k <- 29L
  results[[expid]] <- q[,c('cell_id','experiment','genotype','cell_state',
                          'umap_1','umap_2','cell_state_fold_change',
                          'encoded_state_index','smoothed_index',
                          'nearest_neighbor_k')]
  cat('Fig3B:',geno,nrow(q),'cells, range(smoothed index)=',
      paste(round(range(smoothed),2),collapse=' to '),'\n')
}
b <- do.call(rbind,results)
stopifnot(all(b$smoothed_index >= 1),all(b$smoothed_index <= 100),
          all(is.finite(b$smoothed_index)))
write.table(b, gzfile(file.path(sd_dir,'Fig3B_smoothed_enrichment.tsv.gz')),
            sep='\t',row.names=FALSE,quote=FALSE,na='')
write.table(data.frame(experiment=c('123','127'),
                       genotype=c('DNMT3A','TET2'),
                       n_cells=vapply(results,nrow,integer(1)),
                       smoothing_k=29L,
                       stringsAsFactors=FALSE),
            file.path(val_dir,'Fig3B_export_summary.tsv'),
            sep='\t',row.names=FALSE,quote=FALSE)
cat('Exported Fig3B Source Data. Check original picture for nearest-neighbor tie/order effects.\n')
