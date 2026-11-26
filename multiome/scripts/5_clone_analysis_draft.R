obj.combined <- readRDS('data/processed/multiome_seurat.rds')

obj.combined <- FindClusters(obj.combined, graph.name = "wsnn",
                             cluster.name = 'main_clusters',
                             algorithm = 3, resolution = 1, verbose = FALSE)
table(obj.combined$main_clusters)

table(obj.combined$seurat_clusters)
table(prev.env$obj.combined$main_clusters)

DimPlot(obj.combined, reduction = "wnn.umap", shuffle = T, group.by = 'main_clusters',#'main_clusters', #'detailed_clusters', #'seurat_clusters',
        label.box = T, label = T, label.size = 3,  repel = F)


size.cut <- 50
motif.matrix <- obj.combined@assays$ATAC@motifs@data
motif.names <- obj.combined@assays$ATAC@motifs@motif.names
motif.counts <- t(motif.matrix) %*% obj.combined@assays$ATAC$counts
#
motif.matrix.pks <- Matrix::Matrix(0, nrow = nrow(motif.matrix), ncol = ncol(motif.matrix))
rownames(motif.matrix.pks) <- rownames(motif.matrix); colnames(motif.matrix.pks) <- colnames(motif.matrix)
pks <- da_commitment_peaks[['10']]
str(pks)
pks <- intersect(rownames(motif.matrix), pks)
str(pks)
motif.matrix.pks[pks,] <- motif.matrix[pks,]
motif.counts.pks <- t(motif.matrix.pks) %*% obj.combined@assays$ATAC$counts
#
motif.gc <- t(motif.matrix) %*% meta.feature$GC.percent
motif.gc <- motif.gc[,1] / colSums(motif.matrix)
#
#rownames(motif.counts) <- unlist(motif.names)[rownames(motif.counts)]
pheno.names <- c('1', '4', '5', '8', '10', 'rate', 'size')
geno.names <- c('DNMT3A', 'TET2')

meta.feature <- GetAssayData(obj.combined, assay = "ATAC", slot = "meta.features")
open.peaks <- AccessiblePeaks(obj.combined, idents = c(0, 1, 5))

# generate GC-controlled peak permutation
gc.groups <- as.numeric(cut(meta.feature$GC.percent, 30))
qplot(gc.groups, meta.feature$GC.percent)
gcl <- lapply(1:length(unique(gc.groups)), function(i) which(gc.groups == i))
gc.groups.control <- vector()
for (i in 1:length(gc.groups)){
  if (i %% 1000 == 0) print(i)
  j <- gc.groups[i]
  k <- sample(gcl[[j]], 1)
  gc.groups.control[i] <- k
  gcl[[j]] <- setdiff(gcl[[j]], k)
}
# estimate GC-controled background motif-peak matrix
motif.matrix.gc <- motif.matrix[gc.groups.control,]
motif.matrix.gc <- Matrix::Matrix(motif.matrix.gc, sparse = TRUE)
colnames(motif.matrix.gc) <- colnames(motif.matrix)
rownames(motif.matrix.gc) <- rownames(motif.matrix)
#
motif.gc.control <- t(motif.matrix.gc) %*% meta.feature$GC.percent
motif.gc.control <- motif.gc.control[,1] / colSums(motif.matrix.gc)
# estimate GC-controlled ATAC counts
motif.counts.gc <- t(motif.matrix.gc) %*% obj.combined@assays$ATAC$counts
#
qplot(motif.gc.control, motif.gc)+geom_abline()
i <- 1
ggplot() + geom_density(aes(meta.feature$GC.percent[motif.matrix[,i] == 1], col = 'real')) +
           geom_density(aes(meta.feature$GC.percent[motif.matrix.gc[,i] == 1], col = 'control'))

ind <- which(unlist(motif.names) == 'CEBPD')
ind
inds <- rownames(motif.matrix)[motif.matrix[, ind] == 1]
str(inds)
x <- unlist(lapply(da_commitment_peaks, function(x) length(intersect(x, inds))))
y <- unlist(lapply(da_commitment_peaks, length))
x
y
x/y
x/length(inds)
da_commitment_peaks[['10']]


run.geno <- F
use.hspc <- T
#
clusts <- 1
if (use.hspc) clusts <- c(0, 1, 5)
#
if (run.geno){
  names <- geno.names
}else names <- pheno.names
#
pop <- 'HSC'
if (use.hspc) pop <- 'HSPC'
#
test.list <- list()
run_deseq2 <- TRUE#; test.list.use <- mutant.test.hsc
for (i in 1:length(names)){
  i <- 5
  pheno.name <- names[i]
  print(pheno.name)

  # identify high/low clones
  if (!run.geno){
    batch_name <- 'expgeno'
    batch_name <- 'id'
    pheno.score <- clone.info.multi[, pheno.name]
    #pheno.score <- clone.info.multi[, '10'] - clone.info.multi[, '8']
    med <- median(pheno.score[clone.info.multi$size > size.cut], na.rm = T)
    clones.high <- clone.info.multi$larry[which(pheno.score > med & clone.info.multi$size > size.cut)]
    clones.low <- clone.info.multi$larry[which(pheno.score < med & clone.info.multi$size > size.cut)]
    cells.high <- rownames(cell.clone.multi)[cell.clone.multi$larry %in% clones.high]
    cells.low <- rownames(cell.clone.multi)[cell.clone.multi$larry %in% clones.low]
  }else{
    batch_name <- 'id'
    if (pheno.name == 'DNMT3A'){
      expid <- '123'
      geno <- 'DNMT3A'
    }else{
      expid <- '127'
      geno <- 'TET2'
    }
    cells.high <- rownames(cell.clone.multi)[cell.clone.multi$genotype %in% geno & grepl(expid, cell.clone.multi$id)]
    cells.low <- rownames(cell.clone.multi)[cell.clone.multi$genotype %in% 'WT' & grepl(expid, cell.clone.multi$id)]
  }
  #
  # permute:
  #cells <- c(cells.high, cells.low)
  #cells.high <- sample(cells, length(cells.high))
  #cells.low <- setdiff(cells, cells.high)
  #

  #cells <- rownames(cell.clone.multi)[cell.clone.multi$larry %in% c(clones.high, clones.low)]
  #cells1 <- cells[obj.combined$seurat_clusters[cells] == 0]
  #cells2 <- cells[obj.combined$seurat_clusters[cells] == 5]
  #cells3 <- cells[obj.combined$seurat_clusters[cells] == 1]
  #cells1 <- cells1[!is.na(cells1)]; cells2 <- cells2[!is.na(cells2)]; cells3 <- cells3[!is.na(cells3)]
  #cells.high <- c(sample(cells1, length(cells1) * 0.7),
  #                sample(cells2, length(cells2) * 0.3),
  #                sample(cells3, length(cells3) * 0.5))
  #cells.low <- setdiff(c(cells1, cells2, cells3), cells.high)
  #

  #print(str(cells.high))
  #print(str(cells.low))

  if (run_deseq2){
    cell_conditions <- list(cells.low, cells.high)
    da_analysis <- perform_da(obj.combined,
                              #modality = 'ATAC',
                              count.mat = motif.counts,
                              #count.mat = commitment.counts,
                              do.deseq2 = T,
                              cell_conditions = cell_conditions,
                              cell.info = cell.clone.multi,
                              #batch_name = 'expgeno',
                              #batch_name = 'id',
                              batch_name = batch_name,
                              seurat_class = 'main_clusters', clusts = clusts
                              #seurat_class = 'main_clusters', clusts = c(1)
                              )

    da_analysis0 <- perform_da(obj.combined,
                              modality = 'ATAC',
                              #count.mat = motif.counts,
                              #count.mat = commitment.counts,
                              do.deseq2 = F,
                              cell_conditions = cell_conditions,
                              cell.info = cell.clone.multi,
                              #batch_name = 'expgeno',
                              #batch_name = 'id',
                              batch_name = batch_name,
                              seurat_class = 'main_clusters', clusts = clusts
                              #seurat_class = 'main_clusters', clusts = c(1)
    )
  }else{
    da_analysis <-  test.list.use[[i]]$da_analysis
  }

  #da_analysis <- list(res = res_atac_dnmt3a)
  motif_stat <- estimate_motif_enrichment(obj.combined, da_analysis, motif.matrix, motif.names,
                                          #pcuts = c(0.05, 0.15)
                                          pcuts = c(0.2, 0.4),
                                          #features.match = c('GC.percent'),
                                          run_peaks = F,
                                          permute = F
                                          )


  pll <- make_plots(motif_stat)
  #print(grid.arrange(grobs = pll, nrow = 3, ncol = 3))

  if (run.geno){
    flname <- paste('figures/geno_', pop ,'_', pheno.name, '.pdf', sep = '')
  }else{
    flname <- paste('figures/clone_', pop ,'_', pheno.name, '.pdf', sep = '')
  }
  pdf(flname, width = 15, height = 15)
  print(grid.arrange(grobs = pll, nrow = 3, ncol = 3))
  dev.off()

  test.list[[i]] <- list(cell_conditions = cell_conditions, da_analysis = da_analysis, motif_stat = motif_stat)
}
if (run.geno){
  if (use.hspc){
    mutant.test.hspc <- test.list
  }else{
    mutant.test.hsc <- test.list
  }
}else{
  if (use.hspc){
    clone.test.hspc <- test.list
  }else{
    clone.test.hsc <- test.list
  }
}

#save(mutant.test.hsc, file = 'data/processed/mutant_test_hsc.RData')
#save(mutant.test.hspc, file = 'data/processed/mutant_test_hspc.RData')
#save(clone.test.hspc, file = 'data/processed/clone_test_hspc.RData')
#save(clone.test.hsc, file = 'data/processed/clone_test_hsc.RData')
#saveRDS(clone.test, 'data/processed/clone_test.rds')

pheno.names
i <- 5
x1 <- clone.test.hsc[[i]]$motif_stat$chisql[[1]]
x2 <- clone.test.hspc[[i]]$motif_stat$chisql[[1]]
#x1 <- mutant.test.hspc[[2]]$motif_stat$chisql[[1]]
#x2 <- mutant.test.hsc[[2]]$motif_stat$chisql[[1]]
qplot(x1[,1], x2[,1], label = rownames(x1))+geom_text_repel()+geom_abline()
cor(x1[,1], x2[,1])
pheno.names[i]
#
x1 <- do.call(cbind, lapply(clone.test.hsc, function(x) x$motif_stat$chisql[[1]][,1]))
x2 <- do.call(cbind, lapply(clone.test.hspc, function(x) x$motif_stat$chisql[[1]][,1]))
colnames(x1) <- colnames(x2) <- pheno.names
pheatmap::pheatmap(cor(x1, x1), cluster_rows = F, cluster_cols = F)

i <- 1
pmat <- do.call(cbind, lapply(clone.test.hspc, function(x) x$motif_stat$chisql[[i]][,2]))
#pmat <- cbind(pmat,
#              do.call(cbind, lapply(mutant.test.new, function(x) x$motif_stat$chisql[[1]][,2])))
pmat <- apply(pmat, 2, function(x) p.adjust(x, 'BH'))
colnames(pmat) <- c(pheno.names)#, geno.names)
lfmat <- do.call(cbind, lapply(clone.test.hspc, function(x) x$motif_stat$chisql[[i]][,1]))
#lfmat <- cbind(lfmat,
#               do.call(cbind, lapply(mutant.test.new, function(x) x$motif_stat$chisql[[1]][,1])))
colnames(lfmat) <- c(pheno.names)#, geno.names)
#
pheno.fix <- c('1', '4', '5', '8', '10', 'rate')#, 'size')#, 'DNMT3A', 'TET2')
ind <- apply(pmat[,pheno.fix], 1, min) < 1e-4
pheat <- pheatmap::pheatmap(log(lfmat)[ind, pheno.fix], breaks = seq(-0.4, 0.4, length.out = 100),
                   col = colorRampPalette(c('blue','grey90', 'grey99', 'red'))(100))
pdf('figures/clone_heatmap.pdf', width = 5, height = 16)
pheat
dev.off()

cor(log(lfmat)[ind, pheno.fix])

i <- 6
x1 <- clone.test.hsc[[i]]$motif_stat$chisql[[2]][,1]
x2 <- clone.test.new[[i]]$motif_stat$chisql[[1]][,1]
qplot(x1, x2, label = names(x1), main = pheno.names[i]) + geom_text_repel()
cor(x1, x2)

### plot mutant trends
res1 <- mutant.test.hspc[[1]]$motif_stat$chisql[[1]]
res2 <- mutant.test.hspc[[2]]$motif_stat$chisql[[1]]
pl1 <- ggplot(,aes(res1[,1], res1[,2],
                   label = rownames(res1))) + theme_bw() + scale_y_log10() +
  geom_point() + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() +
  labs(x = paste('enrichment in', 'DNMT3A', 'DA peaks'), y = 'p-value') + theme(text = element_text(size = 15))
#
pl2 <- ggplot(,aes(res2[,1], res2[,2],
                   label = rownames(res2))) + theme_bw() + scale_y_log10() +
  geom_point() + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() +
  labs(x = paste('enrichment in', 'TET2', 'DA peaks'), y = 'p-value') + theme(text = element_text(size = 15))
#
#pl3 <-
pl3 <- ggplot(,aes(res1[,1], res2[,1], label = rownames(res1))) + theme_bw() + geom_point() +
  geom_text_repel() + geom_vline(xintercept = 1, linetype = 2) + geom_hline(yintercept = 1, linetype = 2) +
  scale_x_log10() + scale_y_log10() +
  labs(x = 'DNMT3A', y = 'TET2', title = 'motif enrichment') + theme(text = element_text(size = 20))

pdf('figures/DA_mutant_motifs.pdf', width = 10, height = 10)
#pl1 + pl2
grid.arrange(pl1, pl2, pl3, ncol = 2, nrow = 2)
dev.off()

#
pl1 <- ggplot(,aes(mutant.test.hspc[[1]]$motif_stat$chisql[[1]][,1], mutant.test.hsc[[1]]$motif_stat$chisql[[1]][,1],
              label = rownames(mutant.test.hspc[[1]]$motif_stat$chisql[[1]]))) +
  geom_point() + geom_abline(linetype = 2) + geom_text_repel() + labs(x = 'HSPC', y = 'HSC', title = 'DNMT3A')# + theme_bw()
pl2 <- ggplot(,aes(mutant.test.hspc[[2]]$motif_stat$chisql[[1]][,1], mutant.test.hsc[[2]]$motif_stat$chisql[[1]][,1],
                   label = rownames(mutant.test.hspc[[2]]$motif_stat$chisql[[1]]))) +
  geom_point() + geom_abline(linetype = 2) + geom_text_repel() + labs(x = 'HSPC', y = 'HSC', title = 'TET2')# + theme_bw()
pdf('figures/DA_mutant_motifs_hsc_hspc.pdf', width = 10, height = 5)
pl1 + pl2
dev.off()

# plot clone phenotypes:
pl <- lapply(1:length(clone.test.hspc), function(i){
  x <- clone.test.hspc[[i]]
  res <- x$motif_stat$chisql[[1]]
  ggplot(,aes(res[,1], res[,2],
              label = rownames(res))) + theme_bw() + scale_y_log10() +
    geom_point() + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() +
    labs(x = paste('enrichment in', pheno.names[i], 'DA peaks'), y = 'p-value') + theme(text = element_text(size = 15))
})

pdf('figures/DA_clone_motifs.pdf', width = 15, height = 15)
print(grid.arrange(grobs = pl, nrow = 3, ncol = 3))
dev.off()

### analyse enrichment of fate-specific DA peaks
clust.deriv <- c(1, 4, 8, 5, 9, 10, 3)
clust.prog <- c(1, 2, 6)
da_commitment <- mclapply(clust.deriv, function(cl){
  #cl <- 4
  mat1 <- clust.atac.list[[cl]]
  #mat2 <- clust.atac.list[[10]]
  mat2 <- Reduce('+', clust.atac.list[clust.prog])
  mat <- cbind(mat1, mat2)
  inf <- cbind(hash.match[colnames(mat), , drop = F],
               cond = c(rep('clust', ncol(mat1)), rep('back', ncol(mat2))) )
  colnames(mat) <- rownames(inf)
  dds_mut <- DESeqDataSetFromMatrix(countData = mat,
                                    colData = inf,
                                    design = ~ cond + exp)

  dds_mut <- DESeq(dds_mut)
  resultsNames(dds_mut)
  res_mut <- lfcShrink(dds_mut, coef="cond_clust_vs_back",
                       #type="apeglm"
                       type="normal"
  )
  return(res_mut)
}, mc.cores = length(clust.deriv))
names(da_commitment) <- clust.deriv
da_commitment_max <- apply(do.call(cbind, lapply(da_commitment, function(x) x$log2FoldChange)), 1, max)
da_commitment_peaks <- lapply(da_commitment, function(x){
  #rownames(x)[!is.na(x$padj) & x$padj < 1e-1 &  x$log2FoldChange > 1]
  rownames(x)[!is.na(x$padj) & x$padj < 1e-2 &  x$log2FoldChange >= 1 &
                x$log2FoldChange >= da_commitment_max]
})
unlist(lapply(da_commitment_peaks, length))
commitment.counts <- do.call(rbind, lapply(da_commitment_peaks[-1], function(x){
  x <- intersect(x, rownames(obj.combined@assays$ATAC$counts))
  colSums(obj.combined@assays$ATAC$counts[x,])
}))
da_commitment_peaks1 <- lapply(1:length(da_commitment_peaks), function(i){
  pks <- unlist(da_commitment_peaks[-i])
  setdiff(da_commitment_peaks[[i]], pks)
})
unlist(lapply(da_commitment_peaks1, length))

test.list <- clone.test.hspc
pcut <- 0.1
run_peaks <- F
meta.feature <- GetAssayData(obj.combined, assay = "ATAC", slot = "meta.features")
Idents(obj.combined) <- 'main_clusters'; open.peaks <- AccessiblePeaks(obj.combined, idents = c(0, 1, 5))
mut_bias_scores <- lapply(1:length(test.list), function(igen){
  res <- test.list[[igen]]$da_analysis$res

  pa.scores <- res$log2FoldChange
  names(pa.scores) <- rownames(res)
  #pa.scores <- pa.scores[!is.na(pa.scores)]
  pa.scores <- pa.scores[!is.na(pa.scores) & rank(res$pvalue) <= pcut * sum(!is.na(res$pvalue)) ]

  #
  polars <- table(pa.scores > 0)
  fate_scores <- do.call(rbind, lapply(da_commitment_peaks[-1], function(pks){
    pks <- intersect(names(pa.scores), pks)
    #polars.tf <- table(as.factor(pa.scores > 0)[pks])
    #print(polars.tf)
    #polars <- polars - polars.tf

    peaks.matched <- setdiff(names(pa.scores), pks)
    if (run_peaks){
      peaks.matched <- MatchRegionStats(
        meta.feature = meta.feature[open.peaks,],
        query.feature = meta.feature[pks, ],
        features.match = 'GC.percent', n = 50000)
    }
    polars.tf <- table(as.factor(pa.scores > 0)[pks])
    polars <- table(pa.scores[peaks.matched] > 0)
    tst <- fisher.test(cbind(polars.tf, polars))


    tst <- fisher.test(cbind(polars.tf, polars))
    scores <- (polars.tf[2] / polars.tf[1]) / (polars[2] / polars[1])
    c(scores, tst$p.value, polars.tf)
  }))
  rownames(fate_scores) <- names(da_commitment_peaks[-1])
  #qplot(fate_scores[,1], fate_scores[,2], log = 'y', label = rownames(fate_scores))
  return(fate_scores)
})

#mut_bias_scores <- mut_bias_scores1
i <- 5
pl1 <- ggplot(,aes(mut_bias_scores[[i]][,1], mut_bias_scores[[i]][,2],
                   label = rownames(mut_bias_scores[[i]]))) + geom_point(size = I(3), alpha = I(0.5)) +
  geom_vline(xintercept = 1, linetype = 2) + geom_hline(yintercept = 1e-2, linetype = 2) +
  theme_bw() + geom_text_repel() + scale_x_log10() + scale_y_log10() +
  labs(x = 'enrichment of fate-specific peaks', y = 'p-value',
       )#title = paste(names(fate_mut_scores)[i], 'DA peaks'))
#
j <- 6
pl2 <- ggplot(,aes(mut_bias_scores[[j]][,1], mut_bias_scores[[j]][,2],
                   label = rownames(mut_bias_scores[[j]]))) + geom_point(size = I(3), alpha = I(0.5)) +
  geom_vline(xintercept = 1, linetype = 2) + geom_hline(yintercept = 1e-2, linetype = 2) +
  theme_bw() + geom_text_repel() + scale_x_log10() + scale_y_log10() +
  labs(x = 'enrichment of fate-specific peaks', y = 'p-value',
       )#title = paste(names(fate_mut_scores)[j], 'DA peaks'))
pl1 + pl2


x <- do.call(cbind, lapply(mut_bias_scores, function(x) x[,1]))
#colnames(x) <- geno.names
colnames(x)  <- c('hsc bias', 'lymph bias', 'dend bias', 'my bias', 'ery bias', 'rate', 'size')
pheat <- pheatmap::pheatmap(log(t(x)), breaks = seq(-1, 1, length.out = 100),
                   main = 'fate-specific peaks in clonal DA peaks')
                   #col = colorRampPalette(c('blue', 'grey90', 'grey95', 'darkred'))(100) )
pdf('figures/DA_clone_fate.pdf')#, width = 15, height = 15)
pheat
dev.off()

#################### functions
make_plots <- function(motif_stat){
  p1 <- qplot(motif_stat$fgseal[[1]]$NES, motif_stat$fgseal[[1]]$padj, log = 'y', label = rownames(motif_stat$fgseal[[1]])) + theme_bw() +
    labs(x = 'NES', y = 'p-adj', title = 'GSEA, all peaks') + geom_text_repel() # + theme(text = element_text(size = 20))
  p2 <- qplot(motif_stat$fgseal[[2]]$NES, motif_stat$fgseal[[2]]$padj, log = 'y', label = rownames(motif_stat$fgseal[[2]])) + theme_bw() +
    labs(x = 'NES', y = 'p-adj', title = 'GSEA, peaks p<=0.1') + geom_text_repel() # + theme(text = element_text(size = 20))
  gg <- intersect(rownames(motif_stat$fgseal[[1]]), rownames(motif_stat$fgseal[[2]]))
  p3 <- qplot(motif_stat$fgseal[[1]][match(gg, rownames(motif_stat$fgseal[[1]])),]$NES, motif_stat$fgseal[[2]][match(gg, rownames(motif_stat$fgseal[[2]])),]$NES, label = gg) + theme_bw() +
    labs(x = 'NES, all', y = 'NES, p <= 0.1') + geom_text_repel() # + theme(text = element_text(size = 20))
  #
  p4 <-  qplot(motif_stat$chisql[[1]][,1], p.adjust(motif_stat$chisql[[1]][,2], 'BH'), log = 'y', label = rownames(motif_stat$chisql[[1]])) + theme_bw() +
    labs(x = 'enrichment', y = 'p-adj', title = 'peaks all') + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() # + theme(text = element_text(size = 20))
  p5 <-  qplot(motif_stat$chisql[[2]][,1], p.adjust(motif_stat$chisql[[2]][,2], 'BH'), log = 'y', label = rownames(motif_stat$chisql[[2]])) + theme_bw() +
    labs(x = 'enrichment', y = 'p-adj', title = 'peaks p<=0.1') + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() # + theme(text = element_text(size = 20))
  p6 <- qplot(motif_stat$chisql[[1]][,1], motif_stat$chisql[[2]][,1], label = rownames(motif_stat$chisql[[1]])) + theme_bw() +
    labs(x = 'enrichment, all', y = 'entichment, p<=0.1') + geom_text_repel() # + theme(text = element_text(size = 20))
  #
  gg <- intersect(rownames(motif_stat$fgseal[[1]]), rownames(motif_stat$chisql[[1]]))
  p7 <- qplot(motif_stat$fgseal[[1]][match(gg, rownames(motif_stat$fgseal[[1]])),]$NES, motif_stat$chisql[[1]][gg,1], label = gg) + theme_bw() +
    labs(x = 'GSEA', y = 'entichment', title = 'all peaks') + geom_text_repel() # + theme(text = element_text(size = 20))
  gg <- intersect(rownames(motif_stat$fgseal[[2]]), rownames(motif_stat$chisql[[2]]))
  p8 <- qplot(motif_stat$fgseal[[2]][match(gg, rownames(motif_stat$fgseal[[2]])),]$NES, motif_stat$chisql[[2]][gg,1], label = gg) + theme_bw() +
    labs(x = 'GSEA', y = 'entichment', title = 'peaks p<=0.1') + geom_text_repel() # + theme(text = element_text(size = 20))
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, p7 = p7, p8 = p8))
}

perform_da <- function(obj.combined, cell_conditions = NULL,
                             modality = NULL,  count.mat = NULL,
                             cell.info = NULL, batch_name = NULL,
                            do.deseq2 = TRUE,
                             seurat_class = 'seurat_clusters', clusts = NULL, old.dataset = FALSE){
  if (is.null(clusts)) stop('provide clusters')
  if (!old.dataset){
    if (!is.null(modality)){
      if (modality == 'ATAC'){
        count.mat <- obj.combined@assays$ATAC$counts
      }else if (modality == 'RNA'){
        count.mat <- obj.combined@assays$RNA$counts
      }
    }
    Idents(obj.combined) <- seurat_class
    p2clust1 <- Idents(obj.combined)
  }else{
    p2clust1 <- old.scrna$clust.final
    count.mat <- old.scrna$cd
  }

  print(paste('count.mat:', dim(count.mat)))
  # calculate pseudo-bulk condition-specific matrices
  print('calculate pseudo-bulk matrices')
  clust.clone.list <- lapply(1:2, function(i){
    cells <- cell_conditions[[i]]
    matl <- lapply(clusts, function(clust){
      cls <- names(p2clust1)[!is.na(p2clust1) & p2clust1 %in% clust]
      cls <- intersect(cls, cells)
      #print(paste('cells:',length(cls)))
      smpls <- unique(cell.info[cls, batch_name])
      if (old.dataset) smpls <- unique(cell.info[cls,]$id)
      smpls <- smpls[!is.na(smpls)]

      mat <- do.call(cbind, lapply(smpls, function(smpl){
        #print(smpl)
        if (old.dataset){
          ind <- cell.info[cls,]$id %in% smpl
        }else{
          ind <- cell.info[cls, batch_name] %in% smpl
        }
        if (sum(ind) > 0){
          rowSums(count.mat[, cls][, ind, drop = F])
        }else{
          rep(0, nrow(count.mat))
        }
      }))
      colnames(mat) <- smpls

      return(mat)
    })
    names(matl) <- clusts
    return(matl)
  })

  # combine in a single matrix and record metadata
  mat <- do.call(cbind, do.call(cbind, clust.clone.list))
  #
  hash.match1 <- do.call(rbind, lapply(1:length(clust.clone.list), function(i){
    do.call(rbind, lapply(1:length(clust.clone.list[[i]]), function(j){
      cbind(colnames(clust.clone.list[[i]][[j]]), i, names(clust.clone.list[[i]])[j])
    }))
  }))
  hash.match1 <- data.frame(hash.match1)
  colnames(hash.match1) <- c('batch', 'clone', 'clust')
  hash.match1$clone <- as.factor(hash.match1$clone)
  hash.match1$clust <- as.factor(hash.match1$clust)
  colnames(mat) <- rownames(hash.match1)

  ind <- colSums(mat) > 0
  mat <- mat[, ind]
  hash.match1 <- hash.match1[ind,]

  # perform DESeq2 analysis
  if (do.deseq2){
  print('perform DESeq2 analysis')
  ind1 <- TRUE
  print(hash.match1[ind1,])
  print(dim(mat[,ind1]))
  print(length(length(unique(hash.match1$clust))))
  print(hash.match1$clone)
  if (length(unique(hash.match1$clust)) > 1){
    dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                      design = ~ clone + clust)
  }else{
    dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                      design = ~ clone)
  }
  dds_mut <- DESeq(dds_mut)
  print(resultsNames(dds_mut))
  #res <- lfcShrink(dds_mut, coef="clone_2_vs_1", type="normal")
  res <-  results(dds_mut, name="clone_2_vs_1")#, type="normal")
  }else{
    dds_mut <- NULL; res <- NULL
  }

  return(list(mat = mat, hash.match1 = hash.match1, dds = dds_mut, res = res))
}


estimate_motif_enrichment <- function(obj.combined, da_analysis,  motif.matrix = NULL, motif.names = NULL, pcuts = c(1, 1e-1),
                                      permute = FALSE, run_peaks = TRUE,
                                      features.match = c('GC.percent', "percentile"), run_fgsea = FALSE){
  res <- da_analysis$res
  # run fgsea for motif enrichment
  if (run_fgsea){
  print('perform fgsea for motif enrichment')
  fgseal <- lapply(pcuts, function(pcut){
    print(pcut)
    pa.scores <- res$log2FoldChange
    if (permute) pa.scores <- sample(pa.scores)
    names(pa.scores) <- rownames(res)
    pa.scores <- pa.scores[!is.na(pa.scores) & rank(res$pvalue) <= pcut * sum(!is.na(res$pvalue)) ]
    pks.sets <- lapply(1:ncol(motif.matrix), function(i){
      pks <- names(which(motif.matrix[, i]))
      pks <- intersect(names(pa.scores), pks)
      pks
    })
    names(pks.sets) <- colnames(motif.matrix)
    fgseaRes <- fgsea::fgsea(pathways = pks.sets, pa.scores,
                             minSize=20, maxSize = 30000)#, nperm=1000)
    rownames(fgseaRes) <- unlist(motif.names)[fgseaRes$pathway]

    return(fgseaRes)
  })
  }else{
    fgseal <- NULL
  }

  # run chisq test for motif enrichment
  print('perform fisher test for motif enrichment')
  meta.feature <- GetAssayData(obj.combined, assay = "ATAC", slot = "meta.features")
  open.peaks <- AccessiblePeaks(obj.combined, idents = clusts)
  # TODO: use all 20 cores?
  chisql <- lapply(pcuts, function(pcut){
    scores <- do.call(rbind, mclapply(1:ncol(motif.matrix), function(i){
      if (i %% 100 == 0) print(i)
      pa.scores <- res$log2FoldChange
      if (permute) pa.scores <- sample(pa.scores)
      names(pa.scores) <- rownames(res)
      pa.scores <- pa.scores[!is.na(pa.scores) &  rank(res$pvalue) <= pcut * sum(!is.na(res$pvalue)) ]#res$pvalue <= pcut]
      #
      pks <- names(which(motif.matrix[, i]))
      pks <- intersect(names(pa.scores), pks)
      peaks.matched <- setdiff(names(pa.scores), pks)
      if (run_peaks){
        peaks.matched <- MatchRegionStats(
          meta.feature = meta.feature[open.peaks,],
          query.feature = meta.feature[pks, ],
          features.match = features.match, n = 50000)
      }
      polars.tf <- table(as.factor(pa.scores > 0)[pks])
      polars <- table(pa.scores[peaks.matched] > 0)
      tst <- fisher.test(cbind(polars.tf, polars))
      scores <- (polars.tf[2] / polars.tf[1]) / (polars[2] / polars[1])
      c(scores, tst$p.value, polars.tf, polars)
    }, mc.cores = ifelse(run_peaks, 20, 1)))
    rownames(scores) <- unlist(motif.names)[colnames(motif.matrix)]
    return(scores)
  })


  return(list(fgseal = fgseal, chisql = chisql))
}


