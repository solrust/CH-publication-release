# select high-low clone phenotypes

env.dna <- new.env()
load('../DNA_analysis/data/processed/DNA_clones.RData', env = env.dna)

# TODO:
# assemble matrix of clone phenotypes
# choose 'group1' and 'group2' cells
# choose phenotypes & put them in the matrix of clone.info.multi.
# run through phenotypes

#
run.old <- FALSE
size.cut <- 50
med <- median(clone.info.multi[clone.info.multi$size > size.cut, '8'], na.rm = T)
clones.high <- clone.info.multi$larry[which(clone.info.multi[, '8'] > med & clone.info.multi$size > size.cut)]
clones.low <- clone.info.multi$larry[which(clone.info.multi[, '8'] < med & clone.info.multi$size > size.cut)]
#
do.rate <- F; do.permute <- F
if (do.rate){
  x <- intersect(names(env.dna$rate), clone.info.multi$larry)
  rate <- env.dna$rate[x]
  if (do.permute){
    rate <- sample(rate)
    names(rate) <- names(env.dna$rate[x])
  }
  med <- median(rate)
  clones.high <- clone.info.multi$larry[which(rate[clone.info.multi$larry] > med & clone.info.multi$size > size.cut)]
  clones.low <- clone.info.multi$larry[which(rate[clone.info.multi$larry] < med & clone.info.multi$size > size.cut)]
}
#
if (!run.old){
  cells.high <- rownames(cell.clone.multi)[cell.clone.multi$larry %in% clones.high]
  cells.low <- rownames(cell.clone.multi)[cell.clone.multi$larry %in% clones.low]
}else{
  cells.high <- rownames(cell.clone)[cell.clone$cloneid %in% clones.high]
  cells.low <- rownames(cell.clone)[cell.clone$cloneid %in% clones.low]
}
str(cells.high)
str(cells.low)

# sample/cluster/clone_state -specific profiles
# expr/acc ~ clone_state + cluster (+ genotype/batch) (+ genotype:clone_state)

# build list of two matrices: for clone-high and clone-low.
# aggregaet them as


DimPlot(obj.combined, reduction = "wnn.umap", shuffle = T, group.by = 'projection', #'seurat_clusters',
        label.box = T, label = T, label.size = 3,  repel = F)
par(mfrow=c(1,1))
plot(obj.combined@reductions$wnn.umap@cell.embeddings)
points(obj.combined@reductions$wnn.umap@cell.embeddings[obj.combined$seurat_clusters %in% c(0,1,5),], col = 'red')


# estimate clone phenotype-specific count matrices
Idents(obj.combined) <- "seurat_clusters"
if (!run.old){
  #p2clust1 <- obj.combined$projection
  p2clust1 <- obj.combined$seurat_clusters
  #count.mat <- obj.combined@assays$RNA$counts
  count.mat <- obj.combined@assays$ATAC$counts
}else{
  p2clust1 <- old.scrna$clust.final
  count.mat <- old.scrna$cd
}
n.clust <- length(unique(p2clust1))-1
clusts.prog <- 1#c(1, 6, 2)
clusts <- c(1, 6, 2)
clusts <- c(4,7,8)
#clusts <- 1:10
clusts <- c(0, 1, 5)
clust.clone.list <- lapply(1:2, function(i){
  if (i == 1){
    cells <- cells.low
  }else{
    cells <- cells.high
  }
  matl <- lapply(clusts, function(clusts.prog){
    cls <- names(p2clust1)[!is.na(p2clust1) & p2clust1 %in% clusts.prog]
    cls <- intersect(cls, cells)
    print(paste('cells:',length(cls)))

    run.expgeno <- TRUE
    if (!run.expgeno){
      smpls <- unique(cell.clone.multi[cls,]$id)
      if (run.old) smpls <- unique(cell.clone[cls,]$id)
    }else{
      smpls <- unique(cell.clone.multi[cls,]$expgeno)
    }
    smpls <- smpls[!is.na(smpls)]
    mat <- do.call(cbind, lapply(smpls, function(smpl){
      #print(smpl)
      if (run.old){
        ind <- cell.clone[cls,]$id %in% smpl
      }else{
        if (!run.expgeno){
          ind <- cell.clone.multi[cls,]$id %in% smpl
        }else{
          ind <- cell.clone.multi[cls,]$expgeno %in% smpl

        }
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
### merge clone samples
#mat <- do.call(cbind, clust.clone.list)
#clid <- unlist(lapply(1:length(clust.clone.list), function(i) rep(i, ncol(clust.clone.list[[i]])) ))
#hash.match1 <- hash.match[colnames(mat), , drop = F]
#hash.match1$clid <- clid
mat <- do.call(cbind, do.call(cbind, clust.clone.list))
hash.match1 <- do.call(rbind, lapply(1:length(clust.clone.list), function(i){
  do.call(rbind, lapply(1:length(clust.clone.list[[i]]), function(j){
    rep(i, names(clust.clone.list[[i]])[j], colnames(clust.clone.list[[i]]))
    cbind(colnames(clust.clone.list[[i]][[j]]), i, names(clust.clone.list[[i]])[j])
  }))
}))
hash.match1 <- data.frame(hash.match1)
colnames(hash.match1) <- c('batch', 'clone', 'clust')
colnames(mat) <- rownames(hash.match1)
ind1 <- T#hash.match1$exp == '123' #& hash.match1$clid %in% c(1, 2)#, 6, 7)
if (length(unique(hash.match1$clust)) > 1){
  dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                    design = ~ clone + clust)#clid + exp)
}else{
  dds_mut <- DESeqDataSetFromMatrix(countData = mat[,ind1], colData = hash.match1[ind1,],
                                    design = ~ clone)#clid + exp)
}
dds_mut <- DESeq(dds_mut)
resultsNames(dds_mut)
res_clone <- lfcShrink(dds_mut, coef="clone_2_vs_1", type="normal")
qplot(res_clone$pvalue)
sum(res_clone$padj < 0.1, na.rm = T); sum(res_clone$padj < 0.5, na.rm = T)
lbl <- rownames(res_clone)
lbl[is.na(res_clone$pvalue) | res_clone$pvalue > 1e-4] <- NA
table(is.na(lbl))
qplot(res_clone$log2FoldChange, res_clone$pvalue, log = 'y', label = lbl)+geom_text_repel()

qplot(res_clone$log2FoldChange, res_clone1$log2FoldChange, label = lbl) + geom_text_repel()
cor(res_clone$log2FoldChange, res_clone3$log2FoldChange, use = 'na.or.complete')

# res_clone1 - hsc cluster
# res_clone2 - prog cluters without cluster control
# res_clone3 - prog clusters with cluster control

#rownames(da_commitment[['1']])(da_commitment[['1']]$log2FoldChange > 0)
#res_clone <- res_clone1
pa.scores <- res_clone$log2FoldChange
names(pa.scores) <- rownames(res_clone)
pa.scores <- pa.scores[!is.na(pa.scores) & res_clone$pvalue < 1e-1]
#pa.scores <- pa.scores[!is.na(pa.scores) & abs(res_clone$log2FoldChange) > 0.2]
#pa.scores <- pa.scores[!is.na(pa.scores) & da_commitment[['1']]$log2FoldChange > 0.0]
str(pa.scores)
#
polars <- table(pa.scores > 0)
polars
#
pks.sets <- lapply(1:ncol(motif.matrix), function(i){
  pks <- names(which(motif.matrix[, i]))
  pks <- intersect(names(pa.scores), pks)
  pks
})
names(pks.sets) <- colnames(motif.matrix)
qplot(unlist(lapply(pks.sets, length)))
fgseaRes <- fgsea::fgsea(pathways = pks.sets, pa.scores,
                  minSize=20, maxSize = 30000)#, nperm=1000)
head(fgseaRes[,1:7])
rownames(fgseaRes) <- unlist(motif.object@motif.names)[fgseaRes$pathway]
qplot(fgseaRes$NES, fgseaRes$pval, log = 'y', label = rownames(fgseaRes)) +
  geom_text_repel()
#
gg <- intersect(rownames(fgseaRes), rownames(fgseaRes1))
str(gg)
qplot(fgseaRes[match(gg, rownames(fgseaRes)),]$NES, fgseaRes1[match(gg, rownames(fgseaRes1)),]$NES, label = gg) + geom_text_repel()

open.peaks <- AccessiblePeaks(obj.combined, idents = clusts)
run_peaks <- T
scores <- do.call(rbind,
                  mclapply(1:ncol(motif.matrix), function(i){
                  #lapply(da_commitment_peaks[-c(1)], function(pks){
                    #print(str(pks))
  if (i %% 10 == 0) print(i)
  pks <- names(which(motif.matrix[, i]))
  pks <- intersect(names(pa.scores), pks)
  #
  if (run_peaks){
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[pks, ],
    features.match = c('GC.percent', "percentile"), #c("percentile"),
    n = 50000
  )
  }else{
    peaks.matched <- setdiff(names(pa.scores), pks)
  }
  #
  #pks <- intersect(names(pa.scores), pks)
  #
  polars.tf <- table(as.factor(pa.scores > 0)[pks])
  polars <- table(pa.scores[peaks.matched] > 0)
  #print(polars.tf)
  tst <- fisher.test(cbind(polars.tf, polars))
  scores <- (polars.tf[2] / polars.tf[1]) / (polars[2] / polars[1])
  c(scores, tst$p.value, polars.tf)
}, mc.cores = 20))
#scores
rownames(scores) <- unlist(motif.object@motif.names)[colnames(motif.matrix)]
#
pl1 <- ggplot(,aes(scores[,1], p.adjust(scores[,2], 'BH'),
                   label = rownames(scores))) + #theme_bw() +
  scale_y_log10(limits = c(1e-2, 1e-0)) + #limits = c(1e-7, 1)) + xlim(1, 1.1) +
  geom_point() + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() +
  #labs(x = paste('enrichment in', names(mutant_scores)[1], 'DA peaks'), y = 'p-value') +
  theme(text = element_text(size = 15))
pl1
#
gg <- intersect(rownames(fgseaRes), rownames(scores))
str(gg)
qplot(fgseaRes[match(gg, rownames(fgseaRes)),]$NES, scores[gg,1], label = gg) + geom_text_repel()
#
gg <- intersect(rownames(scores), rownames(scores_small))
qplot(scores[gg,1], scores_small[gg,1], label = gg)+geom_abline()+geom_text_repel()
#
qplot(scores1[gg,1], scores3[gg,1], label = gg)+geom_abline()+geom_text_repel()
qplot(scores_peaks[,1], scores_basic[gg,1], label = rownames(scores_basic))+geom_abline()+geom_text_repel()

ggplot(,aes(scores_small[,1], p.adjust(scores_small[,2], 'BH'),
            label = rownames(scores_small))) + #theme_bw() +
  scale_y_log10()+#limits = c(1e-7, 1e-0)) + #limits = c(1e-7, 1)) + xlim(1, 1.1) +
  geom_point() + geom_vline(xintercept = 1, linetype = 2) + geom_text_repel() +
  #labs(x = paste('enrichment in', names(mutant_scores)[1], 'DA peaks'), y = 'p-value') +
  theme(text = element_text(size = 15))

#
meta.feature <- GetAssayData(obj.combined, assay = "ATAC", slot = "meta.features")
res_mut <- res_clone
open_clust <- rownames(res_mut[!is.na(res_mut$padj) & res_mut$pvalue < 1e-1 &
                                 res_mut$log2FoldChange > 0.0, ])
str(open_clust)
open.peaks <- AccessiblePeaks(obj.combined, idents = 1)
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[open_clust, ],
  #features.match = 'GC.percent', #c("percentile"),
  n = 50000
)
enriched.motifs <- FindMotifs(
  object = obj.combined,
  features = open_clust,
  background = peaks.matched,#peaks.matched,
)
head(enriched.motifs)
pl1 <- MotifPlot(object = obj.combined,
                 motifs = head(rownames(enriched.motifs), 10)
)
pl1
enriched.motifs$motif.name.show <- enriched.motifs$motif.name
enriched.motifs$motif.name.show[enriched.motifs$p.adjust > 1e-2] <- NA
pl2 <- ggplot(enriched.motifs, aes(fold.enrichment, pvalue, label = motif.name.show)) + geom_point(col = old.scrna$clust.final.cols[cl]) + theme_bw() +
  scale_x_log10() + scale_y_log10() + labs(x = 'motif fold enrichment', y = 'p-value', title = paste('cluster', cl)) + theme(text = element_text(size = 20))
pl2 <- pl2 + ggrepel::geom_text_repel()
pl2
#print(grid.arrange(pl2, pl1))
head(enriched.motifs, 10)



