library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(glmGamPoi)
library(parallel)
library(GenomicRanges)

hash.match <- readRDS('../RNA_analysis/data/processed/hashtag_meta.rds')
old.scrna <- new.env()
load(file = '../RNA_analysis/data/processed/scrna.RData', envir = old.scrna)
#
cells.old <- names(old.scrna$sampleid)[old.scrna$sampleid != 'doublets']
str(cells.old)
meta.old <- hash.match[old.scrna$sampleid[cells.old],]
rownames(meta.old) <- cells.old
meta.old$clust.final <- NA
meta.old[names(old.scrna$clust.final),]$clust.final <- old.scrna$clust.final
meta.old$sampleid <- old.scrna$sampleid[rownames(meta.old)]

cells <- cells.old
old.hemo <- CreateSeuratObject(counts = old.scrna$cd[, cells],
                               meta.data = meta.old[cells,])
old.hemo@meta.data <- cbind(old.hemo@meta.data,
                            meta.old[cells, c('exp', 'genotype', 'sampleid', 'clust.final')])
old.hemo@meta.data$dataset <- old.hemo@meta.data$exp
old.hemo[["RNA"]] <- split(old.hemo[["RNA"]], f = meta.old[cells,]$exp)


load(file = 'data/processed/multiome_input.RData')
sample.names <- names(inputdatal)
new.hemo <- mclapply(sample.names, function(sample.name){
  rna_counts <- inputdatal[[sample.name]]$`Gene Expression`
  ind <- colSums(rna_counts) >= 500
  hemo <- CreateSeuratObject(counts = rna_counts[, ind])
  hemo$dataset <- sample.name
  hemo
}, mc.cores = length(sample.names))
names(new.hemo) <- sample.names
#
new.hemo <- merge(
  x = new.hemo[[1]],
  y = new.hemo[-1],
  add.cell.ids = sample.names
)

rna.combined <- merge(
  x = old.hemo,
  y = new.hemo,
)
rna.combined

rna.combined <- SCTransform(rna.combined, verbose = FALSE, vst.flavor = "v2") %>% RunPCA()
rna.combined <- IntegrateLayers(object = rna.combined, method = HarmonyIntegration,
                                orig.reduction = "pca", new.reduction = "harmony",
                                assay = "SCT", verbose = T)
rna.combined[["RNA"]] <- JoinLayers(rna.combined[["RNA"]])
#
rna.combined <- FindNeighbors(rna.combined, graph.name = 'graph.harmony', reduction = "harmony", dims = 1:30)
rna.combined <- RunUMAP(rna.combined, dims = 1:30, reduction = "harmony", reduction.name = 'umap.harmony', reduction.key = 'harmonyUMAP_')
rna.combined <- FindClusters(rna.combined, graph.name = 'graph.harmony', cluster.name = 'harmony.clusters',
                             algorithm = 3, resolution = 0.5, verbose = FALSE)
#
rna.combined$dataset[rna.combined$dataset == '123'] <- 'exp1'
rna.combined$dataset[rna.combined$dataset == '127'] <- 'exp2'
progs <- readRDS("../RNA_analysis/data/genesets/processed/fate_markers.rds")
rna.combined <- CellCycleScoring(rna.combined, s.features = progs[["S phase"]], g2m.features = progs[["G2M"]], set.ident = TRUE)
DimPlot(rna.combined, reduction = "umap.harmony", shuffle = T,
        order = rev(c('exp1', 'exp2', 'C1', 'C2', 'D', 'T')), ncol = 2,
        group.by = c("dataset", 'harmony.clusters', 'clust.final', 'Phase'))

#saveRDS(rna.combined, file = 'data/processed/rna_integration.rds')

# mask multiome-specific clusters,  representing removed 'bad' cells
harmony.clust.mask
ind <- rna.combined$dataset %in% c('exp1', 'exp2')
clust.abund.old <- table(rna.combined$harmony.clusters[ind])
clust.abund <- table(rna.combined$harmony.clusters)
rt <- clust.abund.old / clust.abund[names(clust.abund.old)]
sort(rt)
harmony.clust.mask <- names(rt)[rt < 0.5]
harmony.clust.mask

# project new scRNA-seq cells onto reference clusters
cells.ref <- rownames(rna.combined@meta.data)[rna.combined@meta.data$dataset %in% c('exp1', 'exp2')]
str(cells.ref)
cells.ref <- sample(cells.ref, 30000)
str(cells.ref)
#
cells.query <- rownames(rna.combined@meta.data)[!rna.combined@meta.data$dataset %in% c('exp1', 'exp2')]
str(cells.query)
# estimate cell-cell distance
system.time(
cor.datasets <- cor(
      t(rna.combined@reductions$harmony@cell.embeddings[cells.query,]),
      t(rna.combined@reductions$harmony@cell.embeddings[cells.ref,]))
)

estimate_neighbors <- function(dst, k.nn = 100){
  # find neighbors
  cell.nn <- t(apply(dst, 1, function(x) rank(x) <= k.nn))
  y <- as.factor(old.scrna$clust.final[colnames(cell.nn)])
  #
  cell.clusters <- apply(cell.nn, 1, function(x) as.numeric(table(y[x])) )
  cell.clusters <- t(cell.clusters)
  return(cell.clusters)
}
cell.clusters <- estimate_neighbors(1 - cor.datasets, k.nn = 100)

projection.final <- apply(cell.clusters, 1, which.max)
# NA cells that tend to project to unnanotated regions
projection.final[apply(cell.clusters, 1, sum) < 50] <- NA
projection.final[rna.combined$harmony.clusters[names(projection.final)] %in%
                 harmony.clust.mask] <- NA
rna.combined$projected <- projection.final
rna.combined$projected[!is.na(rna.combined$clust.final)] <- rna.combined$clust.final[!is.na(rna.combined$clust.final)]

DimPlot(rna.combined, reduction = "umap.harmony", shuffle = T,
        #cells = c(rownames(cor.datasets), colnames(cor.datasets)),
        group.by = c("dataset", 'harmony.clusters', 'projected'),
        )


save(rna.combined, old.scrna, old.hemo, new.hemo,
        file = 'data/processed/rna_integration.RData')
saveRDS(list(harmony.clusters = rna.combined$harmony.clusters,
             projection.clusters = projection.final, cell.clusters = cell.clusters),
        file = "data/processed/clusters_annotation.rds")



############# run downstream analysis of RNA-based integrated datasets #################
#load('data/processed/rna_integration.RData') # upload saved object if available

### plot panels
pl1 <- DimPlot(rna.combined, reduction = "umap.harmony", shuffle = T,
               order = rev(c('exp1', 'exp2', 'C1', 'C2', 'D', 'T')),
               group.by = c("dataset"),
)
pl1 <- pl1 +  theme(legend.position = c(0.1,0.8), plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Samples")
pl1
#
pl2 <- DimPlot(rna.combined, reduction = "umap.harmony", shuffle = T,
               label.box = T, label = T, label.size = 3, repel = T,
               group.by = c("harmony.clusters")
)
pl2 <- pl2 + theme(legend.position = 'none', plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Harmony-based clusters")
pl2
#
pl3 <- DimPlot(rna.combined, reduction = "umap.harmony", shuffle = T,
               label.box = T, label = T, label.size = 3, repel = T,
               cols = old.scrna$clust.final.cols, order = names(rna.combined$projected),
               group.by = c("projected")
)
pl3
pl3 <- pl3 + theme(legend.position = 'none', plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Transferred clusters")
pl3
#
pl4 <- DimPlot(rna.combined, reduction = "umap.harmony", shuffle = T,
               #label.box = T, label = T, label.size = 3, repel = T,
               #cols = old.scrna$clust.final.cols,
               #order = names(rna.combined$projected),
               group.by = c("Phase")
)
pl4 <- pl4 + theme(legend.position = c(0.1,0.8), plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Cell cycle phase")
pl4
#

#
png('figures/rna_integration.png')#, width = 900, height = 300)
#pl
grid.arrange(pl1, pl2, pl3, pl4, ncol = 2)
dev.off()


### comparison of RNA patterns of scRNA-seq and multiome
rna.combined <- PrepSCTFindMarkers(rna.combined)
# estimate DE
Idents(rna.combined) <- "dataset"
n.cl <- length(unique(rna.combined$projected[!is.na(rna.combined$projected)]))
res <- mclapply(1:n.cl, function(i){
  print(i)
  x <- lapply(unique(rna.combined$dataset), function(datid){
    print(datid)
    markers <- FindMarkers(rna.combined,
                           ident.1 = i,
                           group.by = 'projected',
                           logfc.threshold = 0, min.pct = 0, min.diff.pct = 0,
                           subset.ident = datid)
    return(markers[, c('avg_log2FC', 'p_val_adj')])
  })
  y <- lapply(x, function(xx) rownames(xx))
  str(y)
  gg <-  Reduce(intersect, y)
  str(gg)
  logFC <- do.call(cbind, lapply(x, function(xx) xx[gg, 1]))
  padj <- do.call(cbind, lapply(x, function(xx) xx[gg, 2]))
  rownames(logFC) <- rownames(padj) <- gg
  colnames(logFC) <- colnames(padj) <- paste(unique(rna.combined$dataset), i, sep ='_')
  return(list(logFC = logFC, padj = padj))
}, mc.cores = n.cl)

#
gg <- Reduce(intersect, lapply(res, function(xx) rownames(xx[[1]])))
gg <- intersect(gg, intersect(rownames(new.hemo@assays$RNA), rownames(old.hemo@assays$RNA)))
str(gg)
logFC <- do.call(cbind, lapply(res, function(xx) xx[[1]][gg, ]))
rownames(logFC) <- gg
str(logFC)
padj <- do.call(cbind, lapply(res, function(xx) xx[[2]][gg, ]))
rownames(padj) <- gg
#
ind <- rna.combined$dataset %in% c('exp1', 'exp2')
rna.coverage.old <- rowSums(rna.combined@assays$RNA$counts[,ind] > 0)
rna.coverage.new <- rowSums(rna.combined@assays$RNA$counts[,!ind] > 0)
# plot RNA coverage
pl <- ggplot() + geom_point(aes(rna.coverage.old[gg] / sum(ind), rna.coverage.new[gg] / sum(!ind)), size = I(1), alpha = I(0.05)) +
  theme_bw() + scale_x_log10() + scale_y_log10() +
  labs(x = 'scRNA-seq', y = 'multiome', title = 'mRNA coverage') + theme(text = element_text(size = 20))
pl
pdf('figures/gene_coverage.pdf', width = 4, height = 4)
pl
dev.off()
#
ind <- rna.coverage.old[gg] > 1e+2 & rna.coverage.new[gg] > 1e+2 &
  apply(padj[gg,], 1, min) < 1e-20 & apply(abs(logFC[gg,]), 1, max) > 2
table(ind)
#
infcl <- do.call(rbind, lapply(colnames(logFC), function(x){
  x <- unlist(strsplit(x, '_'))
  if (x[1] %in% c('exp1', 'exp2')){
    c('RNA', x[2])
  }else{
    c('multiome', x[2])
  }
}))
infcl <- data.frame(infcl)
rownames(infcl) <- colnames(logFC)
colnames(infcl) <- c('type', 'cluster')
pheat <- pheatmap::pheatmap(cor(logFC[gg[ind],]), breaks = seq(-1, 1, length.out = 100),
                            annotation_col = infcl,
                            annotation_colors = list(cluster = old.scrna$clust.final.cols),
                            treeheight_row = 20, treeheight_col = 20,
                            show_rownames = F, show_colnames = F,
                            #annotation_row = infcl,
                            #col = old.scrna$clust.final.cols
)
pheat
#
pdf('figures/rna_integration_heatmap.pdf')#, width = 5, height = 4)
pheat
dev.off()



