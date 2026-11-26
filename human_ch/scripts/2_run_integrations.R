library(ggplot2)
library(ggrepel)
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

#
cells <- cells.old
cells1 <- sample(cells, 20000)
str(cells1)
old.hemo <- CreateSeuratObject(counts = old.scrna$cd[, cells1],
                               meta.data = meta.old[cells1,])
old.hemo@meta.data <- cbind(old.hemo@meta.data,
                            meta.old[cells1, c('exp', 'genotype', 'sampleid', 'clust.final')])
old.hemo@meta.data$dataset <- old.hemo@meta.data$exp
old.hemo[["RNA"]] <- split(old.hemo[["RNA"]], f = meta.old[cells1,]$exp)
old.hemo <- NormalizeData(old.hemo)
old.hemo <- FindVariableFeatures(old.hemo)
old.hemo <- ScaleData(old.hemo)
old.hemo <- RunPCA(old.hemo)


#  joint integration of CH and Vyas
run.combined <- FALSE
if (run.combined){
rna.combined <- merge(x = old.hemo, y = obj.human)
rna.combined <- NormalizeData(rna.combined)
rna.combined <- FindVariableFeatures(rna.combined)
rna.combined <- ScaleData(rna.combined)
rna.combined <- RunPCA(rna.combined)
rna.combined$exp[is.na(rna.combined$exp)] <- 'human'
#
options(future.globals.maxSize = 1 * 2000 * 1024^2)
rna.combined <- IntegrateLayers(object = rna.combined, method = RPCAIntegration,
                                orig.reduction = "pca",
                                new.reduction = "rpca", verbose = FALSE)
rna.combined <- FindNeighbors(rna.combined, reduction = "rpca", dims = 1:30)
rna.combined <- FindClusters(rna.combined, resolution = 2, cluster.name = "rpca_clusters")
rna.combined <- RunUMAP(rna.combined, reduction = "rpca", dims = 1:30, reduction.name = "umap.rpca")
DimPlot(rna.combined, reduction = "umap.rpca", shuffle = T,
        group.by = c('exp', 'clust.final', 'Cluster', 'Sample'), , na.value = "black")
}else{
  rna.combined <- readRDS(file = 'data/integrated_seurat.rds')
}

### annotate clusters
clust_labels <- c('HSC', 'LMPP-1', 'MLP', 'Pro-B', 'pDC', 'LMPP-2', 'early GMP', 'late GMP', 'EMPP/MkP', 'EryP')
names(clust_labels) <- 1:10
#write.csv(data.frame(clust_labels), file = 'data/annotation_clusters.csv')

# match clusters
datx <- as.data.frame(cbind(rna.combined@reductions$umap.rpca@cell.embeddings))
colnames(datx) <- c('emb1', 'emb2')
datx$Annotation <- rna.combined@meta.data[rownames(datx),]$Cluster
datx$cluster <- as.character(rna.combined@meta.data[rownames(datx),]$clust.final)
#
ind <- !is.na(datx$Annotation)
datx_mean <- aggregate(datx[, 1:2], by = list(Annotation = as.factor(datx$Annotation)), FUN = mean, na.rm = TRUE)
pl1 <- ggplot() + geom_point(data = datx[!ind, ], aes(emb1, emb2), size = 0.5, col = 'grey80') +
           geom_point(data = datx[ind,], aes(emb1, emb2, col = Annotation), size = 0.5) + theme_bw() +
           geom_text_repel(data = datx_mean, aes(x = emb1, y = emb2, label = Annotation),
                           size = 4) +
           theme(legend.text = element_text(size = 8), legend.title = element_blank(),
                     legend.key.height = unit(0.5, "cm")) +
          guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) +
          labs(x = '', y = '')
pl1
#
ind <- !is.na(datx$cluster)
table(ind)
datx_mean <- aggregate(datx[, 1:2], by = list(cluster = as.factor(datx$cluster)), FUN = median, na.rm = TRUE)
pl2 <- ggplot() + geom_point(data = datx[!ind, ], aes(emb1, emb2), size = 0.5, col = 'grey80') +
  geom_point(data = datx[ind,], aes(emb1, emb2, col = cluster), size = 0.5) + theme_bw() +
  geom_text_repel(data = datx_mean, aes(x = emb1, y = emb2,
                                        #label = cluster
                                        label = paste(cluster, ": ",clust_labels[as.character(cluster)], sep = '')
                                        ),
                  col = 'black', size = 6) +
  theme(legend.text = element_text(size = 16), legend.title = element_blank(),
        legend.key.height = unit(0.5, "cm")) +
  scale_color_manual(values = old.scrna$clust.final.cols,
                     labels = clust_labels) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) +
  labs(x = '', y = '')
pl2

png("figures/annotation.png", width = 1200, height = 500)
pl1 + pl2
dev.off()

# visualize probabilistic matching of cell types across datasets
DimPlot(rna.combined, reduction = "umap.rpca", shuffle = T,
        group.by = c('exp', 'clust.final', 'rpca_clusters'), ncol = 2)
clust.overlap_xeno <- table(rna.combined$rpca_clusters,
                       clust_labels[rna.combined$clust.final])
clust.overlap_human <- table(rna.combined$rpca_clusters,
                        rna.combined$Cluster)
xx <- cor(clust.overlap_xeno, clust.overlap_human)
pdf('figures/cluster_matching.pdf')
pheatmap::pheatmap(t(xx[order(apply(xx, 1, which.max)),]),
                   cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

# save seurat object
saveRDS(rna.combined, file = 'data/integrated_seurat.rds')

umap_coords <- Embeddings(rna.combined, reduction = "umap.rpca")
metadata_plus_umap <- cbind(rna.combined@meta.data, umap_coords)
saveRDS(metadata_plus_umap, file = 'data/integrated_seurat_metadata.rds')




