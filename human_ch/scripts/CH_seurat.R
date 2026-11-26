### seurat processing of CH data
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


### process CH xenograft data
old.hemo.all <- CreateSeuratObject(counts = old.scrna$cd[,cells.old],
                                   meta.data = meta.old[cells.old,])
old.hemo.all@meta.data <- cbind(old.hemo.all@meta.data,
                                meta.old[cells.old, c('exp', 'genotype', 'sampleid', 'clust.final')])
old.hemo.all@meta.data$dataset <- old.hemo.all@meta.data$exp
old.hemo.all[["RNA"]] <- split(old.hemo.all[["RNA"]], f = meta.old[cells.old,]$exp)
old.hemo.all <- NormalizeData(old.hemo.all)
old.hemo.all <- FindVariableFeatures(old.hemo.all)
old.hemo.all <- ScaleData(old.hemo.all)
old.hemo.all <- RunPCA(old.hemo.all)
#
old.hemo.all <- FindNeighbors(old.hemo.all, reduction = "pca", dims = 1:30)
old.hemo.all <- FindClusters(old.hemo.all, resolution = 2, cluster.name = "pca_clusters")
old.hemo.all <- RunUMAP(old.hemo.all, reduction = "pca", dims = 1:30, reduction.name = "umap.pca")
old.hemo.all[["RNA"]] <- JoinLayers(old.hemo.all[["RNA"]])
#
emb_original <- CreateDimReducObject(embeddings = old.scrna$emb, key = "umap.original", assay = DefaultAssay(old.hemo.all))
old.hemo.all[["umap.original"]] <- emb_original

# add cell cycle score
old.hemo.all <- CellCycleScoring(old.hemo.all, s.features = genes.s, g2m.features = genes.g2m, set.ident = TRUE)
gns.show <- c("PCNA", "TOP2A", "MCM6", "MKI67")
par(mfrow=c(2,2), mar = rep(2,4))
for(gn in gns.show) boxplot(old.hemo.all@assays$RNA$scale.data[gn,] ~ old.hemo.all$Phase, main = gn)

# add Neu and Mono scores
old.hemo.all <- AddModuleScore(object = old.hemo.all,
  features = list(rownames(gs.atlas)[gs.atlas[,'preNeu'] > 0.2],
                  rownames(gs.atlas)[gs.atlas[,'Mono-1'] > 0.2]),
  name = c("preNeu", "Mono"), nbin = 24)

### save seurat object
saveRDS(old.hemo.all, 'data/xeno_ch_seurat.rds')

#
pl_bias <- FeaturePlot(old.hemo.all, reduction = "umap.original",
            features = c('preNeu1', 'Mono2'),
            cols = c("grey95", "red"), pt.size = 0.5, order = TRUE)
pl_bias

# estimate density of Mono-biased and Neu-biase clusters
DimPlot(old.hemo.all, reduction = "umap.original", cells = rownames(old.hemo.all@reductions$umap.original),
        group.by = c('pca_clusters', 'pca_clusters_resolved', 'clust.final'), ncol = 2)

#  reanalyse Myeloid cluster
cells.my = names(which(!is.na(old.hemo.all$clust.final) & old.hemo.all$clust.final == 8))
# manually choose Mono-biased and Neu-biased cells
pl_select <- qplot(subset_hem@reductions$umap.original@cell.embeddings[,1],
      subset_hem@reductions$umap.original@cell.embeddings[,2], size = I(0.2), alpha = I(0.5),
      col = old.hemo.all$Mono2[cells.my] > 0.35) + scale_color_discrete(name = "Mono_biased") + theme_bw()

mono_bias <- data.frame(sid = old.hemo.all$sampleid[cells.my], gid = old.hemo.all$genotype[cells.my],
           mono_bias = as.numeric(old.hemo.all$Mono2[cells.my] > 0.35))
mono_bias_stat <- mono_bias %>% group_by(sid, gid, mono_bias)  %>%
   summarise(count = n()) %>% ungroup()
mono_bias_stat <- data.frame(mono_bias_stat)
mono <- mono_bias_stat[mono_bias_stat$mono_bias == '1',]
neu <- mono_bias_stat[mono_bias_stat$mono_bias == '0',]
rownames(mono) <- mono$sid; rownames(neu) <- neu$sid;
# hsc counts
cells.hsc <- names(old.scrna$clust.final)[which(old.scrna$clust.final %in% '1')]
hsc_stat <- table(old.scrna$sampleid[cells.hsc])

pl_bias_est <- qplot(mono[names(hsc_stat),]$count / hsc_stat, neu[names(hsc_stat),]$count / hsc_stat,
      col = mono[names(hsc_stat),]$gid, shape = hash.match[names(hsc_stat),]$exp,
      size = I(5), alpha = I(1)) +
  scale_color_discrete(name = "Genotype") + scale_shape_discrete(name = "Exp ID") +
  xlab('Mono-biased to HSC') + ylab('Neu-biased to HSC') + theme_bw()

require(gridExtra)
pdf(paste('figures/mono_neu_bias.pdf', sep = ''), width = 8, height = 8)
grid.arrange(pl_bias[[1]], pl_bias[[2]], pl_select, pl_bias_est, ncol = 2)
dev.off()

### analysis of cycling cells
old.hemo.all <- RunPCA(old.hemo.all, features = c(genes.s, genes.g2m),
                       reduction.name = 'pca_cc')
DimPlot(old.hemo.all, reduction = "pca_cc", group.by = c('Phase'),
        size = I(0.02), alpha = I(0.2))
cycle_bias <- data.frame(sid = old.hemo.all$sampleid,
                        gid = old.hemo.all$genotype,
                        clustid = old.hemo.all$clust.final,
                        phase = old.hemo.all$Phase)
#
cycling_bias_stat <- cycle_bias %>% group_by(sid, gid, clustid, phase)  %>%
  summarise(count = n()) %>% ungroup()
cycling_bias_stat <- as.data.frame(cycling_bias_stat)
#
cycling.s <- cycling_bias_stat[cycling_bias_stat$phase %in% c('S'), -4]
rownames(cycling.s) <- paste(cycling.s$sid, cycling.s$clustid, sep = '_')
cycling.g2m <- cycling_bias_stat[cycling_bias_stat$phase %in% c('G2M'), -4]
rownames(cycling.g2m) <- paste(cycling.g2m$sid, cycling.g2m$clustid, sep = '_')
cycling <- cycling.s; cycling$count <-  cycling$count + cycling.g2m[rownames(cycling.s),]$count
noncycling <- cycling_bias_stat[cycling_bias_stat$phase %in% c('G1'), -4]
rownames(noncycling) <- paste(cycling$sid, cycling$clustid, sep = '_')

gn <- 'TET2'
ind <- hash.match[cycling$sid,]$exp == '123'
if (gn == 'TET2') ind <-  hash.match[cycling$sid,]$exp == '127'
pl <- ggplot() + geom_boxplot(aes(as.factor(cycling$clustid[ind]), cycling$count[ind] / (cycling$count[ind]+noncycling[rownames(cycling),]$count[ind]),
                            col = cycling$gid[ind])) + geom_jitter(aes(x = as.factor(cycling$clustid[ind]),
                                                                       y = cycling$count[ind] / (cycling$count[ind] + noncycling[rownames(cycling),]$count[ind]),
                                                                       col = cycling$gid[ind]),
                                                                   width = 0.1,   # Adjusts the amount of horizontal jitter
                                                                   size = 3,      # Adjusts the size of the points
                                                                   alpha = 0.5) + xlab('cluster') + ylab('cycling cells, proportion')
pdf(paste('figures/cycling_prop_',gn,'.pdf', sep = ''))
pl
dev.off()
