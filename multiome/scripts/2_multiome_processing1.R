old.scrna <- new.env()
load(file = '../RNA_analysis/data/processed/scrna.RData', envir = old.scrna)
#
load(file = 'data/processed/multiome_input.RData')
sample.names <- names(inputdatal)
annotation.clusters <- readRDS(file = "data/processed/clusters_annotation.rds")
obj.combined <- readRDS(file = 'data/processed/multiome_seurat.rds')

### plot panels
pl1 <- DimPlot(obj.combined, reduction = "wnn.umap", shuffle = T,
               #order = rev(c('C1', 'C2', 'D', 'T')),
               group.by = c("dataset"),
)
pl1 <- pl1 +  theme(legend.position = c(0.1,0.8), plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Samples")
pl1
#
cll <- unique(obj.combined$projection[!is.na(obj.combined$projection)])
emb.clust <- do.call(rbind, lapply(cll, function(x){
  ind <- !is.na(obj.combined$projection) & obj.combined$projection == x
  c(colMeans(obj.combined@reductions$wnn.umap@cell.embeddings[ind,]))
}))
emb.clust <- data.frame(emb.clust);
emb.clust$clust <- cll
emb.clust$col <- old.scrna$clust.final.cols[cll]
colnames(emb.clust)[3:4] <- c('clust', 'col')
pl2 <- DimPlot(obj.combined, reduction = "wnn.umap", shuffle = T,
               cols = old.scrna$clust.final.cols,
               #cells = c(names(obj.combined$projection)[!is.na(obj.combined$projection)]),
               #order = names(obj.combined$projection)[order(-obj.combined$projection)],
               #order = names(obj.combined$projection),
               na.value = 'grey80',
               group.by = c("projection"),
               #label.box = T, label = T, label.size = 3,  repel = F,

)
pl2 <- pl2 + ggrepel::geom_label_repel(data = emb.clust, aes(x = wnnUMAP_1, y = wnnUMAP_2, label = clust), fill = emb.clust$col)
pl2 <- pl2 + theme(legend.position = 'none', plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Transferred clusters")
pl2

png('figures/multiome_integration.png', width = 600, height = 300)
pl2 + pl1
dev.off()

# plot suppl panels of RNA embedding
### plot panels
pl1 <- DimPlot(obj.combined, reduction = "umap.rna", shuffle = T,
               #order = rev(c('C1', 'C2', 'D', 'T')),
               group.by = c("dataset"),
)
pl1 <- pl1 +  theme(legend.position = c(0.1,0.8), plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Samples")
pl1
#
cll <- unique(obj.combined$projection[!is.na(obj.combined$projection)])
emb.clust <- do.call(rbind, lapply(cll, function(x){
  ind <- !is.na(obj.combined$projection) & obj.combined$projection == x
  c(colMeans(obj.combined@reductions$umap.rna@cell.embeddings[ind,]))
}))
emb.clust <- data.frame(emb.clust); emb.clust$clust <- cll
emb.clust$col <- old.scrna$clust.final.cols[cll]; colnames(emb.clust)[3:4] <- c('clust', 'col')
pl2 <- DimPlot(obj.combined, reduction = "umap.rna", shuffle = T,
               cols = old.scrna$clust.final.cols,
               #cells = c(names(obj.combined$projection)[!is.na(obj.combined$projection)]),
               #order = names(obj.combined$projection)[order(-obj.combined$projection)],
               #order = names(obj.combined$projection),
               na.value = 'grey80',
               group.by = c("projection"),
               #label.box = T, label = T, label.size = 3,  repel = F,

)
pl2 <- pl2 + ggrepel::geom_label_repel(data = emb.clust, aes(x = rnaUMAP_1, y = rnaUMAP_2, label = clust), fill = emb.clust$col)
pl2 <- pl2 + theme(legend.position = 'none', plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Transferred clusters")
pl2
#
png('figures/multiome_integration_rna.png', width = 600, height = 300)
pl2 + pl1
dev.off()

# plot suppl panels of ATAC embedding
### plot panels
pl1 <- DimPlot(obj.combined, reduction = "umap.atac", shuffle = T,
               #order = rev(c('C1', 'C2', 'D', 'T')),
               group.by = c("dataset"),
)
pl1 <- pl1 +  theme(legend.position = c(0.1,0.8), plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Samples")
pl1
#
cll <- unique(obj.combined$projection[!is.na(obj.combined$projection)])
emb.clust <- do.call(rbind, lapply(cll, function(x){
  ind <- !is.na(obj.combined$projection) & obj.combined$projection == x
  c(colMeans(obj.combined@reductions$umap.atac@cell.embeddings[ind,]))
}))
emb.clust <- data.frame(emb.clust); emb.clust$clust <- cll
emb.clust$col <- old.scrna$clust.final.cols[cll]; colnames(emb.clust)[3:4] <- c('clust', 'col')
pl2 <- DimPlot(obj.combined, reduction = "umap.atac", shuffle = T,
               cols = old.scrna$clust.final.cols,
               #cells = c(names(obj.combined$projection)[!is.na(obj.combined$projection)]),
               #order = names(obj.combined$projection)[order(-obj.combined$projection)],
               #order = names(obj.combined$projection),
               na.value = 'grey80',
               group.by = c("projection"),
               #label.box = T, label = T, label.size = 3,  repel = F,

)
pl2 <- pl2 + ggrepel::geom_label_repel(data = emb.clust, aes(x = atacUMAP_1, y = atacUMAP_2, label = clust), fill = emb.clust$col)
pl2 <- pl2 + theme(legend.position = 'none', plot.title = element_text(face = "plain", size = 16, family = "Arial")) +
  labs(x = '', y = '', title = "Transferred clusters")
pl2
#
png('figures/multiome_integration_atac.png', width = 600, height = 300)
pl2 + pl1
dev.off()


### assess correlation between RNA expression and promoter accesibility
clustname <- unique(obj.combined$seurat_clusters)
gene.clust <- do.call(cbind, lapply(clustname, function(cl){
  ind <- obj.combined$seurat_clusters == cl
  rowMeans(obj.combined@assays$SCT$scale.data[,ind], na.rm = T)
}))
colnames(gene.clust) <- clustname
prom.clust <- do.call(cbind, lapply(clustname, function(cl){
  ind <- obj.combined$seurat_clusters == cl
  rowMeans(obj.combined@assays$SCT_prom$scale.data[,ind], na.rm = T)
}))
colnames(prom.clust) <- clustname
#
gg <- intersect(rownames(gene.clust), rownames(prom.clust))
str(gg)
cr <- cor(t(gene.clust[gg,]), t(prom.clust[gg,]), method = 'spearman')
cr <- diag(cr)

#
gn1 <- 'MZB1'
gn2 <- 'SCN3A'
gn1 <- gn2
x <- rbind( data.frame(rna = gene.clust[gn1,], promoter = prom.clust[gn1,], gene = gn1),
            data.frame(rna = gene.clust[gn2,], promoter = prom.clust[gn2,], gene = gn2))
p1 <- ggplot(x) + geom_point(aes(rna, promoter), size = I(3)) + theme_bw() +
  labs(x='mRNA expression', y='Promoter\naccessibility', font = 2) + theme(text = element_text(size = 20), plot.title = element_text(face = "italic")) +
  facet_wrap(facets = ~ gene, scales = 'free') +
  geom_text(data = data.frame(gene = unique(x$gene),
                              #label = round(cr[unique(x$gene)], 2),
                              #label = bquote(paste(rho, ' = ' ) ),
                              #label = 'fd',
                              x = max(x$rna) - 0.0, y = min(x$promoter) + 0.0),
            aes(x = x, y = y),
            label = paste(round(cr[c(gn1, gn2)], 2)),#bquote(paste(italic(rho), ' = ')),#bquote(paste(rho, ' = ' ) ),
            size = 6, vjust = -1, hjust = 1
  )
p1
#
pdf('figures/rna_promoter_scatterplot.pdf', width = 6, height = 3)
p1
dev.off()

p2 <- qplot(cr, bins = 50, xlab = 'RNA vs Promoter acc., Spearman') +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() + theme(text = element_text(size = 20))
p2
#
pdf('figures/rna_promoter_hist.pdf', width = 5, height = 3)
p2
dev.off()


###### analysis of motif activities across cells
motif.object <- obj.combined@assays$ATAC@motifs
chvar <- obj.combined@assays$chromvar@data
chvar[is.na(chvar)] <- 0
rownames(chvar) <- as.character(unlist(motif.object@motif.names)[rownames(chvar)])
#
gg <- intersect(rownames(chvar), rownames(obj.combined@assays$SCT$scale.data))
str(gg)
head(sort(apply(chvar[gg,], 1, sd), decreasing = T), 20)
gn <- 'GATA2' #'GATA2' #'CEBPD'
#gn <- names(sort(diag(xx), decreasing = T))[17]
gn <- names(sort(apply(chvar[gg,], 1, sd), decreasing = T))[8]
pdf('figures/motif_expression_patterns_new.pdf', width = 8, height = 4)
for (gn in names(sort(apply(chvar[gg,], 1, sd), decreasing = T))){
  print(gn)
  #gn <- 'CEBPD'#'MYC'
  motif.name <- ConvertMotifID(obj.combined, name = gn)
  gene_plot <- FeaturePlot(obj.combined, features = paste("sct", gn, sep = '_'), order = T,
                           #raster = TRUE, raster.dpi = c(1000, 1000), pt.size = 6,
                           reduction = 'wnn.umap')#umap.naive.rna 'umap.harmony')
  motif_plot <- FeaturePlot(obj.combined, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), order = T,
                            #raster = TRUE, raster.dpi = c(1000, 1000), pt.size = 6,
                            reduction = 'wnn.umap')
  print(gene_plot | motif_plot)
}
dev.off()

# plot motif activity of cell type-specific TFs:
gn <- c('CEBPD', 'GATA2', 'EBF1', 'ZEB1')
motif.name <- ConvertMotifID(obj.combined, name = gn)
motif_plot <- FeaturePlot(obj.combined, features = motif.name, min.cutoff = 0,
                          cols = c('grey90', "darkred"),
                          reduction = 'wnn.umap', ncol = 2, pt.size = I(0.5),
                          combine = FALSE)
motif_plot
for (i in 1:length(motif_plot)){
  motif_plot[[i]] <- motif_plot[[i]] + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                                             legend.position = 'none', title = element_blank(),
                                             #plot.title = element_text(face = "plain", size = 16, family = "Arial"),
  ) + geom_text(x = -Inf, y = Inf, label = gn[i], #font = 1,
                hjust = -0, vjust = 3, size = 10) +
    labs(x = '', y = '', title = NULL)
}
#grid.arrange(motif_plot, ncol = 2)
#pdf('figures/motif_activity_examples.pdf')#, width = 400, height = 400)
#png('figures/motif_activity_examples.png')#, width = 400, height = 400)
print(do.call(grid.arrange, c(motif_plot, ncol = 2)))
dev.off()



###### cluster-specific motif enrichment analysis
DefaultAssay(obj.combined) <- 'ATAC'
Idents(obj.combined) <- "projection"
meta.feature <- GetAssayData(obj.combined, assay = "ATAC", slot = "meta.features")

enriched.motifs.list <- list()
pdf('figures/motif_cluster_specific_new.pdf', width = 6, height = 8)
for (proj.clust in 1:10){
  print(proj.clust)
  #proj.clust <- '8'
  da_peaks <- FindMarkers(
    object = obj.combined,
    ident.1 = proj.clust,
    test.use = 'wilcox', #'DESeq2', #'LR',
    latent.vars = 'nCount_ATAC'
  )
  open_clust <- rownames(da_peaks[da_peaks$p_val_adj < 1e-2 &
                                    da_peaks$avg_log2FC > 2, ])
  str(open_clust)
  if (length(open_clust) > 0){
    open.peaks <- AccessiblePeaks(obj.combined, idents = c(proj.clust))
    peaks.matched <- MatchRegionStats(
      meta.feature = meta.feature[open.peaks, ],
      query.feature = meta.feature[open_clust, ],
      features.match = c("GC.percent"),
      n = 50000
    )
    enriched.motifs <- FindMotifs(
      object = obj.combined,
      features = open_clust,
      background = peaks.matched,
    )
    enriched.motifs.list[[proj.clust]] <- enriched.motifs
    pl1 <- MotifPlot(object = obj.combined,
                     motifs = head(rownames(enriched.motifs), 10)
    )
    enriched.motifs$motif.name.show <- enriched.motifs$motif.name
    enriched.motifs$motif.name.show[enriched.motifs$p.adjust > 1e-2] <- NA
    pl2 <- ggplot(enriched.motifs, aes(fold.enrichment, pvalue, label = motif.name.show)) + geom_point(col = old.scrna$clust.final.cols[proj.clust]) + theme_bw() +
      scale_x_log10() + scale_y_log10() + labs(x = 'motif fold enrichment', y = 'p-value', title = paste('cluster', proj.clust)) + theme(text = element_text(size = 20))
    pl2 <- pl2 + ggrepel::geom_text_repel()
    print(grid.arrange(pl2, pl1))
  }
}
dev.off()

### plot heatmap of top enriched motifs
enriched.motifs <- enriched.motifs.list[[1]]
fc <- do.call(cbind, lapply(enriched.motifs.list, function(x){
  if (!is.null(x)) {
    x[rownames(enriched.motifs),]$fold.enrichment
  }else{
    rep(1, nrow(enriched.motifs))
  }
}))
rownames(fc) <- enriched.motifs$motif.name
colnames(fc) <- 1:10
padj <- do.call(cbind, lapply(enriched.motifs.list, function(x) if (!is.null(x)) {x[rownames(enriched.motifs),]$p.adj}else{rep(1, nrow(enriched.motifs))} ))
#ind <- apply(padj, 1, min) < 1e-5 & apply(fc, 1, max) > 2.5
ind <- unique(unlist(lapply(1:10, function(i){
  ind1 <- padj[,i] < 1e-3 & fc[,i] > 1.5
  if (sum(ind1) > 0){
    names(sort(fc[ind1,i], decreasing = T))[1:min(sum(ind1), 15)]
  }
})))
length(ind)
pheat <- pheatmap::pheatmap(fc[ind,], main = 'cluster-specific motif enrichment', annotation_col = data.frame(cluster = 1:10),
                            annotation_colors = list(cluster = old.scrna$clust.final.cols), angle_col = 0, annotation_legend = FALSE,
)#, breaks = 2^seq(-3, 3, length.out = 100) )
pheat

pdf('figures/motif_cluster_heatmap.pdf', width = 5, height = 8)#, width = 6, height = 8)
pheat
dev.off()


