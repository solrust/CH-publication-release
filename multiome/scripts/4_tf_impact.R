DefaultAssay(obj.combined) <- "ATAC"

# first compute the GC content for each peak
obj.combined <- RegionStats(obj.combined, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
timing_result <- system.time({
obj.combined <- LinkPeaks(
  object = obj.combined,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  distance = 1e+05,
  score_cutoff = 0.05,
  genes.use = rownames(obj.combined@assays$SCT$counts)#[1:1000]#c("LYZ", "MS4A1")
)
})
timing_result
g2p.seurat <- Links(obj.combined)
g2p.seurat1 <- do.call(rbind, lapply(unique(g2p.seurat$gene), function(gn){
  table(as.factor(g2p.seurat$peak)[which(g2p.seurat$gene == gn)])
}))
rownames(g2p.seurat1) <- unique(g2p.seurat$gene)


### find TSS and gene-peaks
gs <- unique(annotation$gene_name)
tss.info <- do.call(rbind, mclapply(gs, function(gn){
  print(gn)
  #annotation[annotation$gene_name == gn,]
  chr <- unique(as.character(annotation[annotation$gene_name == gn,]@seqnames))[1]
  strand <- unique(as.character(annotation[annotation$gene_name == gn,]@strand))[1]
  if (strand == '+'){
    tss <- min(annotation[annotation$gene_name == gn,]@ranges@start)
  }else{
    tss <- max(annotation[annotation$gene_name == gn,]@ranges@start +
                 annotation[annotation$gene_name == gn,]@ranges@width)
  }
  return((cbind(chr, tss, strand)))
}, mc.cores = 20))
rownames(tss.info) <- gs
#
window <- 1e+5
chrs <- as.character(obj.combined@assays$ATAC@ranges@seqnames)
pos <- obj.combined@assays$ATAC@ranges@ranges@start + obj.combined@assays$ATAC@ranges@ranges@width/2
g2p <- do.call(rbind, lapply(1:nrow(tss.info), function(ii){
  #if (ii %% 100 == 0) print(ii)
  print(ii)
  gn <- rownames(tss.info)[ii]
  #gn <- rownames(tss.info)[2]
  coord <- as.numeric(tss.info[gn, 2])
  indx <- chrs == tss.info[gn,1] & pos > coord - window &
    pos < coord + window
  indx + 0
}))
g2p <- Matrix::Matrix(g2p, sparse = T)
rownames(g2p) <- rownames(tss.info)
colnames(g2p) <- rownames(obj.combined@assays$ATAC@meta.features)
dim(g2p)
qplot(rowSums(g2p))
qplot(colSums(g2p))

### assess correlation between RNA expression and peak accesibility
obj.combined <- FindClusters(obj.combined, graph.name = "wsnn", cluster.name = 'metacells',
                             algorithm = 3, resolution = 3, verbose = FALSE)
sort(table(obj.combined$metacells))
length(table(obj.combined$metacells))
clustname <- unique(obj.combined$seurat_clusters)
gene.clust <- do.call(cbind, lapply(clustname, function(cl){
  ind <- obj.combined$seurat_clusters == cl
  #rowMeans(obj.combined@assays$SCT$scale.data[,ind], na.rm = T)
  rowMeans(obj.combined@assays$RNA$counts[,ind], na.rm = T)
}))
gene.clust <- t(t(gene.clust) / colSums(gene.clust))
colnames(gene.clust) <- clustname
#
atac.clust <- do.call(cbind, lapply(clustname, function(cl){
  ind <- obj.combined$seurat_clusters == cl
  rowSums(obj.combined@assays$ATAC$counts[,ind], na.rm = T)
}))
atac.clust <- t(t(atac.clust) / colSums(atac.clust))
colnames(atac.clust) <- clustname


### estimate gene-peak correlation:
peaks.choose <- rownames(atac.clust)
str(peaks.choose)
gns <- intersect(rownames(g2p), rownames(gene.clust))
str(gns)
g2p.cor <- cor(t(gene.clust[gns,]), t(atac.clust[peaks.choose,]), method = 'spearman')
g2p.cor <- g2p[gns, peaks.choose] * g2p.cor
#
nperm <- 20
g2p.cor.control.list <- mclapply(1:nperm, function(i){
  print(i)
  cr <- cor(t(gene.clust[gns,]), t(atac.clust[peaks.choose, sample(1:ncol(atac.clust))]), method = 'spearman')
  g2p[gns, peaks.choose] * cr
}, mc.cores = 10)
length(g2p.cor.control.list)
g2p.cor.control.mean <- Reduce('+', g2p.cor.control.list)/length(g2p.cor.control.list)
g2p.cor.control.sd <- sqrt(Reduce('+', lapply(g2p.cor.control.list, function(x) (x - g2p.cor.control.mean)^2 ))/length(g2p.cor.control.list))
z.cor <- (g2p.cor[gns, peaks.choose] - g2p.cor.control.mean) / g2p.cor.control.sd
z.cor[is.na(z.cor) | g2p.cor.control.sd == 0] <- 0
znn <- rowSums(abs(z.cor) > 4)
zmax <- apply(z.cor, 1, max)
qplot(g2p[names(sort(znn, decreasing = T))[3], ],
      z.cor[names(sort(znn, decreasing = T))[3],])

y1 <- rowSums(g2p); names(y1) <- rownames(tss.info)
y <- as.numeric(table(x$gene)); names(y) <- names(table(x$gene))
qplot(y, y1[names(y)])
qplot(znn[names(y)], y1[names(y)])
qplot(znn[names(y)], y[names(y)])
cor(znn[names(y)], y[names(y)], use = 'na.or.complete')

#head(sort(table(x$gene), decreasing = T), 30)
Idents(obj.combined) <- 'projection'
idents.plot <- names(table(obj.combined$projection))
#gn <- names(sort(table(x$gene), decreasing = T))[5]
ginter <- rownames(tss.info)[which(rownames(tss.info) %in% rownames(obj.combined@assays$SCT$scale.data))]
#table(unique(annotation$gene_name) %in% rownames(obj.combined@assays$SCT$scale.data))
gn <- ginter[4]
gn <- names(sort(zmax, decreasing = T))[2]
#gn <- names(sort(y, decreasing = T))[2]
p1 <- CoveragePlot(
  object = obj.combined,
  region = gn,
  features = gn,
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 20000,
  extend.downstream = 20000
)
p1
gn
tss.info[gn,]

devtools::install_github("buenrostrolab/FigR")
BiocManager::install("SummarizedExperiment")

