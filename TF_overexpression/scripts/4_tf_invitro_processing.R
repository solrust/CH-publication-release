library(Matrix)
library(pagoda2)
library(ggplot2)
library(parallel)
library(scde)
library(RColorBrewer)


### read out scRNA-seq
cd.tfvitro <- read.10x.matrices(paste('data/raw/rna_filtered_feature_bc_matrix',sep=''))
rownames(cd.tfvitro) <- make.unique(rownames(cd.tfvitro))
#
p2.tfinvitro <- Pagoda2$new(cd.tfvitro, log.scale=TRUE, n.cores=20)
p2.tfinvitro$adjustVariance(plot=T, gam.k=10)
p2.tfinvitro$calculatePcaReduction(nPcs=100, n.odgenes=5e3, maxit = 300)
### make PCA-based analysis
p2.tfinvitro$makeKnnGraph(k=30, type='PCA', center=T, distance='cosine');
p2.tfinvitro$getEmbedding(type='PCA', embeddingType = 'UMAP', name = 'UMAP', distance = 'cosine', #distance = 'euclidean', #spread = 1,
                          n_neighbors = 30)
#p2.tfinvitro$makeGeneKnnGraph()
p2.tfinvitro$getKnnClusters(method = igraph::cluster_leiden, type = "PCA", resolution_parameter = 0.5, objective_function = 'modularity', name = "leiden")
length(unique(p2.tfinvitro$clusters$PCA$leiden))
p2.tfinvitro$plotEmbedding(type = "PCA", embeddingType = 'UMAP', clusterType = 'leiden')
p2.tfinvitro$getDifferentialGenes(type='PCA',verbose=T,clusterType='leiden')
#saveRDS(p2.tfinvitro, file = 'data/processed/p2_tfinvitro.rds')

### upload p2 object
p2.tfinvitro <- readRDS(file = 'data/processed/p2_tfinvitro.rds')

# get DE genes
# build heatmap of DE genes
gns <- vector()
for (i in 1:length(p2.tfinvitro$diffgenes$PCA[[1]])){
  de <- p2.tfinvitro$diffgenes$PCA[[1]][[i]]
  #de <- de[abs(de$Z) > 8, ]
  de <- de[abs(de$Z) > 8 & abs(de$M > 1), ]
  #gns <- c(gns, unique(c(rownames(de)[order(-de$Z)][1:min(10, nrow(de))], rownames(de)[order(de$Z)][1:min(10, nrow(de))])) )
  gns <- c(gns, rownames(de)[order(-de$Z)][1:min(8, nrow(de))])
}
#
gns.add.hsc <- c('HLF', 'PROCR', 'HOPX', 'AVP', 'MEIS1') # Clusters 1 and 6 (HSC)
gns.add.b <-  c('CD79A', 'MS4A1', 'VPREB1/3', 'DNTT', 'RAG1/2') # Cluster 5 (B progenitors)
gns.add.dc <- c('HLA-DRA', 'CD74', 'LST1', 'FCER1A', 'CLEC10A', 'IRF7', 'GZMB') # Cluster 7 9APC/DC)
gns.add <- c(gns.add.hsc, gns.add.b, gns.add.dc)
gns <- c(gns, gns.add)
gns <- unique(gns)
gns <- gns[!is.na(gns) & gns %in% colnames(p2.tfinvitro$counts)]
str(gns)

# order by max expression
y <- do.call(cbind, lapply(p2.tfinvitro$diffgenes$PCA[[1]], function(x) x[gns, 'M']))
rownames(y) <- gns
yy <- apply(y, 1, function(x) which.max(x)[1])
# remove non-informative add-on genes
yy <- yy[!names(yy) %in%
           c(gns.add.hsc[!yy[gns.add.hsc] %in% c(1,6)],
             gns.add.b[!yy[gns.add.b] %in% c(5)],
             gns.add.dc[!yy[gns.add.dc] %in% c(7:9)])]

#
#pdf('figures/invitro/heatmap_de.pdf', height = 8, width = 8)
#p2.tfinvitro$plotGeneHeatmap(genes = gns, groups = p2.tfinvitro$clusters$PCA[[1]], cluster.genes = T, inner.clustering = T)
#dev.off()
grp <- p2.tfinvitro$clusters$PCA[[1]]
min.size <- min(table(grp))
grp.sample <- unlist(lapply(unique(grp), function(x){
  grp1 <- grp[grp %in% x]
  sample(grp1, min.size + (length(grp1)-min.size)^(0.6))
}))
pdf('figures/revision/heatmap_de_1.pdf', height = 6*1, width = 3*1)
#p2.tfinvitro$plotGeneHeatmap(genes = gns, groups = grp.sample, useRaster = F, box = T)
p2.tfinvitro$plotGeneHeatmap(genes = names(sort(yy[!is.na(yy)])), groups = grp.sample, useRaster = F, box = T)

dev.off()


pdf('figures/invitro/emb.pdf', width = 5, height = 5)
p2.tfinvitro$plotEmbedding(type = "PCA", embeddingType = 'UMAP', clusterType = 'leiden')
dev.off()

### show activity of annotate programs on top of embedding
progs_curate <- readRDS("../RNA_analysis/data/genesets/processed/fate_markers.rds")
show_gene <- F
#pdf("/home/solrust/CH/draft/figures/TF_overexpression/invitro/progs_activity.pdf", width = 4 * 4, height = 4 * 4)
png('figures/invitro/map_programs.png', width = 200 * 4, height = 200 * 4)
par(mfrow=c(4,4), mar = rep(0.5, 4))
pp <- 0.3
gnout <- vector()
for (k in setdiff(1:length(progs_curate), 13)){
  print(k)
  inf1 <- project_gset(gset = progs_curate[[k]], mat = p2.tfinvitro$counts, pr = NULL)
  if (!is.null(inf1)){
    if (!show_gene){
      x <- inf1$av.fit
      #x <- sign(x) * log(abs(x) * 1e+8 + 1)^1
      nm <- names(progs_curate)[k]
    }else{
      #nm <- gn <- names(which.max(abs(inf1$gscore)))
      #nm <- gn <- names(which.max(abs(inf1$gscore[setdiff(names(inf1$gscore), gnout)] )))
      nm <- gn <- names(sort(abs(inf1$gscore[setdiff(names(inf1$gscore), gnout)] ), decreasing = T))[1]
      gnout <- c(gnout, gn)
      #nm <- gn <- inf1$gscore[intersect(names(inf1$gscore), gns)]
      x <- p2.tfinvitro$counts[, gn]
    }
    if (!is.null(nrow(x))){
      x <- x[intersect(rownames(x), rownames(p2.tfinvitro$embeddings$PCA$UMAP)),]
    }else{
      x <- x[intersect(names(x), rownames(p2.tfinvitro$embeddings$PCA$UMAP))]
    }
    if (!show_gene){
      #x <- sign(x) * (abs(x) * 1e+8 + 1)^pp
      #x <- sign(x) * log(abs(x) * 1e+8 + 1)^1
    }

    #for (i in 1:ncol(x)){
    cols <- colorRampPalette(c("grey90","grey80", "red" ))(n = 20)[as.numeric(cut(x,20))]
    #cols <- colorRampPalette(c("grey90", "red" ))(n = 20)[as.numeric(cut(x,20))]
    plot(p2.tfinvitro$embeddings$PCA$UMAP[names(x[]),], #main = names(progs_curate)[k],
         #bg = cols, col = 'grey30',
         col = cols,
         #lwd = 0.05, #p2$clusters$cleaned$leiden1[cls],
         xlab = '', ylab = '', pch = 21, cex = 0.2, xaxt = 'n', yaxt = 'n')#, bty = 'n')
    legend('topright', nm[], bty = 'n', cex = 2.5)
    gns.show <- names(sort(abs(inf1$gscore), decreasing = T))[1:9]

    gns.show <- unlist(lapply(gns.show, function(gn){
      paste(substr(gn, 1, 1), tolower(substr(gn, 2, nchar(gn))), sep = '')
    }))

    #unlist(lapply(gns.show, function(gn) substitute(paste(italic(gn)))))
    #legend('bottomright', legend = gns.show, bty = 'n', cex = 1.5)

    #}
  }
}
dev.off()



############### read out scRNA-seq of TF amplicones
cd.tf.ampl <- read.10x.matrices('data/raw/amplicon_raw_feature_bc_matrix/')
cd.ampl <- cd.tf.ampl[(nrow(cd.tfvitro)+1):nrow(cd.tf.ampl), ]
str(cd.ampl)

qplot(colSums(cd.tf.ampl), colSums(cd.ampl), xlab = 'UMIs per cell', ylab = 'barcode UMIs per cells', size = I(1)) + geom_abline(col='red')
qplot(colSums(cd.ampl), apply(cd.ampl, 2, max), xlab = 'TF barcode UMIs, total', ylab = 'TF barcode UMIs, dominant barcode', size = I(1)) +
  theme_bw() + geom_abline(slope=1, linetype = 2, col = 'red')
bid <- rownames(cd.ampl)[apply(cd.ampl, 2, which.max)]
cell.map.tf <- data.frame(cell = colnames(cd.ampl), TFbarcode = bid, TF = barcode2tf$Gene[match(bid, barcode2tf$Sequence)],
                          dominant = apply(cd.ampl, 2, max), total = colSums(cd.ampl))
# annotate perturbation
cell.map.tf$id <- FALSE
cell.map.tf$id[cell.map.tf$dominant > 10 & cell.map.tf$dominant > cell.map.tf$total * 0.9] <- TRUE
TFn <- sort(table(subset(cell.map.tf[rownames(p2.tfinvitro$counts), ], id)$TF), decreasing = T)

pl <- ggplot(cell.map.tf, aes(total, dominant, col = id))  + geom_point(size = I(1)) + theme_bw() +
  scale_color_manual(values = c('TRUE' = 'black', 'FALSE' = 'grey')) + theme(text = element_text(size = 20))
pl
#
pdf('figures/invitro/TFbarcode_stat.pdf')
pl
dev.off()

### upload DNA barcode abundances and compare with RNA-based estimates
dna.tf <- read.table('data/raw/3_Count_161_HSC_filt.txt', header = T)
rna.tf <- table(cell.map.tf$TF[cell.map.tf$id])
pl <- ggplot(,aes(dna.tf$HSC / sum(dna.tf$HSC), rna.tf[dna.tf$TF] / sum(rna.tf[dna.tf$TF], na.rm = T), label = dna.tf$TF)) +
  scale_x_log10() + scale_y_log10() + theme_bw() + theme(text = element_text(size = 20)) +
  geom_point(size = 2) + geom_text_repel() + geom_abline(linetype = 2) + labs(title = 'in vitro TF frequencies', x = 'frequency, DNA', y = 'frequency, RNA')
pl
#
pdf('figures/invitro/comparison_rna_dna.pdf')
pl
dev.off()

##
png('figures/invitro/TF_emb_distribution.png', width = 800, height = 800)
par(mfrow=c(4,4), mar = rep(0.1, 4))
for (tfnm in names(TFn)[1:16]){
  print(tfnm)
  plot(p2.tfinvitro$embeddings$PCA$UMAP, pch = 19, cex = 1, col = 'grey', xaxt = 'n', yaxt = 'n')
  cls <- cell.map.tf$cell[cell.map.tf$cell %in% rownames(p2.tfinvitro$counts) & cell.map.tf$TF == tfnm & cell.map.tf$id]
  points(p2.tfinvitro$embeddings$PCA$UMAP[cls,], pch = 19, cex = 1, col = 'red')
  text(4.5, c(4, 3), labels = c(tfnm, paste(length(cls), 'cells')), cex = 3)
}
dev.off()

pdf('figures/invitro/TF_cells.pdf')
par(mfrow=c(1,1), mar = c(5,5,0.1, 0.1))
barplot(TFn[TFn > 2], #col = ifelse(p.adjust(pv[TFn > 2], 'BH') < 0.5, 'grey30', 'grey'),
        horiz = T, xlab = 'cells', las = 2, cex.names = 1, log = 'x', cex.axis = 2, cex.lab = 2)
dev.off()


# show major TFs
cls <- intersect(rownames(p2.tfinvitro$embeddings$PCA$UMAP), rownames(cell.map.tf))
tfids <- cell.map.tf[cls,]$TF; names(tfids) <- cls
tfids[!cell.map.tf[cls,]$id] <- NA
tfids[!tfids %in% names(sort(table(tfids), decreasing = T))[2:10]] <- NA
tfids <- tfids[!is.na(tfids)]
tfids <- factor(tfids, levels = names(sort(table(tfids), decreasing = T)) )
#cols <- rainbow(length(unique(tfids))); names(cols) <- names(unique(tfids))
cols <- brewer.pal(length(unique(tfids)), "Set1"); names(cols) <- names(unique(tfids))
pl <- ggplot() + geom_point(aes(p2.tfinvitro$embeddings$PCA$UMAP[,1], p2.tfinvitro$embeddings$PCA$UMAP[,2]), col = 'grey85', size = I(1)) +
  geom_point(aes(p2.tfinvitro$embeddings$PCA$UMAP[names(tfids),1], p2.tfinvitro$embeddings$PCA$UMAP[names(tfids),2],
                 col = tfids), size = I(1)) +
  labs(x = '', y = '') + theme_bw() +
  scale_color_manual('Perturbation', values = cols) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text  = element_text(size = 20), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), panel.grid  = element_blank(),
        legend.title  = element_text(size = 20))
pl
#
pdf('figures/invitro/TF_major_emb.pdf',
    width = 5, height = 3)
pl
dev.off()

### estimate TF-specific cell state biases
tf.clust <- do.call(rbind, lapply(names(TFn), function(tf){
  cells <- cell.map.tf$cell[cell.map.tf$id & cell.map.tf$TF == tf & cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
  table(p2.tfinvitro$clusters$PCA$leiden[cells])
}))
rownames(tf.clust) <- names(TFn)
cls <- cell.map.tf$cell[cell.map.tf$id & cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
cls.gfp <- cell.map.tf$cell[cell.map.tf$id & cell.map.tf$TF == 'GFP' & cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
qplot(table(p2.tfinvitro$clusters$PCA$leiden), table(p2.tfinvitro$clusters$PCA$leiden[cls]))
qplot(table(p2.tfinvitro$clusters$PCA$leiden[cls]), table(p2.tfinvitro$clusters$PCA$leiden[cls.gfp]), )

# perform differential cluster abundance testing
p <- tf.clust['GFP',] / sum(tf.clust['GFP',])
pv <- apply(tf.clust, 1, function(x1){
  #tst <- chisq.test(x = x1, p = p)
  #tst <- chisq.test(x = x1, y = tf.clust['GFP',] )
  #tst <- fisher.test(x = x1, y = tf.clust['GFP',])
  tst <- fisher.test( cbind(c(sum(x1[c(1,4,5)]), sum(x1)),
                            c(sum(tf.clust['GFP',c(1,4,5)]), sum(tf.clust['GFP',])) ))

  tst$p.value
})
sort(pv, decreasing = T)
sort(p.adjust(pv, 'BH'), decreasing = T)


### examine expression levels of overexpressed TFs
cls.gfp <- cell.map.tf$cell[cell.map.tf$id & cell.map.tf$TF == 'GFP' &
                              cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
tf_overexpr <- do.call(rbind, lapply(setdiff(names(TFn), 'GFP'), function(tf){
  print(tf)
  cls <- cell.map.tf$cell[cell.map.tf$id & cell.map.tf$TF == tf &
                            cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
  str(cls)
  tst <- wilcox.test(p2.tfinvitro$counts[cls,tf], p2.tfinvitro$counts[cls.gfp,tf])
  c(mean(p2.tfinvitro$counts[cls,tf]), mean(p2.tfinvitro$counts[cls.gfp,tf]),
    tst$p.value)
}))
rownames(tf_overexpr) <- setdiff(names(TFn), 'GFP')
colnames(tf_overexpr) <- c('TF', 'GFP', 'p')
tf_overexpr <- data.frame(tf_overexpr)
tf_overexpr$p[is.na(tf_overexpr$p)] <- FALSE
pp <- ggplot(tf_overexpr, aes(GFP, TF, label = rownames(tf_overexpr), size = -log(p)+10, fill = p.adjust(p, 'BH') < 0.01 )) +
  guides(size = FALSE, fill=guide_legend(title = "FDR < 0.01")) +
  geom_text_repel(size = I(3)) +
  geom_point(shape = 21) + geom_abline(linetype = 2) + theme_bw() +
  labs(title = 'TF abundance', x = 'GFP-overexpressing cells', y = 'TF-overexpressing cells') +
  scale_x_log10(limits = c(1e-4, 1)) + scale_y_log10(limits = c(1e-4, 1))  +
  theme(legend.position = c(.8, .2))#, legend.title  = element_text('FDR < 0.1'))
pp
#
pdf('figures/invitro/TF_GFP_overexpression.pdf',
    width = 4, height = 4)
print(pp)
dev.off()


save(cd.tfvitro, cell.map.tf, TFn, file = 'data/processed/TF_invitro.RData')

write.table(env.tf$cell.map.tf, file = 'data/processed/TF_overexpression_id.txt')
