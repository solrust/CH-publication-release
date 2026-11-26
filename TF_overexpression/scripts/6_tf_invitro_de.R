source('../RNA_analysis/R/draw_go.R')
library(Matrix)
library(pagoda2)
library(ggplot2)
library(ggrepel)
library(parallel)
library(scde)
library(RColorBrewer)
library(DESeq2)
library(data.table)
library(fgsea)
load(file = 'data/processed/TF_invitro.RData')
load(file = '../Association_analysis/data/processed/use/association_summary.RData')
load(file = '../RNA_analysis/data/processed/mutant_de.RData')
p2.tfinvitro <- readRDS(file = 'data/processed/p2_tfinvitro.rds')
gene_sets <- readRDS('../RNA_analysis/data/genesets/processed/gene_sets_use.rds')

### estimate TF-specific DE across cell states
tf.clust.list <- lapply(1:length(unique(p2.tfinvitro$clusters$PCA$leiden)), function(cl){
  print(cl)

  xx <- do.call(cbind, lapply(c('CTL', names(TFn)), function(tf){
    if (tf == 'CTL'){
      cls <- cell.map.tf$cell[!cell.map.tf$id & cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
    }else{
      cls <- cell.map.tf$cell[cell.map.tf$id & cell.map.tf$TF == tf &
                                cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
    }
    cls1 <- cls[p2.tfinvitro$clusters$PCA$leiden[cls] == cl]
    if (length(cls1) > 0){
      rowSums(cd.tfvitro[,cls1, drop = F])
    }else{
      rep(0, nrow(cd.tfvitro))
    }
  }))
  colnames(xx) <- paste(c('CTL', names(TFn)), cl, sep = '_')
  xx
})
names(tf.clust.list) <- 1:length(tf.clust.list)
inds.include <- setdiff(1:length(tf.clust.list), c(2,3,4,10)) # 10 - mito, 3,4 - cycling, 2 - partially cycling
tf.clust.mat <- do.call(cbind, tf.clust.list[inds.include])
tf.features <-  do.call(rbind, lapply(inds.include, function(i){
  #xx <- cbind(names(TFn)[i], 1:ncol(tf.clust.list[[i]]))
  xx <- cbind(c('CTL', names(TFn)), i)
  rownames(xx) <- colnames(tf.clust.list[[i]])
  xx
}))
tf.features <- data.frame(tf.features)
colnames(tf.features) <- c('TF', 'cluster')
tf.features$TF <- factor(tf.features$TF, levels = unique(tf.features$TF))

### run Deseq2
retain.detected <- colSums(tf.clust.mat > 0) > 100
table(retain.detected)
retain.gns <- rowSums(cd.tfvitro > 0) / ncol(cd.tfvitro) > 0.05
table(retain.gns)
dds_tf <- DESeqDataSetFromMatrix(countData = tf.clust.mat[retain.gns, retain.detected],
                                 colData = tf.features[retain.detected,],
                                 design = ~ TF + cluster)
dds_tf <- DESeq(dds_tf)
resultsNames(dds_tf)
res_tf <- results(dds_tf, name = "TF_GFP_vs_CTL")
res_tf <- results(dds_tf, name = "TF_ID1_vs_CTL")
head(data.frame(res_tf[order(res_tf$pvalue),]), 20)

### estimate TF-specific DE across cell states
tf.clust.list <- lapply(1:length(unique(p2.tfinvitro$clusters$PCA$leiden)), function(cl){
  print(cl)

  xx <- do.call(cbind, lapply(c('CTL', names(TFn)), function(tf){
    if (tf == 'CTL'){
      cls <- cell.map.tf$cell[!cell.map.tf$id & cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
    }else{
      cls <- cell.map.tf$cell[cell.map.tf$id & cell.map.tf$TF == tf &
                                cell.map.tf$cell %in% names(p2.tfinvitro$clusters$PCA$leiden)]
    }
    cls1 <- cls[p2.tfinvitro$clusters$PCA$leiden[cls] == cl]
    if (length(cls1) > 0){
      rowSums(cd.tfvitro[,cls1, drop = F])
    }else{
      rep(0, nrow(cd.tfvitro))
    }
  }))
  colnames(xx) <- paste(c('CTL', names(TFn)), cl, sep = '_')
  xx
})
names(tf.clust.list) <- 1:length(tf.clust.list)
inds.include <- setdiff(1:length(tf.clust.list), c(2,3,4,10)) # 10 - mito, 3,4 - cycling, 2 - partially cycling
tf.clust.mat <- do.call(cbind, tf.clust.list[inds.include])
tf.features <-  do.call(rbind, lapply(inds.include, function(i){
  #xx <- cbind(names(TFn)[i], 1:ncol(tf.clust.list[[i]]))
  xx <- cbind(c('CTL', names(TFn)), i)
  rownames(xx) <- colnames(tf.clust.list[[i]])
  xx
}))
tf.features <- data.frame(tf.features)
colnames(tf.features) <- c('TF', 'cluster')
tf.features$TF <- factor(tf.features$TF, levels = unique(tf.features$TF))

## record logFC and p-values
tf.in <- names(TFn)[unlist(lapply(names(TFn), function(nm) sum(grepl(nm, resultsNames(dds_tf))) )) > 0]
logFC.tf <- do.call(cbind, mclapply(tf.in, function(tf){
  print(tf)
  nm <- resultsNames(dds_tf)[grepl(tf, resultsNames(dds_tf))]
  #if (length(nm) > 0){
    #res_tf <- results(dds_tf, name = nm)
    res_tf <- lfcShrink(dds_tf, coef = nm, type="normal")
    res_tf$log2FoldChange
  #}
}, mc.cores = 20))
rownames(logFC.tf) <- rownames(tf.clust.mat)[retain.gns]
colnames(logFC.tf) <- tf.in
#
pv.tf <- do.call(cbind, lapply(tf.in, function(tf){
  print(tf)
  nm <- resultsNames(dds_tf)[grepl(tf, resultsNames(dds_tf))]
  #if (length(nm) > 0){
    res_tf <- results(dds_tf, name = nm)
    res_tf$pvalue
  #}
}))
rownames(pv.tf) <- rownames(tf.clust.mat)[retain.gns]
colnames(pv.tf) <- tf.in
compos.tf <- sign(logFC.tf) * -log(pv.tf)
padj.tf <- apply(pv.tf, 2, function(x) p.adjust(x, 'BH'))
ind <- rowSums(pv.tf < 1e-3) > 0
table(ind)
#save(dds_tf, logFC.tf, pv.tf, file = 'data/processed/TF_invitro_DE.RData')

###
load(file = 'data/processed/TF_invitro_DE.RData')

##
ncell <- 4
tfx <- colnames(logFC.tf)[TFn[colnames(logFC.tf)] >= ncell]
sum(TFn[colnames(logFC.tf)] >= ncell)
ind <- rowSums(padj.tf[,tfx] < 1e-1) > 0 & !rownames(padj.tf) %in% tfx
table(ind)
pheat <- pheatmap::pheatmap(t(logFC.tf[ind,  tfx]),
                            #labels_row = gnms,
                            show_colnames = F,
                            annotation_col_title = 'dsfg',
                            border_color = NA,
                            treeheight_col = 10, treeheight_row = 20,
                            color = colorRampPalette(c("blue", "white", "darkred" ))(n = 100),
                            breaks = seq(-3, 3, length.out = 100))
#
pdf('figures/invitro/TF_response_genes.pdf', width = 4, height = 4)
pheat
dev.off()



### correlation of TF responses
cr <- cor((logFC.tf[ind, tfx]), method = 'spearman')
#cr <- cor(compos.tf[ind,], method = 'spearman')
diag(cr) <- NA
pheat <- pheatmap::pheatmap(cr, main = 'TF overexpression response', color = colorRampPalette(c("darkgreen", "white", "darkorange" ))(n = 100),
                            breaks = seq(-0.5, 0.5, length.out = 100)
                            )
#
pdf('figures/invitro/TF_response_cor.pdf')
pheat
dev.off()


fdr.tf <- apply(pv.tf, 2, function(x) p.adjust(x, 'BH'))
logFC.tf1 <- melt(logFC.tf[,TFn[colnames(logFC.tf)] >= ncell])
pv.tf1 <- melt(pv.tf[,TFn[colnames(logFC.tf)] >= ncell])
fdr.tf1 <- melt(fdr.tf[,TFn[colnames(logFC.tf)] >= ncell])
#
xx <- data.frame(logFC.tf1, pv = pv.tf1$value, fdr = fdr.tf1$value)
colnames(xx) <- c('gene', 'TF', 'logFC', 'p', 'fdr')
#xx$gene[xx$p > 1e-4] <- ''
#xx  <- xx[sample( 1:nrow(xx), 10000), ]
xx$gene <- as.character(xx$gene)
xx$gene[xx$fdr > 2e-1] <- ''
pl <- ggplot(xx,aes(logFC, p, label = gene, col = fdr < 0.3)) + geom_point(size = I(0.2), alpha = I(1)) +
  geom_text_repel(max.overlaps = 15) +
  theme_bw() + scale_y_log10() + labs(x = 'logFC', y = 'p-value, log') +
  scale_color_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme(text = element_text(size = 20), legend.position  = 'none') +
  facet_wrap(~TF,scales='free', ncol = 5)
pl
#
pdf('figures/invitro/TF_response_scatterplots.pdf', width = 14, height = 14 )
pl
dev.off()

### estimate pathway enrichment
gsea.tfover.compose <- mclapply(which(TFn[colnames(logFC.tf)] >= 4) , function(i){
  sc <- compos.tf[,i]
  x <- fgsea(pathways = gene_sets, stats  = sc, minSize = 15, maxSize = 1000)
}, mc.cores = 10)
names(gsea.tfover.compose) <- names(which(TFn[colnames(logFC.tf)] >= 4))

### plot heatmap of enriched pathways
tf.enrich.go <- do.call(cbind, lapply(gsea.tfover.compose, function(x) -log(x$pval) * sign(x$NES) ))
rownames(tf.enrich.go) <- gsea.tfover.compose[[1]]$pathway
tf.enrich.go[is.na(tf.enrich.go)] <- 0
pvmax <- apply(abs(tf.enrich.go), 1, max)
ind1 <- pvmax > 9
table(ind1)
dist.tf.go <-  as.dist( 1-cor(t(tf.enrich.go[ind1,])) )
hclt <- hclust(dist.tf.go)
cl.tf.go <- data.frame(clust = as.factor(cutree(hclt,5)))
col.tf.go <- brewer.pal(length(unique(cl.tf.go$clust)), "Set2")
names(col.tf.go) <- unique(cl.tf.go$clust)
cl.tf.go$col <- col.tf.go[cl.tf.go$clust]
pheat <- pheatmap::pheatmap(t(tf.enrich.go[ind1,]), show_colnames = F, main = 'Pathway enrichment',
                   treeheight_col = 10, treeheight_row = 20,
                   clustering_distance_cols = dist.tf.go,
                   annotation_col = cl.tf.go[,1, drop = F],
                   annotation_colors = list('clust' = col.tf.go),
                   color = colorRampPalette(c("blue", "white", "red" ))(n = 100),
                   breaks = seq(-30, 30, length.out = 100))
#
pdf('figures/invitro/TF_response_GOheatmap.pdf', width = 5, height = 5)
pheat
dev.off()

pl <- draw_go(tf.enrich.go[ind1,], rownames(tf.enrich.go)[ind1],
              cl.tf.go[,1,drop=F], col.tf.go,
              n.terms = 10, only.go = FALSE, term.length = 70, cut.pval = 40)
pl
#
pdf(paste('figures/invitro/TF_response_GOexamples.pdf', sep = ''),
    width = 12, height = 14)
pl
dev.off()


### compare TF-induced changes with clone signatures
gns <- intersect(rownames(logFC.tf), rownames(gsea.clone$gscore.clone))
str(gns)
gns <- gns[rowSums(pv.tf[gns,] < 1e-3) > 0]
str(gns)
#xx <- cor(logFC.tf[gns,], gsea.clone$gscore.clone[gns,])
tf_response_pheno <- cor(compos.tf[gns,], gsea.clone$gscore.clone[gns,])
#tf_response_pheno <- cor(logFC.tf[gns,], gsea.clone$gscore.clone[gns,])
pheat <- pheatmap::pheatmap((tf_response_pheno), angle_col = 45, main = 'TF perturtbation response and clone phenotypes', fontsize = 10,
                            breaks = seq(-max(abs(tf_response_pheno)), max(abs(tf_response_pheno)), length.out = 100))
pdf('figures/invitro/TF_response_pheno_cor_v1.pdf', width = 5, height = 8)
pheat
dev.off()

#
mut_response <- do.call(cbind, lapply(de.mut2wt.shrink[[1]], function(x) sign(x$log2FoldChange) * -log(x$pvalue) ))
rownames(mut_response) <- rownames(de.mut2wt.shrink[[1]]$DNMT3A)
head(mut_response)
#
tf_response_mut <- cor(compos.tf[gns,], mut_response[gns,])
#
pheat <- pheatmap::pheatmap((tf_response_mut), angle_col = 45, main = 'TF perturtbation response and mutation distortion', fontsize = 10,
                   breaks = seq(-max(abs(tf_response_mut)), max(abs(tf_response_mut)), length.out = 100)
)
#
pdf('figures/invitro/TF_response_mutation_cor.pdf', width = 3, height = 8)
pheat
dev.off()

tf_response <- cbind(tf_response_pheno, tf_response_mut)
pheat <- pheatmap::pheatmap(tf_response, angle_col = 45, main = 'TF perturtbation response and signatures', fontsize = 10, cluster_cols = F, gaps_col = 7,
                   breaks = seq(-max(abs(tf_response)), max(abs(tf_response)), length.out = 100))
#
pdf('figures/invitro/TF_response_signatures.pdf', width = 6, height = 8)
pheat
dev.off()


