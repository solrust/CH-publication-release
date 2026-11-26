library(Matrix)
library(pagoda2)
library(ggplot2)
library(ggrepel)
library(parallel)
library(scde)
library(RColorBrewer)
load(file = 'data/processed/TF_invitro.RData')
p2.tfinvitro <- readRDS(file = 'data/processed/p2_tfinvitro.rds')

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
