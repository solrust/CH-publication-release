# TODO:
# define TF-target matrix
  # multiple TF-peak by peak-target matrices
# assess TF impact through DE of target genes.

# estimate TF-target matrix
tf2target0 <- (z.cor[,rownames(motif.matrix)] > 3)  %*% motif.matrix
tf2target1 <- (z.cor[,rownames(motif.matrix)] > 4)  %*% motif.matrix
tf2target2 <- (z.cor[,rownames(motif.matrix)] > 5)  %*% motif.matrix
tf2target3 <- (z.cor[,rownames(motif.matrix)] > 6)  %*% motif.matrix
# seurat-based TF-target matrix
tf2target.seurat <- g2p.seurat1 %*% motif.matrix[colnames(g2p.seurat1), ]
rownames(tf2target.seurat) <- rownames(g2p.seurat1)
# oracle TF-target matrix
library(arrow)
### annotate CellOracle TF-target relationships
xx <- read_parquet("../../hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet")
tf_targets_oracle <- lapply(xx[-c(1,2)], function(y) unique(xx[[2]][which(y > 0)]) )
names(tf_targets_oracle) <- names(xx[-c(1,2)])
gns <- as.factor(unique(unlist(tf_targets_oracle)))
tf2targets.oracle <- do.call(rbind, lapply(tf_targets_oracle, function(targs){
  table(gns[gns %in% targs])
}))
rownames(tf2targets.oracle) <- names(tf_targets_oracle)
tf2targets.oracle <- t(tf2targets.oracle)

# compare my and seurat list of significant tf-peaks
g2p.my <- (z.cor[,rownames(motif.matrix)] > 5) + 0
qplot(rowSums(g2p.seurat1), rowSums(g2p.my[rownames(g2p.seurat1),]))
cor(rowSums(g2p.seurat1), rowSums(g2p.my[rownames(g2p.seurat1),]))
qplot(colSums(g2p.seurat1), colSums(g2p.my[, colnames(g2p.seurat1)]))
cor(colSums(g2p.seurat1), colSums(g2p.my[, colnames(g2p.seurat1)]))

#
tf2target <- tf2target0
#tf2target <- tf2target.seurat
colnames(tf2target) <- unlist(motif.object@motif.names)[colnames(tf2target)]
#tf2target <- tf2targets.oracle
#
tf2target <- tf2target[rowSums(tf2target) > 0, ]
names(de_commitment)
cl <- '4'
scores <- de_commitment_old[[cl]]$log2FoldChange; names(scores) <- rownames(de_commitment_old[[cl]])
#scores <- de_commitment[[cl]]$log2FoldChange; names(scores) <- rownames(de_commitment[[cl]])
scores <- scores[intersect(names(scores), rownames(tf2target))]
#scores <- res_rna_dnmt3a$log2FoldChange; names(scores) <- rownames(res_rna_dnmt3a)
#scores <- res_rna_tet2$log2FoldChange; names(scores) <- rownames(res_rna_tet2)
tf.impact <- do.call(rbind, lapply(1:ncol(tf2target), function(i){
  print(i)
  if (i %% 10 == 0) print(i)
  #i <- 1
  inds <- which(tf2target[,i] > 0)
  gns <- intersect(rownames(tf2target)[which(tf2target[,i] > 0)], names(scores))
  if (length(gns) > 1){
  tst <- t.test(scores[gns], scores)
  shift <- mean(scores[gns] > 0, na.rm=T) /
           mean(scores > 0, na.rm=T)
  c(shift, tst$p.value)
  }else{
    c(NA, NA)
  }
}))
#rownames(tf.impact) <- unlist(motif.object@motif.names)[colnames(tf2target)]
rownames(tf.impact) <- colnames(tf2target)
colnames(tf.impact) <- c('shift', 'p')
tf.impact <- as.data.frame(tf.impact)
sz <- colSums(tf2target > 0); names(sz) <- rownames(tf.impact)
ggplot(tf.impact, aes(shift, p, #p.adjust(p, 'BH'),
                      label = rownames(tf.impact))) +
  geom_point() + geom_vline(xintercept = 1) +
  scale_y_log10() + geom_hline(yintercept = 0.05) + geom_text_repel()

