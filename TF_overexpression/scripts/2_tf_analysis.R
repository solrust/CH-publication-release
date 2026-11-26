library(ggplot2)
library(parallel)
library(ggrepel)
library(data.table)
library(gridExtra)
#
load('data/processed/TF_trends.RData')
load('../Association_analysis/data/processed/use/association_summary.RData')


## examine TF growth and TF expression/activity
tf.trnd <- tf.trends.no25.v1m
#tf.trnd <- tf.trends.v1m
tfexp <- rep(NA, length(tf.trnd$TF)); names(tfexp) <- tf.trnd$TF
gg <- intersect(tf.trnd$TF, rownames(gsea.clone$gscore.clone))
tfexp[gg] <- gsea.clone$gscore.clone[gg, 1]
#
xx <- data.frame(TFgrowth = tf.trnd$t,
                 TFgrowth.p = tf.trnd$p.t,
                 TF = tf.trnd$TF,
                 TFgrowth.expr = tfexp)

p1 <- ggplot(xx, aes(sign(TFgrowth) * -log(TFgrowth.p), TFgrowth.expr, label = TF)) +
  geom_point() + geom_text_repel(max.overlaps = 20) + theme_bw() + theme(text = element_text(size = 15)) +
  geom_hline(yintercept = c(-5, 5), linetype = 2) + geom_vline(xintercept = c(-5, 5), linetype = 2) +
  labs(x = 'TF overexpression effect (0-18w), -log P', y = 'TF expression association, -log P')

###
bias_site1 <- data.frame(bias_site,  TF = barcode2tf$Gene[match(rownames(bias_site), barcode2tf$Sequence)])
tfexp1 <- rep(NA, length(bias_site1$TF)); names(tfexp1) <- bias_site1$TF
gg <- intersect(bias_site1$TF, rownames(gsea.clone$gscore.clone))
tfexp1[gg] <- gsea.clone$gscore.clone[gg, phenm]
xx1 <- data.frame(TFbias = bias_site1[,1],
                  TFbias.p = bias_site1[,2],
                  TF = bias_site1$TF,
                  #bias.activity = tfact1,
                  bias.expr = tfexp1)
p2 <- ggplot(xx1, aes(sign(TFbias) * -log(TFbias.p), bias.expr, label = TF)) +
  geom_point() + geom_text_repel(max.overlaps = 20) + theme_bw() + theme(text = element_text(size = 15)) +
  geom_hline(yintercept = c(-5, 5), linetype = 2) + geom_vline(xintercept = c(-5, 5), linetype = 2) +
  labs(x = 'TF B/My bias, -log P', y = 'TF expression lymphoid bias, -log P')

p2

pdf('figures/TFpert_vs_expr.pdf', width = 10, height = 5)
print(grid.arrange(p1, p2, nrow = 1))
dev.off()

