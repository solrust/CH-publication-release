source('R/project_gset.R')
library(ggrepel)
library(parallel)
library(Matrix)
library(RColorBrewer)
library(cacoa) # only for state enrichment

#
load(file = 'data/processed/scrna.RData')
#
hash.annot <- readRDS(file = "data/processed/hashtag_info.rds")
hash.match <- readRDS(file = "data/processed/hashtag_meta.rds")
#
#p2filter <- readRDS(file = 'data/p2objects/p2filter.rds')


#
emb.clust <- do.call(rbind, lapply(unique(clust.final), function(i){
  cls <- names(clust.final)[!is.na(clust.final) & clust.final == i]
  colMeans(emb[cls,])
}))
rownames(emb.clust) <- unique(clust.final)

#
pdf("figures/emb.pdf", width = 5, height = 5)
par(mfrow=c(1,1), mar = rep(0.2, 4))
plot(emb, cex = 0.3, col = adjustcolor('grey90', 0.2), , xaxt = 'n', yaxt = 'n', bty = 'n')
points(emb[names(clust.final),], cex = 0.3,
       col = adjustcolor(clust.final.cols[clust.final], 0.2))
points(emb.clust, cex = 3, pch = 21, lwd = 1, bg = 'grey95', col = 'black')
text(emb.clust[,1], emb.clust[,2], rownames(emb.clust), font = 2,
     col = clust.final.cols[as.numeric(rownames(emb.clust))])
dev.off()

#
emb.clust.leiden <- do.call(rbind, lapply(unique(p2clusts$PCAncc$leiden), function(i){
  cls <- names(p2clusts$PCAncc$leiden)[!is.na(p2clusts$PCAncc$leiden) & p2clusts$PCAncc$leiden == i]
  colMeans(emb[cls,])
}))
rownames(emb.clust.leiden) <- unique(p2clusts$PCAncc$leiden)
pdf("figures/emb_leiden.pdf", width = 5, height = 5)
par(mfrow=c(1,1), mar = rep(0.2, 4))
cols <- as.numeric(p2clusts$PCAncc$leiden)
plot(emb, cex = 0.5, col = adjustcolor(cols, 0.2), , xaxt = 'n', yaxt = 'n', bty = 'n')
points(emb.clust.leiden, cex = 3, pch = 21, lwd = 1, bg = 'grey95', col = 'black')
text(emb.clust.leiden[,1], emb.clust.leiden[,2], rownames(emb.clust.leiden), font = 2,
     col = as.numeric(unique(p2clusts$PCAncc$leiden)))
dev.off()

png("figures/emb_lymphoid.png")
par(mfrow=c(1,1), mar = rep(0.2, 4))
alpp <- ifelse(clust.final %in% 4, 1, 0.01)
plot(emb[names(clust.final),],  , xaxt = 'n', yaxt = 'n', bty = 'n', cex = 0.3,
     col = colorspace::adjust_transparency(clust.final.cols[clust.final], alpp))
dev.off()

#
png("figures/emb_my.png")
par(mfrow=c(1,1), mar = rep(0.2, 4))
alpp <- ifelse(clust.final %in% 8, 1, 0.01)
plot(emb[names(clust.final),],  , xaxt = 'n', yaxt = 'n', bty = 'n', cex = 0.3,
     col = colorspace::adjust_transparency(clust.final.cols[clust.final], alpp))
dev.off()

png("figures/emb_hsc.png")
par(mfrow=c(1,1), mar = rep(0.2, 4))
alpp <- ifelse(clust.final %in% 1, 1, 0.01)
plot(emb[names(clust.final),],  , xaxt = 'n', yaxt = 'n', bty = 'n', cex = 0.3,
     col = colorspace::adjust_transparency(clust.final.cols[clust.final], alpp))
dev.off()

png("figures/emb_geno_batch.png", width = 300*2, height = 300)
par(mfrow = c(1,2), mar = c(2,2,2,2))
geno <- hash.match[sampleid[rownames(p2counts)],]$genotype
names(geno) <- rownames(p2counts)
smpl <- sample(rownames(emb))
cols <- rep('black', nrow(emb))
cols[geno[rownames(emb)][smpl] == 'DNMT3A'] <- 'red'
cols[geno[rownames(emb)][smpl] == 'TET2'] <- 'green'
plot(emb[smpl,],col = adjustcolor(cols, 0.8), #bg = 'white',
     lwd = 1, main = 'genotypes', cex.main = 2, font.main = 1,
     xlab = '', ylab = '',
     pch = 21, cex = 0.1, xaxt = 'n', yaxt = 'n')#, bty = 'n')
legend('topright',  c('WT', 'DNMT3A', 'TET2'), col = c('black', 'red', 'green'),
     pch = 19, bty = 'n', pt.cex = 2, cex = 1.5)
###
plot(emb[smpl,],col = as.factor(sampleid[smpl]), #bg = 'white',
     lwd = 1,  main = 'samples', cex.main = 2, font.main = 1,
     xlab = '', ylab = '',
     pch = 21, cex = 0.1, xaxt = 'n', yaxt = 'n')#, bty = 'n')
dev.off()

### show cell type proportions across mice:
clust <- as.factor(clust.final)
nms <- names(clust)
clust.info <- data.frame(clust = clust,
                         sample = sampleid[nms],
                         genotype = hash.match[sampleid[nms],]$genotype,
                         exp = hash.match[sampleid[nms],]$exp,
                         count = 1
)
clust.aggr <- aggregate(count ~ clust + sample + genotype + exp, data = clust.info, sum)
#
clust.aggr.tot <- aggregate(count ~ sample, data = clust.aggr, sum)
rownames(clust.aggr.tot) <- clust.aggr.tot$sample
#
clust.aggr.clust <- aggregate(count ~ clust + exp, data = clust.aggr, sum)
clust.aggr.clust.exp <- aggregate(count ~ exp, data = clust.aggr.clust, sum)
rownames(clust.aggr.clust.exp) <- clust.aggr.clust.exp$exp
clust.aggr.clust$freq <- clust.aggr.clust$count / clust.aggr.clust.exp[clust.aggr.clust$exp,]$count
#rownames(clust.aggr.clust) <- clust.aggr.clust$clust
#
clust.aggr$freq <- clust.aggr$count / clust.aggr.tot[clust.aggr$sample,]$count
mtch <- match(paste(clust.aggr$clust, clust.aggr$exp, sep = '_'),
              paste(clust.aggr.clust$clust, clust.aggr.clust$exp, sep = '_'),
)
clust.aggr$freq <- clust.aggr$freq / clust.aggr.clust[mtch,]$freq

# testing differences in proportions
pv.prop <- do.call(cbind, lapply(unique(clust.aggr$exp), function(expp){
  do.call(rbind,lapply(unique(clust.aggr$clust), function(cl){
    ind <- clust.aggr$exp == expp & clust.aggr$clust == cl
    s1 <- clust.aggr[ind & clust.aggr$genotype == 'WT',]$freq
    s2 <- clust.aggr[ind & clust.aggr$genotype != 'WT',]$freq
    tst <- t.test(s1, s2)
    c(mean(s2) / mean(s1), tst$p.value  )
  }))
}))
rownames(pv.prop) <- unique(clust.aggr$clust)
colnames(pv.prop) <- c('fc.DNMT3A', 'p.DNMT3A', 'fc.TET2', 'p.TET2')
pv.prop <- data.frame(pv.prop)
pl <- ggplot(pv.prop, aes(fc.DNMT3A, fc.TET2, label = unique(clust.aggr$clust))) +
  geom_point(size = 5,  color = clust.final.cols[unique(clust.aggr$clust)]) +
  geom_text_repel(size = 10) +
  theme_bw() + labs(x = 'DNMT3A', y = 'TET2', title = 'cell state enrichment')+
  annotate("text", x = min(pv.prop$fc.DNMT3A), y = min(pv.prop$fc.TET2), size = 6,
           label = paste("\u03C1 = ", round(cor(pv.prop[,1], pv.prop[,3]), 2),
                         '\np < ', round(cor.test(pv.prop[,1], pv.prop[,3])$p.value, 2),sep = ''), hjust = 0, vjust = 0) +
  theme(text = element_text(size = 20), legend.position = 'none')
pl
cor(pv.prop[,1], pv.prop[,3])
cor.test(pv.prop[,1], pv.prop[,3])
pdf("figures/cell_state_proportion_scatter.pdf", width = 4, height = 4)
pl
dev.off()

# cacoa-based state enrichment:
cacoa.clust.p <- do.call(cbind, lapply(c('123', '127'), function(ind.exp){
  print(ind.exp)
  sample.groups <- hash.match[hash.match$exp == ind.exp,]$genotype
  names(sample.groups) <- rownames(hash.match)[hash.match$exp == ind.exp]
  #
  cls <- names(clust.final)
  cls <- cls[sampleid[cls] %in% names(sample.groups)]
  cell.groups <- clust.final[cls] # p2filter$clusters$PCAncc$leiden[cls]
  sample.per.cell <- sampleid[cls]
  #
  cao <- Cacoa$new(Matrix(cd[, cls], sparse = T),
                   #t(p2filter$counts)[,cls],  # send raw or normalized matrix?
                   sample.groups = sample.groups,     # +
                   cell.groups = cell.groups,         # +
                   sample.per.cell = sample.per.cell, # +
                   ref.level = "WT",             # +
                   target.level = ifelse(ind.exp == '123', "DNMT3A", "TET2"),       # +
                   embedding = emb[cls,],
                   n.cores = 30, verbose = T)             # +
  cao$estimateCellLoadings()
  cao$test.results$coda$pval[as.character(1:10)]
}))
print(cacoa.clust.p)

pl <- ggplot(clust.aggr,
             aes(x = clust,
                 y = freq, col = genotype)) +
  xlab('clusters') + ylab('cell state enrichment') + theme_bw() +
  scale_y_log10(limits = c(1/2, 2)) +
  scale_color_manual(values = col.genotype) +
  theme(text = element_text(size = 20), legend.title = element_text(size = 0)) +
  geom_boxplot(width = 0.8, outlier.shape=NA, alpha = 1) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             size = 2, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~exp, nrow = 2) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())#+ coord_flip()
pl
#
pdf("figures/cell_state_proportions.pdf", width = 6, height = 4)
pl
dev.off()


### estimate smoothed state fold enrichment for every cell.
enrich.geno <- lapply(c('123', '127'), function(expid){
  print(expid)
  #expid <- '123'
  smpl <- rownames(hash.match[hash.match$exp == expid & hash.match$genotype != 'WT',])
  cls <- intersect(names(sampleid)[sampleid %in% smpl], names(clust.final))
  str(cls)
  #
  i <- 1; if (expid == '127') i <- 3
  x <- log2(pv.prop[clust.final[cls], i])
  #x <- as.integer(x / max(abs(x)) * 50) + 51
  x <- as.integer(x * 50) + 51
  x[x < 1] <- 1; x[x > 100] <- 100
  names(x) <- cls
  y <- unlist(mclapply(cls, function(cl){
    y <- rowSums( (t(t(emb[cls,]) - emb[cl,]))^2)
    clss <- names(y)[rank(y) < 30]
    mean(x[clss])
  }, mc.cores = 20))
  cbind(x, y)
})

png("figures/cell_state_proportion_emb.png", width = 400, height = 800)
par(mfrow=c(2,1), mar = rep(0.5,4))
for(i in 1:2){
  cls <- rownames(enrich.geno[[i]])
  cols <- colorRampPalette(c("blue","grey90", "red" ))(n = 101)[enrich.geno[[i]][,2]]
  plot(emb, col  = 'grey95', cex = 0.05, pch = 19, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  points(emb[cls,], col  = cols, cex = 0.5, pch = 19)
  legend('topright', legend = ifelse(i == 1, 'DNMT3A', 'TET2'), bty = 'n', cex = 4)
}
dev.off()

png("figures/col_heatmap.png", width = 300, height = 300 / 2.5)
par(mfrow=c(1,1), mar = c(5.5,2,0.5,2), mgp = c(3,1.0,0))
image(as.matrix(seq(0.5, 2, length.out = 100)), col = colorRampPalette(c("blue","grey90", "red" ))(n = 100),
      xaxt = 'n',
      breaks = seq(0.5, 2, length.out = 101),
      #breaks = c(0.5, 1, 2),
      yaxt = 'n',# border = 'black', lwd = 5,
      xlab = 'state fold change', cex.lab = 3)
axis(1, at = c(0, 0.5, 1), labels = c(0.5, 1, 2), las = 1, cex.axis = 2)
dev.off()


### estimate cluster-specific count matrix
clust.counts <- do.call(cbind, lapply(1:length(unique(clust.final)), function(cl){
  print(cl)
  cls <- names(clust.final)[clust.final == cl]
  Matrix::colMeans(p2counts[cls,])
}))
rownames(clust.counts) <- colnames(p2counts)
colnames(clust.counts) <- 1:length(unique(clust.final))
clust.counts.norm <- clust.counts / rowMeans(clust.counts)
# choose DE genes
de.genes <- 10 #5 - for heatmap, 10 - for Sup Table
gns <- unlist(lapply(1:length(unique(clust.final)), function(cl){
  print(cl)
  de <- diffGenes$PCAncc$final[[cl]]
  de <- de[order(-de$M),]
  de <- de[de$Z > 4,]
  return(rownames(de)[1:min(de.genes, nrow(de))])
}))
gns <- unique(gns)
str(gns)
coldat <- data.frame(clusts = 1:ncol(clust.counts.norm))
pheat <- pheatmap::pheatmap( log2(clust.counts.norm[gns,] + 1e-2), #breaks = seq()
                         show_colnames = T, cluster_cols = T,
                         annotation_col = coldat, annotation_colors = list(clusts = clust.final.cols), angle_col = 0,
                         annotation_names_col = F, annotation_legend = F, treeheight_col = 10, treeheight_row = 10,
                         labels_row = as.expression(genes.ital[gns]),
                         breaks = seq(-4, 4, length.out = 100),
                         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(n = 100),
                         fontsize = 12
)
#
pdf("figures/markers_heatmap.pdf", width = 4, height = 8*1)
pheat
dev.off()
#
marker.data <- data.frame(marker.cluster = apply(clust.counts.norm[gns,], 1, which.max), clust.counts.norm[gns,])
colnames(marker.data)[-1] <- paste('expr.cluster', 1:(ncol(marker.data)-1), sep = '.')
write.xlsx(marker.data, '../resource/marker_genes.xlsx', row.names = TRUE, col.names = TRUE)



### show activity of annotate programs on top of embedding
show_gene <- FALSE
progs_curate <- readRDS("data/genesets/processed/fate_markers.rds")
#pdf("/home/solrust/CH/draft/figures/rna_analysis/progs_activity_new.pdf", width = 4 * 4, height = 4 * 4)
png(paste("figures/progs_", ifelse(show_gene, 'gene', 'activity'),".png", sep = ''), width = 200 * 6, height = 200 * 6)
par(mfrow=c(4,4), mar = rep(0.5, 4))
pp <- 0.3
gnout <- vector()
for (k in setdiff(1:length(progs_curate), 13)){
  print(k)
  inf1 <- project_gset(gset = progs_curate[[k]], mat = p2counts, pr = NULL)
  if (!is.null(inf1)){
    if (!show_gene){
      x <- inf1$av.fit
      #x <- sign(x) * log(abs(x) * 1e+8 + 1)^1
      nm <- names(progs_curate)[k]
    }else{
      #nm <- gn <- names(which.max(abs(inf1$gscore)))
      #nm <- gn <- names(which.max(abs(inf1$gscore[setdiff(names(inf1$gscore), gnout)] )))
      nm <- gn <- names(sort(abs(inf1$gscore[setdiff(names(inf1$gscore), gnout)] ), decreasing = T))[1:2]
      gnout <- c(gnout, gn)
      #nm <- gn <- inf1$gscore[intersect(names(inf1$gscore), gns)]
      x <- p2init$counts[, gn]
    }
    if (!is.null(nrow(x))){
      x <- x[intersect(rownames(x), rownames(emb)),]
    }else{
      x <- x[intersect(names(x), rownames(emb))]
    }
    if (!show_gene){
      x <- sign(x) * (abs(x) * 1e+8 + 1)^pp
      #x <- sign(x) * log(abs(x) * 1e+8 + 1)^1
    }

    #for (i in 1:ncol(x)){
    cols <- colorRampPalette(c("grey90","grey80", "red" ))(n = 20)[as.numeric(cut(x,20))]
    #cols <- colorRampPalette(c("grey90", "red" ))(n = 20)[as.numeric(cut(x,20))]
    plot(emb[names(x[]),], #main = names(progs_curate)[k],
         #bg = cols, col = 'grey30',
         col = cols,
         #lwd = 0.05, #p2$clusters$cleaned$leiden1[cls],
         xlab = '', ylab = '', pch = 21, cex = 0.2, xaxt = 'n', yaxt = 'n')#, bty = 'n')
    legend('topright', nm[], bty = 'n', cex = 2.5)
    gns.show <- names(sort(abs(inf1$gscore), decreasing = T))[1:5]

    gns.show <- unlist(lapply(gns.show, function(gn){
      paste(substr(gn, 1, 1), tolower(substr(gn, 2, nchar(gn))), sep = '')
    }))

    unlist(lapply(gns.show, function(gn) substitute(paste(italic(gn)))))
    legend('bottomright', legend = gns.show, bty = 'n', cex = 3)

    #legend('bottomright', legend = substitute(paste(italic(gns.show))), bty = 'n', cex = 2)
    #substitute(paste(italic('this text is italic')))

    #}
  }
}
dev.off()

