library(ggplot2)
library(ggrepel)
library(gridExtra)
library(parallel)
library(data.table)
library(dplyr)
load('../RNA_analysis/data/processed/scrna.RData')
load('data/processed/clone_info.RData')
hash.annot <- readRDS(file = "../RNA_analysis/data/processed/hashtag_info.rds")
hash.match <- readRDS(file = "../RNA_analysis/data/processed/hashtag_meta.rds")
expids <- c('exp123', 'exp127')
smpls <- c('RNA1', 'RNA2', 'RNA3', 'RNA4')


clone.dat <- data.frame()
#clone.dat.size <- data.frame()
for (expid in expids){
  for (smpl in smpls){
    print(paste(expid, smpl))
    larry <-readRDS(paste("data/processed/larry_", expid, '_',smpl,".rds", sep = ''))
    cmax <- apply(larry$mat.clones, 2, function(x) max(x))
    cother <- colSums(larry$mat.clones) - apply(larry$mat.clones, 2, function(x) max(x))
    #
    clone.dat <- rbind(clone.dat,
                       data.frame(max = cmax, nonmax = cother, id = paste(expid, smpl, sep = '_'))
    )
    #
    #clone.dat.size <- rbind(clone.dat.size,
    #                        data.frame(size = as.numeric(table(larry$cloneid)), id = paste(expid, smpl, sep = '_'))
    #)
  }
}
#
pl <- ggplot(clone.dat[sample(1:nrow(clone.dat)),], aes(max + nonmax, max)) +
  geom_point(aes(col = id), alpha = I(0.01), size = I(2)) + labs(x = 'all clones, UMIs', y = 'dominant clone, UMIs') +
  lims(x = c(0, 250), y = c(0, 250)) +
  geom_abline(slope = 0.85, col = 'black', lty = 2) + geom_vline(xintercept = 5, col = 'black', lty = 2) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position =  'none')
pl
#
png("figures/dominant_clone_calling.png", width = 300, height = 300)#, width = 1.5 * 3, height = 1.5 * 2)
pl
dev.off()


### plot clone representation in mice
png("figures/clone_hashtag.png", width = 100 * 3, height = 100 * 2)
#pdf("figures/clone_hashtag.pdf", width = 1.5 * 3, height = 1.5 * 2)
cex <- 1.5
doublet.barcodes <- c()
par(mfrow=c(2, 3), mar = c(5,5,1.5,1.5), mgp = c(3,1,0))
for (i in setdiff(1:8, c(2, 7))
){
  hash <- hash.annot[[i]]
  if (i < 5){
    i1 <- i
    expid <- '123'
  }else{
    i1 <- i - 4
    expid <- '127'
  }

  expid1 <- paste('exp', expid, sep = '')
  smpl <- paste('RNA', i1, sep = '')
  larry <-readRDS(paste("data/processed/larry_", expid1, '_',smpl,".rds", sep = ''))

  mat.clones <- larry$mat.clones
  #colnames(mat.clones) <- paste(expid, i1,
  #                        unlist(lapply(strsplit(colnames(mat.clones), 'one_'), function(x) substr(x[2], 1, nchar(x[2])-2))), sep = '_')

  #names(hash) <- paste('one_', names(hash.annot[[i]]), '-1', sep = '')
  cls <- intersect(colnames(mat.clones), names(hash))
  str(cls)
  print(table(hash[cls]))
  #
  brd <-do.call(cbind, lapply(unique(hash[cls]), function(x){
    cls1 <- cls[hash[cls] == x]
    y <- rep(0, nrow(mat.clones))
    names(y) <- rownames(mat.clones)
    z <- table(larry$cloneid[cls1])
    y[names(z)] <- z
    y
  }))
  colnames(brd) <- unique(hash[cls])


  inds <- which(colnames(brd) %in% hash.match$mouse)
  plot(brd[, inds[1]] + 1, brd[, inds[2]] + 1, main = paste('RNA', i, sep = '_'),
       cex.lab = cex, cex.axis = cex, pch = 21, cex = 1, log = 'xy',
       xaxt = 'n', yaxt = 'n',
       bg = 'grey', col = 'black',
       xlab = paste('hashtag', colnames(brd)[inds[1]]), ylab = paste('hashtag', colnames(brd)[inds[2]]))
  axis(1, at = c(1, 1e+1 + 1, 1e+2 + 1), labels = c(0, 1e+1, 1e+2),
       cex.axis = cex)
  axis(2, at = c(1, 1e+1 + 1, 1e+2 + 1), labels = c(0, 1e+1, 1e+2),
       cex.axis = cex)

  indx <- brd[, inds[1]] > 0 & brd[, inds[2]] > 0
  print(paste('doublet barcodes:',sum(indx)))
  if (sum(indx) > 0){
    doublet.barcodes <- c(doublet.barcodes, names(which(indx)))
  }
}
dev.off()

### draw histograms of clones
x <- cbind(cell.clone[,-c(1, 3, 4)], geno = hash.match[cell.clone$id,]$genotype)
#xx <- data.table::dcast(x, id + geno + cloneid + expid ~ ..., length)
xx <- data.table::dcast(data.table::as.data.table(x),
                        id + geno + cloneid + expid ~ ..., length)
colnames(xx)[5] <- 'size'
pl <- ggplot(xx) +  geom_histogram(aes(size, fill = geno), lwd = 0.1, col = 'black', bins = 10) + scale_x_log10() +
  theme_bw() + theme(text = element_text(size = 20))+#, legend.position =  'none') +
  scale_fill_manual(values = col.genotype) +# scale_fill_manual(values = col.genotype)  +
  facet_wrap(~ id + expid, ncol = 2, dir = 'v')  +
  theme(strip.text = element_blank(),
        legend.title = element_blank(), legend.position = 'top')
pl
pdf("/home/solrust/CH/draft/figures/rna_clones/clone_hist.pdf")#, width = 1.5 * 3, height = 1.5 * 2)
pl
dev.off()

### show number of scRNA-seq clones
y <- table(xx[xx$size > 0,]$id)
yy <- data.frame(cbind(y, hash.match[names(y), c('genotype', 'exp')]))
cex <- 1.4
pdf('figures/clone_numbers.pdf', width = 3.5, height = 3.5)
par(mfrow=c(1,1), mar = c(5,5,1,1))
nums0 <- nums <- as.numeric(factor(yy$genotype, levels = c('DNMT3A', 'WT', 'TET2')))
ind1 <- yy$geno == 'WT' & yy$exp == '123'
nums[ind1] <- nums[ind1] - 0.2
ind2 <- yy$geno == 'WT'  & yy$exp == '127'
nums[ind2] <- nums[ind2] + 0.2
plot(nums, yy$Freq, cex = 2, cex.lab = cex, ylab = 'scRNA-seq clones',
     cex.axis = cex, xlab = '', bty = 'n',
     xlim = c(0.8, 3.2), ylim = c(0, max(yy$Freq)+1), bg = col.genotype[yy$geno], col = 'black',  pch = ifelse(yy$exp=='123', 21, 24), xaxt = 'n')
legend('bottomright', c('exp1', 'exp2'), pch = c(21, 24), cex = cex, bty = 'n')
axis(1, labels = c(as.factor(yy$geno)), at = nums0,
     cex.axis = cex)
dev.off()

# estimate p-values
ix <- hash.match[rownames(yy),]$exp %in% '127'
swt <- yy[ix & yy$geno == 'WT',]$Freq
smut <- yy[ix & yy$geno != 'WT',]$Freq
wilcox.test(swt, smut)
#t.test(swt, smut)
swt
smut


### plot all clones
rasterpdf::raster_pdf("figures/clones_patterns.pdf",width = 8, height = 8)
par(mfrow = c(1,1), mar = c(1,1,1,1))
for (cln in order(-clone.info$size)){
  #cl <- clone.info$clone
  #clsid <- names(cloneids)[cloneids==clid]
  cls <- rownames(cell.clone)[cell.clone$cloneid == rownames(clone.info)[cln]]
  cls <- intersect(cls, rownames(emb))
  if (length(cls) > 0){
    plot(emb, pch = 19, cex = 0.1, col = 'grey', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
    points(emb[cls,], pch = 19, cex = 0.1, col = 'blue', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
    legend('topright', paste('clone', cln,
                             #clone.info$clone[cln],
                             '\n',
                             'size:', clone.info$size[cln],
                             '\n',
                             'id:', clone.info$id[cln],
                             '\n',
                             'genotype:', clone.info$genotype[cln]), bty = 'n', cex = 1)
    #legend("bottomleft", paste('size =', clone.states[clone.nm,]$size), bty = "n", cex = 1)
  }
}
rasterpdf::dev.off()


### show clone distribution across cell states
inds <- rowSums(clone.states) >= 100
sum(inds)
mat <- t( abs(clone.states2[inds,]-1)^(0.5) * sign(clone.states2[inds,]-1) )
hcl <- hclust(dist(mat))
pheat <- pheatmap::pheatmap( mat,
   clustering_distance_rows = dist(mat),
  fontsize = 15, #fontsize_row = 20, fontsize_col = 20,
  angle_col = 0, cellwidth = 2.5,
  cluster_cols = T, cluster_rows = T,
  #labels_row = unlist(lapply(strsplit(coarse_ann_new, ' '), function(x) x[1])) ,
  show_colnames = F,
  annotation_col =  data.frame(clone.info[inds, 'genotype', drop = F]),
  annotation_row =  data.frame(clust = factor(1:10)),
  annotation_names_row = F, annotation_names_col = F,
  show_rownames = F, annotation_legend = T,
  annotation_colors = list(genotype = col.genotype,
                           clust = clust.final.cols)
)
#
pdf("figures/rna_clone_cluster_heatmap.pdf", width = 8, height = 5)
pheat
dev.off()

### show clone coupling across cell states
x <- cor(clone.states2[inds,], method = 'spearman'); diag(x) <- NA
pheat <- pheatmap::pheatmap(x, main = 'clone coupling of clusters',
                            clustering_distance_rows = dist(mat),
                            clustering_distance_cols = dist(mat),
                            annotation_row =  data.frame(clust = factor(1:10)),
                            annotation_col =  data.frame(clust = factor(1:10)),
                            annotation_colors = list(clust = clust.final.cols),
                            annotation_names_row = F, annotation_names_col = F, show_colnames = F, show_rownames = F, annotation_legend = F)
#
pdf("figures/rna_clone_coupling_heatmap.pdf", width = 4, height = 4)
pheat
dev.off()


### complexity of clone biases across cell states
clone.states.melt <- data.table::melt(as.data.frame(cbind(clone.states2[inds,],
                                              geno = clone.info[rownames(clone.states2[inds,]),]$genotype)))
#
dt <- data.table::as.data.table(clone.states2[inds, ])
dt[, geno := clone.info[rownames(clone.states2[inds, ]), "genotype"]]
clone.states.melt <- data.table::melt(dt, id.vars = "geno")
pl <- #ggplot(clone.states.melt, aes(x = geno, y = log(value), col = geno)) +
  ggplot(clone.states.melt, aes(x = value, col = geno)) +
  geom_density(bw = 0.2, alpha = I(1), lwd = 1.5) + scale_x_log10() +
  theme_bw() +
  #geom_boxplot(lwd = 0.1) + geom_point(size = I(0.1)) + scale_y_log10() +
  geom_vline(xintercept = 1, lty = 2, lwd = 0.2) +
  facet_wrap(~variable, scales = 'free_x', nrow = 2)+#,
  #labeller = labeller(variable = c("PC1" = "clonal PC1","PC2" = "clonal PC2"))) +
  scale_color_manual(values = col.genotype) + #scale_fill_manual(values = col.genotype)
  labs(x = NULL, y = NULL) + theme(axis.text.y = element_blank(), text = element_text(size = 20),
                                   axis.text = element_text(size = 10),
                                   legend.title = element_blank(), legend.position = "bottom") + #, legend.position = 'none')
  labs(x = 'clone enrichment in a cluster', y = 'density')
pl

#
pdf("figures/rna_clone_clusters.pdf", width = 8, height = 4)
pl
dev.off()
# assess stat significance
do.call(rbind, lapply(1:ncol(clone.states2), function(cl){
  ix0 <- clone.states.melt$variable == cl & clone.states.melt$geno == 'WT'
  ix1 <- clone.states.melt$variable == cl & clone.states.melt$geno == 'TET2'
  ix2 <- clone.states.melt$variable == cl & clone.states.melt$geno == 'DNMT3A'
  #tst1 <- wilcox.test(clone.states.melt[ix0,]$value, clone.states.melt[ix1,]$value)
  #tst2 <- wilcox.test(clone.states.melt[ix0,]$value, clone.states.melt[ix2,]$value)
  tst1 <- ks.test(log(clone.states.melt[ix0,]$value+1), log(clone.states.melt[ix1,]$value+1))
  tst2 <- ks.test(log(clone.states.melt[ix0,]$value+1), log(clone.states.melt[ix2,]$value+1))
  c(tst1$p.value, tst2$p.value)
}))
# 3(2): 0.013
# 9: 0.025
# 10:  0.015

### analyse clonal PCs
#inds <- rowSums(clone.states) > 100
sum(inds)
pr.clone <- prcomp(clone.states2[inds,], scale = T)
pr.clone.loadings <- predict(pr.clone, newdata = clone.states2)
clonevar <- pr.clone$sdev^2 / sum(pr.clone$sdev^2)
#
pl1 <- ggplot(data.frame(pr.clone$rotation[,]), aes(PC1, PC2,
                                                    label = rownames(pr.clone$rotation), col = I(clust.final.cols))) +
  theme_bw() + geom_point(size = I(5), alpha = I(1)) +
  geom_text_repel(size = I(10)) +
  labs(x = paste('clonal PC1, ', round(clonevar[1]*1e+2, 0), '%', sep = ''), y = paste('clonal PC2, ', round(clonevar[2]*1e+2, 0), '%', sep = '')) +
  theme(text = element_text(size = 20), legend.position = 'none',
        axis.title.x = element_text(margin = margin(t = 20, b = 0, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 20, l = 0, unit = "pt")))
pl1
#
pdf("figures/rna_clone_PCA.pdf",
    width = 5, height = 5)
pl1
dev.off()

geno <- clone.info[rownames(pr.clone$x),]$genotype
xx0 <- pr.clone$x[,c(1,2)]; #colnames(xx0) <- c('clonal PC1', 'clonal PC2')
#xx <- data.table::melt(data.frame(xx0, geno = clone.info[rownames(xx0),]$genotype))
dt <- data.table::as.data.table(xx0)
dt[, geno := clone.info[rownames(xx0), "genotype"]]
xx <- data.table::melt(dt, id.vars = "geno")
pl2 <- ggplot(xx,aes(x = value, fill = geno, col = geno)) +
  geom_density(bw = 1, alpha = I(0.2), lwd = 1) + theme_bw() +
  facet_wrap(~variable, scales = 'free', nrow = 2, labeller = labeller(variable = c("PC1" = "clonal PC1","PC2" = "clonal PC2"))) +
  scale_color_manual(values = col.genotype) + scale_fill_manual(values = col.genotype) +
  labs(x = NULL, y = NULL) + theme(axis.text.y = element_blank(), text = element_text(size = 20),
                                   legend.title = element_blank(), legend.position = "bottom")#, legend.position = 'none')
pl2
#
pdf("figures/rna_clonePCA_density_new.pdf", width = 3.5*(1.5+1), height = 3.5*1.5)
#grid.arrange(pl1, pl2, nrow = 1)
grid.arrange(arrangeGrob(pl1), arrangeGrob(pl2), nrow = 1, widths = c(1.5, 1))
dev.off()

### test  differences
sig <- pr.clone$x[,]
mut1 <- clone.info[rownames(sig),]$genotype == 'DNMT3A'
mut2 <- clone.info[rownames(sig),]$genotype == 'TET2'
wt <- clone.info[rownames(sig),]$genotype == 'WT'
i <- 1
tst1 <- t.test(sig[wt,i], sig[mut2,i])
tst2 <- wilcox.test(sig[wt,i], sig[mut2,i])
c(tst1$p.value, tst2$p.value)

# plot selected clones
# 213: lymphoid 217?
# 16: myeloid
# 396: hsc biased  555
titl <- c('myeloid biased', 'lymphoid biased', 'HSC biased')
indxs <- c(16, 213, 396)
#indxs <- match(rownames(clone.states2[inds,])[order(-clone.states2[inds, 3])][1:3], rownames(clone.states2))
#cols <- col.genotype[clone.info[indxs,]$genotype]
cols <- c('#0000FF', '#8000FF', '#FFA500')
#pdf("figures/clones_examples.pdf", width = 12, height = 4)
png("figures/clones_examples.png", width = 900, height = 300)
par(mfrow = c(1,3), mar = c(1,1,4,1))
for (i in 1:length(indxs)){
  cln <- indxs[i]
  cls <- rownames(cell.clone)[cell.clone$cloneid == rownames(clone.info)[cln]]
  cls <- intersect(cls, rownames(emb))
  plot(emb, pch = 19, cex = 0.1, col = 'grey', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
       main = titl[i], font.main = 1, cex.main = 3)
  points(emb[cls,], xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n',
      #   pch = 21, cex = 1, lwd = 1.5, col = cols[i], bg = 'white')
        pch = 21, cex = 1, lwd = 1.5, col = cols[i], bg = 'white')

  legend('topright', paste(#clone.info$genotype[cln],
                           'clone', cln,
                           #clone.info$clone[cln],
                           #'\n',
                           '\nsize =', clone.info$size[cln]
                           #'id:', clone.info$id[cln],
  ), bty = 'n', cex = 3)
}
dev.off()

### Show distribution of individual clone phenotypes
inds <- rowSums(clone.states) >= 100
sum(inds)
mat <- t( abs(clone.states2[inds,]-1)^(0.5) * sign(clone.states2[inds,]-1) )

clones.show <- rownames(clone.states)[rowSums(clone.states) >= 20]
str(clones.show)
clone.states.show <- data.frame('My_bias' = clone.states2[clones.show, 8],
                                'HSC_bias' = clone.states2[clones.show, 8],
                                'Ly_bias' = clone.states2[clones.show, 4],
                                'size' = clone.info[clones.show,]$size,
                                genotype = clone.info[clones.show,]$genotype)
head(clone.states.show)
library(tidyr)
df_long <- clone.states.show %>%
  tidyr::pivot_longer(
    cols = c(size, 'My_bias', 'HSC_bias', 'Ly_bias'),
    names_to = "variable",
    values_to = "value"
  )
df_long$variable <- factor(df_long$variable,
                           levels = c("size", "My_bias", "Ly_bias", "HSC_bias"))
df_long$variable <- dplyr::recode(df_long$variable,
                                  My_bias = "My bias",
                                  HSC_bias = "HSC bias",
                                  Ly_bias = "Ly bias",
                                  size = "Size")
text.scale <- 20
pl <- ggplot(df_long, aes(x = '', y = value)) + geom_boxplot(col = 'black') +
  geom_jitter(aes(bg = genotype), size = I(2), pch = 21, col = 'black') + facet_wrap(~ variable, scales = "free_y", , nrow = 1) +
  theme_bw() +
  labs(y = 'clonal phenotype', x = '') + scale_fill_manual(values = col.genotype) +
  theme(axis.title = element_text(size = text.scale), axis.text = element_text(size = text.scale),
        legend.title = element_blank(),
        legend.text = element_text(size = text.scale),
        strip.text = element_text(size = text.scale),
        #legend.position = c(0.15, 0.9)
        )
pl
pdf('figures/clone_phenotype_scatterplot.pdf', width = 9, height = 4)
pl
dev.off()
