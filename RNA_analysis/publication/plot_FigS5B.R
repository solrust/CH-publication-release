# Plot Supplementary Fig. 5B from exported QC Source Data only.
sdir <- 'RNA_analysis/publication/source_data'
odir <- 'RNA_analysis/publication/figures'
dir.create(odir,recursive=TRUE,showWarnings=FALSE)
cells <- read.delim(gzfile(file.path(sdir,'FigS5B_cell_QC.tsv.gz')))
clusters <- read.delim(file.path(sdir,'FigS5B_cluster_QC.tsv'))
stopifnot(nrow(cells)>0,nrow(clusters)>0)
pdf(file.path(odir,'FigS5B_reproduced.pdf'),width=9,height=9)
par(mfrow=c(2,2),mar=c(4.5,4.5,3,1))
# Mito-high UMAP
pal <- colorRampPalette(c('grey90','grey60','red'))(10)
cols <- pal[as.integer(cut(cells$mito_PC1,breaks=10))]
plot(cells$UMAP_1,cells$UMAP_2,col=cols,pch=19,cex=.07,
     main='mito-score',xlab='',ylab='',xaxt='n',yaxt='n')
legend('bottomleft',c('high','low'),col=c('red','grey90'),pch=19,bty='n')
# Score histogram and historical threshold
hist(cells$mito_PC1,breaks=100,main='',xlab='mito-score',ylab='count',col='grey70')
abline(v=11,col=adjustcolor('red',.5),lwd=3)
# Hashtag doublet UMAP
plot(cells$UMAP_1,cells$UMAP_2,
     col=ifelse(cells$hashtag_doublet,'blue','grey65'),
     pch=19,cex=.07,main='hashtag doublets',xlab='',ylab='',xaxt='n',yaxt='n')
# Proportion of hashtag-defined doublets per cluster
plot(clusters$n_cells,clusters$doublet_fraction, log='x', ylim=c(0,.5),
     pch=21,bg=adjustcolor(ifelse(clusters$doublet_cluster_flag,'blue','grey60'),.5),
     xlab='cluster size',ylab='fraction doublets')
abline(h=.3,col=adjustcolor('blue',.5),lwd=3)
dev.off()
message('Created ',file.path(odir,'FigS5B_reproduced.pdf'))
