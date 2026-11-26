# Reproduce review plots from exported portable TSVs; no original RData required.
# Run from CH_publication_release root.
if (!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2 to make S6 plots')
sd_dir <- 'RNA_analysis/publication/source_data'
fig_dir <- 'RNA_analysis/publication/figures'
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
a <- read.delim(gzfile(file.path(sd_dir,'Fig3A_umap_cells.tsv.gz')),
                stringsAsFactors=FALSE,na.strings=c('','NA'))
a <- a[is.finite(a$umap_1) & is.finite(a$umap_2),]
col <- a$state_color
col[is.na(col) | col==''] <- 'grey90'
centroids <- aggregate(cbind(umap_1,umap_2) ~ cell_state,
                       data=a[!is.na(a$cell_state),],FUN=mean)
# Fig 3A main UMAP
pdf(file.path(fig_dir,'Fig3A_umap_reproduced.pdf'),width=5,height=5,useDingbats=FALSE)
par(mar=rep(0.2,4))
plot(a$umap_1,a$umap_2,pch=19,cex=.12,col=adjustcolor('grey90',alpha.f=.2),
     axes=FALSE,ann=FALSE,bty='n',asp=1)
ind <- !is.na(a$cell_state)
points(a$umap_1[ind],a$umap_2[ind],pch=19,cex=.18,
       col=adjustcolor(col[ind],alpha.f=.25))
points(centroids$umap_1,centroids$umap_2,pch=21,cex=2.5,
       bg='grey95',col='black')
text(centroids$umap_1,centroids$umap_2,labels=centroids$cell_state,
     col='black',cex=.9,font=2)
dev.off()

# Fig3A genotype inset; original randomly shuffled cell rendering order.
# Use reproducible shuffle solely for draw order (not the Source Data).
set.seed(101)
i <- sample(seq_len(nrow(a)))
acol <- c('WT'='black','DNMT3A'='red','TET2'='green')[a$genotype]
acol[is.na(acol)] <- 'grey85'
pdf(file.path(fig_dir,'Fig3A_genotype_inset_reproduced.pdf'),width=4,height=4,
    useDingbats=FALSE)
par(mar=rep(.2,4))
plot(a$umap_1[i],a$umap_2[i],pch=19,cex=.1,axes=FALSE,ann=FALSE,bty='n',
     col=adjustcolor(acol[i],alpha.f=.8),asp=1)
legend('topright',c('WT','DNMT3A','TET2'),pch=19,
       col=c('black','red','green'),bty='n',cex=.85)
dev.off()

# Supplementary Fig 6A: sample-level boxplot + individual samples.
# Deliberately no manually inserted significance asterisks.
sa <- read.delim(file.path(sd_dir,'FigS6A_sample_state_enrichment.tsv'),
                 stringsAsFactors=FALSE)
sa$cell_state <- factor(sa$cell_state,levels=1:10)
sa$experiment <- factor(sa$experiment,levels=c('123','127'))
sa$genotype <- factor(sa$genotype,levels=c('WT','DNMT3A','TET2'))
p6a <- ggplot2::ggplot(sa,ggplot2::aes(x=cell_state,y=normalized_enrichment,
                                      color=genotype)) +
  ggplot2::geom_boxplot(width=.8,outlier.shape=NA,alpha=1) +
  ggplot2::geom_point(position=ggplot2::position_jitterdodge(
    jitter.width=.2,dodge.width=.8),size=1.9,alpha=.65) +
  ggplot2::geom_hline(yintercept=1,linetype=2) +
  ggplot2::facet_wrap(~experiment,nrow=2) +
  ggplot2::scale_y_log10(limits=c(.5,2)) +
  ggplot2::scale_color_manual(values=c('WT'='black','DNMT3A'='red','TET2'='green'),
                              drop=TRUE) +
  ggplot2::labs(x='clusters',y='cell state enrichment') +
  ggplot2::theme_bw(base_size=11)
ggplot2::ggsave(file.path(fig_dir,'FigS6A_reproduced_without_manual_asterisks.pdf'),
                plot=p6a,width=7,height=5)

# Supplementary Fig 6B: Pearson correlation, original code used cor() and
# cor.test() without method argument, hence Pearson r, NOT Spearman rho.
f <- read.delim(file.path(sd_dir,'FigS6B_state_fold_changes.tsv'),
                stringsAsFactors=FALSE)
cstat <- read.delim(file.path(sd_dir,'FigS6B_correlation.tsv'),
                   stringsAsFactors=FALSE)
x <- f[f$experiment=='123',c('cell_state','fold_change')]
y <- f[f$experiment=='127',c('cell_state','fold_change')]
names(x)[2] <- 'DNMT3A'
names(y)[2] <- 'TET2'
xy <- merge(x,y,by='cell_state',sort=TRUE)
xy$cell_state <- factor(xy$cell_state,levels=1:10)
fit_label <- sprintf('Pearson r = %.2f\np = %.3g',cstat$pearson_r,cstat$pvalue)
p6b <- ggplot2::ggplot(xy,ggplot2::aes(x=DNMT3A,y=TET2)) +
  ggplot2::geom_point(ggplot2::aes(color=cell_state),size=3.8) +
  ggplot2::geom_text(ggplot2::aes(label=cell_state),vjust=-.9,size=3.8) +
  ggplot2::scale_color_manual(values=setNames(
     vapply(1:10, function(k) {
       q <- a$state_color[!is.na(a$cell_state) & a$cell_state==k]
       q <- q[!is.na(q) & nzchar(q)]
       if (length(q)) q[1] else 'grey50'
     },character(1)), as.character(1:10))) +
  ggplot2::annotate('text',x=min(xy$DNMT3A),y=min(xy$TET2),
                   label=fit_label,hjust=0,vjust=0,size=3.5) +
  ggplot2::labs(x='DNMT3A',y='TET2',title='cell state enrichment') +
  ggplot2::theme_bw(base_size=12) +
  ggplot2::theme(legend.position='none')
ggplot2::ggsave(file.path(fig_dir,'FigS6B_reproduced.pdf'),
                plot=p6b,width=4.5,height=4.5)

# Fig3B plots are optional until nearest-neighbour export is completed.
bfile <- file.path(sd_dir,'Fig3B_smoothed_enrichment.tsv.gz')
if (file.exists(bfile)) {
  b <- read.delim(gzfile(bfile),stringsAsFactors=FALSE)
  pal <- grDevices::colorRampPalette(c('blue','grey90','red'))(101)
  # Relative palette index follows original integer encoding.
  png(file.path(fig_dir,'Fig3B_reproduced.png'),width=400,height=800)
  par(mfrow=c(2,1),mar=rep(.5,4))
  for(geno in c('DNMT3A','TET2')) {
    q <- b[b$genotype==geno,]
    ci <- pmax(1L,pmin(101L,as.integer(q$smoothed_index)))
    plot(a$umap_1,a$umap_2,col='grey95',cex=.05,pch=19,axes=FALSE,ann=FALSE,bty='n')
    points(q$umap_1,q$umap_2,col=pal[ci],pch=19,cex=.5)
    legend('topright',legend=geno,bty='n',cex=2)
  }
  dev.off()
  cat('Created Fig3B_reproduced.png\n')
} else cat('Fig3B missing; run export_Fig3B.R then rerun this plot script.\n')
cat('Review PDFs/PNG created under ',fig_dir,'\n',sep='')
