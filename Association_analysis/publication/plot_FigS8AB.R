# Review plots for S8A/S8B. Source Data from export_FigS8AB.R.
pub <- 'Association_analysis/publication'
sd <- file.path(pub,'source_data'); fg <- file.path(pub,'figures')
dir.create(fg,recursive=TRUE,showWarnings=FALSE)
for(nm in c('ggplot2','pheatmap')) if(!requireNamespace(nm,quietly=TRUE)) stop('Install ',nm)
cor_tbl <- read.delim(file.path(sd,'FigS8A_correlations.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
needed <- c('growth rate','clone size','HSC bias','lymphoid bias','dendritic bias','myeloid bias','erythroid bias')
a <- matrix(NA_real_,length(needed),length(needed),dimnames=list(needed,needed))
for(i in seq_len(nrow(cor_tbl)))a[cor_tbl$y[i],cor_tbl$x[i]]<-cor_tbl$r_pearson[i]
if(anyNA(a))stop('Incomplete S8A correlation matrix')
# Original: complete linkage of distance 1 - Pearson correlation.
d <- 1-a; diag(d)<-0
cl <- stats::hclust(stats::as.dist(d),method='complete')
diag(a)<-NA_real_
colors <- grDevices::colorRampPalette(c('blue','white','red'))(100)
pdf(file.path(fg,'FigS8A_reproduced_REVIEW.pdf'),width=5.3,height=5)
pheatmap::pheatmap(a,cluster_rows=cl,cluster_cols=cl,color=colors,
   breaks=seq(-0.5,0.5,length.out=101),main='Similarity of clone phenotypes',
   fontsize=10,na_col='white')
dev.off()
x<-read.delim(file.path(sd,'FigS8B_clone_phenotype_pairs.tsv.gz'),check.names=FALSE)
xy <- unique(x[,c('x_phenotype','y_phenotype','pearson_r')]); if(nrow(xy)!=10)stop('Need 10 scatter panels')
order_y <- c('growth rate','HSC bias')
order_x <- c('clone size','HSC bias','erythroid bias','myeloid bias','lymphoid bias','growth rate')
x$panel <- paste(x$y_phenotype,' vs ',x$x_phenotype,sep='')
x$panel<-factor(x$panel,levels=unique(x$panel))
anno<-data.frame(panel=factor(paste(xy$y_phenotype,' vs ',xy$x_phenotype,sep=''),levels=levels(x$panel)),
                 label=sprintf('r = %.2f',xy$pearson_r),stringsAsFactors=FALSE)
p<-ggplot2::ggplot(x,ggplot2::aes(x_value,y_value)) +
  ggplot2::geom_point(size=1.3,color='grey30',alpha=0.5)+
  ggplot2::geom_smooth(method='lm',formula=y~x,color='#2d64ce',fill='grey80',alpha=0.4,se=TRUE)+
  ggplot2::geom_text(data=anno,ggplot2::aes(label=label),x=Inf,y=Inf,
     hjust=1.1,vjust=1.5,color='darkred',size=3.1,inherit.aes=FALSE)+
  ggplot2::facet_wrap(~panel,ncol=5,scales='free')+
  ggplot2::theme_bw(base_size=10)+
  ggplot2::labs(x='Clone phenotype (x)',y='Clone phenotype (y)',title='S8B clone phenotype Pearson correlations (review)')+
  ggplot2::theme(strip.text=ggplot2::element_text(size=8))
ggplot2::ggsave(file.path(fg,'FigS8B_reproduced_REVIEW.pdf'),p,width=15,height=6,device='pdf')
cat('Created S8A/B review PDFs.\n')
