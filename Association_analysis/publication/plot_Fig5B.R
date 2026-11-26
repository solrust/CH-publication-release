# Figure 5B: plot the 21 already-exported cells (no checkpoint loading).
# Use historical first-observed-group sign for the initial visual audit.
# If published direction differs, compare to the explicit up-minus-down column.
if (!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2')
pub <- 'Association_analysis/publication'
input <- file.path(pub,'source_data','Fig5B_heatmap_statistics.tsv')
if (!file.exists(input)) stop('Run export_Fig5B.R first')
dat <- read.delim(input,check.names=FALSE,stringsAsFactors=FALSE)
wanted <- c('growth rate','clone size','HSC bias','lymphoid bias',
            'dendritic bias','myeloid bias','erythroid bias')
if (!identical(dat$phenotype, rep(wanted,each=3L))) stop('Unexpected heatmap row order')
if (!all(c('DNMT3A','TET2','AML') %in% dat$signature)) stop('Missing signatures')
# Discrete placement with a slight visual break between CH and AML.
dat$y <- match(dat$phenotype,rev(wanted))
dat$x <- c(DNMT3A=1,TET2=2,AML=3.25)[dat$signature]
dat$score <- dat$signed_neglog10p_historical
if (any(!is.finite(dat$score))) stop('Non-finite signed values')
lim <- max(abs(dat$score),na.rm=TRUE)
p <- ggplot2::ggplot(dat,ggplot2::aes(x=x,y=y,fill=score)) +
  ggplot2::geom_tile(width=0.96,height=0.96,color='black',linewidth=0.55) +
  ggplot2::geom_text(ggplot2::aes(label=significance_stars),size=8,fontface='bold') +
  ggplot2::scale_fill_gradient2(low='#5150a6',mid='#f3f2f2',high='#961a24',
     limits=c(-lim,lim),midpoint=0,name='Signed -log10(P)') +
  ggplot2::scale_x_continuous(breaks=c(1,2,3.25),labels=c('DNMT3A','TET2','AML'),
     limits=c(0.45,3.85),expand=c(0,0),position='top') +
  ggplot2::scale_y_continuous(breaks=seq_along(wanted),labels=rev(wanted),
     limits=c(0.4,7.6),expand=c(0,0)) +
  ggplot2::coord_fixed(ratio=0.6,clip='off') +
  ggplot2::labs(title='Fig. 5B — checkpoint reproduction (REVIEW)',x=NULL,y=NULL) +
  ggplot2::theme_minimal(base_size=13) +
  ggplot2::theme(panel.grid=ggplot2::element_blank(),
    axis.text.x=ggplot2::element_text(face='italic',size=14),
    axis.text.y=ggplot2::element_text(size=12),
    axis.ticks=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(size=14,face='bold'),
    legend.position='right')
dir.create(file.path(pub,'figures'),recursive=TRUE,showWarnings=FALSE)
output <- file.path(pub,'figures','Fig5B_reproduced_REVIEW.pdf')
ggplot2::ggsave(output,p,width=9.0,height=6.5,units='in',device='pdf')
message('Created ',output)
