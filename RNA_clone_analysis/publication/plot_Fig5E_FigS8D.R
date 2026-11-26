# Replot Figure 5E and Supplementary Figure 8D from checkpoint-based exported tables.
local({
  if (!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2')
  library(ggplot2)
  pub <- 'RNA_clone_analysis/publication'
  src <- file.path(pub,'source_data'); out <- file.path(pub,'figures')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  paths <- file.path(src,c('Fig5E_ROC_coordinates.tsv.gz','Fig5E_FigS8D_clone_scores.tsv.gz'))
  if(any(!file.exists(paths))) stop('Run export_Fig5E_FigS8D.R first.')
  curves <- read.delim(paths[1],stringsAsFactors=FALSE)
  clones <- read.delim(paths[2],stringsAsFactors=FALSE)
  pal <- c(DNMT3A='#e46871',TET2='#7fc357')
  curves$curve <- paste(curves$contrast,curves$split)
  pl1 <- ggplot(curves,aes(x=specificity,y=sensitivity,color=contrast,linetype=split,
                             group=interaction(contrast,split))) +
    geom_path(linewidth=1) + geom_abline(slope=-1,intercept=1,color='grey70') +
    facet_wrap(~level,nrow=1,labeller=as_labeller(c(cell='cell level',clone='clone level'))) +
    scale_x_reverse(limits=c(1,0)) + coord_fixed(ratio=1) +
    scale_color_manual(values=pal) + scale_linetype_manual(values=c(train='dashed',test='solid')) +
    labs(x='Specificity',y='Sensitivity',color=NULL,linetype=NULL) + theme_bw(base_size=12)
  pdf(file.path(out,'Fig5E_reproduced_review.pdf'),width=7.5,height=3.8)
  print(pl1);dev.off()
  data <- clones[clones$clone_size>50 & clones$clone_genotype %in% c('WT','DNMT3A','TET2') &
                   is.finite(clones$score_train) & is.finite(clones$score_test),,drop=FALSE]
  data <- data[data$clone_genotype=='WT' | data$clone_genotype==data$contrast,,drop=FALSE]
  pal1 <- c(WT='#858585',DNMT3A='#ea6971',TET2='#80c35b')
  pl2 <- ggplot(data,aes(score_train,score_test,color=clone_genotype,size=clone_size)) +
    geom_point(alpha=0.8) + geom_smooth(method='lm',formula=y~x,se=TRUE,linewidth=0.6) +
    facet_wrap(~contrast,nrow=1,scales='free') +
    scale_color_manual(values=pal1) + scale_size_continuous(range=c(1,5),guide='none') +
    labs(x='train',y='test',color=NULL) + theme_bw(base_size=12)
  pdf(file.path(out,'FigS8D_reproduced_review.pdf'),width=8,height=3.7)
  print(pl2);dev.off()
  message('Created Fig5E and FigS8D review PDFs. Do not mark validated until compared to accepted artwork.')
})
