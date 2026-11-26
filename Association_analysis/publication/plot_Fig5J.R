# Fig 5J review plot based on exported Spearman means and SDs.
# Run from the repository root after export_Fig5J.R.
if (!requireNamespace('ggplot2', quietly=TRUE)) stop('Install ggplot2')
nsplits <- as.integer(getOption('CH.fig5j.splits',300L))
tag <- if (nsplits == 300L) '' else paste0('_SMOKETEST_',nsplits,'splits')
dir <- 'Association_analysis/publication'
input <- file.path(dir,'source_data',paste0('Fig5J_heritability_summary',tag,'.tsv'))
if (!file.exists(input)) stop('Missing summary: ',input)
x <- read.delim(input, check.names=FALSE)
x <- x[x$phenotype != 'dendritic.bias',,drop=FALSE]
wt <- x[x$genotype=='WT', c('phenotype','mean_rho','sd_rho')]
colnames(wt) <- c('phenotype','WT_mean','WT_sd')
y <- merge(x[x$genotype %in% c('DNMT3A','TET2'),], wt, by='phenotype',sort=FALSE)
y$lineage <- factor(y$lineage,levels=c('B Lymph','GMP','Ery'))
y$genotype <- factor(y$genotype, levels=c('DNMT3A','TET2'))
colors <- c('DNMT3A'='#C94F4F','TET2'='#398A43')
p <- ggplot2::ggplot(y,ggplot2::aes(x=WT_mean,y=mean_rho,color=genotype)) +
  ggplot2::geom_abline(intercept=0,slope=1, linetype='dashed',color='grey60') +
  ggplot2::geom_segment(ggplot2::aes(x=WT_mean-WT_sd,xend=WT_mean+WT_sd,y=mean_rho,yend=mean_rho),linewidth=0.6) +
  ggplot2::geom_segment(ggplot2::aes(x=WT_mean,xend=WT_mean,y=mean_rho-sd_rho,yend=mean_rho+sd_rho),linewidth=0.6) +
  ggplot2::geom_point(size=3) +
  ggplot2::geom_text(ggplot2::aes(label=as.character(lineage)),hjust=-0.12,vjust=-0.7,show.legend=FALSE,size=3.5) +
  ggplot2::scale_color_manual(values=colors) +
  ggplot2::labs(x='WT HSC-to-phenotype Spearman correlation',
       y='Mutant HSC-to-phenotype Spearman correlation',color='Genotype',
       title='Fig. 5J checkpoint reproduction (review)') +
  ggplot2::coord_cartesian(xlim=c(0,0.5),ylim=c(0,0.5)) +
  ggplot2::theme_bw(base_size=11)
dir.create(file.path(dir,'figures'),recursive=TRUE,showWarnings=FALSE)
out <- file.path(dir,'figures',paste0('Fig5J_reproduced_REVIEW',tag,'.pdf'))
ggplot2::ggsave(out,p,width=7,height=5)
message('Created ',out)
