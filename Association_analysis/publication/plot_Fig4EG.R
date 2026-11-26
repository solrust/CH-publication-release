# Reproduce review versions of Fig. 4E and 4G from small exported tables.
# Run after export_Fig4EG.R. Requires ggplot2 only; no large R checkpoint.
stopifnot(requireNamespace('ggplot2',quietly=TRUE))
pub <- 'Association_analysis/publication'
sd <- file.path(pub,'source_data')
fd <- file.path(pub,'figures')
dir.create(fd,showWarnings=FALSE,recursive=TRUE)

for (spec in list(
  list(panel='Fig4E',gene='LSP1',xlab='clone HSC bias',title='LSP1'),
  list(panel='Fig4G',gene='NPM1',xlab='clone growth rate',title='NPM1')
)) {
  p <- file.path(sd,paste0(spec$panel,'_',spec$gene,'_clone_expression.tsv'))
  if (!file.exists(p)) stop('Missing Source Data: ',p)
  x <- read.delim(p,stringsAsFactors=FALSE)
  x$genotype <- factor(x$genotype,levels=c('WT','DNMT3A','TET2'))
  x$experiment <- factor(x$experiment)
  pal <- c(WT='#333333', DNMT3A='#d45b50', TET2='#68a957')
  shapes <- c('exp123'=16,'exp127'=17)
  if (!all(levels(x$experiment) %in% names(shapes))) {
    shapes <- setNames(c(16,17,15,18)[seq_along(levels(x$experiment))],levels(x$experiment))
  }
  g <- ggplot2::ggplot(x,ggplot2::aes(x=phenotype_value,y=relative_gene_expression))+
    ggplot2::geom_point(ggplot2::aes(color=genotype,shape=experiment),alpha=0.8,size=2)+
    ggplot2::geom_smooth(ggplot2::aes(color=genotype),method='lm',formula=y~x,
                         se=TRUE,linewidth=0.6,alpha=0.12)+
    ggplot2::scale_color_manual(values=pal,drop=FALSE)+
    ggplot2::scale_shape_manual(values=shapes,drop=FALSE)+
    ggplot2::labs(x=spec$xlab,y='Expression in HSC clones (fraction)',title=spec$title,
                  color='Genotype',shape='Experiment')+
    ggplot2::theme_bw(base_size=11)+
    ggplot2::theme(plot.title=ggplot2::element_text(face='italic',hjust=0.5))
  outfile <- file.path(fd,paste0(spec$panel,'_reproduced_review.pdf'))
  ggplot2::ggsave(outfile,g,width=5.5,height=4.6)
  message('Created ',outfile)
}
