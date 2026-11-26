# Review plots from exported Fig 4F/4H/S7F numerical Source Data only.
# Run after export_Fig4FH_FigS7F.R; requires ggplot2, gridExtra;
# ggrepel is optional. Does not load the association checkpoint.

if (!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2')
if (!requireNamespace('gridExtra',quietly=TRUE)) stop('Install gridExtra')
pub <- 'Association_analysis/publication'
outdir <- file.path(pub,'figures')
sourcedir <- file.path(pub,'source_data')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

read_one <- function(file) {
  if (!file.exists(file)) stop('Missing export: ',file)
  read.delim(file,check.names=FALSE,stringsAsFactors=FALSE)
}

make_one <- function(panel,phenotype) {
  x <- read_one(file.path(sourcedir,paste0(panel,'_gene_associations.tsv.gz')))
  y <- read_one(file.path(sourcedir,paste0(panel,'_selected_pathways.tsv')))
  x$plot_p <- pmax(x$pvalue,.Machine$double.xmin)
  x$point_alpha <- ifelse(x$label_candidate,1,0.11)
  # Retain unadjusted P on a logarithmic axis as in historical artwork.
  pv <- ggplot2::ggplot(x,ggplot2::aes(x=coefficient_log2,y=plot_p)) +
    ggplot2::geom_point(ggplot2::aes(alpha=point_alpha),color='gray35',size=.65,na.rm=TRUE) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_y_log10() +
    ggplot2::geom_vline(xintercept=0,linetype='dashed',linewidth=.3) +
    ggplot2::labs(title=paste0(phenotype,' gene associations'),
                  x='Association coefficient, log2 scale',y='Gene association P value (unadjusted)') +
    ggplot2::theme_bw(base_size=11)
  labs <- x[x$label_candidate & is.finite(x$coefficient_log2) & is.finite(x$plot_p),,drop=FALSE]
  if (nrow(labs)>60) labs <- labs[order(labs$pvalue),,drop=FALSE][seq_len(60),,drop=FALSE]
  if (nrow(labs) > 0 && requireNamespace('ggrepel',quietly=TRUE)) {
    pv <- pv + ggrepel::geom_text_repel(data=labs,ggplot2::aes(label=gene),
                                       size=2.2,max.overlaps=15,box.padding=.3)
  }

  if (nrow(y)>0) {
    y <- y[is.finite(y$signed_minus_log10_p),,drop=FALSE]
  }
  if (nrow(y)>0) {
    y$direction <- ifelse(y$NES>=0,'positive','negative')
    y <- y[order(y$signed_minus_log10_p),,drop=FALSE]
    y$pathway_label <- factor(y$pathway,levels=y$pathway)
    pb <- ggplot2::ggplot(y,ggplot2::aes(x=pathway_label,y=signed_minus_log10_p,fill=direction)) +
      ggplot2::geom_col(color='gray30',linewidth=.25) + ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values=c(negative='#a9b9d9',positive='#e6aea7')) +
      ggplot2::theme_bw(base_size=9) +
      ggplot2::theme(legend.position='none',axis.title.y=ggplot2::element_blank()) +
      ggplot2::labs(y='Signed -log10(fgsea P value)')
  } else {
    pb <- ggplot2::ggplot() + ggplot2::annotate('text',x=0,y=0,label='No finite pathways selected') +
      ggplot2::theme_void()
  }
  gridExtra::arrangeGrob(pv,pb,nrow=1,widths=c(1.1,1.3))
}

p1 <- make_one('Fig4F','HSC bias')
p2 <- make_one('Fig4H','growth rate')
pdf(file.path(outdir,'Fig4F_reproduced_review.pdf'),width=15,height=7.5,onefile=TRUE);grid::grid.draw(p1);dev.off()
pdf(file.path(outdir,'Fig4H_reproduced_review.pdf'),width=15,height=7.5,onefile=TRUE);grid::grid.draw(p2);dev.off()

s7 <- c('clone_size'='clone size', 'erythroid_bias'='erythroid bias',
        'myeloid_bias'='myeloid bias','lymphoid_bias'='lymphoid bias',
        'dendritic_bias'='dendritic bias')
pdf(file.path(outdir,'FigS7F_reproduced_review.pdf'),width=15,height=7.5,onefile=TRUE)
for (x in names(s7)) {
  grid::grid.newpage()
  grid::grid.draw(make_one(paste0('FigS7F_',x),s7[[x]]))
}
dev.off()
message('Created Fig4F, Fig4H and a 5-page FigS7F review PDF (one phenotype per page).')
