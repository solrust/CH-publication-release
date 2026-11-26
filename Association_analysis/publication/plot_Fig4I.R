# Fig 4I review plot from exported numerical data + optional provisional
# GO-category annotation. No checkpoint required at plotting time.

pub <- 'Association_analysis/publication'
sd_path <- file.path(pub,'source_data','Fig4I_pathway_association_scores.tsv.gz')
cat_path <- file.path(pub,'validation','Fig4I_GO_categories_REVIEW.tsv')
plot_dir <- file.path(pub,'figures')
dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(sd_path)) stop('Run export_Fig4I.R first')
if(!requireNamespace('ggplot2',quietly=TRUE)) stop('ggplot2 required')

x <- read.delim(gzfile(sd_path),check.names=FALSE,stringsAsFactors=FALSE)
x$GO_category <- NA_character_
if(file.exists(cat_path)) {
  cat <- read.delim(cat_path,stringsAsFactors=FALSE)
  if(anyDuplicated(cat$pathway)) stop('Duplicate pathway classifications')
  x$GO_category <- cat$GO_category_review[match(x$pathway,cat$pathway)]
  message('Joined review-only pathway categories. Compare with accepted artwork.')
} else {
  message('No GO category audit yet. Drawing all pathways in grey.')
}
# Plot grey (unclassified) pathways first, highlighted pathways on top.
x <- x[order(!is.na(x$GO_category)),,drop=FALSE]
valid <- is.finite(x$growth_rate_score) & is.finite(x$HSC_bias_score)
message('Fig4I: ',sum(valid),' finite pathways, ',sum(!is.na(x$GO_category)&valid),' colored.')

# The historical figure uses ggplot2 discrete colors. These are the conventional
# 4-group ggplot hues in the alphabetical order of the annotation names.
palette <- c('cell cycle'='#F8766D','immune activation'='#7CAE00',
             'metabolism'='#00BFC4','translation'='#C77CFF')
p <- ggplot2::ggplot(x[valid,,drop=FALSE],
  ggplot2::aes(x=growth_rate_score,y=HSC_bias_score,color=GO_category)) +
  ggplot2::geom_point(size=2.5,alpha=0.75) +
  ggplot2::geom_vline(xintercept=0,linetype='dashed',color='black') +
  ggplot2::geom_hline(yintercept=0,linetype='dashed',color='black') +
  ggplot2::coord_cartesian(xlim=c(-30,30),ylim=c(-13,13)) +
  ggplot2::scale_color_manual(values=palette,na.value='grey60',
                               breaks=names(palette),drop=FALSE) +
  ggplot2::labs(x='Clone growth rate, log pathway significance',
                y='HSC-biased clones, log pathway significance',color=NULL) +
  ggplot2::theme_bw(base_size=12) +
  ggplot2::theme(legend.position=c(0.98,0.98),
       legend.justification=c('right','top'),
       legend.background=ggplot2::element_rect(fill='white',color='grey80'))
out <- file.path(plot_dir,'Fig4I_reproduced_REVIEW.pdf')
ggplot2::ggsave(out,p,width=6.3,height=5.5,units='in',useDingbats=FALSE)
message('Created ',out)
