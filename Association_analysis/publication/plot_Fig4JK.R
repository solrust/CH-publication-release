# Review plots for Fig 4J/K from exported TSVs only.
pub <- 'Association_analysis/publication'
if(!requireNamespace('ggplot2',quietly=TRUE)) stop('ggplot2 required')
figdir <- file.path(pub,'figures');dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
j <- utils::read.delim(file.path(pub,'source_data','Fig4JK_gene_signature_association_scores.tsv.gz'))
tests <- utils::read.delim(file.path(pub,'validation','Fig4J_wilcoxon_summary.tsv'))
if(anyDuplicated(j[,c('gene','signature','phenotype')])) stop('Duplicate gene-signature-phenotype rows')
j <- j[j$phenotype %in% c('HSC bias','growth rate'),,drop=FALSE]
j$phenotype <- factor(j$phenotype,levels=c('HSC bias','growth rate'))
j$signature <- factor(j$signature,levels=c('HSC dormancy','HSC label-retention','HSC output'))
set.seed(723)
p <- ggplot2::ggplot(j,ggplot2::aes(x=group,y=association_score,color=group))+
  ggplot2::geom_hline(yintercept=0,linetype=2,color='grey60')+
  ggplot2::geom_boxplot(outlier.shape=NA,alpha=0.4,width=0.6)+
  ggplot2::geom_jitter(width=0.18,height=0,alpha=0.55,size=0.8)+
  ggplot2::facet_grid(phenotype~signature,scales='free_y')+
  ggplot2::labs(x=NULL,y='Gene-association score: signed -ln(P)',
                title='Fig 4J numerical review; gene-level Wilcoxon p-values in validation table')+
  ggplot2::theme_bw(base_size=11)+ggplot2::theme(legend.position='none',axis.text.x=ggplot2::element_text(angle=25,hjust=1))
ggplot2::ggsave(file.path(figdir,'Fig4J_REVIEW.pdf'),p,width=12,height=7,useDingbats=FALSE)

r <- utils::read.delim(file.path(pub,'source_data','Fig4K_ROC_coordinates.tsv.gz'))
a <- utils::read.delim(file.path(pub,'validation','Fig4K_ROC_AUC_summary.tsv'))
r$signature<-factor(r$signature,levels=c('HSC dormancy','HSC label-retention','HSC output'))
r$phenotype<-factor(r$phenotype,levels=c('growth rate','HSC bias','clone size'))
p2<-ggplot2::ggplot(r,ggplot2::aes(false_positive_rate,sensitivity,color=phenotype,group=phenotype))+
  ggplot2::geom_abline(slope=1,intercept=0,linetype=2,color='grey60')+
  ggplot2::geom_path(linewidth=0.85)+
  ggplot2::facet_wrap(~signature,nrow=1)+
  ggplot2::coord_equal(xlim=c(0,1),ylim=c(0,1))+
  ggplot2::labs(x='1 - specificity',y='Sensitivity',color=NULL,
                title='Fig 4K ROC review; AUCs in validation table')+
  ggplot2::theme_bw(base_size=12)
ggplot2::ggsave(file.path(figdir,'Fig4K_REVIEW.pdf'),p2,width=13,height=4.3,useDingbats=FALSE)
message('Created Fig4J_REVIEW.pdf and Fig4K_REVIEW.pdf; do not approve until accepted-artwork comparison.')
