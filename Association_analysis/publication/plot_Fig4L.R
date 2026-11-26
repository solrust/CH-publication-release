# Fig 4L numerical heatmap review only. Input categorical annotations are not final.
pub<-'Association_analysis/publication';figure_dir<-file.path(pub,'figures')
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
if(!requireNamespace('ggplot2',quietly=TRUE)) stop('ggplot2 required')
x<-utils::read.delim(file.path(pub,'source_data','Fig4L_pathway_score_matrix_REVIEW.tsv.gz'),check.names=FALSE)
ph<-c('lymphoid bias','clone size','growth rate','erythroid bias','myeloid bias','HSC bias','dendritic bias')
if(!all(ph %in% names(x))) stop('Fig4L input lacks seven phenotype columns')
mat<-as.matrix(x[,ph,drop=FALSE]);rownames(mat)<-x$pathway
# Cluster pathways by their seven phenotype scores, retaining panel row order.
# Pathways are rows (132), phenotype scores are columns (7).
# cor(t(mat)) computes a pathway-by-pathway similarity matrix.
if (nrow(mat) < 2L) stop('Need at least two pathways to cluster')
if (anyDuplicated(rownames(mat))) stop('Duplicate pathway IDs')
if (any(!is.finite(mat))) stop('Nonfinite pathway scores; review export')
co <- stats::cor(t(mat), use='pairwise.complete.obs')
if (!is.matrix(co) || nrow(co)!=nrow(mat) || ncol(co)!=nrow(mat))
  stop('Expected a square pathway correlation matrix')
co[!is.finite(co)] <- 0
diag(co) <- 1
# IMPORTANT: pmax(0, matrix) drops dimensions in some R versions.
# Preserve matrix dimensions explicitly for stats::as.dist().
d <- 1 - co
d[!is.finite(d)] <- 1
d[d < 0] <- 0
diag(d) <- 0
ord <- rownames(mat)[stats::hclust(stats::as.dist(d))$order]
message('Fig4L: clustered ', length(ord), ' pathways across ', ncol(mat), ' phenotypes')
long<-do.call(rbind,lapply(ph,function(p) data.frame(pathway=x$pathway,
                  phenotype=p,score=as.numeric(x[[p]]),stringsAsFactors=FALSE)))
long$pathway<-factor(long$pathway,levels=ord)
long$phenotype<-factor(long$phenotype,levels=rev(ph))
p<-ggplot2::ggplot(long,ggplot2::aes(x=pathway,y=phenotype,fill=score))+
 ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low='#507fa9',mid='white',high='#c85244',
 midpoint=0,limits=c(-25,25),oob=scales::squish,name='signed -ln(P)')+
 ggplot2::labs(x='Pathways (clustered; category strip not yet validated)',y=NULL,
 title='Fig 4L numerical review — 132 pathway score profiles')+
 ggplot2::theme_minimal(base_size=12)+
 ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),
 panel.grid=ggplot2::element_blank())
ggplot2::ggsave(file.path(figure_dir,'Fig4L_REVIEW.pdf'),p,width=12,height=4.4,useDingbats=FALSE)
message('Created Fig4L_REVIEW.pdf. Category annotation intentionally omitted pending accepted-color audit.')
