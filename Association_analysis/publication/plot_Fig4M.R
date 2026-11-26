# Fig 4M: review seven HSC module projections from archived cell scores.
pub<-'Association_analysis/publication';figdir<-file.path(pub,'figures')
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
if(!requireNamespace('ggplot2',quietly=TRUE)) stop('ggplot2 required')
x<-utils::read.delim(file.path(pub,'source_data','Fig4M_HSC_module_projection_cells.tsv.gz'),check.names=FALSE)
ph<-setdiff(names(x),c('cell_id','UMAP_1','UMAP_2'))
if(length(ph)!=7L) stop('Expected 7 module columns')
# Match the original HSC embedding figure: within-module rank deciles mapped blue/grey/red.
y<-do.call(rbind,lapply(ph,function(nm){
 sc<-as.numeric(x[[nm]])
 q<-pmin(10L,1L+floor((rank(sc,ties.method='average')-1)/length(sc)*10))
 data.frame(cell_id=x$cell_id,UMAP_1=x$UMAP_1,UMAP_2=x$UMAP_2,
            phenotype=nm,score=sc,score_decile=q,stringsAsFactors=FALSE)
}))
y$phenotype<-factor(y$phenotype,levels=ph)
# One displayed dot is one HSC cell; Source Data is the underlying continuous score.
p<-ggplot2::ggplot(y,ggplot2::aes(UMAP_1,UMAP_2,color=score_decile))+
 ggplot2::geom_point(size=0.32,alpha=0.75,stroke=0)+
 ggplot2::facet_wrap(~phenotype,ncol=4)+
 ggplot2::scale_color_gradientn(colors=c('#295fa2','#edf0ef','#c43a30'),limits=c(1,10),
                                breaks=c(1,5,10),name='Module rank decile')+
 ggplot2::coord_equal()+ggplot2::labs(x='UMAP 1',y='UMAP 2',title='Fig 4M checkpoint module projections (review)')+
 ggplot2::theme_void(base_size=12)+ggplot2::theme(strip.text=ggplot2::element_text(size=11),legend.position='right')
ggplot2::ggsave(file.path(figdir,'Fig4M_REVIEW.pdf'),p,width=11,height=6.1,useDingbats=FALSE)
message('Created Fig4M_REVIEW.pdf (seven HSC-module projection maps).')
