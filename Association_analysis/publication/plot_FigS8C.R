# S8C uses the exact same archived HSC-module projections as Fig4M.
# Reuses Fig4M Source Data without duplicating 6531 x 7 cell-level values.
if(!requireNamespace('ggplot2',quietly=TRUE))stop('ggplot2 required')
pub <- 'Association_analysis/publication'; fg <- file.path(pub,'figures')
dir.create(fg,recursive=TRUE,showWarnings=FALSE)
f <- file.path(pub,'source_data','Fig4M_HSC_module_projection_cells.tsv.gz')
if(!file.exists(f))stop('First export Fig4M Source Data: ',f)
x <- read.delim(f,check.names=FALSE)
ph <- c('growth rate','clone size','HSC bias','erythroid bias','lymphoid bias','dendritic bias','myeloid bias')
if(!all(c('cell_id','UMAP_1','UMAP_2',ph) %in% names(x))) stop('Missing required Fig4M columns')
arr <- lapply(ph,function(nm){
 s <- x[[nm]]
 # Approximate historical rank color bins using 10 equal ranked groups.
 dec <- pmin(10L,1L+floor((rank(s,ties.method='average')-1)/length(s)*10))
 data.frame(UMAP_1=x$UMAP_1,UMAP_2=x$UMAP_2,phenotype=nm,rank_decile=dec)
})
y <- do.call(rbind,arr)
y$phenotype <- factor(y$phenotype,levels=ph)
p<-ggplot2::ggplot(y,ggplot2::aes(UMAP_1,UMAP_2,color=rank_decile))+
 ggplot2::geom_point(size=0.22,alpha=0.55,stroke=0)+
 ggplot2::facet_wrap(~phenotype,ncol=4)+
 ggplot2::scale_color_gradientn(colors=c('#295fa2','#ededed','#c92b2c'),limits=c(1,10),
                                breaks=c(1,5,10),name='Projection rank decile')+
 ggplot2::coord_equal()+ggplot2::theme_void(base_size=10)+
 ggplot2::theme(strip.text=ggplot2::element_text(size=10))+
 ggplot2::labs(title='S8C: clone phenotype projections on HSC cells (review)')
ggplot2::ggsave(file.path(fg,'FigS8C_reproduced_REVIEW.pdf'),p,width=10.5,height=6.2,useDingbats=FALSE)
cat('S8C: reused ',nrow(x),' HSC cells and ',length(ph),' projections from Fig4M Source Data.\n',sep='')
