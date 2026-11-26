# Fig 4L: export checkpoint-derived pathway scores for the HSC gene-module heatmap.
# Uses the 132 variable pathway IDs from the Fig4I review reconstruction.
# WARNING: functional-category assignments and final editorial ordering remain PROVISIONAL.
pub <- 'Association_analysis/publication'
sd_dir <- file.path(pub,'source_data'); val_dir <- file.path(pub,'validation')
dir.create(sd_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(val_dir,recursive=TRUE,showWarnings=FALSE)
checkpoint <- 'Association_analysis/data/processed/default/association_summary.RData'
category_path <- file.path(val_dir,'Fig4I_GO_categories_REVIEW.tsv')
for(f in c(checkpoint,category_path)) if(!file.exists(f)) stop('Missing: ',f,' (run Fig4I first)')
e<-new.env(parent=baseenv());load(checkpoint,envir=e)
mat<-as.matrix(e$gsea.clone$gsea.score)
phenotypes<-names(e$test.basic.list)
if(is.null(phenotypes) || length(phenotypes)!=ncol(mat)) stop('Unexpected phenotype columns')
colnames(mat)<-phenotypes
cat<-utils::read.delim(category_path,stringsAsFactors=FALSE)
if(!all(c('pathway','GO_category_review') %in% names(cat)) || anyDuplicated(cat$pathway)) stop('Invalid Fig4I review-category table')
ids<-as.character(cat$pathway)
if(!length(ids) || !all(ids %in% rownames(mat))) stop('Some reviewed pathway IDs absent from association checkpoint')
# Figure 4L displays the seven clone-phenotype rows; the historical script also
# appended two mutant rows, apparently omitted during final figure assembly.
accepted_order <- c('lymphoid bias','clone size','growth rate','erythroid bias',
                    'myeloid bias','HSC bias','dendritic bias')
if(!all(accepted_order %in% phenotypes)) stop('Expected seven published phenotype labels absent')
out<-data.frame(pathway=ids,GO_category_REVIEW=cat$GO_category_review,
                mat[ids,accepted_order,drop=FALSE],check.names=FALSE,stringsAsFactors=FALSE)
if(!all(is.finite(as.matrix(out[,accepted_order,drop=FALSE])))) stop('Missing/nonfinite pathway score in selected matrix')
con<-gzfile(file.path(sd_dir,'Fig4L_pathway_score_matrix_REVIEW.tsv.gz'),'wt')
tryCatch(utils::write.table(out,con,sep='\t',row.names=FALSE,quote=FALSE),finally=close(con))
utils::write.table(data.frame(checkpoint=checkpoint,pathway_set_source=category_path,
    n_selected_pathways=nrow(out),n_phenotype_rows=length(accepted_order),
    caveat='GO-category membership is reconstructed, not independently validated against accepted artwork',
    score='signed -ln(unadjusted fgsea p), sign from NES',stringsAsFactors=FALSE),
    file.path(val_dir,'Fig4L_export_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
message('Fig4L exported ',nrow(out),' pathway x ',length(accepted_order),' phenotype scores.')
message('GO categories and order remain REVIEW ONLY; do not use as final journal color annotations yet.')
