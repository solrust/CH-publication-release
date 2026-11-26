# Fig 4M: reuse archived HSC axes projections + Fig3A shared UMAP coordinates.
# Does not recalculate PCA, regression projections, or the 7.4 GiB scrna checkpoint.
pub<-'Association_analysis/publication';sd_dir<-file.path(pub,'source_data');val_dir<-file.path(pub,'validation')
dir.create(sd_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(val_dir,recursive=TRUE,showWarnings=FALSE)
axes_file <- 'Association_analysis/data/processed/hsc/hsc_axes.RData'
umap_file <- 'RNA_analysis/publication/source_data/Fig3A_umap_cells.tsv.gz'
for(f in c(axes_file,umap_file)) if(!file.exists(f)) stop('Missing input: ',f)
e<-new.env(parent=baseenv());loaded<-load(axes_file,envir=e)
if(!'pheno.projs' %in% loaded) stop('Saved HSC axes checkpoint lacks pheno.projs')
pr<-as.matrix(e$pheno.projs)
if(!is.numeric(pr) || ncol(pr)!=7L || is.null(rownames(pr)) || anyDuplicated(rownames(pr))) {
 stop('Unexpected HSC projection matrix shape or cell identifiers')
}
if(is.null(colnames(pr)) || anyDuplicated(colnames(pr))) stop('Missing phenotype labels on HSC projections')
ref<-utils::read.delim(umap_file,check.names=FALSE,stringsAsFactors=FALSE)
norm<-tolower(gsub('[^[:alnum:]]','',names(ref)))
pick<-function(choices,label) {
 idx<-match(choices,norm,nomatch=0L);idx<-idx[idx>0L]
 if(!length(idx)) stop('Missing ',label,' in Fig3A Source Data. Available: ',paste(names(ref),collapse=', '))
 names(ref)[idx[1L]]
}
cellcol<-pick(c('cellid','cellbarcode','barcode','cell','cellname'),'cell IDs')
xcol<-pick(c('umap1','umapx','embedding1','emb1','dimension1','dim1','x'),'UMAP 1')
ycol<-pick(c('umap2','umapy','embedding2','emb2','dimension2','dim2','y'),'UMAP 2')
refid<-as.character(ref[[cellcol]])
if(anyNA(refid) || anyDuplicated(refid)) stop('Fig3A cell ID missing or duplicated')
where<-match(rownames(pr),refid)
matched<-which(!is.na(where))
if(length(matched) < 0.9*nrow(pr)) stop('Only ',length(matched),'/',nrow(pr),' HSC projections match Fig3A. Review cell ID conventions.')
pr1<-pr[matched,,drop=FALSE]
sel<-where[matched]
out<-data.frame(cell_id=rownames(pr1),
                UMAP_1=as.numeric(ref[[xcol]][sel]),
                UMAP_2=as.numeric(ref[[ycol]][sel]),
                pr1,check.names=FALSE,stringsAsFactors=FALSE)
if(any(!is.finite(as.matrix(out[,-1,drop=FALSE])))) stop('Nonfinite projection or embedding values')
con<-gzfile(file.path(sd_dir,'Fig4M_HSC_module_projection_cells.tsv.gz'),'wt')
tryCatch(utils::write.table(out,con,sep='\t',row.names=FALSE,quote=FALSE),finally=close(con))
summary<-data.frame(input_checkpoint=axes_file,UMAP_reference=umap_file,
    n_archived_HSC_projection_cells=nrow(pr),n_matched_UMAP_cells=nrow(out),
    n_unmatched=nrow(pr)-nrow(out),n_modules=ncol(pr),
    score_definition='projection of seven clone-phenotype gene signatures onto 10 HSC PCA axes',
    stringsAsFactors=FALSE)
utils::write.table(summary,file.path(val_dir,'Fig4M_export_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
message('Fig4M: ',nrow(out),' HSC cells with UMAP coordinates and ',ncol(pr),' module projection scores; ',summary$n_unmatched,' missing from Fig3A UMAP.')
message('Projection names: ',paste(colnames(pr),collapse=', '))
