# Fig 4J/K: HSC signature gene-level association scores and ROC curves.
# Run from CH_publication_release root; no association-model fitting.
pub <- 'Association_analysis/publication'
sd_dir <- file.path(pub,'source_data'); val_dir <- file.path(pub,'validation')
dir.create(sd_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(val_dir,recursive=TRUE,showWarnings=FALSE)
checkpoint <- 'Association_analysis/data/processed/default/association_summary.RData'
aml_file <- 'RNA_analysis/data/genesets/processed/aml_sets.rds'
output_file <- 'RNA_analysis/data/genesets/processed/HSC_output_RF.xlsx'
for(f in c(checkpoint,aml_file,output_file)) if(!file.exists(f)) stop('Missing: ',f)
if(!requireNamespace('readxl',quietly=TRUE)) stop('Install readxl to parse the HSC-output signature workbook.')
if(!requireNamespace('pROC',quietly=TRUE)) stop('Install pROC to reconstruct the historical Fig 4K ROC curves.')

e <- new.env(parent=baseenv()); loaded <- load(checkpoint,envir=e)
if(!all(c('gsea.clone','test.basic.list') %in% loaded)) stop('Association checkpoint lacks expected objects')
sc <- as.matrix(e$gsea.clone$gscore.clone)
if(is.null(rownames(sc)) || anyDuplicated(rownames(sc))) stop('Invalid gene-score matrix rownames')
ph <- names(e$test.basic.list)
if(is.null(ph) || length(ph)!=ncol(sc) || anyDuplicated(ph)) stop('Unexpected association phenotype names')
colnames(sc) <- ph
wanted <- c('HSC bias','growth rate','clone size')
if(!all(wanted %in% colnames(sc))) stop('Missing gene-score phenotypes: ',paste(setdiff(wanted,colnames(sc)),collapse=','))

aml <- readRDS(aml_file)
if(!all(c('aml_sign','kentsis_sign') %in% names(aml))) stop('aml_sets.rds lacks aml_sign or kentsis_sign')
a <- aml$aml_sign; k <- aml$kentsis_sign
if(is.null(rownames(a)) || !'Zhang' %in% names(a)) stop('aml_sign lacks named genes and Zhang HSC groups')
if(is.null(names(k))) stop('kentsis_sign must be a named vector')
# Reproduce historical read_excel / header-promotion logic without inferring column positions.
h <- as.data.frame(readxl::read_excel(output_file),stringsAsFactors=FALSE,check.names=FALSE)
req <- c('gene_names','logfoldchange','adj.p-value')
if(!all(req %in% names(h))) {
  if(nrow(h)<2) stop('HSC output workbook too short to promote header row')
  names(h) <- as.character(unlist(h[1,,drop=TRUE],use.names=FALSE))
  h <- h[-1,,drop=FALSE]
}
if(!all(req %in% names(h))) stop('HSC output workbook unexpected headers: ',paste(names(h),collapse=', '))
h$gene_names <- as.character(h$gene_names)
h$logfoldchange <- suppressWarnings(as.numeric(h$logfoldchange))
h[['adj.p-value']] <- suppressWarnings(as.numeric(h[['adj.p-value']]))
sig <- !is.na(h$logfoldchange) & !is.na(h[['adj.p-value']]) &
       abs(h$logfoldchange)>2 & h[['adj.p-value']] < 1e-3
# Gene direction agrees with the historical 3_association_analysis.R script.
signatures <- list(
  'HSC dormancy'=list(active=rownames(a)[!is.na(a$Zhang) & a$Zhang=='aHSC'],
                      dormant=rownames(a)[!is.na(a$Zhang) & a$Zhang=='dHSC']),
  'HSC label-retention'=list(negative=names(k)[!is.na(k) & k == -2],
                             positive=names(k)[!is.na(k) & k == 2]),
  'HSC output'=list(high=toupper(h$gene_names[sig & h$logfoldchange<0]),
                    low=toupper(h$gene_names[sig & h$logfoldchange>0]))
)
rows <- list(); audit <- list(); m <- 0L; q <- 0L
for(sname in names(signatures)) {
  pairs <- signatures[[sname]]
  # Historical script assigns group 1 then group 2; group 2 wins any overlap.
  g <- stats::setNames(rep(NA_character_,nrow(sc)),rownames(sc))
  overlap <- length(intersect(toupper(pairs[[1]]),toupper(pairs[[2]])))
  for(group in names(pairs)) {
    ids <- intersect(toupper(as.character(pairs[[group]])),rownames(sc))
    g[ids] <- group
  }
  for(group in names(pairs)) {
    ids <- names(g)[!is.na(g) & g==group]
    q<-q+1L
    audit[[q]] <- data.frame(signature=sname,group=group,genes_in_matrix=length(ids),
                             overlaps_between_original_groups=overlap,stringsAsFactors=FALSE)
    if(length(ids)<2L) stop(sname, '/',group,' has only ',length(ids),' scored genes: inspect inputs')
    for(phenotype in wanted) {
      m<-m+1L
      rows[[m]] <- data.frame(gene=ids,signature=sname,group=group,phenotype=phenotype,
                            association_score=as.numeric(sc[ids,phenotype]),stringsAsFactors=FALSE)
    }
  }
}
dat <- do.call(rbind,rows)
dat <- dat[is.finite(dat$association_score),,drop=FALSE]
if(!nrow(dat)) stop('No finite HSC signature association scores')
con<-gzfile(file.path(sd_dir,'Fig4JK_gene_signature_association_scores.tsv.gz'),'wt')
tryCatch(utils::write.table(dat,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA'),finally=close(con))
utils::write.table(do.call(rbind,audit),file.path(val_dir,'Fig4JK_signature_gene_coverage.tsv'),
                   sep='\t',row.names=FALSE,quote=FALSE)
# Fig J: original code uses a two-sided Wilcoxon rank-sum test of gene scores.
test_rows<-list(); k<-0L
for(phenotype in c('HSC bias','growth rate')) for(sname in names(signatures)) {
  d <- dat[dat$phenotype==phenotype & dat$signature==sname,,drop=FALSE]
  if(length(unique(d$group))!=2L) stop('Expected 2 groups for ',sname)
  wt <- suppressWarnings(stats::wilcox.test(association_score~group,data=d))
  gtab <- table(d$group)
  k<-k+1L
  test_rows[[k]]<-data.frame(phenotype=phenotype,signature=sname,
    group_1=names(gtab)[1],n_genes_1=unname(gtab[1]),
    group_2=names(gtab)[2],n_genes_2=unname(gtab[2]),
    wilcoxon_two_sided_pvalue=wt$p.value,p_adjusted=FALSE,stringsAsFactors=FALSE)
}
utils::write.table(do.call(rbind,test_rows),file.path(val_dir,'Fig4J_wilcoxon_summary.tsv'),
                   sep='\t',row.names=FALSE,quote=FALSE)
# Fig K: reproduce historical pROC::roc(..., direction='auto') for each signature/phenotype.
roc_rows<-list(); auc_rows<-list(); k<-0L
for(sname in names(signatures)) for(phenotype in wanted) {
  d<-dat[dat$signature==sname & dat$phenotype==phenotype,,drop=FALSE]
  d$group<-factor(d$group,levels=sort(unique(d$group)))
  rr<-pROC::roc(d$group,d$association_score,quiet=TRUE,direction='auto')
  k<-k+1L
  auc_rows[[k]]<-data.frame(signature=sname,phenotype=phenotype,auc=as.numeric(pROC::auc(rr)),
    n_genes=nrow(d),class_1=levels(d$group)[1],class_2=levels(d$group)[2],
    roc_direction=rr$direction,stringsAsFactors=FALSE)
  roc_rows[[k]]<-data.frame(signature=sname,phenotype=phenotype,
    false_positive_rate=1-as.numeric(rr$specificities),
    sensitivity=as.numeric(rr$sensitivities),stringsAsFactors=FALSE)
}
utils::write.table(do.call(rbind,auc_rows),file.path(val_dir,'Fig4K_ROC_AUC_summary.tsv'),
                   sep='\t',row.names=FALSE,quote=FALSE)
con<-gzfile(file.path(sd_dir,'Fig4K_ROC_coordinates.tsv.gz'),'wt')
tryCatch(utils::write.table(do.call(rbind,roc_rows),con,sep='\t',row.names=FALSE,quote=FALSE),finally=close(con))
message('Fig4J/K: ',nrow(dat),' gene/phenotype/signature rows; ',length(unique(dat$gene)),' distinct genes.')
message('Fig4J two-sided gene-level Wilcoxon tests:');print(do.call(rbind,test_rows),row.names=FALSE)
message('Fig4K AUCs (classification of gene-signature membership, NOT clone mutation status):')
print(do.call(rbind,auc_rows)[,c('signature','phenotype','auc','n_genes')],row.names=FALSE)
message('NOTE: original y-axis wording "Z-score" needs verification: archived score is signed -ln(gene-association p-value).')
