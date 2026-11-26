# Figure 5E and Supplementary Figure 8D: checkpoint-based source-data export.
# Run from CH_publication_release root. No DESeq2 / mutation prediction rerun.
# This is an INITIAL PASS retaining original ROC cohort definitions; see README.
local({
  pub <- 'RNA_clone_analysis/publication'
  src <- file.path(pub, 'source_data')
  val <- file.path(pub, 'validation')
  dir.create(src, recursive=TRUE, showWarnings=FALSE)
  dir.create(val, recursive=TRUE, showWarnings=FALSE)
  clone_path <- 'RNA_clone_analysis/data/processed/clone_info.RData'
  scores_dir <- 'RNA_clone_analysis/data/processed'
  required <- c(clone_path, file.path(scores_dir,c('DNMT3A_scores.rds','TET2_scores.rds')))
  if (any(!file.exists(required))) stop('Missing input(s): ',paste(required[!file.exists(required)],collapse=', '))
  if (!requireNamespace('pROC',quietly=TRUE)) stop('Install pROC for ROC curve calculations.')
  e <- new.env(parent=baseenv()); load(clone_path,envir=e)
  if (!all(c('clone.info','cell.clone') %in% ls(e))) stop('clone_info.RData missing expected objects.')
  clones <- e$clone.info; assignments <- e$cell.clone
  if (!'cloneid' %in% names(assignments) || !'genotype' %in% names(clones) || !'size' %in% names(clones))
    stop('Unexpected clone_info columns.')
  if (anyDuplicated(rownames(assignments))) stop('Duplicated cell IDs in cell.clone')
  cells_by_clone <- split(rownames(assignments),as.character(assignments$cloneid))
  write_gz <- function(df,path){
    con <- gzfile(path,open='wt')
    tryCatch(write.table(df,con,sep='\t',quote=FALSE,row.names=FALSE,na='NA'),
             finally=close(con))
  }
  cell_all <- list(); clone_all <- list(); roc_all <- list(); summary_all <- list()
  for (gen in c('DNMT3A','TET2')) {
    checkpoint <- readRDS(file.path(scores_dir,paste0(gen,'_scores.rds')))
    if (is.null(checkpoint$mut.sig.de.qq) || length(checkpoint$mut.sig.de.qq)!=10)
      stop('Expected ten state-specific mut.sig.de.qq tables in ',gen,' checkpoint.')
    tbls <- lapply(seq_along(checkpoint$mut.sig.de.qq),function(k){
      z <- checkpoint$mut.sig.de.qq[[k]]
      needed <- c('score','pheno','use','qmut','qwt')
      if (!all(needed %in% colnames(z))) stop(gen,' state ',k,' lacks: ',paste(setdiff(needed,colnames(z)),collapse=','))
      if (is.null(rownames(z))) stop(gen,' state ',k,' has no cell IDs in row names.')
      data.frame(contrast=gen,cell_state=k,cell_id=rownames(z),
                 score=z$score,pheno=as.numeric(z$pheno),
                 split=as.character(z$use),qmut=as.numeric(z$qmut),qwt=as.numeric(z$qwt),
                 stringsAsFactors=FALSE)
    })
    cell <- do.call(rbind,tbls)
    rownames(cell)<-NULL
    if (anyDuplicated(cell$cell_id)) stop('Duplicated cell IDs in ',gen,' score tables.')
    if (!all(cell$split %in% c('train','test'))) stop('Unexpected train/test value in ',gen)
    nmatch <- sum(cell$cell_id %in% rownames(assignments))
    if (nmatch!=nrow(cell)) stop(gen,': only ',nmatch,'/',nrow(cell),' score cells match cell.clone rows.')
    cell_all[[gen]] <- cell
    idx_by_cell <- setNames(seq_len(nrow(cell)),cell$cell_id)
    clone_tbl <- do.call(rbind,lapply(seq_len(nrow(clones)),function(i){
      clid <- rownames(clones)[i]
      cell_ids <- cells_by_clone[[clid]]
      k <- if (is.null(cell_ids)) integer(0) else unname(idx_by_cell[cell_ids])
      k <- k[!is.na(k)]
      take <- cell[k,,drop=FALSE]
      tr <- take$qmut[take$split=='train']; te <- take$qmut[take$split=='test']
      data.frame(contrast=gen,clone_id=clid,
                 clone_genotype=as.character(clones$genotype[i]),
                 clone_size=as.numeric(clones$size[i]),
                 n_train=length(tr),n_test=length(te),
                 score_train=if(length(tr)) mean(tr) else NA_real_,
                 score_test=if(length(te)) mean(te) else NA_real_,
                 stringsAsFactors=FALSE)
    }))
    clone_all[[gen]] <- clone_tbl
    # EXACT historical ROC response: all non-WT clones are treated as cases;
    # DNMT3A and TET2 curves differ by score, not by the clone-level case set.
    for (level in c('cell','clone')) for (split in c('train','test')) {
      if(level=='cell') {
        data <- cell; x <- data$qwt; y <- data$pheno
        keep <- data$split==split
      } else {
        data <- clone_tbl
        x <- if(split=='train') data$score_train else data$score_test
        y <- as.numeric(data$clone_genotype!='WT')
        keep <- data$clone_size>=10
      }
      ok <- keep & is.finite(x) & !is.na(y)
      if(sum(ok)<3 || length(unique(y[ok]))!=2) {
        warning('Skipping ROC with insufficient classes: ',gen,' ',level,' ',split)
        next
      }
      r <- pROC::roc(response=y[ok],predictor=x[ok],quiet=TRUE)
      roc_all[[paste(gen,level,split,sep='_')]] <- data.frame(
        contrast=gen,level=level,split=split,
        specificity=as.numeric(r$specificities),sensitivity=as.numeric(r$sensitivities),
        stringsAsFactors=FALSE)
      summary_all[[paste(gen,level,split,sep='_')]] <- data.frame(
        contrast=gen,level=level,split=split,
        n_total=sum(ok),n_control=sum(y[ok]==0),n_case=sum(y[ok]==1),
        auc=as.numeric(pROC::auc(r)),direction=r$direction,
        note=if(level=='clone') 'historical_cohort_all_non_WT_vs_WT' else 'mutation_specific_cell_label',
        stringsAsFactors=FALSE)
    }
    message(gen,': ',nrow(cell),' scored cells; ',nrow(clone_tbl),' clone rows')
  }
  if(!length(roc_all)) stop('No ROC curves calculated.')
  write_gz(do.call(rbind,cell_all),file.path(src,'Fig5E_cell_scores.tsv.gz'))
  write_gz(do.call(rbind,clone_all),file.path(src,'Fig5E_FigS8D_clone_scores.tsv.gz'))
  write_gz(do.call(rbind,roc_all),file.path(src,'Fig5E_ROC_coordinates.tsv.gz'))
  summary <- do.call(rbind,summary_all)
  write.table(summary,file.path(val,'Fig5E_ROC_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  s8d <- do.call(rbind,lapply(clone_all,function(x){
    sub <- x[x$clone_size>50 & x$clone_genotype %in% c('WT',x$contrast[1]),,drop=FALSE]
    data.frame(contrast=x$contrast[1],genotype=c('WT','mutant'),n_clones=c(sum(sub$clone_genotype=='WT'),sum(sub$clone_genotype==x$contrast[1])),stringsAsFactors=FALSE)
  }))
  write.table(s8d,file.path(val,'FigS8D_clone_count_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  print(summary,row.names=FALSE)
  message('Fig5E/S8D exported. Compare AUCs and visual plots with accepted artwork before validation.')
})
