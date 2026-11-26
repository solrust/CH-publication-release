# Source Data for Supplementary Figure 7B-D
# Run with working directory CH_publication_release (repository root).
# Reuses archived barcode-calling results; does not rerun LARRY processing or load scrna.RData.

local({
  pub <- 'RNA_clone_analysis/publication'
  src <- file.path(pub, 'source_data')
  val <- file.path(pub, 'validation')
  dir.create(src, showWarnings = FALSE, recursive = TRUE)
  dir.create(val, showWarnings = FALSE, recursive = TRUE)

  data_dir <- 'RNA_clone_analysis/data/processed'
  clone_info_file <- file.path(data_dir, 'clone_info.RData')
  hashtag_file <- 'RNA_analysis/data/processed/hashtag_meta.rds'
  if (!file.exists(clone_info_file) || !file.exists(hashtag_file)) {
    stop('Missing clone_info.RData or hashtag_meta.rds; run from repository root.')
  }
  expids <- c('exp123','exp127')
  rnas <- paste0('RNA',1:4)
  inpaths <- unlist(lapply(expids, function(ex) file.path(data_dir,paste0('larry_',ex,'_',rnas,'.rds'))))
  if (any(!file.exists(inpaths))) stop('Missing LARRY RDS: ',paste(inpaths[!file.exists(inpaths)],collapse=', '))

  # S7B: one observation per cell/barcode and RNA library, as in 2_larry_analysis.R.
  barcode_tables <- vector('list',length(inpaths))
  lib_summary <- vector('list',length(inpaths))
  for (i in seq_along(inpaths)) {
    path <- inpaths[i]
    x <- readRDS(path)
    if (is.null(x$mat.clones)) stop('LARRY object lacks mat.clones: ', path)
    m <- x$mat.clones
    if (is.null(ncol(m)) || ncol(m)==0) stop('Invalid mat.clones: ',path)
    # Historical calculation: dominant = column max; total = column sum.
    max_umi <- apply(m, 2L, max)
    total_umi <- if (inherits(m,'Matrix')) Matrix::colSums(m) else colSums(m)
    stopifnot(length(max_umi)==length(total_umi), all(max_umi<=total_umi), all(max_umi>=0))
    lib <- sub('^larry_','',sub('\\.rds$','',basename(path)))
    bc <- colnames(m)
    if (is.null(bc)) bc <- paste0('cell_',seq_along(max_umi))
    barcode_tables[[i]] <- data.frame(
      library=lib, barcode=as.character(bc),
      total_clone_umis=as.numeric(total_umi),
      dominant_clone_umis=as.numeric(max_umi),
      other_clone_umis=as.numeric(total_umi-max_umi),
      stringsAsFactors=FALSE)
    lib_summary[[i]] <- data.frame(library=lib, n_barcodes=ncol(m),
                                   n_total_umi_ge_5=sum(total_umi>=5),
                                   n_dominant_fraction_ge_085=sum(total_umi>0 & max_umi >= 0.85*total_umi),
                                   stringsAsFactors=FALSE)
    message('S7B ',lib,': ',ncol(m),' barcodes')
    rm(x,m); gc(verbose=FALSE)
  }
  barcode_dat <- do.call(rbind,barcode_tables)
  write.table(barcode_dat, gzfile(file.path(src,'FigS7B_dominant_clone_umis.tsv.gz'),'wt'),
              sep='\t',quote=FALSE,row.names=FALSE)
  write.table(do.call(rbind,lib_summary),file.path(val,'FigS7B_library_summary.tsv'),
              sep='\t',quote=FALSE,row.names=FALSE)
  message('S7B exported ',nrow(barcode_dat),' barcode rows.')

  # S7C & S7D: per-clone sizes and per-sample clone numbers from the saved
  # filtered cell.clone table (same observation level as historical script).
  e <- new.env(parent=baseenv()); load(clone_info_file, envir=e)
  if (!exists('cell.clone',envir=e,inherits=FALSE)) stop('clone_info.RData lacks cell.clone')
  cl <- e$cell.clone
  req <- c('id','expid','cloneid')
  if (!all(req %in% names(cl))) stop('Missing cell.clone fields: ',paste(setdiff(req,names(cl)),collapse=', '))
  hash.match <- readRDS(hashtag_file)
  ids <- as.character(cl$id)
  geno <- as.character(hash.match[ids,'genotype'])
  if (length(geno)!=nrow(cl)||anyNA(geno)) stop('Unmatched genotype in hashtag_meta.rds')
  df <- data.frame(sample_id=ids,expid=as.character(cl$expid),
                   genotype=geno,clone_id=as.character(cl$cloneid),
                   stringsAsFactors=FALSE)
  if (anyNA(df)) stop('Missing fields in clone assignments')
  per_clone <- aggregate(rep.int(1L,nrow(df)),
                         by=df[c('sample_id','expid','genotype','clone_id')],FUN=sum)
  names(per_clone)[5] <- 'clone_size'
  per_clone <- per_clone[order(per_clone$expid,per_clone$sample_id,per_clone$clone_id),]
  write.table(per_clone,file.path(src,'FigS7C_clone_sizes_by_sample.tsv'),
              sep='\t',quote=FALSE,row.names=FALSE)

  per_sample <- aggregate(rep.int(1L,nrow(per_clone)),
                          by=per_clone[c('sample_id','expid','genotype')],FUN=sum)
  names(per_sample)[4] <- 'n_clones'
  per_sample <- per_sample[order(per_sample$expid,per_sample$genotype,per_sample$sample_id),]
  write.table(per_sample,file.path(src,'FigS7D_clone_counts_by_sample.tsv'),
              sep='\t',quote=FALSE,row.names=FALSE)
  message('S7C exported ',nrow(per_clone),' clone-by-sample records.')
  message('S7D exported ',nrow(per_sample),' sample-level clone counts:')
  print(per_sample,row.names=FALSE)

  # The original exploratory script tests WT vs mutant within exp127,
  # but the accepted S7D artwork does not show a test p-value.
  wt <- per_sample$n_clones[per_sample$expid=='127' & per_sample$genotype=='WT']
  mut <- per_sample$n_clones[per_sample$expid=='127' & per_sample$genotype=='TET2']
  pv <- if(length(wt)>0 && length(mut)>0) suppressWarnings(stats::wilcox.test(wt,mut)$p.value) else NA_real_
  notes <- data.frame(test='historical_exp127_WT_vs_TET2_unpaired_Wilcoxon',
                      n_wt_samples=length(wt),n_tet2_samples=length(mut),
                      p_value=pv,displayed_in_panel=FALSE)
  write.table(notes,file.path(val,'FigS7D_historical_test.tsv'),
              sep='\t',quote=FALSE,row.names=FALSE)
  message('S7B-D Source Data export finished. Do not claim visual verification until inspecting plots.')
})
