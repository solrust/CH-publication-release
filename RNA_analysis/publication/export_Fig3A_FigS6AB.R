# CH publication export: Fig3A and Supplementary Fig6A/B
# Run from CH_publication_release root. Uses SAVED scRNA checkpoint, not Pagoda2.
# Dependencies: base R, utils (no cacoa, pagoda2, or fgsea required).
stopifnot(file.exists('RNA_analysis/scripts/4_sc_analysis.R'))

resolve_data <- function(...) {
  tail <- file.path(...)
  root <- Sys.getenv('CH_DATA_ROOT', unset = '')
  candidates <- file.path('RNA_analysis', 'data', tail)
  if (nzchar(root)) candidates <- c(
    file.path(root, 'RNA_analysis', 'data', tail),
    file.path(root, 'data', tail), candidates)
  matches <- candidates[file.exists(candidates)]
  if (!length(matches)) stop('Input not found: ', tail, '\nChecked: ', paste(candidates, collapse = '; '))
  matches[1]
}

sd_dir <- 'RNA_analysis/publication/source_data'
val_dir <- 'RNA_analysis/publication/validation'
dir.create(sd_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(val_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- resolve_data('processed', 'scrna.RData')
meta_file <- resolve_data('processed', 'hashtag_meta.rds')
message('Loading archived checkpoint (large): ', input_file)
e <- new.env(parent = emptyenv())
load(input_file, envir = e)
needed <- c('emb', 'clust.final', 'sampleid', 'clust.final.cols')
stopifnot(all(vapply(needed, exists, logical(1), envir = e, inherits = FALSE)))
emb <- e$emb
clust.final <- e$clust.final
sampleid <- e$sampleid
clust.final.cols <- e$clust.final.cols
# Release archived matrices early to avoid unnecessary peak RAM in later processing.
rm(e); gc(verbose = FALSE)
hash.match <- as.data.frame(readRDS(meta_file))
stopifnot(all(c('genotype','exp') %in% names(hash.match)),
          !is.null(rownames(hash.match)),
          is.matrix(emb), ncol(emb) >= 2,
          !is.null(rownames(emb)), !is.null(names(clust.final)),
          !is.null(names(sampleid)), !anyDuplicated(rownames(emb)),
          !anyDuplicated(names(clust.final)), !anyDuplicated(names(sampleid)))

cells <- rownames(emb)
state <- as.integer(clust.final[cells])
samples <- as.character(sampleid[cells])
meta_index <- match(samples, rownames(hash.match))
expid <- as.character(hash.match$exp[meta_index])
genotype <- as.character(hash.match$genotype[meta_index])
stopifnot(all(names(clust.final) %in% cells),
          all(names(clust.final) %in% names(sampleid)))

umap <- data.frame(
  cell_id = cells,
  umap_1 = as.numeric(emb[,1]),
  umap_2 = as.numeric(emb[,2]),
  cell_state = state,
  sample_id = samples,
  experiment = expid,
  genotype = genotype,
  is_annotated = !is.na(state),
  stringsAsFactors = FALSE
)
# Add explicit cell-state palette from the archived R object.
# Use as.character() so state numbers are interpreted as names, not numeric indices.
umap$state_color <- unname(clust.final.cols[as.character(umap$cell_state)])

umap_path <- file.path(sd_dir, 'Fig3A_umap_cells.tsv.gz')
write.table(umap, gzfile(umap_path), sep = '\t', quote = FALSE,
            row.names = FALSE, na = '')

# IMPORTANT: reproduce original cell-state normalization literally:
# (fraction of each state in a sample)/(fraction of that state in experiment).
# Only sample/state pairs with observed cells are included, as in original code.
clust <- factor(clust.final)
nms <- names(clust)
sub_sample <- as.character(sampleid[nms])
sub_idx <- match(sub_sample, rownames(hash.match))
clust.info <- data.frame(
  clust = clust,
  sample = sub_sample,
  genotype = as.character(hash.match$genotype[sub_idx]),
  exp = as.character(hash.match$exp[sub_idx]),
  count = 1L,
  stringsAsFactors = FALSE
)
if (anyNA(clust.info[,c('sample','genotype','exp')])) {
  stop('Missing sample/genotype/experiment metadata for annotated cells')
}
clust.aggr <- aggregate(count ~ clust + sample + genotype + exp,
                         data = clust.info, sum)
clust.aggr.tot <- aggregate(count ~ sample, data = clust.aggr, sum)
rownames(clust.aggr.tot) <- clust.aggr.tot$sample
clust.aggr.clust <- aggregate(count ~ clust + exp, data = clust.aggr, sum)
clust.aggr.clust.exp <- aggregate(count ~ exp,
                                   data = clust.aggr.clust, sum)
rownames(clust.aggr.clust.exp) <- clust.aggr.clust.exp$exp
clust.aggr.clust$freq <- clust.aggr.clust$count /
  clust.aggr.clust.exp[as.character(clust.aggr.clust$exp), 'count']
clust.aggr$sample_total <- clust.aggr.tot[as.character(clust.aggr$sample), 'count']
clust.aggr$within_sample_fraction <- clust.aggr$count / clust.aggr$sample_total
mtch <- match(paste(clust.aggr$clust,clust.aggr$exp,sep='_'),
              paste(clust.aggr.clust$clust,clust.aggr.clust$exp,sep='_'))
stopifnot(!anyNA(mtch))
clust.aggr$experiment_state_baseline <- clust.aggr.clust$freq[mtch]
clust.aggr$freq <- clust.aggr$within_sample_fraction /
  clust.aggr$experiment_state_baseline
stopifnot(all(is.finite(clust.aggr$freq)), all(clust.aggr$freq > 0))

s6a <- data.frame(
  experiment = as.character(clust.aggr$exp),
  sample_id = as.character(clust.aggr$sample),
  genotype = as.character(clust.aggr$genotype),
  cell_state = as.integer(as.character(clust.aggr$clust)),
  n_cells = as.integer(clust.aggr$count),
  n_cells_in_sample = as.integer(clust.aggr$sample_total),
  sample_state_fraction = clust.aggr$within_sample_fraction,
  experiment_state_baseline = clust.aggr$experiment_state_baseline,
  normalized_enrichment = clust.aggr$freq,
  stringsAsFactors = FALSE
)
s6a <- s6a[order(s6a$experiment,s6a$cell_state,s6a$sample_id),]
write.table(s6a, file.path(sd_dir, 'FigS6A_sample_state_enrichment.tsv'),
            sep='\t', quote=FALSE, row.names=FALSE, na='')

# Original unpaired two-sided Welch t-tests and fold-change comparisons
# used by FigS6B. These tests are NOT necessarily the Cacoa results behind
# manually added asterisks in S6A.
exps <- unique(as.character(clust.aggr$exp))
states <- unique(as.character(clust.aggr$clust))
stopifnot(length(exps) == 2L, all(c('123','127') %in% exps),
          length(states) == 10L)
folds <- do.call(rbind, lapply(exps, function(expp) {
  do.call(rbind, lapply(states, function(cl) {
    ind <- as.character(clust.aggr$exp) == expp &
      as.character(clust.aggr$clust) == cl
    s1 <- clust.aggr$freq[ind & clust.aggr$genotype == 'WT']
    s2 <- clust.aggr$freq[ind & clust.aggr$genotype != 'WT']
    if (length(s1) < 2L || length(s2) < 2L) {
      warning('Insufficient sample replication for t.test exp=',expp,' cluster=',cl)
      p <- NA_real_
    } else p <- t.test(s1, s2)$p.value
    data.frame(experiment=expp,cell_state=as.integer(cl),
      genotype=if (expp=='123') 'DNMT3A' else 'TET2',
      n_WT_samples=length(s1),n_mutant_samples=length(s2),
      mean_WT_enrichment=if(length(s1)) mean(s1) else NA_real_,
      mean_mutant_enrichment=if(length(s2)) mean(s2) else NA_real_,
      fold_change=if(length(s1) && length(s2)) mean(s2)/mean(s1) else NA_real_,
      welch_ttest_pvalue=p)
  }))
}))
folds$BH_adjusted_ttest_pvalue <- ave(folds$welch_ttest_pvalue,
                                      folds$experiment,
                                      FUN=function(x) p.adjust(x,'BH'))
folds <- folds[order(folds$experiment,folds$cell_state),]
stopifnot(!anyNA(folds$fold_change), all(is.finite(folds$fold_change)))
write.table(folds, file.path(sd_dir,'FigS6B_state_fold_changes.tsv'),
            sep='\t', quote=FALSE, row.names=FALSE, na='')

# Pearson correlation between the 10 cell-state fold changes
x <- folds$fold_change[folds$experiment=='123'][order(folds$cell_state[folds$experiment=='123'])]
y <- folds$fold_change[folds$experiment=='127'][order(folds$cell_state[folds$experiment=='127'])]
ct <- cor.test(x,y,method='pearson')
corr <- data.frame(method='Pearson', n_cell_states=length(x),
  pearson_r=unname(ct$estimate), t_statistic=unname(ct$statistic),
  degrees_of_freedom=unname(ct$parameter),pvalue=ct$p.value,
  conf_int_lower=unname(ct$conf.int[1]),conf_int_upper=unname(ct$conf.int[2]),
  stringsAsFactors=FALSE)
write.table(corr, file.path(sd_dir,'FigS6B_correlation.tsv'),
            sep='\t', quote=FALSE, row.names=FALSE)

report <- data.frame(
  metric=c('UMAP_cells','annotated_cells','S6A_sample_state_rows',
           'S6B_state_genotype_rows','S6B_pearson_r','S6B_pearson_pvalue',
           'S6A_significance_asterisks_verified'),
  value=c(nrow(umap),sum(umap$is_annotated),nrow(s6a),nrow(folds),
          signif(corr$pearson_r,8),signif(corr$pvalue,8), 'NO'),
  stringsAsFactors=FALSE)
write.table(report,file.path(val_dir,'Fig3A_FigS6AB_export_summary.tsv'),
            quote=FALSE,sep='\t',row.names=FALSE)
cat('Exported Fig3A + S6A/B Source Data:\n'); print(report,row.names=FALSE)
cat('S6A: manual Cacoa/significance asterisks NOT verified by this script.\n')
cat('Fig3B: run export_Fig3B.R next; it uses the exported cell table.\n')
