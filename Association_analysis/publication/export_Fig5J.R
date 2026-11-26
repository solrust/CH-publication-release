# Figure 5J (publication reproduction): repeated cell-level CV Spearman correlations.
# This targets Association_analysis/scripts/7_heritability_updated.R, NOT the
# older 7_heritability.R R-squared checkpoint.
# Run from the CH_publication_release repository root, in a >=64 GB R session.
# A reduced smoke test is available with options(CH.fig5j.splits = 5).

stopifnot(file.exists('RNA_analysis/data/processed/scrna.RData'))
inputs <- list(
  scrna = 'RNA_analysis/data/processed/scrna.RData',
  rna_clones = 'RNA_clone_analysis/data/processed/clone_info.RData',
  dna = 'DNA_analysis/data/processed/DNA_clones.RData'
)
out_dir <- 'Association_analysis/publication'
dir.create(file.path(out_dir, 'source_data'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, 'validation'), recursive = TRUE, showWarnings = FALSE)

message('Loading scRNA checkpoint; allow substantial RAM.')
r <- new.env(parent = emptyenv())
load(inputs$scrna, envir = r)
c <- new.env(parent = emptyenv())
load(inputs$rna_clones, envir = c)
d <- new.env(parent = emptyenv())
load(inputs$dna, envir = d)

for (nm in c('p2pca', 'clust.final')) if (!exists(nm, envir = r, inherits = FALSE)) stop('Missing scrna object: ', nm)
for (nm in c('cell.clone', 'clone.info', 'clone.states2')) if (!exists(nm, envir = c, inherits = FALSE)) stop('Missing clone object: ', nm)
if (!exists('rate', envir = d, inherits = FALSE)) stop('Missing DNA growth-rate object: rate')

# Match the phenotype construction in scripts/1_processing.R.
clons <- rownames(c$clone.states2)
if (!identical(rownames(c$clone.info), rownames(c$clone.states2))) {
  if (!all(clons %in% rownames(c$clone.info))) stop('Clone identifiers differ between clone.info and clone.states2')
}
size_by_sample <- tapply(as.numeric(c$clone.info[clons, 'size']),
                         as.character(c$clone.info[clons, 'id']), sum)
pheno <- cbind(
  'growth.rate' = as.numeric(d$rate[clons]),
  'clone.size' = log(as.numeric(c$clone.info[clons, 'size']) /
                     size_by_sample[as.character(c$clone.info[clons, 'id'])]),
  c$clone.states2[clons, c('1','4','5','8','10'), drop = FALSE]
)
rownames(pheno) <- clons
colnames(pheno) <- c('growth.rate','clone.size','HSC.bias',
                     'lymphoid.bias','dendritic.bias','myeloid.bias','erythroid.bias')

# Match the updated script's data.frame name repair and HSC-cell restriction.
cells <- c$cell.clone
if (is.null(rownames(cells)) || !'cloneid' %in% names(cells)) stop('cell.clone must be rownamed and include cloneid')
cells$genotype <- as.character(c$clone.info[cells$cloneid, 'genotype'])
pheno_cell <- pheno[as.character(cells$cloneid), , drop = FALSE]
colnames(pheno_cell) <- make.names(colnames(pheno_cell))
cell2clone <- data.frame(cells, pheno_cell, check.names = TRUE)
rownames(cell2clone) <- rownames(cells)

# The historical 7_heritability_updated.R uses `!is.na(cell2clone[,9])`.
# We record that prefilter, then require finite target values per phenotype.
if (ncol(cell2clone) < 9) stop('Unexpected structure in original cell2clone data.frame')
prefilter_cells <- rownames(cell2clone)[!is.na(cell2clone$genotype) &
                                          !is.na(cell2clone[[9]])]
pca <- r$p2pca
if (is.null(rownames(pca))) stop('p2pca has no cell rownames')
if (!is.matrix(pca)) pca <- as.matrix(pca)
if (anyDuplicated(rownames(pca))) stop('Duplicated PCA cell identifiers')
state <- r$clust.final
if (is.null(names(state))) stop('clust.final must have cell names')
hsc_cells <- intersect(prefilter_cells, rownames(pca))
hsc_cells <- hsc_cells[!is.na(state[hsc_cells]) & as.character(state[hsc_cells]) == '1']
message('HSC cells with clone metadata and PCA: ', length(hsc_cells))

# In the original script, four lineage phenotypes are fitted, but the
# dendritic-bias phenotype is omitted from the final scatterplot.
targets <- c('lymphoid.bias','dendritic.bias','myeloid.bias','erythroid.bias')
labels <- c('lymphoid.bias'='B Lymph', 'dendritic.bias'='Dn',
            'myeloid.bias'='GMP', 'erythroid.bias'='Ery')
genotypes <- c('WT','TET2','DNMT3A')
nsplits <- getOption('CH.fig5j.splits', 300L)
nsplits <- as.integer(nsplits)
if (length(nsplits) != 1 || is.na(nsplits) || nsplits < 2) stop('CH.fig5j.splits must be >=2')
seed <- as.integer(getOption('CH.fig5j.seed', 20261106L))
set.seed(seed)

# Each fit mirrors r2.cv() in 7_heritability_updated.R:
# 90% of cells for OLS training; Spearman correlation on held-out cells.
records <- vector('list', length(genotypes)*length(targets)*nsplits)
coverage <- list()
k <- 0L
for (geno in genotypes) {
  for (ph in targets) {
    xx <- hsc_cells[as.character(cell2clone[hsc_cells, 'genotype']) == geno]
    y <- suppressWarnings(as.numeric(cell2clone[xx, ph]))
    ok <- is.finite(y) & rowSums(!is.finite(pca[xx, , drop=FALSE])) == 0
    xx <- xx[ok]; y <- y[ok]
    ntrain <- floor(0.9*length(xx))
    if (ntrain < ncol(pca)+3 || (length(xx)-ntrain) < 5) {
      stop('Insufficient ', geno, ' HSC cells for ', ph, ': ', length(xx),
           ' cells; ', ncol(pca), ' PCs')
    }
    nclones <- length(unique(as.character(cell2clone[xx,'cloneid'])))
    nsamples <- if ('id' %in% colnames(cell2clone)) length(unique(as.character(cell2clone[xx,'id']))) else NA_integer_
    coverage[[length(coverage)+1L]] <- data.frame(genotype=geno, phenotype=ph,
                     n_cells=length(xx), n_clones=nclones, n_samples=nsamples,
                     n_PCs=ncol(pca), stringsAsFactors=FALSE)
    message('Fig5J ', geno, ' / ', ph, ': ', length(xx), ' HSC cells, ', nclones,
            ' clones, ', nsplits, ' splits')
    mm <- pca[xx,,drop=FALSE]
    for (i in seq_len(nsplits)) {
      tr <- sample.int(length(xx), size=ntrain, replace=FALSE)
      te <- setdiff(seq_along(xx), tr)
      train <- data.frame(y=y[tr], mm[tr,,drop=FALSE], check.names=TRUE)
      test <- data.frame(mm[te,,drop=FALSE], check.names=TRUE)
      ft <- stats::lm(y ~ ., data=train)
      pred <- as.numeric(stats::predict(ft, newdata=test))
      rho <- suppressWarnings(stats::cor(y[te], pred, method='spearman', use='complete.obs'))
      if (!is.finite(rho)) stop('Non-finite held-out Spearman rho for ', geno, '/', ph, ' split ', i)
      k <- k+1L
      records[[k]] <- data.frame(genotype=geno, phenotype=ph, lineage=unname(labels[ph]),
                      split=i, spearman_rho=rho, n_train=length(tr), n_test=length(te),
                      stringsAsFactors=FALSE)
    }
  }
}
raw <- do.call(rbind, records)
cover <- do.call(rbind, coverage)
split_by <- split(raw, interaction(raw$genotype, raw$phenotype, drop=TRUE))
stat <- do.call(rbind, lapply(split_by, function(z) {
  data.frame(genotype=z$genotype[1], phenotype=z$phenotype[1],
             lineage=z$lineage[1], mean_rho=mean(z$spearman_rho),
             sd_rho=stats::sd(z$spearman_rho), n_splits=nrow(z),
             stringsAsFactors=FALSE)
}))
stat <- merge(stat, cover, by=c('genotype','phenotype'), all.x=TRUE, sort=FALSE)
stat <- stat[order(match(stat$phenotype,targets), match(stat$genotype, genotypes)),]
rownames(stat) <- NULL
if (anyNA(stat)) warning('Some optional sample-count annotations may be NA')

# Distinct filenames for a smoke test: do not overwrite full 300-split outputs.
tag <- if (nsplits == 300L) '' else paste0('_SMOKETEST_',nsplits,'splits')
raw_file <- file.path(out_dir,'source_data',paste0('Fig5J_split_correlations',tag,'.tsv.gz'))
summary_file <- file.path(out_dir,'source_data',paste0('Fig5J_heritability_summary',tag,'.tsv'))
coverage_file <- file.path(out_dir,'validation',paste0('Fig5J_cell_clone_coverage',tag,'.tsv'))
write.table(raw, gzfile(raw_file, 'wt'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(stat, summary_file, sep='\t', quote=FALSE, row.names=FALSE)
write.table(cover, coverage_file, sep='\t', quote=FALSE, row.names=FALSE)
message('Fig5J summary (includes dendritic-bias metric excluded from final panel):')
print(stat[, c('genotype','lineage','mean_rho','sd_rho','n_cells','n_clones','n_samples')], row.names=FALSE)
message('Wrote: ', summary_file)
message('NOTE: original analysis omitted random seed; values may differ modestly from accepted artwork.')
message('NOTE: cell-level splits can contain cells from the same clone in train and test.')
