# Fig. 4E (LSP1 / clone HSC bias) and Fig. 4G (NPM1 / clone growth rate)
# Targeted clone-level expression export; no association-model refitting.
# Run from CH_publication_release root, using a 64-GB R session.
# Mirrors Association_analysis/scripts/1_processing.R:
#   clone.counts[[3]] -> clone filters -> gene prevalence -> fpm normalization.
# Uses sparse matrix multiplication instead of converting all cell counts to dense.

stopifnot(dir.exists('Association_analysis'), dir.exists('RNA_analysis'))
pub <- 'Association_analysis/publication'
sd <- file.path(pub, 'source_data')
vd <- file.path(pub, 'validation')
dir.create(sd, recursive=TRUE, showWarnings=FALSE)
dir.create(vd, recursive=TRUE, showWarnings=FALSE)

sources <- c(
  sc='RNA_analysis/data/processed/scrna.RData',
  clones='RNA_clone_analysis/data/processed/clone_info.RData',
  dna='DNA_analysis/data/processed/DNA_clones.RData'
)
for (f in sources) if (!file.exists(f)) stop('Missing required input: ', f)
if (!requireNamespace('Matrix', quietly=TRUE)) stop('Matrix package required')

message('Loading scrna.RData (~7.4 GiB disk; use sufficient RAM)')
rna <- new.env(parent=emptyenv())
load(sources['sc'], envir=rna)
clone <- new.env(parent=emptyenv())
load(sources['clones'], envir=clone)
dna <- new.env(parent=emptyenv())
load(sources['dna'], envir=dna)

for (nm in c('cd','p2clusts')) if (!exists(nm, envir=rna, inherits=FALSE)) stop('Missing scrna object: ', nm)
for (nm in c('cell.clone','clone.info','clone.states2')) if (!exists(nm, envir=clone, inherits=FALSE)) stop('Missing clone object: ', nm)
if (!exists('rate', envir=dna, inherits=FALSE)) stop('DNA_clones.RData lacks rate vector')

cd <- rna$cd
clusts <- rna$p2clusts$PCAncc$leiden
cell.clone <- clone$cell.clone
clone.info <- clone$clone.info
clone.states2 <- clone$clone.states2
rate <- dna$rate

if (is.null(names(clusts))) stop('Leiden cluster labels must be named by cell barcode')
if (is.null(colnames(cd)) || is.null(rownames(cd))) stop('cd must have named gene rows and cell columns')
if (is.null(names(rate))) stop('rate must have clone IDs as names')
if (!'cloneid' %in% names(cell.clone)) stop('cell.clone lacks cloneid')
if (!all(c('genotype','size','expid','id') %in% colnames(clone.info))) {
  stop('clone.info must have genotype, size, expid, and id columns; inspect object')
}

# The historical script used clone.counts[[3]], with its 3rd Leiden-cluster
# index. We explicitly select cluster label "3" (not final annotated cell state).
cluster_cells <- names(clusts)[!is.na(clusts) & as.character(clusts)=='3']
cells <- intersect(rownames(cell.clone), intersect(cluster_cells, colnames(cd)))
if (!length(cells)) stop('No matching clone-assigned cells in Leiden cluster 3')
cell_clone_ids <- as.character(cell.clone[cells,'cloneid'])
good <- !is.na(cell_clone_ids) & nzchar(cell_clone_ids)
cells <- cells[good]
cell_clone_ids <- cell_clone_ids[good]

# Sparse aggregation: genes-by-cells %*% cells-by-clones.
clone_factor <- factor(cell_clone_ids, levels=unique(cell_clone_ids))
assign_mat <- Matrix::sparseMatrix(
  i=seq_along(cells), j=as.integer(clone_factor), x=rep(1,length(cells)),
  dims=c(length(cells), nlevels(clone_factor)),
  dimnames=list(cells, levels(clone_factor))
)
message('Aggregating ', length(cells), ' Leiden-3 cells into ', ncol(assign_mat), ' clones')
counts <- cd[,cells,drop=FALSE] %*% assign_mat
rm(assign_mat)

# Historical filters in Association_analysis/scripts/1_processing.R
clons <- intersect(colnames(counts), rownames(clone.info))
clone_totals <- as.numeric(Matrix::colSums(counts[,clons,drop=FALSE]))
clons <- clons[!is.na(clone.info[clons,'size']) & clone.info[clons,'size']>50 & clone_totals>0]
clons <- intersect(clons, names(rate))
if (length(clons)<5) stop('Too few clones after coverage/rate filtering')
rate1 <- as.numeric(rate[clons]); names(rate1) <- clons
if (anyNA(rate1) || !is.finite(stats::sd(rate1)) || stats::sd(rate1)==0) stop('Invalid clone rates')
zrate <- (rate1 - mean(rate1))/stats::sd(rate1)
clons <- clons[abs(zrate) < 1]
if (!length(clons)) stop('No clones after growth-rate outlier filter')

# Historical expression prevalence threshold is applied AFTER clone filtering
# and BEFORE fractions are normalized. This affects the denominator.
counts <- counts[, clons, drop=FALSE]
prevalence <- as.numeric(Matrix::rowSums(counts > 0))
names(prevalence) <- rownames(counts)
counts <- counts[prevalence>=20, , drop=FALSE]
if (!all(c('LSP1','NPM1') %in% rownames(counts))) {
  stop('Target gene missing after historical prevalence filter: ',
       paste(setdiff(c('LSP1','NPM1'),rownames(counts)),collapse=', '))
}
library_totals <- as.numeric(Matrix::colSums(counts))
if (any(!is.finite(library_totals)|library_totals<=0)) stop('Invalid pseudobulk totals')

# Extract only plotted genes; fractions correspond to the historical fpm table.
for (task in list(
  list(panel='Fig4E', gene='LSP1', phenotype='HSC bias'),
  list(panel='Fig4G', gene='NPM1', phenotype='growth rate')
)) {
  gene <- task$gene
  gene_counts <- as.numeric(counts[gene,clons])
  phenotype <- if (task$phenotype=='HSC bias') as.numeric(clone.states2[clons,'1']) else as.numeric(rate[clons])
  out <- data.frame(
    clone_id=clons,
    genotype=as.character(clone.info[clons,'genotype']),
    experiment=as.character(clone.info[clons,'expid']),
    sample_id=as.character(clone.info[clons,'id']),
    total_clone_size=as.numeric(clone.info[clons,'size']),
    cells_in_leiden_3=as.integer(table(factor(cell_clone_ids, levels=clons))),
    gene=gene,
    clone_phenotype=task$phenotype,
    phenotype_value=phenotype,
    gene_counts=gene_counts,
    denominator_counts=library_totals,
    relative_gene_expression=gene_counts/library_totals,
    stringsAsFactors=FALSE
  )
  if (anyNA(out$phenotype_value)|anyNA(out$relative_gene_expression)) stop('Missing plotting values: ',task$panel)
  output <- file.path(sd, paste0(task$panel,'_',gene,'_clone_expression.tsv'))
  write.table(out, output, sep='\t', row.names=FALSE, quote=FALSE, na='NA')
  message(task$panel, ': ', nrow(out), ' clones, ',gene,
          ' expression range ', signif(range(out$relative_gene_expression),3)[1],
          ' - ', signif(range(out$relative_gene_expression),3)[2])
}

summ <- do.call(rbind,lapply(c('Fig4E','Fig4G'), function(panel){
  gene <- if (panel=='Fig4E') 'LSP1' else 'NPM1'
  x <- read.delim(file.path(sd, paste0(panel,'_',gene,'_clone_expression.tsv')))
  data.frame(panel=panel,gene=gene,n_clones=nrow(x),
             n_WT=sum(x$genotype=='WT'),n_DNMT3A=sum(x$genotype=='DNMT3A'),
             n_TET2=sum(x$genotype=='TET2'),n_samples=length(unique(x$sample_id)),
             n_experiments=length(unique(x$experiment)),
             min_expression=min(x$relative_gene_expression),
             max_expression=max(x$relative_gene_expression))
}))
write.table(summ,file.path(vd,'Fig4EG_export_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(summ,row.names=FALSE)
message('Fig4E/G Source Data exported. Compare against accepted artwork before final validation.')
