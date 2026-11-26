# S8A/B: reconstruct the published clone-phenotype matrix without rerunning models.
# Root of CH_publication_release is the working directory.
# Reproduces Association_analysis/scripts/1_processing.R phenotype definitions
# and scripts/9_clone_phenotype_similarity.R filtering and statistics.

pub <- 'Association_analysis/publication'
sd <- file.path(pub,'source_data'); val <- file.path(pub,'validation')
dir.create(sd,recursive=TRUE,showWarnings=FALSE)
dir.create(val,recursive=TRUE,showWarnings=FALSE)
cl_file <- 'RNA_clone_analysis/data/processed/clone_info.RData'
dna_file <- 'DNA_analysis/data/processed/DNA_clones.RData'
for(f in c(cl_file,dna_file)) if(!file.exists(f)) stop('Missing required checkpoint: ',f)
e <- new.env(parent=baseenv()); e2 <- new.env(parent=baseenv())
load(cl_file,envir=e); load(dna_file,envir=e2)
if(!all(c('clone.info','clone.states2') %in% ls(e))) stop('Clone checkpoint missing expected objects')
if(!exists('rate',envir=e2,inherits=FALSE)) stop('DNA checkpoint does not contain named growth-rate vector `rate`; inspect names of loaded objects: ',paste(ls(e2),collapse=', '))
inf <- e$clone.info; states <- as.matrix(e$clone.states2); rate <- e2$rate
if(is.null(names(rate))) stop('DNA growth-rate `rate` needs clone ID names')
if(!all(c('size','id') %in% names(inf))) stop('clone.info missing size or id')
if(!all(c('1','4','5','8','10') %in% colnames(states))) stop('Missing expected state columns')
clons <- rownames(states)
if(!all(clons %in% rownames(inf))) stop('Clone metadata missing clone.state IDs')
info <- inf[clons,,drop=FALSE]
ids <- as.character(info$id)
size <- as.numeric(info$size)
if(anyNA(size)|any(size<=0)) stop('Missing or non-positive clone sizes')
# The upstream code uses aggregate(size ~ sample ID) and log(clone size / sample total).
sample_sizes <- tapply(size,ids,sum)
ph <- data.frame(clone_id=clons, sample_id=ids,
  genotype=as.character(info$genotype),n_cells_in_clone=size,
  'growth rate'=as.numeric(rate[clons]),
  'clone size'=log(size/as.numeric(sample_sizes[ids])),
  'HSC bias'=as.numeric(states[clons,'1']),
  'lymphoid bias'=as.numeric(states[clons,'4']),
  'dendritic bias'=as.numeric(states[clons,'5']),
  'myeloid bias'=as.numeric(states[clons,'8']),
  'erythroid bias'=as.numeric(states[clons,'10']),
  check.names=FALSE,stringsAsFactors=FALSE)
needed <- c('growth rate','clone size','HSC bias','lymphoid bias','dendritic bias','myeloid bias','erythroid bias')
# IMPORTANT: DNA_analysis/scripts/1_DNA_processing.R estimates growth rates only
# for DNA clones with >=100 barcode counts across the first five time points.
# Consequently, `rate` intentionally has fewer entries than `clone.info`.
# Original association processing filters to clones with an available rate.
# Match by EXACT clone ID; do not fill missing rates or align by row order.
if (anyDuplicated(names(rate))) stop('Duplicated growth-rate clone IDs in DNA checkpoint')
n_matches <- sum(clons %in% names(rate))
if (n_matches == 0L) {
  stop('No overlapping clone IDs between clone.info and DNA growth rates. ',
       'First RNA IDs: ', paste(head(clons,3),collapse=', '),
       '; first DNA IDs: ', paste(head(names(rate),3),collapse=', '))
}
if (n_matches < 30L) {
  stop('Only ', n_matches, ' clone IDs match the DNA growth-rate vector; ',
       'inspect ID conventions before exporting.')
}
# scripts/9_clone_phenotype_similarity.R excludes clones with growth-rate
# magnitude >1 and <=50 cells. It implicitly drops absent rates from
# complete-case correlations; here we make the omissions explicit.
rate_available <- clons %in% names(rate)
rate_finite <- is.finite(ph[['growth rate']])
other_finite <- rowSums(is.finite(as.matrix(ph[,needed[-1],drop=FALSE]))) == length(needed)-1L
reason <- rep('retained',length(clons))
reason[!rate_available] <- 'no_DNA_growth_rate'
reason[rate_available & !rate_finite] <- 'nonfinite_DNA_growth_rate'
reason[reason=='retained' & ph$n_cells_in_clone <= 50] <- 'clone_size_le_50'
reason[reason=='retained' & abs(ph[['growth rate']]) > 1] <- 'growth_rate_abs_gt_1'
reason[reason=='retained' & !other_finite] <- 'nonfinite_other_phenotype'
# Write an explicit audit even if the retained-clone count is unexpected.
inc <- data.frame(clone_id=clons,sample_id=ids,genotype=as.character(info$genotype),
                  size=size,matched_DNA_rate=rate_available,
                  growth_rate=ph[['growth rate']],inclusion_status=reason,
                  stringsAsFactors=FALSE)
write.table(inc,file.path(val,'FigS8AB_clone_inclusion_audit.tsv'),
            sep='\t',quote=FALSE,row.names=FALSE,na='NA')
cat('S8A/B exact clone ID overlap: ',n_matches,' / ',length(clons),' clones have stored DNA growth rates.\n',sep='')
print(as.data.frame(table(inclusion_status=reason),stringsAsFactors=FALSE),row.names=FALSE)
keep <- reason=='retained'
ph <- ph[keep,,drop=FALSE]
if(nrow(ph)<30L) stop('Unexpectedly few retained clones: ',nrow(ph),
                     '. See validation/FigS8AB_clone_inclusion_audit.tsv')
if(anyDuplicated(ph$clone_id)) stop('Nonunique clone IDs')
con <- gzfile(file.path(sd,'FigS8AB_clone_phenotype_values.tsv.gz'),'wt')
tryCatch(write.table(ph,con,sep='\t',quote=FALSE,row.names=FALSE),finally=close(con))
# S8A correlations use log(bias + 1), while S8B uses the raw values.
mat <- ph[,needed,drop=FALSE]
mat[,3:7] <- log(mat[,3:7]+1)
if(any(!is.finite(as.matrix(mat)))) stop('Nonfinite log-transformed phenotype values')
rmat <- stats::cor(mat,method='pearson',use='na.or.complete')
ord <- expand.grid(y=needed,x=needed,stringsAsFactors=FALSE)
ord$r_pearson <- mapply(function(y,x) rmat[y,x],ord$y,ord$x)
ord$diagonal <- ord$x==ord$y
write.table(ord,file.path(sd,'FigS8A_correlations.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# S8B: ten selected pairs from 9_clone_phenotype_similarity.R original 30-pair grid.
combos <- rbind(data.frame(y='growth rate', x=c('clone size','HSC bias','erythroid bias','myeloid bias','lymphoid bias')),
                data.frame(y='HSC bias', x=c('growth rate','clone size','erythroid bias','myeloid bias','lymphoid bias')))
rows <- lapply(seq_len(nrow(combos)),function(i){
  x <- combos$x[i]; y <- combos$y[i]
  data.frame(clone_id=ph$clone_id,sample_id=ph$sample_id,genotype=ph$genotype,
             x_phenotype=x,y_phenotype=y,x_value=ph[[x]],y_value=ph[[y]],
             pearson_r=unname(stats::cor(ph[[x]],ph[[y]],method='pearson')),
             stringsAsFactors=FALSE)
})
points <- do.call(rbind,rows)
con <- gzfile(file.path(sd,'FigS8B_clone_phenotype_pairs.tsv.gz'),'wt')
tryCatch(write.table(points,con,sep='\t',quote=FALSE,row.names=FALSE),finally=close(con))
summary <- data.frame(input_clones=cl_file,input_DNA=dna_file,
  n_clone_info=length(clons), n_matched_DNA_growth_rates=n_matches,
  n_missing_DNA_growth_rates=sum(reason=='no_DNA_growth_rate'),
  n_excluded_nonfinite_growth_rate=sum(reason=='nonfinite_DNA_growth_rate'),
  n_excluded_size_le_50=sum(reason=='clone_size_le_50'),
  n_excluded_growth_abs_gt_1=sum(reason=='growth_rate_abs_gt_1'),
  n_excluded_other_nonfinite=sum(reason=='nonfinite_other_phenotype'),
  n_clones=nrow(ph), n_samples=length(unique(ph$sample_id)),
  n_WT=sum(ph$genotype=='WT'), n_DNMT3A=sum(ph$genotype=='DNMT3A'),
  n_TET2=sum(ph$genotype=='TET2'), n_pairs=nrow(combos),
  note='S8A correlations use log(bias+1); S8B scatterplots use raw phenotype values',
  stringsAsFactors=FALSE)
write.table(summary,file.path(val,'FigS8AB_export_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
cat('S8A/B: ',nrow(ph),' retained clones from ',length(unique(ph$sample_id)),' samples.\n',sep='')
cat('S8B ten selected Pearson correlations:\n')
print(data.frame(combos,r=vapply(rows,function(x)x$pearson_r[1],numeric(1))),row.names=FALSE)
