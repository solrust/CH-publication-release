# Figure 5B: reproduce genotype/AML signature associations with clone phenotypes.
# Run from CH_publication_release repository root. No model refitting.
# The source code calls wilcox.test then OVERWRITES that result with t.test.
# Therefore the published test is Welch's two-sided two-sample t-test.

pub <- 'Association_analysis/publication'
sd_dir <- file.path(pub, 'source_data')
val_dir <- file.path(pub, 'validation')
dir.create(sd_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(val_dir, recursive=TRUE, showWarnings=FALSE)

ass_file <- 'Association_analysis/data/processed/default/association_summary.RData'
de_file <- 'RNA_analysis/data/processed/mutant_de.RData'
aml_file <- 'RNA_analysis/data/genesets/processed/aml_sets.rds'
for (file in c(ass_file, de_file, aml_file)) if (!file.exists(file)) stop('Missing input: ', file)

a <- new.env(parent=globalenv())
loaded_ass <- load(ass_file, envir=a)
if (!('gsea.clone' %in% loaded_ass) || !('test.basic.list' %in% loaded_ass)) stop('Association checkpoint lacks gsea.clone or test.basic.list')
score <- as.matrix(a$gsea.clone$gscore.clone)
if (is.null(rownames(score)) || anyDuplicated(rownames(score))) stop('Gene-level scores need unique gene names')
phenos <- names(a$test.basic.list)
if (length(phenos)!=7L || anyNA(phenos) || anyDuplicated(phenos)) stop('Unexpected phenotype count/names')
if (ncol(score)!=length(phenos)) stop('Gene score matrix does not match phenotype list')
colnames(score) <- phenos
wanted <- c('growth rate','clone size','HSC bias','lymphoid bias','dendritic bias','myeloid bias','erythroid bias')
if (!setequal(phenos,wanted)) stop('Phenotype labels do not match expected seven: ', paste(phenos, collapse=', '))
score <- score[,wanted,drop=FALSE]
rm(a); gc(verbose=FALSE)

d <- new.env(parent=globalenv())
loaded_de <- load(de_file, envir=d)
if (!('de.mut2wt.shrink' %in% loaded_de)) stop('DE checkpoint lacks de.mut2wt.shrink')
if (!length(d$de.mut2wt.shrink)) stop('Empty HSC DE results')
first <- d$de.mut2wt.shrink[[1]]
if (!all(c('DNMT3A','TET2') %in% names(first))) stop('Unexpected genotype names in HSC DE results')

make_de_membership <- function(res, signature, padj_cut=0.3) {
  vals <- as.data.frame(res)
  if (!all(c('padj','log2FoldChange') %in% names(vals))) stop(signature,' DE lacks padj or log2FoldChange')
  if (is.null(rownames(vals))) stop(signature,' DE lacks gene names')
  ok <- is.finite(vals$padj) & vals$padj < padj_cut & is.finite(vals$log2FoldChange)
  up <- rownames(vals)[ok & vals$log2FoldChange >= 0]
  down <- rownames(vals)[ok & vals$log2FoldChange <= 0]
  list(up=up,down=down)
}
sig <- list(
  'DNMT3A' = make_de_membership(first$DNMT3A,'DNMT3A'),
  'TET2' = make_de_membership(first$TET2,'TET2')
)
rm(d,first); gc(verbose=FALSE)

aml <- readRDS(aml_file)
if (is.null(aml$aml_sign)) stop('aml_sets.rds does not contain aml_sign')
am <- as.data.frame(aml$aml_sign)
if (is.null(rownames(am)) || !all(c('log2FC','logP.sum','n.patients.sign') %in% colnames(am))) {
  stop('Unexpected AML signature columns; inspect aml_sets.rds')
}
keep_aml <- is.finite(am$log2FC) & is.finite(am$logP.sum) & is.finite(am$n.patients.sign) &
  am$logP.sum > -log(1e-5) & am$n.patients.sign >= 2
sig[['AML']] <- list(
  up=rownames(am)[keep_aml & am$log2FC > 0],
  down=rownames(am)[keep_aml & am$log2FC < 0]
)

# The original code assigns 'up' and then 'down'; 'down' wins any overlap.
# We keep the gene-row order of gsea.clone$gscore.clone so the historical
# first-observed group can be audited independently of the explicit up/down sign.
rows <- list(); coverage <- list(); j <- 0L; cidx <- 0L
for (signature in names(sig)) {
  m <- stats::setNames(rep(NA_character_,nrow(score)),rownames(score))
  u <- unique(as.character(sig[[signature]]$up))
  dn <- unique(as.character(sig[[signature]]$down))
  for (group in c('up','down')) {
    ids <- intersect(as.character(sig[[signature]][[group]]), names(m))
    m[ids] <- group
  }
  cidx <- cidx+1L
  coverage[[cidx]] <- data.frame(signature=signature, n_original_up=length(u),
    n_original_down=length(dn),n_overlap_up_down=length(intersect(u,dn)),
    n_scored_up=sum(m=='up', na.rm=TRUE),n_scored_down=sum(m=='down', na.rm=TRUE),stringsAsFactors=FALSE)
  kept <- which(!is.na(m))
  for (ph in wanted) {
    j <- j+1L
    rows[[j]] <- data.frame(gene=names(m)[kept],phenotype=ph,signature=signature,
      perturbation_direction=m[kept],association_score=as.numeric(score[kept,ph]),
      stringsAsFactors=FALSE)
  }
}
genes <- do.call(rbind,rows)
if (!nrow(genes)) stop('No signature genes matched the association score matrix')
write_gz <- function(x, dest) {
  con <- gzfile(dest,'wt')
  tryCatch(write.table(x,file=con,sep='\t',quote=FALSE,row.names=FALSE,na='NA'),
           finally=close(con))
}
write_gz(genes,file.path(sd_dir,'Fig5B_signature_gene_association_scores.tsv.gz'))
write.table(do.call(rbind,coverage),file.path(val_dir,'Fig5B_signature_coverage.tsv'),
            sep='\t',quote=FALSE,row.names=FALSE)

# Both the original first-observed-group direction and explicit up-minus-down
# direction are exported: the historical script depended on row order.
stat_rows <- list(); q <- 0L
for (ph in wanted) for (signature in names(sig)) {
  ss <- genes[genes$phenotype==ph & genes$signature==signature &
                is.finite(genes$association_score),,drop=FALSE]
  u <- ss$association_score[ss$perturbation_direction=='up']
  dn <- ss$association_score[ss$perturbation_direction=='down']
  if (length(u)<2L || length(dn)<2L) stop('Too few scored genes for ',ph,'/',signature)
  test <- stats::t.test(u,dn,alternative='two.sided') # Default Welch variances
  p <- test$p.value
  if (!is.finite(p)) stop('Invalid Welch P for ',ph,'/',signature)
  first_group <- unique(ss$perturbation_direction)[1]
  up_minus_down <- mean(u)-mean(dn)
  sign_up <- sign(up_minus_down)
  sign_historical <- if(first_group=='up') sign_up else -sign_up
  neg_log <- -log10(max(p,.Machine$double.xmin))
  q <- q+1L
  stat_rows[[q]] <- data.frame(phenotype=ph,signature=signature,
    n_genes_up=length(u),n_genes_down=length(dn),
    mean_up=mean(u),mean_down=mean(dn),mean_up_minus_down=up_minus_down,
    welch_t_two_sided_pvalue=p,df_welch=unname(test$parameter),
    signed_neglog10p_up_minus_down=sign_up*neg_log,
    first_group_in_archived_gene_order=first_group,
    signed_neglog10p_historical=sign_historical*neg_log,
    p_adjusted=FALSE,
    significance_stars=if(p < 5e-4) '***' else if(p < 1e-3) '**' else if(p < 0.05) '*' else '',
    stringsAsFactors=FALSE)
}
stats_tbl <- do.call(rbind,stat_rows)
if (nrow(stats_tbl)!=21L) stop('Expected exactly 7 x 3=21 heatmap cells')
write.table(stats_tbl,file.path(sd_dir,'Fig5B_heatmap_statistics.tsv'),
            sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(source=c('association','HSC_DE','AML_signatures'),
  relative_path=c(ass_file,de_file,aml_file)),
  file.path(val_dir,'Fig5B_input_files.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
message('Fig5B: ',nrow(stats_tbl),' heatmap cells; ',length(unique(genes$gene)),' distinct scored genes.')
message('Historical script: Welch t.test overwrites wilcox.test; p-values unadjusted.')
message('NOTE: verify both historical sign convention and stars against accepted Fig5B.')
print(stats_tbl[,c('phenotype','signature','n_genes_up','n_genes_down',
 'welch_t_two_sided_pvalue','first_group_in_archived_gene_order',
 'signed_neglog10p_historical','significance_stars')],row.names=FALSE)
