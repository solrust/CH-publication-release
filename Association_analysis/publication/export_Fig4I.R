# Fig. 4I: pathway association scores for clone growth rate vs HSC bias.
# Run from CH_publication_release root. Uses FINAL-FIGURE-COMPATIBLE
# Association_analysis/data/processed/default/association_summary.RData.
# This script does not recompute fgsea or gene association models.

pub <- 'Association_analysis/publication'
sd_dir <- file.path(pub, 'source_data')
val_dir <- file.path(pub, 'validation')
dir.create(sd_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(val_dir, recursive=TRUE, showWarnings=FALSE)
input <- 'Association_analysis/data/processed/default/association_summary.RData'
if (!file.exists(input)) stop('Missing association checkpoint: ',input)
e <- new.env(parent=emptyenv())
loaded <- load(input, envir=e)
if (!'gsea.clone' %in% loaded || !'test.basic.list' %in% loaded) {
  stop('Checkpoint lacks gsea.clone or test.basic.list')
}
if (is.null(e$gsea.clone$gsea.score)) stop('Missing gsea.clone$gsea.score')
score <- as.matrix(e$gsea.clone$gsea.score)
phenotypes <- names(e$test.basic.list)
if (is.null(phenotypes) || anyNA(phenotypes) || anyDuplicated(phenotypes) ||
    ncol(score) != length(phenotypes) ||
    !all(c('growth rate','HSC bias') %in% phenotypes)) {
  stop('Unexpected phenotype order: inspect default checkpoint before aligning matrix columns')
}
# Historical 3_association_analysis.R explicitly assigns these names by position.
colnames(score) <- phenotypes
if (is.null(rownames(score)) || anyNA(rownames(score)) || anyDuplicated(rownames(score))) {
  stop('Missing/duplicated pathway names')
}

out <- data.frame(pathway=rownames(score),stringsAsFactors=FALSE, check.names=FALSE)
out$growth_rate_score <- as.numeric(score[,'growth rate'])
out$HSC_bias_score <- as.numeric(score[,'HSC bias'])
# Preserve the other scores because the historical GO-category clustering
# selected pathways by variation across ALL seven clone phenotypes.
for (p in setdiff(phenotypes,c('growth rate','HSC bias'))) {
  col <- paste0('score_', gsub('[^A-Za-z0-9]+','_',p))
  out[[col]] <- as.numeric(score[,p])
}
out$finite_for_plot <- is.finite(out$growth_rate_score) & is.finite(out$HSC_bias_score)
out$inside_historical_limits <- out$finite_for_plot &
  out$growth_rate_score >= -30 & out$growth_rate_score <= 30 &
  out$HSC_bias_score >= -13 & out$HSC_bias_score <= 13

out_file <- file.path(sd_dir,'Fig4I_pathway_association_scores.tsv.gz')
con <- gzfile(out_file,'wt')
tryCatch(write.table(out,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA'),
         finally=close(con))

summary <- data.frame(
  input_checkpoint=input,
  n_pathways=nrow(out),
  n_finite=sum(out$finite_for_plot),
  n_in_accepted_axis_range=sum(out$inside_historical_limits),
  n_outside_axis_range=sum(out$finite_for_plot & !out$inside_historical_limits),
  n_phenotypes_used_for_category_selection=ncol(score),
  score_definition='signed -ln(fgsea unadjusted p-value), with sign of ES/NES',
  note='GO categories reconstructed separately; original Leiden labels not stored',
  stringsAsFactors=FALSE
)
write.table(summary,file.path(val_dir,'Fig4I_export_summary.tsv'),
            sep='\t',row.names=FALSE,quote=FALSE)
message('Fig4I: ',nrow(out),' pathways; ',sum(out$finite_for_plot),
        ' finite; ',sum(out$inside_historical_limits),' inside original display range.')
message('Wrote ',out_file)
