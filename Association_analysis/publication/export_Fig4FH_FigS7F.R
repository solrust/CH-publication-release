# Reproduce candidate numerical Source Data for Fig. 4F/H and Supplementary Fig. 7F
# from the historical association checkpoint used for published Fig. 5C.
# Run from CH_publication_release repository root.
# Does not recompute regression coefficients, fgsea, or load scrna.RData.

pub <- 'Association_analysis/publication'
source_dir <- file.path(pub, 'source_data')
val_dir <- file.path(pub, 'validation')
dir.create(source_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(val_dir, showWarnings=FALSE, recursive=TRUE)

checkpoint <- 'Association_analysis/data/processed/default/association_summary.RData'
if (!file.exists(checkpoint)) stop('Missing checkpoint: ', checkpoint)
e <- new.env(parent=globalenv())
obj_names <- load(checkpoint, envir=e)
required <- c('test.basic.list', 'gsea.clone', 'mainPathways.clones')
if (!all(required %in% obj_names)) {
  stop('Checkpoint missing objects: ', paste(setdiff(required, obj_names), collapse=', '))
}

b <- e$test.basic.list
g <- e$gsea.clone$gsea.clone
mp <- e$mainPathways.clones
expected <- c('growth rate', 'clone size', 'HSC bias', 'lymphoid bias',
              'dendritic bias', 'myeloid bias', 'erythroid bias')
if (length(b) != 7 || length(g) != 7 || length(mp) != 7 ||
    !identical(names(b), expected)) {
  stop('Unexpected phenotype order or list sizes; inspect default checkpoint before joining.')
}
# The historical workflow assembled gsea.clone and mainPathways.clones by
# position from test.basic.list; both are therefore joined positionally.
if (!is.null(colnames(e$gsea.clone$gscore.clone)) &&
    length(colnames(e$gsea.clone$gscore.clone)) == 7 &&
    !identical(colnames(e$gsea.clone$gscore.clone), expected)) {
  stop('Stored gsea score column labels conflict with basic phenotype order.')
}

panels <- data.frame(
  panel = c('Fig4H', 'FigS7F_clone_size', 'Fig4F',
            'FigS7F_lymphoid_bias','FigS7F_dendritic_bias',
            'FigS7F_myeloid_bias','FigS7F_erythroid_bias'),
  phenotype = expected,
  stringsAsFactors=FALSE
)

write_tsv <- function(dat, file, compress=FALSE) {
  if (compress) {
    con <- gzfile(file, open='wt')
    on.exit(close(con), add=TRUE)
    write.table(dat, con, sep='\t', row.names=FALSE, col.names=TRUE,
                quote=FALSE, na='NA')
  } else {
    write.table(dat, file, sep='\t', row.names=FALSE, col.names=TRUE,
                quote=FALSE, na='NA')
  }
}

get_main_path_ids <- function(x, nm) {
  if (is.data.frame(x) && 'pathway' %in% colnames(x)) {
    return(as.character(x[['pathway']]))
  }
  if (is.character(x)) return(x)
  stop('Unexpected mainPathways format for ', nm, ': ', paste(class(x),collapse=', '))
}

metrics <- vector('list', length(expected))
for (i in seq_along(expected)) {
  phen <- expected[i]
  panel <- panels$panel[i]
  tab <- as.data.frame(b[[i]])
  if (!all(c('pheno','p.pheno') %in% names(tab))) {
    stop('Missing gene-level coefficients or p values for ', phen)
  }
  if (is.null(rownames(tab)) || anyDuplicated(rownames(tab))) {
    stop('Missing/nonunique gene IDs for ', phen)
  }
  pv <- as.numeric(tab[['p.pheno']]); co <- as.numeric(tab[['pheno']])
  if (any(is.finite(pv) & (pv < 0 | pv > 1))) stop('Invalid gene P values for ', phen)
  genes <- data.frame(
    panel=panel, phenotype=phen, gene=rownames(tab),
    coefficient_log_e=co, coefficient_log2=co/log(2),
    pvalue=pv, padj_BH=p.adjust(pv, method='BH'),
    label_candidate=is.finite(pv) & pv < 0.001,
    stringsAsFactors=FALSE
  )
  write_tsv(genes, file.path(source_dir, paste0(panel,'_gene_associations.tsv.gz')), TRUE)

  fg <- as.data.frame(g[[i]])
  if (!all(c('pathway','pval','padj','NES') %in% names(fg))) {
    stop('Unexpected fgsea columns for ', phen, ': ',paste(names(fg),collapse=', '))
  }
  pathway <- as.character(fg$pathway)
  if (anyDuplicated(pathway)) stop('Duplicate pathway identifiers for ', phen)
  p <- as.numeric(fg$pval)
  q <- as.numeric(fg$padj)
  nes <- as.numeric(fg$NES)
  if (any(is.finite(p) & (p < 0 | p > 1))) stop('Invalid enrichment P values for ',phen)
  main_ids <- unique(get_main_path_ids(mp[[i]],phen))
  if (any(!main_ids %in% pathway)) {
    warning(phen, ': ',sum(!main_ids %in% pathway),' archived main pathway IDs missing from fgsea')
  }
  historical_selected <- rep(FALSE,length(pathway))
  # Reconstruct selection logic in scripts/3_association_analysis.R (c. lines 168-180).
  # HSC bias was intentionally curated differently in the historical script.
  if (identical(phen,'HSC bias')) {
    historical_selected <- !is.na(q) & rank(q, na.last='keep') < 10 & nchar(pathway) < 60
    selection_rule <- 'HSC special: rank(fgsea_padj) < 10 and pathway name length < 60'
  } else {
    ix <- which(pathway %in% main_ids)
    if (length(ix)) {
      rr <- rank(q[ix], na.last='keep')
      historical_selected[ix] <- !is.na(rr) & rr < 15
    }
    selection_rule <- 'Archived mainPathways IDs, retaining first <15 ranks of fgsea_padj'
  }
  enriched <- data.frame(
    panel=panel, phenotype=phen, pathway=pathway,
    pvalue=p, padj=q, NES=nes,
    signed_minus_log10_p= -log10(p)*sign(nes),
    is_in_archived_main_pathways=pathway %in% main_ids,
    selected_for_historical_plot=historical_selected,
    stringsAsFactors=FALSE
  )
  if ('size' %in% names(fg)) enriched$gene_set_size <- as.numeric(fg$size)
  write_tsv(enriched, file.path(source_dir,paste0(panel,'_all_pathways.tsv.gz')),TRUE)
  write_tsv(enriched[historical_selected,,drop=FALSE],
            file.path(source_dir,paste0(panel,'_selected_pathways.tsv')))
  if (!any(historical_selected)) warning(panel, ': no pathways selected. Inspect archival lists.')
  metrics[[i]] <- data.frame(
    panel=panel, phenotype=phen, n_genes=nrow(genes),
    n_label_candidate_genes=sum(genes$label_candidate),
    n_fgsea_pathways=nrow(enriched),
    n_main_pathway_ids=length(main_ids),
    n_selected_pathways=sum(historical_selected),
    n_selected_with_nonfinite_scores=sum(!is.finite(enriched$signed_minus_log10_p[historical_selected])),
    pathway_selection_rule=selection_rule,
    input_checkpoint=checkpoint,
    stringsAsFactors=FALSE
  )
  message(panel,' (',phen,'): ',nrow(genes),' genes, ',nrow(enriched),
          ' pathways, ',sum(historical_selected),' selected bars')
}
summary <- do.call(rbind,metrics)
write_tsv(summary,file.path(val_dir,'Fig4FH_FigS7F_export_summary.tsv'))
write_tsv(panels,file.path(val_dir,'Fig4FH_FigS7F_panel_order.tsv'))
print(summary[,c('panel','phenotype','n_genes','n_selected_pathways')],row.names=FALSE)
message('Export complete. Selected pathway bars are candidates until matched to accepted artwork.')
