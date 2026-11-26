# Export Fig. 5C association counts from the approved `use` results.
# Run from the CH_publication_release repository root.
# Does not recompute the association models or load scrna.RData.

outdir <- 'Association_analysis/publication'
sourcedir <- file.path(outdir, 'source_data')
validir <- file.path(outdir, 'validation')
dir.create(sourcedir, recursive=TRUE, showWarnings=FALSE)
dir.create(validir, recursive=TRUE, showWarnings=FALSE)

input <- 'Association_analysis/data/processed/default/association_summary.RData'
if (!file.exists(input)) stop('Missing final association checkpoint: ', input)
e <- new.env(parent=globalenv())
loaded <- load(input, envir=e)
message("Loaded association checkpoint ", input, ": ", paste(loaded, collapse = ", "))
need <- c('test.basic.list','test.inter.list')
if (!all(need %in% loaded)) stop('Missing required objects: ',paste(setdiff(need,loaded),collapse=', '))
basic <- e$test.basic.list
inter <- e$test.inter.list
if (length(basic) != length(inter)) stop('Basic and interaction list lengths differ')

# The historical scripts/2_processing.R builds both lists from the SAME
# ordered phenotype vector. A subsequent names(test.inter.list) assignment
# mistakenly indexes that already-shortened 7-element vector at 8 and 9,
# producing mislabelled/NA names. Repair ONLY that exact known pattern.
canonical <- c('growth rate', 'clone size', 'HSC bias', 'lymphoid bias',
               'dendritic bias', 'myeloid bias', 'erythroid bias')
historical_bad_names <- canonical[c(1, 2, 3, 5, 6, 8, 9)]
basic_names <- names(basic)
inter_names_original <- names(inter)
repair_reason <- 'none'
if (!identical(basic_names, inter_names_original)) {
  if (identical(basic_names, canonical) &&
      identical(inter_names_original, historical_bad_names)) {
    message('Recognized historical interaction-list name-assignment bug.')
    message('Restoring interaction names by original loop position (same phenotype order).')
    names(inter) <- basic_names
    repair_reason <- 'known_historical_names_bug'
  } else if (!is.null(basic_names) && !is.null(inter_names_original) &&
             !anyNA(basic_names) && !anyNA(inter_names_original) &&
             !anyDuplicated(basic_names) && !anyDuplicated(inter_names_original) &&
             setequal(basic_names, inter_names_original)) {
    message('Both lists have same phenotype names but different order; joining by exact name.')
    inter <- inter[match(basic_names, inter_names_original)]
    repair_reason <- 'exact_name_reorder'
  } else {
    message('Basic phenotype names: ', paste(basic_names, collapse=' | '))
    message('Interaction phenotype names: ', paste(inter_names_original, collapse=' | '))
    stop('Unknown phenotype-name mismatch: inspect before export; nothing was written.')
  }
}
pheno <- names(basic)
if (is.null(pheno) || anyNA(pheno) || any(!nzchar(pheno)) ||
    anyDuplicated(pheno)) stop('Invalid basic phenotype names; inspect checkpoint.')
stopifnot(identical(pheno, names(inter)))

# Keep an explicit audit of the original checkpoint names and the safe alignment.
name_audit <- data.frame(
  list_position = seq_along(pheno),
  basic_name = basic_names,
  interaction_name_original = inter_names_original,
  interaction_name_exported = pheno,
  reason = repair_reason,
  stringsAsFactors = FALSE
)
write.table(name_audit,
            file = file.path(validir, 'Fig5C_phenotype_name_alignment.tsv'),
            sep='\t', quote=FALSE, row.names=FALSE, na='NA')

term_map <- c('phenotype'='p.pheno',
              'phenotype:DNMT3A'='p.pheno:DNMT3A',
              'phenotype:TET2'='p.pheno:TET2')
threshold <- 0.5 # exact published/historical script; do not silently switch to 0.05

get_one <- function(tbl, phenotype, term, p_col) {
  if (is.null(dim(tbl)) || !p_col %in% colnames(tbl)) {
    stop('Missing p-value column ', p_col, ' for phenotype ', phenotype)
  }
  vals <- suppressWarnings(as.numeric(tbl[,p_col]))
  if (length(vals) != nrow(tbl)) stop('Bad p-value vector for ', phenotype, '/', term)
  if (any(is.finite(vals) & (vals < 0 | vals > 1))) stop('P-values outside [0,1]')
  adj <- p.adjust(vals, method='BH')
  genes <- rownames(tbl)
  if (is.null(genes)) genes <- paste0('row_', seq_len(nrow(tbl)))
  data.frame(phenotype=phenotype, term=term,
             gene_id=genes, pvalue=vals, padj_BH=adj,
             pass_FDR_0p5=!is.na(adj) & adj <= threshold,
             stringsAsFactors=FALSE)
}

records <- do.call(rbind, unlist(lapply(seq_along(pheno), function(i) {
  list(get_one(basic[[i]], pheno[i], 'phenotype', term_map['phenotype']),
       get_one(inter[[i]], pheno[i], 'phenotype:DNMT3A', term_map['phenotype:DNMT3A']),
       get_one(inter[[i]], pheno[i], 'phenotype:TET2', term_map['phenotype:TET2']))
}), recursive=FALSE))
rownames(records) <- NULL
z <- gzfile(file.path(sourcedir, 'Fig5C_gene_level_pvalues.tsv.gz'), 'wt')
tryCatch(write.table(records, file=z, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, na='NA'),
         finally=close(z))

counts <- aggregate(records$pass_FDR_0p5,
                    by=list(phenotype=records$phenotype, term=records$term), FUN=sum)
names(counts)[3] <- 'n_genes_BH_FDR_le_0p5'
counts$threshold <- threshold
counts$method <- 'BH adjusted within phenotype and term'
counts <- counts[order(match(counts$phenotype,pheno),match(counts$term,names(term_map))),]
write.table(counts, file=file.path(sourcedir, 'Fig5C_gene_association_counts.tsv'),
            sep='\t',quote=FALSE,row.names=FALSE)

audit <- data.frame(check=c('checkpoint','phenotypes','gene_test_rows','total_significant_entries','FDR_threshold'),
                    value=c(input,length(pheno),nrow(records),sum(records$pass_FDR_0p5),threshold))
write.table(audit,file=file.path(validir,'Fig5C_export_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(counts, row.names=FALSE)
message("Fig5C: exported ", nrow(counts), " summary rows and ", nrow(records), " gene-level records from ", input)
