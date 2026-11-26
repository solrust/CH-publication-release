# Replace historical 14-pathway Fig4H plot selection with the 16 published IDs.
# Read accepted-ID audit, generated from the default/ association checkpoint.
# Run from CH_publication_release root, after:
#   source('Association_analysis/publication/export_Fig4FH_FigS7F.R')
#   source('Association_analysis/publication/audit_Fig4H_accepted_pathways.R')
# Only a small source-data/validation edit; no fgsea recalculation.

pub <- 'Association_analysis/publication'
sd <- file.path(pub, 'source_data')
vd <- file.path(pub, 'validation')
audit_file <- file.path(vd, 'Fig4H_accepted_pathway_checkpoint_audit.tsv')
all_file <- file.path(sd, 'Fig4H_all_pathways.tsv.gz')
selected_file <- file.path(sd, 'Fig4H_selected_pathways.tsv')
old_archive <- file.path(vd, 'Fig4H_historical_14_pathway_selection.tsv')
for (f in c(audit_file, all_file)) if (!file.exists(f)) stop('Missing input: ', f)

x <- read.delim(audit_file, stringsAsFactors=FALSE, check.names=FALSE)
y <- read.delim(all_file, stringsAsFactors=FALSE, check.names=FALSE)
needed <- c('published_order', 'pathway_id', 'checkpoint_pathway',
            'found_in_checkpoint', 'score_finite', 'pvalue', 'padj_BH',
            'NES', 'signed_minus_log10_p')
stopifnot(all(needed %in% names(x)), nrow(x) == 16,
          identical(sort(as.integer(x$published_order)), 1:16),
          !anyDuplicated(x$pathway_id), !anyDuplicated(x$checkpoint_pathway),
          all(x$found_in_checkpoint), all(x$score_finite),
          all(is.finite(x$pvalue)), all(is.finite(x$padj_BH)),
          all(is.finite(x$NES)), all(is.finite(x$signed_minus_log10_p)),
          !anyDuplicated(y$pathway))

x <- x[order(x$published_order), , drop=FALSE]
ix <- match(x$checkpoint_pathway, y$pathway)
if (anyNA(ix)) stop('Accepted pathway absent from full fgsea export')
for (pair in list(c('pvalue','pvalue'), c('padj_BH','padj'),
                  c('NES','NES'), c('signed_minus_log10_p','signed_minus_log10_p'))) {
  a <- as.numeric(x[[pair[1]]]); b <- as.numeric(y[[pair[2]]][ix])
  if (!isTRUE(all.equal(a, b, tolerance=1e-10, check.attributes=FALSE))) {
    stop('Audit and complete fgsea export disagree for ', pair[1])
  }
}

# Preserve the same fields required by plot_Fig4FH_FigS7F.R, and extra
# fields documenting editorial panel selection and original published order.
out <- y[ix, , drop=FALSE]
out$selected_for_published_panel <- TRUE
out$published_order <- as.integer(x$published_order)
out$published_pathway_id <- as.character(x$pathway_id)
out$selection_provenance <- 'Published Fig4H pathway IDs; default association checkpoint fgsea values'

# Save the old automatic selection under validation/, not source_data/.
# Avoid overwriting a conflicting historical archive.
if (file.exists(selected_file)) {
  old <- read.delim(selected_file, stringsAsFactors=FALSE)
  if (nrow(old) != 16 || !('published_pathway_id' %in% names(old))) {
    if (file.exists(old_archive)) {
      if (!identical(unname(tools::md5sum(selected_file)), unname(tools::md5sum(old_archive)))) {
        stop('Historical archive already exists with different contents; inspect before overwriting')
      }
    } else if (!file.copy(selected_file, old_archive, overwrite=FALSE)) {
      stop('Could not archive historical selection at ', old_archive)
    }
    message('Archived historical selection under validation/: ', old_archive)
  }
}
# Atomic replacement on the same filesystem; write a temporary file first.
tmp <- paste0(selected_file, '.tmp')
write.table(out, tmp, sep='\t', row.names=FALSE, quote=FALSE, na='NA')
if (!file.rename(tmp, selected_file)) stop('Could not replace ', selected_file)
message('Canonical Fig4H selected pathway Source Data: ', selected_file)
message('Published pathways: ', nrow(out), ' (from the archived default/ checkpoint)')
message('Note: rerunning export_Fig4FH_FigS7F.R overwrites selection with its historical 14; rerun this curation script afterward.')
