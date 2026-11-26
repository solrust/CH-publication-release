# Compare the 16 Fig. 4H pathway IDs visible in the accepted Illustrator file
# against archived fgsea results exported from default/association_summary.RData.
# Does not modify Fig4H canonical Source Data or rerun fgsea.
# Run from the CH_publication_release repository root.

pub <- 'Association_analysis/publication'
sd <- file.path(pub, 'source_data')
vd <- file.path(pub, 'validation')
fd <- file.path(pub, 'figures')
dir.create(vd, showWarnings = FALSE, recursive = TRUE)
dir.create(fd, showWarnings = FALSE, recursive = TRUE)
ids_file <- file.path(pub, 'resources', 'Fig4H_accepted_pathway_ids.tsv')
all_file <- file.path(sd, 'Fig4H_all_pathways.tsv.gz')
selected_file <- file.path(sd, 'Fig4H_selected_pathways.tsv')
for (p in c(ids_file, all_file, selected_file)) if (!file.exists(p)) stop('Missing: ', p)

ids <- read.delim(ids_file, stringsAsFactors = FALSE, check.names = FALSE)
all <- read.delim(all_file, stringsAsFactors = FALSE, check.names = FALSE)
sel <- read.delim(selected_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(ids) == 16, !anyDuplicated(ids$pathway_id), !anyDuplicated(all$pathway))
required <- c('pathway', 'pvalue', 'padj', 'NES', 'signed_minus_log10_p')
stopifnot(all(required %in% names(all)))

# Match exact pathway IDs, not partial numeric strings or names.
path_ids <- sub("^(GO:[0-9]+|R-HSA-[0-9]+|M[0-9]+).*", "\\1", all$pathway)
if (anyDuplicated(path_ids)) {
  warning('Several fgsea pathways share the same leading ID; requiring exact name/ID review.')
}
match_list <- lapply(ids$pathway_id, function(pid) which(path_ids == pid))
nmatch <- lengths(match_list)
if (any(nmatch > 1)) stop('Ambiguous pathway ID matches: ', paste(ids$pathway_id[nmatch>1],collapse=', '))
ix <- vapply(match_list, function(m) if (length(m)) m else NA_integer_, integer(1))

report <- ids
report$found_in_checkpoint <- !is.na(ix)
report$checkpoint_pathway <- all$pathway[ix]
report$pvalue <- all$pvalue[ix]
report$padj_BH <- all$padj[ix]
report$NES <- all$NES[ix]
report$signed_minus_log10_p <- all$signed_minus_log10_p[ix]
report$selected_by_historical_code <- all$pathway[ix] %in% sel$pathway
report$score_finite <- is.finite(report$signed_minus_log10_p)
report$numeric_match_to_artwork <- 'NOT_TESTED: artwork does not encode a source table'
report_path <- file.path(vd, 'Fig4H_accepted_pathway_checkpoint_audit.tsv')
write.table(report, report_path, sep='\t', row.names=FALSE, quote=FALSE, na='NA')

cat('Accepted Fig4H IDs:', nrow(report), '\n')
cat('Found in archived fgsea:', sum(report$found_in_checkpoint), '\n')
cat('Also selected by original script:', sum(report$selected_by_historical_code,na.rm=TRUE), '\n')
cat('Found with finite plotted scores:', sum(report$score_finite,na.rm=TRUE), '\n')
print(report[,c('published_order','pathway_id','found_in_checkpoint',
                'selected_by_historical_code','pvalue','padj_BH','NES')], row.names=FALSE)
cat('Wrote: ',report_path,'\n')

if (all(report$found_in_checkpoint & report$score_finite)) {
  if (!requireNamespace('ggplot2',quietly=TRUE)) {
    message('Install ggplot2 to render the audit-only pathway review PDF.')
  } else {
    ptab <- report
    ptab$pathway_label <- factor(ptab$checkpoint_pathway, levels=rev(ptab$checkpoint_pathway))
    ptab$direction <- ifelse(ptab$NES > 0, 'positive', 'negative')
    pl <- ggplot2::ggplot(ptab, ggplot2::aes(x=pathway_label,y=signed_minus_log10_p,fill=direction)) +
      ggplot2::geom_col(color='grey35',linewidth=.2) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values=c(negative='#a9b9d9',positive='#e6aea7')) +
      ggplot2::labs(x=NULL,y='Signed -log10(fgsea P value)',
                    title='Fig. 4H: accepted artwork pathway IDs, archived fgsea scores',
                    subtitle='Audit-only review; confirm values against published bars before finalizing Source Data') +
      ggplot2::theme_bw(base_size=10) +
      ggplot2::theme(legend.position='none')
    outfile <- file.path(fd,'Fig4H_accepted_pathway_IDs_review.pdf')
    ggplot2::ggsave(outfile,plot=pl,width=11,height=8.5)
    cat('Created audit-only plot: ',outfile,'\n')
  }
} else {
  message('Some accepted pathways are missing/non-finite: see audit table; no reconstructed figure created.')
}
