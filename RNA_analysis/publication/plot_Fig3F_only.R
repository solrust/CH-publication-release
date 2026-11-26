#!/usr/bin/env Rscript
# Reproduces ONLY the BASE Fig.3F heatmap
# FROM exported TSV files only; no historical checkpoint required.
# Note final Illustrator figure contains added category labels and layout edits.

repo <- normalizePath(getwd(), mustWork = TRUE)
pub <- file.path(repo, 'RNA_analysis', 'publication')
sourcedata <- file.path(pub, 'source_data')
figures <- file.path(pub, 'figures')
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
for (p in c('pheatmap','RColorBrewer')) {
  if (!requireNamespace(p, quietly=TRUE)) stop('Please install R package: ', p)
}

matfile <- file.path(sourcedata, 'Fig3F_pathway_score_matrix.tsv')
annofile <- file.path(sourcedata, 'Fig3F_pathway_annotations.tsv')
statefile <- file.path(sourcedata, 'Fig3F_state_annotations.tsv')
stopifnot(all(file.exists(c(matfile, annofile, statefile))))

wide <- utils::read.delim(matfile, check.names=FALSE, stringsAsFactors=FALSE)
annotations <- utils::read.delim(annofile, stringsAsFactors=FALSE)
states <- utils::read.delim(statefile, stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(wide$pathway), !anyDuplicated(states$state_id),
          setequal(names(wide)[-1], states$state_id),
          setequal(wide$pathway, annotations$pathway))
score <- as.matrix(wide[,-1,drop=FALSE])
rownames(score) <- wide$pathway
stopifnot(is.numeric(score), all(is.finite(score)), ncol(score)==20)

# Restore the historical matrix ordering before re-applying the tree calculations.
score <- score[annotations$pathway,
               states$state_id[order(states$original_matrix_column)], drop=FALSE]
annotations <- annotations[match(rownames(score), annotations$pathway),,drop=FALSE]
states <- states[match(colnames(score), states$state_id),,drop=FALSE]

# Historic hclust, pheatmap parameters from 6_mutant_analysis.R, lines 162-199.
# The accepted Illustrator figure has manually added functional category names.
state_anno <- data.frame(genotype=states$genotype, row.names=states$state_id)
pathway_anno <- data.frame('GO.clusters'=factor(annotations$GO_cluster,levels=1:6),
                          row.names=annotations$pathway,
                          check.names=FALSE)
state_pal <- RColorBrewer::brewer.pal(10, 'Paired')
names(state_pal) <- as.character(1:10)
cluster_pal <- RColorBrewer::brewer.pal(6, 'Set2')
names(cluster_pal) <- as.character(1:6)

geno_pal <- c('DNMT3A'=grDevices::adjustcolor('red', alpha.f=.5),
              'TET2'=grDevices::adjustcolor('green', alpha.f=.5))

pheatmap::pheatmap(
  t(score),
  clustering_distance_rows=stats::as.dist(1-stats::cor(score)),
  clustering_distance_cols=stats::dist(score),
  breaks=seq(-40, 40, length.out=100),
  show_colnames=FALSE,
  treeheight_col=0,
  annotation_col=pathway_anno,
  labels_row=as.character(states$cell_state),
  annotation_row=state_anno,
  annotation_names_row=FALSE,
  annotation_colors=list('GO.clusters'=cluster_pal, genotype=geno_pal),
  filename=file.path(figures, 'Fig3F_reproduced_base.pdf'),
  width=8, height=3
)


capture.output(utils::sessionInfo(), file=file.path(figures, 'Fig3F_sessionInfo.txt'))
cat('Created RNA_analysis/publication/figures/Fig3F_reproduced_base.pdf\n')
cat('Supplementary Fig. 6D remains on hold pending pathway provenance reconciliation.\n')
