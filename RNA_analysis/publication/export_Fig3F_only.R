#!/usr/bin/env Rscript
# Publication Source Data: Fig. 3F only.
# Reads HISTORICAL mutant_de.RData; does not recalculate DESeq2 or fgsea or require S6D.
# Run from CH_publication_release repository root.

repo <- normalizePath(getwd(), mustWork = TRUE)
pub <- file.path(repo, 'RNA_analysis', 'publication')
if (!file.exists(file.path(repo, 'RNA_analysis', 'scripts', '6_mutant_analysis.R'))) {
  stop('Run from the CH_publication_release root (RNA_analysis/scripts not found).')
}

# CH_DATA_ROOT is optional, and should point to a directory containing
# RNA_analysis/data/processed, parallel to the code project layout.
data_root <- Sys.getenv('CH_DATA_ROOT', unset = '')
candidates <- c(
  if (nzchar(data_root)) file.path(data_root, 'RNA_analysis', 'data', 'processed', 'mutant_de.RData'),
  file.path(repo, 'RNA_analysis', 'data', 'processed', 'mutant_de.RData')
)
input <- candidates[file.exists(candidates)][1]
if (is.na(input)) stop('mutant_de.RData not found. Set CH_DATA_ROOT or put data under RNA_analysis/data/processed/')

out <- file.path(pub, 'source_data')
valdir <- file.path(pub, 'validation')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(valdir, recursive = TRUE, showWarnings = FALSE)

# Historical processed checkpoint.
h <- new.env(parent = baseenv())
load(input, envir = h)
if (!exists('fgsea.all.score', envir = h, inherits = FALSE)) {
  stop('Checkpoint is missing fgsea.all.score')
}
a <- as.matrix(h$fgsea.all.score)
stopifnot(is.numeric(a), ncol(a) == 20L, !anyNA(a), all(is.finite(a)),
          !anyDuplicated(rownames(a)), !anyDuplicated(colnames(a)))

# Reproduce the filtering and GO cluster derivation in
# RNA_analysis/scripts/6_mutant_analysis.R, lines 21-23 and 162-169.
fdr_across_pathways <- apply(a, 2, function(x) stats::p.adjust(exp(-abs(x)), method = 'BH'))
selected <- rowSums(fdr_across_pathways < 1e-5, na.rm = TRUE) > 0
stopifnot(any(selected), sum(selected) >= 6L)
score <- a[selected, , drop = FALSE]
cluster_tree <- stats::hclust(stats::dist(score))
go_cluster <- stats::cutree(cluster_tree, k = 6)

category_file <- file.path(pub, 'resources', 'Fig3F_functional_categories.tsv')
category <- utils::read.delim(category_file, stringsAsFactors = FALSE)
stopifnot(identical(category$GO_cluster, 1:6), !anyDuplicated(category$functional_category))
category_names <- setNames(category$functional_category, as.character(category$GO_cluster))

state <- data.frame(
  state_id = colnames(score),
  genotype = sub(' ([0-9]+)$', '', colnames(score)),
  cell_state = as.integer(sub('^.* ', '', colnames(score))),
  stringsAsFactors = FALSE
)
if (anyNA(state$cell_state) ||
    !all(state$genotype %in% c('DNMT3A', 'TET2')) ||
    !all(table(state$genotype) == 10L) ||
    !all(vapply(split(state$cell_state, state$genotype),
                function(x) setequal(x, 1:10), logical(1)))) {
  stop('Unexpected genotype/state column identifiers in fgsea.all.score')
}
state$original_matrix_column <- seq_len(nrow(state))
state$heatmap_row_position <- match(state$state_id,
                                   colnames(score)[stats::hclust(stats::as.dist(1 - stats::cor(score)))$order])

anno <- data.frame(
  pathway = rownames(score),
  GO_cluster = as.integer(go_cluster[rownames(score)]),
  functional_category = unname(category_names[as.character(go_cluster[rownames(score)])]),
  position_in_hierarchical_tree = match(rownames(score), cluster_tree$labels[cluster_tree$order]),
  stringsAsFactors = FALSE
)
stopifnot(!anyNA(anno$functional_category), !anyNA(anno$position_in_hierarchical_tree))

# Numerical matrix: original score = -ln(fgsea raw p-value) * sign(NES).
# No clipping here. Clipping to [-40,40] happens only in the heatmap display.
matrix_out <- data.frame(pathway = rownames(score), score, check.names = FALSE)
utils::write.table(matrix_out, file.path(out, 'Fig3F_pathway_score_matrix.tsv'),
                   sep = '\t', row.names = FALSE, quote = TRUE, na = 'NA')
utils::write.table(anno, file.path(out, 'Fig3F_pathway_annotations.tsv'),
                   sep = '\t', row.names = FALSE, quote = TRUE)
utils::write.table(state, file.path(out, 'Fig3F_state_annotations.tsv'),
                   sep = '\t', row.names = FALSE, quote = TRUE)


# Fig. 3F-only quality-control summary. S6D is deliberately out of scope.
qc <- data.frame(
  metric = c('archived_pathways', 'archived_comparisons',
             'Fig3F_selected_pathways', 'Fig3F_state_genotype_rows',
             'Fig3F_GO_clusters'),
  value = c(nrow(a), ncol(a), nrow(score), ncol(score),
            length(unique(anno$GO_cluster))),
  stringsAsFactors = FALSE
)
utils::write.table(qc, file.path(valdir, 'Fig3F_export_summary.tsv'),
                   sep='\t', row.names=FALSE, quote=TRUE)
message('Fig3F exported three numerical tables to RNA_analysis/publication/source_data/')
print(qc, row.names=FALSE)
message('S6D is NOT exported: see the analysis scope in RNA_analysis/README.md')
